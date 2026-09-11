#!/usr/bin/env python3
"""Post-processing for a pipeline-level vContact3 run's final_assignments.csv.

Given vContact3's final_assignments.csv and the same run's combined
gene2genome.tsv (protein_id, genome_id, keywords — see vcontact3_prep.py),
produce final_assignments_postprocessed.csv where:

  1. Only query genomes are retained (rows whose Genome column contains '||'
     — this pipeline's sample-namespacing separator, as emitted by
     vcontact3_prep.py; reference genomes from the database never contain it).
  2. Each query row's Proteins count is filled in from gene2genome.tsv
     (counting protein rows per genome_id) wherever vContact3's own Proteins
     column is blank. vConTACT3 3.2.4 only populates Proteins (and
     GenomeName, Size_Kb) from its reference-DB metadata table, so every
     query genome's Proteins is blank in final_assignments.csv — comparing
     a blank/NaN Proteins against --min-proteins/--max-proteins is silently
     always False, so this step is required before step 3 can do anything.
  3. A query row's genus_prediction is flagged as uncertain when it is a
     novel-genus call (genus_prediction starts with 'novel_genus') AND the
     genome has too few or too many proteins for that call to be trusted
     (Proteins < --min-proteins or > --max-proteins). Flagged rows get an
     'uncertain_novel_' prefix prepended to genus_prediction.
  4. The result is sorted by Genome (alphabetical) and written to --output.
  5. The remaining novel-taxa calls — genus_prediction still starts with
     'novel_genus' after step 3 (i.e. NOT demoted to 'uncertain_novel_') —
     are additionally written out on their own to --novel-taxa-output, for
     easy downstream review of the genomes confident enough to flag as
     candidate novel genera.

Usage:
    vcontact3_postprocess.py \\
        --assignments       vcontact3_out/exports/final_assignments.csv \\
        --gene2genome       combined_gene2genome.tsv \\
        --output            vcontact3_out/exports/final_assignments_postprocessed.csv \\
        --novel-taxa-output vcontact3_out/exports/final_assignments_noveltaxa.csv \\
        [--min-proteins 5] [--max-proteins 20]
"""

import argparse
import sys

import pandas as pd

# Sample-namespacing separator vcontact3_prep.py builds query genome IDs with
# (`<sample_id>||<bin_or_scaffold>`). Reference genomes from the DB never contain it.
QUERY_SEPARATOR = '||'


def filter_query_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Keep only query genomes (Genome contains '||'); drop reference genomes.

    regex=False matters: '||' as a regex is an empty alternation, which matches
    every row and would silently let the whole reference DB through."""
    return df[df['Genome'].str.contains(QUERY_SEPARATOR, na=False, regex=False)].copy()


def load_protein_counts(gene2genome_tsv: str) -> pd.Series:
    """genome_id -> number of proteins, counted from a gene2genome TSV (protein_id, genome_id, keywords)."""
    g2g = pd.read_csv(gene2genome_tsv, sep='\t')
    return g2g.groupby('genome_id').size()


def fill_protein_counts(df: pd.DataFrame, protein_counts: pd.Series) -> pd.DataFrame:
    """Fill df['Proteins'] from gene2genome-derived counts wherever vContact3's own
    value is missing (query genomes: vConTACT3 3.2.4 leaves Proteins blank for these,
    populating it only for reference genomes from its DB metadata)."""
    df = df.copy()
    df['Proteins'] = pd.to_numeric(df['Proteins'], errors='coerce')
    computed = df['Genome'].map(protein_counts)
    df['Proteins'] = df['Proteins'].where(df['Proteins'].notna(), computed)
    return df


def assert_proteins_known(df: pd.DataFrame) -> None:
    """Fail loudly if any query genome still has no Proteins count after falling
    back to gene2genome.tsv. A silent NaN here is exactly the bug this guards
    against: pandas comparisons against NaN are silently False, so an unnoticed
    gap would make every novel-genus call in that gap pass through unflagged."""
    missing = df.loc[df['Proteins'].isna(), 'Genome'].tolist()
    if missing:
        sys.exit(
            f'ERROR: no protein count available (neither in final_assignments.csv nor '
            f'gene2genome.tsv) for {len(missing)} query genome(s): {", ".join(missing)}'
        )


def flag_uncertain_novel_genus(df: pd.DataFrame, min_proteins: int, max_proteins: int) -> pd.DataFrame:
    """Prefix genus_prediction with 'uncertain_novel_' for novel-genus calls
    on genomes with too few or too many proteins to trust that call."""
    is_novel = df['genus_prediction'].str.startswith('novel_genus', na=False)
    out_of_range = (df['Proteins'] < min_proteins) | (df['Proteins'] > max_proteins)
    flag = is_novel & out_of_range

    df.loc[flag, 'genus_prediction'] = 'uncertain_novel_' + df.loc[flag, 'genus_prediction']
    return df


def extract_novel_taxa(df: pd.DataFrame) -> pd.DataFrame:
    """Rows whose genus_prediction is still a (non-demoted) novel-genus call."""
    return df[df['genus_prediction'].str.startswith('novel_genus', na=False)].copy()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--assignments', required=True, help='final_assignments.csv from a vcontact3 run')
    parser.add_argument('--gene2genome', required=True,
                         help='combined gene2genome.tsv (protein_id, genome_id, keywords) from the same run, '
                              'used to fill in per-genome Proteins counts vConTACT3 leaves blank for query genomes')
    parser.add_argument('--output', required=True, help='where to write final_assignments_postprocessed.csv')
    parser.add_argument('--novel-taxa-output', required=True, help='where to write final_assignments_noveltaxa.csv (the still-novel-genus subset)')
    parser.add_argument('--min-proteins', type=float, default=5, help='below this protein count, a novel-genus call is flagged uncertain (default: 5)')
    parser.add_argument('--max-proteins', type=float, default=20, help='above this protein count, a novel-genus call is flagged uncertain (default: 20)')
    args = parser.parse_args(argv)

    df = pd.read_csv(args.assignments)

    query_df = filter_query_rows(df)
    print(f'{len(df)} total rows -> {len(query_df)} query rows (Genome contains "{QUERY_SEPARATOR}")')

    protein_counts = load_protein_counts(args.gene2genome)
    query_df = fill_protein_counts(query_df, protein_counts)
    assert_proteins_known(query_df)

    query_df = flag_uncertain_novel_genus(query_df, args.min_proteins, args.max_proteins)

    query_df = query_df.sort_values('Genome')
    query_df.to_csv(args.output, index=False)

    novel_taxa_df = extract_novel_taxa(query_df)
    novel_taxa_df.to_csv(args.novel_taxa_output, index=False)

    print(f'\nWrote {len(query_df)} postprocessed rows to {args.output}')
    print(f'Wrote {len(novel_taxa_df)} remaining novel-taxa rows to {args.novel_taxa_output}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
