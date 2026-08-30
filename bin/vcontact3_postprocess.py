#!/usr/bin/env python3
"""Post-processing for a pipeline-level vContact3 run's final_assignments.csv.

Given vContact3's final_assignments.csv (plus the run's own
HMMprofile_vog_support.h5 and the reference-database genome_report.parquet
for the db version used), produce final_assignments_postprocessed.csv where:

  1. Only query genomes are retained (rows whose Genome column contains '#'
     — this pipeline's sample-namespacing separator; reference genomes from
     the database never contain it).
  2. A query row's genus_prediction is flagged as uncertain when it is a
     novel-genus call (genus_prediction starts with 'novel_genus') AND the
     genome has too few or too many proteins for that call to be trusted
     (Proteins < --min-proteins or > --max-proteins). Flagged rows get an
     'uncertain_novel_' prefix prepended to genus_prediction.
  3. Rows where realm_prediction is internally inconsistent with a lower rank
     (kingdom/phylum/class/order/family/genus prediction implies a different
     realm than realm_prediction says) are re-derived from the genome's own
     VOG hits and corrected where that own evidence disagrees — the same
     fix as fix_realm_consistency.py, run here on the query-only subset.
  4. The result is sorted by Genome (alphabetical) and written to --output.
  5. The remaining novel-taxa calls — genus_prediction still starts with
     'novel_genus' after step 2 (i.e. NOT demoted to 'uncertain_novel_') —
     are additionally written out on their own to --novel-taxa-output, for
     easy downstream review of the genomes confident enough to flag as
     candidate novel genera.

Usage:
    vcontact3_postprocess.py \\
        --assignments       vcontact3_out/exports/final_assignments.csv \\
        --vog-support       vcontact3_out/HMMprofile_vog_support.h5 \\
        --genome-report     /path/to/db/RefSeq.<version>.genome_report.parquet \\
        --output            vcontact3_out/exports/final_assignments_postprocessed.csv \\
        --novel-taxa-output vcontact3_out/exports/final_assignments_noveltaxa.csv \\
        [--report           vcontact3_out/exports/postprocess_report.txt]
        [--min-proteins 5] [--max-proteins 20]
"""

import argparse

import pandas as pd

RANKS = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus']

# Same per-genome LCA thresholds vConTACT3 uses internally (config.py: vog_tax_genome_thresholds)
LCA_THRESHOLDS = {
    'kingdom': {'min_hits': 3, 'min_frac': 0.90},
    'phylum': {'min_hits': 3, 'min_frac': 0.90},
    'class': {'min_hits': 3, 'min_frac': 0.90},
    'order': {'min_hits': 2, 'min_frac': 0.80},
    'family': {'min_hits': 2, 'min_frac': 0.90},
    'genus': {'min_hits': 2, 'min_frac': 0.90},
}


def filter_query_rows(df: pd.DataFrame) -> pd.DataFrame:
    """Keep only query genomes (Genome contains '#'); drop reference genomes."""
    return df[df['Genome'].str.contains('#', na=False)].copy()


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


def build_realm_lookup(genome_report: pd.DataFrame) -> dict:
    """rank -> {taxon_label: realm}, built from reference genomes only."""
    lookup = {}
    for rank in RANKS:
        if rank not in genome_report.columns:
            continue
        sub = genome_report.dropna(subset=[rank, 'realm'])
        if sub.empty:
            continue
        mapping = sub.groupby(rank)['realm'].agg(lambda s: s.mode().iat[0])
        lookup[rank] = mapping.to_dict()
    return lookup


def find_inconsistent_rows(df: pd.DataFrame, realm_lookup: dict) -> pd.DataFrame:
    """Flag rows where any rank's prediction implies a realm different from realm_prediction."""
    flags = pd.Series(False, index=df.index)
    reasons = pd.Series('', index=df.index)

    for rank in RANKS:
        pred_col = f'{rank}_prediction'
        if pred_col not in df.columns or rank not in realm_lookup:
            continue

        implied = df[pred_col].map(realm_lookup[rank])
        mismatch = implied.notna() & df['realm_prediction'].notna() & (implied != df['realm_prediction'])
        flags |= mismatch

        for idx in df.index[mismatch]:
            reasons.loc[idx] += (
                f'{rank}={df.at[idx, pred_col]!r} implies realm={implied.at[idx]!r} '
                f'but realm_prediction={df.at[idx, "realm_prediction"]!r}; '
            )

    out = df.loc[flags].copy()
    out['_conflict_reason'] = reasons.loc[flags]
    return out


def recompute_from_own_vog(genome_id: str, vog_support: pd.DataFrame) -> dict:
    """Re-derive each rank purely from this genome's own VOG hits (vConTACT3's own per-genome LCA rule)."""
    sub = vog_support[vog_support['genome_id'] == genome_id]
    corrected = {}
    for rank in RANKS:
        rank_df = sub[sub['rank'] == f'lca_{rank}']
        if rank_df.empty:
            continue

        total_hits = rank_df['vog_hits'].sum()
        top = rank_df.loc[rank_df['vog_hits'].idxmax()]
        top_frac = top['vog_hits'] / total_hits if total_hits else 0

        thr = LCA_THRESHOLDS[rank]
        if total_hits >= thr['min_hits'] and top_frac >= thr['min_frac']:
            corrected[rank] = (top['label'], int(total_hits), top_frac)

    return corrected


def fix_realm_consistency(df: pd.DataFrame, vog_support: pd.DataFrame, genome_report: pd.DataFrame):
    """Correct rows whose realm_prediction conflicts with a lower-rank prediction,
    using each genome's own VOG evidence. Returns (corrected_df, report_lines)."""
    realm_lookup = build_realm_lookup(genome_report)
    inconsistent = find_inconsistent_rows(df, realm_lookup)

    report_lines = [
        f'Checked {len(df)} genomes; found {len(inconsistent)} with a realm/lower-rank taxonomy conflict.'
    ]

    corrected_df = df.copy()
    n_corrected = 0
    n_no_evidence = 0
    n_agreed = 0

    for idx, row in inconsistent.iterrows():
        genome = row['Genome']

        if bool(row.get('Reference')) is True:
            continue  # never touch reference rows (query-only df, so shouldn't trigger anyway)

        corrections = recompute_from_own_vog(genome, vog_support)
        if not corrections:
            n_no_evidence += 1
            report_lines.append(
                f'{genome}: CONFLICT DETECTED but no confident own-VOG evidence to correct with '
                f'(left as-is). {row["_conflict_reason"]}'
            )
            continue

        changes = []
        for rank, (label, hits, frac) in corrections.items():
            pred_col = f'{rank}_prediction'
            ev_col = f'{rank}_evidence'
            old_val = corrected_df.at[idx, pred_col]
            if old_val != label:
                changes.append(f'{rank}: {old_val!r} -> {label!r} ({hits} hits, {frac:.0%} dominant)')
                corrected_df.at[idx, pred_col] = label
                if ev_col in corrected_df.columns:
                    corrected_df.at[idx, ev_col] = 'VOG:postprocess_corrected'

        if changes:
            n_corrected += 1
            report_lines.append(f'{genome}: ' + '; '.join(changes))
        else:
            n_agreed += 1
            report_lines.append(f'{genome}: CONFLICT DETECTED but own-VOG evidence agreed with existing values (left as-is).')

    report_lines.append(
        f'\n{n_corrected} genome(s) corrected, {n_no_evidence} left unresolved (no own evidence), '
        f'{n_agreed} conflicts where own evidence agreed with the existing value.'
    )
    return corrected_df, report_lines


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--assignments', required=True, help='final_assignments.csv from a vcontact3 run')
    parser.add_argument('--vog-support', required=True, help='HMMprofile_vog_support.h5 from the same run')
    parser.add_argument('--genome-report', required=True, help='RefSeq.*.genome_report.parquet for the db version used')
    parser.add_argument('--output', required=True, help='where to write final_assignments_postprocessed.csv')
    parser.add_argument('--novel-taxa-output', required=True, help='where to write final_assignments_noveltaxa.csv (the still-novel-genus subset)')
    parser.add_argument('--report', default=None, help='optional: write a plain-text summary of the realm-consistency check')
    parser.add_argument('--min-proteins', type=float, default=5, help='below this protein count, a novel-genus call is flagged uncertain (default: 5)')
    parser.add_argument('--max-proteins', type=float, default=20, help='above this protein count, a novel-genus call is flagged uncertain (default: 20)')
    args = parser.parse_args(argv)

    df = pd.read_csv(args.assignments)
    vog_support = pd.read_hdf(args.vog_support, key='taxonomy_support')
    genome_report = pd.read_parquet(args.genome_report)

    query_df = filter_query_rows(df)
    print(f'{len(df)} total rows -> {len(query_df)} query rows (Genome contains "#")')

    query_df = flag_uncertain_novel_genus(query_df, args.min_proteins, args.max_proteins)

    corrected_df, report_lines = fix_realm_consistency(query_df, vog_support, genome_report)

    corrected_df = corrected_df.sort_values('Genome')
    corrected_df.to_csv(args.output, index=False)

    novel_taxa_df = extract_novel_taxa(corrected_df)
    novel_taxa_df.to_csv(args.novel_taxa_output, index=False)

    full_report = '\n'.join(report_lines)
    print(full_report)
    if args.report:
        with open(args.report, 'w') as fh:
            fh.write(full_report + '\n')

    print(f'\nWrote {len(corrected_df)} postprocessed rows to {args.output}')
    print(f'Wrote {len(novel_taxa_df)} remaining novel-taxa rows to {args.novel_taxa_output}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
