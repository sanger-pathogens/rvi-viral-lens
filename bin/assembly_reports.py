#!/usr/bin/env python3
"""Report builder for the de novo assembly + viral binning lane.

Produces three run-level CSVs from the lane's own per-sample outputs:

  assembly_sample_summary_report.csv
      One row per sample: the lane's module counts plus the taxonomies geNomad
      and vContact3 assigned anywhere in that sample.

  assembly_scaffold_summary_report.csv
      One row per geNomad viral scaffold: CheckV's per-scaffold QC plus that
      scaffold's own geNomad taxonomy.

  vmag_scaffold_summary_report.csv
      One row per vRhyme bin (vMAG): CheckV's QC for the linked bin plus the
      taxonomy vContact3 assigned to that bin as a query genome.

Runs in two stages, mirroring the lane's own fan-in shape:

  per-sample   reduce one sample's geNomad/vRhyme/CheckV files to a small JSON
               of report rows. Everything here is sample-local.
  merge        concatenate those JSONs into the three CSVs, joining in the
               vContact3 taxonomy. vContact3 runs once for the whole batch, so
               its assignments are only available at this stage.

Usage:
    assembly_reports.py per-sample \\
        --sample-id          SAMPLE \\
        --genomad-summary    SAMPLE_virus_summary.tsv \\
        --vrhyme-membership  vRhyme_best_bins.0.membership.tsv \\
        --scaffold-quality   virus_scaffolds_quality_summary.tsv \\
        [--bin-quality       linked_bins_quality_summary.tsv] \\
        --output             SAMPLE.assembly_report_parts.json

    assembly_reports.py merge \\
        --parts                 *.assembly_report_parts.json \\
        [--vcontact3-assignments final_assignments_postprocessed.csv] \\
        --out-sample            assembly_sample_summary_report.csv \\
        --out-scaffold          assembly_scaffold_summary_report.csv \\
        --out-vmag              vmag_scaffold_summary_report.csv
"""

import argparse
import csv
import json
import os
import re
import sys

# vContact3's ranks, deepest first -- the order the "lowest assigned taxon" walk
# below follows.
VCONTACT3_RANKS = [
    'genus', 'subfamily', 'family', 'order', 'class', 'phylum', 'kingdom', 'realm',
]

# vContact3 fills unassigned ranks with generated placeholder labels rather than
# leaving them blank: 'novel_genus_0_of_...', 'unplaced_kingdom_within_...', and
# (from vcontact3_postprocess.py's protein-range demotion) 'uncertain_novel_...'.
# None of these name a real taxon, so the walk skips them.
PLACEHOLDER_PREFIXES = ('novel_', 'unplaced_', 'uncertain_novel_')

# geNomad writes a ';'-separated lineage whose root element is always 'Viruses'.
# A row classified no further than the root carries no usable taxonomy.
GENOMAD_ROOT = 'Viruses'

# CheckV columns carried through to both scaffold-level reports, in output order.
CHECKV_COLUMNS = [
    'contig_length', 'gene_count', 'checkv_quality', 'completeness', 'completeness_method',
]

SAMPLE_COLUMNS = [
    'Sample_ID',
    'genomad_n_viral_scaffolds',
    'taxonomy_geNomad',
    'vrhyme_n_bins',
    'vrhyme_n_binned_scaffolds',
    'checkv_n_high_quality',
    'checkv_n_medium_quality',
    'vcontact3_n_genomes',
    'vcontact3_n_clusters',
    'taxonomy_vcontact3',
]

SCAFFOLD_COLUMNS = ['Sample_ID', 'scaffold_id'] + CHECKV_COLUMNS + ['taxonomy_geNomad']

VMAG_COLUMNS = ['vMAG_ID', 'Sample_ID'] + CHECKV_COLUMNS + ['vcontact3_taxonomy']

_BIN_SUFFIX_RE = re.compile(r'^(?:vRhyme_)?bin_(\S+)$', re.IGNORECASE)


# --------------------------------------------------------------------------
# shared helpers
# --------------------------------------------------------------------------

def read_tsv(path):
    """Rows of a TSV as dicts. Missing or header-only file -> empty list."""
    if not path or not os.path.exists(path):
        return []
    with open(path, newline='') as fh:
        return list(csv.DictReader(fh, delimiter='\t'))


def normalize_bin_id(raw):
    """Canonicalize a bin identifier to 'bin_<N>'.

    Must match vcontact3_prep.py's _normalize_bin_id exactly: that is what mints
    the '<sample_id>||bin_<N>' genome IDs this report joins vContact3 against.
    CheckV names the linked bins 'vRhyme_bin_<N>', vRhyme's membership.tsv spells
    the same bin '<N>'.
    """
    s = str(raw).strip()
    if not s:
        return ''
    m = _BIN_SUFFIX_RE.match(s)
    if m:
        return f'bin_{m.group(1)}'
    return f'bin_{s}'


def genomad_lowest_taxon(taxonomy):
    """Lowest-ranked element of a geNomad ';'-separated lineage.

    '' when the row has no lineage, or was classified no deeper than the
    'Viruses' root.
    """
    if not taxonomy:
        return ''
    parts = [p.strip() for p in str(taxonomy).split(';') if p.strip()]
    if not parts:
        return ''
    lowest = parts[-1]
    if lowest == GENOMAD_ROOT:
        return ''
    return lowest


def is_placeholder(value):
    """True for vContact3's generated 'not actually assigned' labels."""
    return str(value).startswith(PLACEHOLDER_PREFIXES)


def vcontact3_lowest_taxon(row):
    """'<rank>:<taxon>' for the deepest rank vContact3 assigned a real name to.

    Walks genus -> realm and returns the first prediction that is populated and
    is not one of vContact3's novel_/unplaced_ placeholders, so the value names a
    taxon that exists rather than a generated label. '' when no rank qualifies.
    """
    for rank in VCONTACT3_RANKS:
        value = (row.get(f'{rank}_prediction') or '').strip()
        if value and not is_placeholder(value):
            return f'{rank}:{value}'
    return ''


def format_taxon_set(taxa):
    """Unique taxa as '(a, b, c)'. '' when there are none.

    Sorted so the cell is stable across runs rather than following row order.
    """
    unique = sorted({t for t in taxa if t})
    if not unique:
        return ''
    return '(' + ', '.join(unique) + ')'


# --------------------------------------------------------------------------
# per-sample stage
# --------------------------------------------------------------------------

def build_scaffold_rows(sample_id, genomad_rows, scaffold_quality_rows):
    """One row per CheckV-assessed viral scaffold, with its geNomad taxonomy."""
    taxonomy_by_scaffold = {
        row.get('seq_name', ''): genomad_lowest_taxon(row.get('taxonomy', ''))
        for row in genomad_rows
    }

    rows = []
    for row in scaffold_quality_rows:
        scaffold_id = row.get('contig_id', '')
        entry = {'Sample_ID': sample_id, 'scaffold_id': scaffold_id}
        for col in CHECKV_COLUMNS:
            entry[col] = row.get(col, '')
        entry['taxonomy_geNomad'] = taxonomy_by_scaffold.get(scaffold_id, '')
        rows.append(entry)
    return rows


def build_vmag_rows(sample_id, bin_quality_rows):
    """One row per vRhyme bin CheckV assessed.

    vMAG_ID namespaces the bin by sample ('<sample_id>_<contig_id>'), since
    CheckV's own contig_id ('vRhyme_bin_1') repeats across samples. bin_id is
    kept alongside for the vContact3 join and dropped before writing the CSV.
    """
    rows = []
    for row in bin_quality_rows:
        contig_id = row.get('contig_id', '')
        entry = {
            'vMAG_ID': f'{sample_id}_{contig_id}',
            'Sample_ID': sample_id,
            'bin_id': normalize_bin_id(contig_id),
        }
        for col in CHECKV_COLUMNS:
            entry[col] = row.get(col, '')
        rows.append(entry)
    return rows


def build_sample_row(sample_id, genomad_rows, membership_rows, scaffold_quality_rows):
    """The sample's counts and its set of geNomad taxonomies.

    The vContact3 columns are added at merge time, once its batch-level
    assignments are available.
    """
    bins = [normalize_bin_id(r.get('bin', '')) for r in membership_rows if r.get('bin')]
    qualities = [r.get('checkv_quality', '') for r in scaffold_quality_rows]

    return {
        'Sample_ID': sample_id,
        'genomad_n_viral_scaffolds': len(genomad_rows),
        'taxonomy_geNomad': format_taxon_set(
            genomad_lowest_taxon(r.get('taxonomy', '')) for r in genomad_rows
        ),
        'vrhyme_n_bins': len(set(bins)),
        'vrhyme_n_binned_scaffolds': len(bins),
        'checkv_n_high_quality': qualities.count('High-quality'),
        'checkv_n_medium_quality': qualities.count('Medium-quality'),
    }


def run_per_sample(args):
    genomad_rows = read_tsv(args.genomad_summary)
    membership_rows = read_tsv(args.vrhyme_membership)
    scaffold_quality_rows = read_tsv(args.scaffold_quality)
    # Only samples where vRhyme produced at least one bin have a linked-bins
    # CheckV run, so this input is genuinely optional.
    bin_quality_rows = read_tsv(args.bin_quality)

    parts = {
        'sample_id': args.sample_id,
        'sample': build_sample_row(
            args.sample_id, genomad_rows, membership_rows, scaffold_quality_rows
        ),
        'scaffolds': build_scaffold_rows(args.sample_id, genomad_rows, scaffold_quality_rows),
        'vmags': build_vmag_rows(args.sample_id, bin_quality_rows),
    }

    with open(args.output, 'w') as fh:
        json.dump(parts, fh, indent=2)

    print(
        f'{args.sample_id}: {len(parts["scaffolds"])} scaffold row(s), '
        f'{len(parts["vmags"])} vMAG row(s) -> {args.output}'
    )
    return 0


# --------------------------------------------------------------------------
# merge stage
# --------------------------------------------------------------------------

def load_vcontact3_assignments(path):
    """{genome_id: {'taxon': '<rank>:<taxon>', 'genus': raw genus_prediction}}.

    Keyed by vContact3's own Genome column ('<sample_id>||bin_<N>' or
    '<sample_id>||<scaffold>'), which is what the callers join on. Genomes with
    no real assignment are still present, with an empty 'taxon'.

    'genus' is kept raw, placeholders and all, only to preserve
    vcontact3_n_clusters' original meaning: the number of distinct
    genus_prediction values across the sample's query genomes.
    """
    if not path or not os.path.exists(path):
        return {}
    with open(path, newline='') as fh:
        return {
            row['Genome']: {
                'taxon': vcontact3_lowest_taxon(row),
                'genus': (row.get('genus_prediction') or '').strip(),
            }
            for row in csv.DictReader(fh)
            if row.get('Genome')
        }


def write_csv(path, columns, rows):
    with open(path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, extrasaction='ignore', restval='')
        writer.writeheader()
        writer.writerows(rows)
    print(f'Wrote {len(rows)} row(s) to {path}')


def run_merge(args):
    parts = []
    for path in args.parts:
        with open(path) as fh:
            parts.append(json.load(fh))
    parts.sort(key=lambda p: p.get('sample_id', ''))

    assignments_by_genome = load_vcontact3_assignments(args.vcontact3_assignments)

    sample_rows = []
    scaffold_rows = []
    vmag_rows = []

    for part in parts:
        sample_id = part.get('sample_id', '')
        prefix = f'{sample_id}||'

        # This sample's query genomes, picked back out of the batch-level
        # assignments by the '<sample_id>||' prefix vcontact3_prep.py mints.
        sample_assignments = {
            genome: assignment
            for genome, assignment in assignments_by_genome.items()
            if genome.startswith(prefix)
        }

        sample_row = dict(part.get('sample', {}))
        sample_row['vcontact3_n_genomes'] = len(sample_assignments)
        sample_row['vcontact3_n_clusters'] = len(
            {a['genus'] for a in sample_assignments.values()}
        )
        sample_row['taxonomy_vcontact3'] = format_taxon_set(
            a['taxon'] for a in sample_assignments.values()
        )
        sample_rows.append(sample_row)

        scaffold_rows.extend(part.get('scaffolds', []))

        for row in part.get('vmags', []):
            row = dict(row)
            assignment = assignments_by_genome.get(
                f'{sample_id}||{row.get("bin_id", "")}'
            )
            row['vcontact3_taxonomy'] = assignment['taxon'] if assignment else ''
            vmag_rows.append(row)

    scaffold_rows.sort(key=lambda r: (r.get('Sample_ID', ''), r.get('scaffold_id', '')))
    vmag_rows.sort(key=lambda r: (r.get('Sample_ID', ''), r.get('vMAG_ID', '')))

    write_csv(args.out_sample, SAMPLE_COLUMNS, sample_rows)
    write_csv(args.out_scaffold, SCAFFOLD_COLUMNS, scaffold_rows)
    write_csv(args.out_vmag, VMAG_COLUMNS, vmag_rows)

    if not assignments_by_genome:
        print(
            'NOTE: no vContact3 assignments were supplied or the file was empty -- '
            'the vContact3 taxonomy columns are blank for every row.',
            file=sys.stderr,
        )
    return 0


# --------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = parser.add_subparsers(dest='command', required=True)

    ps = sub.add_parser('per-sample', help="reduce one sample's lane outputs to report rows")
    ps.add_argument('--sample-id', required=True)
    ps.add_argument('--genomad-summary', required=True, help="geNomad's virus_summary.tsv")
    ps.add_argument('--vrhyme-membership', required=True,
                    help="vRhyme's best_bins membership.tsv (scaffold -> bin)")
    ps.add_argument('--scaffold-quality', required=True,
                    help="CheckV's per-scaffold virus_scaffolds_quality_summary.tsv")
    ps.add_argument('--bin-quality', default=None,
                    help="CheckV's linked_bins_quality_summary.tsv; absent for samples with no bins")
    ps.add_argument('--output', required=True)
    ps.set_defaults(func=run_per_sample)

    mg = sub.add_parser('merge', help='concatenate per-sample parts into the three report CSVs')
    mg.add_argument('--parts', nargs='+', required=True, help='*.assembly_report_parts.json')
    mg.add_argument('--vcontact3-assignments', default=None,
                    help='final_assignments_postprocessed.csv for the batch (optional)')
    mg.add_argument('--out-sample', required=True)
    mg.add_argument('--out-scaffold', required=True)
    mg.add_argument('--out-vmag', required=True)
    mg.set_defaults(func=run_merge)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == '__main__':
    raise SystemExit(main())
