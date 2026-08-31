#!/usr/bin/env python3
import argparse


def parse_coverage(path):
    """Parse samtools coverage output (header: #rname startpos endpos numreads covbases
    coverage meandepth meanbaseq meanmapq), keyed by rname (our SEQIDX_<n> tag)."""
    rows = {}
    with open(path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                header = line.lstrip("#").split("\t")
                continue
            fields = line.split("\t")
            row = dict(zip(header, fields))
            rows[row["rname"]] = row
    return rows


def parse_query_lengths(path):
    """rname -> total aligned query (read) bases, from a samtools view | awk pass (samtools
    coverage itself has no such column)."""
    lengths = {}
    with open(path) as fh:
        for line in fh:
            rname, length = line.rstrip("\n").split("\t")
            lengths[rname] = int(length)
    return lengths


def parse_index_label_map(path):
    """SEQIDX_<n> -> label. One sequence per label now (SELECT_REFERENCE_RECORDS picks the
    single longest sequence per species), so this is a strict 1:1 mapping."""
    mapping = {}
    with open(path) as fh:
        for line in fh:
            record_id, label = line.rstrip("\n").split("\t")
            mapping[record_id] = label
    return mapping


def parse_groups(path):
    """label -> relative_abundance, from select_reference_records.py's out-groups TSV."""
    abundances = {}
    with open(path) as fh:
        next(fh)  # header
        for line in fh:
            label, relative_abundance = line.rstrip("\n").split("\t")
            abundances[label] = float(relative_abundance)
    return abundances


def main():
    parser = argparse.ArgumentParser(
        description="Build the per-species mSWEEP-map QC table: one row per species, "
                    "from its single mapped reference sequence."
    )
    parser.add_argument("--coverage", required=True, help="samtools coverage output")
    parser.add_argument("--query-lengths", required=True, help="rname -> total aligned query bases TSV")
    parser.add_argument("--index-label-map", required=True, help="SEQIDX_<n> -> label TSV")
    parser.add_argument("--groups", required=True, help="label -> relative_abundance TSV")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    coverage_rows = parse_coverage(args.coverage)
    query_lengths = parse_query_lengths(args.query_lengths)
    record_id_to_label = parse_index_label_map(args.index_label_map)
    abundances = parse_groups(args.groups)

    with open(args.out, "w") as out:
        out.write(
            "sample_id\tspecies_label\trelative_abundance\treference_length\tquery_length\t"
            "covered_bases\tbreadth_pct\tmean_depth\tmeanbaseq\tmeanmapq\treads_mapped\n"
        )
        for record_id, label in record_id_to_label.items():
            row = coverage_rows.get(record_id)
            if row is None:
                continue

            reference_length = int(row["endpos"]) - int(row["startpos"]) + 1
            out.write(
                f"{args.sample_id}\t{label}\t{abundances.get(label, '')}\t{reference_length}\t"
                f"{query_lengths.get(record_id, 0)}\t{row['covbases']}\t{row['coverage']}\t"
                f"{row['meandepth']}\t{row['meanbaseq']}\t{row['meanmapq']}\t{row['numreads']}\n"
            )


if __name__ == "__main__":
    main()
