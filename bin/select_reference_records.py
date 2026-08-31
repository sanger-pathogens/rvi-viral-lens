#!/usr/bin/env python3
import argparse
import random


def parse_abundances(path, min_abundance):
    """Parse mSWEEP's <sample>_mSWEEP_abundances.txt (label, relative_abundance columns;
    '#'-prefixed lines and any non-numeric second column are treated as header/comment
    lines and skipped). Returns {label: relative_abundance} for labels above min_abundance."""
    selected = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t") if "\t" in line else line.split()
            if len(fields) < 2:
                continue
            label = fields[0]
            try:
                relative_abundance = float(fields[1])
            except ValueError:
                continue  # header row using a non-numeric second column
            if relative_abundance > min_abundance:
                selected[label] = relative_abundance
    return selected


def parse_species_labels(path):
    """Read species_labels.txt: one label per line, Nth line (1-based) = Nth sequence in
    the .thm2 index / reference FASTA. A label may repeat across several lines."""
    labels_by_line = {}
    with open(path) as fh:
        for i, line in enumerate(fh, start=1):
            label = line.strip()
            if label:
                labels_by_line[i] = label
    return labels_by_line


def parse_sequence_lengths(path):
    """Read INDEX_REFERENCE_FASTA's sequence_lengths.tsv: line_no, length (both 1-based /
    positionally aligned with species_labels.txt and the reference FASTA)."""
    lengths_by_line = {}
    with open(path) as fh:
        for line in fh:
            line_no, length = line.rstrip("\n").split("\t")
            lengths_by_line[int(line_no)] = int(length)
    return lengths_by_line


def main():
    parser = argparse.ArgumentParser(
        description="Filter mSWEEP abundances by threshold and, for each surviving "
                    "reference-group label, pick its longest sequence to map against "
                    "(ties broken randomly) rather than every sequence clustered under "
                    "that label."
    )
    parser.add_argument("--abundances", required=True, help="<sample>_mSWEEP_abundances.txt")
    parser.add_argument("--species-labels", required=True, help="species_labels.txt (positionally aligned with the reference FASTA)")
    parser.add_argument("--sequence-lengths", required=True, help="INDEX_REFERENCE_FASTA's sequence_lengths.tsv: line_no, length")
    parser.add_argument("--min-abundance", type=float, required=True)
    parser.add_argument("--seed", type=int, required=True, help="random seed, for reproducible tie-breaking")
    parser.add_argument("--out-groups", required=True, help="output TSV: label, relative_abundance")
    parser.add_argument("--out-record-ids", required=True, help="output: one SEQIDX_<n> token per line, for seqkit grep")
    parser.add_argument("--out-index-label-map", required=True, help="output TSV: SEQIDX_<n>, label")
    args = parser.parse_args()

    selected = parse_abundances(args.abundances, args.min_abundance)
    if not selected:
        # Nothing above threshold for this sample: write no output files at all, so the
        # optional Nextflow outputs are empty and downstream mapping/coverage steps are
        # skipped for this sample rather than run on nothing.
        return

    labels_by_line = parse_species_labels(args.species_labels)
    lengths_by_line = parse_sequence_lengths(args.sequence_lengths)

    matching_lines = {}  # label -> [line numbers]
    for line_no, label in labels_by_line.items():
        if label in selected:
            matching_lines.setdefault(label, []).append(line_no)

    rng = random.Random(args.seed)
    chosen_line = {}
    for label, line_numbers in matching_lines.items():
        max_length = max(lengths_by_line[line_no] for line_no in line_numbers)
        longest = [line_no for line_no in line_numbers if lengths_by_line[line_no] == max_length]
        chosen_line[label] = rng.choice(longest)  # random only among ties

    with open(args.out_groups, "w") as out_groups:
        out_groups.write("label\trelative_abundance\n")
        for label, relative_abundance in selected.items():
            out_groups.write(f"{label}\t{relative_abundance}\n")

    with open(args.out_record_ids, "w") as out_ids, open(args.out_index_label_map, "w") as out_map:
        for label, line_no in chosen_line.items():
            record_id = f"SEQIDX_{line_no}"
            out_ids.write(f"{record_id}\n")
            out_map.write(f"{record_id}\t{label}\n")


if __name__ == "__main__":
    main()
