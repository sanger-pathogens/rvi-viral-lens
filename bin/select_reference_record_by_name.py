#!/usr/bin/env python3
import argparse
import random


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
        description="Given one already-known species name, pick its longest sequence "
                    "(ties broken randomly) in the positionally-aligned reference FASTA "
                    "used to build the Themisto2 index (species_labels.txt line N == "
                    "reference FASTA record N) -- same selection rule as "
                    "select_reference_records.py, without its abundance-threshold "
                    "gating: the caller already knows exactly which species it wants, it "
                    "just needs that species' reference sequence resolved."
    )
    parser.add_argument("--species-name", required=True)
    parser.add_argument("--species-labels", required=True, help="species_labels.txt (positionally aligned with the reference FASTA)")
    parser.add_argument("--sequence-lengths", required=True, help="INDEX_REFERENCE_FASTA's sequence_lengths.tsv: line_no, length")
    parser.add_argument("--seed", type=int, required=True, help="random seed, for reproducible tie-breaking")
    parser.add_argument("--out-record-id", required=True, help="output: one SEQIDX_<n> token, for seqkit grep")
    args = parser.parse_args()

    labels_by_line = parse_species_labels(args.species_labels)
    target = args.species_name.strip().lower()
    matching_lines = [
        line_no for line_no, label in labels_by_line.items()
        if label.strip().lower() == target
    ]
    if not matching_lines:
        # No reference record for this species name in msweep_ref_groups -- leave no
        # output file, so the optional Nextflow output is empty and this (sample,
        # species) pair is naturally dropped downstream (same "optional" convention
        # reference_subset.nf's own processes use throughout).
        return

    lengths_by_line = parse_sequence_lengths(args.sequence_lengths)
    max_length = max(lengths_by_line[line_no] for line_no in matching_lines)
    longest = [line_no for line_no in matching_lines if lengths_by_line[line_no] == max_length]
    rng = random.Random(args.seed)
    chosen_line = rng.choice(longest)  # random only among ties

    with open(args.out_record_id, "w") as out_id:
        out_id.write(f"SEQIDX_{chosen_line}\n")


if __name__ == "__main__":
    main()
