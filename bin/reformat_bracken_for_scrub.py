#!/usr/bin/env python3
import argparse
import csv
import re
import sys

SUBSAMPLE_SUFFIX_RE = re.compile(r"_subsampled\d+[A-Za-z]*-\d+$")


def parse_bracken_summary(path, sample_suffix):
    """Parse bracken_summary_report.tsv: column 1 is a pipe-delimited taxonomy lineage
    (spanning every rank, since kreport2mpa.py was run with --intermediate-ranks), and the
    remaining columns are one per sample (raw Bracken-reestimated counts). Keep only rows
    whose deepest (last '|'-delimited) segment is species-level ('s__'), stripping that
    prefix to get the species name. sample_suffix is stripped from each sample column
    header, since GENERATE_ABUNDANCE_SUMMARY's own cleanup step doesn't match the real
    per-sample filenames and leaves it in place. Samples that were subsampled also carry
    a "_subsampled<depth>-<n>" suffix (e.g. "_subsampled20M-1") from the subsampling step,
    which is stripped too so subsampled and non-subsampled runs of the same sample share
    one ID.

    Returns (sample_ids, {species: {sample_id: count}}).
    """
    species_counts = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        sample_ids = []
        for col in header[1:]:
            if sample_suffix and col.endswith(sample_suffix):
                col = col[: -len(sample_suffix)]
            col = SUBSAMPLE_SUFFIX_RE.sub("", col)
            sample_ids.append(col)

        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            lineage, counts = fields[0], fields[1:]

            last_rank = lineage.split("|")[-1]
            if not last_rank.startswith("s__"):
                continue  # not a species-level row

            species = last_rank[len("s__"):]
            species_counts[species] = dict(zip(sample_ids, counts))

    return sample_ids, species_counts


def read_plate_map(path):
    """Read a SCRuB plate map CSV. Uses utf-8-sig since these files are often exported
    from Excel with a leading BOM.

    Returns {sample_id: is_control} where is_control is a bool parsed from the
    'is_control' column (TRUE/FALSE).
    """
    with open(path, encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh)
        id_field = reader.fieldnames[0]
        return {
            row[id_field]: row["is_control"].strip().upper() == "TRUE"
            for row in reader
            if row[id_field]
        }


def zero_species_in_controls(species_counts, control_sample_ids, species_to_zero):
    """Force the given species' counts to 0 in every control sample, since some
    species (e.g. known lab/reagent contaminants) should never register in a blank."""
    zeroed_count = 0
    not_found = []
    for species in species_to_zero:
        if species not in species_counts:
            not_found.append(species)
            continue
        for sample_id in control_sample_ids:
            if species_counts[species].get(sample_id, "0") != "0":
                zeroed_count += 1
            species_counts[species][sample_id] = "0"

    if not_found:
        print(
            f"WARNING: species requested for zeroing in controls were not found in the "
            f"bracken summary and were skipped: {', '.join(not_found)}",
            file=sys.stderr,
        )
    if zeroed_count:
        print(
            f"INFO: zeroed {zeroed_count} non-zero (species, control sample) count(s) "
            f"per --zero-species-in-controls.",
            file=sys.stderr,
        )


def main():
    parser = argparse.ArgumentParser(
        description="Filter bracken_summary_report.tsv to species-level rows and "
                    "transpose it into SCRuB's expected samples x species orientation."
    )
    parser.add_argument("--bracken-summary", required=True, help="bracken_summary_report.tsv")
    parser.add_argument("--sample-suffix", default="_report_bracken_species",
                         help="suffix to strip from sample column headers to recover clean sample IDs")
    parser.add_argument("--plate-map", required=True,
                         help="SCRuB plate map CSV; samples not listed here are dropped, with a warning")
    parser.add_argument("--zero-species-in-controls", default="",
                         help="comma-separated species names (as they appear after stripping the 's__' "
                              "prefix) whose counts should be forced to 0 in every control sample "
                              "(plate map's is_control=TRUE), e.g. known contaminant species that "
                              "should never appear in a blank")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    sample_ids, species_counts = parse_bracken_summary(args.bracken_summary, args.sample_suffix)
    species_names = list(species_counts.keys())

    plate_map = read_plate_map(args.plate_map)
    missing = [s for s in sample_ids if s not in plate_map]
    if missing:
        print(
            f"WARNING: {len(missing)} sample(s) in {args.bracken_summary} are not present in "
            f"{args.plate_map} and will be dropped from {args.out}: {', '.join(missing)}",
            file=sys.stderr,
        )
    sample_ids = [s for s in sample_ids if s in plate_map]
    # Emit rows in the plate map's own order, not sorted order. SCRuB requires the
    # abundance matrix and the metadata to have equivalent row names and compares
    # them positionally -- with the rows sorted independently it fails outright with
    # "The row names of the `data` and `metadata` inputs must be equivalent!" for any
    # plate map that is not already in sorted order. read_plate_map returns a dict, so
    # (Python 3.7+) it preserves the file's row order.
    plate_map_order = {sample_id: i for i, sample_id in enumerate(plate_map)}
    sample_ids = sorted(sample_ids, key=lambda s: plate_map_order[s])

    species_to_zero = [s for s in args.zero_species_in_controls.split(",") if s]
    if species_to_zero:
        control_sample_ids = [s for s in sample_ids if plate_map[s]]
        zero_species_in_controls(species_counts, control_sample_ids, species_to_zero)

    with open(args.out, "w") as out:
        out.write("\t" + "\t".join(species_names) + "\n")
        for sample_id in sample_ids:
            row = [species_counts[species].get(sample_id, "0") for species in species_names]
            out.write(sample_id + "\t" + "\t".join(row) + "\n")


if __name__ == "__main__":
    main()
