#!/usr/bin/env python3
"""Run-level summary writer shared by the mapping/assembly/abundance report lanes.

Each of these lanes' per-sample step just JSON-dumps that sample's `meta` (already
enriched via `meta.plus([...])` upstream) to a `<sample_id>.properties.json` file --
unlike the taxid lane's report, there's no separate qc/nextclade JSON to merge in, so
no per-lane python script is needed for that step (see the calling `.nf` workflow).

This script does the run-level concatenation: one JSON list of every sample's record,
and a CSV with one row per sample. Column set is the union of keys seen across all
records (rather than a fixed mapping) since these lanes don't have an established
column contract yet -- a sample missing a given key just gets a blank cell for it.
"""

import argparse
import csv
import json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-c", "--concat_files", nargs="+", required=True,
        help="Per-sample *.properties.json files to concatenate",
    )
    parser.add_argument(
        "-o", "--out_prefix", required=True,
        help="Prefix for the two output files: <prefix>_run_summary.json, <prefix>_summary_report.csv",
    )
    args = parser.parse_args()

    records = []
    for f in args.concat_files:
        with open(f) as fh:
            records.append(json.load(fh))

    records.sort(key=lambda r: r.get("id", ""))

    with open(f"{args.out_prefix}_run_summary.json", "w") as out:
        json.dump(records, out, indent=4)

    fieldnames = []
    for r in records:
        for k in r.keys():
            if k not in fieldnames:
                fieldnames.append(k)

    with open(f"{args.out_prefix}_summary_report.csv", "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fieldnames, restval="")
        writer.writeheader()
        writer.writerows(records)


if __name__ == "__main__":
    main()
