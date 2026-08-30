#!/usr/bin/env python3
"""Concatenate per-sample viral fastas into a single pooled fasta, prefixing
each scaffold header with the sample id so the sample-of-origin is recoverable
after binning.

Output header convention:  >sampleId||originalHeader
This is the same convention vRhyme/ViWrap use for renaming geNomad outputs.
"""

import argparse
import sys


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--sample-ids", nargs="+", required=True,
                   help="Sample ids in the same order as --fastas")
    p.add_argument("--fastas", nargs="+", required=True,
                   help="Per-sample virus fasta paths")
    p.add_argument("--output", required=True,
                   help="Pooled fasta to write")
    args = p.parse_args(argv)

    if len(args.sample_ids) != len(args.fastas):
        sys.exit(
            f"--sample-ids ({len(args.sample_ids)}) and --fastas "
            f"({len(args.fastas)}) must have the same length"
        )

    with open(args.output, "w") as out:
        for sid, fna in zip(args.sample_ids, args.fastas):
            with open(fna) as fh:
                for line in fh:
                    if line.startswith(">"):
                        out.write(f">{sid}||{line[1:]}")
                    else:
                        out.write(line)

    return 0


if __name__ == "__main__":
    sys.exit(main())
