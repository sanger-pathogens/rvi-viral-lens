#!/usr/bin/env python3
"""Per-sample prep for the pipeline-level vContact3 step.

For a single sample, given:
  - vRhyme's `vRhyme_best_bins_fasta/` directory of `vRhyme_bin_<N>.faa` files
    (each one holding the proteins predicted on a bin's scaffolds)
  - vRhyme's best_bins membership.tsv (scaffold -> bin)
  - geNomad's virus_proteins.faa (per-scaffold predicted proteins)
  - geNomad's virus_summary.tsv (per-scaffold metadata incl. n_genes, length)
  - CheckV's quality_summary.tsv for geNomad's virus scaffolds (checkv_quality
    per scaffold)

write a sample-namespaced protein FASTA, gene2genome TSV, and lengths TSV
where:
  - each vRhyme bin becomes one "genome" (`<sample_id>||bin_<N>`), with its
    proteins taken from `vRhyme_bin_<N>.faa` — i.e. vRhyme's own gene calls,
    not geNomad's;
  - each unbinned scaffold (in virus_summary with n_genes > 0 but not in the
    membership) becomes a singleton "genome" (`<sample_id>||<scaffold>`),
    with its proteins taken from geNomad's virus_proteins.faa;
  - a BINNED scaffold whose CheckV checkv_quality is High-quality or
    Medium-quality is ALSO emitted as its own standalone singleton "genome"
    (`<sample_id>||<scaffold>||standalone`), on top of contributing to its
    bin — giving vContact3 both the bin-level and scaffold-level signal for
    scaffolds confident enough to stand on their own. Protein/genome IDs for
    this case carry an explicit `||standalone` marker so they never collide
    with that scaffold's proteins already emitted under its bin;
  - protein IDs are `<sample_id>||` prefixed to disambiguate across samples;
  - bin lengths are the sum of constituent scaffold lengths from
    virus_summary.tsv; singleton lengths (unbinned or standalone-promoted)
    are the scaffold's own length.

Usage:
    vcontact3_prep.py --sample-id ID --proteins F.faa --summary S.tsv \\
        --bin-membership M.tsv --bins-fasta BINS_DIR \\
        --checkv-quality Q.tsv \\
        --out-proteins OUT.faa --out-gene2genome OUT.tsv --out-lengths OUT.tsv

Pass an empty/header-only TSV as --bin-membership AND an empty directory as
--bins-fasta when vRhyme produced no bins; in that case all eligible
scaffolds get emitted as singleton genomes via geNomad's proteins. Pass an
empty/header-only TSV as --checkv-quality when CheckV produced nothing for
this sample; in that case no scaffold is promoted to standalone.
"""

import argparse
import csv
import re
import sys
from pathlib import Path

PROMOTE_QUALITIES = {"High-quality", "Medium-quality"}


def load_summary(summary_tsv: str):
    """Return ({scaffolds with n_genes > 0}, {scaffold -> length})."""
    eligible: set = set()
    lengths: dict = {}
    with open(summary_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            seq = row.get("seq_name")
            if not seq:
                continue
            try:
                n_genes = int(row.get("n_genes") or 0)
            except (TypeError, ValueError):
                n_genes = 0
            try:
                length = int(row.get("length") or 0)
            except (TypeError, ValueError):
                length = 0
            if n_genes > 0:
                eligible.add(seq)
            if length > 0:
                lengths[seq] = length
    return eligible, lengths


def _normalize_bin_id(raw) -> str:
    """Canonicalize a bin identifier to `bin_<N>` regardless of whether the
    source spelled it as '1', 'bin_1', 'Bin_1' etc."""
    s = str(raw).strip()
    low = s.lower()
    if low.startswith("bin_"):
        return f"bin_{s[4:]}"
    return f"bin_{s}"


def load_scaffold_to_bin(membership_tsv: str) -> dict:
    """Map scaffold -> canonical 'bin_<N>'. Empty dict if file is missing/empty."""
    mapping: dict = {}
    try:
        with open(membership_tsv) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                scaf = row.get("scaffold")
                bin_raw = row.get("bin")
                if scaf and bin_raw:
                    mapping[scaf] = _normalize_bin_id(bin_raw)
    except FileNotFoundError:
        pass
    return mapping


def load_scaffold_quality(checkv_quality_tsv: str) -> dict:
    """Map scaffold -> checkv_quality. Empty dict if file is missing/empty."""
    quality: dict = {}
    try:
        with open(checkv_quality_tsv) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                scaf = row.get("contig_id")
                q = row.get("checkv_quality")
                if scaf and q:
                    quality[scaf] = q
    except FileNotFoundError:
        pass
    return quality


_FAA_BIN_RE = re.compile(r"^vRhyme_bin_(\S+?)\.faa$")


def discover_bin_faa_files(bins_fasta_dir: str) -> dict:
    """Return {canonical bin_id ('bin_<N>') -> Path to its .faa}.

    Empty dict if the directory is missing or has no matching files.
    """
    p = Path(bins_fasta_dir)
    if not p.is_dir():
        return {}
    out: dict = {}
    for f in sorted(p.iterdir()):
        m = _FAA_BIN_RE.match(f.name)
        if not m:
            continue
        out[f"bin_{m.group(1)}"] = f
    return out


# Prodigal-style protein header: '<scaffold>_<N>'. Anchor on the trailing _<digits>.
PROTEIN_HEADER_RE = re.compile(r"^(.+)_(\d+)$")


def scaffold_of(protein_id: str) -> str:
    match = PROTEIN_HEADER_RE.match(protein_id)
    return match.group(1) if match else protein_id


def write_outputs(
    sample_id: str,
    proteins_faa: str,
    eligible: set,
    scaffold_to_bin: dict,
    scaffold_lengths: dict,
    scaffold_quality: dict,
    bin_faa_files: dict,
    out_proteins: str,
    out_gene2genome: str,
    out_lengths: str,
) -> None:
    # ---- compute genome lengths ---------------------------------------------
    # Binned genome: sum of constituent scaffold lengths.
    # Unbinned eligible scaffold: its own length.
    # Standalone-promoted (binned but high/medium quality) scaffold: its own
    # length, under a distinct genome_id so it doesn't merge with its bin.
    bin_to_scaffolds: dict = {}
    for scaffold, bin_id in scaffold_to_bin.items():
        bin_to_scaffolds.setdefault(bin_id, []).append(scaffold)

    genome_lengths: dict = {}
    for bin_id, scaffolds in bin_to_scaffolds.items():
        total = sum(scaffold_lengths.get(s, 0) for s in scaffolds)
        if total > 0:
            genome_lengths[f"{sample_id}||{bin_id}"] = total

    binned_scaffolds = set(scaffold_to_bin.keys())
    unbinned_scaffolds = eligible - binned_scaffolds
    promoted_scaffolds = {
        s for s in binned_scaffolds
        if s in eligible and scaffold_quality.get(s) in PROMOTE_QUALITIES
    }

    for scaffold in unbinned_scaffolds:
        length = scaffold_lengths.get(scaffold, 0)
        if length > 0:
            genome_lengths[f"{sample_id}||{scaffold}"] = length

    for scaffold in promoted_scaffolds:
        length = scaffold_lengths.get(scaffold, 0)
        if length > 0:
            genome_lengths[f"{sample_id}||{scaffold}||standalone"] = length

    with open(out_proteins, "w") as out_p, \
         open(out_gene2genome, "w") as out_g, \
         open(out_lengths, "w") as out_l:
        out_g.write("protein_id\tgenome_id\tkeywords\n")
        out_l.write("genome_id\tlength\n")
        for genome_id in sorted(genome_lengths):
            out_l.write(f"{genome_id}\t{genome_lengths[genome_id]}\n")

        # ---- emit bin proteins from vRhyme's per-bin .faa files -------------
        for bin_id, bin_faa in sorted(bin_faa_files.items()):
            genome_id = f"{sample_id}||{bin_id}"
            with open(bin_faa) as bf:
                for line in bf:
                    if line.startswith(">"):
                        header = line[1:].split()[0]
                        new_protein_id = f"{sample_id}||{header}"
                        out_p.write(f">{new_protein_id}\n")
                        out_g.write(f"{new_protein_id}\t{genome_id}\tNone\n")
                    else:
                        out_p.write(line)

        # ---- emit unbinned + standalone-promoted proteins from geNomad's
        # proteins.faa. A promoted scaffold's IDs carry a "||standalone" marker
        # so they never collide with the copies already emitted under its bin.
        keep = False
        standalone = False
        with open(proteins_faa) as fh_in:
            for line in fh_in:
                if line.startswith(">"):
                    header = line[1:].split()[0]
                    scaffold = scaffold_of(header)
                    standalone = scaffold in promoted_scaffolds
                    keep = scaffold in unbinned_scaffolds or standalone
                    if keep:
                        suffix = "||standalone" if standalone else ""
                        new_protein_id = f"{sample_id}||{header}{suffix}"
                        genome_id = f"{sample_id}||{scaffold}{suffix}"
                        out_p.write(f">{new_protein_id}\n")
                        out_g.write(f"{new_protein_id}\t{genome_id}\tNone\n")
                elif keep:
                    out_p.write(line)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--proteins", required=True, help="geNomad virus_proteins.faa (used for unbinned scaffolds)")
    parser.add_argument("--summary", required=True, help="geNomad virus_summary.tsv (n_genes, length)")
    parser.add_argument("--bin-membership", required=True,
                        help="vRhyme best_bins.N.membership.tsv (pass header-only TSV if no bins)")
    parser.add_argument("--bins-fasta", required=True,
                        help="vRhyme_best_bins_fasta/ directory holding per-bin .faa files "
                             "(pass an empty directory if vRhyme produced no bins)")
    parser.add_argument("--checkv-quality", required=True,
                        help="CheckV quality_summary.tsv for geNomad's virus scaffolds "
                             "(pass an empty/header-only TSV if CheckV produced nothing)")
    parser.add_argument("--out-proteins", required=True)
    parser.add_argument("--out-gene2genome", required=True)
    parser.add_argument("--out-lengths", required=True)
    args = parser.parse_args(argv)

    eligible, scaffold_lengths = load_summary(args.summary)
    scaffold_to_bin = load_scaffold_to_bin(args.bin_membership)
    scaffold_quality = load_scaffold_quality(args.checkv_quality)
    bin_faa_files = discover_bin_faa_files(args.bins_fasta)

    # Defensive check: every bin referenced in membership should have a .faa.
    # If not, log a warning to stderr but continue — those bins simply won't
    # contribute proteins (they'll still appear in lengths if their scaffolds
    # have summary entries).
    missing = sorted(set(scaffold_to_bin.values()) - set(bin_faa_files))
    if missing:
        sys.stderr.write(
            f"[vcontact3_prep] WARNING: membership references bins with no .faa file: {missing}\n"
        )

    write_outputs(
        sample_id=args.sample_id,
        proteins_faa=args.proteins,
        eligible=eligible,
        scaffold_to_bin=scaffold_to_bin,
        scaffold_lengths=scaffold_lengths,
        scaffold_quality=scaffold_quality,
        bin_faa_files=bin_faa_files,
        out_proteins=args.out_proteins,
        out_gene2genome=args.out_gene2genome,
        out_lengths=args.out_lengths,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
