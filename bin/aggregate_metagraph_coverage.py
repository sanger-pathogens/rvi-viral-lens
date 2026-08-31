#!/usr/bin/env python3
import argparse


def parse_coverage(path):
    """Parse samtools coverage output (header: #rname startpos endpos numreads covbases
    coverage meandepth meanbaseq meanmapq), keyed by rname (the reference accession)."""
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
    """accession -> species, from call_metagraph_species.py's out-index-label-map (one
    most-hit accession chosen per called species, so this is a strict 1:1 mapping)."""
    mapping = {}
    with open(path) as fh:
        for line in fh:
            accession, species = line.rstrip("\n").split("\t")
            mapping[accession] = species
    return mapping


def parse_species_hits(path):
    """species -> hit_count, from call_metagraph_species.py's out-species-hits TSV."""
    hits = {}
    with open(path) as fh:
        next(fh)  # header
        for line in fh:
            _sample_id, species, hit_count, _called = line.rstrip("\n").split("\t")
            hits[species] = int(hit_count)
    return hits


def main():
    parser = argparse.ArgumentParser(
        description="Build the per-species metagraph-map QC table: one row per called "
                    "species, from its single most-hit mapped reference accession."
    )
    parser.add_argument("--coverage", required=True, help="samtools coverage output")
    parser.add_argument("--query-lengths", required=True, help="rname -> total aligned query bases TSV")
    parser.add_argument("--index-label-map", required=True, help="accession -> species TSV")
    parser.add_argument("--species-hits", required=True, help="species -> hit_count TSV")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    coverage_rows = parse_coverage(args.coverage)
    query_lengths = parse_query_lengths(args.query_lengths)
    accession_to_species = parse_index_label_map(args.index_label_map)
    species_hits = parse_species_hits(args.species_hits)

    with open(args.out, "w") as out:
        out.write(
            "sample_id\tspecies\thit_count\treference_accession\treference_length\tquery_length\t"
            "covered_bases\tbreadth_pct\tmean_depth\tmeanbaseq\tmeanmapq\treads_mapped\n"
        )
        for accession, species in accession_to_species.items():
            row = coverage_rows.get(accession)
            if row is None:
                continue

            reference_length = int(row["endpos"]) - int(row["startpos"]) + 1
            out.write(
                f"{args.sample_id}\t{species}\t{species_hits.get(species, '')}\t{accession}\t{reference_length}\t"
                f"{query_lengths.get(accession, 0)}\t{row['covbases']}\t{row['coverage']}\t"
                f"{row['meandepth']}\t{row['meanbaseq']}\t{row['meanmapq']}\t{row['numreads']}\n"
            )


if __name__ == "__main__":
    main()
