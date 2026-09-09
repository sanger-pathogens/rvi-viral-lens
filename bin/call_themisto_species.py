#!/usr/bin/env python3
import argparse
import gzip


def open_maybe_gzip(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_species_labels(path):
    """Read species_labels.txt: one label per line, Nth line (1-based) = Nth sequence in
    the .thm2 index / reference FASTA — the same file mSWEEP uses as its -i ref_groups
    argument. Returned as a list indexed by 0-based Themisto reference index (dense
    array lookup rather than a dict: heavily-clustered viral indexes can pseudoalign one
    read against thousands of reference indices, so this lookup runs hundreds of millions
    of times per sample and a list avoids per-lookup hashing overhead)."""
    labels = []
    with open(path) as fh:
        for line in fh:
            labels.append(line.strip() or None)
    return labels


def parse_pseudoalignments(path, labels, species_hits, index_hits_by_species):
    """Parse one mate's --themisto1-output-format pseudoalignment file: the first token
    of each line is the (0-based) query/read number (unused — only relative order
    matters), every remaining token a 0-based index into the reference sequences the
    .thm2 index was built from (sorted ascending with --sort-output; a read with no
    pseudoalignment hits emits just its read number with no further tokens). A read is
    counted once per distinct species regardless of how many reference indices of that
    species it hits (mirrors call_metagraph_species.py's per-read species dedup). The
    inner loop is inlined (no per-token generator/function call) since it runs hundreds
    of millions of times on a heavily-clustered viral index."""
    num_labels = len(labels)
    with open_maybe_gzip(path) as fh:
        for line in fh:
            fields = line.split()
            if len(fields) < 2:
                continue
            species_in_read = set()
            for field in fields[1:]:
                ref_index = int(field)
                if ref_index >= num_labels:
                    continue
                label = labels[ref_index]
                if label is None:
                    continue
                species_in_read.add(label)
                counts = index_hits_by_species.setdefault(label, {})
                counts[ref_index] = counts.get(ref_index, 0) + 1
            for label in species_in_read:
                species_hits[label] = species_hits.get(label, 0) + 1


def main():
    parser = argparse.ArgumentParser(
        description="Count reads per species from Themisto2 pseudoalignment output (both "
                    "mates) and provisionally call a species present once its read-hit "
                    "count clears --min-hits; for each called species, pick its single "
                    "most-hit reference index to map against downstream."
    )
    parser.add_argument("--pseudoalignment-1", required=True, help="Themisto2 pseudoalignment output for read 1 (--themisto1-output-format, optionally gzipped)")
    parser.add_argument("--pseudoalignment-2", required=True, help="Themisto2 pseudoalignment output for read 2 (--themisto1-output-format, optionally gzipped)")
    parser.add_argument("--species-labels", required=True, help="species_labels.txt: one label per line, Nth line == Nth reference sequence in the .thm2 index")
    parser.add_argument("--min-hits", type=int, required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--out-species-hits", required=True, help="output TSV: sample_id, species, hit_count, provisional_call")
    parser.add_argument("--out-record-ids", required=True, help="output: one SEQIDX_<n> token per line, for seqkit grep -f (see reference_subset.nf)")
    parser.add_argument("--out-index-label-map", required=True, help="output TSV: SEQIDX_<n>, species")
    args = parser.parse_args()

    labels = parse_species_labels(args.species_labels)

    species_hits = {}
    index_hits_by_species = {}
    parse_pseudoalignments(args.pseudoalignment_1, labels, species_hits, index_hits_by_species)
    parse_pseudoalignments(args.pseudoalignment_2, labels, species_hits, index_hits_by_species)

    ranked_species = sorted(species_hits.items(), key=lambda kv: kv[1], reverse=True)

    with open(args.out_species_hits, "w") as out:
        out.write("sample_id\tspecies\thit_count\tprovisional_call\n")
        for species, hit_count in ranked_species:
            out.write(f"{args.sample_id}\t{species}\t{hit_count}\t{hit_count >= args.min_hits}\n")

    called_species = [species for species, hit_count in species_hits.items() if hit_count >= args.min_hits]
    if not called_species:
        # Nothing cleared min-hits: leave record-id/index-label-map files unwritten, so the
        # optional Nextflow outputs are empty and downstream mapping is skipped for this
        # sample rather than run on nothing.
        return

    with open(args.out_record_ids, "w") as out_ids, open(args.out_index_label_map, "w") as out_map:
        for species in called_species:
            most_common_index = max(index_hits_by_species[species], key=index_hits_by_species[species].get)
            record_id = f"SEQIDX_{most_common_index + 1}"  # 1-based, matches INDEX_REFERENCE_FASTA's tagging
            out_ids.write(f"{record_id}\n")
            out_map.write(f"{record_id}\t{species}\n")


if __name__ == "__main__":
    main()
