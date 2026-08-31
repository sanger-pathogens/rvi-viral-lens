#!/usr/bin/env python3
import argparse
import re
from collections import Counter, defaultdict

# metagraph's coordinate indexes append the matched coordinate window to the label,
# and a reference record whose own name already carries one ends up with two:
#   'KM368312.1 | Alphainfluenzavirus influenzae:1612-1762:2314-2464'
# Every window has to come off the species key, or each one becomes its own species.
COORD_SUFFIX_RE = re.compile(r"(?::\d+-\d+)+$")

# `metagraph query --query-mode labels` joins a read's matched labels with ':' rather
# than the ';' used elsewhere, e.g.
#   'OV040211.1 | Betacoronavirus pandemicum:OY127609.1 | Betacoronavirus pandemicum'
# Split on a ':' only when what follows looks like '<accession> | ', so the ':' that
# introduces a coordinate window (':16-166') is left alone.
LABEL_JOIN_RE = re.compile(r":(?=[A-Za-z0-9_.]+ \| )")


def parse_label(label):
    """Split one metagraph annotation label into (species_key, full_label, record_id).
    species_key groups hits into one species-hit bucket; full_label is kept as a fallback
    human-readable display name (overridden by a names.dmp scientific name where one is
    found, see load_scientific_names); record_id is the token used to select this label's
    reference sequence out of a FASTA (for downstream mapping validation, see
    metagraph_map_qc.nf) where one exists. Three label shapes are in use across metagraph
    indexes:
      - 'kraken:taxid|<taxid>|<strain/segment description>' (flu_rbfish index): grouped by
        taxid, so multiple accessions/segments sharing a taxid collapse into one
        species-hit bucket. record_id is the taxid itself: the flu_rbfish reference library
        (a Kraken library.fna) uses this exact 'kraken:taxid|<taxid>|...' header shape too,
        so grepping by taxid reliably finds the matching record regardless of how its
        trailing description differs from this label's.
      - '<accession> | <species_name>' (bigviralindex-rvdbc coordinate index): grouped by
        species_name, as before. In a coordinate index the label also carries the matched
        coordinate window (possibly more than one, see COORD_SUFFIX_RE); those are stripped
        so every window of a species lands in the same bucket rather than its own.
      - a bare '<accession> <description>' (no taxid, no ' | ' separator): grouped by
        accession, one bucket per reference sequence.
    Returns None if the label matches none of these shapes.
    """
    if label.startswith("kraken:taxid|"):
        parts = label.split("|", 2)
        if len(parts) < 3:
            return None
        taxid = parts[1]
        return taxid, label, taxid
    if " | " in label:
        accession, _, species_name = label.partition(" | ")
        return COORD_SUFFIX_RE.sub("", species_name), label, accession
    if " " in label:
        accession, _, _rest = label.partition(" ")
        return accession, label, accession
    return None


def iter_read_labels(line):
    """Extract every label token from one 'metagraph align' output line. Column 1 is the
    read id, column 2 the original query sequence, and every column after that belongs to
    one of a variable number of alternative-alignment groups: (strand, aligned reference
    sequence, score, aligned length, CIGAR, mismatch count, ';'-joined matched labels) for
    an aligned read, or a single '*' placeholder group (no labels field at all) for a read
    that didn't align. Rather than parse that variable-width group structure positionally,
    scan every field after the read id/sequence for anything shaped like a known label (see
    parse_label) — no other field in these groups (strand, score, CIGAR, ...) can match one
    of those shapes, so this is equivalent to and far simpler than tracking group
    boundaries."""
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 2:
        return
    for field in fields[2:]:
        for chunk in field.split(";"):
            for candidate in LABEL_JOIN_RE.split(chunk):
                parsed = parse_label(candidate)
                if parsed is not None:
                    yield parsed


def parse_alignments(path, species_hits, label_hits_by_species, record_id_by_label):
    """Parse one 'metagraph align' output file (no header, see iter_read_labels for its
    shape). A read is counted once per distinct species regardless of how many labels of
    that species it hits."""
    with open(path) as fh:
        for line in fh:
            species_in_read = set()
            for species_key, full_label, record_id in iter_read_labels(line):
                species_in_read.add(species_key)
                label_hits_by_species[species_key][full_label] += 1
                record_id_by_label[full_label] = record_id
            for species_key in species_in_read:
                species_hits[species_key] += 1


def load_scientific_names(names_dmp_path):
    """taxid -> scientific name, from an NCBI-taxonomy-format names.dmp ('field\\t|\\t'
    separated rows). Only 'scientific name' rows are kept, since a taxid can have several
    name rows (acronym, equivalent name, ...). A taxid absent here (or with no scientific
    name row) just isn't in the returned dict — callers fall back to the raw label."""
    names = {}
    with open(names_dmp_path) as fh:
        for line in fh:
            fields = [field.strip() for field in line.rstrip("\n").split("|")]
            if len(fields) < 4:
                continue
            taxid, name, _unique_name, name_class = fields[0], fields[1], fields[2], fields[3]
            if name_class == "scientific name":
                names[taxid] = name
    return names


def build_record_id_pattern(record_id):
    """A grep -n -r (match full header, regex) pattern that selects record_id's reference
    sequence out of a FASTA, tolerant of however its header's trailing text differs from
    what's embedded in the alignment label. record_id is a bare taxid for
    'kraken:taxid|...'-shaped labels (see parse_label) — anchored on the fixed
    'kraken:taxid|<taxid>|' prefix all such headers share, so e.g. taxid 13000336 can't
    accidentally match a header for taxid 130003360. Otherwise record_id is already a
    complete accession, anchored to the start of the header and required to be immediately
    followed by whitespace or end-of-line so it can't match a longer accession sharing it
    as a prefix."""
    if record_id.isdigit():
        return f"^kraken:taxid\\|{record_id}\\|"
    return f"^{re.escape(record_id)}(\\s|$)"


def main():
    parser = argparse.ArgumentParser(
        description="Count reads per species from metagraph align output (both mates) and "
                    "provisionally call a species present once its read-hit count clears "
                    "--min-hits; for each called species, pick its most-hit reference "
                    "record to map against downstream."
    )
    parser.add_argument("--align-1", required=True, help="metagraph align output for read 1")
    parser.add_argument("--align-2", required=True, help="metagraph align output for read 2")
    parser.add_argument("--min-hits", type=int, required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--names-dmp", help="NCBI-taxonomy-format names.dmp: relabels taxid-keyed species with their scientific name instead of the raw alignment label")
    parser.add_argument("--out-species-hits", required=True, help="output TSV: sample_id, species, hit_count, provisional_call")
    parser.add_argument("--out-record-ids", required=True, help="output: one seqkit grep -n -r pattern per called species (see build_record_id_pattern)")
    parser.add_argument("--out-index-label-map", required=True, help="output TSV: record_id (taxid or accession), species")
    args = parser.parse_args()

    species_hits = Counter()
    label_hits_by_species = defaultdict(Counter)
    record_id_by_label = {}
    parse_alignments(args.align_1, species_hits, label_hits_by_species, record_id_by_label)
    parse_alignments(args.align_2, species_hits, label_hits_by_species, record_id_by_label)

    scientific_names = load_scientific_names(args.names_dmp) if args.names_dmp else {}

    def display_name(species_key):
        most_common_label, _ = label_hits_by_species[species_key].most_common(1)[0]
        return scientific_names.get(species_key, most_common_label)

    with open(args.out_species_hits, "w") as out:
        out.write("sample_id\tspecies\thit_count\tprovisional_call\n")
        for species_key, hit_count in species_hits.most_common():
            out.write(f"{args.sample_id}\t{display_name(species_key)}\t{hit_count}\t{hit_count >= args.min_hits}\n")

    called_species = [species_key for species_key, hit_count in species_hits.items() if hit_count >= args.min_hits]
    if not called_species:
        # Nothing cleared min-hits: leave record-id/index-label-map files unwritten, so the
        # optional Nextflow outputs are empty and downstream mapping is skipped for this
        # sample rather than run on nothing.
        return

    # One reference record per called species -- but two species can legitimately pick the
    # same best-hit record (closely related organisms sharing an accession's best label),
    # and emitting it twice puts a duplicate record in the extracted subset FASTA. bowtie2
    # indexes that happily, then samtools sort dies on the resulting header:
    #   [E::sam_hrecs_update_hashes] Duplicate entry "OZ031634.1" in sam header
    # so de-duplicate here. The index-label map keeps the first species to claim a record,
    # matching EXTRACT_METAGRAPH_REFERENCE_SUBSET's own "first match wins" behaviour.
    seen_record_ids = set()
    with open(args.out_record_ids, "w") as out_ids, open(args.out_index_label_map, "w") as out_map:
        for species_key in called_species:
            most_common_label, _ = label_hits_by_species[species_key].most_common(1)[0]
            record_id = record_id_by_label[most_common_label]
            if record_id in seen_record_ids:
                continue
            seen_record_ids.add(record_id)
            out_ids.write(f"{build_record_id_pattern(record_id)}\n")
            out_map.write(f"{record_id}\t{display_name(species_key)}\n")


if __name__ == "__main__":
    main()
