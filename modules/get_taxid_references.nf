// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.
process get_taxid_reference_files{
    /*
    *         Fetch Fasta Sequence Files for a Given Taxid

     The get_taxid_reference_files process is designed to extract
     reference sequences corresponding to a specified taxonomic ID
     (taxid) from a larger FASTA file. This process retrieves the
     relevant sequence, writes it to an output file, and indexes
     it using BWA (Burrows-Wheeler Aligner). This step is essential
     for downstream analysis that requires taxon-specific reference
     sequences. In the context of the pipeline, it is used to
     extract all the reference files from the Kraken database which
     were observed on that input batch.

    * --------------------------------------------------------------

    * Input:
       - taxid: The taxonomic ID for which reference sequences are
           being retrieved.
       - kraken_db_library_path: The path to the source FASTA file
           containing sequences for multiple taxa.

    * Output:
       - A FASTA file (`${taxid}.fa`) containing the sequence
         corresponding to the specified taxid.
       - BWA index files generated from the FASTA file `(${taxid}.fa.*)`

    * NOTE: Output is optional, as sequences may not be found for all
    *      specified taxids.

    * --------------------------------------------------------------
    * > TODO: we need to either remove the optional or, at least,
    *            raise a warning when that happens.
    *
    * --------------------------------------------------------------
    */

    tag "${taxid}"
    label "get_taxid_reference_files"
    label "cpu_1"
    label "mem_100M"

    input:
        val(taxid)
        val(kraken_db_library_path)
    output:
        tuple val(taxid), path("${taxid}.fa"), path("${taxid}.fa.*"), stdout, optional : true

    script:
"""
extract_sequence.py ${taxid} ${taxid}.fa ${kraken_db_library_path}
"""
}

/*
# Script Breakdown

Python Script: The script is written in Python and is responsible for
extracting reference sequences for a given taxid from a source FASTA
file and indexing them with BWA.

- Reading the Source File:
  - The script reads through the `source_fna_path`, which contains
    sequences for multiple taxa.
  - It searches for headers (`>`) that contain the specified taxid
    (`|${taxid}|`).
  - Once a matching sequence is found, it reads the sequence data
    until the next header or end of the file.

- Output Generation:
  - If a matching sequence is found, it is written to a FASTA file named
    ${taxid}.fa.
  - The script prints the number of sequences found for the given taxid
    for logging and verification purposes.

- BWA Indexing:
  - If a sequence is found and written to the output file, BWA is used to
    index the FASTA file (`bwa index ${taxid}.fa`), preparing it for
    alignment and further analysis.
*/
