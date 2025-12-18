// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.
process run_aligner {
    /*
    *              Map reads to reference

    The process is responsible for mapping sequencing reads
    to a reference genome using BWA (Burrows-Wheeler Aligner),
    followed by post-processing steps including conversion to
    BAM format, sorting, and indexing. At current version of
    this pipeline, this process generates high-quality, sorted
    BAM files that are essential for consensus sequence
    generation.

    * --------------------------------------------------------------
    * Input Tuple:
        - `meta`: Metadata associated with the sample. The process
            assumes there is:
            - `taxid` : taxonomic ID
            - `sample_id`: sample ID
        - `fastq`: Path to the FASTQ file(s) containing sequencing
            reads to be aligned.
        - `ref_fasta`: Path to the reference files required for
            alignment

    * Output Tuple:
        - Input tuple, plus:
          - Sorted BAM file (`*.sorted.bam`).
          - Index file for the sorted BAM (`*.sorted.bam.bai`).

    * ---------------------------------------------------------------
    *
    */
    tag "${meta.id}"
    label "aligner"
    label "mem_2"
    label "cpu_1"
    label "time_queue_from_small"

    input:
        tuple val(meta), path(fastq), path(ref_fa)

    output:
        tuple val(meta), path(fastq), path(ref_fa), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai")

    script:
        aligner = "${params.read_aligner}"
        aligner_params = "${params.read_aligner_params}"
        """
        set -e
        set -o pipefail

        if [ "$aligner" = "bwa" ]; then
            echo "[INFO] Running BWA on sample ${meta.id}"
            # bwa needs an index
            bwa index ${ref_fa}

            # Align and convert to BAM
            bwa mem ${aligner_params} ${ref_fa} ${fastq} \
                | samtools view -S -b - \
                | samtools sort -o ${meta.id}.sorted.bam

        elif [ "$aligner" = "minimap2" ]; then
            echo "[INFO] Running minimap2 on sample ${meta.id}"

            # Minimap2 does not need indexing
            minimap2 ${aligner_params} ${ref_fa} ${fastq} \
                | samtools view -S -b - \
                | samtools sort -o ${meta.id}.sorted.bam

        else
            echo "[ERROR] Unknown aligner: $aligner"
            exit 1
        fi

        # Index the sorted BAM
        samtools index ${meta.id}.sorted.bam
        """
}
/*

# Script Breakdown

1. **Error Handling**:

- `set -e`: Exits the script if any command returns a non-zero status,
    ensuring that errors are caught early.
- `set -o pipefail`: Ensures that errors in piped commands are
    correctly propagated.

2. **BWA Alignment**:

- Command: `bwa mem ${ref_fa} ${fastq} | samtools view -S -b - -o ${meta.id}.bam
  - Aligns the input sequencing reads (`${fastq}`) to the specified
    reference genome (`${ref_fa}`), producing a BAM file
    (`${meta.id}.sam`) containing the raw alignments.

3. **Sort BAM File**:

- Command: `samtools sort ${meta.id}.bam -o ${meta.id}.sorted.bam`
  - Sorts the BAM file by the leftmost coordinates of the alignments,
  a necessary step for many downstream analyses.

4. **Generate BAM Index**:

- Command: `samtools index ${meta.id}.sorted.bam`
  - Indexes the sorted BAM file, creating an index file
  (`*.sorted.bam.bai`) that allows for fast access to specific regions
  of the alignments.
*/