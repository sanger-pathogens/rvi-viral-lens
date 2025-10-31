// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

include {bwa_alignment_and_post_processing} from '../modules/bwa_alignment.nf'
include {run_ivar} from '../modules/run_ivar.nf'
include {run_qc_script} from '../modules/run_qc_script.nf'

workflow GENERATE_CONSENSUS {
    /*
    -----------------------------------------------------------------
    Obtain Consensus Sequences

    The GENERATE_CONSENSUS workflow performs read alignment and
    consensus sequence generation for sequencing data. It processes
    paired-end reads by aligning them to reference genomes using BWA,
    followed by consensus calling with iVar. This workflow is designed
    to take in sequencing data for different samples and taxonomic
    IDs, process them, and produce consensus sequences.

    -----------------------------------------------------------------
    # Inputs
    - **Sample Taxid Channel **: A channel containing tuples of
    metadata and paired-end FASTQ files. Metadata (`meta`) must
    include the following keys:
        - `id`: Unique identifier combining sample ID and taxonomic
        ID.
        - `taxid`: Taxonomic ID of the sample.
        - `sample_id`: Sample identifier.

    -----------------------------------------------------------------
    # Key Processes
    - **BWA Alignment**: Aligns sequencing reads to the provided
    reference genome.
    - **Consensus Calling with iVar**: Generates consensus sequences
    from the aligned reads.

    -----------------------------------------------------------------
    # Outputs
    - `run_ivar.out`: A channel containing tuples of metadata and the generated consensus FASTA file.

    */

    take:
        sample_taxid_ch // tuple (meta, reads, ref_files)

    main:
        // align reads to reference
        bwa_alignment_and_post_processing(sample_taxid_ch)
        bams_ch = bwa_alignment_and_post_processing.out

        // set ivar input channel
        bams_ch
            | map {meta, _fastq, ref_fa, _ref_indices, bam, bam_idx ->
                [meta, bam, bam_idx, ref_fa]
            }
            | set {ivar_in_ch}

        run_ivar(ivar_in_ch)

        run_ivar.out
            .set { ivar_out_ch }

        run_qc_script(ivar_out_ch)

        run_qc_script.out
            .map {meta, bam, bam_idx, consensus, variants, qc_json ->
                def json_map = new groovy.json.JsonSlurper().parse(new File(qc_json.toString()))
                def filter_map = [("longest_non_n_subsequence"): json_map["longest_non_n_subsequence"]]
                def new_meta = meta.plus(filter_map)
                [new_meta, bam, bam_idx, consensus, variants, qc_json] }
            .set {qc_out_ch}


        qc_out_ch
        .filter{  it -> (it[0].longest_non_n_subsequence > 0) }
        .set{filtered_consensus_ch}

    emit:
        filtered_consensus_ch
}