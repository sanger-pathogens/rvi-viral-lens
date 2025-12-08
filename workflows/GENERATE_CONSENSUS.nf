// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

include {run_aligner as initial_alignment; run_aligner as re_alignment; run_aligner as final_alignment} from '../modules/alignment.nf'
include {run_ivar as initial_consensus; run_ivar as polish_consensus} from '../modules/run_ivar.nf'
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
        sample_taxid_ch // tuple (meta, reads, ref_genome)

    main:

        // First round of alignment
        alignment_ch = initial_alignment( sample_taxid_ch )

        // set ivar input channel
        consensus_initial_in_ch = alignment_ch.map {
            meta, fastq, _ref_fa, bam, bam_idx ->
                [meta, fastq, bam, bam_idx]
        }
        
        consensus_initial_ch = initial_consensus(
            consensus_initial_in_ch, 
            params.mpileup_max_depth, 
            params.ivar_initial_freq_threshold, 
            params.ivar_initial_min_depth )

        if (params.do_consensus_polishing == 'true') {
            // Realign reads to consensus, and re-call consensus

            realignment_in_ch = consensus_initial_ch.map {
                meta, fastq, _bam, _bam_idx, cons_fa ->
                    [meta, fastq, cons_fa] 
            }

            realignment_ch = re_alignment( realignment_in_ch ) 

            consensus_final_in_ch = realignment_ch.map {
                meta, fastq, _ref_fa, bam, bam_idx ->
                    [meta, fastq, bam, bam_idx]
            }
        
            consensus_final_ch = polish_consensus(
                consensus_final_in_ch, 
                params.mpileup_max_depth, 
                params.ivar_polish_freq_threshold, 
                params.ivar_polish_min_depth )
        } else {
            consensus_final_ch = consensus_initial_ch
        }

        // Realign reads reads back to final consensus
        final_align_in_ch = consensus_final_ch.map {
            meta, fastq, _bam, _bam_idx, cons_fa ->
                [meta, fastq, cons_fa]    
        }

        final_alignment_ch = final_alignment( final_align_in_ch )

        qc_script_in_ch = final_alignment_ch.map {
            meta, _fastq, cons_fa, bam, bam_idx ->
                [meta, bam, bam_idx, cons_fa]
        } 

        qc_ch = run_qc_script(qc_script_in_ch)

        qc_augmented_ch = qc_ch.map {
            meta, bam, bam_idx, consensus, qc_json ->
                def json_map = new groovy.json.JsonSlurper().parse(new File(qc_json.toString()))
                def filter_map = [("longest_non_n_subsequence"): json_map["longest_non_n_subsequence"]]
                def new_meta = meta.plus(filter_map)
                [new_meta, bam, bam_idx, consensus, qc_json] 
        }
        
        filtered_consensus_ch = qc_augmented_ch.filter{  it -> (it[0].longest_non_n_subsequence > 0) }

    emit:
        filtered_consensus_ch
}