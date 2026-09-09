#!/usr/bin/env nextflow

/*
========================================================================================
    viral metagraph-query sub-workflow
    ------------------------
    Pseudoaligns preprocessed paired reads (each mate independently, see
    metagraph_query.nf) against a pre-built metagraph de Bruijn graph + annotation index
    using metagraph query --query-mode labels, counts reads per species from the matched
    labels, and provisionally calls a species present in the sample once its read-hit
    count clears metagraph_align_min_hits. Hit count is the whole calling criterion: no
    reads are mapped here (see the note by the emit block). Same species-calling as
    VIRAL_METAGRAPH_ALIGN.nf, and reuses the same shared index/threshold params -- the two subworkflows are alternative methods against the same reference
    data, not independently configured pipelines. See metagraph_query.nf for why this
    exists as a separate, simpler module rather than the old filter+query pipeline
    metagraph_align.nf itself replaced.
    Reads are capped at metagraph_align_subsample_limit per mate before querying, same
    reasoning as VIRAL_METAGRAPH_ALIGN.nf.
    Per-sample results are published under <outdir>/<sample>/sequenceindex/metagraph_query_hits.
========================================================================================
*/

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/
include { METAGRAPH_QUERY        } from '../modules/metagraph_query.nf'
include { CALL_METAGRAPH_SPECIES } from '../modules/metagraph_species_call.nf'
include { SUBSAMPLE_ITER         } from '../rvi_toolbox/subworkflows/subsample.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow VIRAL_METAGRAPH_QUERY {

    take:
    reads_ch    // tuple( meta, read_1, read_2 ) — preprocessed paired reads

    main:
    // Same pre-built graph/annotation/names_dmp/reference FASTA as
    // VIRAL_METAGRAPH_ALIGN.nf -- one shared index, two query strategies against it.
    graph_ch            = Channel.fromPath(params.metagraph_align_graph).first()
    annotation_ch       = Channel.fromPath(params.metagraph_align_annotation).first()
    annotation_seqs_ch  = Channel.fromPath(params.metagraph_align_annotation_seqs).first()
    names_dmp_ch        = Channel.fromPath(params.metagraph_align_names_dmp).first()

    metagraph_align_subsample_limit_ch = Channel.value( params.metagraph_align_subsample_limit )

    reads_ch.map{ meta, read_1, read_2 ->
        def readCount = read_1.countFastq()
        [meta, read_1, read_2, readCount]
    }.set{ ready_for_subsampling }

    SUBSAMPLE_ITER(ready_for_subsampling, metagraph_align_subsample_limit_ch)
    capped_reads_ch = SUBSAMPLE_ITER.out.final_read_channel

    METAGRAPH_QUERY(capped_reads_ch, graph_ch, annotation_ch, annotation_seqs_ch)

    // 'metagraph_query_hits': a name distinct from VIRAL_METAGRAPH_ALIGN.nf's, so the two
    // methods' outputs cannot overwrite each other when both run for the same sample.
    CALL_METAGRAPH_SPECIES(METAGRAPH_QUERY.out.alignments, names_dmp_ch, 'metagraph_query_hits')

    // NO map-QC here. Species are called on read-hit counts alone
    // (metagraph_align_min_hits); the validation mapping that used to follow -- bowtie2
    // the reads against each called species' reference, then samtools coverage for breadth
    // -- was removed deliberately. Its breadth figure is now obtained downstream instead,
    // from the consensus alignment subworkflows/mapping.nf performs anyway, so a
    // sequence-index species is mapped once rather than twice. METAGRAPH_MAP_QC.nf is kept
    // but unused; see its header.

    emit:
    species_hits    = CALL_METAGRAPH_SPECIES.out.species_hits
    // SEQIDX_<n> -> species for the species that cleared min-hits: the "ideal reference"
    // per call, which map-QC used to consume and MAPPING now uses to extract a consensus
    // reference. Optional per sample (unwritten when nothing cleared min-hits).
    index_label_map = CALL_METAGRAPH_SPECIES.out.index_label_map
}

/*
========================================================================================
    THE END
========================================================================================
*/
