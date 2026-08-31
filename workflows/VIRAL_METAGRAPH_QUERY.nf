#!/usr/bin/env nextflow

/*
========================================================================================
    viral metagraph-query sub-workflow
    ------------------------
    Pseudoaligns preprocessed paired reads (each mate independently, see
    metagraph_query.nf) against a pre-built metagraph de Bruijn graph + annotation index
    using metagraph query --query-mode labels, counts reads per species from the matched
    labels, and provisionally calls a species present in the sample once its read-hit
    count clears metagraph_align_min_hits. Same species-calling and (optional) map-QC
    validation as VIRAL_METAGRAPH_ALIGN.nf, and reuses the same shared index/threshold
    params -- the two subworkflows are alternative methods against the same reference
    data, not independently configured pipelines. See metagraph_query.nf for why this
    exists as a separate, simpler module rather than the old filter+query pipeline
    metagraph_align.nf itself replaced.
    Reads are capped at metagraph_align_subsample_limit per mate before querying, same
    reasoning as VIRAL_METAGRAPH_ALIGN.nf.
    Per-sample results are published under <outdir>/<sample>/sequenceindex/{metagraph_query_hits,metagraph_query_map}.
========================================================================================
*/

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/
include { METAGRAPH_QUERY        } from '../modules/metagraph_query.nf'
include { CALL_METAGRAPH_SPECIES } from '../modules/metagraph_species_call.nf'
include { METAGRAPH_MAP_QC       } from './METAGRAPH_MAP_QC.nf'
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

    // 'metagraph_query_hits'/'metagraph_query_map_summary': distinct from
    // VIRAL_METAGRAPH_ALIGN.nf's names -- see that subworkflow's identical comment.
    CALL_METAGRAPH_SPECIES(METAGRAPH_QUERY.out.alignments, names_dmp_ch, 'metagraph_query_hits')

    if (params.metagraph_align_run_map_qc) {
        reference_fasta_ch = Channel.fromPath(params.metagraph_map_reference_fasta).first()

        METAGRAPH_MAP_QC(
            capped_reads_ch,
            CALL_METAGRAPH_SPECIES.out.record_ids,
            CALL_METAGRAPH_SPECIES.out.index_label_map,
            CALL_METAGRAPH_SPECIES.out.species_hits,
            reference_fasta_ch,
            'metagraph_query_map',
            'metagraph_query_map_summary'
        )
        map_qc_ch = METAGRAPH_MAP_QC.out.qc_table
    } else {
        map_qc_ch = Channel.empty()
    }

    emit:
    species_hits = CALL_METAGRAPH_SPECIES.out.species_hits
    map_qc       = map_qc_ch
}

/*
========================================================================================
    THE END
========================================================================================
*/
