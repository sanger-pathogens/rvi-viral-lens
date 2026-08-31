#!/usr/bin/env nextflow

/*
========================================================================================
    viral metagraph-align sub-workflow
    ------------------------
    Aligns preprocessed paired reads (each mate independently, see metagraph.nf) against a
    pre-built metagraph de Bruijn graph + annotation index, counts reads per species from
    the alignment labels, and provisionally calls a species present in the sample once its
    read-hit count clears metagraph_align_min_hits. When metagraph_align_run_map_qc is set,
    each called species is additionally validated by mapping the same reads directly
    against its single most-hit reference record and recording genome breadth of coverage
    (see metagraph_map_qc.nf) — an independent sanity check on metagraph align's read-hit
    calls.
    Reads are capped at metagraph_align_subsample_limit per mate before alignment (mirrors
    VIRAL_MSWEEP's subsampling): metagraph align's per-query output is one label line per
    matched k-mer window, so it scales with input depth far faster than a normal aligner's
    output would, and species calling only needs enough hits to clear metagraph_align_min_hits
    — feeding it full depth on a deep sample mostly just costs METAGRAPH/CALL_METAGRAPH_SPECIES
    time without changing which species get called.
    Per-sample results are published under <outdir>/<sample>/sequenceindex/{metagraph_hits,metagraph_map}.
========================================================================================
*/

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/
include { METAGRAPH              } from '../modules/metagraph.nf'
include { CALL_METAGRAPH_SPECIES } from '../modules/metagraph_species_call.nf'
include { METAGRAPH_MAP_QC       } from './METAGRAPH_MAP_QC.nf'
include { SUBSAMPLE_ITER         } from '../rvi_toolbox/subworkflows/subsample.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow VIRAL_METAGRAPH_ALIGN {

    take:
    reads_ch    // tuple( meta, read_1, read_2 ) — preprocessed paired reads

    main:
    // Pre-built metagraph graph/annotation and the mapping-validation reference FASTA:
    // collected once and reused as value channels across every sample.
    graph_ch            = Channel.fromPath(params.metagraph_align_graph).first()
    annotation_ch       = Channel.fromPath(params.metagraph_align_annotation).first()
    annotation_seqs_ch  = Channel.fromPath(params.metagraph_align_annotation_seqs).first()
    names_dmp_ch        = Channel.fromPath(params.metagraph_align_names_dmp).first()

    // Cap METAGRAPH input at metagraph_align_subsample_limit reads per mate (mirrors
    // VIRAL_MSWEEP's SUBSAMPLE_ITER step).
    metagraph_align_subsample_limit_ch = Channel.value( params.metagraph_align_subsample_limit )

    reads_ch.map{ meta, read_1, read_2 ->
        def readCount = read_1.countFastq()
        [meta, read_1, read_2, readCount]
    }.set{ ready_for_subsampling }

    SUBSAMPLE_ITER(ready_for_subsampling, metagraph_align_subsample_limit_ch)
    capped_reads_ch = SUBSAMPLE_ITER.out.final_read_channel

    METAGRAPH(capped_reads_ch, graph_ch, annotation_ch, annotation_seqs_ch)

    CALL_METAGRAPH_SPECIES(METAGRAPH.out.alignments, names_dmp_ch)

    // metagraph_align_run_map_qc lets species-hit calling be validated on its own first.
    if (params.metagraph_align_run_map_qc) {
        reference_fasta_ch = Channel.fromPath(params.metagraph_map_reference_fasta).first()

        // Validate against the same (subsampled) reads CALL_METAGRAPH_SPECIES's call was
        // actually based on, not the original deeper reads_ch — otherwise breadth/depth here
        // could look better than what the species call itself saw.
        METAGRAPH_MAP_QC(
            capped_reads_ch,
            CALL_METAGRAPH_SPECIES.out.record_ids,
            CALL_METAGRAPH_SPECIES.out.index_label_map,
            CALL_METAGRAPH_SPECIES.out.species_hits,
            reference_fasta_ch
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
