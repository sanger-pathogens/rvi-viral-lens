#!/usr/bin/env nextflow

/*
========================================================================================
    viral metagraph-align sub-workflow
    ------------------------
    Aligns preprocessed paired reads (each mate independently, see metagraph.nf) against a
    pre-built metagraph de Bruijn graph + annotation index, counts reads per species from
    the alignment labels, and provisionally calls a species present in the sample once its
    read-hit count clears metagraph_align_min_hits. Hit count is the whole calling
    criterion: no reads are mapped here (see the note by the emit block, and
    subworkflows/classifying_index.nf's header for where breadth is measured instead).
    Reads are capped at metagraph_align_subsample_limit per mate before alignment (mirrors
    VIRAL_THEMISTO's subsampling): metagraph align's per-query output is one label line per
    matched k-mer window, so it scales with input depth far faster than a normal aligner's
    output would, and species calling only needs enough hits to clear metagraph_align_min_hits
    — feeding it full depth on a deep sample mostly just costs METAGRAPH/CALL_METAGRAPH_SPECIES
    time without changing which species get called.
    Per-sample results are published under <outdir>/<sample>/sequenceindex/metagraph_hits.
========================================================================================
*/

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/
include { METAGRAPH_ALIGN        } from '../rvi_toolbox/modules/metagraph_align.nf'
include { CALL_METAGRAPH_SPECIES } from '../rvi_toolbox/modules/metagraph_species_call.nf'
include { SUBSAMPLE_ITER         } from '../rvi_toolbox/subworkflows/subsample.nf'
include { INDEX_REFERENCE_LENGTHS } from '../rvi_toolbox/modules/reference_lengths.nf'

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

    // Cap METAGRAPH_ALIGN input at metagraph_align_subsample_limit reads per mate (mirrors
    // VIRAL_THEMISTO's SUBSAMPLE_ITER step).
    metagraph_align_subsample_limit_ch = Channel.value( params.metagraph_align_subsample_limit )

    reads_ch.map{ meta, read_1, read_2 ->
        def readCount = read_1.countFastq()
        [meta, read_1, read_2, readCount]
    }.set{ ready_for_subsampling }

    SUBSAMPLE_ITER(ready_for_subsampling, metagraph_align_subsample_limit_ch)
    capped_reads_ch = SUBSAMPLE_ITER.out.final_read_channel


    // -- Inputs for the two call gates (see rvi_toolbox/modules/themisto_species_call.nf's
    // header). Both are passed to the caller unconditionally -- a process input cannot be
    // conditionally absent -- so when a gate is off an empty placeholder from assets/ is
    // staged instead and the script reads a zero-byte file as "not supplied".
    taxon_table_ch = Channel.fromPath(
        params.run_taxon_filter ? params.taxon_filter_table
                                : "${projectDir}/assets/NO_TAXON_TABLE"
    ).first()

    // Reference lengths come from params.metagraph_map_reference_fasta -- the FASTA
    // metagraph's accession/taxid record ids name records in, and a DIFFERENT file from
    // the Themisto2 side's (see rvi_toolbox/modules/metagraph_species_call.nf).
    // Gated on the threshold rather than run unconditionally: INDEX_REFERENCE_LENGTHS
    // streams a multi-GB FASTA, which is wasted work if nothing will read the result.
    if (params.min_called_reference_length > 0) {
        INDEX_REFERENCE_LENGTHS(Channel.fromPath(params.metagraph_map_reference_fasta))
        reference_lengths_ch = INDEX_REFERENCE_LENGTHS.out.lengths.first()
    } else {
        reference_lengths_ch = Channel.fromPath("${projectDir}/assets/NO_REFERENCE_LENGTHS").first()
    }

    METAGRAPH_ALIGN(capped_reads_ch, graph_ch, annotation_ch, annotation_seqs_ch)

    // 'metagraph_hits': a subdir name distinct from VIRAL_METAGRAPH_QUERY.nf's, so the two
    // methods don't overwrite each other's output when both run for the same sample (see
    // CALL_METAGRAPH_SPECIES's output_subdir).
    CALL_METAGRAPH_SPECIES(
        METAGRAPH_ALIGN.out.alignments, names_dmp_ch, taxon_table_ch, reference_lengths_ch, 'metagraph_hits'
    )

    emit:
    species_hits    = CALL_METAGRAPH_SPECIES.out.species_hits
    // SEQIDX_<n> -> species for the species that cleared min-hits: the "ideal reference"
    // per call, which MAPPING uses to extract a consensus reference. Optional per sample
    // (unwritten when nothing cleared min-hits).
    index_label_map = CALL_METAGRAPH_SPECIES.out.index_label_map
}

/*
========================================================================================
    THE END
========================================================================================
*/
