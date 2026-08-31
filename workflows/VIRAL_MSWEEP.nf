#!/usr/bin/env nextflow

/*
========================================================================================
    viral mSWEEP sub-workflow
    ------------------------
    Pseudoalign preprocessed paired reads (subsampled further here if still above
    msweep_subsample_limit — mSWEEP abundance estimation doesn't need full depth and
    excess reads only cost Themisto2 time/memory) against a pre-built Themisto2 (.thm2)
    index, estimate reference-group abundances with mSWEEP, then validate any
    low-abundance calls by mapping the same (subsampled) reads directly against the
    candidate species' reference sequence and recording genome breadth of coverage (see
    msweep_map_qc.nf) — an independent sanity check on mSWEEP's probabilistic calls, not
    mGEMS-style read binning.
    Per-sample results are published under <outdir>/<sample>/sequenceindex/{msweep,msweep_map}.

    Core processes are adapted from the gemsweep pipeline (themisto2 branch).
========================================================================================
*/

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/
include { THEMISTO_PSEUDOALIGN              } from '../modules/themisto2.nf'
include { MSWEEP                            } from '../modules/msweep.nf'
include { CLEANUP_THEMISTO_PSEUDOALIGNMENTS } from '../modules/cleanup.nf'
include { MSWEEP_MAP_QC                     } from './MSWEEP_MAP_QC.nf'
include { SUBSAMPLE_ITER                    } from '../rvi_toolbox/subworkflows/subsample.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow VIRAL_MSWEEP {

    take:
    reads_ch    // tuple( meta, read_1, read_2 ) — preprocessed/subsampled paired reads

    main:
    // Pre-built Themisto2 index: collected for staging and read in place. The process
    // invokes `-i ${index_prefix}.thm2`, so index_prefix is the index filename minus the
    // .thm2 extension.
    index_files_ch  = Channel.fromPath(params.msweep_themisto_index).collect()
    index_prefix_ch = Channel.value(file(params.msweep_themisto_index).getBaseName())

    // mSWEEP reference grouping / hierarchy file (value channel, reused per sample).
    ref_groups_ch = Channel.fromPath(params.msweep_ref_groups).first()

    // Cap Themisto2 pseudoalignment input at msweep_subsample_limit reads (mirrors
    // ASSEMBLE_META's metaspades_subsample_limit step): mSWEEP is a statistical estimator,
    // not an assembler, so feeding it more than this depth wastes Themisto2 time/memory
    // without improving the abundance estimate.
    msweep_subsample_limit_ch = Channel.value( params.msweep_subsample_limit )

    reads_ch.map{ meta, read_1, read_2 ->
        def readCount = read_1.countFastq()
        [meta, read_1, read_2, readCount]
    }.set{ ready_for_subsampling }

    SUBSAMPLE_ITER(ready_for_subsampling, msweep_subsample_limit_ch)
    capped_reads_ch = SUBSAMPLE_ITER.out.final_read_channel

    pseudoaligned_ch = THEMISTO_PSEUDOALIGN(capped_reads_ch, index_files_ch, index_prefix_ch)
    MSWEEP(pseudoaligned_ch, ref_groups_ch)

    // Validate against the same (subsampled) reads mSWEEP's call was actually based on,
    // not the original deeper reads_ch — otherwise breadth/depth here could look better
    // than what mSWEEP itself saw.
    MSWEEP_MAP_QC(capped_reads_ch, MSWEEP.out.abundances, ref_groups_ch)

    // Pseudoalignment files are never published (see themisto2.nf) and, left unmanaged, accumulate
    // unbounded in the Nextflow work directory. Clean them up only once MSWEEP has consumed them —
    // MSWEEP.out.pseudoalignments passes the same files through as an output, so this cleanup step
    // cannot run until the MSWEEP task for that sample has completed.
    if (params.cleanup_intermediate_files_msweep) {
        CLEANUP_THEMISTO_PSEUDOALIGNMENTS(MSWEEP.out.pseudoalignments)
    }

    emit:
    abundances  = MSWEEP.out.abundances
    map_qc      = MSWEEP_MAP_QC.out.qc_table
}

/*
========================================================================================
    THE END
========================================================================================
*/
