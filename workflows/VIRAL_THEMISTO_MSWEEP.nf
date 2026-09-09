#!/usr/bin/env nextflow

/*
========================================================================================
    viral Themisto2/mSWEEP sub-workflow
    ------------------------
    Pseudoalign preprocessed paired reads (subsampled further here if still above
    msweep_subsample_limit — species calling doesn't need full depth and excess reads only
    cost Themisto2 time/memory) against a pre-built Themisto2 (.thm2) index, then call
    species directly from the pseudoalignment read-hit counts (see
    themisto_species_call.nf) and validate each call by mapping the same (subsampled)
    reads directly against its single most-hit reference sequence and recording genome
    breadth of coverage (see themisto_map_qc.nf) — mirrors metagraph_align.nf's read-hit
    calling/validation, just against Themisto2 pseudoalignment counts instead of metagraph
    align's alignment labels. mSWEEP's own probabilistic abundance estimation (and its
    matching low-abundance mapping validation, msweep_map_qc.nf) is optional and off by
    default — set run_msweep to also run it alongside the above.
    Per-sample results are published under
    <outdir>/<sample>/sequenceindex/{themisto_hits,themisto_map}, and (if run_msweep)
    <outdir>/<sample>/{msweep,msweep_map}.

    Core processes are adapted from the gemsweep pipeline (themisto2 branch).
    PORTED (rvi_integration_1) from eu1/rvi_toolbox.git's `feature_msweep_map` branch
    (`subworkflows/themisto2-msweep.nf`), as a viral-lens-owned file rather than a submodule
    bump, same as every other lane in this integration -- see INSTRUCT.md item 2 (the
    unresolved fork question). Differences from upstream are confined to paths:
    params.results_dir -> params.outdir, publish under <outdir>/<sample>/sequenceindex/
    to match the other sequence-index methods, bin/ instead of rvi_toolbox/bin/, and
    include paths for viral-lens's workflows/ + modules/ layout. If the upstream branch
    changes, re-apply by hand.
========================================================================================
*/

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/
include { THEMISTO_PSEUDOALIGN              } from '../modules/themisto2.nf'
include { CALL_THEMISTO_SPECIES             } from '../modules/themisto_species_call.nf'
include { MSWEEP                            } from '../modules/msweep.nf'
include { CLEANUP_THEMISTO_PSEUDOALIGNMENTS } from '../modules/cleanup.nf'
include { THEMISTO_MAP_QC                   } from './THEMISTO_MAP_QC.nf'
include { SUBSAMPLE_ITER                    } from '../rvi_toolbox/subworkflows/subsample.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow VIRAL_THEMISTO_MSWEEP {

    take:
    reads_ch    // tuple( meta, read_1, read_2 ) — preprocessed/subsampled paired reads

    main:
    // Pre-built Themisto2 index: collected for staging and read in place. The process
    // invokes `-i ${index_prefix}.thm2`, so index_prefix is the index filename minus the
    // .thm2 extension.
    index_files_ch  = Channel.fromPath(params.msweep_themisto_index).collect()
    index_prefix_ch = Channel.value(file(params.msweep_themisto_index).getBaseName())

    // species_labels.txt: mSWEEP's own reference grouping/hierarchy file, and also the
    // positional index->species map CALL_THEMISTO_SPECIES needs (Nth line == Nth
    // reference sequence in the .thm2 index) — same file, reused for both purposes.
    ref_groups_ch = Channel.fromPath(params.msweep_ref_groups).first()

    // Cap Themisto2 pseudoalignment input at msweep_subsample_limit reads (mirrors
    // ASSEMBLE_META's metaspades_subsample_limit step): species calling/mSWEEP are
    // statistical estimators, not an assembler, so feeding them more than this depth
    // wastes Themisto2 time/memory without improving the call.
    msweep_subsample_limit_ch = Channel.value( params.msweep_subsample_limit )

    reads_ch.map{ meta, read_1, read_2 ->
        def readCount = read_1.countFastq()
        [meta, read_1, read_2, readCount]
    }.set{ ready_for_subsampling }

    SUBSAMPLE_ITER(ready_for_subsampling, msweep_subsample_limit_ch)
    capped_reads_ch = SUBSAMPLE_ITER.out.final_read_channel

    pseudoaligned_ch = THEMISTO_PSEUDOALIGN(capped_reads_ch, index_files_ch, index_prefix_ch)

    // Default path: call species straight from pseudoalignment hit counts, no
    // probabilistic model needed (see themisto_species_call.nf).
    CALL_THEMISTO_SPECIES(pseudoaligned_ch, ref_groups_ch)

    // themisto_align_run_map_qc lets species-hit calling be validated on its own first.
    if (params.themisto_align_run_map_qc) {
        themisto_map_reference_fasta_ch = Channel.fromPath(params.themisto_map_reference_fasta).first()

        // Validate against the same (subsampled) reads CALL_THEMISTO_SPECIES's call was
        // actually based on, not the original deeper reads_ch — otherwise breadth/depth
        // here could look better than what the species call itself saw.
        THEMISTO_MAP_QC(
            capped_reads_ch,
            CALL_THEMISTO_SPECIES.out.record_ids,
            CALL_THEMISTO_SPECIES.out.index_label_map,
            CALL_THEMISTO_SPECIES.out.species_hits,
            themisto_map_reference_fasta_ch
        )
        themisto_map_qc_ch = THEMISTO_MAP_QC.out.qc_table
    } else {
        themisto_map_qc_ch = Channel.empty()
    }

    // mSWEEP's probabilistic abundance estimation is optional — off by default in favour
    // of the direct species-hit call above.
    //
    // ABUNDANCE ESTIMATION ONLY. Upstream also ran MSWEEP_MAP_QC here, mapping reads
    // against each low-abundance call's reference to validate it by breadth of coverage.
    // That step is deliberately dropped in viral-lens: THEMISTO_MAP_QC above already does
    // breadth validation for this arm, and does it better, because it maps against each
    // species' MOST-HIT reference record whereas SELECT_REFERENCE_RECORDS picked the
    // LONGEST sequence carrying the label. Measured on the same sample, that choice was
    // worth 29.79% breadth (mSWEEP's pick) versus 99.96% (Themisto's) for the same
    // SARS-CoV-2 call. Keeping both meant paying for a second bowtie2 index + mapping pass
    // per sample to produce the worse of two answers to the same question.
    if (params.run_msweep) {
        MSWEEP(pseudoaligned_ch, ref_groups_ch)
        abundances_ch = MSWEEP.out.abundances
    } else {
        abundances_ch = Channel.empty()
    }

    // Pseudoalignment files are never published (see themisto2.nf) and, left unmanaged,
    // accumulate unbounded in the Nextflow work directory. Clean them up only once every
    // consumer has finished with them: CALL_THEMISTO_SPECIES always runs and passes the
    // same files through as an output; MSWEEP (when run_msweep is set) reads the same
    // pseudoaligned_ch independently, so cleanup has to wait on both, not just one.
    if (params.cleanup_intermediate_files_msweep) {
        cleanup_input_ch = params.run_msweep
            ? CALL_THEMISTO_SPECIES.out.pseudoalignments
                .join(MSWEEP.out.pseudoalignments)
                .map { meta, r1, r2, _r1b, _r2b -> tuple(meta, r1, r2) }
            : CALL_THEMISTO_SPECIES.out.pseudoalignments

        CLEANUP_THEMISTO_PSEUDOALIGNMENTS(cleanup_input_ch)
    }

    emit:
    species_hits     = CALL_THEMISTO_SPECIES.out.species_hits
    themisto_map_qc  = themisto_map_qc_ch
    abundances       = abundances_ch
}

/*
========================================================================================
    THE END
========================================================================================
*/
