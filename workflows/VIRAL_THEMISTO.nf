#!/usr/bin/env nextflow

/*
========================================================================================
    viral Themisto2 sub-workflow
    ------------------------
    Pseudoalign preprocessed paired reads (subsampled further here if still above
    msweep_subsample_limit — species calling doesn't need full depth and excess reads only
    cost Themisto2 time/memory) against a pre-built Themisto2 (.thm2) index, then call
    species directly from the pseudoalignment read-hit counts (see
    themisto_species_call.nf) — mirrors metagraph_align.nf's read-hit calling, just
    against Themisto2 pseudoalignment counts instead of metagraph align's alignment
    labels. Hit count is the whole calling criterion: no reads are mapped here (see the
    note by the emit block, and subworkflows/classifying_index.nf's header). mSWEEP's own
    probabilistic abundance estimation is optional and off by default — set run_msweep to
    also run it alongside the above.
    Per-sample results are published under
    <outdir>/<sample>/sequenceindex/themisto_hits, and (if run_msweep)
    <outdir>/<sample>/msweep.

    Core processes are adapted from the gemsweep pipeline (themisto2 branch).
    PORTED (rvi_integration_1) from eu1/rvi_toolbox.git's `feature_msweep_map` branch
    (`subworkflows/themisto2-msweep.nf`, renamed here since mSWEEP was split out), as a viral-lens-owned file rather than a submodule
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
include { CLEANUP_THEMISTO_PSEUDOALIGNMENTS } from '../modules/cleanup.nf'
include { SUBSAMPLE_ITER                    } from '../rvi_toolbox/subworkflows/subsample.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow VIRAL_THEMISTO {

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

    // Cap Themisto2 pseudoalignment input at msweep_subsample_limit reads PER MATE
    // (mirrors ASSEMBLE_META's metaspades_subsample_limit step): species calling is a
    // statistical estimate, not an assembly, so more depth costs Themisto2 time/memory
    // without changing which species get called. NB the limit is per mate, so the total
    // is twice this -- SUBSAMPLE_ITER's own id suffix reports the doubled figure.
    msweep_subsample_limit_ch = Channel.value( params.msweep_subsample_limit )

    reads_ch.map{ meta, read_1, read_2 ->
        def readCount = read_1.countFastq()
        [meta, read_1, read_2, readCount]
    }.set{ ready_for_subsampling }

    SUBSAMPLE_ITER(ready_for_subsampling, msweep_subsample_limit_ch)

    // SUBSAMPLE_ITER (rvi_toolbox) is LOSSY, and asymmetrically so — restore the meta
    // here or every downstream key silently stops matching.
    //
    // Its two branches emit different meta shapes (see
    // rvi_toolbox/subworkflows/subsample.nf): a sample already below the limit passes
    // through with its meta untouched, but a sample that actually gets subsampled has its
    // meta REBUILT FROM SCRATCH as `meta_new = [:]` carrying a single renamed key,
    // `id = "<original id>_subsampled<limit*2 abbreviated>-<iteration>"`. Everything else,
    // `sample_id` included, is dropped.
    //
    // That cost a real multi-sample run: 2M-read samples subsampled, their ids became
    // `47726_2_53_subsampled2M-1`, and the classifier's report backbone (keyed on the
    // pre-subsample `meta.id`) no longer matched the counts keyed on the post-subsample
    // one. Being `join(..., remainder: true)`, the mismatch did not drop the row — it
    // yielded `meta == null` and aborted the run with
    // "Cannot execute null+{themisto_n_species_considered=115, ...}". Small test samples
    // never tripped it because they sat below the limit and so kept their meta.
    //
    // combine(by: 0), not join: with subsample_iterations > 1 the left side has repeated
    // ids, which join would not tolerate.
    original_meta_ch = reads_ch.map { meta, _r1, _r2 -> [meta.id, meta] }

    capped_reads_ch = SUBSAMPLE_ITER.out.final_read_channel
        .map { meta, read_1, read_2 ->
            def matcher = (meta.id =~ /^(.+)_subsampled\d+[kMG]-\d+$/)
            [matcher ? matcher[0][1] : meta.id, read_1, read_2]
        }
        .combine(original_meta_ch, by: 0)
        .map { _id, read_1, read_2, meta -> [meta, read_1, read_2] }

    pseudoaligned_ch = THEMISTO_PSEUDOALIGN(capped_reads_ch, index_files_ch, index_prefix_ch)

    // Default path: call species straight from pseudoalignment hit counts, no
    // probabilistic model needed (see themisto_species_call.nf).
    CALL_THEMISTO_SPECIES(pseudoaligned_ch, ref_groups_ch)

    // NO map-QC here. Species are called on read-hit counts alone
    // (themisto_align_min_hits); the validation mapping that used to follow -- bowtie2 the
    // reads against each called species' reference, then samtools coverage for breadth --
    // was removed deliberately. Its breadth figure is now obtained downstream instead,
    // from the consensus alignment subworkflows/mapping.nf performs anyway, so a
    // sequence-index species is mapped once rather than twice. THEMISTO_MAP_QC.nf is kept
    // but unused; see its header.

    // mSWEEP is NOT run here any more. It is an abundance estimator, not a classifier,
    // so it now lives behind the abundance lane's own --run_msweep flag
    // (subworkflows/abundance.nf), consuming the `pseudoalignments` emitted below. This
    // workflow is purely Themisto2: pseudoalign, then call species from read-hit counts.
    // Pseudoalignment files are never published (see themisto2.nf) and, left unmanaged,
    // accumulate unbounded in the Nextflow work directory. Clean them up only once every
    // consumer has finished with them: CALL_THEMISTO_SPECIES always runs and passes the
    // same files through as an output; MSWEEP (when run_msweep is set) reads the same
    // pseudoaligned_ch independently, so cleanup has to wait on both, not just one.
    // Only safe to clean up here when nothing else will read them. With --run_msweep set,
    // MSWEEP consumes the same files from the abundance lane, so cleanup is deferred to
    // there (subworkflows/abundance.nf) — deleting them here would race MSWEEP.
    if (params.cleanup_intermediate_files_msweep && !params.run_msweep) {
        CLEANUP_THEMISTO_PSEUDOALIGNMENTS(CALL_THEMISTO_SPECIES.out.pseudoalignments)
    }

    emit:
    species_hits     = CALL_THEMISTO_SPECIES.out.species_hits
    // SEQIDX_<n> -> species for the species that cleared min-hits: the "ideal reference"
    // per call, which map-QC used to consume and MAPPING now uses to extract a consensus
    // reference. Optional per sample (unwritten when nothing cleared min-hits).
    index_label_map  = CALL_THEMISTO_SPECIES.out.index_label_map
    // Passed through by CALL_THEMISTO_SPECIES, so consuming this cannot start before that
    // sample's species call is done. The abundance lane's optional MSWEEP reads these.
    pseudoalignments = CALL_THEMISTO_SPECIES.out.pseudoalignments
    // The species_labels.txt MSWEEP needs as its -i ref_groups, surfaced so the abundance
    // lane does not have to re-derive it from params.
    ref_groups       = ref_groups_ch
}

/*
========================================================================================
    THE END
========================================================================================
*/
