// -- abundance estimation (rvi_integration_1) --------------------------------
// Extracted unchanged from main.nf's inline body. KRAKEN2BRACKEN is a
// viral-lens-owned fork of rvi_toolbox's own subworkflow (same modules, added
// an emit: block -- see ../workflows/KRAKEN2BRACKEN.nf's header for why).
// ABUNDANCE_ESTIMATION is unmodified, included directly from the shared
// submodule. SCRUB_DECONTAM is ported from eu1/rvi_toolbox.git (that fork's
// only copy, same situation as Metagraph -- see INSTRUCT.md's "rvi_toolbox
// fork problem").
include {KRAKEN2BRACKEN} from '../workflows/KRAKEN2BRACKEN.nf'
// mSWEEP moved here from the sequence-index lane: it estimates abundances, it does not
// call species. It is the one process in this lane that is not self-sufficient -- it
// consumes Themisto2's pseudoalignments, so it can only run when CLASSIFYING_INDEX has
// actually produced them (see the run_msweep block below).
include {MSWEEP} from '../modules/msweep.nf'
include {CLEANUP_THEMISTO_PSEUDOALIGNMENTS} from '../modules/cleanup.nf'
include {ABUNDANCE_ESTIMATION} from '../rvi_toolbox/subworkflows/abundance_estimation.nf'
include {SCRUB_DECONTAM} from '../workflows/SCRUB_DECONTAM.nf'
include {GENERATE_ABUNDANCE_REPORT} from '../workflows/GENERATE_ABUNDANCE_REPORT.nf'
include {publish_lane_json as publish_abundance_lane_json} from '../modules/publish_lane_report.nf'
include {publish_run_files as publish_abundance_run_files} from '../modules/publish_lite.nf'

workflow ABUNDANCE {
    /*
    Kraken2+Bracken and ABUNDANCE_ESTIMATION run in parallel off the same
    preprocessed reads (not downstream of each other, and not downstream of
    the assembly/mapping lanes either). SCRuB is a final, whole-run step on
    Kraken2+Bracken's output only (see ../workflows/SCRUB_DECONTAM.nf) --
    ABUNDANCE_ESTIMATION doesn't go through it. Same join-backbone pattern as
    the sequence-index lane.
    */

    take:
        preprocessed_3tuple_ch     // tuple (meta, read1, read2)
        themisto_pseudoaln_ch      // CLASSIFYING_INDEX.out.themisto_pseudoalignments -- empty unless run_themisto
        themisto_ref_groups_ch     // CLASSIFYING_INDEX.out.themisto_ref_groups -- empty unless run_themisto

    main:
        abundance_sample_ch = preprocessed_3tuple_ch
            .map { meta, _r1, _r2 -> [meta.id, meta] }

        if (params.run_kraken2bracken) {
            KRAKEN2BRACKEN(preprocessed_3tuple_ch)

            bracken_counts_ch = KRAKEN2BRACKEN.out.mpa_abundance_report
                .map { meta, mpa -> [meta.id, count_bracken_species(mpa)] }

            if (params.run_scrub) {
                SCRUB_DECONTAM(KRAKEN2BRACKEN.out.abundance_summary)

                // SCRuB runs once per run against the whole-run abundance_summary, not
                // per sample -- its output (an RDS result object) has no natural
                // per-sample split, so every sample in this run just records that it
                // went through the step, not a per-sample decontamination magnitude.
                scrub_ran_ch = abundance_sample_ch
                    .map { id, _meta -> [id, [scrub_ran: true]] }
            } else {
                scrub_ran_ch = Channel.empty()
            }
        } else {
            bracken_counts_ch = Channel.empty()
            scrub_ran_ch = Channel.empty()
        }

        if (params.run_abundance_estimation) {
            ABUNDANCE_ESTIMATION(preprocessed_3tuple_ch)

            // ABUNDANCE_ESTIMATION (rvi_toolbox, unmodified -- included directly, not
            // forked) has no emit: block, so nothing per-sample is reachable from it
            // the way KRAKEN2BRACKEN's fork now is. Left as a pass-through call: its own
            // outputs still publish normally under outdir, but the report only records
            // that it ran, not per-sample metrics. Deepen this into a real wrapper (same
            // pattern as ../workflows/KRAKEN2BRACKEN.nf) only once it's clear the flag
            // alone isn't enough for the report -- same "don't build ahead of need"
            // reasoning as the assembly lane's deferred per-scaffold granularity.
            abund_est_ran_ch = abundance_sample_ch
                .map { id, _meta -> [id, [abundance_estimation_ran: true]] }
        } else {
            abund_est_ran_ch = Channel.empty()
        }

        // -- mSWEEP probabilistic abundance estimation (opt-in, --run_msweep) -----------
        // Off by default. Unlike everything else in this lane, mSWEEP does not start from
        // the reads: it reads Themisto2's pseudoalignments, so it is only meaningful when
        // the sequence-index lane actually ran Themisto2. Checked up front on the params
        // rather than left to fail late: without run_themisto the handover channels are
        // Channel.empty(), so MSWEEP would simply never be scheduled and --run_msweep
        // would look like it had been honoured while producing nothing at all.
        if (params.run_msweep) {
            if (!params.do_sequence_index || !params.run_themisto) {
                error("--run_msweep needs Themisto2's pseudoalignments, which are only " +
                      "produced by the sequence-index lane. Enable both " +
                      "--do_sequence_index true and --run_themisto true (the default), or " +
                      "drop --run_msweep. " +
                      "(do_sequence_index=${params.do_sequence_index}, " +
                      "run_themisto=${params.run_themisto})")
            }

            MSWEEP(themisto_pseudoaln_ch, themisto_ref_groups_ch.first())

            msweep_counts_ch = MSWEEP.out.abundances
                .map { meta, abundances, _probs -> [meta.id, count_msweep_abundances(abundances)] }

            // Themisto2's pseudoalignments are unpublished intermediates that accumulate
            // in the work directory. VIRAL_THEMISTO deliberately does NOT clean them up
            // when run_msweep is set (that would race this MSWEEP), so the cleanup lands
            // here instead, after MSWEEP has passed the same files through as an output.
            if (params.cleanup_intermediate_files_msweep) {
                CLEANUP_THEMISTO_PSEUDOALIGNMENTS(MSWEEP.out.pseudoalignments)
            }
        } else {
            msweep_counts_ch = Channel.empty()
        }

        abundance_sample_ch
            .join(bracken_counts_ch, remainder: true)
            .join(scrub_ran_ch, remainder: true)
            .join(abund_est_ran_ch, remainder: true)
            .join(msweep_counts_ch, remainder: true)
            .map { id, meta, b_counts, s_ran, ae_ran, m_counts ->
                def new_meta = meta + (b_counts ?: EMPTY_BRACKEN_COUNTS) +
                    (s_ran ?: [scrub_ran: false]) + (ae_ran ?: [abundance_estimation_ran: false]) +
                    (m_counts ?: EMPTY_MSWEEP_COUNTS)
                [id, new_meta]
            }
            .set { abundance_report_prep_ch }

        GENERATE_ABUNDANCE_REPORT(abundance_report_prep_ch)

        // PUBLISH (abundance lane)
        publish_abundance_lane_json(GENERATE_ABUNDANCE_REPORT.out.publish_seq_level_ch)
        publish_abundance_run_files(GENERATE_ABUNDANCE_REPORT.out.publish_run_level_summaries_ch)
}

// --- rvi_integration_1: sample-level count helper for the abundance report ---

EMPTY_BRACKEN_COUNTS = [bracken_n_species_called: 0]

// Moved here with mSWEEP from subworkflows/classifying_index.nf.
EMPTY_MSWEEP_COUNTS = [msweep_n_groups: 0, msweep_top_group: '', msweep_top_abundance: 0.0]

def count_msweep_abundances(txt) {
    // <sample>_mSWEEP_abundances.txt: "<label>\t<relative_abundance>", '#'-prefixed and
    // non-numeric-second-column lines skipped -- same rule bin/select_reference_records.py
    // applies in parse_abundances(), so these counts describe the same set of groups the
    // downstream map-QC step actually considered.
    if (txt == null || !txt.exists()) return EMPTY_MSWEEP_COUNTS
    def rows = []
    txt.readLines().each { line ->
        def trimmed = line.trim()
        if (!trimmed || trimmed.startsWith('#')) return
        def cols = trimmed.split('\t')
        if (cols.size() < 2) return
        try {
            rows << [cols[0], cols[1] as Double]
        } catch (NumberFormatException ignored) {
            // header or malformed row -- skipped, as parse_abundances does
        }
    }
    if (!rows) return EMPTY_MSWEEP_COUNTS
    def above = rows.findAll { row -> row[1] >= params.msweep_map_min_abundance }
    def top = rows.max { row -> row[1] }
    return [
        msweep_n_groups:      above.size(),
        msweep_top_group:     top[0],
        msweep_top_abundance: top[1]
    ]
}


def count_bracken_species(mpa) {
    // <sample>_report_bracken_species.mpa.txt (rvi_toolbox/modules/krakentools.nf's
    // KREPORT2MPA): tab-separated "<pipe-delimited lineage>\t<count>", one row per
    // taxonomic rank (kreport2mpa.py run with --intermediate-ranks, so every rank
    // appears, not just leaves). Counts species-level rows (last lineage segment starts
    // 's__', same rule bin/reformat_bracken_for_scrub.py's parse_bracken_summary applies)
    // with a non-zero read count.
    if (mpa == null || !mpa.exists()) return EMPTY_BRACKEN_COUNTS
    def n_called = 0
    mpa.readLines().each { String line ->
        def cols = line.split('\t')
        if (cols.size() < 2) return
        def last_rank = cols[0].split('\\|')[-1]
        if (!last_rank.startsWith('s__')) return
        try {
            if ((cols[1] as Double) > 0) n_called++
        } catch (NumberFormatException ignored) {
            // header or malformed row -- skipped
        }
    }
    return [bracken_n_species_called: n_called]
}
