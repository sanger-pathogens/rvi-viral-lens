include {KRAKEN2BRACKEN} from '../rvi_toolbox/subworkflows/kraken2bracken.nf'
include {MSWEEP} from '../rvi_toolbox/modules/msweep.nf'
include {CLEANUP_THEMISTO_PSEUDOALIGNMENTS} from '../rvi_toolbox/modules/cleanup.nf'
include {ABUNDANCE_ESTIMATION} from '../rvi_toolbox/subworkflows/abundance_estimation.nf'
include {SCRUB_DECONTAM} from '../rvi_toolbox/subworkflows/scrub.nf'
include {GENERATE_ABUNDANCE_REPORT} from '../workflows/GENERATE_ABUNDANCE_REPORT.nf'
include {publish_lane_json as publish_abundance_lane_json} from '../modules/publish_lane_report.nf'
include {publish_run_files as publish_abundance_run_files} from '../modules/publish_lite.nf'

workflow ABUNDANCE {

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

                // SCRuB runs once per run against the whole-run abundance_summary
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

            // ABUNDANCE_ESTIMATION has no emit: block, so nothing per-sample is reachable from it the way KRAKEN2BRACKEN's fork now is. 
            abund_est_ran_ch = abundance_sample_ch
                .map { id, _meta -> [id, [abundance_estimation_ran: true]] }
        } else {
            abund_est_ran_ch = Channel.empty()
        }

        // -- mSWEEP probabilistic abundance estimation (opt-in, --run_msweep) -----------
        // Off by default. mSWEEP reads Themisto2's pseudoalignments, so it is only meaningful when
        // the sequence-index lane actually ran Themisto2. 
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

            // 2-tuple, not 3: MSWEEP no longer writes the per-read probability matrix
            // (--write-probs), which cost tens of GB per sample and nothing consumed.
            msweep_counts_ch = MSWEEP.out.abundances
                .map { meta, abundances -> [meta.id, count_msweep_abundances(abundances)] }

            // Themisto2's pseudoalignments clean them up
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
                def new_meta = meta + (b_counts ?: empty_bracken_counts()) +
                    (s_ran ?: [scrub_ran: false]) + (ae_ran ?: [abundance_estimation_ran: false]) +
                    (m_counts ?: empty_msweep_counts())
                [id, new_meta]
            }
            .set { abundance_report_prep_ch }

        GENERATE_ABUNDANCE_REPORT(abundance_report_prep_ch)

        // PUBLISH (abundance lane)
        publish_abundance_lane_json(GENERATE_ABUNDANCE_REPORT.out.publish_seq_level_ch)
        publish_abundance_run_files(GENERATE_ABUNDANCE_REPORT.out.publish_run_level_summaries_ch)
}

// --- sample-level count helper for the abundance report ---

// A step that never ran and a step that ran and found nothing are different facts; the
// report spells the first NA and the second 0, rather than conflating them (same rule as
// subworkflows/classifying_index.nf).
NOT_RUN = 'NA'

def empty_bracken_counts() {
    return [bracken_n_species_called: params.run_kraken2bracken ? 0 : NOT_RUN]
}

def empty_msweep_counts() {
    if (params.run_msweep) return EMPTY_MSWEEP_COUNTS
    return [msweep_n_groups: NOT_RUN, msweep_top_group: NOT_RUN, msweep_top_abundance: NOT_RUN]
}

EMPTY_BRACKEN_COUNTS = [bracken_n_species_called: 0]
EMPTY_MSWEEP_COUNTS = [msweep_n_groups: 0, msweep_top_group: '', msweep_top_abundance: 0.0]

def count_msweep_abundances(txt) {
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
    // Counts species-level rows with a non-zero read count.
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
