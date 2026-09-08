// --- map reads to sequence indexes (rvi_integration_1) ----------------------
// Extracted unchanged from main.nf's inline body. Up to three methods run in
// parallel off the same preprocessed reads (not downstream of one another),
// each independently gated, all feeding ONE GENERATE_MAPPING_REPORT call --
// see INSTRUCT.md item 3: "feed its per-sample counts into the same
// mapping_report_prep_ch rather than building a second report path".
// sequence_index_sample_ch is the join backbone (every sample that reaches
// this lane) so a sample report row exists even if only one method ran for
// it, or (with more than one flag on) one row carries every enabled method's
// counts.
include {VIRAL_MSWEEP} from '../workflows/VIRAL_MSWEEP.nf'
include {VIRAL_METAGRAPH_ALIGN} from '../workflows/VIRAL_METAGRAPH_ALIGN.nf'
include {VIRAL_METAGRAPH_QUERY} from '../workflows/VIRAL_METAGRAPH_QUERY.nf'
include {GENERATE_MAPPING_REPORT} from '../workflows/GENERATE_MAPPING_REPORT.nf'
include {publish_lane_json as publish_mapping_lane_json} from '../modules/publish_lane_report.nf'
include {publish_run_files as publish_mapping_run_files} from '../modules/publish_lite.nf'

workflow SEQUENCE_INDEX {
    take:
        preprocessed_3tuple_ch // tuple (meta, read1, read2)

    main:
        sequence_index_sample_ch = preprocessed_3tuple_ch
            .map { meta, _r1, _r2 -> [meta.id, meta] }

        // -- Themisto2 pseudoalignment + mSWEEP abundance, then breadth-of-coverage
        // validation of the low-abundance calls.
        if (params.run_msweep) {
            VIRAL_MSWEEP(preprocessed_3tuple_ch)

            msweep_counts_ch = VIRAL_MSWEEP.out.abundances
                .map { meta, abundances, _probs -> [meta.id, count_msweep_abundances(abundances)] }

            // map_qc drops samples with nothing above msweep_map_min_abundance (the
            // MSWEEP_MAP_QC processes emit optional outputs); joined with remainder
            // below rather than losing those samples from the report entirely.
            mapqc_counts_ch = VIRAL_MSWEEP.out.map_qc
                .map { meta, qc_tsv -> [meta.id, count_msweep_map_qc(qc_tsv)] }
        } else {
            msweep_counts_ch = Channel.empty()
            mapqc_counts_ch = Channel.empty()
        }

        // -- Sequence-to-graph alignment via Metagraph (metagraph align), then its own
        // map-QC validation.
        if (params.run_metagraph_align) {
            VIRAL_METAGRAPH_ALIGN(preprocessed_3tuple_ch)

            metagraph_align_counts_ch = VIRAL_METAGRAPH_ALIGN.out.species_hits
                .map { meta, tsv -> [meta.id, count_metagraph_species_hits(tsv, 'metagraph_align')] }

            // map_qc drops samples with nothing above metagraph_align_min_hits, same
            // reasoning as mSWEEP's map_qc above.
            metagraph_align_mapqc_counts_ch = VIRAL_METAGRAPH_ALIGN.out.map_qc
                .map { meta, tsv -> [meta.id, count_metagraph_map_qc(tsv, 'metagraph_align')] }
        } else {
            metagraph_align_counts_ch = Channel.empty()
            metagraph_align_mapqc_counts_ch = Channel.empty()
        }

        // -- Pseudoalignment via Metagraph (metagraph query --query-mode labels), same
        // shared index, alternative method to the alignment above -- see
        // ../workflows/VIRAL_METAGRAPH_QUERY.nf / ../modules/metagraph_query.nf for why
        // this exists as its own module rather than the old filter+query pipeline that
        // metagraph_align.nf itself replaced (found ~zero real hits, see git history).
        if (params.run_metagraph_query) {
            VIRAL_METAGRAPH_QUERY(preprocessed_3tuple_ch)

            metagraph_query_counts_ch = VIRAL_METAGRAPH_QUERY.out.species_hits
                .map { meta, tsv -> [meta.id, count_metagraph_species_hits(tsv, 'metagraph_query')] }

            metagraph_query_mapqc_counts_ch = VIRAL_METAGRAPH_QUERY.out.map_qc
                .map { meta, tsv -> [meta.id, count_metagraph_map_qc(tsv, 'metagraph_query')] }
        } else {
            metagraph_query_counts_ch = Channel.empty()
            metagraph_query_mapqc_counts_ch = Channel.empty()
        }

        sequence_index_sample_ch
            .join(msweep_counts_ch, remainder: true)
            .join(mapqc_counts_ch, remainder: true)
            .join(metagraph_align_counts_ch, remainder: true)
            .join(metagraph_align_mapqc_counts_ch, remainder: true)
            .join(metagraph_query_counts_ch, remainder: true)
            .join(metagraph_query_mapqc_counts_ch, remainder: true)
            .map { id, meta, m_counts, qc_counts, mga_counts, mga_qc_counts, mgq_counts, mgq_qc_counts ->
                def new_meta = meta + (m_counts ?: EMPTY_MSWEEP_COUNTS) + (qc_counts ?: EMPTY_MAP_QC_COUNTS) +
                    (mga_counts ?: EMPTY_METAGRAPH_ALIGN_COUNTS) + (mga_qc_counts ?: EMPTY_METAGRAPH_ALIGN_MAPQC_COUNTS) +
                    (mgq_counts ?: EMPTY_METAGRAPH_QUERY_COUNTS) + (mgq_qc_counts ?: EMPTY_METAGRAPH_QUERY_MAPQC_COUNTS)
                [id, new_meta]
            }
            .set { mapping_report_prep_ch }

        GENERATE_MAPPING_REPORT(mapping_report_prep_ch)

        // PUBLISH (mapping/sequence-index lane)
        publish_mapping_lane_json(GENERATE_MAPPING_REPORT.out.publish_seq_level_ch)
        publish_mapping_run_files(GENERATE_MAPPING_REPORT.out.publish_run_level_summaries_ch)
}

// --- rvi_integration_1: sample-level count helpers for the mapping report ---
// Default counts for a sample a given optional step produced no output for -- either
// because that method didn't run at all (join remainder is null) or because the method
// ran but the step's own output is itself optional per-sample (e.g. nothing above a
// min-abundance/min-hits threshold). Named constants (not inline [:]) so every sample
// still gets the same report columns regardless of which method(s) actually ran for it.
EMPTY_MSWEEP_COUNTS = [msweep_n_groups: 0, msweep_top_group: '', msweep_top_abundance: 0.0]
EMPTY_MAP_QC_COUNTS = [mapqc_n_species: 0, mapqc_max_breadth_pct: 0.0]
// One pair per Metagraph method (align, query) -- both call count_metagraph_species_hits()/
// count_metagraph_map_qc() with a distinct prefix, since both methods' counts can merge
// into the same per-sample meta and would otherwise collide on field name.
EMPTY_METAGRAPH_ALIGN_COUNTS = [metagraph_align_n_species_considered: 0, metagraph_align_n_species_called: 0]
EMPTY_METAGRAPH_ALIGN_MAPQC_COUNTS = [metagraph_align_mapqc_n_species: 0, metagraph_align_mapqc_max_breadth_pct: 0.0]
EMPTY_METAGRAPH_QUERY_COUNTS = [metagraph_query_n_species_considered: 0, metagraph_query_n_species_called: 0]
EMPTY_METAGRAPH_QUERY_MAPQC_COUNTS = [metagraph_query_mapqc_n_species: 0, metagraph_query_mapqc_max_breadth_pct: 0.0]

def empty_metagraph_counts(prefix) {
    return prefix == 'metagraph_align' ? EMPTY_METAGRAPH_ALIGN_COUNTS : EMPTY_METAGRAPH_QUERY_COUNTS
}

def empty_metagraph_mapqc_counts(prefix) {
    return prefix == 'metagraph_align' ? EMPTY_METAGRAPH_ALIGN_MAPQC_COUNTS : EMPTY_METAGRAPH_QUERY_MAPQC_COUNTS
}

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

def count_metagraph_species_hits(tsv, prefix) {
    // <sample>_species_hits.tsv (bin/call_metagraph_species.py): sample_id, species,
    // hit_count, provisional_call -- the last written by Python's str(bool), so "True"/
    // "False", not lowercase. Shared by both Metagraph methods (see
    // ../workflows/VIRAL_METAGRAPH_ALIGN.nf / VIRAL_METAGRAPH_QUERY.nf, both call
    // CALL_METAGRAPH_SPECIES) -- prefix ('metagraph_align' or 'metagraph_query') keeps
    // their counts from colliding when both merge into the same per-sample meta.
    if (tsv == null || !tsv.exists()) return empty_metagraph_counts(prefix)
    def lines = tsv.readLines()
    if (lines.size() < 2) return empty_metagraph_counts(prefix)
    def header = lines[0].split('\t')
    def called_idx = header.findIndexOf { String col -> col == 'provisional_call' }
    def n_considered = lines.size() - 1
    def n_called = 0
    if (called_idx >= 0) {
        n_called = lines[1..-1].count { String line ->
            def cols = line.split('\t')
            called_idx < cols.size() && cols[called_idx] == 'True'
        }
    }
    return ["${prefix}_n_species_considered": n_considered, "${prefix}_n_species_called": n_called]
}

def count_metagraph_map_qc(tsv, prefix) {
    // <sample>_metagraph_map_qc.tsv (bin/aggregate_metagraph_coverage.py): sample_id,
    // species, hit_count, reference_accession, reference_length, query_length,
    // covered_bases, breadth_pct, mean_depth, meanbaseq, meanmapq, reads_mapped. Same
    // prefix reasoning as count_metagraph_species_hits() above.
    if (tsv == null || !tsv.exists()) return empty_metagraph_mapqc_counts(prefix)
    def lines = tsv.readLines()
    if (lines.size() < 2) return empty_metagraph_mapqc_counts(prefix)
    def header = lines[0].split('\t')
    def breadth_idx = header.findIndexOf { String col -> col == 'breadth_pct' }
    def breadths = lines[1..-1].collect { String line ->
        def cols = line.split('\t')
        if (breadth_idx < 0 || breadth_idx >= cols.size()) return 0.0
        try { return cols[breadth_idx] as Double } catch (NumberFormatException ignored) { return 0.0 }
    }
    return [
        "${prefix}_mapqc_n_species":       lines.size() - 1,
        "${prefix}_mapqc_max_breadth_pct": breadths ? breadths.max() : 0.0
    ]
}

def count_msweep_map_qc(tsv) {
    // <sample>_msweep_map_qc.tsv, written by bin/aggregate_species_coverage.py:
    // sample_id, species_label, relative_abundance, reference_length, query_length,
    // covered_bases, breadth_pct, mean_depth, meanbaseq, meanmapq, reads_mapped
    if (tsv == null || !tsv.exists()) return EMPTY_MAP_QC_COUNTS
    def lines = tsv.readLines()
    if (lines.size() < 2) return EMPTY_MAP_QC_COUNTS
    def header = lines[0].split('\t')
    def breadth_idx = header.findIndexOf { String col -> col == 'breadth_pct' }
    def breadths = lines[1..-1].collect { line ->
        def cols = line.split('\t')
        if (breadth_idx < 0 || breadth_idx >= cols.size()) return 0.0
        try { return cols[breadth_idx] as Double } catch (NumberFormatException ignored) { return 0.0 }
    }
    return [
        mapqc_n_species:       lines.size() - 1,
        mapqc_max_breadth_pct: breadths ? breadths.max() : 0.0
    ]
}
