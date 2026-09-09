// --- classify reads against pre-built sequence indexes (rvi_integration_1) ----
// The sequence-index counterpart to subworkflows/classifying_kraken2.nf: up to three
// methods run in parallel off the same preprocessed reads (not downstream of one
// another), each independently gated, all feeding ONE GENERATE_MAPPING_REPORT call --
// see INSTRUCT.md item 3: "feed its per-sample counts into the same
// mapping_report_prep_ch rather than building a second report path".
// sequence_index_sample_ch is the join backbone (every sample that reaches this lane) so
// a sample report row exists even if only one method ran for it, or (with more than one
// flag on) one row carries every enabled method's counts.
//
// This classifier CALLS AND REPORTS ONLY. It does no reference extraction, no read
// mapping and no consensus generation: it emits the species each method called, the
// reference record that method validated each one against, and the supporting hit
// counts, and subworkflows/mapping.nf takes it from there -- preferring Kraken2's calls
// where the two classifiers agree, and mapping the genuinely new ones itself. Doing the
// reference/read work here as well was the previous arrangement, and meant resolving
// references for species that were about to be discarded as already-known.
include {VIRAL_THEMISTO_MSWEEP} from '../workflows/VIRAL_THEMISTO_MSWEEP.nf'
include {VIRAL_METAGRAPH_ALIGN} from '../workflows/VIRAL_METAGRAPH_ALIGN.nf'
include {VIRAL_METAGRAPH_QUERY} from '../workflows/VIRAL_METAGRAPH_QUERY.nf'
include {GENERATE_MAPPING_REPORT} from '../workflows/GENERATE_MAPPING_REPORT.nf'
include {publish_lane_json as publish_mapping_lane_json} from '../modules/publish_lane_report.nf'
include {publish_run_files as publish_mapping_run_files} from '../modules/publish_lite.nf'

workflow CLASSIFYING_INDEX {
    take:
        preprocessed_3tuple_ch  // tuple (meta, read1, read2)

    main:
        sequence_index_sample_ch = preprocessed_3tuple_ch
            .map { meta, _r1, _r2 -> [meta.id, meta] }

        // -- Themisto2 pseudoalignment + mSWEEP abundance, then breadth-of-coverage
        // validation of the low-abundance calls.
        //
        // THE LANE'S DEFAULT METHOD (run_themisto defaults true). Species are called
        // directly from Themisto2 pseudoalignment read-hit counts (CALL_THEMISTO_SPECIES)
        // and validated by THEMISTO_MAP_QC -- no probabilistic model involved. mSWEEP's
        // abundance estimate is an optional add-on *inside* this arm, gated by run_msweep
        // (default false) inside VIRAL_THEMISTO_MSWEEP itself; see nextflow.config's note
        // on run_msweep's changed meaning.
        if (params.run_themisto) {
            VIRAL_THEMISTO_MSWEEP(preprocessed_3tuple_ch)

            themisto_counts_ch = VIRAL_THEMISTO_MSWEEP.out.species_hits
                .map { meta, tsv -> [meta.id, count_species_hits(tsv, 'themisto')] }

            // themisto_map_qc is Channel.empty() when themisto_align_run_map_qc is off,
            // and drops samples with nothing above themisto_align_min_hits; joined with
            // remainder below rather than losing those samples from the report.
            themisto_mapqc_counts_ch = VIRAL_THEMISTO_MSWEEP.out.themisto_map_qc
                .map { meta, tsv -> [meta.id, count_map_qc_breadth(tsv, 'themisto')] }

            // Both of these are Channel.empty() unless run_msweep is set (see
            // ../workflows/VIRAL_THEMISTO_MSWEEP.nf), so they need no gate of their own --
            // an empty channel simply contributes no counts and the joins below fill in
            // EMPTY_*_COUNTS.
            msweep_counts_ch = VIRAL_THEMISTO_MSWEEP.out.abundances
                .map { meta, abundances, _probs -> [meta.id, count_msweep_abundances(abundances)] }

            // Kept as its own variable so the new-species block below can consume it
            // without touching VIRAL_THEMISTO_MSWEEP.out, which is undefined unless the
            // subworkflow was actually invoked.
            themisto_map_qc_ch = VIRAL_THEMISTO_MSWEEP.out.themisto_map_qc
        } else {
            themisto_counts_ch       = Channel.empty()
            themisto_mapqc_counts_ch = Channel.empty()
            msweep_counts_ch         = Channel.empty()
            themisto_map_qc_ch       = Channel.empty()
        }

        // -- Sequence-to-graph alignment via Metagraph (metagraph align), then its own
        // map-QC validation.
        if (params.run_metagraph_align) {
            VIRAL_METAGRAPH_ALIGN(preprocessed_3tuple_ch)

            metagraph_align_counts_ch = VIRAL_METAGRAPH_ALIGN.out.species_hits
                .map { meta, tsv -> [meta.id, count_species_hits(tsv, 'metagraph_align')] }

            // map_qc drops samples with nothing above metagraph_align_min_hits, same
            // reasoning as mSWEEP's map_qc above.
            metagraph_align_mapqc_counts_ch = VIRAL_METAGRAPH_ALIGN.out.map_qc
                .map { meta, tsv -> [meta.id, count_map_qc_breadth(tsv, 'metagraph_align')] }

            metagraph_align_map_qc_ch = VIRAL_METAGRAPH_ALIGN.out.map_qc
        } else {
            metagraph_align_counts_ch = Channel.empty()
            metagraph_align_mapqc_counts_ch = Channel.empty()
            metagraph_align_map_qc_ch = Channel.empty()
        }

        // -- Pseudoalignment via Metagraph (metagraph query --query-mode labels), same
        // shared index, alternative method to the alignment above -- see
        // ../workflows/VIRAL_METAGRAPH_QUERY.nf / ../modules/metagraph_query.nf for why
        // this exists as its own module rather than the old filter+query pipeline that
        // metagraph_align.nf itself replaced (found ~zero real hits, see git history).
        if (params.run_metagraph_query) {
            VIRAL_METAGRAPH_QUERY(preprocessed_3tuple_ch)

            metagraph_query_counts_ch = VIRAL_METAGRAPH_QUERY.out.species_hits
                .map { meta, tsv -> [meta.id, count_species_hits(tsv, 'metagraph_query')] }

            metagraph_query_mapqc_counts_ch = VIRAL_METAGRAPH_QUERY.out.map_qc
                .map { meta, tsv -> [meta.id, count_map_qc_breadth(tsv, 'metagraph_query')] }

            metagraph_query_map_qc_ch = VIRAL_METAGRAPH_QUERY.out.map_qc
        } else {
            metagraph_query_counts_ch = Channel.empty()
            metagraph_query_mapqc_counts_ch = Channel.empty()
            metagraph_query_map_qc_ch = Channel.empty()
        }

        // -- Species calls handed to MAPPING (opt-in): every species a sequence-index
        // method called with real breadth of coverage, plus the reference record that
        // method validated it against and its read-hit count.
        //
        // Reporting ONLY. This classifier does no reference extraction and no read
        // mapping: MAPPING decides which of these species are actually new (Kraken2's
        // calls win) and resolves/maps just those. See MAPPING's own comment for why the
        // filter lives there.
        if (params.call_consensus_for_new_species) {
            // Parse each enabled method's own already-computed map-QC table. These consume
            // the *_map_qc_ch variables set in each method's if/else above, NOT
            // VIRAL_*.out.map_qc directly: a subworkflow that was never invoked has no
            // .out at all, so reaching for it aborts the run with "Access to
            // 'VIRAL_METAGRAPH_ALIGN.out' is undefined" the moment this feature is enabled
            // with any subset of the three methods. flatMap over an empty channel emits
            // nothing, which is what "that method is off" should mean here.
            themisto_calls_ch = themisto_map_qc_ch
                .flatMap { meta, tsv -> parse_species_calls(tsv, 'themisto').collect { call -> [meta.id, call] } }

            metagraph_align_calls_ch = metagraph_align_map_qc_ch
                .flatMap { meta, tsv -> parse_species_calls(tsv, 'metagraph_align').collect { call -> [meta.id, call] } }

            metagraph_query_calls_ch = metagraph_query_map_qc_ch
                .flatMap { meta, tsv -> parse_species_calls(tsv, 'metagraph_query').collect { call -> [meta.id, call] } }

            // One call per (sample, species). .unique() streams -- it emits each
            // non-duplicate immediately as it passes, it does not need to see the whole
            // channel close first (unlike groupTuple()).
            //
            // CAVEAT with more than one method enabled: which method's row wins here, and
            // therefore which reference record and hit count get attributed to a species
            // both methods called, is whichever arrives first -- task completion order, so
            // not reproducible run to run. Harmless on the defaults (only run_themisto is
            // on, so there is nothing to race), and the species itself is unaffected --
            // only the record and counts reported for it. If a multi-method run ever needs
            // determinism here, rank by method instead of by arrival.
            //
            // No mSWEEP calls: mSWEEP now estimates abundance only and produces no breadth
            // table to threshold on (see ../workflows/VIRAL_THEMISTO_MSWEEP.nf).
            species_calls_ch = themisto_calls_ch
                .mix(metagraph_align_calls_ch, metagraph_query_calls_ch)
                .unique { sample_id, call -> [sample_id, call.species_name.trim().toLowerCase()] }

            // Counts what this classifier can honestly measure: species it called above the
            // breadth threshold and reported to MAPPING, BEFORE MAPPING drops the ones
            // Kraken2 already found. For the post-filter truth, count
            // classification-report rows carrying `discovered_by: 'sequence_index'`.
            new_species_counts_ch = species_calls_ch
                .map { sample_id, _call -> [sample_id, 1] }
                .groupTuple()
                .map { sample_id, ones -> [sample_id, [new_species_candidates_n: ones.size()]] }
        } else {
            species_calls_ch = Channel.empty()
            new_species_counts_ch = Channel.empty()
        }

        sequence_index_sample_ch
            .join(themisto_counts_ch, remainder: true)
            .join(themisto_mapqc_counts_ch, remainder: true)
            .join(msweep_counts_ch, remainder: true)
            .join(metagraph_align_counts_ch, remainder: true)
            .join(metagraph_align_mapqc_counts_ch, remainder: true)
            .join(metagraph_query_counts_ch, remainder: true)
            .join(metagraph_query_mapqc_counts_ch, remainder: true)
            .join(new_species_counts_ch, remainder: true)
            .map { id, meta, t_counts, t_qc_counts, m_counts, mga_counts, mga_qc_counts, mgq_counts, mgq_qc_counts, ns_counts ->
                def new_meta = meta + (t_counts ?: EMPTY_THEMISTO_COUNTS) + (t_qc_counts ?: EMPTY_THEMISTO_MAPQC_COUNTS) +
                    (m_counts ?: EMPTY_MSWEEP_COUNTS) +
                    (mga_counts ?: EMPTY_METAGRAPH_ALIGN_COUNTS) + (mga_qc_counts ?: EMPTY_METAGRAPH_ALIGN_MAPQC_COUNTS) +
                    (mgq_counts ?: EMPTY_METAGRAPH_QUERY_COUNTS) + (mgq_qc_counts ?: EMPTY_METAGRAPH_QUERY_MAPQC_COUNTS) +
                    (ns_counts ?: EMPTY_NEW_SPECIES_COUNTS)
                [id, new_meta]
            }
            .set { mapping_report_prep_ch }

        GENERATE_MAPPING_REPORT(mapping_report_prep_ch)

        // PUBLISH (mapping/sequence-index lane)
        publish_mapping_lane_json(GENERATE_MAPPING_REPORT.out.publish_seq_level_ch)
        publish_mapping_run_files(GENERATE_MAPPING_REPORT.out.publish_run_level_summaries_ch)

    emit:
        // Handover to subworkflows/mapping.nf: one entry per (sample, species) this lane
        // called above the breadth threshold, as [sample_id, call] where call is a Map of
        // species_name, reference_record (a SEQIDX_<n> token from the method's own map-QC
        // table, the id INDEX_REFERENCE_FASTA mints and seqkit grep matches), hit_count,
        // breadth_pct and method.
        //
        // Note this is NOT the shape classifying_kraken2.nf hands over -- that one emits
        // reads-plus-reference ready for consensus, because SORT_READS_BY_REF resolves its
        // own references upstream. The asymmetry is deliberate: it is what lets MAPPING
        // resolve references only for species that survive its filter. MAPPING builds the
        // consensus-ready shape for this side itself.
        //
        // Channel.empty() unless --call_consensus_for_new_species is set, so MAPPING can
        // consume it unconditionally.
        species_calls_ch // [sample_id, [species_name:, reference_record:, hit_count:, breadth_pct:, method:]]
}

// --- rvi_integration_1: sample-level count helpers for the mapping report ---
// Default counts for a sample a given optional step produced no output for -- either
// because that method didn't run at all (join remainder is null) or because the method
// ran but the step's own output is itself optional per-sample (e.g. nothing above a
// min-abundance/min-hits threshold). Named constants (not inline [:]) so every sample
// still gets the same report columns regardless of which method(s) actually ran for it.
EMPTY_THEMISTO_COUNTS = [themisto_n_species_considered: 0, themisto_n_species_called: 0]
EMPTY_THEMISTO_MAPQC_COUNTS = [themisto_mapqc_n_species: 0, themisto_mapqc_max_breadth_pct: 0.0]
EMPTY_MSWEEP_COUNTS = [msweep_n_groups: 0, msweep_top_group: '', msweep_top_abundance: 0.0]
// One pair per Metagraph method (align, query) -- both call count_species_hits()/
// count_map_qc_breadth() with a distinct prefix, since both methods' counts can merge
// into the same per-sample meta and would otherwise collide on field name.
EMPTY_METAGRAPH_ALIGN_COUNTS = [metagraph_align_n_species_considered: 0, metagraph_align_n_species_called: 0]
EMPTY_METAGRAPH_ALIGN_MAPQC_COUNTS = [metagraph_align_mapqc_n_species: 0, metagraph_align_mapqc_max_breadth_pct: 0.0]
EMPTY_METAGRAPH_QUERY_COUNTS = [metagraph_query_n_species_considered: 0, metagraph_query_n_species_called: 0]
EMPTY_METAGRAPH_QUERY_MAPQC_COUNTS = [metagraph_query_mapqc_n_species: 0, metagraph_query_mapqc_max_breadth_pct: 0.0]
// New-species candidates (--call_consensus_for_new_species): how many species a
// sequence-index method called with real breadth and got a reference resolved for, i.e.
// how many were handed to MAPPING as consensus candidates -- 0 whenever the feature is
// off, or on but nothing cleared the threshold for this sample. Whether MAPPING then
// kept them (Kraken2 hadn't already found them) is not visible from here, by design:
// see the counts block above.
EMPTY_NEW_SPECIES_COUNTS = [new_species_candidates_n: 0]

def empty_species_hits_counts(prefix) {
    return prefix == 'metagraph_align' ? EMPTY_METAGRAPH_ALIGN_COUNTS : EMPTY_METAGRAPH_QUERY_COUNTS
}

def empty_map_qc_counts(prefix) {
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

def count_species_hits(tsv, prefix) {
    // <sample>_species_hits.tsv: sample_id, species, hit_count, provisional_call -- the
    // last written by Python's str(bool), so "True"/"False", not lowercase. All three
    // read-hit methods emit this identical schema (bin/call_metagraph_species.py for both
    // Metagraph methods, bin/call_themisto_species.py for Themisto2), so one parser serves
    // them all -- prefix ('themisto', 'metagraph_align' or 'metagraph_query') keeps their
    // counts from colliding when several merge into the same per-sample meta.
    if (tsv == null || !tsv.exists()) return empty_species_hits_counts(prefix)
    def lines = tsv.readLines()
    if (lines.size() < 2) return empty_species_hits_counts(prefix)
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

def count_map_qc_breadth(tsv, prefix) {
    // <sample>_metagraph_map_qc.tsv (bin/aggregate_metagraph_coverage.py): sample_id,
    // species, hit_count, reference_accession, reference_length, query_length,
    // covered_bases, breadth_pct, mean_depth, meanbaseq, meanmapq, reads_mapped. Same
    // prefix reasoning as count_species_hits() above.
    if (tsv == null || !tsv.exists()) return empty_map_qc_counts(prefix)
    def lines = tsv.readLines()
    if (lines.size() < 2) return empty_map_qc_counts(prefix)
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


def parse_species_calls(tsv, method) {
    // A method's own *_map_qc.tsv -- one row per species that method validated by mapping
    // reads against a single chosen reference record. Returns one call Map per row whose
    // breadth_pct clears new_species_min_breadth_pct, carrying everything MAPPING needs to
    // map that species without re-deriving anything: the species name (original case), the
    // reference record, and the supporting counts.
    //
    // The reference column is named differently by the two aggregators that write these
    // tables -- `reference_record` (bin/aggregate_themisto_coverage.py) vs
    // `reference_accession` (bin/aggregate_metagraph_coverage.py) -- for the same 4th
    // column holding the same SEQIDX_<n> token. Both spellings are accepted; a table with
    // neither fails loudly rather than silently reporting species with no reference, since
    // MAPPING cannot map those and the sample would just quietly lose them.
    if (tsv == null || !tsv.exists()) return []
    def lines = tsv.readLines()
    if (lines.size() < 2) return []
    def header = lines[0].split('\t')
    def species_idx = header.findIndexOf { String col -> col == 'species' }
    def breadth_idx = header.findIndexOf { String col -> col == 'breadth_pct' }
    def hits_idx    = header.findIndexOf { String col -> col == 'hit_count' }
    def record_idx  = header.findIndexOf { String col -> col == 'reference_record' }
    if (record_idx < 0) {
        record_idx = header.findIndexOf { String col -> col == 'reference_accession' }
    }
    if (species_idx < 0 || breadth_idx < 0) return []
    if (record_idx < 0) {
        error("map-QC table ${tsv} has neither a reference_record nor a reference_accession " +
              "column (header: ${header}). One of bin/aggregate_{themisto,metagraph}_coverage.py " +
              "changed its output -- update parse_species_calls() in " +
              "subworkflows/classifying_index.nf.")
    }
    def max_idx = [species_idx, breadth_idx, record_idx].max()
    def calls = []
    lines[1..-1].each { line ->
        def cols = line.split('\t')
        if (max_idx >= cols.size()) return
        try {
            def breadth = cols[breadth_idx] as Double
            if (breadth <= params.new_species_min_breadth_pct) return
            calls << [
                species_name:     cols[species_idx],
                reference_record: cols[record_idx],
                breadth_pct:      breadth,
                hit_count:        (hits_idx >= 0 && hits_idx < cols.size()) ? cols[hits_idx] : '',
                method:           method,
            ]
        } catch (NumberFormatException ignored) {
            // header or malformed row -- skipped
        }
    }
    return calls
}
