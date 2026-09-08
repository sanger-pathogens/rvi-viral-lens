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
// -- new-species consensus (rvi_integration_1, opt-in via --call_consensus_for_new_species):
// when a sequence-index method calls a species above new_species_min_breadth_pct breadth
// that MAPPING's Kraken2/SORT_READS_BY_REF did NOT already find for that sample, run it
// through GENERATE_CONSENSUS too, same as any Kraken2-found taxid.
include {INDEX_REFERENCE_FASTA; EXTRACT_REFERENCE_SUBSET} from '../modules/reference_subset.nf'
include {SELECT_REFERENCE_RECORD_BY_NAME} from '../modules/select_reference_record_by_name.nf'
include {GENERATE_CONSENSUS} from '../workflows/GENERATE_CONSENSUS.nf'
include {publish_consensus_files as publish_new_species_consensus_files} from '../modules/publish_lite.nf'

workflow SEQUENCE_INDEX {
    take:
        preprocessed_3tuple_ch  // tuple (meta, read1, read2)
        identified_species_ch   // [sample_id, [normalized_species_name, ...]] -- MAPPING.out.identified_species_ch

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

            // Kept as its own variable so the new-species block below can consume it
            // without touching VIRAL_MSWEEP.out, which is undefined unless the
            // subworkflow was actually invoked.
            msweep_map_qc_ch = VIRAL_MSWEEP.out.map_qc
        } else {
            msweep_counts_ch = Channel.empty()
            mapqc_counts_ch = Channel.empty()
            msweep_map_qc_ch = Channel.empty()
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
                .map { meta, tsv -> [meta.id, count_metagraph_species_hits(tsv, 'metagraph_query')] }

            metagraph_query_mapqc_counts_ch = VIRAL_METAGRAPH_QUERY.out.map_qc
                .map { meta, tsv -> [meta.id, count_metagraph_map_qc(tsv, 'metagraph_query')] }

            metagraph_query_map_qc_ch = VIRAL_METAGRAPH_QUERY.out.map_qc
        } else {
            metagraph_query_counts_ch = Channel.empty()
            metagraph_query_mapqc_counts_ch = Channel.empty()
            metagraph_query_map_qc_ch = Channel.empty()
        }

        // -- New-species consensus (opt-in): species a sequence-index method called with
        // real breadth of coverage that MAPPING's Kraken2 pass never found for that
        // sample. See subworkflows/mapping.nf's identified_species_ch header comment for
        // why this compares by free-text species name (the only thing Kraken2 taxids and
        // the mSWEEP/Metagraph reference indexes' own labels have in common) and why the
        // per-sample "already identified" set needs no batch-wide wait.
        if (params.call_consensus_for_new_species) {
            // Union every method's own already-computed map_qc breadth table, filtered to
            // real hits (breadth_pct > new_species_min_breadth_pct). These consume the
            // *_map_qc_ch variables set in each method's if/else above, NOT
            // VIRAL_*.out.map_qc directly: a subworkflow that was never invoked has no
            // .out at all, so reaching for it aborts the run with "Access to
            // 'VIRAL_METAGRAPH_ALIGN.out' is undefined" the moment this feature is enabled
            // with any subset of the three methods. flatMap over an empty channel emits
            // nothing, which is what "that method is off" should mean here.
            msweep_candidates_ch = msweep_map_qc_ch
                .flatMap { meta, tsv -> parse_new_species_candidates(tsv, 'species_label').collect { name -> [meta.id, name] } }

            metagraph_align_candidates_ch = metagraph_align_map_qc_ch
                .flatMap { meta, tsv -> parse_new_species_candidates(tsv, 'species').collect { name -> [meta.id, name] } }

            metagraph_query_candidates_ch = metagraph_query_map_qc_ch
                .flatMap { meta, tsv -> parse_new_species_candidates(tsv, 'species').collect { name -> [meta.id, name] } }

            // .unique() streams -- it emits each non-duplicate immediately as it passes,
            // it does not need to see the whole channel close first (unlike groupTuple()).
            candidate_new_species_ch = msweep_candidates_ch
                .mix(metagraph_align_candidates_ch, metagraph_query_candidates_ch)
                .unique { sample_id, name -> [sample_id, name.trim().toLowerCase()] }

            // Drop anything MAPPING already found for that sample. remainder:true so a
            // sample with no MAPPING entry at all (e.g. every fastq filtered as empty
            // upstream) still passes its candidates through -- treated as "nothing
            // already identified", not as "drop everything".
            candidate_new_species_ch
                .map { sample_id, name -> [sample_id, name, name.trim().toLowerCase()] }
                .join(identified_species_ch, remainder: true)
                .filter { _sample_id, _name, name_norm, identified -> !(identified != null && identified.contains(name_norm)) }
                .map { sample_id, name, _name_norm, _identified -> [sample_id, name] }
                .set { new_species_ch }

            // Build the synthetic per-(sample,species) meta GENERATE_CONSENSUS/its publish
            // step need: taxid here is a filesystem-safe slug of the species name, NOT a
            // real Kraken taxid -- there isn't one, these species were never Kraken2-sorted.
            new_species_meta_ch = new_species_ch
                .map { sample_id, species_name ->
                    def slug = species_name.replaceAll(/[^A-Za-z0-9]+/, '_').replaceAll(/^_+|_+$/, '')
                    def meta = [
                        id: "${sample_id}.${slug}",
                        sample_id: sample_id,
                        taxid: slug,
                        species_name: species_name,
                        discovered_by: 'sequence_index',
                    ]
                    [meta, species_name]
                }

            // Same reference (msweep_ref_groups/msweep_map_reference_fasta) mSWEEP's own
            // map_qc already indexes, reused regardless of run_msweep -- both params always
            // have real defaults (see nextflow.config).
            INDEX_REFERENCE_FASTA(Channel.fromPath(params.msweep_map_reference_fasta))
            new_species_indexed_reference_ch = INDEX_REFERENCE_FASTA.out.fasta.first()
            new_species_sequence_lengths_ch  = INDEX_REFERENCE_FASTA.out.lengths.first()
            new_species_labels_ch            = Channel.fromPath(params.msweep_ref_groups).first()

            SELECT_REFERENCE_RECORD_BY_NAME(new_species_meta_ch, new_species_labels_ch, new_species_sequence_lengths_ch)
            EXTRACT_REFERENCE_SUBSET(SELECT_REFERENCE_RECORD_BY_NAME.out.record_id, new_species_indexed_reference_ch)

            // Deliberately no cross-sample de-duplication here: a species found "new" in
            // many samples gets its reference extracted once per sample rather than once
            // per run -- keeps everything scoped per-sample (no batch-wide wait) at the
            // cost of some redundant extraction work.
            reads_by_sample_ch = preprocessed_3tuple_ch
                .map { meta, r1, r2 -> [meta.id, r1, r2] } // meta.id == sample_id here (pre-lane)

            EXTRACT_REFERENCE_SUBSET.out.subset_fasta
                .map { meta, ref_fa -> [meta.sample_id, meta, ref_fa] }
                .combine(reads_by_sample_ch, by: 0)
                .map { _sample_id, meta, ref_fa, r1, r2 -> [meta, [r1, r2], ref_fa] }
                .set { new_species_consensus_in_ch }

            GENERATE_CONSENSUS(new_species_consensus_in_ch)

            GENERATE_CONSENSUS.out.filtered_consensus_ch
                .map { meta, bam, bam_idx, consensus, _qc_json -> [meta, [bam, bam_idx, consensus]] }
                .set { new_species_aln_publish_ch }

            publish_new_species_consensus_files(new_species_aln_publish_ch)

            new_species_counts_ch = GENERATE_CONSENSUS.out.filtered_consensus_ch
                .map { meta, _bam, _bam_idx, _consensus, _qc_json -> [meta.sample_id, 1] }
                .groupTuple()
                .map { sample_id, ones -> [sample_id, [new_species_consensus_n: ones.size()]] }
        } else {
            new_species_counts_ch = Channel.empty()
        }

        sequence_index_sample_ch
            .join(msweep_counts_ch, remainder: true)
            .join(mapqc_counts_ch, remainder: true)
            .join(metagraph_align_counts_ch, remainder: true)
            .join(metagraph_align_mapqc_counts_ch, remainder: true)
            .join(metagraph_query_counts_ch, remainder: true)
            .join(metagraph_query_mapqc_counts_ch, remainder: true)
            .join(new_species_counts_ch, remainder: true)
            .map { id, meta, m_counts, qc_counts, mga_counts, mga_qc_counts, mgq_counts, mgq_qc_counts, ns_counts ->
                def new_meta = meta + (m_counts ?: EMPTY_MSWEEP_COUNTS) + (qc_counts ?: EMPTY_MAP_QC_COUNTS) +
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
// New-species consensus (--call_consensus_for_new_species): how many species a
// sequence-index method called, with real breadth, that MAPPING's Kraken2 pass didn't
// already find for that sample -- 0 whenever the feature is off, or on but nothing new
// was found for this sample.
EMPTY_NEW_SPECIES_COUNTS = [new_species_consensus_n: 0]

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

def parse_new_species_candidates(tsv, species_col) {
    // A method's own *_map_qc.tsv (msweep_map_qc.tsv: species_label; metagraph_map_qc.tsv,
    // both align and query: species) -- one row per species that method already validated
    // by mapping. Returns the species names (original case, for reference-record lookup)
    // whose breadth_pct clears new_species_min_breadth_pct.
    if (tsv == null || !tsv.exists()) return []
    def lines = tsv.readLines()
    if (lines.size() < 2) return []
    def header = lines[0].split('\t')
    def species_idx = header.findIndexOf { String col -> col == species_col }
    def breadth_idx = header.findIndexOf { String col -> col == 'breadth_pct' }
    if (species_idx < 0 || breadth_idx < 0) return []
    def candidates = []
    lines[1..-1].each { line ->
        def cols = line.split('\t')
        if (species_idx >= cols.size() || breadth_idx >= cols.size()) return
        try {
            if ((cols[breadth_idx] as Double) > params.new_species_min_breadth_pct) {
                candidates << cols[species_idx]
            }
        } catch (NumberFormatException ignored) {
            // header or malformed row -- skipped
        }
    }
    return candidates
}
