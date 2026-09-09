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
// What this classifier produces is CALLS, not consensus sequences: per sample, the species
// each method found, the reference record its own index points at for that species, and
// the read-hit count supporting it. subworkflows/mapping.nf takes it from there -- preferring Kraken2's calls where
// the two classifiers agree, and building consensus for the genuinely new ones.
//
// It maps NO reads. Species are called on read-hit counts alone (themisto_align_min_hits
// / metagraph_align_min_hits). The map-QC step that used to run here -- bowtie2 the reads
// against each called species' reference, then samtools coverage for breadth -- was
// removed: it meant a surviving species got mapped twice, once to measure breadth and
// again for consensus. Breadth is now measured once, downstream, from the consensus
// alignment MAPPING performs anyway, and MAPPING applies the breadth threshold there (see
// params.new_species_min_breadth_pct). THEMISTO_MAP_QC.nf / METAGRAPH_MAP_QC.nf are kept
// but unused.
//
// The consequence to be aware of: calls leaving here are hit-count-only, so they are
// less filtered than they used to be. Index noise that breadth would have rejected now
// reaches MAPPING and is rejected after its consensus alignment instead.
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

        // -- Themisto2 pseudoalignment, species called from read-hit counts.
        //
        // THE LANE'S DEFAULT METHOD (run_themisto defaults true). Species are called
        // directly from Themisto2 pseudoalignment read-hit counts (CALL_THEMISTO_SPECIES)
        // -- no probabilistic model and no validation mapping involved. mSWEEP's
        // abundance estimate is an optional add-on *inside* this arm, gated by run_msweep
        // (default false) inside VIRAL_THEMISTO_MSWEEP itself; see nextflow.config's note
        // on run_msweep's changed meaning.
        if (params.run_themisto) {
            VIRAL_THEMISTO_MSWEEP(preprocessed_3tuple_ch)

            themisto_counts_ch = VIRAL_THEMISTO_MSWEEP.out.species_hits
                .map { meta, tsv -> [meta.id, count_species_hits(tsv, 'themisto')] }


            // Both of these are Channel.empty() unless run_msweep is set (see
            // ../workflows/VIRAL_THEMISTO_MSWEEP.nf), so they need no gate of their own --
            // an empty channel simply contributes no counts and the joins below fill in
            // EMPTY_*_COUNTS.
            msweep_counts_ch = VIRAL_THEMISTO_MSWEEP.out.abundances
                .map { meta, abundances, _probs -> [meta.id, count_msweep_abundances(abundances)] }

            // Kept as its own variable so the species-calls block below can consume it
            // without touching VIRAL_THEMISTO_MSWEEP.out, which is undefined unless the
            // subworkflow was actually invoked. Optional per sample: unwritten when
            // nothing cleared min-hits.
            themisto_hits_ch     = VIRAL_THEMISTO_MSWEEP.out.species_hits
            themisto_labels_ch   = VIRAL_THEMISTO_MSWEEP.out.index_label_map
        } else {
            themisto_counts_ch   = Channel.empty()
            msweep_counts_ch     = Channel.empty()
            themisto_hits_ch     = Channel.empty()
            themisto_labels_ch   = Channel.empty()
        }

        // -- Sequence-to-graph alignment via Metagraph (metagraph align).
        if (params.run_metagraph_align) {
            VIRAL_METAGRAPH_ALIGN(preprocessed_3tuple_ch)

            metagraph_align_counts_ch = VIRAL_METAGRAPH_ALIGN.out.species_hits
                .map { meta, tsv -> [meta.id, count_species_hits(tsv, 'metagraph_align')] }

            metagraph_align_hits_ch   = VIRAL_METAGRAPH_ALIGN.out.species_hits
            metagraph_align_labels_ch = VIRAL_METAGRAPH_ALIGN.out.index_label_map
        } else {
            metagraph_align_counts_ch = Channel.empty()
            metagraph_align_hits_ch   = Channel.empty()
            metagraph_align_labels_ch = Channel.empty()
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

            metagraph_query_hits_ch   = VIRAL_METAGRAPH_QUERY.out.species_hits
            metagraph_query_labels_ch = VIRAL_METAGRAPH_QUERY.out.index_label_map
        } else {
            metagraph_query_counts_ch = Channel.empty()
            metagraph_query_hits_ch   = Channel.empty()
            metagraph_query_labels_ch = Channel.empty()
        }

        // -- Species calls handed to MAPPING (opt-in): every species a sequence-index
        // method called on read-hit count alone, plus the reference record the index
        // points at for it.
        //
        // Built from two files the species-calling step already writes -- species_hits.tsv
        // (species, hit_count, provisional_call) and index_label_map.tsv
        // (record_id -> species) -- joined per sample. No mapping and no breadth here:
        // hit count is the whole calling criterion now, and breadth is measured downstream
        // by MAPPING off the consensus alignment (see this file's header).
        //
        // Both files are per-sample and index_label_map is optional (unwritten when
        // nothing cleared min-hits), so join() -- 1:1 per sample -- naturally drops
        // samples with no calls, which is correct: there is nothing to hand over for them.
        if (params.call_consensus_for_new_species) {
            // These consume the *_hits_ch/*_labels_ch variables set in each method's
            // if/else above, NOT VIRAL_*.out directly: a subworkflow that was never invoked
            // has no .out at all, so reaching for it aborts the run with "Access to
            // 'VIRAL_METAGRAPH_ALIGN.out' is undefined" the moment this feature is enabled
            // with any subset of the three methods. An empty channel simply contributes
            // nothing, which is what "that method is off" should mean here.
            themisto_calls_ch = themisto_hits_ch
                .join(themisto_labels_ch)
                .flatMap { meta, hits, labels -> parse_species_calls(hits, labels, 'themisto').collect { call -> [meta.id, call] } }

            metagraph_align_calls_ch = metagraph_align_hits_ch
                .join(metagraph_align_labels_ch)
                .flatMap { meta, hits, labels -> parse_species_calls(hits, labels, 'metagraph_align').collect { call -> [meta.id, call] } }

            metagraph_query_calls_ch = metagraph_query_hits_ch
                .join(metagraph_query_labels_ch)
                .flatMap { meta, hits, labels -> parse_species_calls(hits, labels, 'metagraph_query').collect { call -> [meta.id, call] } }

            // One call per (sample, species). .unique() streams -- it emits each
            // non-duplicate immediately as it passes, it does not need to see the whole
            // channel close first (unlike groupTuple()).
            //
            // CAVEAT with more than one method enabled: which method's row wins here, and
            // therefore which reference record and hit count get attributed to a species
            // both methods called, is whichever arrives first -- task completion order, so
            // not reproducible run to run. Harmless on the defaults (only run_themisto is
            // on, so there is nothing to race), and the species itself is unaffected --
            // only the record and count reported for it. If a multi-method run ever needs
            // determinism here, rank by method instead of by arrival.
            //
            // No mSWEEP calls: mSWEEP estimates abundance only and calls no species (see
            // ../workflows/VIRAL_THEMISTO_MSWEEP.nf).
            species_calls_ch = themisto_calls_ch
                .mix(metagraph_align_calls_ch, metagraph_query_calls_ch)
                .unique { sample_id, call -> [sample_id, call.species_name.trim().toLowerCase()] }

            // Counts what this classifier can honestly measure: species it called on hit
            // count and reported to MAPPING, BEFORE MAPPING drops the ones Kraken2 already
            // found and before MAPPING's post-consensus breadth gate. For what actually
            // survived, count classification-report rows carrying
            // `discovered_by: 'sequence_index'`.
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
            .join(msweep_counts_ch, remainder: true)
            .join(metagraph_align_counts_ch, remainder: true)
            .join(metagraph_query_counts_ch, remainder: true)
            .join(new_species_counts_ch, remainder: true)
            // Parameter count matches the tuple width: five joins onto the backbone, so
            // six values after the key. The three *_mapqc_* slots that used to sit in here
            // went with map-QC -- breadth is no longer measured at this stage.
            .map { id, meta, t_counts, m_counts, mga_counts, mgq_counts, ns_counts ->
                def new_meta = meta + (t_counts ?: EMPTY_THEMISTO_COUNTS) +
                    (m_counts ?: EMPTY_MSWEEP_COUNTS) +
                    (mga_counts ?: EMPTY_METAGRAPH_ALIGN_COUNTS) +
                    (mgq_counts ?: EMPTY_METAGRAPH_QUERY_COUNTS) +
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
        // called on read-hit count, as [sample_id, call] where call is a Map of
        // species_name, reference_record, reference_source, hit_count and method.
        //
        // reference_source ('seqidx' or 'metagraph') says which reference FASTA
        // reference_record indexes into -- the two families of method report ids in
        // different namespaces (see parse_species_calls()), and MAPPING extracts from the
        // matching file. There is no breadth here: MAPPING measures it from its own
        // consensus alignment and applies params.new_species_min_breadth_pct there.
        //
        // Note this is NOT the shape classifying_kraken2.nf hands over -- that one emits
        // reads-plus-reference ready for consensus, because SORT_READS_BY_REF resolves its
        // own references upstream. The asymmetry is deliberate: it is what lets MAPPING
        // resolve references only for species that survive its filter. MAPPING builds the
        // consensus-ready shape for this side itself.
        //
        // Channel.empty() unless --call_consensus_for_new_species is set, so MAPPING can
        // consume it unconditionally.
        species_calls_ch // [sample_id, [species_name:, reference_record:, reference_source:, hit_count:, method:]]
}

// --- rvi_integration_1: sample-level count helpers for the mapping report ---
// Default counts for a sample a given optional step produced no output for -- either
// because that method didn't run at all (join remainder is null) or because the method
// ran but the step's own output is itself optional per-sample (e.g. nothing above a
// min-abundance/min-hits threshold). Named constants (not inline [:]) so every sample
// still gets the same report columns regardless of which method(s) actually ran for it.
EMPTY_THEMISTO_COUNTS = [themisto_n_species_considered: 0, themisto_n_species_called: 0]
EMPTY_MSWEEP_COUNTS = [msweep_n_groups: 0, msweep_top_group: '', msweep_top_abundance: 0.0]
// One per Metagraph method (align, query) -- both call count_species_hits() with a
// distinct prefix, since both methods' counts can merge into the same per-sample meta and
// would otherwise collide on field name.
EMPTY_METAGRAPH_ALIGN_COUNTS = [metagraph_align_n_species_considered: 0, metagraph_align_n_species_called: 0]
EMPTY_METAGRAPH_QUERY_COUNTS = [metagraph_query_n_species_considered: 0, metagraph_query_n_species_called: 0]
// New-species candidates (--call_consensus_for_new_species): how many species a
// sequence-index method called on read-hit count and has a reference record for, i.e. how
// many were handed to MAPPING as consensus candidates -- 0 whenever the feature is off, or
// on but nothing cleared min-hits for this sample. Whether MAPPING then kept them
// (Kraken2 hadn't already found them, and the consensus cleared
// new_species_min_breadth_pct) is not visible from here, by design: see the counts block
// above.
EMPTY_NEW_SPECIES_COUNTS = [new_species_candidates_n: 0]

def empty_species_hits_counts(prefix) {
    return prefix == 'metagraph_align' ? EMPTY_METAGRAPH_ALIGN_COUNTS : EMPTY_METAGRAPH_QUERY_COUNTS
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



def parse_species_calls(hits_tsv, label_map_tsv, method) {
    // The two files a method's species-calling step already writes, joined into one call
    // per species MAPPING should consider for consensus:
    //
    //   hits_tsv      <sample>_species_hits.tsv -- sample_id, species, hit_count,
    //                 provisional_call. provisional_call is `hit_count >= min_hits`
    //                 rendered by Python's str(bool), so "True"/"False", not lowercase.
    //   label_map_tsv <sample>_index_label_map.tsv -- headerless "<record_id>\t<species>",
    //                 written only for the species that cleared min_hits. This is the
    //                 "ideal reference" per call.
    //
    // Hit count is the entire calling criterion here: the map-QC step that used to measure
    // breadth for these species is gone, so there is no breadth_pct to threshold on yet.
    // MAPPING applies params.new_species_min_breadth_pct after its consensus alignment
    // instead -- see this file's header and subworkflows/mapping.nf.
    //
    // Species are matched between the two files on the name as written, which is the same
    // string in both (both come from one display_name()/species key in the same Python
    // loop), normalized only for whitespace/case so a trailing-space difference cannot
    // silently drop a call.
    //
    // reference_record is only meaningful against the reference FASTA the calling method's
    // own index was built from, and the two families of method disagree about that:
    // Themisto2 reports positional SEQIDX_<n> ids into params.msweep_map_reference_fasta,
    // Metagraph reports a bare taxid or an accession into
    // params.metagraph_map_reference_fasta. reference_source carries that distinction to
    // MAPPING, which needs it to extract the record from the right file.
    if (hits_tsv == null || !hits_tsv.exists()) return []
    if (label_map_tsv == null || !label_map_tsv.exists()) return []

    def record_by_species = [:]
    label_map_tsv.readLines().each { String line ->
        def trimmed = line.trim()
        if (!trimmed) return
        def cols = trimmed.split('\t')
        if (cols.size() < 2) return
        // First writer of a species wins, matching how the calling scripts de-duplicate.
        def key = cols[1].trim().toLowerCase()
        if (!record_by_species.containsKey(key)) record_by_species[key] = cols[0].trim()
    }

    def lines = hits_tsv.readLines()
    if (lines.size() < 2) return []
    def header = lines[0].split('\t')
    def species_idx = header.findIndexOf { String col -> col == 'species' }
    def hits_idx    = header.findIndexOf { String col -> col == 'hit_count' }
    def called_idx  = header.findIndexOf { String col -> col == 'provisional_call' }
    if (species_idx < 0 || hits_idx < 0 || called_idx < 0) {
        error("species-hits table ${hits_tsv} is missing one of the species/hit_count/" +
              "provisional_call columns (header: ${header}). bin/call_themisto_species.py " +
              "or bin/call_metagraph_species.py changed its output -- update " +
              "parse_species_calls() in subworkflows/classifying_index.nf.")
    }

    def source = method == 'themisto' ? 'seqidx' : 'metagraph'
    def max_idx = [species_idx, hits_idx, called_idx].max()
    def calls = []
    lines[1..-1].each { String line ->
        def cols = line.split('\t')
        if (max_idx >= cols.size()) return
        if (cols[called_idx].trim() != 'True') return
        def species_name = cols[species_idx].trim()
        def record = record_by_species[species_name.toLowerCase()]
        if (record == null) {
            // The label map is written from the same called_species list, so a called
            // species with no record means the two files disagree -- except for Metagraph,
            // which legitimately drops a species whose best-hit record another species
            // already claimed (see call_metagraph_species.py's seen_record_ids). Skipped
            // either way: MAPPING has no reference to map it against.
            log.warn("${method}: called species '${species_name}' has no reference record in " +
                     "${label_map_tsv.name}; not offered to MAPPING for consensus")
            return
        }
        calls << [
            species_name:     species_name,
            reference_record: record,
            reference_source: source,
            hit_count:        cols[hits_idx].trim(),
            method:           method,
        ]
    }
    return calls
}
