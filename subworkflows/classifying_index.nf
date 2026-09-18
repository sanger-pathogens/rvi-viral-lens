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
// It maps NO reads. Breadth is measured downstream, from the consensus alignment MAPPING
// performs anyway, and MAPPING applies the breadth threshold there (see
// params.new_species_min_breadth_pct).
//
// THREE GATES decide a call, all of them inside the species callers and all folded into
// the one `provisional_call` column, so everything here and downstream respects them
// without knowing they exist (see rvi_toolbox/modules/themisto_species_call.nf):
//
//   READ HITS   themisto_align_min_hits / metagraph_align_min_hits.
//   TAXONOMY    the species' lineage must sit under params.taxon_filter_whitelist (by
//               default the eight respiratory virus families) and not under
//               taxon_filter_blacklist, resolved through params.taxon_filter_table.
//   REFERENCE   the record the call resolved to must be at least
//               params.min_called_reference_length bases.
//
// The last two exist because these methods query a WHOLE-VIROME index, unlike Kraken2's
// curated database: read hits there answer "is this sequence in the index and did reads
// match it", not "is this a virus we report, with enough genome behind it to be worth a
// consensus". Without them, most of what cleared min-hits on a respiratory sample was
// neither -- phage, plant and insect viruses sharing k-mers, and partial-CDS records that
// pass a percentage-of-reference breadth gate trivially precisely because they are short.
//
// Breadth is still measured only downstream, so the remaining noise these gates do not
// catch (a whitelisted family's species with a full-length reference and no real coverage)
// still reaches MAPPING and is rejected there, after its consensus alignment.
include {VIRAL_THEMISTO} from '../workflows/VIRAL_THEMISTO.nf'
include {VIRAL_METAGRAPH_ALIGN} from '../workflows/VIRAL_METAGRAPH_ALIGN.nf'
include {VIRAL_METAGRAPH_QUERY} from '../workflows/VIRAL_METAGRAPH_QUERY.nf'
include {GENERATE_MAPPING_REPORT} from '../workflows/GENERATE_MAPPING_REPORT.nf'
include {publish_lane_json as publish_mapping_lane_json} from '../modules/publish_lane_report.nf'
include {publish_run_files as publish_mapping_run_files} from '../modules/publish_lite.nf'

workflow CLASSIFYING_INDEX {
    take:
        preprocessed_3tuple_ch  // tuple (meta, read1, read2)
        identified_species_ch   // [sample_id, [normalized_species_name, ...]] -- CLASSIFYING_KRAKEN2

    main:
        sequence_index_sample_ch = preprocessed_3tuple_ch
            .map { meta, _r1, _r2 -> [meta.id, meta] }

        // -- Themisto2 pseudoalignment, species called from read-hit counts.
        //
        // THE LANE'S DEFAULT METHOD (run_themisto defaults true). Species are called
        // directly from Themisto2 pseudoalignment read-hit counts (CALL_THEMISTO_SPECIES)
        // -- no probabilistic model and no validation mapping involved. mSWEEP is NOT part
        // of this lane: it estimates abundance rather than calling species, so it lives
        // behind the abundance lane's --run_msweep (subworkflows/abundance.nf), consuming
        // the pseudoalignments VIRAL_THEMISTO emits.
        if (params.run_themisto) {
            VIRAL_THEMISTO(preprocessed_3tuple_ch)

            themisto_counts_ch = VIRAL_THEMISTO.out.species_hits
                .map { meta, tsv -> [meta.id, count_species_hits(tsv, 'themisto')] }


            // Kept as its own variable so the species-calls block below can consume it
            // without touching VIRAL_THEMISTO.out, which is undefined unless the
            // subworkflow was actually invoked. Optional per sample: unwritten when
            // nothing cleared min-hits.
            themisto_hits_ch     = VIRAL_THEMISTO.out.species_hits
            themisto_labels_ch   = VIRAL_THEMISTO.out.index_label_map
            // Handover for the abundance lane's optional MSWEEP (--run_msweep). Empty
            // whenever run_themisto is off, which is exactly what "no Themisto2 result
            // available" has to mean there.
            themisto_pseudoaln_ch = VIRAL_THEMISTO.out.pseudoalignments
            themisto_ref_groups_ch = VIRAL_THEMISTO.out.ref_groups
        } else {
            themisto_counts_ch   = Channel.empty()
            themisto_hits_ch     = Channel.empty()
            themisto_pseudoaln_ch = Channel.empty()
            themisto_ref_groups_ch = Channel.empty()
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
        // ../workflows/VIRAL_METAGRAPH_QUERY.nf / ../rvi_toolbox/modules/metagraph_query.nf for why
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
        // method called, plus the reference record the index points at for it.
        //
        // Built from two files the species-calling step already writes -- species_hits.tsv
        // (species, hit_count, provisional_call, + the gate columns) and index_label_map.tsv
        // (record_id -> species) -- joined per sample. No mapping and no breadth here: the
        // callers' three gates have already decided provisional_call, and breadth is
        // measured downstream by MAPPING off the consensus alignment (see this file's
        // header). Nothing in here re-checks the gates, and nothing needs to: an excluded
        // species simply never carries provisional_call True.
        //
        // Both files are per-sample and index_label_map is optional (unwritten when
        // nothing cleared min-hits), so join() -- 1:1 per sample -- naturally drops
        // samples with no calls, which is correct: there is nothing to hand over for them.
        // -- Which species each method CALLED, independent of every downstream gate -----
        // Deliberately NOT derived from species_calls_ch below: that lives behind
        // --call_consensus_for_new_species and is empty by default, which would make both
        // overlapping_n_species and the report's Discovered_By silently wrong rather than
        // absent. This reads the same species-hits tables the counts come from, so it is
        // available whenever the method ran at all.
        //
        // [sample_id, species_lower, method] -- one entry per (sample, species, method).
        called_species_ch = themisto_hits_ch
            .flatMap { meta, tsv -> called_species_names(tsv).collect { sp -> [meta.id, sp, 'themisto2'] } }
            .mix(
                metagraph_align_hits_ch
                    .flatMap { meta, tsv -> called_species_names(tsv).collect { sp -> [meta.id, sp, 'metagraph_align'] } },
                metagraph_query_hits_ch
                    .flatMap { meta, tsv -> called_species_names(tsv).collect { sp -> [meta.id, sp, 'metagraph_query'] } }
            )

        // How many of the species Kraken2 SELECTED for this sample were also called by a
        // sequence-index method. Kraken2's side is identified_species_ch -- the species it
        // acted on (k2r pre-report's virus_name + ref_selected, already normalized), not
        // every row of the raw Kraken2 report. Counted per sample over the distinct
        // species names, so a species two methods both called counts once.
        //
        // Skipped entirely with --do_mapping false, so the column fills as NA rather than
        // 0. identified_species_ch is empty then (Kraken2 never ran), and the
        // remainder: true join below would happily pair every sample against a null
        // Kraken2 side and count an overlap of 0 -- which reads as "nothing Kraken2 found
        // was corroborated" when the truth is that Kraken2 was never asked. Same
        // NA-vs-0 distinction the count helpers below exist for.
        if (params.do_mapping) {
            overlap_counts_ch = called_species_ch
                .map { sample_id, species, _method -> [sample_id, species] }
                .unique()
                .groupTuple()
                .join(identified_species_ch, remainder: true)
                .filter { _sample_id, called, _identified -> called != null }
                .map { sample_id, called, identified ->
                    def kraken2_set = (identified ?: []) as Set
                    [sample_id, [overlapping_n_species: called.count { sp -> kraken2_set.contains(sp) }]]
                }
        } else {
            overlap_counts_ch = Channel.empty()
        }

        if (params.call_consensus_for_new_species) {
            // These consume the *_hits_ch/*_labels_ch variables set in each method's
            // if/else above, NOT VIRAL_*.out directly: a subworkflow that was never invoked
            // has no .out at all, so reaching for it aborts the run with "Access to
            // 'VIRAL_METAGRAPH_ALIGN.out' is undefined" the moment this feature is enabled
            // with any subset of the three methods. An empty channel simply contributes
            // nothing, which is what "that method is off" should mean here.
            themisto_calls_ch = themisto_hits_ch
                .join(themisto_labels_ch)
                .flatMap { meta, hits, labels -> parse_species_calls(hits, labels, 'themisto2').collect { call -> [meta.id, call] } }

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
            // ../workflows/VIRAL_THEMISTO.nf).
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
            .join(metagraph_align_counts_ch, remainder: true)
            .join(metagraph_query_counts_ch, remainder: true)
            .join(new_species_counts_ch, remainder: true)
            .join(overlap_counts_ch, remainder: true)
            // Parameter count matches the tuple width: four joins onto the backbone, so
            // five values after the key.
            .map { id, meta, t_counts, mga_counts, mgq_counts, ns_counts, ov_counts ->
                def new_meta = meta + (t_counts ?: empty_species_hits_counts('themisto')) +
                    (mga_counts ?: empty_species_hits_counts('metagraph_align')) +
                    (mgq_counts ?: empty_species_hits_counts('metagraph_query')) +
                    (ns_counts ?: empty_new_species_counts()) +
                    (ov_counts ?: empty_overlap_counts())
                [id, new_meta]
            }
            .set { mapping_report_prep_ch }

        GENERATE_MAPPING_REPORT(mapping_report_prep_ch)

        // PUBLISH (mapping/sequence-index lane)
        publish_mapping_lane_json(GENERATE_MAPPING_REPORT.out.publish_seq_level_ch)
        publish_mapping_run_files(GENERATE_MAPPING_REPORT.out.publish_run_level_summaries_ch)

    emit:
        // Handover to subworkflows/mapping.nf: one entry per (sample, species) this lane
        // called -- read hits, taxonomy and reference length all cleared -- as
        // [sample_id, call] where call is a Map of species_name, reference_record,
        // reference_source, hit_count and method.
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

        // Every (sample, species, method) this lane called, with no gate on it -- unlike
        // species_calls_ch above, which is empty unless --call_consensus_for_new_species.
        // MAPPING uses it to record ALL the methods that found a species in the report's
        // Discovered_By, including for species Kraken2 also found and therefore won.
        called_species_ch // [sample_id, species_lower, method]

        // Themisto2's pseudoalignments and the species_labels.txt that goes with them,
        // for the abundance lane's optional MSWEEP (--run_msweep). Both Channel.empty()
        // when run_themisto is off; abundance.nf treats that as "Themisto2 has produced
        // nothing to estimate from" and refuses to run mSWEEP.
        themisto_pseudoalignments = themisto_pseudoaln_ch
        themisto_ref_groups       = themisto_ref_groups_ch
}

// --- rvi_integration_1: sample-level count helpers for the mapping report ---
// Default counts for a sample a given optional step produced no output for -- either
// because that method didn't run at all (join remainder is null) or because the method
// ran but the step's own output is itself optional per-sample (e.g. nothing above a
// min-abundance/min-hits threshold). Named constants (not inline [:]) so every sample
// still gets the same report columns regardless of which method(s) actually ran for it.
// A count of 0 and a method that never ran are different facts, and a report that spells
// both "0" cannot be read correctly -- that is exactly how new_species_candidates_n's 0
// was misread as "the index found nothing new" when the feature was simply off. So the
// fill value is NA when the corresponding flag is off, and 0 only when the step genuinely
// ran and found nothing.
NOT_RUN = 'NA'

def empty_overlap_counts() {
    // NA when the mapping lane did not run: there were no Kraken2 calls to overlap with,
    // which is not the same fact as an overlap of zero.
    return [overlapping_n_species: params.do_mapping ? 0 : NOT_RUN]
}

def empty_new_species_counts() {
    return [new_species_candidates_n: params.call_consensus_for_new_species ? 0 : NOT_RUN]
}

// New-species candidates (--call_consensus_for_new_species): how many species a
// sequence-index method called and has a reference record for, i.e. how many were handed
// to MAPPING as consensus candidates -- 0 whenever the feature is off, or on but nothing
// cleared the callers' gates for this sample. Whether MAPPING then kept them
// (Kraken2 hadn't already found them, and the consensus cleared
// new_species_min_breadth_pct) is not visible from here, by design: see the counts block
// above.
EMPTY_NEW_SPECIES_COUNTS = [new_species_candidates_n: 0]

def empty_species_hits_counts(prefix) {
    def ran = prefix == 'metagraph_align' ? params.run_metagraph_align
            : prefix == 'metagraph_query' ? params.run_metagraph_query
            : params.run_themisto
    def v = ran ? 0 : NOT_RUN
    // The gate counts get their own fill value, on the same NA-vs-0 reasoning as the
    // others but one level down: a method that never ran reports NA, and a method that ran
    // with a gate switched off reports NA for THAT gate specifically -- 0 there would read
    // as "the gate ran and rejected nothing", which is a different fact from "the gate was
    // never applied". Same trap as new_species_candidates_n's 0 (see NOT_RUN above).
    def taxon_v = (ran && params.run_taxon_filter) ? 0 : NOT_RUN
    def ref_v = (ran && params.min_called_reference_length > 0) ? 0 : NOT_RUN
    // Parentheses are required: a GString map key is a computed key, and Groovy parses a
    // bare one as a label instead.
    return [
        ("${prefix}_n_species_considered".toString()):      v,
        ("${prefix}_n_species_called".toString()):          v,
        ("${prefix}_n_species_taxon_filtered".toString()):  taxon_v,
        ("${prefix}_n_species_short_reference".toString()): ref_v,
    ]
}


def count_species_hits(tsv, prefix) {
    // <sample>_species_hits.tsv: sample_id, species, hit_count, provisional_call, then the
    // call-gate columns (taxon_filter, taxonomy_id, family, family_taxon_id,
    // reference_record, reference_length, reference_filter). provisional_call is written by
    // Python's str(bool), so "True"/"False", not lowercase. All three read-hit methods emit
    // this identical schema (rvi_toolbox/bin/call_metagraph_species.py for both
    // Metagraph methods, rvi_toolbox/bin/call_themisto_species.py for Themisto2), so one parser serves
    // them all -- prefix ('themisto', 'metagraph_align' or 'metagraph_query') keeps their
    // counts from colliding when several merge into the same per-sample meta.
    //
    // Columns are found BY NAME, never by position, so the gate columns could be appended
    // without touching this -- and so a future column added in the middle cannot silently
    // shift what gets counted.
    if (tsv == null || !tsv.exists()) return empty_species_hits_counts(prefix)
    def lines = tsv.readLines()
    if (lines.size() < 2) return empty_species_hits_counts(prefix)
    def header = lines[0].split('\t')
    def called_idx = header.findIndexOf { String col -> col == 'provisional_call' }
    def taxon_idx  = header.findIndexOf { String col -> col == 'taxon_filter' }
    def ref_idx    = header.findIndexOf { String col -> col == 'reference_filter' }
    def n_considered = lines.size() - 1
    def n_called = 0
    // Why these two are counted separately rather than lumped into one "filtered" figure:
    // they answer different questions about a run. A large taxon-filtered count is the gate
    // doing its job on a whole-virome index and says nothing is wrong. A large
    // short-reference count says the index's representative records for otherwise-wanted
    // families are fragments, which is a property of the reference set worth noticing. And
    // a taxon-filtered count of ~everything, with the callers' "not one label resolved"
    // warning in the log, means the taxonomy table does not match the index at all.
    def n_taxon_filtered = 0
    def n_short_reference = 0
    lines[1..-1].each { String line ->
        def cols = line.split('\t')
        if (called_idx >= 0 && called_idx < cols.size() && cols[called_idx] == 'True') n_called += 1
        // NB `+= 1` rather than `++` on these three counters is for tooling, not style:
        // `nextflow lint` (25.10.3) cannot parse a postfix `++` here -- it crashes on a
        // braceless `if x` body ("Range [70, 71) out of bounds") and rejects a braced one
        // inside a closure ("Unexpected input: '}'"). Groovy accepts every form; only the
        // linter does not, and `+= 1` is the one it reads cleanly.
        //
        // Any reason other than a pass or a disabled gate is a rejection by that gate.
        // Counted this way (rather than matching each reason string) so a new rejection
        // reason is included the day it is added, instead of quietly counting as zero.
        if (taxon_idx >= 0 && taxon_idx < cols.size()) {
            def reason = cols[taxon_idx].trim()
            if (reason && reason != 'pass' && reason != 'off') n_taxon_filtered += 1
        }
        if (ref_idx >= 0 && ref_idx < cols.size()) {
            def reason = cols[ref_idx].trim()
            // 'not_evaluated' is not a reference rejection: it means the taxonomy gate had
            // already rejected the species, so no reference was ever resolved for it.
            // Counting it here would double-count every taxon-filtered species.
            if (reason && !(reason in ['pass', 'off', 'not_evaluated'])) n_short_reference += 1
        }
    }
    // A gate that is switched off reports NA, not the 0 it counted -- and it has to be
    // decided from the params here, exactly as empty_species_hits_counts() decides it for a
    // sample with no hits file at all. Reading it off the file instead would make the two
    // paths disagree within one run: with a gate off, every row's reason is 'off' so
    // nothing counts as a rejection, and this would report 0 for samples that produced a
    // hits table while empty_species_hits_counts() reported NA for samples that did not.
    // Same 0-is-not-NA trap as new_species_candidates_n (see NOT_RUN above), one level in.
    return [
        ("${prefix}_n_species_considered".toString()):      n_considered,
        ("${prefix}_n_species_called".toString()):          n_called,
        ("${prefix}_n_species_taxon_filtered".toString()):  params.run_taxon_filter ? n_taxon_filtered : NOT_RUN,
        ("${prefix}_n_species_short_reference".toString()): params.min_called_reference_length > 0 ? n_short_reference : NOT_RUN,
    ]
}



def called_species_names(hits_tsv) {
    // Normalized names of the species a method CALLED (provisional_call == "True") in its
    // <sample>_species_hits.tsv. Same table count_species_hits() counts, same "True"/"False"
    // spelling from Python's str(bool). Lowercased/trimmed to match identified_species_ch,
    // which normalizes the Kraken2 side the same way.
    if (hits_tsv == null || !hits_tsv.exists()) return []
    def lines = hits_tsv.readLines()
    if (lines.size() < 2) return []
    def header = lines[0].split('\t')
    def species_idx = header.findIndexOf { String col -> col == 'species' }
    def called_idx  = header.findIndexOf { String col -> col == 'provisional_call' }
    if (species_idx < 0 || called_idx < 0) return []
    def max_idx = [species_idx, called_idx].max()
    return lines[1..-1].collect { String line -> line.split('\t') }
        .findAll { cols -> max_idx < cols.size() && cols[called_idx].trim() == 'True' }
        .collect { cols -> cols[species_idx].trim().toLowerCase() }
        .findAll { name -> name }
        .unique()
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
    // provisional_call is the callers' verdict across all three gates (read hits, taxonomy
    // whitelist/blacklist, reference length), so filtering on it here is all that is needed
    // -- the taxon_filter/reference_filter columns are for reading the table, not for
    // re-deciding. There is no breadth_pct at this stage: MAPPING applies
    // params.new_species_min_breadth_pct after its consensus alignment instead -- see this
    // file's header and subworkflows/mapping.nf.
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
              "provisional_call columns (header: ${header}). rvi_toolbox/bin/call_themisto_species.py " +
              "or rvi_toolbox/bin/call_metagraph_species.py changed its output -- update " +
              "parse_species_calls() in subworkflows/classifying_index.nf.")
    }

    def source = method == 'themisto2' ? 'seqidx' : 'metagraph'
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
