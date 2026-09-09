// --- consensus generation, lineage calling and classification reporting --------
// The shared half of what mapping.nf used to do end to end: everything from
// GENERATE_CONSENSUS onward. The Kraken2 half moved to
// subworkflows/classifying_kraken2.nf, and this is now driven by EITHER classifier:
//
//   subworkflows/classifying_kraken2.nf  Kraken2 + Kraken2Ref taxid selection
//   subworkflows/classifying_index.nf    Themisto2/Metagraph species calls Kraken2 missed
//
// One consensus/Nextclade/subtyping/report pass runs over the union of both, rather than
// each classifier growing its own parallel copy of it. That union is why a species only
// Themisto2/Metagraph found also gets Nextclade, SARS-CoV-2 subtyping and a row in the
// classification report -- it once had its consensus published on its own with none of
// that.
//
// The two classifiers hand over DIFFERENT shapes, deliberately:
//   - CLASSIFYING_KRAKEN2 arrives consensus-ready (reads + reference), because
//     SORT_READS_BY_REF resolves its references as part of classifying;
//   - CLASSIFYING_INDEX arrives as species calls plus the reference record its own index
//     points at for each, and no reads. Extracting that reference and pairing reads for
//     consensus is this subworkflow's job, done only for the calls that survive the
//     "Kraken2 already found it" filter.
//
// This is the ONLY place either classifier's reads get mapped. CLASSIFYING_INDEX calls
// species on read-hit counts alone and maps nothing, so there is no breadth figure
// attached to its calls when they arrive -- which is why the breadth threshold that
// decides whether a sequence-index-only species is worth reporting is applied here, AFTER
// the consensus alignment, rather than up in the classifier (see
// params.new_species_min_breadth_pct below). The cost of that ordering is that a noise
// call's consensus is computed and then discarded; the saving is that a real call is
// mapped once instead of twice.
include {INDEX_REFERENCE_FASTA; EXTRACT_REFERENCE_RECORD} from '../modules/reference_subset.nf'
include {EXTRACT_METAGRAPH_REFERENCE_RECORD} from '../modules/metagraph_reference_subset.nf'
include {GENERATE_CONSENSUS} from '../workflows/GENERATE_CONSENSUS.nf'
include {SCOV2_SUBTYPING} from '../workflows/SCOV2_SUBTYPING.nf'
include {GENERATE_CLASSIFICATION_REPORT} from '../workflows/GENERATE_CLASSIFICATION_REPORT.nf'
include {RUN_NEXTCLADE} from '../workflows/RUN_NEXTCLADE.nf'
include {publish_consensus_files as publish_aln_files; publish_consensus_files as publish_nc_files; publish_consensus_files as publish_per_sample_json} from '../modules/publish_lite.nf'
include {publish_run_files} from '../modules/publish_lite.nf'

workflow MAPPING {
    /*
    -----------------------------------------------------------------
    Generates a consensus sequence per (sample, reference) pair
    (GENERATE_CONSENSUS), optionally runs Nextclade and SARS-CoV-2
    subtyping, and writes the final per-sample classification report.

    Also decides which species are worth a consensus at all: where
    both classifiers found the same species for a sample, Kraken2's
    call and reference win, and only the species Kraken2 missed are
    resolved and mapped off the sequence indexes' calls.
    -----------------------------------------------------------------
    # Inputs

    From CLASSIFYING_KRAKEN2, consensus-ready:

    - **kraken2_sample_taxid_ch**: tuple(meta, [read_1, read_2],
      reference_fasta). meta carries `id` ("<sample_id>.<taxid>"),
      `sample_id`, `taxid` and `reference_header`.
    - **kraken2_report_ch**: [join_key, report_meta], join_key == the
      matching consensus's meta.id. report_meta holds the per-(sample,
      reference) descriptive fields the classification report writes
      out (sample_id, virus_name, ref_selected, ...).
    - **identified_species_ch**: [sample_id, [species, ...]],
      normalized -- what Kraken2 already found, i.e. what the
      sequence-index side gets filtered against.

    From CLASSIFYING_INDEX, calls only (Channel.empty() when it didn't
    run):

    - **index_species_calls_ch**: [sample_id, call], call being a Map
      of species_name, reference_record, reference_source, hit_count
      and method. No breadth: nothing has mapped these reads yet.

    Plus **reads_ch**, tuple(meta, read_1, read_2) per sample, needed
    to map the index-side species this subworkflow resolves itself.
    -----------------------------------------------------------------
    */

    take:
        kraken2_sample_taxid_ch  // tuple (meta, [read_1, read_2], reference_fasta) -- CLASSIFYING_KRAKEN2
        kraken2_report_ch        // [join_key, report_meta]                         -- CLASSIFYING_KRAKEN2
        index_species_calls_ch   // [sample_id, call]                               -- CLASSIFYING_INDEX, or Channel.empty()
        identified_species_ch    // [sample_id, [normalized_species_name, ...]]     -- CLASSIFYING_KRAKEN2
        reads_ch                 // tuple (meta, read_1, read_2) -- preprocessed, one per sample

    main:
        // --- Kraken2's calls win; sequence-index species are mapped only if new --------
        // CLASSIFYING_INDEX reports species and the ideal reference record for each, but
        // generates no consensus, so both halves of "prefer Kraken2, map what it missed"
        // are decided here, where the consensus actually gets spent:
        //
        //   1. drop index calls for species Kraken2 already found for that sample --
        //      Kraken2's own reference selection is kept in preference to the index's;
        //   2. resolve a reference and pair reads for whatever survives, so it reaches
        //      GENERATE_CONSENSUS in exactly the shape the Kraken2 side arrives in.
        //
        // Resolving references only after step 1 is the point of doing it here: the
        // previous arrangement resolved them up in the classifier, i.e. also for species
        // that were about to be discarded as already-known.
        //
        // Complete the "already identified" side over the samples that actually have index
        // calls: exactly one entry each, empty list where CLASSIFYING_KRAKEN2 produced no
        // pre-report for that sample (its emit only covers samples with a non-empty one,
        // and a sample Kraken2 found nothing in is exactly the interesting case here).
        // remainder:true fills those with null. The third shape it yields -- a sample
        // Kraken2 found species in but no index method called anything for,
        // [key, null, list] -- is dropped: there is nothing to filter.
        identified_by_sample_ch = index_species_calls_ch
            .map { sample_id, _call -> [sample_id, sample_id] }
            .unique()
            .join(identified_species_ch, remainder: true)
            .filter { _sample_id, has_calls, _identified -> has_calls != null }
            .map { sample_id, _has_calls, identified -> [sample_id, identified ?: []] }

        // combine(by: 0), NOT join: many calls per sample against one identified-species
        // list, and join() pairs keys one-to-one instead of broadcasting -- the bug fixed
        // in 4973f35, which this must not reintroduce. Note the trailing throwaway
        // parameter on the .map: combine() leaves the joined element on the tuple, and a
        // closure's parameter count has to match the tuple's width.
        index_new_calls_ch = index_species_calls_ch
            .combine(identified_by_sample_ch, by: 0)
            .filter { _sample_id, call, identified ->
                !identified.contains(call.species_name.trim().toLowerCase())
            }
            .map { sample_id, call, _identified -> [sample_id, call] }

        // Build the synthetic per-(sample, species) meta. Three things to know:
        //
        // `taxid`/`selected_taxid` are a filesystem-safe slug of the species name, NOT a
        // real Kraken taxid -- there isn't one, these species were never Kraken2-sorted.
        // They still have to be *something*, since `id` (and therefore every publish path
        // and report join key) is built from them, exactly as "<sample_id>.<taxid>" is on
        // the Kraken2 side.
        //
        // The descriptive fields mirror bin/k2r_report.py's pre-report columns (sample_id,
        // virus, virus_name, selected_taxid, ref_selected, sample_subtype, flu_segment,
        // virus_subtype, parent_selected, num_reads, report_name), because this map is
        // carried through as the base of the Nextclade/report meta for either classifier.
        // Fields the sequence-index side has no honest equivalent for are left empty
        // rather than faked: `virus` (Kraken2's source species taxid) and
        // `flu_segment`/`virus_subtype`/`sample_subtype` (k2r_report.py's
        // influenza-specific parsing, which never ran for these). `num_reads` is empty
        // too: Kraken2's per-taxon read count is not the same measurement as a
        // pseudoalignment hit count, so the index's own figures go in their own fields.
        //
        // `ref_selected` matters most: SARS-CoV-2 subtyping branches on it below, so a
        // SARS-CoV-2 infection Kraken2 missed but a sequence index caught gets subtyped.
        //
        // discovered_by / discovered_by_method / index_* record that this species came
        // from a sequence index and what supported it. They ride along into the
        // classification report (which dumps meta per consensus), which is therefore where
        // to see what the sequence indexes actually contributed.
        index_new_meta_ch = index_new_calls_ch
            .map { sample_id, call ->
                def slug = call.species_name.replaceAll(/[^A-Za-z0-9]+/, '_').replaceAll(/^_+|_+$/, '')
                def meta = [
                    id: "${sample_id}.${slug}".toString(),
                    sample_id: sample_id,
                    taxid: slug,
                    selected_taxid: slug,
                    species_name: call.species_name,
                    virus: '',
                    virus_name: call.species_name,
                    ref_selected: call.species_name,
                    report_name: call.species_name,
                    sample_subtype: '',
                    flu_segment: '',
                    virus_subtype: '',
                    parent_selected: false,
                    num_reads: '',
                    // The Kraken2 side gets this from get_taxid_reference_files; here the
                    // reference is the single record the calling method validated the
                    // species against, so the species name is the honest label.
                    reference_header: call.species_name,
                    discovered_by: 'sequence_index',
                    discovered_by_method: call.method,
                    index_reference_record: call.reference_record,
                    index_reference_source: call.reference_source,
                    index_hit_count: call.hit_count,
                ]
                [meta, call]
            }

        // A record id only means something against the reference FASTA the calling
        // method's own index was built from, and the methods disagree about which that is:
        //
        //   'seqidx'    Themisto2 -- positional SEQIDX_<n> into msweep_map_reference_fasta
        //   'metagraph' Metagraph -- a taxid or accession into metagraph_map_reference_fasta
        //
        // so the calls are split and extracted separately, then re-merged. Sending them all
        // through one extractor would quietly produce no reference for the other family's
        // ids (an unmatched grep, an absent optional output, a species silently gone),
        // even though both FASTAs happen to be builds of the same RVDB release today.
        index_new_meta_ch
            .branch { _meta, call ->
                seqidx_ch:    call.reference_source == 'seqidx'
                metagraph_ch: call.reference_source == 'metagraph'
                unknown_ch:   true
            }
            .set { index_new_by_source_ch }

        // branch drops anything no arm matched, so the third arm exists purely to turn a
        // new/typo'd reference_source into a loud failure instead of a species that
        // vanishes between the classifier and the report.
        index_new_by_source_ch.unknown_ch
            .map { meta, call ->
                error("${meta.id}: sequence-index call for '${call.species_name}' has " +
                      "reference_source '${call.reference_source}', which MAPPING has no " +
                      "reference FASTA for. Add an arm to the branch above (and an " +
                      "extractor) or fix parse_species_calls() in " +
                      "subworkflows/classifying_index.nf.")
            }

        // Gated even though CLASSIFYING_INDEX already emits nothing when the feature is
        // off: without this, INDEX_REFERENCE_FASTA would still run its (mem_16) pass over
        // the reference FASTA on every default run, and every run would then require
        // msweep_map_reference_fasta to exist.
        if (params.call_consensus_for_new_species) {
            // The same reference FASTA the Themisto2 index was built from, re-tagged here
            // with the same positional SEQIDX_<n> ids -- a deterministic pass over the same
            // file, which is what makes the record ids the classifier reported valid to
            // grep for. Costs one extra pass over that FASTA per run.
            INDEX_REFERENCE_FASTA(Channel.fromPath(params.msweep_map_reference_fasta))

            EXTRACT_REFERENCE_RECORD(
                index_new_by_source_ch.seqidx_ch.map { meta, call -> [meta, call.reference_record] },
                INDEX_REFERENCE_FASTA.out.fasta.first()
            )

            // metagraph_record_pattern() builds the grep pattern in Groovy rather than in
            // the process's shell -- see modules/metagraph_reference_subset.nf.
            EXTRACT_METAGRAPH_REFERENCE_RECORD(
                index_new_by_source_ch.metagraph_ch.map { meta, call ->
                    [meta, metagraph_record_pattern(call.reference_record)]
                },
                Channel.fromPath(params.metagraph_map_reference_fasta).first()
            )

            // Deliberately no cross-sample de-duplication: a species found new in many
            // samples gets its reference extracted once per sample rather than once per
            // run -- keeps everything scoped per-sample (no batch-wide wait) at the cost of
            // some redundant extraction.
            reads_by_sample_ch = reads_ch
                .map { meta, r1, r2 -> [meta.id, r1, r2] } // meta.id == sample_id pre-classifier

            EXTRACT_REFERENCE_RECORD.out.subset_fasta
                .mix(EXTRACT_METAGRAPH_REFERENCE_RECORD.out.subset_fasta)
                .map { meta, ref_fa -> [meta.sample_id, meta, ref_fa] }
                .combine(reads_by_sample_ch, by: 0)
                .map { _sample_id, meta, ref_fa, r1, r2 -> [meta, [r1, r2], ref_fa] }
                .set { index_new_sample_taxid_ch }

            // Keyed the same way CLASSIFYING_KRAKEN2 keys its report rows: by the matching
            // consensus's meta.id.
            index_new_report_ch = index_new_sample_taxid_ch
                .map { meta, _reads, _ref_fa -> [meta.id, meta] }
        } else {
            index_new_sample_taxid_ch = Channel.empty()
            index_new_report_ch = Channel.empty()
        }

        // The union both classifiers feed. mix() (not join/combine) because post-filter
        // these are disjoint sets of (sample, reference) pairs rather than two views of the
        // same pair -- disjoint because of the filter above, not by assumption.
        sample_taxid_ch = kraken2_sample_taxid_ch.mix(index_new_sample_taxid_ch)
        sample_report_with_join_key_ch = kraken2_report_ch.mix(index_new_report_ch)

        GENERATE_CONSENSUS(sample_taxid_ch)

        // --- breadth gate on sequence-index-only species -------------------------------
        // The consensus alignment is the first (and now only) time these species' reads
        // are mapped, so this is the earliest point their breadth of coverage can be
        // judged -- hence a filter after GENERATE_CONSENSUS rather than a threshold in the
        // classifier. Below params.new_species_min_breadth_pct the species is dropped
        // completely: no consensus published, no Nextclade, no subtyping, no report row.
        // The wasted work is the consensus itself, which is the price of not mapping every
        // real call twice.
        //
        // WHICH breadth: percent_non_n_bases, the share of the consensus that is a real
        // base rather than N. `samtools mpileup -aa` (modules/run_ivar.nf) emits every
        // reference position including zero-coverage ones, and `ivar consensus -n N` pads
        // those with N, so the consensus is exactly reference-length and this figure is
        // genome breadth at iVar's minimum depth (params.ivar_polish_min_depth, 10x) over
        // the whole reference. It is deliberately NOT depth>=1 breadth, which the QC JSON
        // does not carry (bin/qc.py buckets depth in steps of 5 from 0, so "positions with
        // any coverage" is not among them) and which would be the wrong gate anyway: the
        // index noise this is meant to reject sits at 3-14% breadth and 0.09-0.46x mean
        // depth, so at 10% a depth>=1 threshold would admit the top of that range while a
        // 10x one cannot.
        //
        // Kraken2-side consensuses pass through untouched -- they are gated by Kraken2's
        // own read-count selection upstream, and retro-fitting a breadth threshold onto
        // them would silently change what the pipeline has always reported.
        GENERATE_CONSENSUS.out.filtered_consensus_ch
            .map { meta, bam, bam_idx, consensus, qc_json ->
                def json_map = new groovy.json.JsonSlurper().parse(new File(qc_json.toString()))
                def breadth = (json_map['percent_non_n_bases'] ?: 0) as Double
                // Recorded on meta for every consensus, both classifiers', so the report
                // shows the number the gate was applied to (or would have been).
                [meta + [consensus_breadth_pct: breadth], bam, bam_idx, consensus, qc_json]
            }
            .filter { meta, _bam, _bam_idx, _consensus, _qc_json ->
                if (meta.discovered_by != 'sequence_index') return true
                if (meta.consensus_breadth_pct >= params.new_species_min_breadth_pct) return true
                log.info("Dropping sequence-index-only species '${meta.species_name}' for " +
                         "${meta.sample_id}: consensus breadth ${meta.consensus_breadth_pct}% " +
                         "< new_species_min_breadth_pct (${params.new_species_min_breadth_pct}%)")
                return false
            }
            .set { consensus_ch }

        consensus_ch
            .map { meta, _bam, _bam_idx, consensus, _qc_json -> [meta.id, meta, consensus] }
            .set { consensus_fa_ch }

        // Nextclade's input meta is the report row (descriptive fields) widened with the
        // two fields only the consensus side knows: which reference record was used, and
        // the taxid the consensus was actually built against.
        sample_report_with_join_key_ch
            .combine(consensus_fa_ch, by: 0)
            .map { _id, report_meta, fa_meta, fa ->
                def final_meta = report_meta + [reference_header: "${fa_meta.reference_header}", taxid: "${fa_meta.taxid}"]
                [final_meta, fa]
            }
            .set { nextclade_In_ch }

        // TODO add check parameters
        if (params.nextclade_index_json == null) {
            log.warn("No nextclade_index_json provided, skipping nextclade analysis step")
            publish_nextclade_outputs_ch = channel.empty()
            per_consensus_nextclade_json_ch = channel.empty()
        } else {
            RUN_NEXTCLADE(nextclade_In_ch)
            RUN_NEXTCLADE.out
                .map { meta, _agg_json, tar_gz -> [meta, tar_gz] }
                .set { publish_nextclade_outputs_ch }

            RUN_NEXTCLADE.out
                .map { meta, json, _tarball -> [meta.id, json] }
                .set { per_consensus_nextclade_json_ch }
        }

        // add report info to out qc metric channel and branch for SCOV2 subtyping
        consensus_ch
            .map { meta, _bam, _bam_idx, consensus, _qc -> [meta.id, meta, consensus] }
            .join(sample_report_with_join_key_ch)
            .map { _id, meta, fasta, report ->
                def new_meta = meta.plus(report)
                [new_meta, fasta]
            }
            .branch { it ->
                scv2_subtyping_workflow_in_ch: it[0].ref_selected.contains("${params.scv2_keyword}")
                no_subtyping_ch: true
            }
            .set { filtered_consensus_by_type_ch }

        if (params.do_scov2_subtyping == true) {
            SCOV2_SUBTYPING(filtered_consensus_by_type_ch.scv2_subtyping_workflow_in_ch)
            scov2_subtyped_ch = SCOV2_SUBTYPING.out
        } else {
            scov2_subtyped_ch = channel.empty()
        }

        // write final classification reports
        filtered_consensus_by_type_ch.no_subtyping_ch.concat(scov2_subtyped_ch)
            .map { meta, _fasta -> [meta.id, meta] }
            .set { report_in_ch }

        consensus_ch
            .map { meta, _bam, _bam_idx, _consensus, qc_json -> [meta.id, qc_json] }
            .set { qc_json_simplified_ch }

        consensus_ch
            .map { meta, bam, bam_idx, consensus, _qc_json -> [meta, [bam, bam_idx, consensus]] }
            .set { aln_publish_ch }

        report_in_ch // [meta.id, meta]
            .join(qc_json_simplified_ch, remainder: true) // [meta.id, meta, qc_json]
            .join(per_consensus_nextclade_json_ch, remainder: true) // [meta.id, meta, qc_json, nc_json]
            .set { report_prep_ch }

        GENERATE_CLASSIFICATION_REPORT(report_prep_ch)

        // PUBLISH (mapping lane)
        publish_aln_files(aln_publish_ch)
        publish_nc_files(publish_nextclade_outputs_ch)
        publish_per_sample_json(GENERATE_CLASSIFICATION_REPORT.out.publish_seq_level_ch)
        publish_run_files(GENERATE_CLASSIFICATION_REPORT.out.publish_run_level_summaries_ch)
}

// The Groovy half of a rule that also exists in Python: bin/call_metagraph_species.py's
// build_record_id_pattern(). Metagraph reports a species' reference either as a bare taxid
// (for 'kraken:taxid|<taxid>|...'-shaped index labels) or as a complete accession, and the
// two need different anchoring to grep out of metagraph_map_reference_fasta:
//
//   taxid     anchored on the fixed 'kraken:taxid|<taxid>|' prefix every such header
//             shares, so taxid 13000336 cannot match a header for 130003360;
//   accession anchored to the start of the header and required to be followed by
//             whitespace or end-of-line, so it cannot match a longer accession that has it
//             as a prefix.
//
// Kept here rather than in the process's shell because the pattern would otherwise have to
// survive Nextflow's string interpolation and then the shell's quoting on its way to
// seqkit; the Python side is the same rule for the file-of-patterns path that
// EXTRACT_METAGRAPH_REFERENCE_SUBSET still uses. Change one, change the other.
def metagraph_record_pattern(record_id) {
    def id = record_id.toString().trim()
    if (id ==~ /^[0-9]+$/) {
        return "^kraken:taxid\\|${id}\\|".toString()
    }
    // Metacharacters escaped one by one rather than wrapped in a \Q...\E literal span:
    // seqkit's regex engine is Go's, which does accept \Q...\E, but this cannot be
    // exercised without seqkit installed, and an accession's '.' silently acting as a
    // wildcard (matching a near-identical accession) is a worse failure than none. The
    // class below is the set Python's re.escape would touch in this input; '.' is the only
    // one accessions actually contain.
    def escaped = id.replaceAll(/([.^$*+?()\[\]{}|\\])/, '\\\\$1')
    return "^${escaped}(\\s|\$)".toString()
}
