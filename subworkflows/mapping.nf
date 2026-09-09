// --- consensus generation, lineage calling and classification reporting --------
// The shared half of what mapping.nf used to do end to end: everything from
// GENERATE_CONSENSUS onward. The Kraken2 half moved to
// subworkflows/classifying_kraken2.nf, and this is now driven by EITHER classifier:
//
//   subworkflows/classifying_kraken2.nf  Kraken2 + Kraken2Ref taxid selection
//   subworkflows/classifying_index.nf    Themisto2/Metagraph species calls Kraken2 missed
//
// Both hand over the same two channel shapes (see either file's emit block), so this
// runs one consensus/Nextclade/subtyping/report pass over the union rather than each
// classifier growing its own parallel copy of it. That union is why a species only
// Themisto2/Metagraph found now also gets Nextclade, SARS-CoV-2 subtyping and a row in
// the classification report -- previously its consensus was published on its own with
// none of that.
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

    Reference selection itself happens upstream, in whichever
    classifier produced the input -- this subworkflow is deliberately
    agnostic about which one that was.
    -----------------------------------------------------------------
    # Inputs

    Two channels per classifier, both already in their final shape:

    - **sample_taxid_ch**: tuple(meta, [read_1, read_2], reference_fasta).
      meta carries `id` ("<sample_id>.<taxid>"), `sample_id`, `taxid` and
      `reference_header`.
    - **report_with_join_key_ch**: [join_key, report_meta], join_key ==
      the matching consensus's meta.id. report_meta holds the
      per-(sample, reference) descriptive fields the classification
      report writes out (sample_id, virus_name, ref_selected, ...).

    Pass Channel.empty() for a classifier that isn't running.
    -----------------------------------------------------------------
    */

    take:
        kraken2_sample_taxid_ch // tuple (meta, [read_1, read_2], reference_fasta) -- CLASSIFYING_KRAKEN2
        kraken2_report_ch       // [join_key, report_meta]                         -- CLASSIFYING_KRAKEN2
        index_sample_taxid_ch   // tuple (meta, [read_1, read_2], reference_fasta) -- CLASSIFYING_INDEX, or Channel.empty()
        index_report_ch         // [join_key, report_meta]                         -- CLASSIFYING_INDEX, or Channel.empty()
        identified_species_ch   // [sample_id, [normalized_species_name, ...]]     -- CLASSIFYING_KRAKEN2

    main:
        // --- "one consensus per (sample, species)" is enforced HERE ---------------
        // Deliberately in MAPPING rather than in the classifier that produces the
        // candidates. MAPPING is what actually spends a consensus on a species, so the
        // invariant holds structurally for anything that reaches it, instead of resting
        // on each classifier remembering to filter itself -- and a classifier needing to
        // know what a *different* classifier found was an odd coupling to begin with.
        // CLASSIFYING_INDEX consequently emits every species it called above the breadth
        // threshold and knows nothing about Kraken2.
        //
        // The cost of moving it here: CLASSIFYING_INDEX resolves a reference record for
        // every candidate, including ones dropped just below. That is wasted
        // SELECT_REFERENCE_RECORD_BY_NAME + EXTRACT_REFERENCE_SUBSET work, and on a
        // sample where a method calls several species and Kraken2 already found most of
        // them, most of that work is wasted. The clean fix if it ever matters is to move
        // reference resolution in here too, behind this filter -- at the price of the
        // symmetric "both classifiers hand over identical shapes" interface, since the
        // Kraken2 side resolves its own references upstream inside SORT_READS_BY_REF.
        //
        // Complete the "already identified" side over the samples that actually have
        // index candidates: exactly one entry each, empty list where CLASSIFYING_KRAKEN2
        // produced no pre-report for that sample (its emit only covers samples with a
        // non-empty one, and a sample Kraken2 found nothing in is exactly the
        // interesting case here). remainder:true fills those with null. The third shape
        // it yields -- a sample Kraken2 found species in but no index method called
        // anything for, [key, null, list] -- is dropped: there is nothing to filter.
        identified_by_sample_ch = index_sample_taxid_ch
            .map { meta, _reads, _ref_fa -> [meta.sample_id, meta.sample_id] }
            .unique()
            .join(identified_species_ch, remainder: true)
            .filter { _sample_id, has_candidates, _identified -> has_candidates != null }
            .map { sample_id, _has_candidates, identified -> [sample_id, identified ?: []] }

        // combine(by: 0), NOT join: many candidates per sample against one
        // identified-species list, and join() pairs keys one-to-one instead of
        // broadcasting. That was the bug fixed in 4973f35 -- moving the filter here must
        // not reintroduce it.
        index_sample_taxid_ch
            .map { meta, reads, ref_fa -> [meta.sample_id, meta, reads, ref_fa] }
            .combine(identified_by_sample_ch, by: 0)
            .filter { _sample_id, meta, _reads, _ref_fa, identified ->
                !identified.contains(meta.species_name.trim().toLowerCase())
            }
            // The trailing _identified matters: combine() left it on the tuple and this
            // closure's parameter count has to match the tuple's width, or Nextflow
            // aborts with "Invalid method invocation `call` with arguments ...".
            .map { _sample_id, meta, reads, ref_fa, _identified -> [meta, reads, ref_fa] }
            .set { index_new_sample_taxid_ch }

        // Same predicate over the report side. Filtered independently rather than
        // semi-joined against the surviving taxid channel: both carry sample_id and
        // species_name, and the predicate is pure, so both necessarily reach the same
        // verdict for a given (sample, species).
        index_report_ch
            .map { join_key, report_meta -> [report_meta.sample_id, join_key, report_meta] }
            .combine(identified_by_sample_ch, by: 0)
            .filter { _sample_id, _join_key, report_meta, identified ->
                !identified.contains(report_meta.species_name.trim().toLowerCase())
            }
            .map { _sample_id, join_key, report_meta, _identified -> [join_key, report_meta] }
            .set { index_new_report_ch }

        // The union both classifiers feed. mix() (not join/combine) because these are
        // disjoint sets of (sample, reference) pairs, not two views of the same one --
        // disjoint because of the filter immediately above, not by assumption.
        sample_taxid_ch = kraken2_sample_taxid_ch.mix(index_new_sample_taxid_ch)
        sample_report_with_join_key_ch = kraken2_report_ch.mix(index_new_report_ch)

        GENERATE_CONSENSUS(sample_taxid_ch)

        GENERATE_CONSENSUS.out.filtered_consensus_ch
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
        GENERATE_CONSENSUS.out.filtered_consensus_ch
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

        GENERATE_CONSENSUS.out.filtered_consensus_ch
            .map { meta, _bam, _bam_idx, _consensus, qc_json -> [meta.id, qc_json] }
            .set { qc_json_simplified_ch }

        GENERATE_CONSENSUS.out.filtered_consensus_ch
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
