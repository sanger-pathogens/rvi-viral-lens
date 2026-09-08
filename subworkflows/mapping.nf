// --- map reads to taxid, generate consensus, classify (rvi_integration_1) ---
// Extracted unchanged from main.nf's inline body (was mirrored by
// mapping_pipeline_main.nf, the legacy standalone entry-point for this lane).
include {SORT_READS_BY_REF} from '../workflows/SORT_READS_BY_REF.nf'
include {GENERATE_CONSENSUS} from '../workflows/GENERATE_CONSENSUS.nf'
include {SCOV2_SUBTYPING} from '../workflows/SCOV2_SUBTYPING.nf'
include {GENERATE_CLASSIFICATION_REPORT} from '../workflows/GENERATE_CLASSIFICATION_REPORT.nf'
include {RUN_NEXTCLADE} from '../workflows/RUN_NEXTCLADE.nf'
include {publish_consensus_files as publish_aln_files; publish_consensus_files as publish_nc_files; publish_consensus_files as publish_per_sample_json} from '../modules/publish_lite.nf'
include {publish_run_files} from '../modules/publish_lite.nf'

workflow MAPPING {
    /*
    -----------------------------------------------------------------
    Maps preprocessed reads to a per-sample reference taxid
    (SORT_READS_BY_REF), generates a consensus sequence
    (GENERATE_CONSENSUS), optionally runs Nextclade and SARS-CoV-2
    subtyping, and writes the final per-sample classification report.
    -----------------------------------------------------------------
    */

    take:
        preprocessed_3tuple_ch // tuple (meta, read1, read2)

    main:
        // reconstruct the tuple(meta, [read1, read2]) shape SORT_READS_BY_REF expects
        preprocessed_3tuple_ch
            .map { meta, read1, read2 -> [meta, [read1, read2]] }
            .set { sort_reads_in_ch }

        SORT_READS_BY_REF(sort_reads_in_ch)
        GENERATE_CONSENSUS(SORT_READS_BY_REF.out.sample_taxid_ch)

        GENERATE_CONSENSUS.out.filtered_consensus_ch
            .map { meta, _bam, _bam_idx, consensus, _qc_json -> [meta.id, meta, consensus] }
            .set { consensus_fa_ch }

        SORT_READS_BY_REF.out.sample_pre_report_ch
            .map { meta ->
                def new_meta = meta + [id: "${meta.sample_id}.${meta.selected_taxid}"]
                [new_meta.id, new_meta]
            }
            .combine(consensus_fa_ch, by: 0)
            .map { _id, pre_report_meta, fa_meta, fa ->
                def final_meta = pre_report_meta + [reference_header: "${fa_meta.reference_header}", taxid: "${fa_meta.taxid}"]
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

        // branching output from generate_consensus for viral specific subtyping
        SORT_READS_BY_REF.out.sample_pre_report_ch
            .map { it -> ["${it.sample_id}.${it.selected_taxid}".toString(), it] }
            .set { sample_report_with_join_key_ch }

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

        // --- rvi_integration_1: per-sample "species already identified by MAPPING" ---
        // Consumed by SEQUENCE_INDEX to decide which of its own species calls are
        // genuinely new. Built straight off SORT_READS_BY_REF's raw per-sample
        // pre-report FILE (one element per sample, available as soon as THAT sample's
        // Kraken2/k2r pass finishes) rather than the exploded sample_pre_report_ch +
        // groupTuple() -- groupTuple() can't emit a group until its whole upstream
        // channel closes, which would mean waiting for every sample in the run, not
        // just this one.
        //
        // Match on virus_name, NOT ref_selected. Both are free-text names lifted from the
        // Kraken2 report by bin/k2r_report.py, but they sit at different ranks and only
        // virus_name shares a vocabulary with the sequence-index methods'
        // species_label/species fields:
        //   virus_name   <- kraken2ref's source_taxid, rank S      e.g. "Betacoronavirus pandemicum"
        //   ref_selected <- the selected reference, rank S1/S2/S3  e.g. "Severe acute
        //                                                          respiratory syndrome coronavirus 2"
        // mSWEEP/Metagraph label the RVDB index with ICTV species binomials, so comparing
        // against ref_selected alone never matched: every species MAPPING had already
        // found looked "new", and the feature called a redundant consensus for it rather
        // than a new-species one. Proven on the farm -- see INSTRUCT.md's
        // "new-species consensus: matched on the wrong column" section.
        //
        // virus_name is reliably species rank: kraken2ref decomposes species -> below-species
        // by construction, and its decomposed JSON carries the rank code in its own `source`
        // field ('S' for all 68 reference selections across the 10 samples run so far).
        //
        // ref_selected is still folded into the same set, as a widening rather than a
        // replacement: a candidate whose label happens to match a strain-level reference
        // Kraken2 already selected is also genuinely already covered, so suppressing it is
        // correct too. Adding it can only ever suppress a candidate, never invent one.
        SORT_READS_BY_REF.out.raw_sample_pre_report_ch
            .filter { it -> it.size() > 1 } // mirror SORT_READS_BY_REF's own empty-file guard
            .map { report_file ->
                def lines = report_file.readLines()
                def header = lines[0].split('\t')
                def sample_id_idx    = header.findIndexOf { String col -> col == 'sample_id' }
                def virus_name_idx   = header.findIndexOf { String col -> col == 'virus_name' }
                def ref_selected_idx = header.findIndexOf { String col -> col == 'ref_selected' }
                // Fail loudly rather than silently mis-matching: Groovy's row[-1] returns the
                // LAST field, so a renamed/removed column would quietly compare against
                // report_name instead of erroring.
                if (sample_id_idx < 0 || virus_name_idx < 0 || ref_selected_idx < 0) {
                    error("pre-report ${report_file} lacks one of the sample_id/virus_name/" +
                          "ref_selected columns (header: ${header}). bin/k2r_report.py's output " +
                          "format has changed -- update identified_species_ch in subworkflows/mapping.nf.")
                }
                // split('\t') drops trailing empty fields, so a row is only usable if it
                // actually reaches the columns being read.
                def max_idx = [sample_id_idx, virus_name_idx, ref_selected_idx].max()
                def rows = lines[1..-1]
                    .collect { line -> line.split('\t') }
                    .findAll { row -> row.size() > max_idx }
                if (!rows) {
                    // Not reachable via bin/k2r_report.py today: a sample with nothing
                    // selected yields a column-less 1-byte file, already dropped by the
                    // filter above, and any real row carries all 11 columns. Fail loudly
                    // instead of letting rows[0] throw an opaque IndexOutOfBounds.
                    error("pre-report ${report_file} has a header but no parseable data rows " +
                          "(need >${max_idx} tab-separated fields). Check bin/k2r_report.py's output.")
                }
                def sample_id = rows[0][sample_id_idx]
                def species = rows
                    .collectMany { row -> [row[virus_name_idx], row[ref_selected_idx]] }
                    .collect { name -> name.trim().toLowerCase() }
                    .findAll { name -> name }
                    .unique()
                [sample_id, species]
            }
            .set { identified_species_ch }

    emit:
        identified_species_ch // [sample_id, [normalized_species_name, ...]]
}
