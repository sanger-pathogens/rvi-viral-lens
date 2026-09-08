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
}
