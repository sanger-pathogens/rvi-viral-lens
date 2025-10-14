// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

include { write_sequence_level_summaries; write_run_level_summaries } from '../modules/write_summary_report.nf'


workflow GENERATE_CLASSIFICATION_REPORT {
    /*
    -----------------------------------------------------------------
    Write Classification Report

    The GENERATE_CLASSIFICATION_REPORT workflow generates some reports
    based on metadata collected during the pipeline:
    - A JSON properties file for each consensus sequence
    - A collated JSON properties files for all consensus sequennce
    - A collated TSV "classification report" file for all consensii

    -----------------------------------------------------------------
    # Inputs

    - **Metadata Channel**: A channel containing metadata for each
    sample.

    # Outputs
        - Per-consensus properties.json file channel
        - Collated properties.json file channel
        - Classification report CSV channel (see code below for columns)
    */

    take:
        report_in_ch // [meta.id, meta]
        qc_json_simplified_ch // [meta.id, per-consensus-qc-json]
        per_con_nexclade_collated_json_in_ch // [meta.id, per-consensus-nextclade-json]

    main:
        report_in_ch // [meta.id, meta]
            .join(qc_json_simplified_ch, remainder: true) // [meta.id, meta, qc_json]
            .join(per_con_nexclade_collated_json_in_ch, remainder: true) //  [meta.id, meta, qc_json, nc_json]
            .set{ report_prep_ch }


        // qc_json_simplified_ch.view{ "QC_JSON_CH: $it \n" }
        // report_prep_ch.view{ "REPORT_PREP_CH: $it \n" }
        write_sequence_level_summaries(report_prep_ch)
        write_sequence_level_summaries.out.set{publish_seq_level_ch}

        write_run_level_summaries(write_sequence_level_summaries.out.collect())
        write_run_level_summaries.out.set{publish_run_level_summaries_ch}

    emit:
        publish_seq_level_ch
        publish_run_level_summaries_ch

}


workflow {
    manifest_channel = Channel.fromPath(params.manifest_file)
    | splitCsv(header: true, sep: ',')
    | map { row ->
        def meta = [[id:row.sample_id,
            taxid:row.taxid,
            ref_selected:row.ref_selected,
            virus_name:row.virus_name,
            virus_subtype:row.virus_subtype,
            flu_segment:row.flu_segment,
            percent_positions_exceeding_depth_10:row.percent_positions_exceeding_depth_10,
            reads_mapped:row.reads_mapped,
            bases_mapped:row.bases_mapped,
            longest_non_n_subsequence:row.longest_non_n_subsequence,
            percent_non_n_bases:row.percent_non_n_bases,
            qc_status:row.qc_status,
            mutations:row.mutations,
            insertions:row.insertions,
            deletions:row.deletions,
            snps:row.snps,
            ti_tv_ratio:row.ti_tv_ratio
        ]]
    }

    GENERATE_CLASSIFICATION_REPORT(manifest_channel)
}

