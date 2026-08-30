include { write_lane_sequence_summary; write_lane_run_summary } from '../modules/write_lane_report.nf'

workflow GENERATE_MAPPING_REPORT {
    /*
    -----------------------------------------------------------------
    Write Mapping Report

    Same shape as GENERATE_ASSEMBLY_REPORT.nf, for the "map reads to
    sequence indexes" lane (pseudoalignment/sequence-to-graph alignment
    + QC Mapping). `meta` arrives already carrying that lane's QC
    Mapping fields.
    -----------------------------------------------------------------
    # Inputs

    - **report_prep_ch**: tuple(sample_id, meta) -- one per sample.

    # Outputs
        - Per-sample properties.json file channel
        - Run-level summary JSON + CSV channel
    -----------------------------------------------------------------
    */

    take:
        report_prep_ch // tuple(sample_id, meta)

    main:
        write_lane_sequence_summary(report_prep_ch)
        write_lane_sequence_summary.out.set { publish_seq_level_ch }

        write_lane_sequence_summary.out
            .map { meta, per_sample_json -> per_sample_json }
            .collect()
            .set { all_summaries_pre_ch }

        write_lane_run_summary(all_summaries_pre_ch, "mapping")
        write_lane_run_summary.out.set { publish_run_level_summaries_ch }

    emit:
        publish_seq_level_ch
        publish_run_level_summaries_ch
}
