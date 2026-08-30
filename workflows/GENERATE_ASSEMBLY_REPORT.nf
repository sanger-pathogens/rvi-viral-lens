include { write_lane_sequence_summary; write_lane_run_summary } from '../modules/write_lane_report.nf'

workflow GENERATE_ASSEMBLY_REPORT {
    /*
    -----------------------------------------------------------------
    Write Assembly Report

    Mirrors GENERATE_CLASSIFICATION_REPORT.nf's shape for the de novo
    assembly + viral binning lane. Unlike that report, there's no
    separate qc/nextclade JSON to merge in here: `meta` arrives already
    carrying every module's sample-level counts (genomad_*, vrhyme_*,
    checkv_*, vcontact3_*), added upstream in main.nf via meta.plus().
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

        write_lane_run_summary(all_summaries_pre_ch, "assembly")
        write_lane_run_summary.out.set { publish_run_level_summaries_ch }

    emit:
        publish_seq_level_ch
        publish_run_level_summaries_ch
}
