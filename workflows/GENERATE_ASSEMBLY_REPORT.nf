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

        // Publish only the run-level JSON. write_lane_run_summary also writes
        // <lane>_summary_report.csv, but the assembly lane's sample-level CSV is
        // now assembly_sample_summary_report.csv, built by ASSEMBLY_REPORTS from
        // the modules' own outputs (subworkflows/assembly.nf). Publishing both
        // would ship two sample-level CSVs whose columns had already diverged.
        // The shared writer is left alone -- the mapping and abundance lanes
        // still use its CSV.
        write_lane_run_summary.out
            .map { run_summary_json, _summary_report_csv -> run_summary_json }
            .set { publish_run_level_summaries_ch }

    emit:
        publish_seq_level_ch
        publish_run_level_summaries_ch
}
