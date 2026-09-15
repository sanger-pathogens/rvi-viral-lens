include { write_lane_sequence_summary } from '../modules/write_lane_report.nf'

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

        // No run-level report from this lane. assembly_run_summary.json used to be
        // written here and is gone: ASSEMBLY_REPORTS (subworkflows/assembly.nf) now
        // builds the lane's run-level CSVs straight from the modules' own outputs,
        // and the JSON was not merely a duplicate of those -- it was wrong. It came
        // off the inner-join chain below, which drops any sample vRhyme never ran
        // for, so a 95-sample run produced a 35-record JSON while
        // assembly_sample_summary_report.csv correctly held all 95.
        //
        // The per-sample properties.json above is still published, unchanged.

    emit:
        publish_seq_level_ch
}
