include { write_lane_sequence_summary; write_lane_run_summary } from '../modules/write_lane_report.nf'

workflow GENERATE_ABUNDANCE_REPORT {
    /*
    -----------------------------------------------------------------
    Write Abundance Report

    Same shape as GENERATE_ASSEMBLY_REPORT.nf, for the abundance
    estimation lane (Kraken2+Bracken, optionally followed by SCRuB
    decontamination, and ABUNDANCE_ESTIMATION). Bracken's own species
    table already *is* the detailed report; this just carries the
    sample-level counts and stays consistent with the other three
    lanes' reporting shape rather than being a special case.
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

        write_lane_run_summary(all_summaries_pre_ch, "abundance")
        write_lane_run_summary.out.set { publish_run_level_summaries_ch }

    emit:
        publish_seq_level_ch
        publish_run_level_summaries_ch
}
