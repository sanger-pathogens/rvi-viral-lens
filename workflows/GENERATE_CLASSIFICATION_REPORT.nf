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
        report_prep_ch //  [meta.id, meta, qc_json, nc_json]

    main:
        write_sequence_level_summaries(report_prep_ch)
        write_sequence_level_summaries.out.set{publish_seq_level_ch}

        write_sequence_level_summaries.out
            .map{meta, per_con_json -> [per_con_json]}
            .collect()
            .set{all_summaries_pre_ch}

        write_run_level_summaries(all_summaries_pre_ch)
        write_run_level_summaries.out.set{publish_run_level_summaries_ch}

    emit:
        publish_seq_level_ch
        publish_run_level_summaries_ch

}

