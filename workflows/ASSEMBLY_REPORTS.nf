include { ASSEMBLY_REPORT_PER_SAMPLE; ASSEMBLY_REPORTS_MERGE } from '../modules/assembly_reports.nf'

workflow ASSEMBLY_REPORTS {
    /*
    -----------------------------------------------------------------
    Assembly lane reports

    Builds the de novo assembly + viral binning lane's three run-level
    CSVs from the lane's own outputs:

      assembly_sample_summary_report.csv    one row per sample
      assembly_scaffold_summary_report.csv  one row per viral scaffold
      vmag_scaffold_summary_report.csv      one row per vRhyme bin (vMAG)

    Everything sample-local is reduced per sample first; the merge stage
    then concatenates those and joins in vContact3's taxonomy, which only
    exists at batch level (vContact3 runs once for the whole run, not per
    sample).
    -----------------------------------------------------------------
    # Inputs

    - **genomad_summary_ch**:   tuple(meta, virus_summary.tsv)                  -- GENOMAD_CLASSIFY.out.virus_summary
    - **membership_ch**:        tuple(meta, membership.tsv)                     -- VRHYME_BIN.out.membership
    - **scaffold_quality_ch**:  tuple(meta, virus_scaffolds_quality_summary)    -- CHECKV_QC.out.virus_scaffolds_quality_summary
    - **bin_quality_ch**:       tuple(meta, linked_bins_quality_summary)        -- CHECKV_QC.out.linked_bins_quality_summary
                                (only samples where vRhyme produced a bin)
    - **vcontact3_assignments_ch**: final_assignments_postprocessed.csv         -- VCONTACT3_RUN.out.postprocessed_assignments
                                (may be empty if vContact3 did not run)

    # Outputs
        - The three report CSVs, published to params.results_dir
    -----------------------------------------------------------------
    */

    take:
        genomad_summary_ch
        membership_ch
        scaffold_quality_ch
        bin_quality_ch
        vcontact3_assignments_ch

    main:
        // Placeholder standing in for the two genuinely optional inputs. Nextflow
        // has no null `path`, so "absent" has to be spelled as a real file whose
        // name the process recognises and drops the corresponding flag for.
        def no_file = file("${projectDir}/assets/NO_FILE")

        // CheckV only assesses linked bins for samples that have at least one bin,
        // so bin_quality_ch is missing entries the other three channels have.
        // remainder: true keeps those samples in rather than silently dropping them
        // from the reports; their vMAG rows are simply empty.
        genomad_summary_ch
            .join(membership_ch)
            .join(scaffold_quality_ch)
            .join(bin_quality_ch, remainder: true)
            .map { meta, genomad, membership, scaffold_q, bin_q ->
                [meta, genomad, membership, scaffold_q, bin_q ?: no_file]
            }
            .dump(tag: 'assembly_reports_per_sample_inputs')
            .set { ch_per_sample_inputs }

        reports_script_ch = Channel.value(file("${projectDir}/bin/assembly_reports.py"))

        ASSEMBLY_REPORT_PER_SAMPLE(ch_per_sample_inputs, reports_script_ch)

        // Run-level barrier: every sample's parts must be in before the merge.
        ASSEMBLY_REPORT_PER_SAMPLE.out.parts
            .map { _meta, json -> json }
            .collect()
            .dump(tag: 'assembly_reports_collected_parts')
            .set { all_parts_ch }

        // vContact3 is optional for the lane; fall back to the placeholder so the
        // reports are still produced (with blank vContact3 columns) without it.
        vcontact3_assignments_ch
            .ifEmpty { no_file }
            .set { ch_vcontact3_assignments }

        ASSEMBLY_REPORTS_MERGE(all_parts_ch, ch_vcontact3_assignments, reports_script_ch)

    emit:
        sample_summary   = ASSEMBLY_REPORTS_MERGE.out.sample_summary
        scaffold_summary = ASSEMBLY_REPORTS_MERGE.out.scaffold_summary
        vmag_summary     = ASSEMBLY_REPORTS_MERGE.out.vmag_summary
}
