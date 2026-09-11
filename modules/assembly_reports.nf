// Report builders for the de novo assembly + viral binning lane.
//
// Two stages, matching the lane's own fan-in shape: a per-sample reduction of
// each sample's geNomad/vRhyme/CheckV outputs, then a run-level merge that
// concatenates those and joins in vContact3's batch-level assignments.
//
// Both processes run the same rvi_toolbox/bin/assembly_reports.py; see that
// script's docstring for the exact column semantics.

// Reduce one sample's lane outputs to a small JSON of report rows. Purely
// sample-local -- no vContact3 involvement, since vContact3 runs once for the
// whole batch and is only available to the merge stage below.
//
// --bin-quality is genuinely optional: CheckV only runs on linked bins for
// samples where vRhyme produced at least one bin. Samples without one are
// handed the assets/NO_FILE placeholder, and the flag is dropped.
process ASSEMBLY_REPORT_PER_SAMPLE {
    tag "${meta.id}"

    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    container "quay.io/gsu-pipelines/rvi-vp-basecontainer"

    input:
    tuple val(meta), path(genomad_summary), path(vrhyme_membership), path(scaffold_quality), path(bin_quality)
    path(reports_script)

    output:
    tuple val(meta), path("${meta.id}.assembly_report_parts.json"), emit: parts

    script:
    def bin_quality_arg = bin_quality.name != 'NO_FILE' ? "--bin-quality ${bin_quality}" : ''
    """
    python3 ${reports_script} per-sample \\
        --sample-id         ${meta.id} \\
        --genomad-summary   ${genomad_summary} \\
        --vrhyme-membership ${vrhyme_membership} \\
        --scaffold-quality  ${scaffold_quality} \\
        ${bin_quality_arg} \\
        --output            ${meta.id}.assembly_report_parts.json
    """
}

// Concatenate every sample's parts into the lane's three run-level CSVs,
// joining in the taxonomy vContact3 assigned to each query genome.
//
// vcontact3_assignments is the postprocessed (query-genomes-only) CSV. It is
// optional: the lane still produces all three reports when vContact3 did not
// run, with the vContact3 taxonomy columns left blank.
process ASSEMBLY_REPORTS_MERGE {
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    container "quay.io/gsu-pipelines/rvi-vp-basecontainer"

    publishDir "${params.outdir}", mode: 'copy', overwrite: true, pattern: "assembly_sample_summary_report.csv"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, pattern: "assembly_scaffold_summary_report.csv"
    publishDir "${params.outdir}", mode: 'copy', overwrite: true, pattern: "vmag_scaffold_summary_report.csv"

    input:
    path(parts_json_files)
    path(vcontact3_assignments)
    path(reports_script)

    output:
    path 'assembly_sample_summary_report.csv',   emit: sample_summary
    path 'assembly_scaffold_summary_report.csv', emit: scaffold_summary
    path 'vmag_scaffold_summary_report.csv',     emit: vmag_summary

    script:
    def vcontact3_arg = vcontact3_assignments.name != 'NO_FILE' ? "--vcontact3-assignments ${vcontact3_assignments}" : ''
    """
    python3 ${reports_script} merge \\
        --parts        ${parts_json_files} \\
        ${vcontact3_arg} \\
        --out-sample   assembly_sample_summary_report.csv \\
        --out-scaffold assembly_scaffold_summary_report.csv \\
        --out-vmag     vmag_scaffold_summary_report.csv
    """
}
