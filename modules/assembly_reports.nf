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
// --vrhyme-membership and --bin-quality are both genuinely optional: vRhyme emits
// nothing for a sample it cannot bin, and CheckV only assesses linked bins where a
// bin exists. Samples without either are handed the empty assets/NO_FILE placeholder,
// which the script reads as zero rows.
process ASSEMBLY_REPORT_PER_SAMPLE {
    tag "${meta.id}"

    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    container "quay.io/gsu-pipelines/rvi-vp-basecontainer"

    input:
    // Every input is staged under a fixed name of its own, for two reasons that would
    // otherwise each cause a filename collision:
    //   - CheckV writes both its runs as '<source>/quality_summary.tsv', so the
    //     per-scaffold and linked-bins files share a basename.
    //   - the NO_FILE placeholder stands in for more than one absent input, so a sample
    //     missing both its membership and its linked-bins QC would stage NO_FILE twice.
    tuple val(meta), \
          path(genomad_summary,   stageAs: 'genomad_virus_summary.tsv'), \
          path(vrhyme_membership, stageAs: 'vrhyme_membership.tsv'), \
          path(scaffold_quality,  stageAs: 'scaffold_quality_summary.tsv'), \
          path(bin_quality,       stageAs: 'bin_quality_summary.tsv')
    path(reports_script)

    output:
    tuple val(meta), path("${meta.id}.assembly_report_parts.json"), emit: parts

    script:
    // Every flag is passed unconditionally: an absent input arrives as the empty NO_FILE
    // placeholder, which the script reads as zero rows -- the same result as omitting it,
    // without the module having to decide which spelling means "absent".
    """
    python3 ${reports_script} per-sample \\
        --sample-id         ${meta.id} \\
        --genomad-summary   ${genomad_summary} \\
        --vrhyme-membership ${vrhyme_membership} \\
        --scaffold-quality  ${scaffold_quality} \\
        --bin-quality       ${bin_quality} \\
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
