// Heatmap figure comparing relative abundance before/after SCRuB decontamination, plus
// the per-species change, for a quick visual QC of the SCRuB decontamination subworkflow.
// Runs in the same container as SCRUB itself (already has R/ggplot2/dplyr/tidyr/scales),
// no additional image to pull.

params.script_src_path = "${projectDir}/bin/"

process SCRUB_HEATMAP {
    label 'cpu_1'
    label 'mem_4'
    label 'time_1'

    container 'quay.io/sangerpathogens/scrub:fcbb852'

    publishDir "${params.outdir}/abundance/scrub", mode: 'copy', overwrite: true

    input:
    path(formatted_bracken)
    path(scrub_output)

    output:
    path("decontam_heatmap.png"), emit: heatmap

    script:
    """
    ${params.script_src_path}plot_scrub_heatmap.R \\
        --bracken-tsv ${formatted_bracken} \\
        --scrub-rds ${scrub_output}/scrub_result.rds \\
        --top-n ${params.scrub_heatmap_top_n} \\
        --out decontam_heatmap.png
    """
}
