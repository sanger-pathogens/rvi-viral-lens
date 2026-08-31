// SCRuB (Shenhav & Korem labs; Austin et al., Nat Biotechnol 2023) cross-contamination
// detection/correction, run once per pipeline run against the whole-run species-abundance
// matrix and a user-supplied plate map. run_scrub.R's CLI args map directly onto the
// public SCRuB() R function's arguments (input_data, metadata, control_order,
// dist_threshold); --output-dir is added by this wrapper to serialize results to disk.

process SCRUB {
    label 'cpu_1'
    label 'mem_4'
    label 'time_1'

    container 'quay.io/sangerpathogens/scrub:fcbb852'

    publishDir "${params.outdir}/abundance/scrub", mode: 'copy', overwrite: true

    input:
    path(formatted_bracken)
    path(plate_map)

    output:
    path("scrub_out"), emit: scrub_output

    script:
    """
    run_scrub.R \\
        --input-data ${formatted_bracken} \\
        --metadata ${plate_map} \\
        --control-order "${params.scrub_control_order}" \\
        --dist-threshold ${params.scrub_dist_threshold} \\
        --output-dir scrub_out
    """
}
