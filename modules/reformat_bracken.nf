// Reformats bracken_summary_report.tsv (KRAKEN2BRACKEN's whole-run summary, taxa as rows
// spanning every rank, samples as columns) into SCRuB's expected input: species-level
// rows only, transposed to samples x species. Samples absent from the plate map are
// dropped (with a warning), since SCRuB requires the data and metadata row names to match.
// Optionally (params.scrub_zero_species_in_controls), forces named species' counts to 0
// in every control sample as a pre-cleaning step ahead of SCRuB itself.

params.script_src_path = "${projectDir}/bin/"

process REFORMAT_BRACKEN_FOR_SCRUB {
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    container "quay.io/gsu-pipelines/rvi-vp-basecontainer"

    input:
    path(bracken_summary)
    path(plate_map)

    output:
    path("bracken_for_scrub.tsv"), emit: formatted

    script:
    """
    ${params.script_src_path}reformat_bracken_for_scrub.py \\
        --bracken-summary ${bracken_summary} \\
        --sample-suffix "${params.scrub_sample_suffix}" \\
        --plate-map ${plate_map} \\
        --zero-species-in-controls "${params.scrub_zero_species_in_controls}" \\
        --out bracken_for_scrub.tsv
    """
}
