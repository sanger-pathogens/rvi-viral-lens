// Genome breadth/depth QC for Themisto-hit species-call validation: reuses
// SAMTOOLS_COVERAGE (samtools_coverage.nf) for the raw per-reference metrics, then joins
// them against the per-species hit counts and index->species map from
// themisto_species_call.nf into one row per called species.

params.script_src_path = "${projectDir}/bin/"

process AGGREGATE_THEMISTO_COVERAGE {
    tag "${meta.id}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    container "quay.io/gsu-pipelines/rvi-vp-basecontainer"

    publishDir "${params.outdir}/${meta.id}/sequenceindex/themisto_map", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(coverage), path(query_lengths), path(index_label_map), path(species_hits)

    output:
    tuple val(meta), path("${meta.id}_themisto_map_qc.tsv"), emit: qc_table

    script:
    """
    ${params.script_src_path}aggregate_themisto_coverage.py \\
        --coverage ${coverage} \\
        --query-lengths ${query_lengths} \\
        --index-label-map ${index_label_map} \\
        --species-hits ${species_hits} \\
        --sample-id ${meta.id} \\
        --out ${meta.id}_themisto_map_qc.tsv
    """
}

process GENERATE_THEMISTO_MAP_SUMMARY {
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'

    publishDir "${params.outdir}/themisto_map_summary", mode: 'copy', overwrite: true

    input:
    path(qc_tables)

    output:
    path("themisto_map_summary.tsv")

    script:
    """
    head -n1 \$(ls *_themisto_map_qc.tsv | head -n1) > themisto_map_summary.tsv
    tail -q -n +2 *_themisto_map_qc.tsv >> themisto_map_summary.tsv
    """
}
