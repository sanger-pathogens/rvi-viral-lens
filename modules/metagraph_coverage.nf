// Genome breadth/depth QC for metagraph-align species-call validation: reuses
// SAMTOOLS_COVERAGE (samtools_coverage.nf) for the raw per-reference metrics, then joins
// them against the per-species hit counts and accession->species map from
// metagraph_species_call.nf into one row per called species.

params.script_src_path = "${projectDir}/bin/"

process AGGREGATE_METAGRAPH_COVERAGE {
    tag "${meta.id}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    container "quay.io/gsu-pipelines/rvi-vp-basecontainer"

    // output_subdir: see the same parameter on CALL_METAGRAPH_SPECIES
    // (metagraph_species_call.nf) -- same reason, same two callers.
    publishDir "${params.outdir}/${meta.id}/sequenceindex/${output_subdir}", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(coverage), path(query_lengths), path(index_label_map), path(species_hits)
    val(output_subdir)

    output:
    tuple val(meta), path("${meta.id}_metagraph_map_qc.tsv"), emit: qc_table

    script:
    """
    ${params.script_src_path}aggregate_metagraph_coverage.py \\
        --coverage ${coverage} \\
        --query-lengths ${query_lengths} \\
        --index-label-map ${index_label_map} \\
        --species-hits ${species_hits} \\
        --sample-id ${meta.id} \\
        --out ${meta.id}_metagraph_map_qc.tsv
    """
}

process GENERATE_METAGRAPH_MAP_SUMMARY {
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'

    // summary_name: 'metagraph_map_summary' (align) or 'metagraph_query_map_summary'
    // (query) -- distinct run-level filenames for the same reason as output_subdir above.
    publishDir "${params.outdir}/${summary_name}", mode: 'copy', overwrite: true

    input:
    path(qc_tables)
    val(summary_name)

    output:
    path("${summary_name}.tsv")

    script:
    """
    head -n1 \$(ls *_metagraph_map_qc.tsv | head -n1) > ${summary_name}.tsv
    tail -q -n +2 *_metagraph_map_qc.tsv >> ${summary_name}.tsv
    """
}
