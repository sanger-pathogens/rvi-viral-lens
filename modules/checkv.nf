process CHECKV {
    tag "${meta.id}:${source}"
    label 'cpu_4'
    label 'mem_16'
    label 'time_12'

    container 'quay.io/biocontainers/checkv:1.1.1--pyh106432d_0'

    publishDir "${params.results_dir}/${meta.id}/checkv", mode: 'copy', overwrite: true, pattern: "${source}/quality_summary.tsv",  saveAs: { f -> "${source}_quality_summary.tsv" }
    publishDir "${params.results_dir}/${meta.id}/checkv", mode: 'copy', overwrite: true, pattern: "${source}/completeness.tsv",     saveAs: { f -> "${source}_completeness.tsv" }
    publishDir "${params.results_dir}/${meta.id}/checkv", mode: 'copy', overwrite: true, pattern: "${source}/contamination.tsv",    saveAs: { f -> "${source}_contamination.tsv" }
    publishDir "${params.results_dir}/${meta.id}/checkv", mode: 'copy', overwrite: true, pattern: "${source}/complete_genomes.tsv", saveAs: { f -> "${source}_complete_genomes.tsv" }

    input:
    tuple val(meta), val(source), path(fasta)

    // completeness/contamination/complete_genomes stay optional: scaffolds may
    // pass the upstream length filter yet yield no bins (a normal outcome when
    // binning signal is too sparse), in which case CheckV simply doesn't write
    // some of these (so not a failure), and nothing downstream needs a
    // guaranteed real path for those three.
    // quality_summary is NOT optional: the script below guarantees a real
    // (possibly empty) file even when CheckV wrote nothing, so consumers (e.g.
    // vContact3's standalone-scaffold-promotion check) never have to handle a
    // missing quality_summary as a special case.
    output:
    tuple val(meta), val(source), path("${source}/quality_summary.tsv"),  emit: quality_summary
    tuple val(meta), val(source), path("${source}/completeness.tsv"),     emit: completeness,     optional: true
    tuple val(meta), val(source), path("${source}/contamination.tsv"),    emit: contamination,    optional: true
    tuple val(meta), val(source), path("${source}/complete_genomes.tsv"), emit: complete_genomes, optional: true

    script:
    """
    checkv end_to_end \\
        ${fasta} \\
        ${source} \\
        -t ${task.cpus} \\
        -d ${params.checkv_db}

    mkdir -p ${source}
    touch ${source}/quality_summary.tsv
    """
}
