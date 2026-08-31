process GENOMAD {
    tag "${meta.id}"
    label 'cpu_8'
    label 'mem_32'
    label 'time_12'

    publishDir "${params.outdir}/${meta.id}/assembly/genomad", mode: 'copy', overwrite: true, pattern: "*_virus_summary.tsv"
    publishDir "${params.outdir}/${meta.id}/assembly/genomad", mode: 'copy', overwrite: true, pattern: "*_virus.fna"
    publishDir "${params.outdir}/${meta.id}/assembly/genomad", mode: 'copy', overwrite: true, pattern: "*_summary.json"
    publishDir "${params.outdir}/${meta.id}/assembly/genomad", mode: 'copy', overwrite: true, pattern: "*_virus_genes.tsv"
    publishDir "${params.outdir}/${meta.id}/assembly/genomad", mode: 'copy', overwrite: true, pattern: "*_virus_proteins.faa"

    container 'quay.io/biocontainers/genomad:1.12.0--pyhdfd78af_0'

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path(virus_summary),  emit: virus_summary
    tuple val(meta), path(virus_fna),      emit: virus_fna
    tuple val(meta), path(summary_json),   emit: summary_json
    tuple val(meta), path(virus_genes),    emit: virus_genes
    tuple val(meta), path(virus_proteins), emit: virus_proteins

    script:
    prefix          = assembly.getBaseName()
    outdir          = "${meta.id}_genomad_output"
    virus_summary   = "${meta.id}_virus_summary.tsv"
    virus_fna       = "${meta.id}_virus.fna"
    summary_json    = "${meta.id}_summary.json"
    virus_genes     = "${meta.id}_virus_genes.tsv"
    virus_proteins  = "${meta.id}_virus_proteins.faa"
    """
    genomad end-to-end \
        --threads ${task.cpus} \
        --cleanup --enable-score-calibration \
        ${assembly} \
        ${outdir} \
        ${params.genomad_db}

    mv ${outdir}/${prefix}_summary/${prefix}_virus_summary.tsv  ${virus_summary}
    mv ${outdir}/${prefix}_summary/${prefix}_virus.fna          ${virus_fna}
    mv ${outdir}/${prefix}_summary/${prefix}_summary.json       ${summary_json}
    mv ${outdir}/${prefix}_summary/${prefix}_virus_genes.tsv    ${virus_genes}
    mv ${outdir}/${prefix}_summary/${prefix}_virus_proteins.faa ${virus_proteins}
    """
}
