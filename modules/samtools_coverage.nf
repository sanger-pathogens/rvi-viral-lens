// Genome breadth/depth QC for mSWEEP low-abundance hit validation: samtools coverage
// gives, per reference sequence, breadth (% covered), mean depth, read count, mean base
// quality and mean mapping quality in one command; a second samtools view + awk pass sums
// aligned query (read) bases per reference sequence, since samtools coverage doesn't
// report that. Since SELECT_REFERENCE_RECORDS now picks exactly one sequence per species,
// aggregate_species_coverage.py joins these straight through to one row per species —
// no more rolling up multiple sequences per label.

params.script_src_path = "${projectDir}/bin/"

process SAMTOOLS_COVERAGE {
    tag "${meta.id}"
    label 'cpu_1'
    label 'mem_2'
    label 'time_1'

    container 'quay.io/sangerpathogens/bowtie2-samtools:1.1-c1'

    input:
    tuple val(meta), path(sorted_bam)

    output:
    tuple val(meta), path("${meta.id}_coverage.tsv"), emit: coverage
    tuple val(meta), path("${meta.id}_query_lengths.tsv"), emit: query_lengths

    script:
    """
    samtools coverage ${sorted_bam} > ${meta.id}_coverage.tsv
    samtools view -F 4 ${sorted_bam} | awk -F'\\t' '{sum[\$3]+=length(\$10)} END{for (r in sum) print r"\\t"sum[r]}' > ${meta.id}_query_lengths.tsv
    """
}

process AGGREGATE_SPECIES_COVERAGE {
    tag "${meta.id}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    container "quay.io/gsu-pipelines/rvi-vp-basecontainer"

    publishDir "${params.outdir}/${meta.id}/sequenceindex/msweep_map", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(coverage), path(query_lengths), path(index_label_map), path(groups)

    output:
    tuple val(meta), path("${meta.id}_msweep_map_qc.tsv"), emit: qc_table

    script:
    """
    ${params.script_src_path}aggregate_species_coverage.py \\
        --coverage ${coverage} \\
        --query-lengths ${query_lengths} \\
        --index-label-map ${index_label_map} \\
        --groups ${groups} \\
        --sample-id ${meta.id} \\
        --out ${meta.id}_msweep_map_qc.tsv
    """
}

process GENERATE_MSWEEP_MAP_SUMMARY {
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'

    publishDir "${params.outdir}/msweep_map_summary", mode: 'copy', overwrite: true

    input:
    path(qc_tables)

    output:
    path("msweep_map_summary.tsv")

    script:
    """
    head -n1 \$(ls *_msweep_map_qc.tsv | head -n1) > msweep_map_summary.tsv
    tail -q -n +2 *_msweep_map_qc.tsv >> msweep_map_summary.tsv
    """
}
