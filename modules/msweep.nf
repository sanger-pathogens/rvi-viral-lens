// mSWEEP group abundance estimation from Themisto2 pseudoalignments.
// Adapted from the gemsweep pipeline (themisto2 branch). Results are published per
// sample under <outdir>/<sample>/sequenceindex/msweep.

process MSWEEP {
    tag "${meta.id}"
    label 'cpu_16'

    container 'quay.io/biocontainers/msweep:2.2.1--h503566f_1'

    memory { 40.GB * task.attempt }
    // Only escalate time after an LSF runlimit/OOM kill (exit 140); queue is left to the
    // global time-aware selector in nextflow-commons.config rather than overridden here,
    // so it advances (normal->long->week) in step with the escalated time.
    time { task.attempt > 1 ? ( task.previousTrace?.exit == 140 ? task.previousTrace.time * 2 : (task.previousTrace?.time ?: 12.h)) : (12.h) }

    // pattern restricts publishing to the actual mSWEEP results; without it, the
    // pseudoalignments passthrough emit below (kept only so CLEANUP_THEMISTO_PSEUDOALIGNMENTS
    // can consume it downstream) would get published too.
    publishDir mode: 'copy', path: "${params.outdir}/${meta.id}/sequenceindex/msweep", pattern: "*_mSWEEP_*"

    input:
    tuple val(meta), path(pseudoalignment_1), path(pseudoalignment_2)
    path(ref_groups)

    output:
    tuple val(meta),
          path("${meta.id}_mSWEEP_abundances.txt"),
          path("${meta.id}_mSWEEP_probs.tsv"), emit: abundances
    tuple val(meta), path(pseudoalignment_1), path(pseudoalignment_2), emit: pseudoalignments

    script:
    output_prefix = "${meta.id}_mSWEEP"

    """
    mSWEEP --themisto-1 ${pseudoalignment_1} --themisto-2 ${pseudoalignment_2} -o ${output_prefix} -i ${ref_groups} --write-probs -t ${task.cpus} --min-hits 100
    """
}
