// Themisto2 pseudoalignment against a pre-built (.thm2) viral index.
// Adapted from the gemsweep pipeline (themisto2 branch); only the pseudoalignment
// process is required for the prebuilt-index mSWEEP workflow.

process THEMISTO_PSEUDOALIGN {
    tag "${meta.id}"
    label 'cpu_16'
    label 'time_queue_from_normal'

    // The upstream gemsweep module ships without a container (placeholder TODO).
    container 'quay.io/sangerpathogens/themisto2:0.0.1'

    // 25GB on the first attempt, escalating on retry. The upstream module asked for
    // 100GB up front, which reserves far more than the run needs and makes the job
    // queue behind large-memory slots for no benefit.
    memory { 25.GB * task.attempt }

    input:
    tuple val(meta), path(reads_1), path(reads_2)
    path index_files    // For staging
    val index_prefix    // For use in command

    output:
    tuple val(meta), path("pseudoalignments_1.aln.gz"), path("pseudoalignments_2.aln.gz")

    script:
    pseudoalignment_params = "-i ${index_prefix}.thm2 --n-threads ${task.cpus}"

    """
    themisto2 intersection-pseudoalign -q ${reads_1} ${pseudoalignment_params} --sort-output --themisto1-output-format | gzip > pseudoalignments_1.aln.gz
    themisto2 intersection-pseudoalign -q ${reads_2} ${pseudoalignment_params} --sort-output --themisto1-output-format | gzip > pseudoalignments_2.aln.gz
    """
}
