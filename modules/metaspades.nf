process METASPADES {
    tag "${meta.id}"
    label 'cpu_16'
    label "mem_${params.metaspades_base_mem_gb}"
    label 'time_12'

    container 'quay.io/biocontainers/spades:3.15.5--h95f258a_1'

    publishDir "${params.results_dir}/${basemetaID}/metaspades/", mode: 'copy', overwrite: true, pattern: 'contigs.fasta', saveAs: { filename -> "${meta.id}_contigs.fa" }
    publishDir "${params.results_dir}/${basemetaID}/metaspades/", mode: 'copy', overwrite: true, pattern: 'scaffolds.fasta', saveAs: { filename -> "${meta.id}_scaffolds.fa" }

    input:
    tuple val(meta), path(R1), path(R2)

    output:
    tuple val(meta), path("contigs.fasta"), path("scaffolds.fasta"), emit: contigs_channel
    path("spades.log"), emit: log_channel

    script:
    splitmeta = "${meta.id}".split("_subsampled-")
    basemetaID = "${splitmeta[0]}"
    """
    metaspades.py -1 ${R1} -2 ${R2} -o . --tmp-dir tmp/
    status=\${?}
    if [ \${status} -gt 0 ] ; then
        # try and catch known errors from spades.log
        ## empty input read file - exit 7
        grep '== Error ==  file is empty' spades.log && exit 7
        ## segmentation fault, possibly due to farm environment and spades not being compiled against the machine/in the singularity container - exit 3
        grep '== Error ==  system call for:.\\+/usr/local/bin/spades-hammer.\\+finished abnormally' spades.log 1>&2 && exit 3
    else
        # try and catch known warnings from the spades.log when not causing non-zero exit code
        ## empty output contigs.fasta file and no scaffold file, often meaning low read input - exit 7
        grep '======= SPAdes pipeline finished WITH WARNINGS!' spades.log 1>&2 \
          && grep ' * Assembled contigs are in .\\+contigs.fasta' spades.log 1>&2 \
          && [ ! -s contigs.fasta ] && exit 7 
        # NB: in case scaffolds.fasta file is missing not due to the above, nextflow will error out as expecting it as ouput file
    fi
    # if not caught known exception, process should not have exited yet - do it now with stored metaspades exit status
    exit \${status}
    """
}  