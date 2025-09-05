process run_nextclade {

    tag "${meta.tag_id}"
    label "nextclade"

    input:
    tuple val(meta), path(input_fa), path(data_dir)

    output:
    tuple val(meta), path("*.csv"), path("*.tar.gz")

    script:
    def nextclade_output_verbosity = params.nextclade_output_verbosity
    //def data_dir = "${meta.data_dir}"
    def ref_fasta = "${data_dir}/reference.fasta"
    def ref_tree = "${data_dir}/tree.json"
    def ref_gff = "${data_dir}/genome_annotation.gff3"
    //def input_fa = "${meta.consensus_path}"

    if (file(ref_tree).exists())
        """
        #!/bin/bash
        set -e

        nextclade run \
            -r ${ref_fasta} \
            -a ${ref_tree} \
            -m ${ref_gff} \
            -O ${meta.tag_id}_nextclade \
            -s ${nextclade_output_verbosity} \
            -n ${meta.tag_id} \
            --include-reference true \
            ${input_fa}

        cp ./${meta.tag_id}_nextclade/${meta.tag_id}.csv .
        tar -czf ${meta.tag_id}_nextclade.tar.gz ./${meta.tag_id}_nextclade
        """

    else
        """
        #!/bin/bash
        set -e

        nextclade run \
            -r ${ref_fasta} \
            -m ${ref_gff} \
            -O ${meta.tag_id}_nextclade \
            -s ${nextclade_output_verbosity} \
            -n ${meta.tag_id} \
            --include-reference true \
            ${input_fa}

        cp ./${meta.tag_id}_nextclade/${meta.tag_id}.csv .
        tar -czf ${meta.tag_id}_nextclade.tar.gz ./${meta.tag_id}_nextclade
        """
}