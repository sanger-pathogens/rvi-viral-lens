process run_nextclade {

    tag "${meta.id}"
    label "nextclade"

    input:
    tuple val(meta), path(input_fa), val(data_dir_list)

    output:
    tuple val(meta), path("*.json"), path("*_nextclade.tar.gz")

    script:
    def agg_label = "${meta.sample_id}.${meta.selected_taxid}"
    def cmd_lines = data_dir_list.collect { assembly_path ->
        def assembly_name = assembly_path.toString().split('/')[-1]

        def ref_fasta = "${assembly_path}/reference.fasta"
        def ref_tree = "${assembly_path}/tree.json"
        def ref_gff = "${assembly_path}/genome_annotation.CDS.gff3"

        def data_label = ""
        if (meta.flu_segment != "") {
            data_label = "${meta.sample_id}.${assembly_name}.segment${meta.flu_segment}"
        } else {
            data_label = "${meta.sample_id}.${assembly_name}"
        }

        if (file(ref_tree).exists()) {
            """
            nextclade run \
                -r ${ref_fasta} \
                -a ${ref_tree} \
                -m ${ref_gff} \
                -O ${data_label}_nextclade \
                -s "all" \
                -n ${data_label} \
                --include-reference true \
                ${input_fa}
            cp ${data_label}_nextclade/${data_label}.json .
            """.stripIndent()

        } else {
            """
            nextclade run \
                -r ${ref_fasta} \
                -m ${ref_gff} \
                -O ${data_label}_nextclade \
                -s "all" \
                -n ${data_label} \
                --include-reference true \
                ${input_fa}
            cp ${data_label}_nextclade/${data_label}.json .
            """.stripIndent()
        }
    }

    """
    #!/bin/bash
    set -e

    ${cmd_lines.join('\n\n')}

    mkdir $agg_label && mv ./*_nextclade/ $agg_label
    tar -czf ${agg_label}_nextclade.tar.gz $agg_label
    """
}

process collate_nextclade_jsons {
    tag "${meta.id}"

    input:
    tuple val(meta), val(json_path_list), path(nextclade_tarball)

    output:
    tuple val(meta), path("${meta.sample_id}.${meta.selected_taxid}.nextclade.json"), path(nextclade_tarball)

    script:
    def json_file_args = ([json_path_list].flatten()*.toString()).join(" ")
    """
    aggregate_nextclade_json.py ${json_file_args} -o ${meta.sample_id}.${meta.selected_taxid}.nextclade.json
    """
}