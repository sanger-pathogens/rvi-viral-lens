process run_nextclade {

    tag "${meta.id}"
    label "nextclade"
    label "cpu_1"
    label "time_queue_from_small"

    input:
    tuple val(meta), path(input_fa), val(data_dir_list)

    output:
    tuple val(meta), path("*.json"), path("*.nextclade.tar.gz")

    script:
    def agg_label = "${meta.sample_id}.${meta.selected_taxid}"
    def cmd_lines = data_dir_list.collect { assembly_path ->
        def assembly_name = assembly_path.toString().split('/')[-1]

        def ref_fasta = "${assembly_path}/reference.fasta"

        def data_label = ""
        if (meta.flu_segment != "") {
            data_label = "${meta.sample_id}.${meta.selected_taxid}.${assembly_name}.segment${meta.flu_segment}"
        } else {
            data_label = "${meta.sample_id}.${meta.selected_taxid}.${assembly_name}"
        }


        """
        nextclade run \
            --input-dataset ${assembly_path} \
            -r ${ref_fasta} \
            -O ${data_label}.nextclade \
            -s "all" \
            -n ${data_label} \
            --include-reference true \
            ${input_fa}
        cp ${data_label}.nextclade/${data_label}.json .
        """.stripIndent()
        }


    """
    #!/bin/bash
    set -e

    ${cmd_lines.join('\n\n')}

    mkdir $agg_label && mv ./*.nextclade/ $agg_label
    mv $agg_label ${agg_label}.nextclade
    tar -czf ${agg_label}.nextclade.tar.gz ${agg_label}.nextclade
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