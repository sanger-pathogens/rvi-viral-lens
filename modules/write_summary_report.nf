process write_sequence_level_summaries {
    input:
    tuple val(sample_id), val(meta), val(qc_json), val(nc_json)

    output:
    path("${sample_id}.properties.json")

    script:
    def metaJson = groovy.json.JsonOutput.toJson(meta)
    def formatted = groovy.json.JsonOutput.prettyPrint(metaJson)
    if (nc_json == null )
        """
        echo '${formatted}' > ${sample_id}.meta.json
        write_all_summaries.py -rm "sample" -i "${sample_id}" -m  ${sample_id}.meta.json -p "${qc_json}"
        """
    else
        """
        echo '${formatted}' > ${sample_id}.meta.json
        write_all_summaries.py -rm "sample" -i "${sample_id}" -m ${sample_id}.meta.json -p "${qc_json}" -n "${nc_json}"
        """
}


process write_run_level_summaries {
    input:
    val(json_file_list)

    output:
    tuple path("consensus_sequence_properties.json"), path("summary_report.csv")

    script:
    json_arg = json_file_list.join(" ")
    """
    write_all_summaries.py -rm "run" -c ${json_arg} -r 1
    """
}