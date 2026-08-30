// Shared per-sample / run-level report writers for the mapping/assembly/abundance
// lanes. Mirrors modules/write_summary_report.nf's two-step shape (per-sample dump,
// then run-level concatenation), but simpler: these lanes have no separate qc/nextclade
// JSON to merge in, so the per-sample step is a plain `meta` dump with no python
// involved, and the run-level step uses the shared bin/write_lane_summary.py rather
// than a lane-specific script (see that script's docstring for why).

process write_lane_sequence_summary {
    tag "${sample_id}"

    input:
    tuple val(sample_id), val(meta)

    output:
    tuple val(meta), path("${sample_id}.properties.json")

    script:
    def metaJson = groovy.json.JsonOutput.toJson(meta)
    def formatted = groovy.json.JsonOutput.prettyPrint(metaJson)
    """
    echo '${formatted}' > ${sample_id}.properties.json
    """
}

process write_lane_run_summary {
    input:
    val(json_file_list)
    val(lane_name)

    output:
    tuple path("${lane_name}_run_summary.json"), path("${lane_name}_summary_report.csv")

    script:
    json_arg = json_file_list.join(" ")
    """
    write_lane_summary.py -c ${json_arg} -o ${lane_name}
    """
}
