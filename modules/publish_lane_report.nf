// Per-sample publish step for the mapping/assembly/abundance lanes' properties.json.
// Can't reuse modules/publish_lite.nf's publish_consensus_files as-is: its
// `consensus_output` label publishes under `${meta.sample_id}/${meta.taxid}/`, and
// these lanes' meta has no `.taxid`. `lane_output` publishes under `${meta.sample_id}/`
// only (see nextflow.config).

process publish_lane_json {
    label "lane_output"

    input:
        tuple val(meta), path(file)

    output:
        tuple val(meta), path(file)

    script:
    """
    echo "Published: ${file}"
    """
}
