// Reference-sequence selection for one already-known species name (used by the
// sequence-index lane's new-species-consensus feature, subworkflows/sequence_index.nf,
// to resolve a reference for a species mSWEEP/Metagraph called but Kraken2 didn't).
// Same "longest sequence per label" rule as reference_subset.nf's
// SELECT_REFERENCE_RECORDS, minus its abundance-threshold gating -- the caller already
// knows exactly which species it wants.

params.script_src_path = "${projectDir}/bin/"

process SELECT_REFERENCE_RECORD_BY_NAME {
    tag "${meta.id}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    container "quay.io/gsu-pipelines/rvi-vp-basecontainer"

    input:
    tuple val(meta), val(species_name)
    path(species_labels)
    path(sequence_lengths)

    output:
    tuple val(meta), path("${meta.id}_record_id.txt"), emit: record_id, optional: true

    script:
    """
    ${params.script_src_path}select_reference_record_by_name.py \\
        --species-name "${species_name}" \\
        --species-labels ${species_labels} \\
        --sequence-lengths ${sequence_lengths} \\
        --seed ${params.msweep_map_reference_seed} \\
        --out-record-id ${meta.id}_record_id.txt
    """
}
