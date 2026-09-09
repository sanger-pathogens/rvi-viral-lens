// Reference-sequence selection for one already-known species name: the "longest sequence
// per label" rule from reference_subset.nf's SELECT_REFERENCE_RECORDS, minus its
// abundance-threshold gating, for a caller that already knows which species it wants.
//
// UNUSED since the classifier/mapping split moved reference resolution into
// subworkflows/mapping.nf. That resolution now reuses the reference record the calling
// method's own map-QC table already names (EXTRACT_REFERENCE_RECORD in
// reference_subset.nf), so nothing needs to re-derive one from the species name. Kept, not
// deleted, for two reasons: it is the only deterministic species -> record rule in the
// codebase (seeded longest-sequence, independent of which method called the species),
// which is the fix if the multi-method non-determinism noted in
// subworkflows/classifying_index.nf ever matters; and it is the way to resolve a reference
// for a species that was never map-QC'd at all.

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
