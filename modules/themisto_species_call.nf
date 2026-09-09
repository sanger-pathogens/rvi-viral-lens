// Provisional species calling straight from Themisto2 pseudoalignment counts: counts
// reads per species (deduped per read across any multi-mapped reference indices of the
// same species) from both mates' pseudoalignment output, calls a species present once
// its hit count clears themisto_align_min_hits, and — for each called species — records
// its single most-hit reference index (see call_themisto_species.py). Unlike mSWEEP,
// this needs no probabilistic abundance model: species_labels.txt already maps each
// 0-based Themisto reference index straight to a label, so the read-hit count itself is
// enough to call a species present and to pick which one of its sequences to validate.
// Reference indices are also 1:1 positional with the msweep_map_reference_fasta/
// INDEX_REFERENCE_FASTA SEQIDX_<n> tagging (reference_subset.nf), so no separate
// reference-extraction module is needed downstream (see themisto_map_qc.nf).

// VIRAL-LENS DEVIATION from upstream: bin/call_themisto_species.py skips placeholder
// labels (currently the literal "NA", 7500 lines of rvdb_clustered_virome_species_labels.txt)
// instead of treating them as a species. See UNUSABLE_LABELS in that script for why.

params.script_src_path = "${projectDir}/bin/"

process CALL_THEMISTO_SPECIES {
    tag "${meta.id}"
    label 'cpu_1'
    label 'mem_4'
    // Single-threaded pure-Python parse over --themisto1-output-format output: on a
    // heavily-clustered viral index (e.g. rvdb_clustered_virome) one read can still
    // pseudoalign against thousands of near-identical reference indices, so this scales
    // with pseudoalignment fan-out, not just read count — time_1 undersells it badly.
    label 'time_queue_from_normal'

    container "quay.io/gsu-pipelines/rvi-vp-basecontainer"

    publishDir "${params.outdir}/${meta.id}/sequenceindex/themisto_hits", mode: 'copy', overwrite: true, pattern: "*_species_hits.tsv"

    input:
    tuple val(meta), path(pseudoalignment_1), path(pseudoalignment_2)
    path(species_labels)

    output:
    tuple val(meta), path("${meta.id}_species_hits.tsv"), emit: species_hits
    tuple val(meta), path("${meta.id}_record_ids.txt"), emit: record_ids, optional: true
    tuple val(meta), path("${meta.id}_index_label_map.tsv"), emit: index_label_map, optional: true
    tuple val(meta), path(pseudoalignment_1), path(pseudoalignment_2), emit: pseudoalignments

    script:
    """
    ${params.script_src_path}call_themisto_species.py \\
        --pseudoalignment-1 ${pseudoalignment_1} \\
        --pseudoalignment-2 ${pseudoalignment_2} \\
        --species-labels ${species_labels} \\
        --min-hits ${params.themisto_align_min_hits} \\
        --sample-id ${meta.id} \\
        --out-species-hits ${meta.id}_species_hits.tsv \\
        --out-record-ids ${meta.id}_record_ids.txt \\
        --out-index-label-map ${meta.id}_index_label_map.tsv
    """
}
