// Provisional species calling from metagraph align output: counts reads per species
// (deduped per read across any multi-mapped records of the same species) from both
// mates' alignment output, calls a species present once its hit count clears
// metagraph_align_min_hits, and — for each called species — records its single most-hit
// reference record (a taxid, or a bare accession for the few non-taxid labels — see
// call_metagraph_species.py). Species keyed by taxid get relabeled with their names_dmp
// scientific name (strain-level, not species-level, for this index) in place of the raw
// alignment label. Unlike mSWEEP's mapping validation, no positional species-labels
// indirection is needed here: that record can be extracted straight from the reference
// FASTA in metagraph_map_qc.nf.

params.script_src_path = "${projectDir}/bin/"

process CALL_METAGRAPH_SPECIES {
    tag "${meta.id}"
    label 'cpu_1'
    label 'mem_2'
    label 'time_1'

    container "quay.io/gsu-pipelines/rvi-vp-basecontainer"

    // output_subdir lets VIRAL_METAGRAPH_ALIGN.nf ('metagraph_hits') and
    // VIRAL_METAGRAPH_QUERY.nf ('metagraph_query_hits') share this process without
    // overwriting each other's output when both methods run for the same sample.
    publishDir "${params.outdir}/${meta.id}/sequenceindex/${output_subdir}", mode: 'copy', overwrite: true, pattern: "*_species_hits.tsv"

    input:
    tuple val(meta), path(align_1), path(align_2)
    path(names_dmp)
    val(output_subdir)

    output:
    tuple val(meta), path("${meta.id}_species_hits.tsv"), emit: species_hits
    tuple val(meta), path("${meta.id}_record_ids.txt"), emit: record_ids, optional: true
    tuple val(meta), path("${meta.id}_index_label_map.tsv"), emit: index_label_map, optional: true

    script:
    """
    ${params.script_src_path}call_metagraph_species.py \\
        --align-1 ${align_1} \\
        --align-2 ${align_2} \\
        --min-hits ${params.metagraph_align_min_hits} \\
        --sample-id ${meta.id} \\
        --names-dmp ${names_dmp} \\
        --out-species-hits ${meta.id}_species_hits.tsv \\
        --out-record-ids ${meta.id}_record_ids.txt \\
        --out-index-label-map ${meta.id}_index_label_map.tsv

    # align_1/align_2 are METAGRAPH's raw per-mate query output (can reach tens
    # to hundreds of GB for large samples) and this is their only consumer downstream
    # of METAGRAPH. Free the real bytes in METAGRAPH's own work dir now
    # rather than leaving them until the whole run finishes — readlink resolves through
    # the symlinked stage-in so this removes the source file, not just the local link.
    # Trade-off: a future `-resume` that needs to re-run this task specifically (e.g.
    # after changing this script or --min-hits) will have to re-run METAGRAPH
    # too, since its cached output will be gone.
    rm -f "\$(readlink -f ${align_1})" "\$(readlink -f ${align_2})"
    """
}
