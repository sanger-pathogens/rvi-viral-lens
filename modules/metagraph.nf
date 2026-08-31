// Runs metagraph align (each mate independently, see call_metagraph_species.py) directly
// against the pre-built de Bruijn graph + annotation, label/trace-consistent alignment on
// (mirrors themisto2.nf's pattern for paired reads: both mates run once each, rather than a
// single joint-mate invocation). --min-kmers-fraction-label / --align-min-exact-match
// require true alignment tolerant of mismatches to find a match at all — a plain
// query-presence prefilter (this module's earlier approach) discards any read that isn't a
// 100%-identical k-mer path through the graph, which throws away most real hits before they
// ever reach label matching. Output shape (no header): col 1 the read id, col 2 the query
// sequence, then a variable number of 7-field alternative-alignment groups (strand, aligned
// reference sequence, score, aligned length, CIGAR, mismatch count, ';'-joined matched
// labels) for reads that aligned, or a single '*' placeholder group (no labels field at
// all) for reads that didn't — see call_metagraph_species.py, which parses this exact shape.

process METAGRAPH {
    tag "${meta.id}"
    label 'cpu_16'
    time { task.attempt > 1 ? ( task.previousTrace?.exit == 140 ? task.previousTrace.time * 2 : (task.previousTrace?.time ?: 12.h)) : (12.h) }

    container 'quay.io/biocontainers/metagraph:0.5.2--h5c1d0b5_0'

    memory { 25.GB * task.attempt }

    input:
    tuple val(meta), path(reads_1), path(reads_2)
    path(graph_dbg)
    path(annotation_annodbg)
    // metagraph derives this sidecar's name from annotation_annodbg's basename and looks
    // for it in the working directory; without it, alignment falls back to bogus
    // 'noid:<coord-ranges>' labels instead of real ones (see metagraph_align_annotation_seqs
    // in metagraph_align.config) -- only relevant for coordinate-type annotations.
    path(annotation_seqs)

    output:
    tuple val(meta), path("${meta.id}_metagraph_align_1.tsv"), path("${meta.id}_metagraph_align_2.tsv"), emit: alignments

    script:
    align_params = "-i ${graph_dbg} -a ${annotation_annodbg} -p ${task.cpus} " +
        "--min-kmers-fraction-label ${params.metagraph_align_min_kmers_fraction_label} " +
        "--align-min-exact-match ${params.metagraph_align_min_exact_match}"

    """
    bash -c 'metagraph align ${align_params} -o ${meta.id}_metagraph_align_1.tsv ${reads_1}'
    bash -c 'metagraph align ${align_params} -o ${meta.id}_metagraph_align_2.tsv ${reads_2}'
    """
}
