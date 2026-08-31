// Runs metagraph query (each mate independently, same per-mate pattern as
// metagraph_align.nf) directly against the pre-built de Bruijn graph + annotation,
// using --query-mode labels to print just the matched labels per query -- no alignment
// (no CIGAR/score/mismatch scoring), so this is a presence/pseudoalignment call rather
// than metagraph_align.nf's true alignment. Requested explicitly instead of reviving the
// old filter+query two-stage pipeline (git history: commit ccae756, "Switch METAGRAPH to
// true alignment instead of plain k-mer query") -- that approach's harsh unannotated
// presence-filter pre-stage discarded almost every real hit on real test data, and this
// module skips that stage entirely, querying the annotated index directly.
//
// Output shape, CONFIRMED against a real run: tab-separated, column 1 a per-query index,
// column 2 a query descriptor, and every column after that carrying the matched labels.
// The one thing the original assumption got wrong is the separator -- query mode joins a
// read's labels with ':' , not the ';' used elsewhere, e.g.
//   'OV040211.1 | Betacoronavirus pandemicum:OY127609.1 | Betacoronavirus pandemicum'
// On the first real run 6364 of 6367 rows were such multi-label fields, and treating each
// as a single label made every one its own "species". bin/call_metagraph_species.py splits
// on that separator (LABEL_JOIN_RE); it is careful not to split the ':' that introduces a
// coordinate window, so align-mode labels are unaffected.

process METAGRAPH_QUERY {
    tag "${meta.id}"
    label 'cpu_16'
    time { task.attempt > 1 ? ( task.previousTrace?.exit == 140 ? task.previousTrace.time * 2 : (task.previousTrace?.time ?: 12.h)) : (12.h) }

    container 'quay.io/biocontainers/metagraph:0.5.2--h5c1d0b5_0'

    memory { 25.GB * task.attempt }

    input:
    tuple val(meta), path(reads_1), path(reads_2)
    path(graph_dbg)
    path(annotation_annodbg)
    // See metagraph_align.nf's identical input for why this sidecar is needed.
    path(annotation_seqs)

    output:
    tuple val(meta), path("${meta.id}_metagraph_query_1.tsv"), path("${meta.id}_metagraph_query_2.tsv"), emit: alignments

    script:
    query_params = "-i ${graph_dbg} -a ${annotation_annodbg} -p ${task.cpus} " +
        "--query-mode labels " +
        "--min-kmers-fraction-label ${params.metagraph_align_min_kmers_fraction_label}"

    """
    bash -c 'metagraph query ${query_params} ${reads_1} > ${meta.id}_metagraph_query_1.tsv'
    bash -c 'metagraph query ${query_params} ${reads_2} > ${meta.id}_metagraph_query_2.tsv'
    """
}
