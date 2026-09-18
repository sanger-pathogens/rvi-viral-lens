// Reference-record selection and extraction out of a positionally-indexed reference FASTA.
// INDEX_REFERENCE_FASTA tags every record with its 1-based position (SEQIDX_<n>) and
// records its length in the same pass; the EXTRACT_* processes pull records back out by
// that tag. The abundance-driven selection this file used to open with
// (SELECT_REFERENCE_RECORDS) went with mSWEEP map-QC, which nothing invoked.

params.script_src_path = "${projectDir}/bin/"

process INDEX_REFERENCE_FASTA {
    label 'cpu_1'
    label 'mem_16'
    label 'time_queue_from_normal'

    // Runs once per pipeline run (not per-sample): tags every record in the reference
    // FASTA with its 1-based positional index so later per-sample extraction can match on
    // a stable exact token (SEQIDX_<n>) instead of fragile numeric record-slicing.
    //
    // The `lengths` emit is a by-product of the same pass. Its only consumer was
    // SELECT_REFERENCE_RECORDS, and the sequence-index callers get record lengths from
    // INDEX_REFERENCE_LENGTHS instead -- which streams the FASTA without writing a
    // tagged copy of it. Left emitted rather than removed: it costs nothing here, and
    // splitting the two passes is the caller's choice, not this process's.
    stageInMode 'symlink'

    input:
    path(reference_fasta)

    output:
    path("indexed_reference.fasta"), emit: fasta
    path("sequence_lengths.tsv"), emit: lengths

    script:
    """
    awk '
        /^>/{
            if (n > 0) print n"\\t"len > "sequence_lengths.tsv"
            n++; len = 0
            print ">SEQIDX_" n " " substr(\$0,2) > "indexed_reference.fasta"
            next
        }
        { len += length(\$0); print > "indexed_reference.fasta" }
        END{ if (n > 0) print n"\\t"len > "sequence_lengths.tsv" }
    ' ${reference_fasta}
    """
}

// Single-record variant of EXTRACT_REFERENCE_SUBSET, taking the record id as a VALUE
// rather than a file of ids. Used by subworkflows/mapping.nf for sequence-index species:
// the calling method's own index-label map already names the one reference record it
// picked for that species (a SEQIDX_<n> token, same id space INDEX_REFERENCE_FASTA
// mints), so there is nothing to write to a file first -- seqkit grep -p matches the
// sequence ID directly.
process EXTRACT_REFERENCE_RECORD {
    tag "${meta.id}"
    label 'cpu_2'
    label 'mem_4'
    label 'time_queue_from_normal'

    container 'quay.io/biocontainers/seqkit:2.10.0--h9ee0642_0'

    input:
    tuple val(meta), val(record_id)
    path(indexed_reference_fasta)

    output:
    tuple val(meta), path("${meta.id}_subset.fasta"), emit: subset_fasta, optional: true

    script:
    """
    seqkit grep -p "${record_id}" ${indexed_reference_fasta} > ${meta.id}_subset.fasta

    # An id that matched nothing still leaves a 0-byte FASTA behind, since `>` creates the
    # file either way. Remove it so the optional output is genuinely absent and this
    # species is dropped, rather than reaching GENERATE_CONSENSUS with no reference.
    if [ ! -s ${meta.id}_subset.fasta ]; then
        rm -f ${meta.id}_subset.fasta
    fi
    """
}

process EXTRACT_REFERENCE_SUBSET {
    tag "${meta.id}"
    label 'cpu_2'
    label 'mem_4'
    label 'time_queue_from_normal'

    container 'quay.io/biocontainers/seqkit:2.10.0--h9ee0642_0'

    input:
    tuple val(meta), path(record_ids)
    path(indexed_reference_fasta)

    output:
    tuple val(meta), path("${meta.id}_subset.fasta"), emit: subset_fasta, optional: true

    script:
    """
    seqkit grep -f ${record_ids} ${indexed_reference_fasta} > ${meta.id}_subset.fasta
    """
}
