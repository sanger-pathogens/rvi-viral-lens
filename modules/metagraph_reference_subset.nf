// Reference-sequence extraction for metagraph-align species-call validation.
// Unlike reference_subset.nf's EXTRACT_REFERENCE_SUBSET (mSWEEP's own mapping validation),
// this can't do one exact-ID grep -f pass against a positionally-indexed reference: several
// distinct reference sequences in the flu_rbfish library can legitimately share one taxid
// (multiple strains/segments under the same species), so each called species' record_ids.txt
// pattern (see call_metagraph_species.py) is grep'd separately and capped at its first
// match, keeping the same "one mapping per species" contract mSWEEP's validation uses.
// Headers are then rewritten to just their record_id (the awk pass below) so the samtools
// rname downstream exactly matches metagraph_species_call.nf's index-label-map keys,
// regardless of how a header's original trailing description text reads. Plain awk rather
// than a bin/ script since this process's container is seqkit-only (no python3).

process EXTRACT_METAGRAPH_REFERENCE_SUBSET {
    tag "${meta.id}"
    label 'cpu_2'
    label 'mem_4'
    label 'time_queue_from_normal'

    container 'quay.io/biocontainers/seqkit:2.10.0--h9ee0642_0'

    input:
    tuple val(meta), path(record_ids)
    path(reference_fasta)

    output:
    tuple val(meta), path("${meta.id}_subset.fasta"), emit: subset_fasta, optional: true

    script:
    """
    : > ${meta.id}_subset_raw.fasta
    while IFS= read -r pattern; do
        # `seqkit head -n 1` closing its stdin after the first record can SIGPIPE the grep
        # still writing further matches (some called taxids match many library records) --
        # under this script's `set -o pipefail` that would abort the whole task even though
        # the pipe did exactly what we wanted (one record kept). The `|| true` swallows
        # grep's side of that race; head's own exit code is what pipefail actually checks.
        (seqkit grep -n -r -p "\${pattern}" ${reference_fasta} || true) | seqkit head -n 1 >> ${meta.id}_subset_raw.fasta
    done < ${record_ids}

    awk '
        /^>/{
            header = substr(\$0, 2)
            if (match(header, /^kraken:taxid\\|[0-9]+\\|/)) {
                taxid_part = substr(header, RSTART, RLENGTH)
                gsub(/^kraken:taxid\\|/, "", taxid_part)
                gsub(/\\|\$/, "", taxid_part)
                print ">" taxid_part
            } else {
                split(header, arr, " ")
                print ">" arr[1]
            }
            next
        }
        { print }
    ' ${meta.id}_subset_raw.fasta > ${meta.id}_subset.fasta
    """
}
