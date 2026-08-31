// Reference-sequence selection for mSWEEP low-abundance hit validation.
// Given mSWEEP's per-sample relative-abundance calls, for every reference group above a
// threshold pick one (randomly, seeded for reproducibility) of its sequences in the
// positionally-aligned reference FASTA used to build the Themisto2 index
// (species_labels.txt line N == reference FASTA record N) — one mapping per species,
// not an aggregate across every sequence clustered under that label.

params.script_src_path = "${projectDir}/bin/"

process SELECT_REFERENCE_RECORDS {
    tag "${meta.id}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    container "quay.io/gsu-pipelines/rvi-vp-basecontainer"

    input:
    tuple val(meta), path(abundances)
    path(species_labels)
    path(sequence_lengths)

    output:
    tuple val(meta), path("${meta.id}_abundant_groups.tsv"), emit: groups, optional: true
    tuple val(meta), path("${meta.id}_record_ids.txt"), emit: record_ids, optional: true
    tuple val(meta), path("${meta.id}_index_label_map.tsv"), emit: index_label_map, optional: true

    script:
    """
    ${params.script_src_path}select_reference_records.py \\
        --abundances ${abundances} \\
        --species-labels ${species_labels} \\
        --sequence-lengths ${sequence_lengths} \\
        --min-abundance ${params.msweep_map_min_abundance} \\
        --seed ${params.msweep_map_reference_seed} \\
        --out-groups ${meta.id}_abundant_groups.tsv \\
        --out-record-ids ${meta.id}_record_ids.txt \\
        --out-index-label-map ${meta.id}_index_label_map.tsv
    """
}

process INDEX_REFERENCE_FASTA {
    label 'cpu_1'
    label 'mem_16'
    label 'time_queue_from_normal'

    // Runs once per pipeline run (not per-sample): tags every record in the reference
    // FASTA with its 1-based positional index so later per-sample extraction can match on
    // a stable exact token (SEQIDX_<n>) instead of fragile numeric record-slicing. Also
    // records each record's sequence length in the same pass, so SELECT_REFERENCE_RECORDS
    // can pick the longest sequence per species label.
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
