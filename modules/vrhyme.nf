// Coverage-table-driven vRhyme: take the per-sample subset of the pooled
// coverm metabat coverage TSV (rows are this sample's scaffolds, columns are
// avg/var pairs across all probe samples) and run vRhyme via -c. No -b/-r:
// the coverage signal is fully encoded in the TSV.
process VRHYME_WITH_COVERAGE {
    tag "${meta.id}"
    label 'cpu_8'
    label 'mem_16'
    label 'time_12'

    publishDir "${params.outdir}/${meta.id}/assembly/binning/vrhyme", mode: 'copy', overwrite: true, pattern: "vrhyme_out/log_vRhyme_*.log",                  saveAs: { f -> file(f).getName() }
    publishDir "${params.outdir}/${meta.id}/assembly/binning/vrhyme", mode: 'copy', overwrite: true, pattern: "vrhyme_out/vRhyme_best_bins.*.membership.tsv", saveAs: { f -> file(f).getName() }
    publishDir "${params.outdir}/${meta.id}/assembly/binning/vrhyme", mode: 'copy', overwrite: true, pattern: "vrhyme_out/vRhyme_best_bins.*.summary.tsv",    saveAs: { f -> file(f).getName() }
    publishDir "${params.outdir}/${meta.id}/assembly/binning/vrhyme", mode: 'copy', overwrite: true, pattern: "vrhyme_out/vRhyme_best_bins_fasta",            saveAs: { f -> file(f).getName() }

    container 'quay.io/biocontainers/vrhyme:1.1.0--pyhdfd78af_1'

    input:
    tuple val(meta), path(virus_fna), path(coverage_tsv)

    // log/summary stay optional: scaffolds may pass the upstream length filter
    // yet yield no bins (a normal outcome when binning signal is too sparse),
    // and nothing downstream needs a guaranteed real path for those two.
    // membership/bins_fasta are NOT optional: the script below guarantees a
    // real (possibly empty) file/dir for both even when vRhyme found no bins,
    // so every sample carries a real path — no optional-output handling needed
    // in consuming subworkflows (VCONTACT3_RUN, CHECKV_QC).
    output:
    tuple val(meta), path("vrhyme_out/log_vRhyme_*.log"),                  emit: log,        optional: true
    tuple val(meta), path("vrhyme_out/vRhyme_best_bins.*.membership.tsv"), emit: membership
    tuple val(meta), path("vrhyme_out/vRhyme_best_bins.*.summary.tsv"),    emit: summary,    optional: true
    tuple val(meta), path("vrhyme_out/vRhyme_best_bins_fasta"),            emit: bins_fasta
    val(meta),                                                             emit: done

    script:
    """
    vRhyme \\
        -i ${virus_fna} \\
        -c ${coverage_tsv} \\
        -t ${task.cpus} \\
        -o vrhyme_out

    # vRhyme writes neither the membership table nor the bins directory when
    # it finds no bins. Guarantee both exist (empty) so downstream consumers
    # never have to distinguish "no bins" from "no output".
    mkdir -p vrhyme_out/vRhyme_best_bins_fasta
    if ! ls vrhyme_out/vRhyme_best_bins.*.membership.tsv >/dev/null 2>&1; then
        printf 'scaffold\\tbin\\n' > vrhyme_out/vRhyme_best_bins.0.membership.tsv
    fi
    """
}

// Concatenate every qualified sample's virus.fna into one pooled fasta,
// prefixing each scaffold header with `<sampleId>||` so sample-of-origin
// survives bowtie/coverm and can be recovered when subsetting coverage back
// per sample for vRhyme.
process POOL_VIRAL_SCAFFOLDS {
    label 'cpu_1'
    label 'mem_4'
    label 'time_1'

    container 'quay.io/biocontainers/genomad:1.12.0--pyhdfd78af_0'

    input:
    tuple val(sample_ids), path(virus_fnas)
    path(pool_script)

    output:
    path "pooled_viral_scaffolds.fna", emit: pooled

    script:
    """
    python3 ${pool_script} \\
        --sample-ids ${sample_ids.join(' ')} \\
        --fastas ${virus_fnas.collect { it.name }.join(' ')} \\
        --output pooled_viral_scaffolds.fna
    """
}

// Filter a sorted BAM to only cross-sample-discriminating alignments: drops
// reads that mapped below `--min-read-percent-identity` or with less than
// `--min-read-aligned-length` aligned bases. Without this filter, a sample's
// reads spuriously stick to other samples' near-identical viral contigs and
// inflate cross-sample coverage covariance.
process COVERM_FILTER {
    tag "${meta.id}"
    label 'cpu_4'
    label 'mem_8'
    label 'time_1'

    container 'quay.io/biocontainers/coverm:0.7.0--hcb7b614_4'

    input:
    tuple val(meta), path(sorted_bam)

    output:
    tuple val(meta), path("${meta.id}.filtered.bam"), emit: bam

    script:
    """
    coverm filter \\
        --bam-files ${sorted_bam} \\
        --output-bam-files ${meta.id}.filtered.bam \\
        --min-read-aligned-length ${params.vrhyme_pool_min_aligned_length} \\
        --min-read-percent-identity ${params.vrhyme_pool_min_identity} \\
        --threads ${task.cpus}
    """
}

// Aggregate all filtered BAMs into one MetaBAT-style coverage TSV via coverm's
// `contig --methods metabat`. Each row is one (pooled, sample-prefixed) scaffold
// and there is one `<bam>` / `<bam>-var` column pair per input BAM. vRhyme
// accepts this format via `-c` after we drop contigLen/totalAvgDepth in the
// per-sample subset step.
process COVERM_DEPTH {
    label 'cpu_8'
    label 'mem_16'
    label 'time_12'

    container 'quay.io/biocontainers/coverm:0.7.0--hcb7b614_4'

    input:
    path(filtered_bams)

    output:
    path "master_coverage.tsv", emit: coverage

    script:
    """
    coverm contig \\
        --methods metabat \\
        --bam-files ${filtered_bams} \\
        --threads ${task.cpus} \\
        > master_coverage.tsv
    """
}

// Per-qualified-sample subset of the master coverage table: keep only rows
// whose scaffold name starts with `<sampleId>||`, strip that prefix back to
// the original scaffold name (so it matches geNomad's per-sample virus.fna),
// and drop coverm's contigLen + totalAvgDepth meta columns so the layout is
// the vRhyme-`-c`-compatible `contigName <bam1> <bam1>-var <bam2> ...` shape.
process SUBSET_COVERAGE_FOR_SAMPLE {
    tag "${meta.id}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    input:
    tuple val(meta), path(master_tsv)

    output:
    tuple val(meta), path("${meta.id}_vrhyme_coverage.tsv"), emit: coverage

    script:
    """
    awk -F'\\t' -v sid="${meta.id}||" '
        NR==1 {
            printf "%s", \$1
            for (i=4; i<=NF; i++) printf "\\t%s", \$i
            printf "\\n"
            next
        }
        index(\$1, sid) == 1 {
            scaffold = substr(\$1, length(sid)+1)
            printf "%s", scaffold
            for (i=4; i<=NF; i++) printf "\\t%s", \$i
            printf "\\n"
        }
    ' ${master_tsv} > ${meta.id}_vrhyme_coverage.tsv
    """
}

// Reclaim disk space from the per-sample BAMs (sorted + filtered) produced by
// the binning path. Only safe to run after COVERM_DEPTH has consumed them —
// the subworkflow gates this process on COVERM_DEPTH's output so timing is
// correct.
process CLEANUP_BOWTIE_BAMS_VRHYME {
    tag "${meta.id}"
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    input:
    tuple val(meta), path(bams), val(vrhyme_done)

    script:
    """
    # `path()` stages the BAMs as symlinks into this task's workdir; we resolve
    # them to their actual paths in the producing BOWTIE2SAMTOOLS workdirs and
    # delete those, freeing the real disk space. The symlinks themselves go
    # away when this task's workdir is cleaned up by Nextflow.
    for bam in ${bams}; do
        if [ -L "\$bam" ]; then
            real=\$(readlink -f "\$bam")
            rm -f "\$real"
        else
            rm -f "\$bam"
        fi
    done
    """
}


// Link scaffolds within each vRhyme bin into a single sequence per bin using vRhyme's
// link_bin_sequences.py helper. Each scaffold concatenation is interleaved with a run
// of Ns (default 1000) so downstream gene-callers won't predict ORFs across joins.
// Generic helper for any process downstream of vRhyme that needs one record per bin.
// See: https://github.com/AnantharamanLab/vRhyme/issues/4
process LINK_BIN_SCAFFOLDS {
    tag "${meta.id}"
    label 'cpu_1'
    label 'mem_4'
    label 'time_1'

    container 'quay.io/biocontainers/vrhyme:1.1.0--pyhdfd78af_1'

    input:
    tuple val(meta), path(bins_dir)

    output:
    tuple val(meta), path("${meta.id}_linked_bins.fna"), emit: linked

    script:
    """
    link_bin_sequences.py \\
        -i ${bins_dir} \\
        -o linked_bins_dir \\
        -e fasta \\
        -n ${params.vrhyme_link_n}
    cat linked_bins_dir/*.fasta > ${meta.id}_linked_bins.fna
    """
}
