/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/

include { BOWTIE_INDEX; BOWTIE2SAMTOOLS                                                  } from '../modules/bowtie.nf'
include { VRHYME_WITH_COVERAGE; POOL_VIRAL_SCAFFOLDS; COVERM_FILTER; COVERM_DEPTH;
          SUBSET_COVERAGE_FOR_SAMPLE; CLEANUP_BOWTIE_BAMS_VRHYME                         } from '../modules/vrhyme.nf'

/*
========================================================================================
    HELPERS
========================================================================================
*/

// Count viral scaffolds in a geNomad virus_summary.tsv whose `length` exceeds a threshold.
// Returns 0 if the file is empty, has only a header, or lacks a `length` column.
def count_long_scaffolds(tsv, int min_length) {
    def lines = tsv.readLines()
    if (lines.size() < 2) return 0
    def header = lines[0].split('\t')
    def length_idx = header.findIndexOf { it == 'length' }
    if (length_idx < 0) return 0
    return lines[1..-1].count { line ->
        def cols = line.split('\t')
        length_idx < cols.size() && (cols[length_idx] as Integer) > min_length
    }
}

/*
========================================================================================
    VRHYME_BIN — ViWrap-style pooled cross-sample binning

    1. Gate samples on >=2 viral scaffolds longer than `vrhyme_min_scaffold_length`.
    2. Pool every qualified sample's virus.fna into one fasta, prefixing each
       scaffold header with `<sampleId>||` so sample-of-origin is recoverable.
    3. Build ONE bowtie2 index on the pooled fasta.
    4. Map every sample's reads against that single index (N x 1 bowtie tasks).
    5. Filter each BAM with `coverm filter` for identity + aligned-length cutoffs
       (drops cross-sample spurious mappings).
    6. Aggregate all filtered BAMs into one master coverage TSV via
       `coverm contig --methods metabat`.
    7. Per qualified sample: subset the master TSV by the `<sampleId>||` prefix
       (strip it back) and feed to vRhyme via -c.
    8. Cleanup: drop the per-sample BAMs once coverm has consumed them.
========================================================================================
*/

workflow VRHYME_BIN {

    take:
    virus_fna_ch     // tuple val(meta), path(virus_fna)         — from GENOMAD_CLASSIFY.out.virus_fna
    virus_summary_ch // tuple val(meta), path(virus_summary_tsv) — from GENOMAD_CLASSIFY.out.virus_summary
    reads_ch         // tuple val(meta), path(R1), path(R2)      — post-preprocessing reads (ALL samples — used as probes)

    main:

    // Gate: only samples with >=2 viral scaffolds long enough to bin pass.
    virus_fna_ch
        .join(virus_summary_ch)
        .map { meta, fna, summary ->
            def n_long = count_long_scaffolds(summary, params.vrhyme_min_scaffold_length)
            def verdict = n_long >= 2 ? 'pass' : 'skip'
            log.info "[vRhyme gate] ${meta.id}: ${n_long} scaffold(s) > ${params.vrhyme_min_scaffold_length}bp — ${verdict}"
            tuple(meta, fna, n_long)
        }
        .filter { meta, fna, n_long -> n_long >= 2 }
        .map { meta, fna, n_long -> tuple(meta, fna) }
        .set { ch_virus_qualified }

    // Pool qualified virus.fnas into one fasta with `<sampleId>||` headers.
    pool_script_ch = Channel.value(file("${projectDir}/bin/pool_viral_scaffolds.py"))

    ch_virus_qualified
        .toList()
        .map { items ->
            def ids  = items.collect { it[0].id }
            def fnas = items.collect { it[1] }
            tuple(ids, fnas)
        }
        .dump(tag: 'vrhyme_pool_input')
        .set { ch_pool_input }

    POOL_VIRAL_SCAFFOLDS(ch_pool_input, pool_script_ch)

    // One bowtie2 index on the pooled fasta. Synthetic meta `[id: 'pooled']`
    // so BOWTIE_INDEX names files `pooled_index*` for the downstream lookup.
    POOL_VIRAL_SCAFFOLDS.out.pooled
        .map { fna -> tuple([id: 'pooled'], fna) }
        .set { ch_pooled_index_input }

    BOWTIE_INDEX(ch_pooled_index_input)

    // Broadcast the single index to every sample's reads (N x 1 bowtie).
    reads_ch
        .combine(BOWTIE_INDEX.out.bowtie_index)
        .map { meta_reads, r1, r2, meta_index, bt2_files ->
            tuple(meta_reads, r1, r2, bt2_files, "${meta_index.id}_index")
        }
        .dump(tag: 'vrhyme_pool_bowtie_input')
        .set { ch_bowtie_input }

    BOWTIE2SAMTOOLS(ch_bowtie_input, params.vrhyme_bowtie_threads)

    // Identity / aligned-length filter on each sorted BAM.
    COVERM_FILTER(BOWTIE2SAMTOOLS.out.bam_file)

    // All filtered BAMs -> one master coverage TSV (rows = pooled scaffolds,
    // cols = per-sample avg/var pairs).
    COVERM_FILTER.out.bam
        .map { meta, bam -> bam }
        .collect()
        .dump(tag: 'vrhyme_pool_filtered_bams')
        .set { ch_all_filtered_bams }

    COVERM_DEPTH(ch_all_filtered_bams)

    // Per qualified sample: subset the master TSV to rows starting with
    // `<sampleId>||`, strip the prefix.
    ch_virus_qualified
        .map { meta, fna -> meta }
        .combine(COVERM_DEPTH.out.coverage)
        .set { ch_subset_input }

    SUBSET_COVERAGE_FOR_SAMPLE(ch_subset_input)

    // Pair the subset coverage TSV with the sample's virus.fna and run vRhyme via -c.
    ch_virus_qualified
        .join(SUBSET_COVERAGE_FOR_SAMPLE.out.coverage)
        .dump(tag: 'vrhyme_pool_vrhyme_input')
        .set { ch_vrhyme_input }

    VRHYME_WITH_COVERAGE(ch_vrhyme_input)

    // Cleanup: once coverm has produced master_coverage.tsv, the filtered and
    // sorted BAMs can be deleted. Use a keyed `.join()` to gate strictly on
    // COVERM_DEPTH having emitted (chained `.combine()` on value channels
    // carrying lists mis-pairs).
    all_bams_to_clean = BOWTIE2SAMTOOLS.out.bam_file
        .mix(COVERM_FILTER.out.bam)
        .map { meta, bam -> bam }
        .collect()
        .map { all_bams -> tuple('pooled', all_bams) }

    coverm_done = COVERM_DEPTH.out.coverage
        .map { coverage -> tuple('pooled', coverage) }

    all_bams_to_clean
        .join(coverm_done)
        .map { key, bams, coverage -> tuple([id: key], bams, 'done') }
        .set { ch_cleanup_input }

    CLEANUP_BOWTIE_BAMS_VRHYME(ch_cleanup_input)

    emit:
    log_file   = VRHYME_WITH_COVERAGE.out.log
    membership = VRHYME_WITH_COVERAGE.out.membership
    summary    = VRHYME_WITH_COVERAGE.out.summary
    bins_fasta = VRHYME_WITH_COVERAGE.out.bins_fasta
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
        -o ${meta.id}_linked_bins.fna \\
        -e fasta \\
        -n ${params.vrhyme_link_n}
    """
}
