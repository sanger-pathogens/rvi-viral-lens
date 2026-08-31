#!/usr/bin/env nextflow

/*
========================================================================================
    mSWEEP low-abundance hit validation sub-workflow
    ------------------------
    For every reference group mSWEEP calls above msweep_map_min_abundance, map the
    sample's own preprocessed reads directly against that species' longest reference
    sequence (pulled positionally from a reference FASTA aligned to species_labels.txt)
    and record genome breadth of coverage, mean depth, mapping/base quality, and reads
    mapped — one mapping per species. This is a sanity check on mSWEEP's probabilistic
    calls — independent of, and not the same as, mGEMS binning.
    Per-sample results are published under <outdir>/<sample>/sequenceindex/msweep_map; a
    run-wide summary is published under <outdir>/msweep_map_summary.
========================================================================================
*/

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/
include { SELECT_REFERENCE_RECORDS
          INDEX_REFERENCE_FASTA
          EXTRACT_REFERENCE_SUBSET       } from '../modules/reference_subset.nf'
include { BOWTIE_INDEX; BOWTIE2SAMTOOLS  } from '../modules/bowtie.nf'
include { SAMTOOLS_COVERAGE
          AGGREGATE_SPECIES_COVERAGE
          GENERATE_MSWEEP_MAP_SUMMARY    } from '../modules/samtools_coverage.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow MSWEEP_MAP_QC {

    take:
    reads_ch           // tuple( meta, read_1, read_2 ) — same preprocessed reads VIRAL_MSWEEP receives
    abundances_ch      // MSWEEP.out.abundances: tuple( meta, mSWEEP_abundances.txt, mSWEEP_probs.tsv )
    species_labels_ch  // value channel: species_labels.txt (same file mSWEEP used as ref_groups)

    main:
    // Reference FASTA whose Nth record corresponds to the Nth line of species_labels.txt
    // (the same order the .thm2 index was built from). Tagged once with positional
    // SEQIDX_<n> IDs (+ per-record lengths, for picking the longest sequence per species
    // below) and reused as value channels across every sample.
    INDEX_REFERENCE_FASTA(Channel.fromPath(params.msweep_map_reference_fasta))
    indexed_reference_ch  = INDEX_REFERENCE_FASTA.out.fasta.first()
    sequence_lengths_ch   = INDEX_REFERENCE_FASTA.out.lengths.first()

    abundances_only_ch = abundances_ch.map { meta, abundances, _probs -> tuple(meta, abundances) }

    SELECT_REFERENCE_RECORDS(abundances_only_ch, species_labels_ch, sequence_lengths_ch)

    // Samples with nothing above msweep_map_min_abundance emit nothing here (optional
    // outputs) and are naturally dropped by every join() below — no branching needed.
    EXTRACT_REFERENCE_SUBSET(SELECT_REFERENCE_RECORDS.out.record_ids, indexed_reference_ch)

    BOWTIE_INDEX(EXTRACT_REFERENCE_SUBSET.out.subset_fasta)

    // BOWTIE2SAMTOOLS wants (meta, r1, r2, bt2_files, index_prefix); BOWTIE_INDEX names its
    // index files "${meta.id}_index*", so the prefix string is derived rather than tracked
    // through the channel.
    mapping_input_ch = reads_ch
        .join(BOWTIE_INDEX.out.bowtie_index)
        .map { meta, r1, r2, bt2_files -> tuple(meta, r1, r2, bt2_files, "${meta.id}_index") }

    BOWTIE2SAMTOOLS(mapping_input_ch, params.msweep_map_bowtie_threads)

    SAMTOOLS_COVERAGE(BOWTIE2SAMTOOLS.out.bam_file)

    qc_input_ch = SAMTOOLS_COVERAGE.out.coverage
        .join(SAMTOOLS_COVERAGE.out.query_lengths)
        .join(SELECT_REFERENCE_RECORDS.out.index_label_map)
        .join(SELECT_REFERENCE_RECORDS.out.groups)

    AGGREGATE_SPECIES_COVERAGE(qc_input_ch)

    GENERATE_MSWEEP_MAP_SUMMARY(
        AGGREGATE_SPECIES_COVERAGE.out.qc_table.map { _meta, qc_table -> qc_table }.collect()
    )

    emit:
    qc_table = AGGREGATE_SPECIES_COVERAGE.out.qc_table
}

/*
========================================================================================
    THE END
========================================================================================
*/
