#!/usr/bin/env nextflow

/*
========================================================================================
    UNUSED (rvi_integration_1). Nothing invokes this subworkflow: sequence-index species
    are called on read-hit counts alone, and the breadth figure this used to produce is
    now measured downstream from the consensus alignment subworkflows/mapping.nf performs
    anyway (params.new_species_min_breadth_pct is applied there). Keeping a validation
    mapping here meant every real call was mapped twice, by two different aligners.

    Kept, not deleted, because it is the only thing that can measure breadth for a species
    BEFORE deciding to spend a consensus on it. Restore it if that ordering ever matters --
    e.g. if index noise volume makes the discarded consensuses expensive, or if breadth is
    wanted for calls that never get a consensus at all (ones Kraken2 already found). Its
    params are still defined in nextflow.config, marked UNUSED.
========================================================================================
*/

/*
========================================================================================
    Themisto-hit species-call validation sub-workflow
    ------------------------
    For every species CALL_THEMISTO_SPECIES calls above themisto_align_min_hits, map the
    sample's own preprocessed reads directly against that species' single most-hit
    reference sequence (pulled positionally from a reference FASTA aligned to
    species_labels.txt, tagged SEQIDX_<n> by INDEX_REFERENCE_FASTA — same convention
    mSWEEP's own mapping validation uses, see msweep_map_qc.nf) and record genome breadth
    of coverage, mean depth, mapping/base quality, and reads mapped — one mapping per
    species. This is a sanity check on Themisto's own pseudoalignment read-hit calls,
    independent of mSWEEP's probabilistic abundance estimation.
    Per-sample results are published under <results_dir>/<sample>/themisto_map; a run-wide
    summary is published under <results_dir>/themisto_map_summary.
    PORTED (rvi_integration_1) from eu1/rvi_toolbox.git's `feature_msweep_map` branch
    (`subworkflows/themisto_map_qc.nf`), as a viral-lens-owned file rather than a submodule
    bump, same as every other lane in this integration -- see INSTRUCT.md item 2 (the
    unresolved fork question). Differences from upstream are confined to paths:
    params.results_dir -> params.outdir, publish under <outdir>/<sample>/sequenceindex/
    to match the other sequence-index methods, bin/ instead of rvi_toolbox/bin/, and
    include paths for viral-lens's workflows/ + modules/ layout. If the upstream branch
    changes, re-apply by hand.
========================================================================================
*/

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/
include { INDEX_REFERENCE_FASTA
          EXTRACT_REFERENCE_SUBSET     } from '../modules/reference_subset.nf'
include { BOWTIE_INDEX; BOWTIE2SAMTOOLS  } from '../modules/bowtie.nf'
include { SAMTOOLS_COVERAGE               } from '../modules/samtools_coverage.nf'
include { AGGREGATE_THEMISTO_COVERAGE
          GENERATE_THEMISTO_MAP_SUMMARY } from '../modules/themisto_coverage.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow THEMISTO_MAP_QC {

    take:
    reads_ch            // tuple( meta, read_1, read_2 ) — same (capped) reads CALL_THEMISTO_SPECIES received
    record_ids_ch        // CALL_THEMISTO_SPECIES.out.record_ids: tuple( meta, record_ids.txt ), optional per sample
    index_label_map_ch  // CALL_THEMISTO_SPECIES.out.index_label_map: tuple( meta, index_label_map.tsv ), optional per sample
    species_hits_ch      // CALL_THEMISTO_SPECIES.out.species_hits: tuple( meta, species_hits.tsv )
    reference_fasta_ch  // value channel: reference multi-FASTA whose Nth record == Nth line of species_labels.txt

    main:
    // Reference FASTA whose Nth record corresponds to the Nth line of species_labels.txt
    // (the same order the .thm2 index was built from) — tagged once with positional
    // SEQIDX_<n> IDs and reused as a value channel across every sample.
    INDEX_REFERENCE_FASTA(reference_fasta_ch)
    indexed_reference_ch = INDEX_REFERENCE_FASTA.out.fasta.first()

    // Samples with nothing above themisto_align_min_hits emit nothing here (optional
    // outputs) and are naturally dropped by every join() below — no branching needed.
    EXTRACT_REFERENCE_SUBSET(record_ids_ch, indexed_reference_ch)

    BOWTIE_INDEX(EXTRACT_REFERENCE_SUBSET.out.subset_fasta)

    // BOWTIE2SAMTOOLS wants (meta, r1, r2, bt2_files, index_prefix); BOWTIE_INDEX names its
    // index files "${meta.id}_index*", so the prefix string is derived rather than tracked
    // through the channel.
    mapping_input_ch = reads_ch
        .join(BOWTIE_INDEX.out.bowtie_index)
        .map { meta, r1, r2, bt2_files -> tuple(meta, r1, r2, bt2_files, "${meta.id}_index") }

    BOWTIE2SAMTOOLS(mapping_input_ch, params.themisto_map_bowtie_threads)

    SAMTOOLS_COVERAGE(BOWTIE2SAMTOOLS.out.bam_file)

    qc_input_ch = SAMTOOLS_COVERAGE.out.coverage
        .join(SAMTOOLS_COVERAGE.out.query_lengths)
        .join(index_label_map_ch)
        .join(species_hits_ch)

    AGGREGATE_THEMISTO_COVERAGE(qc_input_ch)

    GENERATE_THEMISTO_MAP_SUMMARY(
        AGGREGATE_THEMISTO_COVERAGE.out.qc_table.map { _meta, qc_table -> qc_table }.collect()
    )

    emit:
    qc_table = AGGREGATE_THEMISTO_COVERAGE.out.qc_table
}

/*
========================================================================================
    THE END
========================================================================================
*/
