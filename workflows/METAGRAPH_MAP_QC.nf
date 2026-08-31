#!/usr/bin/env nextflow

/*
========================================================================================
    metagraph-align species-call validation sub-workflow
    ------------------------
    For every species CALL_METAGRAPH_SPECIES calls above metagraph_align_min_hits, map the
    sample's own preprocessed reads directly against that species' single most-hit
    reference record (a taxid or an accession, see call_metagraph_species.py) — pulled
    straight from metagraph's own alignment labels, no positional species-labels
    indirection needed (unlike mSWEEP's mapping validation, see msweep_map_qc.nf) — and
    record genome breadth of coverage, mean depth, mapping/base quality, and reads mapped:
    one mapping per species. This is a sanity check on metagraph align's read-hit calls.
    Per-sample results are published under <outdir>/<sample>/sequenceindex/metagraph_map; a
    run-wide summary is published under <outdir>/metagraph_map_summary.
========================================================================================
*/

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/
include { EXTRACT_METAGRAPH_REFERENCE_SUBSET                            } from '../modules/metagraph_reference_subset.nf'
include { BOWTIE_INDEX; BOWTIE2SAMTOOLS                                 } from '../modules/bowtie.nf'
include { SAMTOOLS_COVERAGE                                             } from '../modules/samtools_coverage.nf'
include { AGGREGATE_METAGRAPH_COVERAGE; GENERATE_METAGRAPH_MAP_SUMMARY  } from '../modules/metagraph_coverage.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow METAGRAPH_MAP_QC {

    take:
    reads_ch            // tuple( meta, read_1, read_2 ) — same preprocessed reads METAGRAPH received
    record_ids_ch       // CALL_METAGRAPH_SPECIES.out.record_ids: tuple( meta, record_ids.txt ), optional per sample
    index_label_map_ch  // CALL_METAGRAPH_SPECIES.out.index_label_map: tuple( meta, index_label_map.tsv ), optional per sample
    species_hits_ch     // CALL_METAGRAPH_SPECIES.out.species_hits: tuple( meta, species_hits.tsv )
    reference_fasta_ch  // value channel: reference multi-FASTA whose headers match metagraph's alignment accessions

    main:
    // Samples with nothing above metagraph_align_min_hits emit nothing here (optional
    // outputs) and are naturally dropped by every join() below — no branching needed.
    EXTRACT_METAGRAPH_REFERENCE_SUBSET(record_ids_ch, reference_fasta_ch)

    BOWTIE_INDEX(EXTRACT_METAGRAPH_REFERENCE_SUBSET.out.subset_fasta)

    // BOWTIE2SAMTOOLS wants (meta, r1, r2, bt2_files, index_prefix); BOWTIE_INDEX names its
    // index files "${meta.id}_index*", so the prefix string is derived rather than tracked
    // through the channel.
    mapping_input_ch = reads_ch
        .join(BOWTIE_INDEX.out.bowtie_index)
        .map { meta, r1, r2, bt2_files -> tuple(meta, r1, r2, bt2_files, "${meta.id}_index") }

    BOWTIE2SAMTOOLS(mapping_input_ch, params.metagraph_map_bowtie_threads)

    SAMTOOLS_COVERAGE(BOWTIE2SAMTOOLS.out.bam_file)

    qc_input_ch = SAMTOOLS_COVERAGE.out.coverage
        .join(SAMTOOLS_COVERAGE.out.query_lengths)
        .join(index_label_map_ch)
        .join(species_hits_ch)

    AGGREGATE_METAGRAPH_COVERAGE(qc_input_ch)

    GENERATE_METAGRAPH_MAP_SUMMARY(
        AGGREGATE_METAGRAPH_COVERAGE.out.qc_table.map { _meta, qc_table -> qc_table }.collect()
    )

    emit:
    qc_table = AGGREGATE_METAGRAPH_COVERAGE.out.qc_table
}

/*
========================================================================================
    THE END
========================================================================================
*/
