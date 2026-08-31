#!/usr/bin/env nextflow

/*
========================================================================================
    KRAKEN2BRACKEN (viral-lens fork, with emits)
    ------------------------
    Identical orchestration to rvi_toolbox/subworkflows/kraken2bracken.nf (same shared
    modules, same control flow) -- that subworkflow has no `emit:` block at all, so
    nothing it produces is reachable from outside it. Reporting (see
    GENERATE_ABUNDANCE_REPORT.nf) and SCRuB decontamination (see SCRUB_DECONTAM.nf,
    which needs the whole-run abundance_summary) both need something out of this, so
    this is a viral-lens-owned fork adding an emit: block rather than a modification to
    the shared rvi_toolbox submodule (see INSTRUCT.md's "rvi_toolbox fork problem").
    Keep this in sync with rvi_toolbox's version if that one changes.
========================================================================================
*/

/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/
include { KRAKEN2; KRAKEN2_GET_CLASSIFIED; COMPRESS_READS } from '../rvi_toolbox/modules/kraken2.nf'
include { BRACKEN } from '../rvi_toolbox/modules/bracken.nf'
include { KREPORT2MPA; GENERATE_ABUNDANCE_SUMMARY } from '../rvi_toolbox/modules/krakentools.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow KRAKEN2BRACKEN {

    take:
    ch_reads // meta, paired_reads

    main:

    Channel.fromPath(params.kraken2bracken_kraken2_db)
        .set { ch_kraken2_db }

    //
    // CLASSIFICATION
    //
    ch_reads
        .combine(ch_kraken2_db)
        .set { ch_reads_and_kraken2_db }
    if (params.kraken2bracken_get_classified_reads) {
        KRAKEN2_GET_CLASSIFIED(ch_reads_and_kraken2_db)

        KRAKEN2_GET_CLASSIFIED.out.kraken2_sample_report.set { ch_kraken2_sample_report }
        KRAKEN2_GET_CLASSIFIED.out.classified_reads.set { ch_kraken2_classified_reads }
        KRAKEN2_GET_CLASSIFIED.out.unclassified_reads.set { ch_kraken2_unclassified_reads }
        ch_kraken2_classified_reads.join(ch_kraken2_unclassified_reads)
            .map { meta, classified, unclassified -> [meta, classified + unclassified] }
            .set { ch_reads_to_compress }

        COMPRESS_READS(
            ch_reads_to_compress
        )
    } else {
        KRAKEN2(ch_reads_and_kraken2_db)
        KRAKEN2.out.kraken2_sample_report.set { ch_kraken2_sample_report }
    }

    //
    // ABUNDANCE ESTIMATION
    //
    kraken2_db_dir = file(params.kraken2bracken_kraken2_db, checkIfExists: true)
    // Assume a pre-built bracken database file has been generated from the given kraken2 database and moved into this database directory
    required_kmer_distrib = file("${kraken2_db_dir}/database${params.kraken2bracken_read_len}mers.kmer_distrib")
    if (!required_kmer_distrib.exists()) {
        log.error("The required bracken kmer distribution database file cannot be found in the kraken database directory ${kraken2_db_dir}")
        exit 1
    } else {
        Channel.fromPath(required_kmer_distrib)
            .set { ch_kmer_distrib }
    }

    ch_kraken2_sample_report
        .combine(ch_kmer_distrib)
        .set { ch_kraken2_report_and_kmer_distrib }
    BRACKEN(
        ch_kraken2_report_and_kmer_distrib
    )

    KREPORT2MPA(BRACKEN.out.kraken_style_bracken_report)

    //
    // SUMMARISE ABUNDANCE
    //
    KREPORT2MPA.out.mpa_abundance_report
        .map { meta, report -> report }
        .collect()
        .set { ch_mpa_abundance_reports }
    GENERATE_ABUNDANCE_SUMMARY(
        ch_mpa_abundance_reports
    )

    emit:
    bracken_report       = BRACKEN.out.bracken_report             // tuple(meta, *.bracken) -- per sample
    mpa_abundance_report = KREPORT2MPA.out.mpa_abundance_report    // tuple(meta, *.mpa.txt) -- per sample
    abundance_summary    = GENERATE_ABUNDANCE_SUMMARY.out.abundance_summary  // path: bracken_summary_report.tsv -- whole run
}

/*
========================================================================================
    THE END
========================================================================================
*/
