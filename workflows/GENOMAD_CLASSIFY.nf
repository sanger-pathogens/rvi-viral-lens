#!/usr/bin/env nextflow

/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULES
//
include { GENOMAD } from '../modules/genomad.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow GENOMAD_CLASSIFY {

    take:
    contigs_channel // tuple val(meta), path(contigs), path(scaffolds)

    main:


    contigs_channel
        .map { meta, contigs, scaffolds -> [meta, scaffolds] }
        .dump(tag: 'scaffolds')
        .set { ch_scaffolds }

    GENOMAD(ch_scaffolds)

    emit:
    virus_summary  = GENOMAD.out.virus_summary
    virus_fna      = GENOMAD.out.virus_fna
    summary_json   = GENOMAD.out.summary_json
    virus_genes    = GENOMAD.out.virus_genes
    virus_proteins = GENOMAD.out.virus_proteins
}

/*
========================================================================================
    THE END
========================================================================================
*/
