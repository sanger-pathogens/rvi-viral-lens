#!/usr/bin/env nextflow

/*
========================================================================================
    SCRuB cross-contamination decontamination sub-workflow (optional, off by default)
    ------------------------
    Reformats the pipeline's whole-run Bracken species-abundance summary into SCRuB's
    expected samples x species orientation, then runs SCRuB against it and a
    user-supplied plate map to detect/correct cross-sample contamination (well-to-well
    leakage, shared control-sample contamination). Runs once per pipeline run, not
    per-sample, since decontamination inherently needs the whole batch (samples +
    controls) together. Also renders a before/after/change relative-abundance heatmap
    for visual QC. Results are published under <outdir>/abundance/scrub.
========================================================================================
*/

/*
========================================================================================
    IMPORT MODULES
========================================================================================
*/
include { REFORMAT_BRACKEN_FOR_SCRUB } from '../modules/reformat_bracken.nf'
include { SCRUB                      } from '../modules/scrub.nf'
include { SCRUB_HEATMAP              } from '../modules/scrub_heatmap.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/
workflow SCRUB_DECONTAM {

    take:
    bracken_summary_ch   // path: bracken_summary_report.tsv (single whole-run file)

    main:
    plate_map_ch = Channel.fromPath(params.scrub_plate_map)

    REFORMAT_BRACKEN_FOR_SCRUB(bracken_summary_ch, plate_map_ch)
    SCRUB(REFORMAT_BRACKEN_FOR_SCRUB.out.formatted, plate_map_ch)
    SCRUB_HEATMAP(REFORMAT_BRACKEN_FOR_SCRUB.out.formatted, SCRUB.out.scrub_output)

    emit:
    scrub_output = SCRUB.out.scrub_output
    heatmap      = SCRUB_HEATMAP.out.heatmap
}

/*
========================================================================================
    THE END
========================================================================================
*/
