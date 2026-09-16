// --- De novo assembly + viral binning (rvi_integration_1) --------------------
// Extracted unchanged from main.nf's inline body. Ported from
// rvi-viral-metagenomics-pipeline/main.nf, wiring unchanged (see
// docs/nf-metro/route_map.mmd "De novo assembly" + "Viral binning" sections).
include {ASSEMBLE_META} from '../rvi_toolbox/subworkflows/assemble.nf'
include {GENOMAD_CLASSIFY} from '../rvi_toolbox/subworkflows/genomad.nf'
include {VRHYME_BIN} from '../rvi_toolbox/subworkflows/vrhyme.nf'
include {CHECKV_QC} from '../rvi_toolbox/subworkflows/checkv.nf'
include {VCONTACT3_RUN} from '../workflows/VCONTACT3_RUN.nf'
include {ASSEMBLY_REPORTS} from '../rvi_toolbox/subworkflows/assembly_reports.nf'

workflow ASSEMBLY {
    take:
        preprocessed_3tuple_ch // tuple (meta, read1, read2)

    main:
        ASSEMBLE_META(preprocessed_3tuple_ch)
        GENOMAD_CLASSIFY(ASSEMBLE_META.out.contigs_channel)
        VRHYME_BIN(
            GENOMAD_CLASSIFY.out.virus_fna,
            GENOMAD_CLASSIFY.out.virus_summary,
            preprocessed_3tuple_ch
        )
        CHECKV_QC(
            GENOMAD_CLASSIFY.out.virus_fna,
            VRHYME_BIN.out.bins_fasta
        )
        VCONTACT3_RUN(
            GENOMAD_CLASSIFY.out.virus_proteins,
            GENOMAD_CLASSIFY.out.virus_summary,
            VRHYME_BIN.out.membership,
            VRHYME_BIN.out.bins_fasta,
            CHECKV_QC.out.virus_scaffolds_quality_summary
        )

        // The lane's three run-level CSVs (sample / scaffold / vMAG level), built by
        // ASSEMBLY_REPORTS from the modules' own output files. Downstream of vContact3
        // for its taxonomy.
        //
        // This lane deliberately writes no per-sample properties.json. It used to, via
        // GENERATE_ASSEMBLY_REPORT, off a `meta` enriched with counts derived here --
        // but assembly_reports.py never read that file, so it was a second, parallel
        // derivation of the same numbers feeding nothing. Worse, it came off a chain of
        // inner joins that dropped any sample vRhyme never ran for: on a 95-sample run
        // it covered 35. The CSVs are the lane's report.
        ASSEMBLY_REPORTS(
            GENOMAD_CLASSIFY.out.virus_summary,
            VRHYME_BIN.out.membership,
            CHECKV_QC.out.virus_scaffolds_quality_summary,
            CHECKV_QC.out.linked_bins_quality_summary,
            VCONTACT3_RUN.out.postprocessed_assignments
        )
}
