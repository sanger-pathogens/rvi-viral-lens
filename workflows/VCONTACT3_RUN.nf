include { VCONTACT3_PER_SAMPLE_PREP; VCONTACT3; VCONTACT3_POSTPROCESS } from '../modules/vcontact3.nf'

workflow VCONTACT3_RUN {

    take:
    virus_proteins_ch    // tuple val(meta), path(virus_proteins.faa) — from GENOMAD_CLASSIFY.out.virus_proteins
    virus_summary_ch     // tuple val(meta), path(virus_summary.tsv)  — from GENOMAD_CLASSIFY.out.virus_summary
    bin_membership_ch    // tuple val(meta), path(membership.tsv)     — from VRHYME_BIN.out.membership (always emitted, real even with 0 bins)
    bins_fasta_ch        // tuple val(meta), path(bins_fasta_dir)     — from VRHYME_BIN.out.bins_fasta (always emitted, real even with 0 bins)
    checkv_quality_ch    // tuple val(meta), path(quality_summary.tsv)— from CHECKV_QC.out.virus_scaffolds_quality_summary (always emitted, real even if empty)

    main:

    // Per-sample prep input: [meta, virus_proteins, virus_summary, membership, bins_fasta, checkv_quality]
    virus_proteins_ch
        .join(virus_summary_ch)
        .join(bin_membership_ch)
        .join(bins_fasta_ch)
        .join(checkv_quality_ch)
        .dump(tag: 'vcontact3_per_sample_inputs')
        .set { ch_per_sample_inputs }

    // Stage the prep script through the channel instead of referencing
    // ${projectDir} inside the process — the latter doesn't get staged into
    // the task workdir and is flagged as a Nextflow lint warning.
    // Use Channel.value(...) so the single script broadcasts to every task
    // rather than being consumed by the first one only.
    prep_script_ch = Channel.value(file("${projectDir}/bin/vcontact3_prep.py"))

    VCONTACT3_PER_SAMPLE_PREP(ch_per_sample_inputs, prep_script_ch)

    // Pipeline-level fan-in: every sample's prep output must be in before
    // VCONTACT3 starts. `.collect()` is the barrier — it waits for the
    // upstream channel to be exhausted, which only happens once all upstream
    // vRhyme processes (and therefore all samples) have finished.
    VCONTACT3_PER_SAMPLE_PREP.out.proteins
        .map { meta, faa -> faa }
        .collect()
        .dump(tag: 'vcontact3_collected_proteins')
        .set { all_proteins }

    VCONTACT3_PER_SAMPLE_PREP.out.gene2genome
        .map { meta, tsv -> tsv }
        .collect()
        .dump(tag: 'vcontact3_collected_gene2genome')
        .set { all_gene2genome }

    VCONTACT3_PER_SAMPLE_PREP.out.lengths
        .map { meta, tsv -> tsv }
        .collect()
        .dump(tag: 'vcontact3_collected_lengths')
        .set { all_lengths }

    VCONTACT3(all_proteins, all_gene2genome, all_lengths)

    // Stage the postprocess script the same way as the prep script above.
    postprocess_script_ch = Channel.value(file("${projectDir}/bin/vcontact3_postprocess.py"))

    VCONTACT3_POSTPROCESS(
        VCONTACT3.out.final_assignments,
        VCONTACT3.out.vog_support,
        postprocess_script_ch
    )

    emit:
    final_assignments         = VCONTACT3.out.final_assignments
    performance_metrics       = VCONTACT3.out.performance_metrics
    postprocessed_assignments = VCONTACT3_POSTPROCESS.out.postprocessed
    novel_taxa_assignments    = VCONTACT3_POSTPROCESS.out.novel_taxa
    postprocess_report        = VCONTACT3_POSTPROCESS.out.report
}
