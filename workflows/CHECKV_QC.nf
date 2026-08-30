include { LINK_BIN_SCAFFOLDS } from '../modules/vrhyme.nf'
include { CHECKV             } from '../modules/checkv.nf'

workflow CHECKV_QC {

    take:
    virus_fna_ch  // tuple val(meta), path(virus_fna)   — from GENOMAD_CLASSIFY.out.virus_fna
    bins_fasta_ch // tuple val(meta), path(bins_dir)    — from VRHYME_BIN.out.bins_fasta (always emitted; empty dir when vRhyme produced no bins)

    main:

    // 1) CheckV on geNomad's virus.fna — runs for every sample where geNomad produced
    //    an output, irrespective of whether vRhyme later ran or produced a bin.
    virus_fna_ch
        .map { meta, fna -> tuple(meta, 'virus_scaffolds', fna) }
        .dump(tag: 'checkv_virus_input')
        .set { ch_checkv_virus_input }

    // 2) CheckV on linked vRhyme bins — only for samples where vRhyme produced at
    //    least one bin. bins_fasta_ch now always carries a real directory (vRhyme
    //    guarantees it, empty or not, when it finds no bins), so filter out the
    //    empty ones here rather than handing link_bin_sequences.py nothing to link.
    //    LINK_BIN_SCAFFOLDS joins each bin's scaffolds with N-spacers into a single
    //    sequence per bin (per the vRhyme-recommended approach for feeding
    //    multi-scaffold bins to scaffold-level QC tools).
    bins_fasta_ch
        .filter { meta, bins_dir -> (bins_dir.toFile().list() ?: []).size() > 0 }
        .set { ch_bins_fasta_present }

    LINK_BIN_SCAFFOLDS(ch_bins_fasta_present)

    LINK_BIN_SCAFFOLDS.out.linked
        .filter { meta, fna -> fna.size() > 0 }
        .map { meta, fna -> tuple(meta, 'linked_bins', fna) }
        .dump(tag: 'checkv_linked_bins_input')
        .set { ch_checkv_linked_input }

    CHECKV(ch_checkv_virus_input.mix(ch_checkv_linked_input))

    // Per-scaffold quality calls (source == 'virus_scaffolds' only, dropping
    // the 'linked_bins' rows) — one guaranteed-real row set per sample, since
    // quality_summary is no longer an optional CHECKV output. Consumed by
    // VCONTACT3_RUN to promote high/medium-quality binned scaffolds.
    CHECKV.out.quality_summary
        .filter { meta, source, tsv -> source == 'virus_scaffolds' }
        .map { meta, source, tsv -> tuple(meta, tsv) }
        .set { ch_virus_scaffolds_quality }

    emit:
    quality_summary                 = CHECKV.out.quality_summary
    completeness                    = CHECKV.out.completeness
    contamination                   = CHECKV.out.contamination
    complete_genomes                = CHECKV.out.complete_genomes
    virus_scaffolds_quality_summary = ch_virus_scaffolds_quality
}
