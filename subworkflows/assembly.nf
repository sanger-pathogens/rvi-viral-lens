// --- De novo assembly + viral binning (rvi_integration_1) --------------------
// Extracted unchanged from main.nf's inline body. Ported from
// rvi-viral-metagenomics-pipeline/main.nf, wiring unchanged (see
// docs/nf-metro/route_map.mmd "De novo assembly" + "Viral binning" sections).
include {ASSEMBLE_META} from '../workflows/ASSEMBLE_META.nf'
include {GENOMAD_CLASSIFY} from '../workflows/GENOMAD_CLASSIFY.nf'
include {VRHYME_BIN} from '../workflows/VRHYME_BIN.nf'
include {CHECKV_QC} from '../workflows/CHECKV_QC.nf'
include {VCONTACT3_RUN} from '../workflows/VCONTACT3_RUN.nf'
include {GENERATE_ASSEMBLY_REPORT} from '../workflows/GENERATE_ASSEMBLY_REPORT.nf'
include {ASSEMBLY_REPORTS} from '../workflows/ASSEMBLY_REPORTS.nf'
include {publish_lane_json} from '../modules/publish_lane_report.nf'
include {publish_run_files as publish_assembly_run_files} from '../modules/publish_lite.nf'

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

        // --- sample-level meta enrichment for the assembly report ---
        // Each module's own file outputs still publish unchanged below; this only
        // adds small counts (see docs/nf-metro/README.md's linked plan) so
        // GENERATE_ASSEMBLY_REPORT can populate its report straight from meta.
        sample_ids_ch = ASSEMBLE_META.out.contigs_channel
            .map { meta, _contigs, _scaffolds -> meta.id }

        genomad_counts_ch = GENOMAD_CLASSIFY.out.virus_summary
            .map { meta, tsv -> [meta.id, meta, count_genomad_summary(tsv)] }

        vrhyme_counts_ch = VRHYME_BIN.out.membership
            .map { meta, tsv -> [meta.id, count_vrhyme_membership(tsv)] }

        checkv_counts_ch = CHECKV_QC.out.virus_scaffolds_quality_summary
            .map { meta, tsv -> [meta.id, count_checkv_quality(tsv)] }

        // vContact3 runs once for the whole batch (final_assignments.csv is not
        // per-sample), so fan its single output out to every sample and let each
        // sample pick its own rows back out by the `<sample_id>||` genome_id prefix
        // vcontact3_prep.py mints (see bin/vcontact3_prep.py).
        vcontact3_counts_ch = sample_ids_ch
            .combine(VCONTACT3_RUN.out.final_assignments.ifEmpty(null))
            .map { sample_id, csv -> [sample_id, count_vcontact3_for_sample(csv, sample_id)] }

        genomad_counts_ch
            .join(vrhyme_counts_ch)
            .join(checkv_counts_ch)
            .join(vcontact3_counts_ch)
            .map { id, meta, g_counts, v_counts, c_counts, vc3_counts ->
                def new_meta = meta + g_counts + v_counts + c_counts + vc3_counts
                [id, new_meta]
            }
            .set { assembly_report_prep_ch }

        GENERATE_ASSEMBLY_REPORT(assembly_report_prep_ch)

        // The lane's three CSV reports (sample / scaffold / vMAG level). Separate
        // from GENERATE_ASSEMBLY_REPORT above, which still writes the per-sample
        // properties.json + assembly_run_summary.json from `meta`: these are
        // multi-row-per-sample tables built from the modules' own output files,
        // which `meta` cannot carry. Downstream of vContact3 for its taxonomy.
        ASSEMBLY_REPORTS(
            GENOMAD_CLASSIFY.out.virus_summary,
            VRHYME_BIN.out.membership,
            CHECKV_QC.out.virus_scaffolds_quality_summary,
            CHECKV_QC.out.linked_bins_quality_summary,
            VCONTACT3_RUN.out.postprocessed_assignments
        )

        // PUBLISH (assembly lane)
        publish_lane_json(GENERATE_ASSEMBLY_REPORT.out.publish_seq_level_ch)
        publish_assembly_run_files(GENERATE_ASSEMBLY_REPORT.out.publish_run_level_summaries_ch)
}

// --- rvi_integration_1: sample-level count helpers for the assembly report ---
// Same pattern as rvi_toolbox's own vrhyme.nf count_long_scaffolds() helper: read
// a small per-sample TSV/CSV straight from the resolved output path in a channel
// .map{} closure. Naive split() parsing (not a real CSV/TSV reader) -- fine for
// these known, script-generated files, but would need replacing if any of these
// columns can carry embedded delimiters.

def count_genomad_summary(tsv) {
    // virus_summary.tsv: seq_name, n_genes, length, ... (bin/vcontact3_prep.py's load_summary)
    def lines = tsv.readLines()
    if (lines.size() < 2) return [genomad_n_scaffolds: 0, genomad_n_eligible: 0]
    def header = lines[0].split('\t')
    def n_genes_idx = header.findIndexOf { String col -> col == 'n_genes' }
    def n_total = lines.size() - 1
    def n_eligible = 0
    if (n_genes_idx >= 0) {
        n_eligible = lines[1..-1].count { String line ->
            def cols = line.split('\t')
            n_genes_idx < cols.size() && (cols[n_genes_idx] as Integer) > 0
        }
    }
    return [genomad_n_scaffolds: n_total, genomad_n_eligible: n_eligible]
}

def count_vrhyme_membership(tsv) {
    // vRhyme_best_bins.*.membership.tsv: header "scaffold\tbin" (modules/vrhyme.nf)
    def lines = tsv.readLines()
    if (lines.size() < 2) return [vrhyme_n_bins: 0, vrhyme_n_binned_scaffolds: 0]
    def bins = lines[1..-1].collect { line -> line.split('\t')[1] }
    // Groovy's List.unique() mutates in place and returns the same list, so
    // reading bins.size() after bins.unique().size() yields the DISTINCT count,
    // not the row count. Take the row count first.
    def n_binned_scaffolds = bins.size()
    return [vrhyme_n_bins: bins.unique().size(), vrhyme_n_binned_scaffolds: n_binned_scaffolds]
}

def count_checkv_quality(tsv) {
    // quality_summary.tsv: has contig_id, checkv_quality (bin/vcontact3_prep.py's load_scaffold_quality)
    def lines = tsv.readLines()
    if (lines.size() < 2) return [checkv_n_high_quality: 0, checkv_n_medium_quality: 0]
    def header = lines[0].split('\t')
    def q_idx = header.findIndexOf { String col -> col == 'checkv_quality' }
    if (q_idx < 0) return [checkv_n_high_quality: 0, checkv_n_medium_quality: 0]
    def qualities = lines[1..-1].collect { line -> line.split('\t')[q_idx] }
    return [
        checkv_n_high_quality:   qualities.count { String q -> q == 'High-quality' },
        checkv_n_medium_quality: qualities.count { String q -> q == 'Medium-quality' }
    ]
}

def count_vcontact3_for_sample(csv, sample_id) {
    // final_assignments.csv: has Genome, genus_prediction (bin/vcontact3_postprocess.py).
    // Genome IDs are "<sample_id>||bin_N" / "<sample_id>||scaffold[||standalone]"
    // (bin/vcontact3_prep.py), so this sample's rows are picked out by that prefix.
    if (csv == null || !csv.exists()) {
        return [vcontact3_n_genomes: 0, vcontact3_n_clusters: 0]
    }
    def lines = csv.readLines()
    if (lines.size() < 2) return [vcontact3_n_genomes: 0, vcontact3_n_clusters: 0]
    def header = lines[0].split(',')
    def genome_idx = header.findIndexOf { String col -> col == 'Genome' }
    def genus_idx = header.findIndexOf { String col -> col == 'genus_prediction' }
    if (genome_idx < 0) return [vcontact3_n_genomes: 0, vcontact3_n_clusters: 0]
    def sample_rows = lines[1..-1].findAll { line ->
        def cols = line.split(',')
        genome_idx < cols.size() && cols[genome_idx].contains("${sample_id}||")
    }
    def clusters = (genus_idx >= 0)
        ? sample_rows.collect { line -> line.split(',')[genus_idx] }.unique()
        : []
    return [vcontact3_n_genomes: sample_rows.size(), vcontact3_n_clusters: clusters.size()]
}
