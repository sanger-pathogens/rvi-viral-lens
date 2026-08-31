#!/usr/bin/env nextflow
// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

// enable dsl2
nextflow.enable.dsl = 2

// --- import modules ---------------------------------------------------------
include {check_sort_reads_params} from './workflows/SORT_READS_BY_REF.nf'
include {validateParameters; paramsSummaryLog} from 'plugin/nf-schema'

include {PREPROCESSING} from "./rvi_toolbox/subworkflows/preprocessing.nf"
include {SORT_READS_BY_REF} from './workflows/SORT_READS_BY_REF.nf'
include {GENERATE_CONSENSUS} from './workflows/GENERATE_CONSENSUS.nf'
include {SCOV2_SUBTYPING} from './workflows/SCOV2_SUBTYPING.nf'
include {GENERATE_CLASSIFICATION_REPORT} from './workflows/GENERATE_CLASSIFICATION_REPORT.nf'
include {RUN_NEXTCLADE} from './workflows/RUN_NEXTCLADE.nf'
include {publish_consensus_files as publish_aln_files; publish_consensus_files as publish_nc_files; publish_consensus_files as publish_per_sample_json} from './modules/publish_lite.nf'
include {publish_run_files; publish_run_files as publish_assembly_run_files; publish_run_files as publish_mapping_run_files} from './modules/publish_lite.nf'

// --- De novo assembly + viral binning (rvi_integration_1) --------------------
// Ported from rvi-viral-metagenomics-pipeline/main.nf, wiring unchanged (see
// docs/nf-metro/route_map.mmd "De novo assembly" + "Viral binning" sections).
// Off by default (params.do_assembly) until decisions B/C/D in that route map are
// settled and the mapping/abundance lanes exist to feed the same Reporting block.
include {ASSEMBLE_META} from './workflows/ASSEMBLE_META.nf'
include {GENOMAD_CLASSIFY} from './workflows/GENOMAD_CLASSIFY.nf'
include {VRHYME_BIN} from './workflows/VRHYME_BIN.nf'
include {CHECKV_QC} from './workflows/CHECKV_QC.nf'
include {VCONTACT3_RUN} from './workflows/VCONTACT3_RUN.nf'
include {GENERATE_ASSEMBLY_REPORT} from './workflows/GENERATE_ASSEMBLY_REPORT.nf'
include {publish_lane_json; publish_lane_json as publish_mapping_lane_json} from './modules/publish_lane_report.nf'

// -- map reads to sequence indexes (rvi_integration_1)
include {VIRAL_MSWEEP} from './workflows/VIRAL_MSWEEP.nf'
include {VIRAL_METAGRAPH_ALIGN} from './workflows/VIRAL_METAGRAPH_ALIGN.nf'
include {VIRAL_METAGRAPH_QUERY} from './workflows/VIRAL_METAGRAPH_QUERY.nf'
include {GENERATE_MAPPING_REPORT} from './workflows/GENERATE_MAPPING_REPORT.nf'
// GENERATE_ABUNDANCE_REPORT already exists (workflows/GENERATE_ABUNDANCE_REPORT.nf) and is
// covered by its own nf-test workflow test, but isn't wired in below yet: its upstream
// lane (abundance estimation + SCRuB) isn't ported yet -- decisions A/D on the route map
// are still open.

// Main entry-point workflow
workflow {
  /*
  * ANSI escape codes to color output messages
  */
  ANSI_GREEN = "\033[1;32m"
  ANSI_RED = "\033[1;31m"
  ANSI_RESET = "\033[0m"
  _ANSI_BOLD = "\033[1m"


  log.info """${ANSI_RESET}
  ===========================================
  Viral Lens [v1.5.1]
  Used parameters:
  -------------------------------------------
  --> general pipeline parameters:
    --outdir                   : ${params.outdir}

  --> SORT_READS_BY_REF workflow parameters:
    --manifest                   : ${params.manifest}
    --db_path                    : ${params.db_path}
    --db_library_fa_path         : ${params.db_library_fa_path}
    --min_reads_for_taxid        : ${params.min_reads_for_taxid}
    --k2r_max_total_reads_per_fq : ${params.k2r_max_total_reads_per_fq}

  --> GENERATE_CONSENSUS workflow parameters:
    --do_consensus_polishing      : ${params.do_consensus_polishing }
    --read_aligner                : ${params.read_aligner}
    --read_aligner_params         : ${params.read_aligner_params}
    --mpileup_max_depth           : ${params.mpileup_max_depth }
    --ivar_initial_min_depth      : ${params.ivar_initial_min_depth}
    --ivar_initial_freq_threshold : ${params.ivar_initial_freq_threshold}
    --ivar_polish_min_depth       : ${params.ivar_polish_min_depth}
    --ivar_polish_freq_threshold  : ${params.ivar_polish_freq_threshold}

  --> viral subtyping branching parameters:
    --scv2_keyword             : ${params.scv2_keyword}

  --> Nextclade parameters:
    --nextclade_index_json      : ${params.nextclade_index_json}

  --> De novo assembly + viral binning parameters (rvi_integration_1):
    --do_assembly               : ${params.do_assembly}

  --> resource management:
    --default_error_strategy   : ${params.default_error_strategy}
    --mem_k2r_b0_offset        : ${params.mem_k2r_b0_offset}
    --mem_k2r_b0               : ${params.mem_k2r_b0}
    --mem_k2r_b0_final         : ${params.mem_k2r_b0_final}
    --mem_k2r_b1               : ${params.mem_k2r_b1}
    --mem_k2r_f1               : ${params.mem_k2r_f1}
    --mem_k2r_a2               : ${params.mem_k2r_a2}
    --max_attempts             : ${params.max_attempts}
  ------------------------------------------
  Runtime data:
  -------------------------------------------
  Running with profile:   ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
  Running as user:        ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
  Launch dir:             ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
  Base dir:               ${ANSI_GREEN}${baseDir}${ANSI_RESET}
  ------------------------------------------
""".stripIndent()

    // Validate input parameters
    validateParameters()
    // Print summary of supplied parameters
    log.info paramsSummaryLog(workflow)

    // === 1 - Process input ===
    check_main_params()
    // ==========================
    reads_ch = parse_mnf(params.manifest) // tuple(meta, [fastq_1, fastq_2])

    // === Preprocessing ===
    // preprocessed_3tuple_ch (meta, read1, read2) mirrors PREPROCESSING's own
    // output shape either way, and feeds both the existing taxid lane (reshaped
    // below into sort_reads_in_ch) and the new de novo assembly lane.
    if (params.do_preprocessing) {
        reads_ch.map{ meta, fastqs ->
            return [meta, fastqs[0], fastqs[1]]
        }.set{preproc_in_ch}

        PREPROCESSING(preproc_in_ch)
        PREPROCESSING.out.out_ch.set{ preprocessed_3tuple_ch }

        preprocessed_3tuple_ch.map{meta, read1, read2 ->
            return [meta, [read1, read2]]
        }.set{sort_reads_in_ch}

    } else {
        sort_reads_in_ch = reads_ch
        // file() matters: parse_mnf yields the manifest's raw strings, while
        // PREPROCESSING emits real paths. Consumers that only declare `path`
        // inputs coerce either, but ASSEMBLE_META calls R1.countFastq() in a
        // map closure, which needs a Path. Coerce here so both branches really
        // do emit the same shape, as the comment above claims.
        reads_ch.map{ meta, fastqs ->
            return [meta, file(fastqs[0]), file(fastqs[1])]
        }.set{ preprocessed_3tuple_ch }
    }

    // ==========================
    // === 2 - Map reads to taxid
    SORT_READS_BY_REF(sort_reads_in_ch)
    // === 3 - Generate consensus ==
    GENERATE_CONSENSUS( SORT_READS_BY_REF.out.sample_taxid_ch )

    // === Run Nextclade
    GENERATE_CONSENSUS.out.filtered_consensus_ch
    .map{meta, _bam, _bam_idx, consensus, _qc_json -> [meta.id, meta, consensus]}
    .set{consensus_fa_ch}

    SORT_READS_BY_REF.out.sample_pre_report_ch
    .map{meta ->
        def new_meta = meta + [id:"${meta.sample_id}.${meta.selected_taxid}"]
        return [new_meta.id, new_meta ]
    }
    .combine(consensus_fa_ch, by:0)
    .map{_id, pre_report_meta, fa_meta, fa ->
        def final_meta = pre_report_meta + [reference_header:"${fa_meta.reference_header}", taxid:"${fa_meta.taxid}"]
        [final_meta, fa]
    }
    .set {nextclade_In_ch}

    // TODO add check parameters
    if (params.nextclade_index_json == null){
        log.warn("No nextclade_index_json provided, skipping nextclade analysis step")
        publish_nextclade_outputs_ch = channel.empty()
        per_consensus_nextclade_json_ch = channel.empty()
    } else {
        RUN_NEXTCLADE(nextclade_In_ch)
        RUN_NEXTCLADE.out
            .map{meta, _agg_json, tar_gz -> [meta, tar_gz] }
            .set { publish_nextclade_outputs_ch }

        RUN_NEXTCLADE.out
            .map{meta, json, _tarball ->
                    def join_key = meta.id
                    [join_key, json]
                }
            .set{ per_consensus_nextclade_json_ch }
    }

    // === 5 - branching output from generate_consensus for viral specific subtyping

    SORT_READS_BY_REF.out.sample_pre_report_ch
        .map{ it ->
            def id="${it.sample_id}.${it.selected_taxid}"
            [id, it]
        }
        .set{sample_report_with_join_key_ch}

    // 5.1 - add report info to out qc metric chanel and branch for SCOV2 subtyping
    GENERATE_CONSENSUS.out.filtered_consensus_ch
        .map { meta, _bam, _bam_idx, consensus, _qc ->
            [meta.id, meta, consensus]
        }
        .join(sample_report_with_join_key_ch)
        .map {_id, meta, fasta, report ->
            def new_meta = meta.plus(report)
            [new_meta, fasta]
        }
        .branch{ it ->
            scv2_subtyping_workflow_in_ch: it[0].ref_selected.contains("${params.scv2_keyword}")
            no_subtyping_ch: true
        }
        .set { filtered_consensus_by_type_ch }

    // 5.2 - do SCOV2 subtyping
    if (params.do_scov2_subtyping == true){
        SCOV2_SUBTYPING(filtered_consensus_by_type_ch.scv2_subtyping_workflow_in_ch)
        SCOV2_SUBTYPING.out.set{scov2_subtyped_ch}
    }

    // === 6 - write final classification reports

    if (!params.do_scov2_subtyping == true){
        scov2_subtyped_ch = channel.empty()
    }
    filtered_consensus_by_type_ch.no_subtyping_ch.concat(scov2_subtyped_ch)
        .map{ meta, _fasta ->  [meta.id, meta] }
        .set{ report_in_ch }

    GENERATE_CONSENSUS.out.filtered_consensus_ch.map { meta, _bam, _bam_idx, _consensus, qc_json ->
                    [meta.id, qc_json]
                    }
                .set{ qc_json_simplified_ch }

    GENERATE_CONSENSUS.out.filtered_consensus_ch.map { meta, bam, bam_idx, consensus, _qc_json ->
                    [meta, [bam, bam_idx, consensus]]
                    }
                .set{ aln_publish_ch }

    report_in_ch // [meta.id, meta]
            .join(qc_json_simplified_ch, remainder: true) // [meta.id, meta, qc_json]
            .join(per_consensus_nextclade_json_ch, remainder: true) //  [meta.id, meta, qc_json, nc_json]
            .set{ report_prep_ch }

    GENERATE_CLASSIFICATION_REPORT(report_prep_ch)

    // === 7 - De novo assembly + viral binning (rvi_integration_1, opt-in) ===
    if (params.do_assembly) {
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

        // PUBLISH (assembly lane)
        publish_lane_json(GENERATE_ASSEMBLY_REPORT.out.publish_seq_level_ch)
        publish_assembly_run_files(GENERATE_ASSEMBLY_REPORT.out.publish_run_level_summaries_ch)
    }

    // === 8 - Map reads to sequence indexes (rvi_integration_1, opt-in) ===
    // Up to three methods run in parallel off the same preprocessed reads the assembly
    // lane uses (not downstream of it), each independently gated, all feeding ONE
    // GENERATE_MAPPING_REPORT call -- see INSTRUCT.md item 3: "feed its per-sample
    // counts into the same mapping_report_prep_ch rather than building a second report
    // path." sequence_index_sample_ch is the join backbone (every sample that reaches
    // this lane) so a sample report row exists even if only one method ran for it, or
    // (with more than one flag on) one row carries every enabled method's counts.
    if (params.do_sequence_index) {
        sequence_index_sample_ch = preprocessed_3tuple_ch
            .map { meta, _r1, _r2 -> [meta.id, meta] }

        // -- Themisto2 pseudoalignment + mSWEEP abundance, then breadth-of-coverage
        // validation of the low-abundance calls.
        if (params.run_msweep) {
            VIRAL_MSWEEP(preprocessed_3tuple_ch)

            msweep_counts_ch = VIRAL_MSWEEP.out.abundances
                .map { meta, abundances, _probs -> [meta.id, count_msweep_abundances(abundances)] }

            // map_qc drops samples with nothing above msweep_map_min_abundance (the
            // MSWEEP_MAP_QC processes emit optional outputs); joined with remainder
            // below rather than losing those samples from the report entirely.
            mapqc_counts_ch = VIRAL_MSWEEP.out.map_qc
                .map { meta, qc_tsv -> [meta.id, count_msweep_map_qc(qc_tsv)] }
        } else {
            msweep_counts_ch = Channel.empty()
            mapqc_counts_ch = Channel.empty()
        }

        // -- Sequence-to-graph alignment via Metagraph (metagraph align), then its own
        // map-QC validation.
        if (params.run_metagraph_align) {
            VIRAL_METAGRAPH_ALIGN(preprocessed_3tuple_ch)

            metagraph_align_counts_ch = VIRAL_METAGRAPH_ALIGN.out.species_hits
                .map { meta, tsv -> [meta.id, count_metagraph_species_hits(tsv, 'metagraph_align')] }

            // map_qc drops samples with nothing above metagraph_align_min_hits, same
            // reasoning as mSWEEP's map_qc above.
            metagraph_align_mapqc_counts_ch = VIRAL_METAGRAPH_ALIGN.out.map_qc
                .map { meta, tsv -> [meta.id, count_metagraph_map_qc(tsv, 'metagraph_align')] }
        } else {
            metagraph_align_counts_ch = Channel.empty()
            metagraph_align_mapqc_counts_ch = Channel.empty()
        }

        // -- Pseudoalignment via Metagraph (metagraph query --query-mode labels), same
        // shared index, alternative method to the alignment above -- see
        // workflows/VIRAL_METAGRAPH_QUERY.nf and modules/metagraph_query.nf for why this
        // exists as its own module rather than the old filter+query pipeline that
        // metagraph_align.nf itself replaced (found ~zero real hits, see git history).
        if (params.run_metagraph_query) {
            VIRAL_METAGRAPH_QUERY(preprocessed_3tuple_ch)

            metagraph_query_counts_ch = VIRAL_METAGRAPH_QUERY.out.species_hits
                .map { meta, tsv -> [meta.id, count_metagraph_species_hits(tsv, 'metagraph_query')] }

            metagraph_query_mapqc_counts_ch = VIRAL_METAGRAPH_QUERY.out.map_qc
                .map { meta, tsv -> [meta.id, count_metagraph_map_qc(tsv, 'metagraph_query')] }
        } else {
            metagraph_query_counts_ch = Channel.empty()
            metagraph_query_mapqc_counts_ch = Channel.empty()
        }

        sequence_index_sample_ch
            .join(msweep_counts_ch, remainder: true)
            .join(mapqc_counts_ch, remainder: true)
            .join(metagraph_align_counts_ch, remainder: true)
            .join(metagraph_align_mapqc_counts_ch, remainder: true)
            .join(metagraph_query_counts_ch, remainder: true)
            .join(metagraph_query_mapqc_counts_ch, remainder: true)
            .map { id, meta, m_counts, qc_counts, mga_counts, mga_qc_counts, mgq_counts, mgq_qc_counts ->
                def new_meta = meta + (m_counts ?: EMPTY_MSWEEP_COUNTS) + (qc_counts ?: EMPTY_MAP_QC_COUNTS) +
                    (mga_counts ?: EMPTY_METAGRAPH_ALIGN_COUNTS) + (mga_qc_counts ?: EMPTY_METAGRAPH_ALIGN_MAPQC_COUNTS) +
                    (mgq_counts ?: EMPTY_METAGRAPH_QUERY_COUNTS) + (mgq_qc_counts ?: EMPTY_METAGRAPH_QUERY_MAPQC_COUNTS)
                [id, new_meta]
            }
            .set { mapping_report_prep_ch }

        GENERATE_MAPPING_REPORT(mapping_report_prep_ch)

        // PUBLISH (mapping lane)
        publish_mapping_lane_json(GENERATE_MAPPING_REPORT.out.publish_seq_level_ch)
        publish_mapping_run_files(GENERATE_MAPPING_REPORT.out.publish_run_level_summaries_ch)
    }

    // PUBLISH
    publish_aln_files(aln_publish_ch)
    publish_nc_files(publish_nextclade_outputs_ch)
    // publish_per_sample_json(per_seq_json_publish_ch)
    publish_per_sample_json(GENERATE_CLASSIFICATION_REPORT.out.publish_seq_level_ch)
    publish_run_files(GENERATE_CLASSIFICATION_REPORT.out.publish_run_level_summaries_ch)

    workflow.onComplete = {
        // Log colors ANSI codes
        /*
        * ANSI escape codes to color output messages
        */

        println """
        Pipeline execution summary
        ---------------------------
        Completed at : ${ANSI_GREEN}${workflow.complete}${ANSI_RESET}
        Duration     : ${ANSI_GREEN}${workflow.duration}${ANSI_RESET}
        Success      : ${workflow.success ? ANSI_GREEN : ANSI_RED}${workflow.success}${ANSI_RESET}
        Results Dir  : ${ANSI_GREEN}${file(params.outdir)}${ANSI_RESET}
        Work Dir     : ${ANSI_GREEN}${workflow.workDir}${ANSI_RESET}
        Exit status  : ${ANSI_GREEN}${workflow.exitStatus}${ANSI_RESET}
        Error report : ${ANSI_GREEN}${workflow.errorReport ?: '-'}${ANSI_RESET}
        """.stripIndent()
    }
}

def __check_if_params_file_exist(param_name, param_value){
    def error = 0

    if (!(param_value==null)){
        def param_file = file(param_value)
        if (!param_file.exists()){
            log.error("${param_file} does not exist")
            error +=1
        }
    }

    if (param_value==null){
        log.error("${param_name} must be provided")
        error +=1
    }
    return error
}

def check_main_params(){

    def errors = 0

    errors += check_sort_reads_params()

    if (errors > 0) {
        log.error("Parameter errors were found, the pipeline will not run.")
        exit 1
    }
}
/* Introspection
 *
 * https://www.nextflow.io/docs/latest/metadata.html
 */

def parse_mnf(mnf) {
    /*
    -----------------------------------------------------------------
    Parses the manifest file to create a channel of metadata and
    FASTQ file pairs.

    Also, checks if there are sample_id duplicated and/or containing
    non alphanumeric characters. Only exception accepted is "_", as
    long as it is not two consecutives "__".

    -----------------------------------------------------------------

    - **Input**:
        mnf (path to the manifest file)

    - **Output**:
        Channel with tuples of metadata and FASTQ file pairs.

    -----------------------------------------------------------------
    */
    // Read manifest file into a list of rows
    def mnf_rows = channel.fromPath(mnf).splitCsv(header: true, sep: ',')

    // Collect sample IDs and validate
    def sample_ids = []
    def errors = 0

    def _errors_ch = mnf_rows.map { row ->
        def sample_id = row.sample_id

        // Check if sample_id is empty
        if (!sample_id) {
            log.error("Empty sample_id detected.")
            errors += 1
        } else {
            // Check for unique sample IDs
            if (sample_ids.contains(sample_id)) {
                log.error("${sample_id} is duplicated")
                errors += 1
            } else {
                sample_ids << sample_id
            }

            // Check if sample_id is alphanumeric, allows underscores but not consecutive
            if (!sample_id.matches(/^(?!.*__)[A-Za-z0-9_]+$/)) {
                log.error("Non alphanumeric sample id ${sample_id} ['_' is permitted]")
                errors += 1
            }
            return errors
        }
        }
        // be sure that the number of errors is evaluated after all rows are processed
        .collect()
        // kill the pipeline if errors are found
        .subscribe{ _v ->
        if (errors > 0) {
            log.error("${errors} critical errors in the manifest were detected. Please check README for more details.")
            exit 1
        }
    }

    // If validation passed, create the channel as before
    def mnf_ch = mnf_rows.map { row ->
                    // set meta
                    def meta = [
                      // id is internal to the pipeline and taxid
                      // is added to it latter
                      id: row.sample_id,
                      // sample_id is explictily used on the
                      // publishing of files paths
                      sample_id: row.sample_id
                    ]
                    // set files
                    def reads = [row.reads_1, row.reads_2]
                    // declare channel shape
                    [meta, reads]
                 }

    return mnf_ch // tuple(meta, [fastq_pairs])
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

// Default counts for a sample a given optional step produced no output for -- either
// because that method didn't run at all (join remainder is null) or because the method
// ran but the step's own output is itself optional per-sample (e.g. nothing above a
// min-abundance/min-hits threshold). Named constants (not inline [:]) so every sample
// still gets the same report columns regardless of which method(s) actually ran for it.
EMPTY_MSWEEP_COUNTS = [msweep_n_groups: 0, msweep_top_group: '', msweep_top_abundance: 0.0]
EMPTY_MAP_QC_COUNTS = [mapqc_n_species: 0, mapqc_max_breadth_pct: 0.0]
// One pair per Metagraph method (align, query) -- both call count_metagraph_species_hits()/
// count_metagraph_map_qc() with a distinct prefix, since both methods' counts can merge
// into the same per-sample meta and would otherwise collide on field name.
EMPTY_METAGRAPH_ALIGN_COUNTS = [metagraph_align_n_species_considered: 0, metagraph_align_n_species_called: 0]
EMPTY_METAGRAPH_ALIGN_MAPQC_COUNTS = [metagraph_align_mapqc_n_species: 0, metagraph_align_mapqc_max_breadth_pct: 0.0]
EMPTY_METAGRAPH_QUERY_COUNTS = [metagraph_query_n_species_considered: 0, metagraph_query_n_species_called: 0]
EMPTY_METAGRAPH_QUERY_MAPQC_COUNTS = [metagraph_query_mapqc_n_species: 0, metagraph_query_mapqc_max_breadth_pct: 0.0]

def empty_metagraph_counts(prefix) {
    return prefix == 'metagraph_align' ? EMPTY_METAGRAPH_ALIGN_COUNTS : EMPTY_METAGRAPH_QUERY_COUNTS
}

def empty_metagraph_mapqc_counts(prefix) {
    return prefix == 'metagraph_align' ? EMPTY_METAGRAPH_ALIGN_MAPQC_COUNTS : EMPTY_METAGRAPH_QUERY_MAPQC_COUNTS
}

def count_msweep_abundances(txt) {
    // <sample>_mSWEEP_abundances.txt: "<label>\t<relative_abundance>", '#'-prefixed and
    // non-numeric-second-column lines skipped -- same rule bin/select_reference_records.py
    // applies in parse_abundances(), so these counts describe the same set of groups the
    // downstream map-QC step actually considered.
    if (txt == null || !txt.exists()) return EMPTY_MSWEEP_COUNTS
    def rows = []
    txt.readLines().each { line ->
        def trimmed = line.trim()
        if (!trimmed || trimmed.startsWith('#')) return
        def cols = trimmed.split('\t')
        if (cols.size() < 2) return
        try {
            rows << [cols[0], cols[1] as Double]
        } catch (NumberFormatException ignored) {
            // header or malformed row -- skipped, as parse_abundances does
        }
    }
    if (!rows) return EMPTY_MSWEEP_COUNTS
    def above = rows.findAll { row -> row[1] >= params.msweep_map_min_abundance }
    def top = rows.max { row -> row[1] }
    return [
        msweep_n_groups:      above.size(),
        msweep_top_group:     top[0],
        msweep_top_abundance: top[1]
    ]
}

def count_metagraph_species_hits(tsv, prefix) {
    // <sample>_species_hits.tsv (bin/call_metagraph_species.py): sample_id, species,
    // hit_count, provisional_call -- the last written by Python's str(bool), so "True"/
    // "False", not lowercase. Shared by both Metagraph methods (see
    // VIRAL_METAGRAPH_ALIGN.nf / VIRAL_METAGRAPH_QUERY.nf, both call
    // CALL_METAGRAPH_SPECIES) -- prefix ('metagraph_align' or 'metagraph_query') keeps
    // their counts from colliding when both merge into the same per-sample meta.
    if (tsv == null || !tsv.exists()) return empty_metagraph_counts(prefix)
    def lines = tsv.readLines()
    if (lines.size() < 2) return empty_metagraph_counts(prefix)
    def header = lines[0].split('\t')
    def called_idx = header.findIndexOf { String col -> col == 'provisional_call' }
    def n_considered = lines.size() - 1
    def n_called = 0
    if (called_idx >= 0) {
        n_called = lines[1..-1].count { String line ->
            def cols = line.split('\t')
            called_idx < cols.size() && cols[called_idx] == 'True'
        }
    }
    return ["${prefix}_n_species_considered": n_considered, "${prefix}_n_species_called": n_called]
}

def count_metagraph_map_qc(tsv, prefix) {
    // <sample>_metagraph_map_qc.tsv (bin/aggregate_metagraph_coverage.py): sample_id,
    // species, hit_count, reference_accession, reference_length, query_length,
    // covered_bases, breadth_pct, mean_depth, meanbaseq, meanmapq, reads_mapped. Same
    // prefix reasoning as count_metagraph_species_hits() above.
    if (tsv == null || !tsv.exists()) return empty_metagraph_mapqc_counts(prefix)
    def lines = tsv.readLines()
    if (lines.size() < 2) return empty_metagraph_mapqc_counts(prefix)
    def header = lines[0].split('\t')
    def breadth_idx = header.findIndexOf { String col -> col == 'breadth_pct' }
    def breadths = lines[1..-1].collect { String line ->
        def cols = line.split('\t')
        if (breadth_idx < 0 || breadth_idx >= cols.size()) return 0.0
        try { return cols[breadth_idx] as Double } catch (NumberFormatException ignored) { return 0.0 }
    }
    return [
        "${prefix}_mapqc_n_species":       lines.size() - 1,
        "${prefix}_mapqc_max_breadth_pct": breadths ? breadths.max() : 0.0
    ]
}

def count_msweep_map_qc(tsv) {
    // <sample>_msweep_map_qc.tsv, written by bin/aggregate_species_coverage.py:
    // sample_id, species_label, relative_abundance, reference_length, query_length,
    // covered_bases, breadth_pct, mean_depth, meanbaseq, meanmapq, reads_mapped
    if (tsv == null || !tsv.exists()) return EMPTY_MAP_QC_COUNTS
    def lines = tsv.readLines()
    if (lines.size() < 2) return EMPTY_MAP_QC_COUNTS
    def header = lines[0].split('\t')
    def breadth_idx = header.findIndexOf { String col -> col == 'breadth_pct' }
    def breadths = lines[1..-1].collect { line ->
        def cols = line.split('\t')
        if (breadth_idx < 0 || breadth_idx >= cols.size()) return 0.0
        try { return cols[breadth_idx] as Double } catch (NumberFormatException ignored) { return 0.0 }
    }
    return [
        mapqc_n_species:       lines.size() - 1,
        mapqc_max_breadth_pct: breadths ? breadths.max() : 0.0
    ]
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
