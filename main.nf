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
include {publish_run_files} from './modules/publish_lite.nf'

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
  Viral Lens [v1.5.0]
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
    --ivar_initial_freq_threshold : ${params.ivar_initial_min_depth}
    --ivar_polish_min_depth       : ${params.ivar_polish_min_depth}
    --ivar_polish_freq_threshold  : ${params.ivar_polish_freq_threshold}
    
  --> viral subtyping branching parameters:
    --scv2_keyword             : ${params.scv2_keyword}

  --> Nextclade parameters:
    --nextclade_index_json      : ${params.nextclade_index_json}

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
    if (params.do_preprocessing) {
        reads_ch.map{ meta, fastqs ->
            return [meta, fastqs[0], fastqs[1]]
        }.set{preproc_in_ch}

        PREPROCESSING(preproc_in_ch).out_ch.map{meta, read1, read2 ->
            return [meta, [read1, read2]]
        }.set{sort_reads_in_ch}

    } else {
        sort_reads_in_ch = reads_ch
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
