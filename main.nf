#!/usr/bin/env nextflow
// Copyright (C) 2023 Genome Surveillance Ltd/Genome Research Ltd.

// enable dsl2
nextflow.enable.dsl = 2

// --- import modules ---------------------------------------------------------
include {check_sort_reads_params} from './workflows/SORT_READS_BY_REF.nf'
include {validateParameters; paramsSummaryLog} from 'plugin/nf-schema'

include {PREPROCESSING} from "./rvi_toolbox/subworkflows/preprocessing.nf"
// Widened input handling (rvi_integration_1) -- already lives in viral-lens's own
// rvi_toolbox (rvi/rvi_toolbox.git), unlike every other new lane this integration; no
// fork/port needed.
include {MIXED_INPUT} from "./rvi_toolbox/subworkflows/mixed_input.nf"

// --- rvi_integration_1 lanes -------------------------------------------------
// Each lane is its own subworkflow under subworkflows/, self-contained (including
// its own sample-level report-count helpers and PUBLISH calls). main.nf just wires
// preprocessed reads into whichever lanes are enabled.
include {MAPPING} from './subworkflows/mapping.nf'
include {ASSEMBLY} from './subworkflows/assembly.nf'
include {SEQUENCE_INDEX} from './subworkflows/sequence_index.nf'
include {ABUNDANCE} from './subworkflows/abundance.nf'

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
    --do_mixed_input           : ${params.do_mixed_input}
    --do_preprocessing         : ${params.do_preprocessing}
    --do_assembly              : ${params.do_assembly}
    --do_sequence_index        : ${params.do_sequence_index}
    --do_abundance             : ${params.do_abundance}
    --default_error_strategy   : ${params.default_error_strategy}
    --max_attempts             : ${params.max_attempts}

  --> MAPPING workflow parameters (subworkflows/mapping.nf: taxid mapping, consensus, Nextclade, SCOV2 subtyping):
    --manifest                    : ${params.manifest}
    --db_path                     : ${params.db_path}
    --db_library_fa_path          : ${params.db_library_fa_path}
    --min_reads_for_taxid         : ${params.min_reads_for_taxid}
    --k2r_max_total_reads_per_fq  : ${params.k2r_max_total_reads_per_fq}
    --mem_k2r_b0_offset           : ${params.mem_k2r_b0_offset}
    --mem_k2r_b0                  : ${params.mem_k2r_b0}
    --mem_k2r_b0_final            : ${params.mem_k2r_b0_final}
    --mem_k2r_b1                  : ${params.mem_k2r_b1}
    --mem_k2r_f1                  : ${params.mem_k2r_f1}
    --mem_k2r_a2                  : ${params.mem_k2r_a2}
    --do_consensus_polishing      : ${params.do_consensus_polishing }
    --read_aligner                : ${params.read_aligner}
    --read_aligner_params         : ${params.read_aligner_params}
    --mpileup_max_depth           : ${params.mpileup_max_depth }
    --ivar_initial_min_depth      : ${params.ivar_initial_min_depth}
    --ivar_initial_freq_threshold : ${params.ivar_initial_freq_threshold}
    --ivar_polish_min_depth       : ${params.ivar_polish_min_depth}
    --ivar_polish_freq_threshold  : ${params.ivar_polish_freq_threshold}
    --scv2_keyword                : ${params.scv2_keyword}
    --nextclade_index_json        : ${params.nextclade_index_json}

  --> ASSEMBLY workflow parameters (subworkflows/assembly.nf; only used if --do_assembly true):
    --genomad_db                  : ${params.genomad_db}
    --checkv_db                   : ${params.checkv_db}
    --vcontact3_db_path           : ${params.vcontact3_db_path}
    --vcontact3_db_version        : ${params.vcontact3_db_version}
    --metaspades_subsample_limit  : ${params.metaspades_subsample_limit}
    --vrhyme_min_scaffold_length  : ${params.vrhyme_min_scaffold_length}

  --> SEQUENCE_INDEX workflow parameters (subworkflows/sequence_index.nf; only used if --do_sequence_index true):
    --run_msweep                  : ${params.run_msweep}
    --run_metagraph_align         : ${params.run_metagraph_align}
    --run_metagraph_query         : ${params.run_metagraph_query}
    --msweep_themisto_index       : ${params.msweep_themisto_index}
    --msweep_ref_groups           : ${params.msweep_ref_groups}
    --msweep_map_min_abundance    : ${params.msweep_map_min_abundance}
    --metagraph_align_graph       : ${params.metagraph_align_graph}
    --metagraph_align_annotation  : ${params.metagraph_align_annotation}
    --metagraph_align_min_hits    : ${params.metagraph_align_min_hits}

  --> ABUNDANCE workflow parameters (subworkflows/abundance.nf; only used if --do_abundance true):
    --run_kraken2bracken             : ${params.run_kraken2bracken}
    --run_abundance_estimation       : ${params.run_abundance_estimation}
    --run_scrub                      : ${params.run_scrub}
    --kraken2bracken_kraken2_db      : ${params.kraken2bracken_kraken2_db}
    --kraken2bracken_classification_level : ${params.kraken2bracken_classification_level}
    --scrub_plate_map                : ${params.scrub_plate_map}
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
    // Widened input handling (rvi_integration_1, opt-in via --do_mixed_input): existing
    // --manifest usage (parse_mnf()) is completely unchanged when this is off (the
    // default). When on, MIXED_INPUT merges a local reads manifest (its OWN id/R1/R2
    // columns -- not parse_mnf()'s sample_id/reads_1/reads_2), ENA download, and/or iRODS
    // retrieval into the same downstream shape.
    if (params.do_mixed_input) {
        MIXED_INPUT()

        reads_ch = MIXED_INPUT.out.all_reads_ready_ch
            .map { meta, r1, r2 ->
                // MIXED_INPUT's sources (INPUT_CHECK/ENA_DOWNLOAD/DOWNLOAD_FROM_IRODS,
                // rvi_toolbox/subworkflows/{input_check,ena_input,irods}.nf) only ever set
                // meta.id; everything downstream (publishDir paths, report columns) keys
                // off meta.sample_id.
                def new_meta = meta + [sample_id: meta.id]
                [new_meta, [r1, r2]]
            }
    } else {
        reads_ch = parse_mnf(params.manifest) // tuple(meta, [fastq_1, fastq_2])
    }

    // === Preprocessing ===
    // preprocessed_3tuple_ch (meta, read1, read2) is the single shared input every
    // lane below (mapping, assembly, sequence-index, abundance) consumes.
    if (params.do_preprocessing) {
        reads_ch.map{ meta, fastqs ->
            return [meta, fastqs[0], fastqs[1]]
        }.set{preproc_in_ch}

        PREPROCESSING(preproc_in_ch)
        PREPROCESSING.out.out_ch.set{ preprocessed_3tuple_ch }

    } else {
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
    // === 2 - Map to taxid, generate consensus, classify (see subworkflows/mapping.nf)
    MAPPING(preprocessed_3tuple_ch)

    // === 3 - De novo assembly + viral binning (rvi_integration_1, opt-in) ===
    if (params.do_assembly) {
        ASSEMBLY(preprocessed_3tuple_ch)
    }

    // === 4 - Map reads to sequence indexes (rvi_integration_1, opt-in) ===
    // MAPPING.out.identified_species_ch is always available (MAPPING runs
    // unconditionally above) -- SEQUENCE_INDEX uses it to tell which of its own species
    // calls are genuinely new (see subworkflows/mapping.nf / sequence_index.nf).
    if (params.do_sequence_index) {
        SEQUENCE_INDEX(preprocessed_3tuple_ch, MAPPING.out.identified_species_ch)
    }

    // === 5 - Abundance estimation (rvi_integration_1, opt-in) ===
    if (params.do_abundance) {
        ABUNDANCE(preprocessed_3tuple_ch)
    }

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
