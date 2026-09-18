#!/usr/bin/env nextflow
// Copyright (C) 2023 Genome Surveillance Ltd/Genome Research Ltd.

// enable dsl2
nextflow.enable.dsl = 2

include {check_sort_reads_params} from './workflows/SORT_READS_BY_REF.nf'
include {validateParameters; paramsSummaryLog} from 'plugin/nf-schema'
include {PREPROCESSING} from "./rvi_toolbox/subworkflows/preprocessing.nf"
include {MIXED_INPUT} from "./rvi_toolbox/subworkflows/mixed_input.nf"
include {CLASSIFYING_KRAKEN2} from './subworkflows/classifying_kraken2.nf'
include {CLASSIFYING_INDEX} from './subworkflows/classifying_index.nf'
include {MAPPING} from './subworkflows/mapping.nf'
include {ASSEMBLY} from './subworkflows/assembly.nf'
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
  Viral Lens [v2.0.0-beta]
  Used parameters:
  -------------------------------------------
  --> general pipeline parameters:
    --outdir                   : ${params.outdir}
    --do_mixed_input           : ${params.do_mixed_input}
    --do_preprocessing         : ${params.do_preprocessing}
    --do_assembly              : ${params.do_assembly}
    --do_mapping               : ${params.do_mapping}
    --do_sequence_index        : ${params.do_sequence_index}
    --do_abundance             : ${params.do_abundance}
    --default_error_strategy   : ${params.default_error_strategy}
    --max_attempts             : ${params.max_attempts}

  --> CLASSIFYING_KRAKEN2 + MAPPING workflow parameters (subworkflows/classifying_kraken2.nf: taxid classification; subworkflows/mapping.nf: consensus, Nextclade, SCOV2 subtyping, classification report):
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

  --> CLASSIFYING_INDEX workflow parameters (subworkflows/classifying_index.nf; only used if --do_sequence_index true):
    --run_themisto                : ${params.run_themisto}
    --run_metagraph_align         : ${params.run_metagraph_align}
    --run_metagraph_query         : ${params.run_metagraph_query}
    --themisto_align_min_hits     : ${params.themisto_align_min_hits}
    --msweep_themisto_index       : ${params.msweep_themisto_index}
    --msweep_ref_groups           : ${params.msweep_ref_groups}
    --msweep_map_min_abundance    : ${params.msweep_map_min_abundance}
    --metagraph_align_graph       : ${params.metagraph_align_graph}
    --metagraph_align_annotation  : ${params.metagraph_align_annotation}
    --metagraph_align_min_hits    : ${params.metagraph_align_min_hits}

  --> Sequence-index call gates (applied by all three methods' species callers):
    --run_taxon_filter            : ${params.run_taxon_filter}
    --taxon_filter_table          : ${params.taxon_filter_table}
    --taxon_filter_whitelist      : ${params.taxon_filter_whitelist}
    --taxon_filter_blacklist      : ${params.taxon_filter_blacklist ?: '(none)'}
    --min_called_reference_length : ${params.min_called_reference_length}

  --> New-species consensus (subworkflows/mapping.nf; sequence-index calls Kraken2 missed):
    --call_consensus_for_new_species : ${params.call_consensus_for_new_species}
    --new_species_min_breadth_pct    : ${params.new_species_min_breadth_pct}
    --msweep_map_reference_fasta     : ${params.msweep_map_reference_fasta}
    --metagraph_map_reference_fasta  : ${params.metagraph_map_reference_fasta}

  --> ABUNDANCE workflow parameters (subworkflows/abundance.nf; only used if --do_abundance true):
    --run_kraken2bracken             : ${params.run_kraken2bracken}
    --run_abundance_estimation       : ${params.run_abundance_estimation}
    --run_scrub                      : ${params.run_scrub}
    --run_msweep                     : ${params.run_msweep} (needs --do_sequence_index + --run_themisto)
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
    if (params.do_mixed_input) {
        MIXED_INPUT()

        reads_ch = MIXED_INPUT.out.all_reads_ready_ch
            .map { meta, r1, r2 ->
                // meta.sample_id > meta.id 
                def new_meta = meta + [sample_id: meta.id]
                [new_meta, [r1, r2]]
            }
    } else {
        reads_ch = parse_mnf(params.manifest) // tuple(meta, [fastq_1, fastq_2])
    }

    // === Preprocessing ===
    // preprocessed_3tuple_ch (meta, read1, read2) is the single shared input 
    if (params.do_preprocessing) {
        reads_ch.map{ meta, fastqs ->
            return [meta, fastqs[0], fastqs[1]]
        }.set{preproc_in_ch}

        PREPROCESSING(preproc_in_ch)
        PREPROCESSING.out.out_ch.set{ preprocessed_3tuple_ch }

    } else {
        reads_ch.map{ meta, fastqs ->
            return [meta, file(fastqs[0]), file(fastqs[1])]
        }.set{ preprocessed_3tuple_ch }
    }

    // ==========================
    // === 2 - Classify reads by Kraken2
    // Kraken2 classification and the consensus pass below are the two halves of one
    // pipeline, so --do_mapping switches both. Outputs are hoisted into variables because
    // a subworkflow that was never invoked has no .out at all -- reaching for
    // CLASSIFYING_KRAKEN2.out with the flag off aborts the run.
    if (params.do_mapping) {
        CLASSIFYING_KRAKEN2(preprocessed_3tuple_ch)
        kraken2_sample_taxid_ch = CLASSIFYING_KRAKEN2.out.sample_taxid_ch
        kraken2_report_ch       = CLASSIFYING_KRAKEN2.out.sample_report_with_join_key_ch
        identified_species_ch   = CLASSIFYING_KRAKEN2.out.identified_species_ch
    } else {
        kraken2_sample_taxid_ch = Channel.empty()
        kraken2_report_ch       = Channel.empty()
        identified_species_ch   = Channel.empty()
    }

    // === 3 - De novo assembly + viral binning ===
    if (params.do_assembly) {
        ASSEMBLY(preprocessed_3tuple_ch)
    }

    // === 4 - Classify reads against sequence indexes ===
    if (params.do_sequence_index) {
        CLASSIFYING_INDEX(preprocessed_3tuple_ch, identified_species_ch)
        index_species_calls_ch     = CLASSIFYING_INDEX.out.species_calls_ch
        index_called_species_ch    = CLASSIFYING_INDEX.out.called_species_ch
        // Handover for the abundance lane's optional mSWEEP abundances estimation
        themisto_pseudoaln_ch      = CLASSIFYING_INDEX.out.themisto_pseudoalignments
        themisto_ref_groups_ch     = CLASSIFYING_INDEX.out.themisto_ref_groups
    } else {
        index_species_calls_ch     = Channel.empty()
        index_called_species_ch    = Channel.empty()
        themisto_pseudoaln_ch      = Channel.empty()
        themisto_ref_groups_ch     = Channel.empty()
    }

    // === 5 - Consensus, lineage calling and classification report ===
    if (params.do_mapping) {
        MAPPING(
            kraken2_sample_taxid_ch,
            kraken2_report_ch,
            index_species_calls_ch,
            identified_species_ch,
            index_called_species_ch,
            preprocessed_3tuple_ch
        )
    }

    // === 6 - Abundance estimation ===
    if (params.do_abundance) {
        ABUNDANCE(preprocessed_3tuple_ch, themisto_pseudoaln_ch, themisto_ref_groups_ch)
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

    // Called whatever --do_mapping is: it validates the manifest, which every lane reads,
    // as well as the kraken2 database, which only the mapping lane reads. The db check is
    // the one conditioned on --do_mapping, inside the function itself.
    errors += check_sort_reads_params()

    // Every lane off would preprocess the reads and then stop, producing no result at all.
    // Caught here rather than left to finish "successfully" with an empty outdir.
    if (!(params.do_mapping || params.do_sequence_index || params.do_assembly || params.do_abundance)) {
        log.error("No lane is enabled: --do_mapping, --do_sequence_index, --do_assembly " +
                  "and --do_abundance are all false, so the run would produce nothing. " +
                  "Enable at least one.")
        errors += 1
    }

    // The new-species consensus path IS the mapping lane -- MAPPING is what resolves those
    // references, maps them and applies the breadth gate. With --do_mapping false the
    // sequence-index lane still reports its calls; there is just nothing to hand them to.
    if (params.call_consensus_for_new_species && !params.do_mapping) {
        log.error("--call_consensus_for_new_species needs --do_mapping true: consensus for " +
                  "sequence-index species is built by the MAPPING lane, which is switched off.")
        errors += 1
    }

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
