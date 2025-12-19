// Copyright (C) 2023 Genome Surveillance Unit/Genome Research Ltd.

include {run_kraken} from '../modules/run_kraken.nf'
include {get_taxid_reference_files} from '../modules/get_taxid_references.nf'
include {run_k2r_sort_reads; run_k2r_dump_fastqs_and_pre_report} from '../modules/run_kraken2ref_and_pre_report.nf'

workflow SORT_READS_BY_REF {
    /*
    -----------------------------------------------------------------
    Sort Reads to A Given Reference

    The SORT_READS_BY_REF workflow processes paired-end sequencing
    reads by sorting them according to taxonomic classifications
    obtained from Kraken2. This workflow uses a manifest file to
    process multiple samples and produces sorted by taxid FASTQ files
    for each sample and classification reports.

    -----------------------------------------------------------------
    # Inputs

    - **Manifest File**: A CSV file containing sample metadata and
    paths to paired-end FASTQ files.
    - **Kraken Database Path**: Path to the Kraken database.
    - **Kraken2Ref Library Fasta**: Optional. Path to the Kraken2Ref
    library fasta file. If none is provided, it assumes there
    is a `${params.db_path}/library/library.fna`.

    -----------------------------------------------------------------
    # Key Processes

    1. **Run Kraken**: Classifies reads using Kraken2 against a
    specified database.

    2. **Sort Reads**: Uses Kraken2Ref to sort reads by taxonomic ID
    and extracts them into separate FASTQ files.

    3. **Merge FASTQ Parts**: Merges split FASTQ parts if necessary.

    4. **Collect Reference Files**: Retrieves reference sequences
    based on taxonomic IDs for downstream analysis.

    -----------------------------------------------------------------
    # Outputs
    - `sample_taxid_ch`: Channel containing tuples of metadata and
    sorted reads per taxonomic ID.
    - `sample_pre_report_ch`: Channel containing pre-reports with
    sample-level summaries.

    */

    take:
        mnf_ch // tuple(meta, [reads_1, reads_2])
    main:

        // drop channel tuples where 1 or more FASTQ files are empty
        filtered_mnf_ch = mnf_ch.branch {
            empty: file(it[1][0]).size() < 5000 || file(it[1][1]).size() < 5000
                log.warn("Empty fastq file(s) for ${it[0].id}")
            not_empty: true
        }.not_empty
        // -------------------------------------------------//

        // 1 - run kraken and get outputs
        run_kraken(filtered_mnf_ch, params.db_path)
        // -------------------------------------------------//

        // 2 - obtain unique taxid reference fasta file and metadata from kraken db

        // check if library fasta was set,
        // if not, assume is at library/library.fna on the kraken db dir
        if (params.db_library_fa_path == null){
            library_fa_path = "${params.db_path}/library/library.fna"
            log.warn("No db_library_fa_path set, assuming ${library_fa_path} exists")
        } else {
            library_fa_path = params.db_library_fa_path
        }
        // -------------------------------------------------//

        // 3 - run Kraken2Ref

        // 3.1 - sort reads
        // drop unclassified fq filepair and store kraken2 files on meta
        sort_reads_in_ch = run_kraken.out 
            .map { it -> // tuple (meta, kraken_output, [classified_fq_filepair], [unclassified_fq_filepair], kraken_report)
                tuple(it[0], it[1], it[2], it[4]) 
            }
            .map {  meta, kraken_output, classified_fq_filepair, kraken_report ->  
                // get input file sizes to estimate resources usage (cpu and mem)
                def fq_1_size = classified_fq_filepair[0].size() // byte
                def fq_2_size = classified_fq_filepair[1].size() // byte
                def _new_meta = meta.plus([
                    classified_fq_filepair: classified_fq_filepair, // store classified fq filepair on meta
                    // store kraken2 outputs on meta to simplify channel gymnastics
                    kraken_report: kraken_report,
                    kraken_output: kraken_output,
                    fqs_total_size: fq_1_size + fq_2_size, // store total size of fq files on meta
                ])

                //_new_meta.fqs_total_size = fq_1_size + fq_2_size
                return [_new_meta, kraken_output, kraken_report]
            }
     
        // 3.2 run sort reads (meta, kraken_report, kraken_output)
        run_k2r_sort_reads(sort_reads_in_ch)

        // 3.3 - prepare chanel for k2r dump fqs
        k2r_dump_fq_In_ch = run_k2r_sort_reads.out.json_files 
            .map { meta, tax_to_reads_json, decomposed_json ->
                [ meta, meta.classified_fq_filepair, tax_to_reads_json, decomposed_json, meta.kraken_report]
            }
              
        run_k2r_dump_fastqs_and_pre_report(k2r_dump_fq_In_ch)

        // -------------------------------------------------//

        // 4 - run kraken2ref and collect all per-sample per-taxon fq filepairs
        raw_sample_pre_report_ch = run_k2r_dump_fastqs_and_pre_report.out.report_file
        per_taxid_fqs_Ch = run_k2r_dump_fastqs_and_pre_report.out.fq_files 

        // 4.1 - remove empty pre_report files
        sample_pre_report_ch = raw_sample_pre_report_ch
            .filter{ it -> (it.size() > 1) } // remove empty pre_reports
            .splitCsv(header: true, sep:"\t")

        // raise warning for sample taxids which had empty pre_reports
        raw_sample_pre_report_ch
            .filter{it -> (it.size() <= 1)}
            .view{it -> log.warn("Excluding ${it} as input due to small size ( < 1 byte)")}

        // prepare channel to be emitted
        per_sample_taxid_ch = per_taxid_fqs_Ch
            .map {_meta, reads ->
                // group pairs of fastqs based on file names, and add new info to meta
                reads
                    .groupBy { filePath -> filePath.getName().tokenize("_")[0..-2].join(".")}
                    .collect { identifier, paths ->[[identifier, paths]]}
                // this map returns a channel with a single value: a list with all the ids and files.
            }
            .flatten() // flat the list [id_a, fqa1, fqa2, idb, fqb1, fqb2, ...]
            .collate(3) // tuple (id.taxid, fq_1, fq_2)
            .map { it ->
                // rebuild meta and reads structure
                def meta = [sample_id:it[0].tokenize(".")[0..-2].join("_"), //run.lane.sample_info
                        taxid:it[0].tokenize(".")[-1]] //taxid
                meta.id = "${meta.sample_id}.${meta.taxid}"
                def reads = [it[1], it[2]]
                [meta, reads] // new meta object created for taxid FASTQ pair - tuple(meta, [fq_1, fq_2])
            }

        // 5 - get references files
        unq_taxid_ch = per_sample_taxid_ch
            .map{meta, _reads -> [meta.taxid]}
            .collect()
            .flatten()
            .unique()

        get_taxid_reference_files(unq_taxid_ch, library_fa_path)
        get_taxid_reference_files.out.set{taxid_ref_files_map_ch}

        sample_taxid_ch = per_sample_taxid_ch
            .map {meta, reads ->
                [meta.taxid, meta, reads]
            }
            .combine(taxid_ref_files_map_ch, by:0) // [taxid, meta, reads, ref_files, ref_header]
            .map {_taxid, meta, reads, ref_fa, ref_header ->
                def new_meta = meta.plus([reference_header: ref_header])
                [new_meta, reads, ref_fa]
            }

    emit:
        sample_taxid_ch // tuple (meta, reads, ref_files)
        sample_pre_report_ch // pre_report

}

def check_sort_reads_params(){
    /*
    -----------------------------------------------------------------
    Checks for necessary parameters and validates paths to ensure
    they exist. Logs errors if any required parameters are missing.
    -----------------------------------------------------------------

    - **Output**: Number of errors encountered during the checks.

    -----------------------------------------------------------------

    */
    def errors = 0
    // was the kraken database provided?
    if (params.db_path == null){
        log.error("No kraken database path provided")
        errors +=1
    }

    // TODO: if provided, check if it is a valid path

    // was the manifest provided?
    if (params.manifest == null){
        log.error("No manifest provided")
        errors +=1
    }
    // if yes, is it a file which exists?
    if (params.manifest){
        def manifest_file = file(params.manifest)
        if (!manifest_file.exists()){
            log.error("The manifest provided (${params.manifest}) does not exist.")
            errors += 1
        }
        //TODO: validate manifest
    }

    return errors
}
