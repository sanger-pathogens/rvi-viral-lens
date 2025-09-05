#!/usr/bin/env nextflow

// import workflow
nextflow.enable.dsl = 2
include {run_nextclade} from '../modules/run_nextclade.nf'

workflow RUN_NEXTCLADE {
    take:
    input_ch // tuple(Sample_ID, Virus_Taxon_ID, Flu_Segment, Reference_Subtype)

    main:
    split_ch = input_ch.branch {meta, it ->
        def (data_dir, tag_id) = check_data_dir(it)
        meta.data_dir = data_dir
        meta.tag_id = tag_id
        found: data_dir != false
        not_found: data_dir == false
    }
    if (split_ch.found.count() == 0) {
        log.warn "No data found for samples in given manifest!"
        exit 0
    }
    //TODO: explore using buffer for this process
    run_nextclade(split_ch.found).set { out_ch } //.buffer(size: params.buffer_size, remainder: true)

    channel.of().set { out_ch }
    emit:
    out_ch
}

def _count_dirs(in_path) {
    def dir = new File(in_path);
    def dirs = [];
    dir.traverse(type: groovy.io.FileType.DIRECTORIES, maxDepth: 1) { d ->
        dirs.add(d)
    };
    return dirs.size
}

def check_data_dir(input_ch_tuple) {
    try {
        def (sample_id, virus_taxid, flu_seg, ref_subtype) = input_ch_tuple
        def path_to_check = "${params.nextclade_data_dir}/${virus_taxid}"

        if ( file(path_to_check).isDirectory() ) {
            def subdirs = _count_dirs(path_to_check)
            def next_path_to_check = "${path_to_check}/${ref_subtype}"
            if (subdirs > 0 && !file(next_path_to_check).isDirectory()) {
                return [false, null]
            } else if (subdirs > 0 && file(next_path_to_check).isDirectory()) {
                def last_subdirs = _count_dirs(next_path_to_check)
                def last_path_to_check = "${path_to_check}/${ref_subtype}/${flu_seg}"
                if (last_subdirs == 0) {
                    def tag_id = "${sample_id}.${virus_taxid}.${ref_subtype}"
                    return [next_path_to_check, tag_id]
                } else if (last_subdirs > 0 && !file(last_path_to_check).isDirectory()) {
                    return [false, null]
                } else if (last_subdirs > 0 && file(last_path_to_check).isDirectory()) {
                    def tag_id = "${sample_id}.${virus_taxid}.${ref_subtype}.segment${flu_seg}"
                    return [last_path_to_check, tag_id]
                }
            } else if (subdirs == 0) {
                def tag_id = "${sample_id}.${virus_taxid}"
                return [path_to_check, tag_id]
            }
        } else {return [false, null]}

    } catch (Exception e) {
        log.warn "Error checking data directory: ${e.message}"
        return [false, null]
    }
}
