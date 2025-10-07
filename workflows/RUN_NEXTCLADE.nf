#!/usr/bin/env nextflow

// import workflow
nextflow.enable.dsl = 2
include {run_nextclade} from '../modules/run_nextclade.nf'
include {publish_consensus_files} from '../modules/publish_lite.nf'

workflow RUN_NEXTCLADE {
    take:
    // meta, fa
    input_ch // tuple(Sample_ID, Virus_Taxon_ID, Flu_Segment, Reference_Subtype)

    main:

    input_ch
    .map{meta, consensus_fa ->
        def ref_dirs_map = buildReferenceTags(
                params.nextclade_data_dir,
                meta.sample_id,
                meta.virus,
                meta.sample_subtype,
                meta.flu_segment)
        return [meta, consensus_fa, ref_dirs_map] // [[meta, fa, [[ref_dir_a], [ref_dir_b], ...]]
    }
    .branch{ meta, fa, ref_dirs_map ->
        found : ref_dirs_map != []
            [meta, fa, ref_dirs_map]

        not_found: ref_dirs_map == []
            [meta, fa]
    }
    .set{ data_dir_ch}

    //data_dir_ch.not_found.view{"No nextclade data on nextclade_data_dir found for ${it[0].id} ${it}"}

    data_dir_ch.found
    .flatMap{list -> expandList(list)} // [[meta, fa, ref_dir_map_a], [meta, fa, ref_dir_map_b]]
    .map{meta, fa, ref_dirs_map ->
        def new_meta = meta + ref_dirs_map

        [new_meta, fa, new_meta.dir]
    }
    .set {nextclade_In_ch}

    // publish nextclade output files
    run_nextclade(nextclade_In_ch, params.nextclade_output_verbosity)
    run_nextclade.out
    .map{meta, csv, tar_gz -> [meta, [csv,tar_gz]] }
    .set {publish_files_In_ch}

    publish_consensus_files{publish_files_In_ch}

}


def expandList(List input) {
    if (input.isEmpty()) return []

    def prefix = input[0..-2] // everything except the last element
    def last   = input[-1] // last element (may be a list or not)

    if (last instanceof List) {
        // Expand each element of the last list with the prefix
        return last.collect { item -> prefix + item }
    } else {
        // If last isn't a list, just wrap the input
        return [input]
    }
}


List<File> findReferenceDirs(dataDir, virusTaxid,subtype = null, segNumber = null) {
    /**
    * Find directories that contain a `reference.fasta` exactly one level below a base path.
    *
    * Case 1: base = <dataDir>/<virusTaxid>
    * Case 2: base = <dataDir>/<virusTaxid>/<subtype>/<segNumber>
    *
    * Returns a List<File> of the matching directories (sorted), or [] if none.
    */
    // Build the base path depending on whether subtype/segNumber are provided
    def parts = [dataDir, virusTaxid]

    if (segNumber != null) {
        parts << segNumber
    }

    if (subtype != null) {
        parts << subtype
    }

    def baseDir = new File(parts.join(File.separator))

    if (!baseDir.isDirectory()) {
        return []
    }

    // Look only at immediate subdirectories under baseDir,
    // and keep those that contain a file named "reference.fasta"
    def output = (baseDir.listFiles()) //?: [])
        .findAll { it.isDirectory() && new File(it, "reference.fasta").isFile() }
        .collect { it.canonicalFile }   // normalize paths
        .sort { it.path }               // deterministic order
    return output
}

List<Map> buildReferenceTags(dataDir, sampleId,taxid,subtype = null, segNumber = null) {
    /**
    * Build tag IDs for all reference dirs found under the Nextclade hierarchy.
    *
    * Case 1: <taxid>/<assembly>/reference.fasta
    *   tag_id = <sample_id>.<taxid>.<assembly>
    *
    * Case 2: <taxid>/<subtype>/<seg_number>/<assembly>/reference.fasta
    *   tag_id = <sample_id>.<taxid>.<subtype>.segment<seg_number>.<assembly>
    *
    * @return List of [File referenceDir, String tagId] pairs
    */
    def dirs = findReferenceDirs(dataDir, taxid, subtype, segNumber)

    if (dirs == []){
        return []
    } else {
        return dirs.collect { dir ->
            def assembly = dir.name  // directory itself is the assembly name
            def tagId
            if (subtype && segNumber) {
                // Case 2
                tagId = "${sampleId}.${taxid}.${subtype}.segment${segNumber}.${assembly}"
            } else {
                // Case 1
                tagId = "${sampleId}.${taxid}.${assembly}"
            }
            [dir: dir.toString(), tag_id: tagId]
        }
    }
}
