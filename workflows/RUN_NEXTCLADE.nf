#!/usr/bin/env nextflow

// import workflow
nextflow.enable.dsl = 2
include {publish_consensus_files} from '../modules/publish_lite.nf'
include {run_nextclade; collate_nextclade_jsons} from '../modules/run_nextclade.nf'

workflow RUN_NEXTCLADE {
    take:
    input_ch // [meta, fa]

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
    //data_dir_ch.not_found.view{"No nextclade data on nextclade_data_dir found for ${it[0].id}"}

    data_dir_ch.found.map{meta, fa, ref_dirs_map ->
        def dir_list = ref_dirs_map.collect { it.dir }
        [meta, fa, dir_list]
    }
    .set { nextclade_In_ch }
    run_nextclade(nextclade_In_ch)
    collate_nextclade_jsons(run_nextclade.out)

    emit:
    collate_nextclade_jsons.out
}


List<File> findReferenceDirs(nc_index_file, virusTaxid, segNumber = null) {
    def nc_index = new groovy.json.JsonSlurper().parse(new File(nc_index_file))

    if (!segNumber) {
        segNumber = "ALL"
    }

    def data_dirs = []
    def virus_data = nc_index.get(virusTaxid)
    if (virus_data != null) {
        data_dirs = virus_data.get(segNumber)
    } else {
        log.warn("No nextclade data for $virusTaxid")
    }

    if (!data_dirs) {
        log.warn("No nextclade data for $virusTaxid, $segNumber")
        return []
    }

    def output = data_dirs
                    .collect { new File(it).canonicalFile }
                    .sort { it.path }
    return output
}

List<Map> buildReferenceTags(dataDir, sampleId, species_taxid, subtype = null, segNumber = null) {
    /**
    * Build tag IDs for all reference dirs found under the Nextclade hierarchy.
    *
    * Case 1: <species_taxid>/<assembly>/reference.fasta
    *   tag_id = <sample_id>.<species_taxid>.<assembly>
    *
    * Case 2: <species_taxid>/<subtype>/<seg_number>/<assembly>/reference.fasta
    *   tag_id = <sample_id>.<species_taxid>.<subtype>.segment<seg_number>.<assembly>
    *
    * @return List of [File referenceDir, String tagId] pairs
    */
    def dirs = findReferenceDirs(params.nextclade_index_json, species_taxid, segNumber)

    if (dirs == []){
        return []
    } else {
        return dirs.collect { dir ->
            def assembly = dir.name  // directory itself is the assembly name
            def tagId
            if (subtype && segNumber) {
                // Case 2
                tagId = "${sampleId}.${species_taxid}.${subtype}.segment${segNumber}.${assembly}"
            } else {
                // Case 1
                tagId = "${sampleId}.${species_taxid}.${assembly}"
            }
            [dir: dir.toString(), tag_id: tagId]
        }
    }
}
