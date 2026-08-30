include { SUBSAMPLE_ITER } from './subsample.nf'
include { METASPADES } from '../modules/metaspades.nf'

workflow ASSEMBLE_META {
    // run metaspades on data set with at most X reads to avoid loss of quality when too much data fed in (statistical inconsistency)
    take:
    meta_reads // tuple val(meta), path(R1), path(R2)

    main:
    metaspades_subsample_limit_ch = Channel.value( params.metaspades_subsample_limit ) 
    
    meta_reads.map{ meta, R1, R2 -> 
        def readCount = R1.countFastq() //as we made these files already they have been verified once
        [meta, R1, R2, readCount]
    }
    | set{ready_for_subsampling}

    SUBSAMPLE_ITER(ready_for_subsampling, metaspades_subsample_limit_ch)
    | METASPADES


    emit:
    contigs_channel = METASPADES.out.contigs_channel
}