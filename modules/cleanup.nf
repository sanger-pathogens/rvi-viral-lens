process CLEANUP_SORTED_BAM_FILES {
    /**
    * Cleanup intermediate files
    */

    input:
        path(sorted_bam_file)

    script:
        """
        # Remove sorted bam files
        orig_bam=\$(readlink -f ${sorted_bam_file})
        rm ${sorted_bam_file} # Remove symlink
        rm \${orig_bam} # Remove original files
        """
}

process CLEANUP_TRIMMED_FASTQ_FILES {
    /**
    * Cleanup intermediate files
    */

    input:
        tuple val(meta), path(first_read), path(second_read)

    script:
        """
        # Remove trimmed fastq files
        orig_fastqs=\$(readlink -f ${first_read} ${second_read})
        rm ${first_read} ${second_read} # Remove symlink
        rm \${orig_fastqs} # Remove original files
        """
}

process CLEANUP_INSTRAIN_OUTPUT {
    /**
    * Cleanup unused output
    */

    input:
         tuple val(meta), path(workdir)

    script:
        """
        # Remove instrain results
        instrain_dir=\$(cat $workdir)
        cd \$instrain_dir
        rm -rf ${meta.id}*
        """
}

process CLEANUP_THEMISTO_PSEUDOALIGNMENTS {
    /**
    * Cleanup intermediate files
    */

    input:
        tuple val(meta), path(pseudoalignment_1), path(pseudoalignment_2)

    script:
        """
        # Remove themisto pseudoalignment files
        orig_files=\$(readlink -f ${pseudoalignment_1} ${pseudoalignment_2})
        rm ${pseudoalignment_1} ${pseudoalignment_2} # Remove symlink
        rm \${orig_files} # Remove original files
        """
}
