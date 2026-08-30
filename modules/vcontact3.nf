// Sample-namespace and assemble a single sample's vContact3-shaped inputs.
// Proteins come from two sources:
//   - vRhyme bin .faa files (one per bin) — each bin becomes one "genome"
//   - geNomad's virus_proteins.faa — for scaffolds vRhyme didn't bin (those
//     become singleton "genomes"), AND for binned scaffolds CheckV rates
//     High-/Medium-quality (those become an ADDITIONAL standalone "genome"
//     on top of their bin, giving vContact3 both signals)
// See rvi_toolbox/bin/vcontact3_prep.py for the exact semantics.
process VCONTACT3_PER_SAMPLE_PREP {
    tag "${meta.id}"
    label 'cpu_1'
    label 'mem_2'
    label 'time_1'

    // The Sanger vcontact3 image doesn't expose python3 on PATH (conda env only
    // activated for the vcontact3 entrypoint). This script only needs stdlib,
    // so use the minimal Sanger pysam image rather than pulling in a heavier
    // biocontainer just for python3.
    container 'quay.io/sangerpathogens/pysam:0.0.2'

    input:
    tuple val(meta), path(virus_proteins), path(virus_summary), path(bin_membership), path(bins_fasta_dir, stageAs: 'bins_fasta'), path(checkv_quality)
    path(prep_script)

    output:
    tuple val(meta), path("${meta.id}_proteins.faa"),    emit: proteins
    tuple val(meta), path("${meta.id}_gene2genome.tsv"), emit: gene2genome
    tuple val(meta), path("${meta.id}_lengths.tsv"),     emit: lengths

    script:
    """
    python3 ${prep_script} \\
        --sample-id        ${meta.id} \\
        --proteins         ${virus_proteins} \\
        --summary          ${virus_summary} \\
        --bin-membership   ${bin_membership} \\
        --bins-fasta       bins_fasta \\
        --checkv-quality   ${checkv_quality} \\
        --out-proteins     ${meta.id}_proteins.faa \\
        --out-gene2genome  ${meta.id}_gene2genome.tsv \\
        --out-lengths      ${meta.id}_lengths.tsv
    """
}

// Pipeline-level vContact3 run. Takes the collected per-sample prep outputs
// (one file per sample), concatenates them into a single proteins.faa and
// gene2genome.tsv, then runs `vcontact3 run`. vcontact3 3.2.0 emits result
// CSVs into `vcontact3_out/exports/`, so the publishDir patterns and output
// declarations target that nested path.
process VCONTACT3 {
    label 'cpu_16'
    label 'mem_8'
    label 'time_12'

    container 'quay.io/sangerpathogens/vcontact3:3.2.0'

    publishDir "${params.results_dir}/vcontact3", mode: 'copy', overwrite: true, pattern: "vcontact3_out/exports/final_assignments.csv",   saveAs: { f -> file(f).getName() }
    publishDir "${params.results_dir}/vcontact3", mode: 'copy', overwrite: true, pattern: "vcontact3_out/exports/performance_metrics.csv", saveAs: { f -> file(f).getName() }

    input:
    path(all_proteins_faa)
    path(all_gene2genome_tsv)
    path(all_lengths_tsv)

    output:
    path 'vcontact3_out/exports/final_assignments.csv',   emit: final_assignments,   optional: true
    path 'vcontact3_out/exports/performance_metrics.csv', emit: performance_metrics, optional: true
    path 'vcontact3_out/HMMprofile_vog_support.h5',       emit: vog_support,         optional: true

    script:
    """
    # Concatenate per-sample protein FASTAs
    cat ${all_proteins_faa} > combined_proteins.faa

    # Concatenate per-sample tables: keep the header from the first file in each
    # set, drop the headers from the rest.
    concat_with_header () {
        local out=\$1; shift
        local first; first=\$(ls -1 "\$@" | head -1)
        head -n 1 "\$first" > "\$out"
        for f in "\$@"; do
            tail -n +2 "\$f" >> "\$out"
        done
    }

    concat_with_header combined_gene2genome.tsv ${all_gene2genome_tsv}
    concat_with_header combined_lengths.tsv     ${all_lengths_tsv}

    # vcontact3 calls pd.read_parquet directly on --gene2genome (markers.py:249)
    # and on --len-nucleotide (profiles.py:158), despite docs claiming both
    # accept TSV. Convert both combined TSVs to parquet before invoking it.
    # Container's own conda Python has pandas + pyarrow (vcontact3 deps).
    /opt/conda/bin/python - <<'PYEOF'
import pandas as pd
pd.read_csv("combined_gene2genome.tsv", sep="\t") \\
    .to_parquet("combined_gene2genome.parquet")
pd.read_csv("combined_lengths.tsv", sep="\t") \\
    .to_parquet("combined_lengths.parquet")
PYEOF

    vcontact3 run \\
        --proteins        combined_proteins.faa \\
        --gene2genome     combined_gene2genome.parquet \\
        --len-nucleotide  combined_lengths.parquet \\
        --db-path         ${params.vcontact3_db_path} \\
        --db-version      ${params.vcontact3_db_version} \\
        --db-domain       ${params.vcontact3_db_domain} \\
        --output          vcontact3_out \\
        -t                ${task.cpus}
    """
}

// Post-process VCONTACT3's final_assignments.csv: keep query genomes only,
// flag out-of-protein-range novel-genus calls as uncertain, fix realm/lower-rank
// taxonomy conflicts using each genome's own VOG evidence, and split off the
// remaining (still-confident) novel-genus calls into their own file.
// --genome-report points at the same reference database VCONTACT3 itself used
// (params.vcontact3_db_path/-version), same convention as --db-path/-version
// above: referenced directly by path, not staged, since it lives in the shared
// reference database directory. Note the "v" prefix on the version directory:
// vcontact3 resolves --db-path/--db-version to <db_path>/v<version>/, and the
// genome report lives in there beside the rest of that version's files.
// See rvi_toolbox/bin/vcontact3_postprocess.py for the exact semantics.
process VCONTACT3_POSTPROCESS {
    label 'cpu_1'
    label 'mem_4'
    label 'time_1'

    container 'quay.io/sangerpathogens/vcontact3:3.2.0'

    publishDir "${params.results_dir}/vcontact3", mode: 'copy', overwrite: true, pattern: "final_assignments_postprocessed.csv"
    publishDir "${params.results_dir}/vcontact3", mode: 'copy', overwrite: true, pattern: "final_assignments_noveltaxa.csv"
    publishDir "${params.results_dir}/vcontact3", mode: 'copy', overwrite: true, pattern: "postprocess_report.txt"

    input:
    path(final_assignments)
    path(vog_support)
    path(postprocess_script)

    output:
    path 'final_assignments_postprocessed.csv', emit: postprocessed
    path 'final_assignments_noveltaxa.csv',     emit: novel_taxa
    path 'postprocess_report.txt',              emit: report

    script:
    """
    /opt/conda/bin/python ${postprocess_script} \\
        --assignments       ${final_assignments} \\
        --vog-support       ${vog_support} \\
        --genome-report     ${params.vcontact3_db_path}/v${params.vcontact3_db_version}/RefSeq.${params.vcontact3_db_version}.genome_report.parquet \\
        --output            final_assignments_postprocessed.csv \\
        --novel-taxa-output final_assignments_noveltaxa.csv \\
        --report            postprocess_report.txt \\
        --min-proteins      ${params.vcontact3_postprocess_min_proteins} \\
        --max-proteins      ${params.vcontact3_postprocess_max_proteins}
    """
}
