# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

- **[removed]**: the two map-QC steps `workflows/THEMISTO_MAP_QC.nf` and `workflows/METAGRAPH_MAP_QC.nf`, kept unused since breadth moved downstream to the consensus alignment `subworkflows/mapping.nf` performs anyway. Nothing invoked either, and keeping them meant every comment describing the sequence-index lane had to explain a stage that never ran
- **[removed]**: `modules/select_reference_record_by_name.nf` and `bin/select_reference_record_by_name.py`, unused since reference resolution moved into `subworkflows/mapping.nf`, which reuses the record id the calling method's own index-label map already names
- **[removed]**: four parameters left defined by the above and referenced by nothing — `themisto_map_bowtie_threads`, `metagraph_map_bowtie_threads`, `msweep_map_reference_seed`, `themisto_map_publish_subdir` and, with them, their `nextflow_schema.json` entries. A command still passing one will now be rejected by schema validation
- **[removed]**: five `withName` executor entries in `nextflow.config` for processes no longer reachable from `main.nf` — `SAMTOOLS_COVERAGE`, `AGGREGATE_THEMISTO_COVERAGE`, `GENERATE_THEMISTO_MAP_SUMMARY`, `AGGREGATE_METAGRAPH_COVERAGE` and `GENERATE_METAGRAPH_MAP_SUMMARY`
- **[change]**: the `docs/nf-metro/route_map.mmd` station that matched `.*_MAP_QC` now matches `CALL_.*_SPECIES` and reads "Call species from read hits" — the stage the lane actually runs. `route_map.svg` re-rendered

## [2.0.0]

A major release: viral-lens was one pipeline — Kraken2 classification into consensus
building — and is now four independent lanes over the same preprocessed reads, selected by
`--do_*` switches. The 1.x behaviour is the `--do_mapping` lane and stays **on by default**,
so an existing command still runs and still produces the same consensuses; the version bump
is for the two renamed run-level reports below, and for the restructuring around them.

> ### ⚠️ Breaking: both run-level report files are renamed
>
> | before (1.x) | after (2.0.0) | content |
> | --- | --- | --- |
> | `summary_report.csv` | `consensus_summary_report.csv` | per-consensus results (unchanged content) |
> | `mapping_summary_report.csv` | `sequenceindex_summary_report.csv` | sequence-index lane per-sample results (unchanged content) |
>
> **No report is called `mapping_*` any more**, and that is the point: anything reading
> either 1.x name now fails to find the file instead of quietly reading the wrong one.
> `summary_report.csv` never said what it held, and `mapping_summary_report.csv` was
> actively misleading -- it is written by the sequence-index lane
> (`subworkflows/classifying_index.nf`), not the mapping lane. Each now carries the name of
> what produces it.

- **[added]**: **three new optional lanes**, each with its own master switch, all running off the same preprocessed reads and independent of one another — `--do_sequence_index` (classify against pre-built Themisto2/Metagraph sequence indexes), `--do_assembly` (de novo assembly, viral identification, binning, QC and vContact3 taxonomy) and `--do_abundance` (Kraken2+Bracken, sourmash/inStrain, SCRuB decontamination, mSWEEP). All default off, so a 1.x command runs unchanged
- **[added]**: `--do_mapping` (default `true`), the master switch for what viral-lens was in 1.x: Kraken2 taxid classification plus the consensus/Nextclade/subtyping/classification-report pass. On by default so existing commands are unaffected. With it off the run is reference-free, no Kraken2 database is required, `--call_consensus_for_new_species` is rejected, and the sequence-index lane reports `overlapping_n_species` as `NA`. At least one lane must be enabled or the run is rejected at startup
- **[change]**: `main.nf` no longer holds the pipeline body inline. Each lane is a self-contained subworkflow under [`subworkflows/`](subworkflows/) owning its own report counts and publishing, and `main.nf` only wires preprocessed reads into whichever lanes are enabled. Verified output-identical to the pre-split pipeline on a farm run
- **[change]**: the old end-to-end mapping workflow is split at the consensus boundary — `subworkflows/classifying_kraken2.nf` (Kraken2 + Kraken2Ref taxid selection) and `subworkflows/mapping.nf` (everything from `GENERATE_CONSENSUS` on). `mapping.nf` is driven by **either** classifier, so one consensus/Nextclade/subtyping/report pass covers both rather than each growing a parallel copy. A species only the sequence-index lane found therefore now gets Nextclade, SARS-CoV-2 subtyping and a row in the classification report; it previously had its consensus published on its own with none of that
- **[added]**: `--call_consensus_for_new_species` (default `false`) — hand species that a sequence-index method called but Kraken2 missed to the same shared consensus path. Where both classifiers agree, Kraken2's call and reference win. Because nothing maps these reads before consensus, the breadth threshold deciding whether such a species is reported at all (`--new_species_min_breadth_pct`, default `10.0`) is applied *after* the consensus alignment rather than in the classifier
- **[added]**: `--do_mixed_input` (default `false`) — widens input beyond `--manifest` to a local reads manifest, ENA download and/or iRODS retrieval, merged into the same downstream shape. Existing `--manifest` usage is untouched when off
- **[added]**: **taxonomy whitelist/blacklist on sequence-index species calls** (`--run_taxon_filter`, default on). A call is kept only if a whitelisted taxon appears anywhere in its lineage; `--taxon_filter_whitelist` defaults to the eight respiratory virus families (`10508,10780,11118,2560066,12058,11158,11244,11308`), and `--taxon_filter_blacklist` (default empty) rejects even within an accepted family. Lineages are resolved through `--taxon_filter_table`, the `RVDB_Taxon_Current.tab` shipped with the index, whose per-accession lineage is already flattened — so no taxdump and no parent-walk is involved. A species the table cannot resolve is rejected, and the callers log the resolved fraction so a stale or mismatched table is visible rather than silently rejecting everything. Applies to all three methods (Themisto2, `metagraph_align`, `metagraph_query`)
- **[added]**: **minimum reference length on sequence-index species calls** (`--min_called_reference_length`, default `1000`; `0` disables). RVDB carries partial-CDS and single-gene records alongside complete genomes, and on a clustered index those take read-hits just as readily — but they are not usable as a consensus reference, and they clear the downstream `new_species_min_breadth_pct` gate trivially *because* they are short (breadth is a percentage of reference length, so a 400bp reference needs 40 covered bases to reach 10%). The verdict is on the call's most-hit record with no re-pick, so a reference is never silently swapped. Adds `INDEX_REFERENCE_LENGTHS`, which prices every record once per run per reference FASTA
- **[added]**: seven columns in the per-sample `*_species_hits.tsv`: `taxon_filter`, `taxonomy_id`, `family`, `family_taxon_id`, `reference_record`, `reference_length`, `reference_filter`. Filtered species keep their row with `provisional_call` `False` and the reason recorded, so a rejected call stays readable against its own hit count instead of vanishing from the table. Appended after the existing four columns, and every consumer reads them by name, so nothing shifts
- **[added]**: `*_n_species_taxon_filtered` and `*_n_species_short_reference` per method in `sequenceindex_summary_report.csv`, following the `NA`-vs-`0` convention below: `NA` when the method did not run *or* ran with that gate off, `0` when the gate ran and rejected nothing
- **[added]**: `Discovered_By` column in `consensus_summary_report.csv`, listing **every** method that called the species as a `;`-separated list (`kraken2`, `themisto2`, `metagraph_align`, `metagraph_query`) rather than only the one whose reference won. Previously the report gave no way to tell which classifier found a species, and the corroboration was discarded before it reached either the CSV or the JSON. Appended as the last column so existing column positions do not shift
- **[added]**: `overlapping_n_species` column in `sequenceindex_summary_report.csv` — how many of the species Kraken2 selected were also called by a sequence-index method
- **[added]**: `assembly_sample_summary_report.csv` (one row per sample), `assembly_scaffold_summary_report.csv` (one row per geNomad viral scaffold) and `assembly_vmag_summary_report.csv` (one row per vRhyme bin), replacing `assembly_summary_report.csv`. All three add the taxonomy geNomad and vContact3 assigned; see README
- **[change]**: summary reports now write `NA`, not `0`, for a step that did not run in that execution. `0` now means "ran and found nothing". Affects `themisto_*`, `metagraph_align_*`, `metagraph_query_*`, `new_species_candidates_n`, `bracken_n_species_called` and `msweep_*`. `new_species_candidates_n` in particular read `0` whenever `--call_consensus_for_new_species` was off, which was indistinguishable from the lane genuinely finding nothing new
- **[change]**: `summary_report.csv` now called `consensus_summary_report.csv` (see breaking note above)
- **[change]**: `mapping_summary_report.csv` now called `sequenceindex_summary_report.csv` (see breaking note above)
- **[change]**: mSWEEP moved out of the sequence-index lane into the abundance lane — it estimates abundance rather than calling species. `--run_msweep` now requires `--do_sequence_index true --run_themisto true` and errors up front if they are missing
- **[change]**: shared modules and subworkflows are now included from the `rvi_toolbox` submodule instead of vendored as forked copies, removing ~4700 lines of duplicated Nextflow that had begun to drift from upstream
- **[change]**: per-sample outputs are grouped by lane, and everything — including Nextflow's own execution report, timeline and DAG — is published under `--outdir` rather than `$launchDir/results`
- **[change]**: the pipeline's Nextflow code is now documented per directory — [`subworkflows/README.md`](subworkflows/README.md) for the lanes and [`workflows/README.md`](workflows/README.md) for the steps they are built from — replacing the main README's `Sub-workflows` section, which had drifted (it still documented `COMPUTE_QC_METRICS`, gone; a viral-lens-owned `KRAKEN2BRACKEN` fork, now taken from `rvi_toolbox`; and `ASSEMBLE_META`/`GENOMAD_CLASSIFY`/`VRHYME_BIN`/`CHECKV_QC` as if local). The main README keeps the user-facing material and links out
- **[change]**: the SCRuB decontamination step and the reference-subset module are now taken from `rvi_toolbox` too. `workflows/SCRUB_DECONTAM.nf`, `modules/reformat_bracken.nf`, `bin/reformat_bracken_for_scrub.py` and `modules/reference_subset.nf` are removed. Two viral-lens fixes went upstream first so nothing regresses: SCRuB rows are ordered by the plate map (sorting them independently makes SCRuB fail outright), and the plate map is read as `utf-8-sig` for Excel-exported BOMs. `EXTRACT_REFERENCE_RECORD` moved upstream as well
- **[removed]**: `scrub_zero_species_in_controls`. Zeroing named species in control samples ahead of SCRuB is no longer supported; the parameter and its `--zero-species-in-controls` implementation are gone rather than carried into `rvi_toolbox`
- **[change]**: the vContact3 step is now taken from `rvi_toolbox` rather than kept locally. `workflows/VCONTACT3_RUN.nf`, `modules/vcontact3.nf`, `bin/vcontact3_prep.py` and `bin/vcontact3_postprocess.py` are removed; both scripts were byte-identical to the toolbox's and the subworkflow's interface was the same, so this is a drop-in. The `publishDir saveAs` fix viral-lens was carrying its own copy of the module for went upstream first, so nothing regresses
- **[removed]**: `bin/pool_viral_scaffolds.py`, dead since the binning helpers moved into `rvi_toolbox` — `binning_helper_processes.nf` uses the toolbox's copy
- **[removed]**: `SELECT_REFERENCE_RECORDS` and `bin/select_reference_records.py`, the input step of mSWEEP low-abundance hit validation, from both viral-lens and `rvi_toolbox`. Nothing invoked either; `rvi/rvi_toolbox` master had already dropped both
- **[removed]**: four parameters left defined after map-QC was dropped and referenced by nothing — `themisto_align_run_map_qc`, `metagraph_align_run_map_qc`, `themisto_map_reference_fasta` and `msweep_map_bowtie_threads`. They were kept so an existing `--flag` would still pass schema validation; a command still passing one will now be rejected
- **[removed]**: mSWEEP map-QC (`MSWEEP_MAP_QC` and its per-species coverage roll-up) from `rvi_toolbox`. Nothing invoked it: breadth now comes from the consensus alignment the mapping lane performs anyway, so validating a call by mapping it first meant mapping every real call twice with two different aligners
- **[change]**: mSWEEP no longer writes the per-read probability matrix (`--write-probs`), adopted from `rvi/rvi_toolbox` master. It was tens of GB per sample and nothing consumed it. `<sample>_mSWEEP_probs.tsv` is no longer produced
- **[change]**: the SCRuB heatmap's read-change threshold is now the `scrub_heatmap_min_read_change` parameter (default `20`) instead of a constant inside the R script, adopted from `rvi/rvi_toolbox` master
- **[removed]**: `mapping_pipeline_main.nf`, the copy of the 1.x pipeline kept as a second entry point while `main.nf` was being rebuilt into lanes. `--do_mapping` (on by default) now gives that behaviour from `main.nf` itself, so the copy had become a second definition of one lane, free to drift from the real one. Anyone invoking `nextflow run mapping_pipeline_main.nf` directly should run `main.nf` with the other `--do_*` lanes left off
- **[removed]**: the assembly lane's per-sample `<sample_id>.properties.json`. `assembly_reports.py` never read it, so it was a second derivation of the same counts feeding nothing, and it came off a chain of inner joins that dropped any sample vRhyme never ran for (35 of 95 on a full run). The lane's three CSVs are its report
- **[removed]**: `mapping_run_summary.json` and `abundance_run_summary.json` — each duplicated its CSV exactly (same keys, same records, no nested values), so neither carried anything the CSV did not. `consensus_sequence_properties.json` is **kept**: it holds nested `nextclade_results` and per-position depth that a CSV cannot represent
- **[removed]**: `assembly_run_summary.json` — superseded by the assembly CSVs below. It was also under-reporting: built from a chain of inner joins that dropped any sample vRhyme never ran for, it listed 35 of 95 samples on a full run
- **[fix]**: a sample that produced no vRhyme bins aborted the entire run while publishing (`No signature of method: ScriptBinding.file()`); a missing optional output no longer takes the run down
- **[fix]**: vContact3 post-processing selected query genomes on the wrong separator and so matched none, leaving `final_assignments_postprocessed.csv` and `final_assignments_noveltaxa.csv` empty on every run
- **[fix]**: vContact3 post-processing now fills in the `Proteins` count vContact3 leaves blank for query genomes, without which the novel-genus protein-range check silently never fired

## [1.5.2]

-- **[changed]**: configuration changes to use current queues supported on Sanger HPC when running with `sanger_standard` profile.

## [1.5.1]

-- **[changed]**: updated to version 2.2.1 of Kraken2ref, fixing an edge-case in the polling algorithm that led to valid reference candidates being discarded

## [1.5.0]

- **[added]**: add Nextclade to compute QC metrics (optional)
- **[added]**: added ability to run pre-processing workflow up front (optional)
- **[added]**: consensus workflow now (by default) performs a second round of iVar consensus, after realigning reads to the initial consensus (original behaviour can be achieved with parameters; see README)
- **[added]**: new script `write_all_summaries.py` which produces per-consensus properties files (JSON) and per-run properties files (JSON and CSV). Run-level JSON/CSV contains ONLY "best fit" nextclade dataset outputs
- **[change]**: read alignment now uses minimap2 by default (bwa still available as an option)
- **[change]**: samtools mpileup max depth (-d) can now be controlled with a parameter (default 2000)
- **[change]**: moved QC metrics calling into GENERATE_CONSENSUS; this workflow now emits ONLY data where consensus is not all-N. COMPUTE_QC_METRICS workflow now removed.
- **[change]**: QC JSON, and Nextclade JSON now aggregated at per-consensus level to create single properties file
- **[change]**: `classification_report.csv` now called `summary_report.csv`
- **[change]**: tweaks to sanger_standard execution profile
- **[removed]**: No longer run iVar variants (superseded by nextClade analysis); iVar variants properties removed from `summary_report.csv`

## [1.4.1]

- **[improvement]**: promote kraken report per sample (*.kraken_report.txt) to a primary output
- **[fix]**: restrict filter non-simple types from QC output (their inclusion in meta resulted in unstable behaviour)
- **[fix]**: run_pangolin reverted to local executor in sanger_standard profile (due to instability when run under lsf)

## [1.4.0]

- **[fix]**: correct ivar_min_quality_threshold to use integer value
- **[added]**: rewrite of qc.py to compute more properties and produce properties.json per consensus sequence
- **[added]**: a consensus_sequence_properties.json file is produced for the run
- **[added]**: additional column has been added to classification_report.csv: Total_Mapped_Bases
- **[change]**: classification report columns Percentage_N_bases changes to Percentage_non_N_bases
- **[change]**: reversion to a single classification report and removal of parameters for asserting "invalid" sequences
- **[change]**: filtering out of all consensus sequences (and associated files) that comprise 100% Ns
- **[change]**: changes to workflow and process interfaces for clarity
- **[change]**: sanger_standard execution profile now submits bwa and ivar tasks to LSF

## [1.3.1]

- **[fix]**: additional fix to `developer_publish` sorted reads FASTQ files

## [1.3.0]

- **[fix]**: fix `developer_publish` (did not publish some files previously)

## [1.2.0]

- **[fix]**: remove singularity cachedir pointing to container_dir
- **[fix]**: change default queue from "long" to "normal"
- **[fix]**: refactor workflows to no longer modify meta object in-place; adding safety to the workflows

## [1.1.1]

- **[fix]**: fix sanger local

## [1.1.0]

- **[added]**: add report output with non-valid sequences
- **[removed]**: remove entry points
- **[removed]**: `PREPROCESSING` workflow removed
- **[fix]**: tests don't use `sanger_standard` profile by default
- **[added]**: reimplement ivar variants and add mutation statistics
- **[added]**: `min_coverage_percent` parameter added

## [1.0]

- **[fix]**: Fix unknown config attribute `CONTAINER_DIR` for tests
- **[improvement]**: update samtools to `1.21` for the qc container
- **[fix]**: unmapped reads extracted from flagstat by qc script fixed
- **[improvement]**: run nf-test on local containers support added
- **[fix]**: ivar unit test emtpy consensus fix
- **[improvement]**: update samtools to `1.21` for the base and ivar container
- **[fix]**: Remove parens from reference fasta header to prevent propagation to BAM header

## [0.4.1]

- **[hotfix]**: `outdir` default value set as the same in `nextflow-commons.config`
- **[added]**: `containers_dir` as the default value for singularity `cachedir`

## [0.4.0]

- **[added]**: Implement PAM's nextflow commons retry strategy
- **[added]**: add sanger specific settings on sanger profile
- **[added]**: use quay.io containers by default
- **[added]**: add preprocessing subworkflow (under rvi_toolbox)
- **[added]**: add docker recipes for all containers
- **[improvement]**: new container for ivar, without conda.

## [0.3.2]

- **[added]**: obtain flu B segment number
- **[added]**: `sample_id` collumns now is checked for duplicated, non-alphanumeric or empty row values.

## [0.3.1] - 2024-10-24

- **[hotfix]**: fix k2r dump fq bash syntax and add kraken2 memory request set by a parameter

## [0.3.0] - 2024-10-21

### Changed

- **[improvement]**: Mpileup output retained by run_ivar & used by the QC script for calculating % genome coverage.
- **[improvement]**: Removed unnecessary code from qc.py and run_qc.nf including the plot generation.
- **[improvement]**: Modified qc.py to read input files from command line including using samtools flagstat for read counts.
- **[improvement]**: Unit test files and snapshot files for run_ivar, run_qc_script, and GENERATE_CONSENSUS to account for changes

### Added

- **[added]**: update and added extensive documentation
- **[improvement]**: Update container of kraken2ref from v2.0.0 to v2.1.0
- **[added]**: Add a parameter to set the polling method for kraken2ref (default method set to max)
- **[added]**: Container for the run_qc process
- **[added]**: Unit test for COMPUTE_QC_METRICS workflow
- **[added]**: Mpileup test data
- **[improvement]**: implement k2r release new features
- **[improvement]**: split fastq files if higher than a set numbers of reads per fq

## [0.2.2] - 2024-08-02

### Changed

- **[improvement]**: The ivar module has been updated to adhere to the ARTIC pipeline standards

### Added

- **[added]**: LSF memory escalation strategy for kraken2ref
- **[added]**: Columns Virus_Taxon_ID, Virus, Species, Reference_Taxon_ID, Selected_Reference added/populated to classification report

### Added
- **[added]**: add LSF memory escalation strategy for kraken2ref
- **[added]**: Columns Virus_Taxon_ID, Virus, Species, Reference_Taxon_ID, Selected_Reference added/populated to classification report

## [0.2.1] - 2024-06-20

### Fixed

- **[bug]**: Classification report generation would crash if ' was present in output report file lines
- **[bug]**: Independent workflow stanza for GENERATE_CLASSIFICATION_REPORT.nf was outdated / broken

## [0.2.0] - 2024-05-29

### Fixed

- **[bug]**: Classification report and pre report parsing errors fixed

### Changed

- **[improvement]**: Remove mpileup repeated command calls on ivar process.
- **[improvement]**: Remove redundant processes, rewiring and tiding up code base.
- **[improvement]**: Qc metrics using the same method of the Artic pipeline
- **[improvement]**: add kraken2ref as the new reference selection tool
- **[updated]**: Unit tests adapted to new channel and processes structure

### Added

- **[added]**: A script (`k2r_report.py`) was added to generate a pre report file from k2r software
- **[added]**: unit tests for new `run_kraken2ref_and_pre_report.nf`

## [0.1.0] - 2024-02-05

### Added

- **[added]**: Viral subtyping and classification reports routines integrated to pipeline
- **[added]**: `Percentage Coverage` and `number of mapped reads` are now computed at a new QC metrics workflow
- **[added]**: new workflow (`SUBTYPE_AND_SEGMENT_FLU.nf`) attempts to retrieve the flu subtype and segment from kraken report file and populates the meta with these values
- **[added]**: new module (`retrieve_flu_subtype_and_segment.nf`) attempts to parse out the flu subtype and segment from the kraken report file and sets these values to Null if nothing retrieved
- **[added]**: QC metrics workflow, currently computes reads depth and percentage genome coverage
- **[added]**: SARS-CoV-2 sequences subtyping via pangolin
- **[added]**: branching `GENERATE_CONSENSUS` workflow output for viral subtyping routines
- **[added]**: new parameter (`min_reads_for_taxid`)to set a treshold for minimum number of reads assigned for a taxid to be considered
- **[added]**: new workflow and module (`GENERATE_CLASSIFICATION_REPORT and write_classification_report`) to generate a classifcation report file
- **[added]**: unit tests written in `nf-test` covering modules, workflows and pipeline

### Changed

- **[improvement]**: `taxid` respective `rank` and `name` are available on meta
- **[improvement]**: Taxid reference fasta files for consensus sequence are obtained from kraken database
- **[improvement]**: Channels now rely on Meta Mapping
- **[improvement]**: Output folder now have the following structure `output/<sample_id>/<taxid>`
- **[improvement]**: `write_manifest.py` relies on glob expression

### Fixed

- **[bug]**: Samples with no taxids to process raises a warning and are now ignored

### Removed

- **[Removed]**: writing manifest process removed from `SORT_READS_BY_REF.nf`
- **[Removed]**: json resource files and fasta files provided on the repo

---
## [0.0.1] - 2023-12-01

This is the first prototype versioning. This pipeline 1) sort reads via Kraken and 2) generate consensus sequences using ivar.
