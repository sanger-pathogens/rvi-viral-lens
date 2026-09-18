# Changelog
## [2.0.0-beta]
This is a major release as we bring in three subworkflows that was initially developed for rvi-viralmetagenomics pipeline: de novo assembly + binning (assembly.nf), read classification with sequence indexes (classifying_index.nf) and abundance estimation (abundance.nf). Functionality of viral lens 1.5.2 is wrapped in classifying_kraken2.nf and mapping.nf subworkflows.

- **[breaking-change]**: summary_report.csv is renamed as consensus_summary_report.csv
- **[breaking-change]**: per-sample outputs of the 1.5.2 path now sit under a `mapping/` subdirectory, so that the three new lanes can publish alongside them without collision. Consensus output moves from `<outdir>/<sample_id>/<taxid>/` to `<outdir>/<sample_id>/mapping/<taxid>/`, and the Kraken2 report from `<outdir>/<sample_id>/` to `<outdir>/<sample_id>/mapping/`. Anything reading those paths needs updating; `--outdir` itself still defaults to `$launchDir/results`

- **[added]**: `--do_assembly` (default `false`) — de novo assembly and viral binning: metaSPAdes, geNomad viral identification, vRhyme binning, CheckV QC and vContact3 taxonomy. Reports per sample, per viral scaffold and per vMAG in `assembly_sample_summary_report.csv`, `assembly_scaffold_summary_report.csv` and `assembly_vmag_summary_report.csv`
- **[added]**: `--do_sequence_index` (default `false`) — read classification against pre-built sequence indexes, by Themisto2 pseudoalignment (`--run_themisto`, on by default within the lane), `metagraph align` (`--run_metagraph_align`) and `metagraph query` (`--run_metagraph_query`). Species are called from read-hit counts; per-sample calls land in `*_species_hits.tsv` and the run-level report is `sequenceindex_summary_report.csv`
- **[added]**: `--do_abundance` (default `false`) — abundance estimation: Kraken2+Bracken, sourmash/inStrain, SCRuB cross-contamination decontamination (`--run_scrub`, needs a `--scrub_plate_map`) and mSWEEP (`--run_msweep`, which needs `--do_sequence_index true --run_themisto true` and errors up front if they are missing)
- **[added]**: `--do_mapping` (default `true`), the master switch for the 1.5.2 behaviour. On by default, so an existing command runs unchanged. With it off the run is reference-free and no Kraken2 database is required. At least one lane must be enabled or the run is rejected at startup
- **[added]**: `--do_mixed_input` (default `false`) — widens input beyond `--manifest` to a local reads manifest, ENA download and/or iRODS retrieval, merged into the same downstream shape. Existing `--manifest` usage is untouched when off
- **[added]**: `--call_consensus_for_new_species` (default `false`) — species that a sequence-index method called but Kraken2 missed are handed to the same consensus path, so they get Nextclade, SARS-CoV-2 subtyping and a row in the classification report rather than a consensus published on its own. Where both classifiers agree, Kraken2's call and reference win. Because nothing maps these reads before consensus, the breadth threshold deciding whether such a species is reported at all (`--new_species_min_breadth_pct`, default `10.0`) is applied *after* the consensus alignment rather than in the classifier. Requires `--do_mapping`
- **[added]**: two call gates on sequence-index species calls, both folded into the one `provisional_call` verdict. `--run_taxon_filter` (default on) keeps a call only if a whitelisted taxon appears in its lineage — `--taxon_filter_whitelist` defaults to the eight respiratory virus families (`10508,10780,11118,2560066,12058,11158,11244,11308`), resolved through the `RVDB_Taxon_Current.tab` shipped with the index — and `--min_called_reference_length` (default `1000`, `0` disables) rejects calls whose most-hit record is too short to serve as a consensus reference. A rejected species keeps its row with the reason recorded, so it stays readable against its own hit count
- **[added]**: `Discovered_By` column in `consensus_summary_report.csv`, listing **every** method that called the species as a `;`-separated list (`kraken2`, `themisto2`, `metagraph_align`, `metagraph_query`) rather than only the one whose reference won. Appended as the last column, so existing column positions do not shift
- **[added]**: Nextflow's own execution report, timeline, trace and DAG are now written for every run, under `<outdir>/pipeline_info/`

- **[change]**: `main.nf` no longer holds the pipeline body inline. Each lane is a self-contained subworkflow under [`subworkflows/`](subworkflows/) owning its own report counts and publishing, and `main.nf` only wires preprocessed reads into whichever lanes are enabled. Verified output-identical to 1.5.2 on a farm run
- **[change]**: the 1.5.2 end-to-end path is split at the consensus boundary — `subworkflows/classifying_kraken2.nf` (Kraken2 + Kraken2Ref taxid selection) and `subworkflows/mapping.nf` (everything from `GENERATE_CONSENSUS` on). `mapping.nf` is driven by **either** classifier, so one consensus/Nextclade/subtyping/report pass covers both rather than each growing a parallel copy
- **[change]**: summary reports write `NA`, not `0`, for a step that did not run in that execution; `0` means "ran and found nothing". This matters most for the new lanes' columns, which are `NA` whenever their lane is switched off
- **[change]**: the pipeline's Nextflow code is documented per directory — [`subworkflows/README.md`](subworkflows/README.md) for the lanes and [`workflows/README.md`](workflows/README.md) for the steps they are built from — replacing the main README's `Sub-workflows` section, which had drifted. The main README keeps the user-facing material and links out, and `docs/nf-metro/route_map.svg` draws the lanes as they now run

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
