[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.10.1-23aa62.svg)](https://www.nextflow.io/) [![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/) ![Nf-test](https://img.shields.io/badge/NFtest-%E2%89%A50.8.4-23aa62.svg?labelColor=0000)

# Viral Lens

**Viral Lens** (also known as viral-lens) is a bioinformatics pipeline for reconstructing and classifying viral genomes from short-read sequencing data. It was developed for the [RVI project](https://www.sanger.ac.uk/group/respiratory-virus-and-microbiome-initiative/) at the Wellcome Sanger Institute, and has been validated using data generated using the RVI target enrichment protocol.

---
## Contents
- [Pipeline Summary](#pipeline-summary)
- [How to Cite](#how-to-cite)
- [Basic usage](#basic-usage)
- [Installation and dependencies](#installation-and-dependencies)
  - [Software](#software)
  - [Containers](#containers)
- [Inputs](#inputs)
  - [Manifest](#manifest)
  - [Kraken2 Database](#kraken2-database)
  - [NextClade config](#nextclade-index-json)
- [Outputs](#outputs)
  - [Primary outputs](#primary-outputs)
  - [Secondary outputs](#secondary-outputs)
  - [De novo assembly + viral binning outputs](#de-novo-assembly--viral-binning-outputs)
- [Configuration](#configuration)
  - [Parameters](#parameters)
  - [Parameter switchboard](#parameter-switchboard)
  - [Profiles](#profiles)
- [Unit Tests](#unit-tests)
- [Pipeline components documentation](#pipeline-components-documentation)
  - [Processes](#processes)
  - [Lanes and steps](#lanes-and-steps)
- [Licence](#licence)

---

## Pipeline Summary

The pipeline takes as input (a) a manifest containing  **fastq pairs file** paths and (b) a **kraken detabase**  (see [Inputs section](#inputs) for more details) and outputs a collection of sequences reconstucted from the reads (using an alignment-and-consensus approach; see below), along with a number of reports. The broad steps followed by the pipeline are as follows:

0. **Preprocessing** (Optional): An optional preprocessing workflow (activated by `--do_preprocessing true`). This remove adapters (via `trimmomatic`), tandem repeats (via `TRF`) and human reads (via `sra-human-scrubber`) from the input fastq files. Each of those steps can be set on/off (`--run_trimmomatic`, `--run_trf`, `--run_hrr`).

1. **Classify Reads and select references** (`subworkflows/classifying_kraken2.nf`): `kraken2` is initially used to classify the reads in the input fastq, using the input Kraken database. The resulting Kraken2 report is used to select partition the reads into groups, each associated with a selected reference sequence that will be used to guide the reconstruction of the viral genome. 

2. **Generate Consensus** (`subworkflows/mapping.nf`): The reads sets produced in the previous step are aligned to their respective references (via `bwa`
or minimap2), with the resulting pileup being provided to `ivar` to determine the sequence by consensus (in either one or two rounds). This step is deliberately *shared*: reads can reach it from either classifier — Kraken2 (step 1) or, when `--call_consensus_for_new_species true`, the sequence-index methods for species Kraken2 missed (step 5). Both hand over the same channel shapes, so steps 2-4 run once over the union rather than once per classifier.

3. **NextClade analysis**: (Optional) NextClade is run on the resulting viral genomes.

4. **Pangolin analysis**: (Optional) For SARS-CoV-2 genomes, Pangolin is run to sub-type the genome. 

5. **Classify reads against sequence indexes** (Optional, activated by `--do_sequence_index true`) (`subworkflows/classifying_index.nf`): a second, independent classifier lane running alongside step 1. Reads are pseudoaligned against a prebuilt sequence index — Themisto2 (the lane's default, `--run_themisto`, on unless turned off) and/or Metagraph (`--run_metagraph_align`, `--run_metagraph_query`) — and species are called directly from read-hit counts, with no probabilistic model. Each call must then clear two further gates — a **taxonomy whitelist** (`--taxon_filter_whitelist`, by default the eight respiratory virus families) and a **minimum reference length** (`--min_called_reference_length`, default 1000bp) — see [Sequence-index call gates](#sequence-index-call-gates-taxonomy-and-reference-length). It does not replace Kraken2: where both classifiers find the same species, Kraken2's call and reference win. Species **only** this lane found are handed to step 2 for consensus when `--call_consensus_for_new_species true` (default `false`); with the default the lane reports its calls and nothing more. Outputs `sequenceindex_summary_report.csv` plus per-sample hit tables.

6. **De novo assembly + viral binning** (Optional, activated by `--do_assembly true`): the same preprocessed reads are also assembled de novo (`metaSPAdes`), classified for viral content (`geNomad`), binned into putative genomes (`vRhyme`), quality-checked (`CheckV`) and clustered/taxonomically assigned (`vContact3`) — a parallel lane alongside steps 1-4, not a replacement for them.

7. **Abundance estimation** (Optional, activated by `--do_abundance true`) (`subworkflows/abundance.nf`): Kraken2 + Bracken species abundance (`--run_kraken2bracken`), sourmash/inStrain genome-level profiling (`--run_abundance_estimation`), SCRuB cross-contamination decontamination (`--run_scrub`), and mSWEEP probabilistic abundance (`--run_msweep`). mSWEEP is the one step here that does not start from the reads: it consumes Themisto2's pseudoalignments, so it additionally requires step 5 (`--do_sequence_index true --run_themisto true`).

The diagram below (rendered with [nf-metro](https://github.com/seqeralabs/nf-metro) from [`docs/nf-metro/route_map.mmd`](docs/nf-metro/route_map.mmd)) shows how these lanes relate. Each optional lane is its own line, and they run in parallel over the same preprocessed reads:

| line | lane | steps |
| --- | --- | --- |
| grey | every sample (input + preprocessing) | 0 |
| green | classify by Kraken2 taxid, then consensus | 1-4 |
| orange | classify against sequence indexes | 5 |
| purple | de novo assembly + viral binning | 6 |
| blue | abundance estimation | 7 |


![viral-lens route map](docs/nf-metro/route_map.svg)

[**(&uarr;)**](#contents)
---

## How to Cite

This software will be published soon. Until it is, please provide the URL to this GitHub repository when you use the software in your own work.

[**(&uarr;)**](#contents)
---

## Basic usage

```bash
PIPELINE_CODE=<path to viral lens repo>
MANIFEST=<path to my manifest>
KRAKEN_DB_PATH=<path to my kraken DB>
PIPELINE_CONTAINERS=<path to my containers dir>
NEXTCLADE_INDEX_JSON=<path to my nextclade_index.json>

## nextclade_index_json is optional -- required if Nextclade output is required
nextflow run ${PIPELINE_CODE}/main.nf \
    --manifest ${MANIFEST} \
    --db_path ${KRAKEN_DB_PATH} \
    --nextclade_index_json ${NEXTCLADE_INDEX_JSON} \
    --outdir ./output/ \
    -with-trace -with-report -with-timeline \
    -profile sanger_local \
    -resume
```

[**(&uarr;)**](#contents)

---
## Installation and dependencies

### Software

The following software in required to run the pipeline (with versions used for testing and validation listed):

- Minimal requirements
  - Nextflow = `23.10.1`
  - Python3
    - pytest = `6.2.2`
    - importlib-resources = `5.1.0`
    - flake8 = `7.0.0`
    - pandas = `2.1.4`
    - cached-property = `1.5.2`
    - scipy = `1.12.0`
  - kraken2 = `v2.1.3`
  - Samtools/htslib = `1.21`
  - BWA = `0.7.17`
  - minimap2 = `2.30`
  - iVar = `1.4.3`
  - [kraken2ref](https://github.com/genomic-surveillance/kraken2ref) = `v2.1.0`
- If requiring NextClade analysis
  - NextClade (CLI) = `3.16.0`
- If requiring Pangolin analysis
  - Pangolin `4.3.1`
- When using containers
  - [Singularity](https://docs.sylabs.io/guides/latest/user-guide/) is required to use Singularity Containers, tested on ``ce version 3.11.4``

> NOTE: Running the pipeline under an environment which has the above software installed but with different versions may work, but has not been tested and validated.

### Containers

Out-of-the-box, the pipeline uses containers for these dependencies. Custom container images required by the pipeline have been deposited at quay.io at [quay.io gsu-pipeline](https://quay.io/organization/gsu-pipelines), and the pipeline is configured to use these containers. 

Source code for these can be found in `"$projectDir/containers"`. To build the containers, run the commands bellow.

```bash
cd containers/
sudo singularity build alignment.sif alignmentContainer.sing
sudo singularity build ivar.sif ivarContainer.sing
sudo singularity build pangolin.sif pangolinContainer.sing
sudo singularity build kraken2ref.sif kraken2ref.sing
```

> NOTE: To use local containers on the pipeline set the parameters, you will need to edit `"$projectDir/conf/containers.config"` and change the locations for each container. 

[**(&uarr;)**](#contents)


## Inputs

This pipeline relies on three **main inputs**:

- **`manifest`** : CSV Manifest of input fastq file pairs.
  - Must have `sample_id`,`reads_1` and `reads_2` columns
  - If you have your set of fastq pairs in a single dir, a script (`write_manifest.py`) is provided to facilitate this process.
- **`db_path`** : Path of a valid [kraken2 database](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#kraken-2-databases)
- **`nextclade_index_json`** : (Optional, if wanting NextClade analysis) Path to a JSON file mapping sequence products from the pipeline to NextClade datasets

### Manifest

The pipeline require as input a manifest containing a unique sample id (`sample_id`) and paths to each of the fastq pair file (`reads_1` and `reads_3`)

```csv
sample_id,reads_1,reads_2
sample1,/path/to/output/sample1_R1.fq,/path/to/output/sample1_R2.fq
sample2,/path/to/output/sample2_R1.fq,/path/to/output/sample2_R2.fq
sample3,/path/to/output/sample3_R1.fq,/path/to/output/sample3_R2.fq
```

> NOTE: The collumn **`sample_id`** must be **unique, alphanumerics (non consecutive "_" are accepted) and cannot be empty**. Pipeline will fail if any of these conditions are not met.

### Kraken2 Database

This pipeline essentially works with a Kraken2 database that has been built from Viral RefSeq in the standard way (e.g. see [here](https://benlangmead.github.io/aws-indexes/k2) for publicly available Kraken3 DBs). However, for respiratory viruses, we have improved the performance of the pipeline using a Kraken2 database that has been prepared in a specific way, namely:

- The Influenza taxonomy is modified below the level of the species such that each segment is represented as its own distinct branch. 
- Further, for specifically segments 4 (HA) and 6 (NA), those branches of the custom taxonomy have an additional level so that nodes representing subtypes H1, H2, H3... are directly below the segment 4 node, and similarly, nodes representing subtypes N1, N2, N3... fall directly under the segment 6 node. Essentially, the segments of the flu genome are considered as distinct sequences for the purposes of the kraken2 classification step
- The number of Influenza, RSV (A and B), Rhinovirus and Metapneumovirus in the database is also expanded as compared to those found in viral RefSeq

An example database prepared in this way can be found [here](https://rvi_kraken2_dbs.cog.sanger.ac.uk/refseq_ncbiFlu_kfv2_20241027.tar.gz). 

A complimentary tool to viral-lens, [vl-kraken-prep](https://github.com/genomic-surveillance/krakenflu), can be used to prepare a Kraken2 database for viral-lens. 

### NextClade index JSON

The index JSON file maps informs the pipeline of the locations of relevent NextClade datasets for the viral genomes it has reconstructed. It is organised by NCBI taxonomy id and segment id (where "ALL" is used as the segment id for monopartite viruses). Here is an example containing configuration for Influenza B (tax id 2955465) and Gammapapillomavirus 11 (tax id 1513256)
```json
{
    "2955465": {
        "1": [
            "/path/to/nextclade_data/data/nextstrain/flu/vic/pb1"
        ],
        "2": [
            "/path/to/nextclade_data/data/nextstrain/flu/vic/pb2"
        ],
        "3": [
            "/path/to/nextclade_data/data/nextstrain/flu/vic/pa"
        ],
        "4": [
            "/path/to/nextclade_data/data/nextstrain/flu/vic/ha/KX058884"
        ],
        "5": [
            "/path/to/nextclade_data/data/nextstrain/flu/vic/np"
        ],
        "6": [
            "/path/to/nextclade_data/data/nextstrain/flu/vic/na/CY073894"
        ],
        "7": [
            "/path/to/nextclade_data/data/nextstrain/flu/vic/mp"
        ],
        "8": [
            "/path/to/nextclade_data/data/nextstrain/flu/vic/ns"
        ]
    },
    "1513256": {
        "ALL": [
            "/path/to/local_data/1513256/GCF_000898855.1",
            "/path/to/local_data/1513256/GCF_000896435.1",
            "/path/to/local_data/1513256/GCF_000908695.1",
            "/path/to/local_data/1513256/GCF_000899695.1"
        ]
    }
}
```

Note in this example that (a) the 8 segments of Influenza B each point to a different data set, and (b) Gammapapillomavirus 11 has 4 datasets (the pipeline will run NextClade against all 4 in this case).

A script is bundled with the pipeline to help with the preparation of this file. Usage:

```bash
python <path/to/viral-lens>/bin/create_index.py /path/to/cloned/github/nextclade_data/data /path/to/local/custom/nextclade/datasets nextstrain,enpen
```

The third positional argument refers to the subdirectories under `path/to/cloned/github/nextclade_data/data` which should be included in the index JSON file.

The expected structure of `/path/to/local/custom/nextclade/datasets` is as follows:
```bash
path/to/local/custom/nextclade/datasets
|--- virus_species_taxID
      |--- assembly_ID
            ├── genome_annotation.gff3
            ├── pathogen.json
            └── reference.fasta
```

[**(&uarr;)**](#contents)

## Outputs

### Primary outputs

The output file tree should look like the tree bellow:

```bash
<output_dir>/
├── mapping_summary_report.csv
├── consensus_sequence_properties.json
├── <sample_id>
│   ├── <sample_id>.kraken_report.txt
│   ├── <ref_id>
│   │   ├── <sample_id>.<ref_id>.consensus.fa
│   │   ├── <sample_id>.<ref_id>.properties.json
│   │   ├── <sample_id>.<ref_id>.nextclade.tar.gz
│   │   ├── <sample_id>.<ref_id>.sorted.bam
│   │   └── <sample_id>.<ref_id>.sorted.bam.bai
│   ├── [...]
├── [...]
```

...where <sample_id> is the identifier provided in the manifest, and <ref_id> is the taxonomy ID of the reference sequence used to create the consensus sequence (note: this may often not correspond to a real NCBI taxonomy id, because the custom kraken2 database used by viral-lens will usually contain many "artificial" nodes introduced to represent segments of multi-partite viruses; see above)

#### <sample_id>/<ref_id>/<sample_id>.<ref_id>.consensus.fa

Inferred consensus sequence for reference `ref_id` in sample `sample_id`
 
#### <sample_id>/<ref_id>/<sample_id>.<ref_id>.properties.json

A collection of observed and computed properties for the inferred consensus sequence. 

- `id` 
  - A unique identifier for the consensus sequence across the run (formed from sample id and reference id)  
  - Example: `50213_1_67.8120647`
- `sample_id` 
  - ID of the sample (as provided in the manifest)
  - Example: `50213_1_67`
- `tax_id`  
  - ID of the reference used to build the consensus sequence
  - Example: `8120647`
- `selected_taxid`
  - Tax ID (in the custom database) of the reference used to construct the consensus
  - Example: `8120647`
- `ref_selected` 
  - Description of the reference used to construct the consensus
  - Example: `"A/swine/Guangxi/NS2394/2012(H3N2) segment 4"`
- `reference_length` 
  - Length of the selected reference
- `virus_subtype`
  - Overall subtype of the selected reference
  - Example: `"H3N2"`
- `virus_name`
  - NCBI species name of the selected reference
  - Example: `Alphainfluenzavirus influenzae`
- `report_name`
  - "Common" name for the virus species for reporting
- `virus` 
  - NCBI taxonomy ID of the selected virus (species level node)   
  - Example: `2955291"`
- `sample_subtype`
  - (Where sub-typing has been possible): Inferred sub-type of the sequence
  - Example: `"H3"`
- `flu_segment`
  - (For segmented / multi-partite viruses) Inferred segment number of the consensus sequence
  - Example: `4`
- `longest_non_n_subsequence`
  - Length of the longest stretch of non-N sequence in the consensus 
- `num_reads`
  - Number of read associated with the reference (by kraken2ref) 
- `reads_mapped` 
  - Number of reads successfully mapped back to the consensus
- `reads_unmapped`
  - Number of reads unmapped
- `bases_mapped`
  - Number of bases mapped mapped back to the consensus
- `reads_mapped_in_proper_pairs` 
  - Reads mapped in proper pairs (expected orientation and distance)  
- `positions_exceeding_depth`
  - Histogram containing the number of positions exceeding depths from 0 to 100 
- `percent_positions_exceeding_depth_10` 
  - Percentage of position exceeding depth 10 
- `percent_non_n_bases`
  - Percentage of bases in the final consensus that are non-M
- `mean_depth_per_position`  
  - Total number of mapped positions divided by length of consensus
- `consensus_length`  
  - Length of final consensus sequence
- `nc.selected_dataset`
  - The NextClade dataset used for the nc.qc.* properties. If NextClade was run on multiple datasets, this is the dataset that resulted in the lowest overall score
  - Example: `nextclade_data/data/nextstrain/flu/h3n2/pb1`
- `nc.coverage` 
  - Coverage when aligning the consensus sequence to the reference in the selected NextClade dataset 
  - Example: `0.9726612558735583`
- `nc.qc.{overallScore,overallStaus,missingData,mixedSites,privateMutations,snpClusters,frameShifts,stopCodons}` - 
  - QC Properties from the NextClade analysis using the selected dataset (see NextClade documentation for details)
- `num_nextclade_datasets`  
  - Number of NextClade datasets used for analysis (and correspondingly number of entries in the `nextclade_results` list)
- `nextclade_results` 
  - Full NextClade results (extracted from the JSON file produced by NextClade; see NextClade documentation for details)

Note: the values for all NextClade properties are set to `NextCladeNotRun` if NextClade was not configured to run. 

#### <sample_id>/<ref_id>/<sample_id>.<ref_id>.nextclade.tar.gz

The raw output of NextClade for the consensus sequence. Only present if NextClade was configured to be run (see above). See NextClade documentation for details of the files in this tarball. 

#### <sample_id>/<ref_id>/<sample_id>.<ref_id>.sorted.bam

Result of re-aligning the reads to the final consensus sequence (associated index also included for convenience)

#### <sample_id>/<sample_id>.kraken_report.txt

A `tsv` file sumarizing the number of reads associated to a given item in the taxonomic tree of the kraken database. For more details, check [this file format kraken2 documentation](https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown#sample-report-output-format)

  - Here is an example of what the content of this file should look like:

```tsv
81.03	3954901	3954901	U	0	unclassified
18.97	926041	0	R	1	root
18.97	926041	29	D	10239	  Viruses
18.76	915890	0	D1	2559587	    Riboviria
18.76	915890	0	K	2732396	      Orthornavirae
17.42	850471	0	P	2732408	        Pisuviricota
17.42	850445	0	C	2732506	          Pisoniviricetes
17.42	850445	0	O	76804	            Nidovirales
17.42	850445	0	O1	2499399	              Cornidovirineae
17.42	850445	0	F	11118	                Coronaviridae
[...]
```

#### consensus_sequence_properties.json

Collation of all of `<sample_id>/<ref_id>/<sample_id>.<ref_id>.properties.json` files for all consequence sequences in the entire run (for convenience)

#### mapping_summary_report.csv

> Renamed from `summary_report.csv`; see the changelog's breaking-change note. The name
> `mapping_summary_report.csv` previously belonged to the sequence-index lane's report,
> which is now `sequenceindex_summary_report.csv`.

A csv file with selected properties (per sequence) from the properties.json files above. Note that column names are different fron the JSON property names for legacy / backwards compatibility readsons. Columns:

- Sample_ID (correponds to `sample_id` in JSON)
- Virus_Taxon_ID (`virus` in JSON)
- Virus (`report_name` in JSON)
- Species (`virus_name` in JSON)
- Reference_Taxon_ID (`selected_taxid` in JSON)
- Selected_Reference (`selected_ref` in JSON)
- Flu_Segment (`flu_segment` in JSON)
- Reference_Subtype (`virus_subtype` in JSON)
- Sample_Subtype (`sample_subtype` in JSON)
- Percentage_of_Genome_Covered (`percent_positions_exceeding_depth_10` in JSON)
- Total_Mapped_Reads (`reads_mapped` in JSON)
- Total_Mapped_Bases (`bases_mapped` in JSON)
- Longest_non_N_segment (`longest_non_n_subsequence` in JSON)
- Percentage_non_N_bases (`percent_non_n_bases` in JSON)
- nc.selected_dataset (identical in JSON)
- nc.{coverage,overallScore,overallStatus,missingData,mixedSites,privateMutation,snpClusters,frameShifts,stopCodons} (identical in JSON)
- file_prefix (`id` in JSON)
- Discovered_By (`discovered_by` in JSON) — **every** method that called this species, as a
  sorted `;`-separated list: `kraken2`, `themisto2`, `metagraph_align`, `metagraph_query`.
  A row reading `kraken2;themisto2` means both found the species; Kraken2 still wins the
  reference where they agree, so the row is otherwise Kraken2's. A row with no `kraken2`
  is one only the sequence-index lane found, which requires
  `--call_consensus_for_new_species true` (default `false`).


### Secondary outputs

> If the `--developer_puplish` parameter is set to `true`, the following additional files will appear in the output folder:

```bash
<output_dir>/
├── developer_publish
│   ├── <sample_id>
│   │   └── reads_by_taxon
│   │       ├── <sample_id>_<ref_id>_R1.fq
│   │       ├── <sample_id>_<ref_id>_R2.fq
│   │       ├── [...]
│   │       ├── <sample_id>_decomposed.json
│   │       ├── <sample_id>_pre_report.tsv
│   │       └── <sample_id>_tax_to_reads.json
│   └── reference_files
│       ├── <ref_id>.fa
│       ├── [...]
```

**Kraken Output files** generated by run_kraken process

#### <sample_id>/reads_by_taxon/<sample_id>_decomposed.json 

TBD

#### <sample_id>/reads_by_taxon/<sample_id>_tax_to_reads.json

TBD

#### <sample_id>/reads_by_taxon/<sample_id>_pre_report.tsv

TBD

#### <sample_id>/reads_by_taxon/<sample_id>.<ref_id>_{R1,R2}.fq 

Pair of fastq files containing all reads which were associated to the reference with id `ref_id` n the database.

[**(&uarr;)**](#contents)

### De novo assembly + viral binning outputs

> Only produced if `--do_assembly true`.

```bash
<output_dir>/
├── assembly_sample_summary_report.csv
├── assembly_scaffold_summary_report.csv
├── assembly_vmag_summary_report.csv
├── vcontact3/
│   ├── final_assignments.csv
│   ├── final_assignments_postprocessed.csv
│   ├── final_assignments_noveltaxa.csv
│   └── postprocess_report.txt
├── <sample_id>/
│   ├── metaspades/
│   │   ├── <sample_id>_contigs.fa
│   │   └── <sample_id>_scaffolds.fa
│   ├── genomad/
│   │   ├── <sample_id>_virus_summary.tsv
│   │   ├── <sample_id>_virus.fna
│   │   └── <sample_id>_virus_proteins.faa
│   ├── vrhyme/
│   │   ├── vRhyme_best_bins.*.membership.tsv
│   │   └── vRhyme_best_bins_fasta/
│   └── checkv/
│       ├── virus_scaffolds_quality_summary.tsv
│       └── linked_bins_quality_summary.tsv
├── [...]
```

The lane's three run-level CSVs are built by `ASSEMBLY_REPORTS`
(`rvi_toolbox/bin/assembly_reports.py`) straight from the modules' own outputs,
**not** by collating the per-sample properties.json:

| file | one row per | carries |
| --- | --- | --- |
| `assembly_sample_summary_report.csv` | sample | the module counts above, plus `taxonomy_geNomad` and `taxonomy_vcontact3` — the set of taxa each classifier assigned anywhere in that sample |
| `assembly_scaffold_summary_report.csv` | geNomad viral scaffold | CheckV per-scaffold QC (`contig_length`, `gene_count`, `checkv_quality`, `completeness`, `completeness_method`) plus that scaffold's own `taxonomy_geNomad` |
| `assembly_vmag_summary_report.csv` | vRhyme bin (vMAG) | `vMAG_ID` (`<sample_id>_vRhyme_bin_<N>`), the same CheckV QC for the linked bin, plus the `vcontact3_taxonomy` assigned to that bin |

geNomad taxonomy is the lowest rank of its `;`-separated lineage. vContact3
taxonomy is reported as `<rank>:<taxon>` for the deepest rank carrying a real
name, skipping the generated `novel_*` / `unplaced_*` placeholders that mean
"not actually assigned"; it is blank where no rank qualifies.

[**(&uarr;)**](#contents)

## Configuration

### Parameters

The following command-line parameters can be used to modify the behaviour of the pipeline.  

#### Input and output 
- `manifest`: Path to the manifest file 
- `db_path`: Path to the Kraken database.
- `db_library_fa_path` (OPTIONAL): Path to the Kraken database library FASTA file.
  - By default, it assumes there is a `${params.db_path}/library/library.fna`.
- `nextclade_index_json` (OPTIONAL) : JSON file specifiying locations of datasets for nextclade analysis (see later section for how to create this file)
  - If not provided, NextClade analysis will not be performed. 
- `outdir` : folder where output files should be published. By default, it will create a subfolder called `results` in the pipeline launch directory. 

#### Kraken2Ref Handling

- `k2r_fq_load_mode`: Loading mode for Kraken2 fastq files (either `full` or `chunks`).
  - Default: `"full"`.
- `k2r_max_total_reads_per_fq`: Maximum number of reads to process per fastq file.
  - Default: `10,000,000`.

#### Kraken2ref Report Filter

- `min_reads_for_taxid`: Minimum number of reads required to assign a taxonomic ID.
  - Default: `100`.

#### Consensus building Parameters

- `do_consensus_polishing` : "polish" consensus by re-aligning the reads to the initial consensus and re-calling the consensus 
  - Default: `"true"`.
- `read_aligner` : Use bwa or minimap2 for read alignment
  - Default: `"minimap2"`
- `read_aligner_params`: Parameters to supply to the read aligner
  - Default: `"-ax sr -k11 -w 4"` (assumes minimap2)
- `mpileup_max_depth` : max depth for samtools mpileup input to ivar
  - Default: `2000`
- `ivar_initial_min_depth` : minimum depth for initial round of ivar consensus 
  - Default: `1`
- `ivar_initial_freq_threshold` : frequency threshold for initial round of ivar consensus
  - Default: `0.60`
- `ivar_polish_min_depth` : minimum depth for second (final) round of ivar consensus
  - Default: `10`
- `ivar_polish_freq_threshold` : frequency threshold for initial round of ivar consensus
  - Default: `0.75`

#### Virus Subtyping

- `scv2_keyword`: Keyword to identify SARS-CoV-2 sequences. Any taxid name equal to the string set by this parameter will be considered as SCOV2 and subjected to specific SARS-CoV-2 subtyping.
  - Default: `"Severe acute respiratory syndrome coronavirus 2"`.
- `do_scov2_subtyping`: Boolean flag to enable or disable SARS-CoV-2 subtyping via Pangolin.
  - Default: `true`.

#### De novo assembly + viral binning (`--do_assembly`)

Off by default; ported from `rvi-viral-metagenomics-pipeline`.
Requires three reference databases with no bundled default — the pipeline will
fail validation if `--do_assembly true` is set without all three:

- `genomad_db`: path to a [geNomad](https://github.com/apcamargo/genomad) database.
- `checkv_db`: path to a [CheckV](https://bitbucket.org/berkeleylab/checkv/) reference database.
- `vcontact3_db_path` (+ `vcontact3_db_version`, `vcontact3_db_domain`): path to a
  [vContact3](https://bitbucket.org/MAVERICLab/vcontact3) reference database.

Everything else has a default carried over from the source pipeline unchanged:
`metaspades_base_mem_gb`, `metaspades_subsample_limit`,
`vrhyme_min_scaffold_length`, `vrhyme_bowtie_threads`,
`vrhyme_pool_min_identity`, `vrhyme_pool_min_aligned_length`, `vrhyme_link_n`,
`vcontact3_postprocess_min_proteins`, `vcontact3_postprocess_max_proteins`,
`subsample_iterations`, `subsample_seed` — see `nextflow.config` and
`nextflow_schema.json` for their values, or `rvi-viral-metagenomics-pipeline`'s
`rvi_toolbox/subworkflows/{assemble,genomad,checkv,vrhyme,vcontact3}.json` for
the per-parameter rationale.

#### Sequence-index call gates (taxonomy and reference length)

Applied by all three sequence-index methods' species callers
(`rvi_toolbox/bin/call_themisto_species.py`, `call_metagraph_species.py`, via the shared
`bin/taxon_filter.py` and `bin/reference_lengths.py`) and folded into the same
`provisional_call` column the read-hit threshold already decides — so the lane's report
counts and the new-species handover to `MAPPING` both respect them automatically.

These gates exist because these methods query a **whole-virome** index
(`rvdb_clustered_virome`: 1,321,608 reference sequences covering every viral family RVDB
carries, and RVDB is a sequence collection rather than a genome collection). A read-hit
count therefore answers "is this sequence in the index and did reads match it", not "is
this a virus this pipeline reports, backed by enough genome to be worth a consensus".
Kraken2's lane needs neither gate: its database is already curated to the viruses of
interest, so the index itself does this job there.

**1. Taxonomy whitelist/blacklist**

- `run_taxon_filter`: master switch. Default: `true`.
- `taxon_filter_table`: `RVDB_Taxon_Current.tab`, shipped in the index's own `data/`
  directory. One row per RVDB accession carrying that accession's **already-flattened**
  NCBI lineage (realm…strain, each as a name and a taxon id), which is why no taxdump and
  no parent-walk is needed — the family taxid a whitelist entry is compared against is a
  column. Keep it in step with the index: it is regenerated per RVDB release, and a table
  from a different release resolves fewer labels, each then rejected as `unresolved`.
  Default: `/data/pam/software/themisto2/viromeindex/1.0/data/RVDB_Taxon_Current.tab`.
- `taxon_filter_whitelist`: comma-separated NCBI taxon ids. A species is kept only if one
  of them appears **anywhere in its lineage**, so a family taxid accepts every species
  beneath it, and a genus or species taxid accepts just that clade. Default is the eight
  respiratory virus families:

  | taxid | family |
  |---|---|
  | 10508 | Adenoviridae |
  | 10780 | Parvoviridae |
  | 11118 | Coronaviridae |
  | 2560066 | Sedoreoviridae |
  | 12058 | Picornaviridae |
  | 11158 | Paramyxoviridae |
  | 11244 | Pneumoviridae |
  | 11308 | Orthomyxoviridae |

- `taxon_filter_blacklist`: taxon ids rejected **even when whitelisted**, for carving an
  exception out of an accepted family (e.g. one genus within Coronaviridae). Tested before
  the whitelist, which is the only ordering that makes it useful. Default: empty.

A species the table cannot resolve at all is **rejected**, not admitted: "everything under
these families, nothing else" cannot be satisfied by a species that fails to demonstrate
it sits under one. Because that also describes a mismatched table, the callers log the
resolved fraction and name each whitelisted taxid with the family the table gives it, so a
wrong or stale table shows up in `.nextflow.log` as `unknown to the table` / `not one
label resolved` rather than as a quietly empty run.

**2. Minimum reference length**

- `min_called_reference_length`: reject a call whose resolved reference record is shorter
  than this. `0` disables the gate. Default: `1000`.

RVDB carries partial-CDS, single-gene and mRNA records alongside complete genomes, and on
a clustered index those take read-hits just as readily. A consensus against a few hundred
bases is not usable, and worse, it clears the downstream breadth gate trivially —
`new_species_min_breadth_pct` is a *percentage* of the reference length, so a 400bp
reference needs only 40 covered bases to reach 10%. Short references are the one class of
call that gets *easier* to accept the less genome there is behind it, which is why they are
cut here rather than there.

The verdict is on the record the call actually resolved to — its **most-hit** record — with
no re-pick. A species whose most-hit record is a fragment is rejected even where the index
also holds a longer record for it that took hits. That keeps the reported reference exactly
"the most-hit record", so a call's reference is never silently swapped, and the reason a
species disappeared is always visible in its own row.

When this gate is on, `INDEX_REFERENCE_LENGTHS`
(`rvi_toolbox/modules/reference_lengths.nf`) runs **once per run per reference FASTA** to
price every record — from `msweep_map_reference_fasta` for Themisto2 and
`metagraph_map_reference_fasta` for the Metagraph methods, i.e. whichever FASTA that
method's own record ids point into.

**What the outputs show**

A rejected species **keeps its row** in `<sample>/sequenceindex/*/[sample]_species_hits.tsv`
with `provisional_call` `False` and the reason recorded, so a filtered call stays readable
against its own hit count. That table gains seven columns:

| column | meaning |
|---|---|
| `taxon_filter` | `pass`, `not_whitelisted`, `blacklisted`, `unresolved`, or `off` |
| `taxonomy_id` | the accession's own NCBI taxon id, from the table |
| `family` / `family_taxon_id` | the resolved family (`;`-joined if the name is ambiguous across RVDB rows) |
| `reference_record` | the record the call resolved to (`SEQIDX_<n>` for Themisto2, an accession or taxid for Metagraph) |
| `reference_length` | that record's length in bases |
| `reference_filter` | `pass`, `short_reference`, `unknown_length`, `off`, or `not_evaluated` (the taxonomy gate already rejected it, so no reference was resolved) |

`sequenceindex_summary_report.csv` gains two per-method counts alongside
`*_n_species_considered` / `*_n_species_called`:
`*_n_species_taxon_filtered` and `*_n_species_short_reference` (`NA` when that method did
not run, or ran with that gate switched off — distinct from `0`, which means the gate ran
and rejected nothing).

#### Abundance estimation (`--do_abundance`)

Off by default, three independent sub-flags (`run_kraken2bracken`,
`run_abundance_estimation`, `run_scrub`) so any combination can run:

- `kraken2bracken_kraken2_db`: path to a Kraken2 database with a matching pre-built
  Bracken kmer-distribution file (`database<kraken2bracken_read_len>mers.kmer_distrib`)
  in the same directory. No bundled default.
- `run_scrub` requires `scrub_plate_map`: a user-supplied metadata CSV
  (`is_control`, `sample_type`, optionally `sample_well`) — mandatory if set, no default.
- `run_abundance_estimation` requires `genome_file_abundance_estimation`,
  `precomputed_index_abundance_estimation`, and `stb_file_abundance_estimation` — none
  have a working default upstream either (same "must be supplied" situation as
  `genomad_db`).

Everything else carries an upstream default unchanged, except
`cleanup_intermediate_files_abundance_estimation` (see the note above) — see
`nextflow.config`/`nextflow_schema.json`, or `rvi-viral-metagenomics-pipeline`'s
`rvi_toolbox/subworkflows/kraken2bracken.json`/`abundance_estimation.json` and
`eu1/rvi_toolbox.git`'s `subworkflows/scrub.json` for the per-parameter rationale.

### Parameter switchboard

`docs/switchboard.html` is a standalone page for answering "if I turn these flags on,
what actually runs?". Open it in a browser (it needs network access — the fonts and the
diagram renderer come from a CDN), tick the lane switches, and it redraws the process
graph for that combination, tallies the processes per lane, and writes out the matching
`nextflow run` command. It also flags a switch that is **inert** — set, but missing a
prerequisite, so Nextflow accepts the run and silently skips the step. `--run_scrub`
without `--run_kraken2bracken` and `--run_msweep` without `--run_themisto` are both
that shape.

The graph is not hand-drawn. It is assembled from one `-preview -with-dag` run per
parameter combination, so the wiring is whatever the DSL actually resolves to. Rebuild
it in three steps after changing a lane switch or the `if (params.*)` gating around a
subworkflow:

```bash
# 1. one preview DAG per case, ~23 runs. Under bsub: main.nf never runs on the head
#    node, not even -preview. No containers or real data needed -- -preview resolves
#    the graph without executing a task, on the fixtures in tests/test_data.
bsub -q normal -n 2 -M 6000 -R "select[mem>6000] rusage[mem=6000]" \
     -o preview.%J.out -e preview.%J.err bash docs/switchboard/preview_dags.sh
#    -> nf_runs/param_dags/<case>/dag.mmd (untracked scratch; pass a directory to move it)

# 2. contract those to a process-level graph, attributing each process to its switch
python3 docs/switchboard/build_graph.py        # -> docs/switchboard/graph.json

# 3. inline the graph into the page
python3 docs/switchboard/build_switchboard.py  # -> docs/switchboard.html
```

Step 2 is the part that needs care. Nextflow numbers DAG nodes by position, so ids mean
nothing across runs; processes are identified by subworkflow path plus name, and channel
and operator nodes are contracted away. Attribution comes from the cases differing by
**one** switch each: a process present with the switch on and absent with it off is gated
on that switch. A batch of all-flags-on runs would carry no attribution at all.

So adding a switch means adding an isolating case, in three places:

1. `docs/switchboard/preview_dags.sh` — a case that flips only the new switch relative to
   an existing case. If the switch is nested inside another (`run_scrub` lives inside
   `if (params.run_kraken2bracken)`), its baseline is the enclosing switch turned on,
   not the pipeline default.
2. `docs/switchboard/build_graph.py` — a `PAIRS` entry naming that baseline, the case and
   the switch, plus a `FLAGS` entry for the page's provenance table.
3. `docs/switchboard/template.html` — the switch in `SWITCHES`, with `pre:` listing the
   flags it needs to do anything, and `grp:` putting it in a lane.

`build_graph.py` prints a warning for any non-baseline process left gated on no switch,
and for any `PAIRS` entry whose cases are missing, so a forgotten step is visible rather
than silently producing a graph that ignores the new flag.

### Profiles

The pipeline is bundled with a pre-defined computation profile, `sanger_standard`. This has been crafted for use on the Wellcome Sanger Institute's HPC infrastructure. A `standard` profile for use with Singularity has also been provided, but has not been tested.  User should consider writing a profile that is compatible with their own compute infrastructure (see [profiles Nextflow documentation](https://www.nextflow.io/docs/latest/config.html#config-profiles) for more details),

[**(&uarr;)**](#contents)

---

## Unit Tests

The workflow & process unit tests for this pipeline are written in the [nf-test](https://www.nf-test.com/) (`v0.8.4`) Nextflow testing framework.

[**(&uarr;)**](#contents)

---
## Pipeline components documentation

### Processes

##### run_kraken

Executes Kraken2 on paired-end FASTQ files, producing outputs that include the classification results, classified and unclassified reads, and a summary report.

##### run_k2r_sort_reads

This process runs Kraken2Ref to parse Kraken reports and sort reads by taxonomic classification. It generates JSON files that map taxonomic IDs to read IDs and performs sorting based on the decomposed taxonomy tree if available.

##### run_k2r_dump_fastqs_and_pre_report

Extract classified reads into FASTQ files and generate a preliminary report based on the taxonomic classification data. It processes classified reads and produces a detailed report for further analysis.

##### concatenate_fqs_parts
This process concatenates FASTQ files from multiple parts into final combined FASTQ files for each taxonomic classification. This process ensures that all parts corresponding to the same taxonomic ID are merged into single files.

##### get_taxid_references

Retrieves sequences for a given taxid from a source FASTA file and indexes them for further analysis.

##### run_aligner

Mapping sequencing reads to a reference genome using BWA or minimap2, followed by post-processing steps including conversion to BAM format, sorting, and indexing. 

##### run_ivar

Generates a consensus sequence from a sorted BAM file using samtools mpileup and ivar consensus.

##### run_pangolin

The process runs the Pangolin tool on a consensus FASTA file to determine the SARS-CoV-2 lineage and extracts relevant metadata from the output.

##### run_qc_script

This process runs a QC analysis on the input BAM, FASTA, and reference files, outputting QC metrics.

##### run_nextclade

This process runs nextClade on the reconstructed sequences, recording the results in JSON (see outputs section)

[**(&uarr;)**](#contents)

### Lanes and steps

The pipeline's own Nextflow code is documented next to it, one README per directory, so a
description sits beside the file it describes rather than drifting from it here:

| directory | holds | README |
| --- | --- | --- |
| [`subworkflows/`](subworkflows/) | the **lanes** — one per `--do_*` flag, wired straight into `main.nf` | [`subworkflows/README.md`](subworkflows/README.md) |
| [`workflows/`](workflows/) | the **steps** those lanes are built from | [`workflows/README.md`](workflows/README.md) |
| [`modules/`](modules/) | individual processes — see [Processes](#processes) above | — |
| [`rvi_toolbox/subworkflows/`](rvi_toolbox/subworkflows/) | steps shared with the other RVI pipelines | [`rvi_toolbox/README.md`](rvi_toolbox/README.md) |

Start from [`subworkflows/README.md`](subworkflows/README.md): it shows how the lanes fit
together and which are on by default. Nothing in `workflows/` runs on its own — it runs
because a lane calls it — so to find out whether a given step executes, start from the lane
that calls it.

[**(&uarr;)**](#contents)

---

## Licence

[GPL-3](https://www.gnu.org/licenses/gpl-3.0.en.html)

[**(&uarr;)**](#contents)