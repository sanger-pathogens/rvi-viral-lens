# `workflows/` — the steps lanes are built from

**[⇦ main README](../README.md)** · **[⇦ subworkflows](../subworkflows/README.md)** · **[workflows](./README.md)**

One file per **step**: a reusable unit that one or more
[lanes](../subworkflows/README.md) call. Nothing here is wired into
[`main.nf`](../main.nf) directly — a step runs because a lane calls it, so to find out
*whether* something runs, start from the lane.

| | [`subworkflows/`](../subworkflows/README.md) | `workflows/` |
| --- | --- | --- |
| is | a lane, wired into `main.nf` | a step a lane calls |
| invoked by | `main.nf` only | one or more lanes |
| naming | `lower_snake_case.nf` | `UPPER_SNAKE_CASE.nf` |

Steps shared across repositories live in
[`rvi_toolbox/subworkflows/`](../rvi_toolbox/subworkflows/) instead; individual processes live
in [`modules/`](../modules/).

A step that moves out to `rvi_toolbox` stops being documented here and is documented under
the lane that composes it, in [`subworkflows/README.md`](../subworkflows/README.md) —
`VCONTACT3_RUN` and `SCRUB_DECONTAM` went that way, and are covered under
[`assembly.nf`](../subworkflows/README.md#assemblynf) and
[`abundance.nf`](../subworkflows/README.md#abundancenf).

## Index, by the lane that calls it

| step | called by | what it does |
| --- | --- | --- |
| [`SORT_READS_BY_REF`](#sort_reads_by_refnf) | [`classifying_kraken2`](../subworkflows/README.md#classifying_kraken2nf) | Kraken2 + Kraken2Ref; sorts reads by taxid and picks a reference |
| [`GENERATE_CONSENSUS`](#generate_consensusnf) | [`mapping`](../subworkflows/README.md#mappingnf) | aligns reads to a reference, calls a consensus with iVar |
| [`RUN_NEXTCLADE`](#run_nextcladenf) | [`mapping`](../subworkflows/README.md#mappingnf) | Nextclade QC/clade assignment on the consensus |
| [`SCOV2_SUBTYPING`](#scov2_subtypingnf) | [`mapping`](../subworkflows/README.md#mappingnf) | Pangolin SARS-CoV-2 lineage assignment |
| [`GENERATE_CLASSIFICATION_REPORT`](#generate_classification_reportnf) | [`mapping`](../subworkflows/README.md#mappingnf) | the pipeline's primary per-consensus report |
| [`VIRAL_THEMISTO`](#viral_themistonf) | [`classifying_index`](../subworkflows/README.md#classifying_indexnf) | Themisto2 pseudoalignment + species calling |
| [`VIRAL_METAGRAPH_ALIGN`](#viral_metagraph_alignnf) | [`classifying_index`](../subworkflows/README.md#classifying_indexnf) | `metagraph align` + species calling |
| [`VIRAL_METAGRAPH_QUERY`](#viral_metagraph_querynf) | [`classifying_index`](../subworkflows/README.md#classifying_indexnf) | `metagraph query` + species calling |
| [`GENERATE_MAPPING_REPORT`](#generate_mapping_reportnf) | [`classifying_index`](../subworkflows/README.md#classifying_indexnf) | that lane's per-sample + run-level report |
| [`GENERATE_ABUNDANCE_REPORT`](#generate_abundance_reportnf) | [`abundance`](../subworkflows/README.md#abundancenf) | that lane's per-sample + run-level report |
| [`THEMISTO_MAP_QC`](#themisto_map_qcnf--metagraph_map_qcnf) | *nothing* | **unused** — kept deliberately, see below |
| [`METAGRAPH_MAP_QC`](#themisto_map_qcnf--metagraph_map_qcnf) | *nothing* | **unused** — kept deliberately, see below |

---

## The Kraken2 → consensus path

### `SORT_READS_BY_REF.nf`

| | |
| --- | --- |
| **take** | `mnf_ch` — `(meta, [reads_1, reads_2])` |
| **emit** | `sample_taxid_ch` — `(meta, reads, ref_files)`<br>`sample_pre_report_ch`, `raw_sample_pre_report_ch` |

Runs Kraken2, then Kraken2Ref to parse the report and sort reads by taxonomic
classification, producing per-taxid FASTQs and a pre-report per sample. It also **resolves
each taxid's reference** out of the Kraken2 database's `library/library.fna`, which is why
the Kraken2 lane hands `mapping.nf` a consensus-ready tuple while the sequence-index lane
cannot.

Also exports `check_sort_reads_params()`, which `main.nf` calls during parameter validation —
the one thing in this directory reached from `main.nf`, and only as a function.

Reads are processed in `full` or `chunks` mode (`k2r_fq_load_mode`); see
[`docs/k2r_memory_escalation.md`](../docs/k2r_memory_escalation.md) for the memory-retry
behaviour.

### `GENERATE_CONSENSUS.nf`

| | |
| --- | --- |
| **take** | `sample_taxid_ch` — `(meta, reads, ref_genome)` |
| **emit** | `filtered_consensus_ch` |

Aligns paired reads to the reference (`bwa` or `minimap2`, per `read_aligner`) and calls a
consensus with `samtools mpileup` + `ivar consensus`. Two rounds when
`do_consensus_polishing` is set: an initial consensus at `ivar_initial_min_depth`, then reads
re-aligned to it and a final call at `ivar_polish_min_depth`.

`ivar consensus -n N` over `samtools mpileup -aa` makes the consensus **reference-length**,
with `N` wherever coverage was too thin. That is what makes `percent_non_n_bases` a genome
breadth figure, and it is the measurement `mapping.nf`'s
`new_species_min_breadth_pct` gate reads.

### `RUN_NEXTCLADE.nf`

| | |
| --- | --- |
| **take** | `input_ch` — `[meta, fa]` |
| **emit** | collated Nextclade JSONs |

Runs Nextclade for sequences that have a matching dataset, resolved through the
`nextclade_index_json` map of taxid → segment → dataset directory (`"ALL"` as the segment for
monopartite viruses). **If `nextclade_index_json` is not provided, this step does not run** —
it is optional, not required. Dataset layout and the index JSON's format:
[NextClade index JSON](../README.md#nextclade-index-json).

### `SCOV2_SUBTYPING.nf`

| | |
| --- | --- |
| **take** | `consensus_seq_ch` — `(meta, consensus_seq)` |
| **emit** | `scov2_subtype_out_ch` — same tuple, `meta` carrying the assigned lineage |

Runs Pangolin on consensus sequences whose taxid name matches `scv2_keyword`, adding the
SARS-CoV-2 lineage to `meta`. Disable with `--do_scov2_subtyping false`.

### `GENERATE_CLASSIFICATION_REPORT.nf`

| | |
| --- | --- |
| **take** | `report_prep_ch` — `[meta.id, meta, qc_json, nc_json]` |
| **emit** | `publish_seq_level_ch`, `publish_run_level_summaries_ch` |

The pipeline's primary report. Collects each consensus's metadata, formats it into a report
line and aggregates them into `mapping_summary_report.csv` — plus a second file listing the
sequences that were **filtered out**, so an absent row is explicable rather than just missing.
Columns: [`mapping_summary_report.csv`](../README.md#mapping_summary_reportcsv).

## The sequence-index methods

All three take the same input and emit the same two things, so
[`classifying_index`](../subworkflows/README.md#classifying_indexnf) can treat them
interchangeably:

| | |
| --- | --- |
| **take** | `reads_ch` — `(meta, read_1, read_2)`, preprocessed |
| **emit** | `species_hits` — `<sample>_species_hits.tsv`<br>`index_label_map` — record id → species, **optional** per sample (unwritten when nothing cleared the gates) |

Each caps its input depth before querying (`msweep_subsample_limit` /
`metagraph_align_subsample_limit`, **per mate**). Species calling is a statistical estimate,
not an assembly: more depth costs time and memory without changing which species get called.

None of them map any reads — see
[`classifying_index`](../subworkflows/README.md#classifying_indexnf) for the three gates that
decide a call, and [`mapping`](../subworkflows/README.md#mappingnf) for where breadth is
measured instead.

### `VIRAL_THEMISTO.nf`

The lane's **default** method (`--run_themisto`, on unless turned off). Pseudoaligns against a
pre-built Themisto2 `.thm2` index, then calls species directly from pseudoalignment read-hit
counts — no probabilistic model. Additionally emits `pseudoalignments` and `ref_groups` for
the abundance lane's optional mSWEEP.

Species are resolved positionally: line *N* of `species_labels.txt` is reference sequence *N*
in the index, which is also record *N* of `msweep_map_reference_fasta`. That 1:1 alignment is
what lets a call name its reference as a bare `SEQIDX_<n>` token and lets `mapping.nf` extract
it by exact id with nothing in between. Core processes are adapted from the gemsweep pipeline
(themisto2 branch).

> `SUBSAMPLE_ITER` is **lossy and asymmetrically so**: a sample below the limit passes through
> with `meta` untouched, but one that is actually subsampled has `meta` rebuilt from scratch
> carrying only a renamed `id`. This workflow restores the original `meta` afterwards, and
> must — without it, report counts key off the post-subsample id while the report backbone
> keys off the pre-subsample one, and the run aborts on a null `meta`. Small test samples sit
> below the limit and never trip it.

### `VIRAL_METAGRAPH_ALIGN.nf`

Aligns each mate independently against a pre-built metagraph de Bruijn graph + annotation
index, counts reads per species from the alignment labels, and calls species above
`metagraph_align_min_hits`. `--run_metagraph_align`, off by default.

Subsampling matters more here than elsewhere: `metagraph align` emits one label line per
matched k-mer window, so its output scales with input depth far faster than a normal aligner's
would.

### `VIRAL_METAGRAPH_QUERY.nf`

The same species calling against the **same** graph, annotation and thresholds, but reaching
them through `metagraph query --query-mode labels` instead of `metagraph align`.
`--run_metagraph_query`, off by default.

The two Metagraph steps are **alternative methods against the same reference data**, not
independently configured pipelines — which is why this one reuses every `metagraph_align_*`
param. They publish to different subdirectories (`metagraph_hits` vs `metagraph_query_hits`)
so running both for one sample cannot overwrite either result.

## Lane reports

### `GENERATE_MAPPING_REPORT.nf`
### `GENERATE_ABUNDANCE_REPORT.nf`

| | |
| --- | --- |
| **take** | `report_prep_ch` — `(sample_id, meta)`, one per sample |
| **emit** | `publish_seq_level_ch`, `publish_run_level_summaries_ch` |

The same shape, one per lane. `meta` arrives already carrying that lane's counts, so each step
just JSON-dumps it per sample and then concatenates into one run-level CSV with a row per
sample. Column set is the **union of keys seen across all records** rather than a fixed
mapping — a sample missing a key gets a blank cell — so a lane adding a count needs no change
here.

Despite the name, `GENERATE_MAPPING_REPORT` belongs to the **sequence-index** lane and writes
`sequenceindex_summary_report.csv`. The `mapping_summary_report.csv` name belongs to
[`GENERATE_CLASSIFICATION_REPORT`](#generate_classification_reportnf); see the changelog's
breaking-change note.

Neither writes a run-level JSON any more: it duplicated the CSV exactly, so nothing published
it.

## Kept but unused

### `THEMISTO_MAP_QC.nf` / `METAGRAPH_MAP_QC.nf`

**Nothing invokes either.** For every called species they mapped the sample's own reads
against that species' most-hit reference record and recorded breadth of coverage, mean depth,
mapping/base quality and reads mapped — one mapping per species, as a sanity check on the
read-hit calls.

They were removed from the lane because they meant every real call was mapped **twice**, by
two different aligners: once here to measure breadth, then again for consensus. Breadth is now
measured once, downstream, from the consensus alignment
[`mapping.nf`](../subworkflows/README.md#mappingnf) performs anyway.

They are **kept rather than deleted** because they are the only thing that can measure breadth
for a species *before* deciding to spend a consensus on it. Restore them if that ordering ever
matters — if index noise volume makes the discarded consensuses expensive, or if breadth is
wanted for calls that never get a consensus at all, such as ones Kraken2 already found. Their
params are still defined in `nextflow.config`, marked UNUSED.
