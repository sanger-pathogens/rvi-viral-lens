# `subworkflows/` — the pipeline's lanes

**[⇦ main README](../README.md)** · **[subworkflows](./README.md)** · **[workflows ⇨](../workflows/README.md)**

One file per **lane**: a self-contained branch of the pipeline that [`main.nf`](../main.nf)
switches on with a `--do_*` flag and feeds the same preprocessed reads. A lane owns its own
report counts and its own `publish*` calls, so adding or removing one touches `main.nf` and
nothing else.

The reusable steps a lane composes live one directory over, in
[`workflows/`](../workflows/README.md). The distinction is worth keeping straight:

| | `subworkflows/` | [`workflows/`](../workflows/README.md) |
| --- | --- | --- |
| is | a lane, wired directly into `main.nf` | a step a lane calls |
| invoked by | `main.nf` only | one or more lanes |
| how many run | whichever flags are on | as many times as the lanes need |
| naming | `lower_snake_case.nf` | `UPPER_SNAKE_CASE.nf` |

## The lanes

| lane | flag | default | run-level report |
| --- | --- | --- | --- |
| [`classifying_kraken2.nf`](#classifying_kraken2nf) | `--do_mapping` | **on** | (feeds `mapping.nf`) |
| [`mapping.nf`](#mappingnf) | `--do_mapping` | **on** | `mapping_summary_report.csv` |
| [`classifying_index.nf`](#classifying_indexnf) | `--do_sequence_index` | off | `sequenceindex_summary_report.csv` |
| [`assembly.nf`](#assemblynf) | `--do_assembly` | off | three `assembly_*_summary_report.csv` |
| [`abundance.nf`](#abundancenf) | `--do_abundance` | off | `abundance_summary_report.csv` |

`--do_mapping` covers **two** files, because Kraken2 classification and the consensus pass
are halves of one pipeline — the original viral-lens, before the other lanes existed. It is
the only lane on by default, for that reason: defaulting it off would silently change every
existing command. At least one lane must be enabled, or the run is rejected at startup.

### How they fit together

Two lanes **classify** and one **builds consensus for both of them**:

```
preprocessed reads ─┬─→ classifying_kraken2.nf ─┐   (--do_mapping)
                    │                            ├─→ mapping.nf ─→ consensus, Nextclade,
                    ├─→ classifying_index.nf ───┘   (--do_mapping)  subtyping, classification
                    │        (--do_sequence_index)                  report
                    ├─→ assembly.nf      (--do_assembly)
                    └─→ abundance.nf     (--do_abundance)
                              ↑
                    Themisto2 pseudoalignments, if the
                    sequence-index lane produced them
```

`mapping.nf` is deliberately **not** per-lane. Two classifiers can each find a species worth
a consensus, and both hand over to the one shared consensus/Nextclade/subtyping/report pass
rather than growing parallel copies of it. That union is why a species only Themisto2 or
Metagraph found still gets Nextclade, SARS-CoV-2 subtyping and a row in the classification
report — it once had its consensus published on its own with none of that.

Cross-lane dependencies, all rejected at startup rather than mid-run:

| this | needs | because |
| --- | --- | --- |
| `--run_msweep` | `--do_sequence_index --run_themisto` | mSWEEP estimates from Themisto2's pseudoalignments, not from reads |
| `--call_consensus_for_new_species` | `--do_mapping` | `mapping.nf` is what resolves those references, maps them and applies the breadth gate |

**Running with `--do_mapping false`** makes the run reference-free — there is nowhere for a
consensus to come from. Two knock-on effects worth knowing:

- Kraken2 database parameters stop being required, so an assembly-only or abundance-only run
  needs no Kraken2 database. The **manifest** is still validated; every lane reads it.
- `classifying_index.nf` still runs and still reports its calls, but reports
  `overlapping_n_species` as `NA` — there are no Kraken2 calls to overlap with, which is not
  the same fact as an overlap of zero.

### Report convention shared by every lane

Run-level CSVs write **`NA`, not `0`, for a step that did not run** in that execution; `0`
means "ran and found nothing". The two are different facts and a report that spells both `0`
cannot be read correctly — `new_species_candidates_n` read `0` whenever
`--call_consensus_for_new_species` was off, indistinguishable from the lane genuinely finding
nothing new. The sequence-index gate counts take this one level further: `NA` when the method
did not run *or* ran with that gate switched off.

---

## `classifying_kraken2.nf`

Classifies reads by Kraken2 taxid and selects a reference per taxid — everything up to, but
not including, consensus. Runs unless `--do_mapping false`.

| | |
| --- | --- |
| **take** | `preprocessed_3tuple_ch` — `(meta, read1, read2)` |
| **emit** | `sample_taxid_ch` — `(meta, [read_1, read_2], reference_fasta)`<br>`sample_report_with_join_key_ch` — `[join_key, report_meta]`<br>`identified_species_ch` — `[sample_id, [normalized_species_name, ...]]` |
| **composes** | [`SORT_READS_BY_REF`](../workflows/README.md#sort_reads_by_refnf) |

This is the Kraken2/Kraken2Ref half of what `mapping.nf` used to do end to end. The interface
it emits is **deliberately identical in shape** to what the sequence-index lane provides, so
`mapping.nf` can be driven by either.

`identified_species_ch` is the set of species Kraken2 actually *acted on* (the pre-report's
`virus_name` + `ref_selected`, normalised), not every row of the raw Kraken2 report. Two
places consume it: `classifying_index.nf` counts how many of its own calls overlap, and
`mapping.nf` uses it to decide which sequence-index calls are genuinely new.

## `mapping.nf`

Consensus generation, lineage calling and classification reporting, over the union of **both**
classifiers' findings. Runs unless `--do_mapping false`. The only place in the pipeline where
either classifier's reads get mapped.

| | |
| --- | --- |
| **take** | `kraken2_sample_taxid_ch`, `kraken2_report_ch`, `identified_species_ch` — from `classifying_kraken2.nf`<br>`index_species_calls_ch`, `index_called_species_ch` — from `classifying_index.nf`, or `Channel.empty()`<br>`reads_ch` — `(meta, read1, read2)` |
| **composes** | [`GENERATE_CONSENSUS`](../workflows/README.md#generate_consensusnf) · [`RUN_NEXTCLADE`](../workflows/README.md#run_nextcladenf) · [`SCOV2_SUBTYPING`](../workflows/README.md#scov2_subtypingnf) · [`GENERATE_CLASSIFICATION_REPORT`](../workflows/README.md#generate_classification_reportnf) · `INDEX_REFERENCE_FASTA` / `EXTRACT_REFERENCE_RECORD` / `EXTRACT_METAGRAPH_REFERENCE_RECORD` ([`modules/reference_subset.nf`](../modules/reference_subset.nf)) |

**The two classifiers hand over different shapes, on purpose.** Kraken2 arrives
*consensus-ready* (reads + reference), because `SORT_READS_BY_REF` resolves its references as
part of classifying. The sequence-index lane arrives as *species calls plus a reference
record id*, and no reads. Extracting that reference and pairing reads for consensus is this
lane's job — done only for the calls that survive the "Kraken2 already found it" filter,
which is the whole point of the asymmetry.

**Where Kraken2 wins.** Where both classifiers found the same species, Kraken2's call and its
reference win; only the species Kraken2 missed are resolved and mapped off the index's calls.

**Why the breadth gate lives here.** `classifying_index.nf` calls species without mapping
anything, so its calls arrive with no breadth figure attached. `params.new_species_min_breadth_pct`
is therefore applied *after* the consensus alignment rather than up in the classifier. The cost
of that ordering is that a noise call's consensus is computed and then discarded; the saving is
that a real call is mapped once instead of twice. The measurement is `percent_non_n_bases` from
the consensus QC JSON — breadth at iVar's minimum depth over the full reference length.
Kraken2-found species are not subject to it.

Two reference-extraction processes, not one, because the two index families report record ids
in different namespaces: positional `SEQIDX_<n>` into `msweep_map_reference_fasta` for
Themisto2, a bare taxid or accession into `metagraph_map_reference_fasta` for Metagraph. The
`reference_source` field on each call is what keeps them straight.

## `classifying_index.nf`

The sequence-index counterpart to `classifying_kraken2.nf`: classifies reads against
**pre-built sequence indexes** and produces species **calls**, not consensus sequences.
Opt-in with `--do_sequence_index`.

| | |
| --- | --- |
| **take** | `preprocessed_3tuple_ch` — `(meta, read1, read2)`<br>`identified_species_ch` — from `classifying_kraken2.nf` |
| **emit** | `species_calls_ch` — `[sample_id, {species_name, reference_record, reference_source, hit_count, method}]`; empty unless `--call_consensus_for_new_species`<br>`called_species_ch` — `[sample_id, species_lower, method]`, ungated<br>`themisto_pseudoalignments`, `themisto_ref_groups` — handover for the abundance lane's mSWEEP |
| **composes** | [`VIRAL_THEMISTO`](../workflows/README.md#viral_themistonf) · [`VIRAL_METAGRAPH_ALIGN`](../workflows/README.md#viral_metagraph_alignnf) · [`VIRAL_METAGRAPH_QUERY`](../workflows/README.md#viral_metagraph_querynf) · [`GENERATE_MAPPING_REPORT`](../workflows/README.md#generate_mapping_reportnf) |

Up to three methods run **in parallel off the same reads** — not downstream of one another —
each independently flagged, all feeding **one** `GENERATE_MAPPING_REPORT` call rather than
building a second report path. `sequence_index_sample_ch` is the join backbone (every sample
that reaches the lane), so a sample gets a report row even if only one method ran for it, and
with several flags on, one row carries every enabled method's counts.

It maps **no reads**. See [`mapping.nf`](#mappingnf) for where breadth is measured instead.

**Three gates decide a call**, all inside the species callers and all folded into the single
`provisional_call` column — so this lane and everything downstream respect them without
knowing they exist:

| gate | threshold |
| --- | --- |
| read hits | `themisto_align_min_hits` / `metagraph_align_min_hits` |
| taxonomy | lineage must sit under `taxon_filter_whitelist` and not under `taxon_filter_blacklist`, resolved via `taxon_filter_table` |
| reference | the resolved record must be at least `min_called_reference_length` bases |

The last two exist because these methods query a **whole-virome** index, unlike Kraken2's
curated database. Read hits there answer "is this sequence in the index and did reads match
it", not "is this a virus we report, with enough genome behind it to be worth a consensus".
Without them, most of what cleared min-hits on a respiratory sample was neither — phage, plant
and insect viruses sharing k-mers, and partial-CDS records that pass a
percentage-of-reference breadth gate trivially *because* they are short. Full configuration:
[Sequence-index call gates](../README.md#sequence-index-call-gates-taxonomy-and-reference-length).

`called_species_ch` is deliberately **not** derived from `species_calls_ch`: the latter sits
behind `--call_consensus_for_new_species` and is empty by default, which would make both
`overlapping_n_species` and the report's `Discovered_By` silently wrong rather than absent.

> With more than one method enabled, which method's row wins for a species both called — and
> therefore which record and hit count get reported for it — is task completion order, so not
> reproducible run to run. Harmless on the defaults (only `run_themisto` is on). The species
> itself is unaffected; only its reported record and count.

## `assembly.nf`

Reference-free lane: de novo assembly, viral identification, binning, QC and taxonomic
assignment. Opt-in with `--do_assembly`. Runs alongside the classifier lanes, not instead of
them.

| | |
| --- | --- |
| **take** | `preprocessed_3tuple_ch` — `(meta, read1, read2)` |
| **composes** | `ASSEMBLE_META` · `GENOMAD_CLASSIFY` · `VRHYME_BIN` · `CHECKV_QC` (all [`rvi_toolbox/subworkflows/`](../rvi_toolbox/subworkflows/)) · [`VCONTACT3_RUN`](../workflows/README.md#vcontact3_runnf) · `ASSEMBLY_REPORTS` |

Almost every step is per-sample, with two exceptions that need the whole batch: `VRHYME_BIN`
pools every qualifying sample's scaffolds into one bowtie2 index so coverage covariance across
the batch informs each sample's binning, and `VCONTACT3_RUN` runs vContact3 once across all
samples' combined input. Requires three reference databases with no bundled default —
`genomad_db`, `checkv_db`, `vcontact3_db_path`.

Reports: three run-level CSVs at different grains — one row per sample, per geNomad viral
scaffold, and per vRhyme bin. The lane writes **no** per-sample `properties.json`;
`ASSEMBLY_REPORTS` builds the CSVs from the modules' own output files.

## `abundance.nf`

Species-abundance estimation by four independent, separately flagged methods. Opt-in with
`--do_abundance`.

| | |
| --- | --- |
| **take** | `preprocessed_3tuple_ch` — `(meta, read1, read2)`<br>`themisto_pseudoaln_ch`, `themisto_ref_groups_ch` — from `classifying_index.nf`, empty unless `--run_themisto` |
| **composes** | `KRAKEN2BRACKEN` · `ABUNDANCE_ESTIMATION` · `MSWEEP` ([`rvi_toolbox`](../rvi_toolbox/subworkflows/)) · [`SCRUB_DECONTAM`](../workflows/README.md#scrub_decontamnf) · [`GENERATE_ABUNDANCE_REPORT`](../workflows/README.md#generate_abundance_reportnf) |

| method | flag | reads from |
| --- | --- | --- |
| Kraken2 + Bracken | `--run_kraken2bracken` | preprocessed reads |
| sourmash / inStrain | `--run_abundance_estimation` | preprocessed reads |
| SCRuB decontamination | `--run_scrub` | Kraken2+Bracken's **whole-run** output |
| mSWEEP | `--run_msweep` | Themisto2's **pseudoalignments** |

Kraken2+Bracken and `ABUNDANCE_ESTIMATION` run in parallel off the same reads, and not
downstream of any other lane. Two are shaped differently:

- **SCRuB** is a final whole-run step on Kraken2+Bracken's output only — decontamination
  inherently needs the whole batch (samples + controls) together, so it cannot be per-sample,
  and `ABUNDANCE_ESTIMATION` does not go through it.
- **mSWEEP** is the one step here that does not start from reads. It estimates abundance from
  Themisto2's pseudoalignments, so it also requires `--do_sequence_index --run_themisto`.
  mSWEEP lives here rather than in the sequence-index lane because it *estimates abundance*
  rather than calling species — it calls no species at all.

Pseudoalignment cleanup is owned by whichever lane reads them last: with `--run_msweep` set,
this lane deletes them after mSWEEP, because deleting them in the sequence-index lane would
race it.
