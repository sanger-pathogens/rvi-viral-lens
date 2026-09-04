# rvi_integration_1 — handoff to continue on HPC

You're picking up mid-flight on branch `rvi_integration_1` of `viral-lens`. This branch is
integrating functionality from a sibling pipeline, `rvi-viral-metagenomics-pipeline`, into
`viral-lens/main.nf`. This file has been updated five times now — first written at
commit `d759968`, updated through `bf71c67` after a farm run (item 1, and the Themisto2
method of item 3), updated through `d0760aa` after a farm-less round ported Metagraph
sequence-to-graph alignment (item 3's second method), updated through `0241641` after a
third, also farm-less round added Metagraph pseudoalignment (item 3's third and last
method), updated through `1536267` after a fourth, still farm-less round built the
entire abundance estimation lane (item 4), updated again through `71c15cb` (`HEAD` at
the time of this edit) after a fifth round wired up wider input handling (item 5).
**Every round since the first has had no `singularity`/reference-data/network access at
all** — items 3, 4, and 5 are all fully wired and DSL-checked, but only item 1 and
Themisto2/mSWEEP have actually executed. `git log` is the source of truth if it's since
moved further.

Read this whole file before touching anything. It front-loads facts (repo layout, remotes,
commit hashes, file conventions) discovered the hard way in the prior session, specifically
so you don't have to re-derive them.

## Why this exists

The previous session ran in a sandbox with no `singularity`/`docker` and none of the
reference databases this pipeline needs (geNomad/CheckV/vContact3 DBs). Everything was
built and reasoned through against real source code, but the actual containerized
processes have never executed. **Your first job on HPC, where those things exist, is to
actually run this for real and fix whatever real execution surfaces** — treat everything
below as "believed correct, not yet proven."

## Repo layout

Two sibling repos under the same parent directory (the parent itself is not a git repo):

- `viral-lens/` — the pipeline being extended. You're working here. Branch
  `rvi_integration_1`, forked from `main`.
- `rvi-viral-metagenomics-pipeline/` — the source of the functionality being ported. Treat
  as read-only reference; don't need to change it. It's not a copy for this task, it's the
  actual sibling repo, so `git log`/`git show`/`git diff` against it directly.

  **On the farm it is NOT next to `viral-lens/`.** There are ~15 checkouts, on different
  branches, under
  `/lustre/scratch126/pam/projects/rvidata/personal/eu1/rvi-viral-metagenomics/*/rvi-viral-metagenomics-pipeline`.
  They differ in what they contain, so pick deliberately:
  `metagraph-workflow/` (branch `main`) is the richest — it has `msweep.nf`,
  `metagraph_align.nf` and `metagraph_map_qc.nf` together, and is what the Themisto2 port
  was taken from. `broad-pipeline/` (`feature/assisted-assembly-genomad-bins`) has
  *dropped* `VIRAL_MSWEEP`. `msweep-map-sourmash-ref/` shows `VIRAL_MSWEEP` actually
  invoked from `main.nf`. Checking these branches is genuinely useful for seeing how
  something was implemented upstream — the user asked for it explicitly.

Both have their own `rvi_toolbox/` git submodule, but pointed at **two different, diverged
forks** (this matters a lot, see below):

- `viral-lens/rvi_toolbox` → remote `rvi/rvi_toolbox.git` (shared org fork)
- `rvi-viral-metagenomics-pipeline/rvi_toolbox` → remote `eu1/rvi_toolbox.git` (personal fork)

## What's already done (commits `af41a75`..`d759968` on `rvi_integration_1`)

In order:

1. `155dbc5` — renamed old `main.nf` → `mapping_pipeline_main.nf` (kept intact, runnable),
   rebuilt `main.nf` as the entry point being extended.
2. `af41a75` — added `docs/nf-metro/route_map.{mmd,svg}` + `docs/nf-metro/README.md`: a
   hand-authored [nf-metro](https://github.com/seqeralabs/nf-metro) diagram of the target
   end-state pipeline, and the method for regenerating/extending it. **Update this diagram
   whenever you change the pipeline's shape** — it's meant to stay current, not be a
   one-off sketch. `docs/nf-metro/README.md` has the full authoring workflow
   (install/validate/render, including the layout-engine gotchas we hit).
3. `c545c24` — ported the de novo assembly + viral binning lane (`ASSEMBLE_META` →
   `GENOMAD_CLASSIFY` → `VRHYME_BIN` + `CHECKV_QC` → `VCONTACT3_RUN`) from
   `rvi-viral-metagenomics-pipeline`, byte-for-byte except two hardcoded
   `${projectDir}/rvi_toolbox/bin/...` script paths fixed to `${projectDir}/bin/...`.
   Landed as **viral-lens-owned files** (`workflows/`, `modules/`, `bin/`), not inside the
   `rvi_toolbox` submodule — see "The rvi_toolbox fork problem" below for why, and don't
   "fix" this without reading that section first.
4. `d727cbf` — added `GENERATE_ASSEMBLY_REPORT.nf` / `GENERATE_MAPPING_REPORT.nf` /
   `GENERATE_ABUNDANCE_REPORT.nf`, all mirroring the existing
   `GENERATE_CLASSIFICATION_REPORT.nf` shape (per-sample `meta` JSON dump, then run-level
   concatenation), sharing new `modules/write_lane_report.nf` +
   `bin/write_lane_summary.py`. Each has an nf-test workflow test
   (`tests/workflows/GENERATE_*_REPORT.nf.test`).
5. `375485a` — wired the assembly lane into `main.nf` behind `params.do_assembly` (default
   `false`), including sample-level `meta.plus()` counting glue (see "Meta propagation" below).
6. `de76027` — fixed two real bugs a `nextflow run main.nf -preview --do_assembly true ...`
   dry-run surfaced (stale include path, duplicate process name — see commit message).
   **This is the level of verification so far: DSL wiring resolves, nothing has actually
   executed.**
7. `d759968` — documented all of the above in `README.md`, including a
   "rvi_integration_1: work in progress" section that's the shorter, in-repo version of
   this file — keep it in sync as you go.

`GENERATE_MAPPING_REPORT.nf` / `GENERATE_ABUNDANCE_REPORT.nf` exist and pass their own
standalone tests but are **not called from `main.nf` yet** — their upstream lanes don't
exist yet. That's most of what's left.

## Item 1 is DONE (see commits after `bf71c67`)

The assembly lane has been run for real on the farm and completes end to end
(`Success: true`). Everything in this section is now established fact, not belief.

### Environment that actually works

```bash
module load ISG/experimental/irods/4.2.7
module load ISG/singularity/3.11.4
module load nextflow/24.10.6      # NOT 23.10.1: nextflow.config requires >=24.10.3
module load bsub.py
module load cellgen/nf-test/0.9.5 # nf-test is not in the pam module namespace

# The shared library is read-only and lacks several images; keep a personal cache.
export NXF_SINGULARITY_LIBRARYDIR=/data/pam/installs/custom_installs/nextflow_singularity_library/
export NXF_SINGULARITY_CACHEDIR=/lustre/scratch126/pam/projects/rvidata/personal/eu1/pipeline-integration/.singularity_cache
# singularity's own blob cache defaults to $HOME/.singularity, whose quota is tiny.
export SINGULARITY_CACHEDIR=/lustre/scratch126/pam/projects/rvidata/personal/eu1/pipeline-integration/.singularity_blobcache
export SINGULARITY_TMPDIR=$SINGULARITY_CACHEDIR/tmp
export LSB_DEFAULTGROUP=rvidata
```

Reference data, all verified present:

| param | path |
|---|---|
| `--db_path` | `/lustre/scratch126/pam/projects/rvidata/pipeline_resources/kraken_databases/production/viral_lens_kdb_v1.5.2` |
| `--genomad_db` | `/data/pam/software/genomad/genomad_db/genomad_db` |
| `--checkv_db` | `/data/pam/software/ViWrap/CheckV_db/` |
| `--vcontact3_db_path` | `/data/pam/software/vcontact3/` (holds `v232`, `v236`; default version is now 236) |

Working run directories (launch scripts, logs, outputs, LSF driver output) live at
`/lustre/scratch126/pam/projects/rvidata/personal/eu1/pipeline-integration/nf_runs/`:
`run1/` (assembly lane) and `run2_seqindex/` (Themisto2 lane). Each has a `launch.sh`
that is the whole reproducible invocation — read it first, then
`bsub -q normal -n 2 -M 8000 -R "select[mem>8000] rusage[mem=8000]" -o driver.out bash launch.sh`.
Run the nextflow driver itself under `bsub` too; the user asked for real runs on LSF, not
local. Cheap `-preview` DAG checks locally are fine.

Traps that cost real time, so you don't repeat them:

- **Run from lustre, not the agent scratchpad.** `/tmp/claude-*` is local to the submit
  host; an LSF job cannot see it and dies with no output at all.
- **Primary unix group is `team230f`, quota ~1 MB.** Setgid dirs make normal writes
  inherit `rvidata`, but singularity's unprivileged image build runs in a user namespace
  that uses the primary gid and fails with "disk quota exceeded". Pull under
  `sg rvidata -c '...'`, or better, reuse an existing image (below).
- **`$HOME` was full** (51146M of 51200M, 12G of it `~/.singularity`). Worth clearing.
- **A partial pull leaves a corrupt SIF** that fails later with "bad superblock for
  squashfs image". Compare byte sizes against a known-good copy.
- Prebuilt images live in
  `/data/pam/installs/custom_installs/nextflow_singularity_library/` (assembly-lane
  tools) and
  `/lustre/scratch126/pam/projects/rvidata/personal/eu1/rvi-viral-lens/provisional_v1.5.3_49643_1/work/singularity`
  (taxid-lane `rvi-vp-*` images). Symlinking these into `NXF_SINGULARITY_CACHEDIR`
  avoids every pull.

### What the real run proved and fixed

The `count_*` helpers parse real files correctly -- column names confirmed against real
output (`seq_name`/`n_genes`, `scaffold`/`bin`, `contig_id`/`checkv_quality`,
`Genome`/`genus_prediction`). One was still wrong for a different reason:
`count_vrhyme_membership()` reported the bin count for both of its fields, because
Groovy's `List.unique()` de-duplicates **in place**. Fixed.

Also fixed: strings-not-Paths into `ASSEMBLE_META`, the whole lane running on the submit
host under `sanger_standard`, an undeclared param, a nonexistent default vContact3 DB
version, and the vContact3 genome report resolving outside the versioned DB directory.
See the git log for each.

Still open from this section:

- `count_vcontact3_for_sample()` still does naive `.split(',')` CSV parsing. The real
  `final_assignments.csv` did not trip it, but nothing guarantees that for other inputs.
- `checkv_n_high_quality`/`checkv_n_medium_quality` count only
  `virus_scaffolds_quality_summary.tsv`. The first real sample had 8 Low-quality
  scaffolds there but 1 High-quality *bin* in `linked_bins_quality_summary.tsv`, which
  no report field currently reflects. Decide whether bin quality belongs in the report.
- `sanger_standard` caps `max_time` at 6h, clamping the lane's `time_12` labels.
- Only single-sample has been run. **`VRHYME_BIN` pools scaffolds across samples**
  (`POOL_VIRAL_SCAFFOLDS` -> one bowtie2 index -> `COVERM_DEPTH` -> per-sample subset),
  so a multi-sample run is still needed to exercise that machinery.

### Output layout (agreed with the user, implemented)

Per-sample outputs are grouped by lane; see README's "Output layout". `sequenceindex/`
and `abundance/` are reserved for items 3 and 4 below. All publishing now goes through
`params.outdir`; the ported `params.results_dir` is gone from viral-lens-owned files.

## FARM-RUN ROUND: items 3, 4 and 5 have now been executed for real

The three lanes below were merged "pending a farm run". They have now been run on
LSF with real reference data. Read this section before the three that follow it —
those were written before any of it had executed, and this supersedes them.

### Status after this round

| lane | flags | farm status |
|---|---|---|
| assembly (single + **multi-sample**) | `--do_assembly` | **PASS** |
| Themisto2/mSWEEP | `--do_sequence_index --run_msweep` | **PASS** |
| pseudoalign via Metagraph | `--do_sequence_index --run_metagraph_query` | **PASS** |
| Kraken2+Bracken | `--do_abundance --run_kraken2bracken` | **PASS** |
| SCRuB | `+ --run_scrub --scrub_plate_map` | **PASS** |
| ABUNDANCE_ESTIMATION | `--do_abundance --run_abundance_estimation` | **PASS** |
| Metagraph **align** | `--do_sequence_index --run_metagraph_align` | **PASS** |
| MIXED_INPUT (local reads manifest) | `--do_mixed_input --manifest_of_reads` | **PASS** |
| MIXED_INPUT (ENA / iRODS sources) | `--manifest_ena` / `--studyid` etc. | **NOT run** — needs ENA/iRODS access |

**Every lane now has a passing farm run.** The only untested paths left are
MIXED_INPUT's ENA and iRODS *sources*; its local-reads-manifest source is proven.

### What the farm runs fixed

Every one of these was invisible to `-preview`:

- **All three lanes shipped with `null`/`""` reference-data params**, so none could
  start: `Channel.fromPath(null)` fails with an opaque "Missing `fromPath`
  parameter". Real values are now defaulted (see `nextflow.config`), taken from the
  source pipeline's configs and confirmed to exist. Note `kraken2bracken_kraken2_db`:
  rvi_toolbox names `viral_kraken/`, which does not exist here.
- **`instrain_profile_options = "--database-mode"`** — inStrain 1.9.0 spells it
  `--database_mode`; the hyphenated form exits 2. Straight from rvi_toolbox's config,
  so this lane cannot ever have completed upstream against this container.
- **`genomad_db`/`checkv_db`/`vcontact3_db_path` were still null**, so the assembly
  lane failed deep into a run with a literal `null` as geNomad's database argument.
- **Repeated `withName:` selectors clobber each other.** `sanger_standard` selects
  KRAKEN2/BRACKEN/INSTRAIN_PROFILE by bare name for `executor = 'lsf'`, and that
  block *replaced* the top-level one wholesale rather than merging. This silently
  discarded the `shell = ['/bin/bash','-u']` settings those processes need, and
  swallowed a first attempt at the output-layout fix. Both now use qualified regex
  selectors (`'.*:KRAKEN2'`). **Check `nextflow config -profile sanger_standard`
  after touching any `withName:` block** — the resolved output is the only reliable
  way to see what survived.
- **Abundance outputs were not lane-grouped**: kraken2/bracken/instrain published to
  `<sample>/<tool>/`. Re-homed under `<sample>/abundance/<tool>/` via config
  overrides rather than editing the shared submodule modules.
- **SCRuB row order** (multi-sample only): the abundance matrix was natural-sorted
  while the metadata kept plate-map order; SCRuB compares them positionally and
  refused to run.
- **vRhyme aborted the whole run** on a sparse sample (a negative control) that its
  own internal screen rejected. The `VRHYME_BIN.nf` gate counts scaffolds in the
  fasta, which is necessarily an over-estimate of what vRhyme keeps, so the gate
  cannot prevent this; that one message is now treated as "no bins".

### Memory, measured rather than guessed

`bjobs -l <id>` gives `TERM_MEMLIMIT` and `MAXMEM` — use it instead of guessing.

| process | peak | setting |
|---|---|---|
| THEMISTO_PSEUDOALIGN | 6.6GB | `25.GB * task.attempt` — ample |
| METAGRAPH_QUERY | 15.6GB | `25.GB * task.attempt` — ample |
| METAGRAPH_ALIGN | **52.7GB** | `32.GB * task.attempt * 2` (64GB first) |

METAGRAPH_ALIGN is driven by the 15GB coordinate annotation; it was OOM-killed at
both 25GB and 50GB. Note its `time`/`memory` directives are raw (not `time_*`/`mem_*`
labels), so they **bypass `check_max`** and `sanger_standard`'s `max_time = 6.h` /
`max_memory = 128.GB` do not clamp them — that is why it asks for 12:00.

### Reproducing

`nf_runs/` holds one directory per configuration, each with the `driver.out` its run
produced: `regress`, `mg_align`, `mg_query`, `abund_k2b`, `abund_est`, `scrub_multi`,
`assembly_multi`, `mixed_mgquery` (plus `run1`/`run2_seqindex` from the first round).

Three launchers, all reading `RUN_NAME` and `EXTRA_ARGS` from the environment, all
sourcing `env.sh`:

| launcher | input | manifest |
|---|---|---|
| `launch_generic.sh` | 1 sample, 10k reads | `one_sample_manifest.csv` |
| `launch_multi.sh` | 3 samples incl. a control | `multi_sample_manifest.csv` + `scrub_plate_map.csv` |
| `launch_mixed.sh` | 2 samples, real depth, via MIXED_INPUT | `two_sample_fixed.csv` |

Submit with
`RUN_NAME=<x> EXTRA_ARGS="<flags>" bsub -q normal -n 2 -M 8000 -R "select[mem>8000] rusage[mem=8000]" -env "all" -o <run>/driver.out bash <launcher>`.
The `-env "all"` matters: without it neither variable reaches the job.

Two traps that cost time, both mine, both avoidable:

- **Give each concurrent run its own launchDir.** Nextflow keeps `.nextflow/` history
  in the launch directory; runs sharing one die with "Unable to acquire lock on
  session". Both launchers `cd $RUN` first.
- **Do not edit a launcher while a job is running it** — bash re-reads the file
  mid-execution and the job dies with exit 127.

### The metagraph species-calling parser was broken in two separate ways

Both were invisible until the numbers were actually read. `metagraph align` reported
**155,181 "species" considered and 0 called** on a sample that is unambiguously
SARS-CoV-2 + influenza A. It was neither a depth problem nor a crash — the counts were
being split across the wrong keys.

1. **Coordinate windows became part of the species key.** In a coordinate index
   metagraph appends the matched window to the label, and a reference record whose own
   name already carries one ends up with two:

       AB847956.1 | Alphainfluenzavirus influenzae:16-166
       KM368312.1 | Alphainfluenzavirus influenzae:1612-1762:2314-2464

   `parse_label` took everything after `" | "` as the species, so every window was its
   own species. Stripped by `COORD_SUFFIX_RE`.

2. **`metagraph query` joins a read's labels with `:` , not `;`.** This is the output
   shape `modules/metagraph_query.nf` had flagged UNVERIFIED. 6364 of 6367 rows were
   multi-label fields being treated as one label. Split by `LABEL_JOIN_RE`, which only
   splits a `:` followed by `<accession> | ` so coordinate windows are left alone.

Effect, same reads each time:

| | before | after |
|---|---|---|
| align: considered / called | 155181 / **0** | 4 / **2** |
| query: considered / called | 6367 / 1 | 22 / **10** |

and the calls are now the right organisms, agreeing with Themisto2/mSWEEP on the same
reads (*Betacoronavirus pandemicum* + *Alphainfluenzavirus influenzae*), each validated
by mapping at 99.4% / 99.96% breadth.

**`min_hits` context:** `metagraph_align_min_hits = 100`. Before the fix the best
evidence for any single key was 60 hits, which is why nothing was called. After
aggregation the same sample gives *Betacoronavirus pandemicum* 572,162. If you ever see
0 called again, check the *number of keys* first — a huge `n_species_considered` is the
signature of a key-fragmentation bug, not of a threshold that needs lowering.

Fixing this exposed a latent one: two species can pick the same best-hit reference
record, which puts a duplicate in the extracted subset FASTA and makes `samtools sort`
die on the header (`Duplicate entry ... in sam header`). Now de-duplicated.

### Real-depth verification

`nf_runs/two_sample_fixed.csv` is a 2-sample, real-depth (19-27MB/mate) MIXED_INPUT
manifest. Run via `launch_mixed.sh`. Results are per-sample distinct and biologically
coherent — influenza for `50376_2_16`, seasonal coronaviruses for `50376_2_78` — with
all called species validated by mapping (96.9% / 100% breadth).

> The manifest originally supplied for this,
> `rvi-viral-metagenomics/minimal_test/two_sample.csv`, cannot be used as-is: all four
> of its FASTQs live under a nextflow work directory that has since been cleaned
> (`.../rvi-viral-metagenomics/50376_2/` no longer exists). `#78` was recovered from
> `10-run-comparative-space/50376_themisto_run/`; `#76` was not found, and `#16`
> substitutes for it. Note also that its `id,R1,R2` header is MIXED_INPUT's format, not
> `parse_mnf`'s `sample_id,reads_1,reads_2` — and `parse_mnf` would reject these sample
> ids anyway, since it forbids the `#` in `50376_2#78`.

### Still open

- **MIXED_INPUT's ENA and iRODS sources** have never been executed (needs credentials
  and real accessions). Its local-reads-manifest source is proven.
- The **`metagraph query` module comment is now confirmed**, but note both metagraph
  methods share `bin/call_metagraph_species.py`; a change there affects both.

---

## Item 3's Metagraph methods are BOTH ported and wired, but only DSL-checked — no farm access this or the last round

Two rounds now with no farm access (sandboxed, no `singularity`, no reference data — same
constraint as the very first handoff). Between them, all three route-map mapping methods
now exist: Themisto2/mSWEEP (**farm-verified**, see the section above) and both Metagraph
methods (**DSL-checked only**, this round added the second one). Treat the Metagraph half
exactly like item 1 was treated before its farm run: believed correct, not proven.

**`metagraph align` (sequence-to-graph alignment)** — `workflows/VIRAL_METAGRAPH_ALIGN.nf`
+ `modules/metagraph_align.nf` (renamed this round from `metagraph.nf`, process `METAGRAPH`
→ `METAGRAPH_ALIGN`, for symmetry with the sibling below). Ported from `eu1/rvi_toolbox.git`'s
`feature_metagraph_align`, merge commit `f23d592` — the fork problem in item 2 still
applies. Gated by `--run_metagraph_align`.

**`metagraph query` (pseudoalignment) — the third method, previously flagged as "not a
real module anywhere," now built.** Read this carefully before touching it further — it
was built deliberately differently from what you'd get by reviving old history:

- The user gave the exact spec: a **new, simpler** module using `metagraph query
  --query-mode labels`, not a revival of the old approach.
- The old approach (recovered from commit `ccae756`'s parent, `eu1/rvi_toolbox.git`) was a
  **two-stage pipeline**: `metagraph align --query-presence --filter-present` (cheap,
  unannotated presence filter) piping into `metagraph query` (annotated, on the filtered
  reads only). Its own commit message says manual testing found it discarded almost every
  real hit — any read with a single sequencing error or SNP never survives the exact-match
  presence filter — and it "found zero" real hits on a real test sample that the
  align-based replacement later found plenty of. **Do not resurrect that two-stage design.**
  If you're ever tempted to "improve" `metagraph_query.nf` by adding a presence-filter
  pre-stage back in, re-read this paragraph first.
- What's actually built (`workflows/VIRAL_METAGRAPH_QUERY.nf` + `modules/metagraph_query.nf`):
  a single `metagraph query --query-mode labels -i <graph> -a <annotation>
  --min-kmers-fraction-label <threshold>` call per mate, directly against the full
  (subsampled) reads — no presence-filter stage. Mirrors `VIRAL_METAGRAPH_ALIGN.nf`'s
  shape otherwise (subsample → run → `CALL_METAGRAPH_SPECIES` → optional
  `METAGRAPH_MAP_QC`), and shares the same graph/annotation/threshold params — same
  reference data, alternative method against it.
- **The exact output shape of `metagraph query --query-mode labels` is unverified** — no
  `metagraph` binary available where this was written, so it couldn't be checked against a
  real invocation. `bin/call_metagraph_species.py` (already proven against `metagraph
  align`'s output) is reused unchanged on the assumption that it scans every
  tab-separated field for a recognizable label shape rather than depending on a fixed
  column count — check this holds for query's actual output on your first real run; if
  `--query-mode labels` output doesn't parse cleanly, that script's `parse_label()`/
  `iter_read_labels()` (see their docstrings) is where to fix it, not
  `metagraph_query.nf` itself.
- A real collision this surfaced, now fixed: `CALL_METAGRAPH_SPECIES` and the map-QC
  processes (`metagraph_species_call.nf`, `metagraph_coverage.nf`) are shared by both
  methods and publish to fixed per-sample paths. If both methods run for the same sample
  they'd have silently overwritten each other's output. Both now take an
  `output_subdir`/`summary_name` value, threaded from each caller (`'metagraph_hits'`/
  `'metagraph_map'` vs `'metagraph_query_hits'`/`'metagraph_query_map'`) — if you add a
  fourth thing that shares these processes, thread a new distinct name the same way rather
  than hardcoding one.
- `GENERATE_MAPPING_REPORT` is now fed by all three methods through one join chain in
  `main.nf`'s `if (params.do_sequence_index) { ... }` block —
  `sequence_index_sample_ch` (built from `preprocessed_3tuple_ch`) is the join backbone,
  each method's counts left-joined with `remainder: true`, defaulted via a named
  `EMPTY_*_COUNTS` constant when that method didn't run or produced no optional output for
  a sample. `count_metagraph_species_hits()`/`count_metagraph_map_qc()` now take a
  `prefix` arg (`'metagraph_align'`/`'metagraph_query'`) so both methods' fields
  (`metagraph_align_n_species_called` vs `metagraph_query_n_species_called`, etc.) can
  merge into one meta without colliding. If a fourth method lands, extend this same chain.

**Reference-data paths are still `null` by default and genuinely unverified** — unlike
every other path param in this file, none of `metagraph_align_graph`/
`metagraph_align_annotation`/`metagraph_align_annotation_seqs`/`metagraph_map_reference_fasta`
were confirmed against real files, across either round. `rvi_toolbox`'s own
`metagraph_align.config` (comments, not necessarily current) suggests starting from eu1's
personal scratch —
`/lustre/scratch126/pam/projects/rvidata/personal/eu1/metagraph/bigviralindex-rvdbc/{viromeindex_clustered_graph.dbg,column_annotation_coordinates.column_coord.annodbg,column_annotation_coordinates.seqs}`
— check it still exists before assuming it does; the equivalent mSWEEP paths in
`rvi_toolbox`'s config turned out to be stale (see item 3's Themisto2 section above), so
don't assume this one is current either. `metagraph_map_reference_fasta` has no candidate
path suggested anywhere — you'll need to work out what reference FASTA actually
corresponds to whichever graph/annotation you end up using. Both methods share these same
four params (one index, two query strategies), so fixing them once fixes both methods.

Your first move on the farm should be exactly what item 1's writeup already describes:
one sample, real containers, `-profile sanger_standard` (or your farm's equivalent) under
`bsub`, confirm the reference-data paths above (correct them in `nextflow.config` once you
know the real ones — don't leave `null` defaults that happen to work only because you
passed the right value on the CLI once), and fix whatever real execution surfaces. Test
`--run_metagraph_align` and `--run_metagraph_query` separately before together — if only
one works, you want to know which. Update this file and the README's "Output layout"
bullet the same way item 1's run did.

---

## Item 4 (abundance estimation lane) is ported and wired, but only DSL-checked — same round as Metagraph pseudoalign, still no farm access

`--do_abundance` (master switch) gates three independently-flagged pieces, all feeding
one `GENERATE_ABUNDANCE_REPORT` call the same join-backbone way the mapping lane's three
methods do:

- **`run_kraken2bracken`** → `workflows/KRAKEN2BRACKEN.nf`, a **viral-lens-owned fork**
  of `rvi_toolbox`'s own `kraken2bracken.nf` (same shared modules, unmodified, included
  directly — `rvi_toolbox/modules/{kraken2,bracken,krakentools}.nf`). The fork exists
  purely because the upstream subworkflow has no `emit:` block at all — nothing it
  produces was reachable from outside it, and both the report and SCRuB need something
  out of it. If `rvi_toolbox`'s `kraken2bracken.nf` ever changes upstream, this file needs
  the same change applied by hand — it's not a re-export, it's a parallel copy.
- **`run_scrub`** → `workflows/SCRUB_DECONTAM.nf` + modules, ported from
  `eu1/rvi_toolbox.git`'s `feature_scrub_decontam` (merge commit `7111c02`) — exists only
  on that fork, same situation as Metagraph (item 2). Runs once per pipeline run against
  `KRAKEN2BRACKEN`'s whole-run `abundance_summary`, requires `--scrub_plate_map` (no
  default, a real user-supplied CSV).
- **`run_abundance_estimation`** → `rvi_toolbox`'s `ABUNDANCE_ESTIMATION`, included
  **unmodified** (not forked) — sourmash/inStrain genome-level profiling against a
  GTDB-style reference set, heavier and less viral-specific than the rest of this lane.
  It also has no `emit:` block, but unlike `KRAKEN2BRACKEN` this one was **not** forked to
  add one — it's wired as a pass-through call (runs, publishes its own files under
  `outdir`, but only contributes an `abundance_estimation_ran: true/false` flag to the
  report, not real per-sample metrics). Deepen this into a real wrapper (same pattern as
  `KRAKEN2BRACKEN.nf`) only once the flag alone is proven insufficient — don't build it
  speculatively.

**A real, previously-latent bug found and worked around, not fixed at the source:**
`rvi_toolbox/subworkflows/abundance_estimation.nf`'s cleanup branch references an
undefined `INSTRAIN` (only `INSTRAIN_PROFILE`/`INSTRAIN_QUICKPROFILE` are ever included)
whenever `cleanup_intermediate_files_abundance_estimation=true` **and**
`bowtie2_samtools_only_abundance_estimation=false` — both the upstream defaults, so this
bug hits every default-configured run of that subworkflow, in both pipelines, and
apparently always has (nothing exercised it before this). `nextflow.config` now defaults
`cleanup_intermediate_files_abundance_estimation` to `false` (not upstream's `true`)
specifically to avoid the branch. If you ever need that cleanup step for real, either
fix it in a viral-lens-owned fork (same pattern as `KRAKEN2BRACKEN.nf`) or fix it upstream
and re-point once the fork situation (item 2) is resolved — don't just flip the flag back
to `true` without one of those, it will crash.

**`results_dir` is now aliased, not left undeclared:** every shared-submodule module this
lane touches (`kraken2.nf`, `bracken.nf`, `krakentools.nf`, `instrain.nf`, `bowtie.nf`,
`sourmash.nf`, `subset_fasta.nf`, `merge_fastq.nf`, `cleanup.nf`) publishes under
`params.results_dir`, an `rvi_toolbox` default viral-lens never declared because nothing
called them before now. `nextflow.config` sets `results_dir = params.outdir` rather than
editing any of those files — everything viral-lens-owned in this lane
(`KRAKEN2BRACKEN.nf`, `SCRUB_DECONTAM.nf` and their modules) uses `outdir` directly and
ignores this alias.

**Also fixed, a gap from the Metagraph round before this one:** `METAGRAPH_ALIGN`/
`METAGRAPH_QUERY` and their downstream processes had no explicit LSF executor override
under `sanger_standard` (which defaults to `executor='local'`) — they'd have silently run
on the submit host. Added alongside this lane's own new overrides.

**Reference-data paths are `null`/unverified, same pattern as every other lane this
round:** `kraken2bracken_kraken2_db` (needs a matching pre-built Bracken kmer-distribution
file alongside it), `genome_file_abundance_estimation`/
`precomputed_index_abundance_estimation`/`stb_file_abundance_estimation` (all three
required together if `run_abundance_estimation` is enabled — none had a working upstream
default either, both pipelines' configs used `""` as their own "must be supplied"
placeholder), `genome_dir_abundance_estimation`/`sourmash_db_abundance_estimation`
(only matter if `sourmash_subset_abundance_estimation=true`, not the default),
`bmtagger_db_abundance_estimation`.

Verified via `-preview`: each of the three abundance sub-flags alone, `run_kraken2bracken`
+ `run_scrub` together, and everything across all four lanes (taxid + assembly + all
three mapping methods + full abundance lane) enabled at once. **None of it has executed a
single real task.**

---

## Item 5 (wider input handling) is wired, but only DSL-checked — same round as items 3/4, still no farm access

Good news first: **this one needed no porting or forking at all.** Unlike items 2-4,
`MIXED_INPUT`, `ENA_DOWNLOAD`, and `DOWNLOAD_FROM_IRODS` already live in viral-lens's own
`rvi_toolbox` submodule (`rvi/rvi_toolbox.git`) — check `ls rvi_toolbox/subworkflows/ |
grep -iE 'mixed|ena|irods'` if you want to confirm this still holds. So the fork problem
(item 2) genuinely doesn't apply here; only wiring was needed.

`--do_mixed_input` (default `false`) gates the whole thing. When off, `main.nf` behaves
exactly as before — `parse_mnf()`, `--manifest`, `sample_id`/`reads_1`/`reads_2` columns,
byte-for-byte unchanged. When on, `MIXED_INPUT()` (no `take:` — it reads `params.*`
directly) replaces it entirely, and its own internal `validate_parameters()`
(`rvi_toolbox/modules/validate_parameters.nf`) decides which of up to three sources
activate, purely from which params are set:

- a local reads manifest via `--manifest_of_reads` (or bare `--manifest`, treated as an
  alias) — **but in `MIXED_INPUT`'s own `id`/`R1`/`R2` column format, not
  `parse_mnf()`'s** `sample_id`/`reads_1`/`reads_2`. These are NOT interchangeable
  manifests; a user switching `--do_mixed_input` on has to reformat their manifest.
- ENA download via `--manifest_ena` (a TSV of run accessions)
- iRODS retrieval via `--studyid`/`--runid`/`--laneid`/`--plexid` (CLI) or
  `--manifest_of_lanes` (a manifest of the same)

Two real gaps found and fixed, both in `main.nf`'s wiring / viral-lens's own code, not the
shared subworkflows themselves:

- **`meta.sample_id` was never set.** All three of `MIXED_INPUT`'s sources
  (`rvi_toolbox/subworkflows/{input_check,ena_input,irods}.nf`) only ever populate
  `meta.id` — none of them know about viral-lens's own `meta.sample_id` convention, which
  every `publishDir` path and report column downstream keys off. Fixed with a `.map{}`
  right after `MIXED_INPUT.out.all_reads_ready_ch` that adds
  `sample_id: meta.id`. If you extend any of the three source subworkflows, remember
  this mapping happens *after* them, not inside — don't duplicate it there.
- **`check_sort_reads_params()` (`workflows/SORT_READS_BY_REF.nf`) unconditionally
  required `--manifest`.** This would have thrown "No manifest provided" on any
  ENA-only or iRODS-only run that never sets `--manifest` at all — even though
  `MIXED_INPUT`'s own `validate_parameters()` already enforces "at least one input
  source" on its own terms. Now skipped entirely when `do_mixed_input` is set. If you
  touch this function again, keep that guard — it's a genuine, easy-to-reintroduce
  regression.

One `includeConfig "./rvi_toolbox/subworkflows/mixed_input.config"` line in
`nextflow.config` pulls in everything this needs (`studyid`/`runid`/`laneid`/`plexid`/
`manifest_ena`/`manifest_of_lanes`/`manifest_of_reads`, plus its own nested
`includeConfig`s of `irods.config` — Sanger-specific `REF_PATH` env var and LSF
`clusterOptions` for `RETRIEVE_CRAM` already baked in, since this config is already
tailored for this exact farm, unlike every path param flagged unverified elsewhere in
this file — and `ena_downloader.config`). All corresponding params were still added to
`nextflow_schema.json` individually (that file has no `includeConfig` equivalent).

Verified via `-preview`: local `--manifest_of_reads`, ENA-only, iRODS-only (via
`--studyid`), the "nothing specified" case (confirms `validate_parameters()`'s own clear
error fires correctly, not a confusing one), and everything across all four lanes at
once. **None of it has executed for real** — ENA needs live network access, iRODS needs
`iinit` auth (interactive login, so check whether your farm session already has a valid
one — it's not something a pipeline run can establish itself) plus the `baton` binary.
Test the local-manifest path first (cheapest to verify), then ENA (network only, no
farm-specific auth), then iRODS last (the most infrastructure-dependent of the three).

---

## Your remaining work, roughly in priority/dependency order

### 1. Prove the assembly lane actually runs -- DONE, see above

Nothing in commits 3-6 has executed for real. On HPC, with real containers and (if
available) real `genomad_db`/`checkv_db`/`vcontact3_db_path` reference data:

```bash
nf-test test tests/workflows/GENERATE_ASSEMBLY_REPORT.nf.test  # should already pass; confirms the environment gap was purely local
nextflow run main.nf --do_assembly true --manifest <a small real/test manifest> \
    --db_path <kraken db> --genomad_db <path> --checkv_db <path> --vcontact3_db_path <path> \
    -profile sanger_standard -resume
```

Follow the plan's original staging advice: **one sample with 2-3 viral scaffolds before
any multi-sample run.** This matters specifically because `VRHYME_BIN` pools scaffolds
*across samples* for its coverage step (see `workflows/VRHYME_BIN.nf` — `POOL_VIRAL_SCAFFOLDS`
→ one bowtie2 index → `COVERM_DEPTH` → per-sample subset) — a bug there won't show up on a
single sample, and a bug elsewhere will be easier to isolate before that cross-sample
machinery is in play.

Specifically check:
- Do `main.nf`'s new `count_genomad_summary()` / `count_vrhyme_membership()` /
  `count_checkv_quality()` / `count_vcontact3_for_sample()` helpers (bottom of `main.nf`)
  actually parse the real TSV/CSV files correctly? They were written against column names
  confirmed by reading `bin/vcontact3_prep.py`/`bin/vcontact3_postprocess.py`
  (`seq_name`/`n_genes` in geNomad's `virus_summary.tsv`; `scaffold`/`bin` in vRhyme's
  `membership.tsv`; `contig_id`/`checkv_quality` in CheckV's `quality_summary.tsv`;
  `Genome`/`genus_prediction` in vContact3's `final_assignments.csv`) but never run against
  real output.
- `count_vcontact3_for_sample()` does naive `.split(',')` CSV parsing (no quoting support).
  If any real column contains an embedded comma, this breaks — swap in a real CSV
  read if so.
- Container/resource labels: everything ported uses `cpu_N`/`mem_N`/`time_N` labels from
  `rvi_toolbox/nextflow-commons.config` (already included). Check these resolve sanely for
  your HPC's actual hardware/queues — they were tuned for the source pipeline's
  environment, not verified against viral-lens's.

### 2. The `rvi_toolbox` fork problem — resolve before building the mapping/abundance lanes

This blocks item 3 below and is worth doing first. Facts as of the prior session (verify
they still hold, forks move):

- `viral-lens/rvi_toolbox` → `rvi/rvi_toolbox.git`, pinned at `4ed291e`d... check
  `git submodule status` for the current pin — was ~39-51 commits behind that remote's own
  `master` (missing vRhyme/vContact3/geNomad-calibration work already merged there — some
  of which commit 3 above just re-ported by hand from the *other* fork instead).
- `rvi-viral-metagenomics-pipeline/rvi_toolbox` → `eu1/rvi_toolbox.git`, was 14 commits
  behind *its* `master`.
- **The SCRuB decontamination subworkflow and the Metagraph alignment subworkflow — both
  needed for the abundance and mapping lanes below — exist ONLY on `eu1/rvi_toolbox.git`**,
  merged to that fork's `master`:
  - SCRuB: branch `feature_scrub_decontam`, merge commit `7111c02`. Adds
    `subworkflows/scrub.nf` (`SCRUB_DECONTAM`, takes a whole-run Bracken summary +
    `params.scrub_plate_map`), `modules/{reformat_bracken,scrub,scrub_heatmap}.nf`,
    `bin/{reformat_bracken_for_scrub.py,plot_scrub_heatmap.R}`.
  - Metagraph: branch `feature_metagraph_align`, merge commit `f23d592`. Adds
    `subworkflows/metagraph_{align,map_qc}.nf`,
    `modules/metagraph{,_coverage,_reference_subset,_species_call}.nf`,
    `bin/{call_metagraph_species.py,aggregate_metagraph_coverage.py}`.
  - Neither is on `rvi/rvi_toolbox.git` at all (checked at the time: not on that fork's
    `master`, not on any of its branches).
- Themisto2/mSWEEP pseudoalignment (needed as one of the mapping lane's three methods) is
  on `rvi/rvi_toolbox.git`'s `origin/feature-msweep` branch — not yet merged to *that*
  fork's `master` either. (`VIRAL_MSWEEP`/`msweep.nf` — using `THEMISTO_PSEUDOALIGN` +
  `MSWEEP` — already exists and runs today in `rvi-viral-metagenomics-pipeline`; the
  `feature-msweep` branch is a rename/cleanup of the same thing on the other fork.)
- Other unmerged branches on `eu1/rvi_toolbox.git` not investigated at all yet — check
  whether any matter before treating the abundance/assembly lanes as complete:
  `feature_mGEMS`, `feature_msweep_map`, `msweep_map_sourmash_ref`,
  `fix/abundance-estimation-bowtie2-index-tuple`, `wip/vcontact3-and-label-fixes`.

You need to decide (this wasn't resolved, just documented — it's a real decision, loop in
a human if the call isn't obvious): does `viral-lens/rvi_toolbox` start tracking
`eu1/rvi_toolbox.git` instead, does `eu1`'s SCRuB/Metagraph work get merged into
`rvi/rvi_toolbox.git` first, or does this integration keep landing everything as
viral-lens-owned files (as commit `c545c24` did) until the forks are reconciled someday
separately? Whichever way, be consistent with whatever the assembly lane already did.

### 3. Build the "map reads to sequence indexes" lane — all three methods now exist; both Metagraph ones need a farm run

Per the route map (`docs/nf-metro/route_map.mmd`, section `seq_index_mapping`): three
methods converging on one QC Mapping step, feeding `GENERATE_MAPPING_REPORT.nf` — which
is now actually wired and running, no longer an orphan. **Your next action here is to run
both Metagraph methods on the farm** — see "Item 3's Metagraph methods..." above for
exactly what that involves and what's unverified.

- **Pseudoalign via Themisto2 — DONE.** `VIRAL_MSWEEP` + `MSWEEP_MAP_QC` are ported as
  viral-lens-owned files (`workflows/VIRAL_MSWEEP.nf`, `workflows/MSWEEP_MAP_QC.nf`,
  `modules/{themisto2,msweep,cleanup,reference_subset}.nf`,
  `bin/{select_reference_records,aggregate_species_coverage}.py`), wired into `main.nf`,
  and **verified end to end on LSF**: THEMISTO_PSEUDOALIGN → MSWEEP → MSWEEP_MAP_QC →
  GENERATE_MAPPING_REPORT all green, with mSWEEP calling *Betacoronavirus pandemicum* at
  0.999 on a SARS-CoV-2 sample and map-QC confirming it at 29.8% breadth.

  Gated by **two** flags, both default `false`: `--do_sequence_index` (lane master
  switch) and `--run_msweep` (this method). One boolean per method, so the Metagraph
  options below can be enabled independently rather than fighting over a single
  method string.

  Reference-data defaults are **corrected, not copied** — `rvi_toolbox`'s
  `msweep.config` points at flat paths under `viromeindex/` that no longer exist (the
  index now lives under a versioned `1.0/` directory), and its `msweep_map_qc.config`
  leaves `msweep_map_reference_fasta` empty ("path TBD"). Now:

  | param | value |
  |---|---|
  | `msweep_themisto_index` | `/data/pam/software/themisto2/viromeindex/1.0/rvdb_clustered_virome.thm2` |
  | `msweep_ref_groups` | `/data/pam/software/themisto2/viromeindex/1.0/rvdb_clustered_virome_species_labels.txt` |
  | `msweep_map_reference_fasta` | `/data/pam/software/themisto2/viromeindex/1.0/data/C-RVDBv32.0.fasta` |

  That FASTA must stay **positionally aligned** with the labels file (record N == line N);
  both currently hold exactly 1321608 entries. If you swap either, re-check that count.

  `THEMISTO_PSEUDOALIGN`'s memory was cut from the upstream `100.GB * task.attempt` to
  `25.GB * task.attempt` at the user's request. Metagraph's own `METAGRAPH` process
  already requests `25.GB * task.attempt` upstream — checked when porting it, no change
  needed. A broader request-vs-actual audit across every module is still wanted later.
- **Sequence-to-graph alignment via Metagraph (`metagraph align`) — ported and wired, not
  farm-run.** `VIRAL_METAGRAPH_ALIGN` + `METAGRAPH_MAP_QC` are in as viral-lens-owned
  files (`workflows/VIRAL_METAGRAPH_ALIGN.nf`, `workflows/METAGRAPH_MAP_QC.nf`,
  `modules/metagraph_align.nf` — renamed from `metagraph.nf` this round —,
  `modules/metagraph{_coverage,_reference_subset,_species_call}.nf`,
  `bin/{call_metagraph_species,aggregate_metagraph_coverage}.py`), gated by
  `--run_metagraph_align`.
- **Pseudoalignment via Metagraph (`metagraph query --query-mode labels`) — ported and
  wired, not farm-run.** The route map used to flag this as *not currently a real,
  maintained module*; it's now built as `workflows/VIRAL_METAGRAPH_QUERY.nf` +
  `modules/metagraph_query.nf`, gated by `--run_metagraph_query`. **Read "Item 3's
  Metagraph methods..." above in full before changing this one** — it's a deliberately
  new, simpler design, specifically NOT a revival of an old two-stage approach that was
  already tried and found to produce ~zero real hits.
- **QC Mapping**: both Metagraph methods reuse the same `METAGRAPH_MAP_QC` subworkflow
  (mSWEEP has its own, `MSWEEP_MAP_QC`) — don't invent a new unified QC step unless
  there's a concrete reason to. `CALL_METAGRAPH_SPECIES` and `METAGRAPH_MAP_QC`'s
  downstream processes now take an `output_subdir`/`summary_name` value precisely so this
  sharing doesn't make align and query overwrite each other's published output — see
  "Item 3's Metagraph methods..." above.

`GENERATE_MAPPING_REPORT` is called once, fed by whichever of the three methods are
enabled (see `main.nf`'s `if (params.do_sequence_index) { ... }` block,
`sequence_index_sample_ch` as the join backbone). If a fourth method ever lands, extend
that same chain rather than building a second report path — that was the whole point of
restructuring it earlier this round. The helpers that populate the meta
(`count_msweep_abundances()`, `count_msweep_map_qc()`, `count_metagraph_species_hits()`,
`count_metagraph_map_qc()` — the last two now take a `prefix` arg so align's and query's
counts don't collide) are at the bottom of `main.nf` beside the assembly ones, each
paired with a named `EMPTY_*_COUNTS` constant for the "this method didn't run, or ran but
produced no optional output for this sample" case — every `.join(..., remainder: true)`
in that block depends on one of those.

### 4. Build the abundance estimation lane — DONE (ported + wired), needs a farm run

See "Item 4 (abundance estimation lane) is ported and wired..." above for the full
picture: `KRAKEN2BRACKEN` (viral-lens fork, adds emits), `SCRUB_DECONTAM` (ported from
`eu1/rvi_toolbox.git`), and `ABUNDANCE_ESTIMATION` (included unmodified, pass-through
report contribution only) all run behind `--do_abundance` + their own sub-flag, feeding
`GENERATE_ABUNDANCE_REPORT`. **Your next action here is to run all three on the farm** —
`run_kraken2bracken` and `run_scrub` first (lighter, more central to this pipeline's
actual purpose), `run_abundance_estimation` second (heavier GTDB/sourmash/inStrain
dependency chain, less proven relevant) — and fix the reference-data paths, none of which
are verified yet.

### 5. Widen input handling — DONE (wired), needs a farm run

See "Item 5 (wider input handling) is wired..." above for the full picture. **Your next
action here is to run all three sources on the farm**, cheapest/most-independent first:
local manifest (via `--do_mixed_input --manifest_of_reads`) to prove the wiring itself,
then ENA (network access only), then iRODS last (needs `iinit` auth + `baton`).
`rvi_toolbox/subworkflows/mixed_input_README.md` documents the exact activation rules if
anything about which param triggers which source is unclear.

### 6. (Deferred, lower priority) Per-scaffold meta granularity

Today, `main.nf`'s assembly-report meta enrichment is **sample-level only**
(`genomad_n_scaffolds`, `vrhyme_n_bins`, etc. — counts, not per-scaffold detail). This was
a deliberate scope decision for the first pass, not an oversight: one report row per
scaffold, carrying that scaffold's own bin/quality/cluster assignment, would need exploding
each sample's TSVs into one channel element per scaffold (mirroring
`workflows/SORT_READS_BY_REF.nf`'s `per_sample_taxid_ch` pattern — synthetic
`meta.id = "${sample_id}.${scaffold}"`), keyed off a scaffold-to-bin/quality/cluster join
that `bin/vcontact3_prep.py` already computes internally (`eligible`, `scaffold_to_bin`,
`scaffold_lengths`, `scaffold_quality`, `promoted_scaffolds` — see that script) but doesn't
currently expose as its own output file. Only take this on once it's clear the sample-level
counts are actually insufficient for whatever report the assembly lane needs to produce —
don't build it speculatively.

## Working conventions established so far — follow these, don't reinvent

- **Before assuming something needs porting or forking from `rvi-viral-metagenomics-pipeline`,
  check whether it already exists in viral-lens's own `rvi_toolbox`** (`ls
  rvi_toolbox/subworkflows/`, `rvi_toolbox/modules/`). Items 2-4 all needed real ports
  because their pieces genuinely only existed on the other fork; item 5 needed none at
  all — `MIXED_INPUT`/`ENA_DOWNLOAD`/`DOWNLOAD_FROM_IRODS` were sitting in viral-lens's
  own submodule the whole time, just never wired up. Don't repeat item 5's initial
  assumption (from the very first version of this file) that it needed the same
  fork-and-port treatment as everything else.
- **New subworkflows/modules land in `viral-lens/workflows/`, `modules/`, `bin/` directly**,
  not inside the `rvi_toolbox` submodule, until the fork situation (item 2) is resolved.
  Module `include` paths stay relative (`../modules/x.nf` from `workflows/`) since both
  repos use the same `workflows/` + `modules/` sibling layout.
- **A shared `rvi_toolbox` subworkflow that already exists but has no `emit:` you need**
  gets forked as a viral-lens-owned file with the same name, same orchestration, same
  shared modules included unmodified — just with an added `emit:` block (see
  `workflows/KRAKEN2BRACKEN.nf`'s header comment for the exact reasoning). This is a
  parallel copy, not a re-export: if the upstream subworkflow changes, apply the same
  change here by hand. Don't fork one just to fix a bug or add a param, though (see
  `ABUNDANCE_ESTIMATION`'s `INSTRAIN` bug, worked around via a `nextflow.config` default
  instead) — forking is specifically for widening the `take`/`emit` contract.
- **After adding any new heavy process** (label `cpu_8`+ or similar), add a matching
  `withName:<PROCESS> { executor = "lsf" }` under the `sanger_standard` profile in
  `nextflow.config` — it defaults to `executor='local'`, so a process without an explicit
  entry runs on the submit host. This was missed once already for Metagraph and had to be
  fixed retroactively (see item 4's writeup above) — check this before considering any
  new lane done, not just when the farm run fails.
- **Report subworkflows** (`GENERATE_*_REPORT.nf`) all reuse
  `modules/write_lane_report.nf` (`write_lane_sequence_summary` + `write_lane_run_summary`)
  and `bin/write_lane_summary.py` — don't write a fourth near-duplicate script; extend the
  shared one if it's missing something.
- **Per-sample publish** for these new lanes uses `modules/publish_lane_report.nf`'s
  `publish_lane_json` (label `lane_output`, publishes under `<outdir>/<sample_id>/` — no
  `.taxid`, unlike the existing `consensus_output` label). Reuse it; don't add a fourth
  near-identical publish process.
- **Multiple calls to the same process** (e.g. `publish_run_files` used for both the
  taxid lane and a new lane) need a distinct `include ... as` alias each — Nextflow DSL2
  rejects invoking one process twice under the same name in one workflow scope. See
  `main.nf`'s `publish_assembly_run_files` alias for the pattern.
- **New params** go in both `nextflow.config` (with a real default, or `null` if there
  isn't a sensible one — see `genomad_db`/`checkv_db`/`vcontact3_db_path`) *and*
  `nextflow_schema.json` (or `validateParameters()` rejects them at runtime).
- **Testing**: one nf-test workflow test per new `GENERATE_*_REPORT.nf`, modeled on
  `tests/workflows/GENERATE_CLASSIFICATION_REPORT.nf.test` (literal `params.meta` map,
  `snapshot(workflow.out).match()`). Stage single-sample fixtures before multi-sample ones,
  specifically because of `VRHYME_BIN`'s cross-sample pooling (see item 1).
- **Sanity-check DSL wiring before assuming it's done**: `nextflow run main.nf -preview
  --do_assembly true --genomad_db /tmp/x --checkv_db /tmp/x --vcontact3_db_path /tmp/x
  --manifest tests/test_data/test_manifests/test_input_manifest.csv --db_path
  tests/test_data/test_kraken_databases/minimal` builds and validates the full DAG without
  executing any task — cheap, catches include-path/duplicate-process/channel-shape bugs
  fast, no containers or real data needed. Do this after any `main.nf` change before a real
  run.
- **Keep the nf-metro route map current** (`docs/nf-metro/route_map.mmd`) as the pipeline's
  shape changes — it's referenced from `README.md` and is the map everyone (human or agent)
  orients from. `docs/nf-metro/README.md` has the authoring/rendering method, including a
  real layout-engine bug we hit and worked around (many-lines-converging-on-one-section) —
  read that before fighting the renderer from scratch again.
- **Commit in small, reviewable increments** with why-focused messages (see `git log` on
  this branch for the tone/format) — this was an explicit ask partway through the prior
  session and is worth continuing.

## Don't do this

- Don't touch `mapping_pipeline_main.nf` — it's the frozen, still-runnable original.
- Don't push/commit into the `rvi_toolbox` submodule's own history without a clear decision
  on item 2 above (which fork it should even go to) — it's shared with other consumers of
  that repo.
- Don't build the per-scaffold meta explosion (item 6) before the sample-level version has
  actually been used and found wanting.
- Don't mark the mapping lane "done" without farm-running both Metagraph methods (item 3)
  — all three exist now, but only Themisto2/mSWEEP has actually executed.
- Don't flip `cleanup_intermediate_files_abundance_estimation` back to `true` without
  either forking `abundance_estimation.nf` to fix the undefined-`INSTRAIN` bug or fixing
  it upstream first (item 4) — it will crash with the upstream default.
- Don't treat `ABUNDANCE_ESTIMATION`'s `abundance_estimation_ran` pass-through flag as
  "good enough forever" without checking — same "don't build ahead of need" reasoning as
  item 6, but also don't let it become permanent by default; revisit once the lane's
  actually been run.
- Don't mix up `--manifest` (parse_mnf()'s `sample_id`/`reads_1`/`reads_2` columns) with
  `--manifest_of_reads` (`MIXED_INPUT`'s `id`/`R1`/`R2` columns, only active when
  `--do_mixed_input true`) — they're different formats for a similarly-named param, not
  aliases of each other, despite `validate_parameters()` treating bare `--manifest` as a
  fallback for `--manifest_of_reads` when `do_mixed_input` is on.
