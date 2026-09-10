# rvi_integration_1 — handoff to continue on HPC

You're picking up mid-flight on branch `rvi_integration_1` of `viral-lens`. This branch is
integrating functionality from a sibling pipeline, `rvi-viral-metagenomics-pipeline`, into
`viral-lens/main.nf`. This file has been updated five times now — first written at
commit `d759968`, updated through `bf71c67` after a farm run (item 1, and the Themisto2
method of item 3), updated through `d0760aa` after a farm-less round ported Metagraph
sequence-to-graph alignment (item 3's second method), updated through `0241641` after a
third, also farm-less round added Metagraph pseudoalignment (item 3's third and last
method), updated through `1536267` after a fourth, still farm-less round built the
entire abundance estimation lane (item 4), updated through `71c15cb` after a fifth round
wired up wider input handling (item 5), and updated again after a farm-less round split
`MAPPING` at the consensus boundary so both classifiers share it (see
"Classifier/consensus split" below — `84d292b`, `HEAD` at the time of this edit).

Note the rounds have alternated between farm-capable and farm-less sessions, and the
sections below are written per round rather than merged — where two sections disagree,
the later one wins and says so explicitly. Items 3, 4 and 5 have since been farm-run (see
"FARM-RUN ROUND"); the classifier/consensus split has not. `git log` is the source of
truth if it's since moved further.

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

A `VERIFY_FASTQ` + `SUBSAMPLE_ITER(initial_subsample_limit)` first-step was tried and then
deliberately reverted (see "Assembly lane: VERIFY_FASTQ + initial-subsample step tried and
reverted" further down) — the user wants to test the lane without it first. Don't
re-add it without checking that section.

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

Reference data, all verified present. Every one of these is now a `nextflow.config`
default, so none has to be passed explicitly for a normal farm run -- pass one only to
point at something else:

| param | path |
|---|---|
| `--db_path` | `/lustre/scratch126/pam/projects/rvidata/pipeline_resources/kraken_databases/production/viral_lens_kdb_v1.5.2` (**now the default**, see below) |
| `--genomad_db` | `/data/pam/software/genomad/genomad_db/genomad_db` |
| `--checkv_db` | `/data/pam/software/ViWrap/CheckV_db/` |
| `--vcontact3_db_path` | `/data/pam/software/vcontact3/` (holds `v232`, `v236`; default version is now 236) |

`db_path` used to default to `null`, which made `--db_path` mandatory and produced only
`ERROR ~ No kraken database path provided` (`check_sort_reads_params()` in
`workflows/SORT_READS_BY_REF.nf`) when it was forgotten. It now defaults to the production
database above, since every real run has used that one. `manifest` is now the only input
with no usable default. `db_library_fa_path` stays `null` on purpose: `SORT_READS_BY_REF`
derives `${params.db_path}/library/library.fna` when it is unset, which is exactly where it
lives under the default database (564 MB, verified). Verified with a `-preview` that passes
**no** `--db_path` at all (`nf_runs/dbdefault/`): `Success: true`, the banner reports the
default path, and the only remaining note is the expected "No db_library_fa_path set,
assuming .../library/library.fna exists" warning.

Two traps that this error tends to travel with, both seen in a real invocation:

- **`--results_dir` is NOT the output flag; `--outdir` is.** `results_dir = params.outdir`
  in `nextflow.config` is an internal alias that exists only so the `rvi_toolbox` submodule's
  processes publish somewhere sensible without editing the shared submodule. Passing
  `--results_dir` overrides the alias but leaves `outdir` at its `$launchDir/results/`
  default, so output **splits**: the 15 `rvi_toolbox` files that publish via `results_dir`
  go where you asked, the 14 viral-lens files that publish via `outdir` do not.
- **nextflow 23.10.1 silently disables parameter validation.** It logs `Nextflow
  self-contained distribution allows only core plugins -- User config plugins will be
  ignored: nf-schema@2.2.0`, which means `validateParameters()` does nothing and a typo'd
  param is accepted in silence; the `>=24.10.3` manifest gate only warns. Use 24.10.6.
  Likewise, without `-profile sanger_standard` the executor stays `local` and every process
  runs inside the single driver job instead of being submitted to LSF.

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

## `main.nf` split into `subworkflows/` — pure refactor, now FARM-VERIFIED as output-identical

Same no-farm-access constraint as every round above. `main.nf` had grown to ~850 lines
because four of the five rvi_integration_1 lanes (mapping/taxid, assembly, sequence-index,
abundance) were inlined directly in its top-level `workflow { ... }` block, each dragging
its own sample-level `count_*()` report helpers along as bottom-of-file `def`s. Only
`PREPROCESSING` was already a proper subworkflow (`rvi_toolbox/subworkflows/preprocessing.nf`).

Reorganized so each lane is its own `take:`/`main:` subworkflow file under the new
`viral-lens/subworkflows/` directory, self-contained (own includes, own `count_*()`
helpers, own PUBLISH calls):

- `subworkflows/classifying_kraken2.nf` — `CLASSIFYING_KRAKEN2`: `SORT_READS_BY_REF`
  (Kraken2 + Kraken2Ref taxid selection and reference resolution), stopping short of
  consensus. Also owns `identified_species_ch`. **Split out of `mapping.nf` later than
  the rest of this reorganization — see "Classifier/consensus split" below.**
- `subworkflows/mapping.nf` — `MAPPING`: `GENERATE_CONSENSUS` → Nextclade → SCOV2
  subtyping → `GENERATE_CLASSIFICATION_REPORT`, driven by *either* classifier. Originally
  this file held the Kraken2 half too, mirroring the still-frozen
  `mapping_pipeline_main.nf` lane byte-for-byte; that half is now
  `classifying_kraken2.nf`.
- `subworkflows/classifying_index.nf` — `CLASSIFYING_INDEX` (was `sequence_index.nf` /
  `SEQUENCE_INDEX`): the three mapping methods (Themisto2/mSWEEP, Metagraph align,
  Metagraph query) → `GENERATE_MAPPING_REPORT`, plus the
  `count_msweep_*()`/`count_metagraph_*()` helpers and `EMPTY_*_COUNTS` constants
  referenced in item 3 above. Produces species **calls**, not consensus sequences, and
  hands them to `MAPPING`. It does still map reads for validation inside its map-QC step —
  see the note under "Classifier reports, MAPPING maps".
- `subworkflows/assembly.nf` — `ASSEMBLY`: `ASSEMBLE_META` → `GENOMAD_CLASSIFY` →
  `VRHYME_BIN`/`CHECKV_QC` → `VCONTACT3_RUN` → `GENERATE_ASSEMBLY_REPORT`, plus the
  `count_genomad_summary()`/`count_vrhyme_membership()`/`count_checkv_quality()`/
  `count_vcontact3_for_sample()` helpers referenced in item 1 above.
- `subworkflows/abundance.nf` — `ABUNDANCE`: `KRAKEN2BRACKEN`/`SCRUB_DECONTAM`/
  `ABUNDANCE_ESTIMATION` → `GENERATE_ABUNDANCE_REPORT`, plus `count_bracken_species()` and
  `EMPTY_BRACKEN_COUNTS`, referenced in item 4 above.

`main.nf` itself is now ~300 lines: the log banner, `--do_mixed_input`/`parse_mnf()` input
handling, the `PREPROCESSING` call, and one gated call per lane
(`CLASSIFYING_KRAKEN2` and `MAPPING` unconditionally,
`ASSEMBLY`/`CLASSIFYING_INDEX`/`ABUNDANCE` behind their existing `params.do_*` flags). At
the time of *this* reorganization no process, channel-shape, or param-gating logic changed
— it was a pure move, verified by grepping `main.nf` for direct calls to any of the moved
processes (`SORT_READS_BY_REF(`, `ASSEMBLE_META(`, `KRAKEN2BRACKEN(`, etc.) and confirming
none remain. The classifier/consensus split below came afterwards and did change wiring.

One real (behavior-neutral) simplification made along the way: previously, when
`params.do_preprocessing` was `false`, the taxid-mapping lane's input
(`sort_reads_in_ch`) came straight from the un-coerced `reads_ch` (raw manifest strings),
while every other lane already used the `file()`-coerced `preprocessed_3tuple_ch`. That
asymmetry is gone — `MAPPING` now derives its input from `preprocessed_3tuple_ch` like
every other lane, so all four lanes share exactly one upstream channel. Nextflow's `path`
input coercion means this should be a no-op for `SORT_READS_BY_REF`, but call this out
specifically if anything in the taxid/consensus lane behaves differently after a real run.

### Verification: DONE, farm-run and diffed against pre-refactor baselines

This refactor has now been verified for real, not just DSL-checked. Two things were
checked: that the moved code is textually the same logic, and that real runs produce the
same output as runs from before the split.

**Static check.** Each old inlined lane block in `fbe08d3:main.nf` was diffed against the
corresponding `fde04c7:subworkflows/*.nf` file, normalizing indentation, blank lines and
comments. All four lanes are identical in logic: the only per-lane body change is the
`if (params.do_*) {` wrapper becoming the subworkflow header (gating moved to `main.nf`),
plus the `count_*()` helpers moving along with their lane. Note the mapping lane's four
`publish_*` calls only *look* new in such a diff — they sat at the bottom of the old
top-level `workflow {}` (old lines 494-498), after the abundance block, not inside the
mapping section. Beyond the documented `preprocessed_3tuple_ch` change, exactly two other
real differences exist, both benign:

- A `.toString()` was added to `sample_report_with_join_key_ch`'s join key. This looked
  dangerous: `meta.id` is built by interpolation at `workflows/SORT_READS_BY_REF.nf:145`,
  so it is a **GString**, and at the Groovy level `GString.equals(String)` is `false` with
  a different `hashCode()` (`37 + s.hashCode()`), which would silently empty the join that
  feeds `GENERATE_CLASSIFICATION_REPORT`. Probed directly on LSF
  (`nf_runs/gstring_probe/probe.nf`, kept for future reference): **Nextflow normalizes
  join keys**, and GString-keyed and String-keyed joins both match. Harmless — but don't
  "tidy" GString/String key handling elsewhere assuming raw Groovy semantics apply.
- `if (!params.do_scov2_subtyping == true)` became a proper `else`. Only differs when
  `do_scov2_subtyping` is truthy but not `true`, where the old code left
  `scov2_subtyped_ch` undefined. A fix, not a regression.

**Also checked: the new-species consensus feature is genuinely inert when off.** Its
consumer is gated on `params.call_consensus_for_new_species`, but `identified_species_ch`
is built **unconditionally** in `MAPPING`, and its `.map{}` would throw on a
header-with-no-rows pre-report. That case cannot arise: `bin/k2r_report.py` (lines
~200-206) builds the frame from an empty dict when nothing is selected, so a zero-hit
sample yields a 1-byte, column-less file that the `size() > 1` filter drops before the map
ever sees it. If `k2r_report.py` is ever changed to emit real column headers in the empty
case, that map needs a row-count guard.

**Real farm runs.** All submitted under `bsub` (the user's standing rule: never run
`main.nf` on the head node, not even `-preview`). Driver harness: `nf_runs/env.sh` plus the
new `nf_runs/launch_case.sh`, which takes `RUN_NAME`/`MANIFEST`/`EXTRA_ARGS` and gives each
run its own launchDir (concurrent runs sharing one launchDir collide on the session lock).

| case | run dir | result |
|---|---|---|
| all four lanes + `--run_kraken2bracken`, `-preview -with-dag` | `nf_runs/preview_dag/` | `Success: true` |
| `--do_preprocessing false`, single sample | `nf_runs/nopreproc/` | `Success: true`, 2m57s |
| `--do_assembly --do_sequence_index --run_msweep`, single sample | `nf_runs/regress_post/` | `Success: true`, 16m57s |
| `--do_assembly`, 3 samples | `nf_runs/assembly_multi_post/` | `Success: true`, 14m44s |

The last two deliberately re-ran commands **parameter-identical** to the pre-refactor
baselines already sitting in `nf_runs/regress/` and `nf_runs/assembly_multi/` (verified:
the `genomad_db`/`checkv_db`/`vcontact3_db_path` values passed explicitly match the
`nextflow.config` defaults those baselines picked up), so outputs could be diffed rather
than merely exit-code-checked:

- **`regress_post` vs `regress`**: identical file inventory (70/70), all 15 sequence
  outputs byte-identical (8 per-taxid consensus FASTAs, metaspades contigs+scaffolds,
  geNomad `virus.fna`/proteins, vRhyme bin `.fasta`/`.faa`/`.ffn`), 35/38 text reports
  byte-identical. The 3 diffs: `mapping_run_summary.json` and `mapping_summary_report.csv`
  differ *only* by the added `new_species_candidates_n: 0` column (the opt-in feature
  contributing its count even when off — every pre-existing column matched exactly,
  including `msweep_top_group`/`0.999273` and `mapqc_max_breadth_pct`/`29.7918`); and
  `_mSWEEP_probs.tsv` differs on 6 of 7328 lines by a last-significant-digit float jitter
  (`2.56821e-06` vs `2.56822e-06`) from mSWEEP's own EM solver — its `abundances` output,
  which is what actually feeds the report, is identical.
- **`assembly_multi_post` vs `assembly_multi`**: identical file inventory (142/142), all
  sequence outputs byte-identical, `assembly_summary_report.csv` and
  `assembly_run_summary.json` **identical**, 65/67 text files identical. The 2 diffs are
  both tool-internal float jitter on the Zeptometrix sample (CheckV completeness
  `93.06734216838821` vs `...18`, ~1e-14; geNomad virus score `0.9720` vs `0.9721` on one
  365 bp scaffold) and neither perturbs any `count_*()` value.

This multi-sample case is the one that exercises `VRHYME_BIN`'s cross-sample scaffold
pooling (item 1's specific worry), and it came out clean.

**Side benefit — item 1's open question about the assembly `count_*()` helpers is now
answered.** They parse real geNomad/vRhyme/CheckV/vContact3 output correctly, across three
samples including the awkward case (Zeptometrix: `vrhyme_n_bins=0` but
`checkv_n_high_quality=1`, i.e. the unbinned-scaffold path):
`genomad_n_scaffolds` 12/8/7, `genomad_n_eligible` 12/8/7, `vrhyme_n_bins` 1/1/0,
`vrhyme_n_binned_scaffolds` 6/6/0, `vcontact3_n_genomes` 7/3/7. The
`count_vcontact3_for_sample()` naive-`.split(',')` caveat still stands in principle — no
embedded comma has turned up in real `final_assignments.csv` output yet.

Also worth noting: `nextflow lint subworkflows/abundance.nf` crashes with an internal
parser error (`Range [51, 52) out of bounds for length 51`), traced to a `try`/`catch`
inside a `.each{}` closure sitting in a top-level `def` (copied verbatim from the original
`main.nf`, so pre-existing, not introduced here) — that's a bug in the separate `nextflow
lint` static-analysis tool, not in the pipeline; `nextflow run`/`-preview` (the actual
compiler) handle it fine. Don't spend time trying to "fix" that code to satisfy the linter.

## Assembly lane: VERIFY_FASTQ + initial-subsample step tried and reverted

**Current status: NOT in the code.** This was added, DSL-checked (see "Verification so
far" below), then explicitly reverted at the user's request — they want to test the
assembly lane without it first. `subworkflows/assembly.nf` is back to calling
`ASSEMBLE_META(preprocessed_3tuple_ch)` / `VRHYME_BIN(..., preprocessed_3tuple_ch)`
directly; the three new params (`initial_subsample_limit`/`minimum_fastq_reads`/
`fastq_error_handling_mode`) were removed from `nextflow.config` and
`nextflow_schema.json` again. The rest of this section is kept as a record of what was
tried and why, in case a real farm run later shows it's actually needed (see the gap it
was meant to fix, below) — **don't re-add it without checking with the user first.**

The gap this would have fixed is specifically about large/misbehaving inputs, which is
exactly what a `-preview` DAG check cannot exercise (it never runs a real process, so it
can't observe metaspades/preprocessing actually choking on an oversized fastq) — worth
revisiting if a real farm run hits that problem.

`rvi-viral-metagenomics-pipeline/main.nf` runs `VERIFY_FASTQ` then `SUBSAMPLE_ITER`
(capped at `initial_subsample_limit`, 10M read pairs by default) as the very first thing
in its workflow, before `PREPROCESSING` even starts — a blanket safety cap so nothing
downstream (preprocessing, assembly, anything) chokes on a pathologically large or
unreadable input:

```groovy
MIXED_INPUT
| VERIFY_FASTQ

initial_subsample_limit_ch = Channel.value( params.initial_subsample_limit )
SUBSAMPLE_ITER(VERIFY_FASTQ.out.verified_fastq_ch, initial_subsample_limit_ch)

SUBSAMPLE_ITER.out.final_read_channel
.set{ capped_reads_ch }
```

`viral-lens/workflows/ASSEMBLE_META.nf` (ported in item 1) already has its own, separate,
much tighter `SUBSAMPLE_ITER` call internally (`metaspades_subsample_limit`, 500k reads —
this one exists specifically because metaspades assembly quality/tractability needs a much
lower depth than general use), but the *earlier*, blanket 10M-cap pass from the reference
pipeline's top-level workflow was never ported at all — nobody had hit a real input large
enough to need it yet.

**What was tried** (`subworkflows/assembly.nf`, top of `main:`, before `ASSEMBLE_META` is
called): `VERIFY_FASTQ(preprocessed_3tuple_ch)` → `SUBSAMPLE_ITER(...,
initial_subsample_limit_ch)` → `capped_reads_ch`, which now feeds both `ASSEMBLE_META`
(which still does its own further 500k-cap internally, unchanged) and `VRHYME_BIN`'s raw
reads argument (previously `preprocessed_3tuple_ch` directly) — matching the reference
pipeline, where `VRHYME_BIN`'s equivalent (`ready_reads_ch`) is also post-10M-cap, not
post-500k-cap.

**Deliberately scoped to the assembly lane only**, not moved to the top of `main.nf`
ahead of `PREPROCESSING` the way the reference pipeline has it: `preprocessed_3tuple_ch`
in viral-lens is shared by three other, already-farm-verified lanes (mapping,
sequence-index, abundance — see items 1/3/4/5 above) that were never built expecting this
cap. Applying it globally would silently change their input the moment `--do_assembly` is
turned on for an unrelated run. If a real farm run later shows the other lanes need the
same blanket protection (e.g. a genuinely oversized or corrupt input reaching
`PREPROCESSING` itself, before assembly even starts), that's a deliberate follow-up
decision, not something to add reflexively — loop in a human, same as the `rvi_toolbox`
fork-problem decision (item 2 below).

**New params** (`nextflow.config` + `nextflow_schema.json`, both updated, defaults copied
verbatim from the reference pipeline's own `nextflow.config`): `initial_subsample_limit`
(10000000), `minimum_fastq_reads` (0), `fastq_error_handling_mode` (`'exit_on_error'` —
one of `exit_on_error`/`only_unreadable`/`ignore`, see `VERIFY_FASTQ`'s own header comment
in `rvi_toolbox/subworkflows/verify_fastq.nf` for exact semantics). `subsample_iterations`/
`subsample_seed` already existed (pre-staged for this exact purpose, per the comment
already sitting above them in `nextflow.config` before this round touched it) and needed
no change. `VERIFY_FASTQ`/`SUBSAMPLE_ITER` both already lived in viral-lens's own
`rvi_toolbox` submodule (same "check before porting" pattern as item 5) — no fork/port
needed, only wiring.

**Verification so far**: `nextflow run main.nf -preview --do_assembly true --genomad_db
/tmp/x --checkv_db /tmp/x --vcontact3_db_path /tmp/x --manifest
tests/test_data/test_manifests/test_input_manifest.csv --db_path
tests/test_data/test_kraken_databases/minimal` — DAG builds (`Success: true`), and the
process list visibly shows **two** distinct `SUBSAMPLE_ITER:SUBSAMPLE_SEQTK` entries (the
new blanket one, and `ASSEMBLE_META`'s pre-existing tighter one), confirming both stages
are wired as separate, sequential caps rather than one accidentally shadowing the other.
Also re-ran with all four lanes enabled together (`--do_assembly --do_sequence_index
--run_msweep --do_abundance`, `run_kraken2bracken`/`run_abundance_estimation` left off to
dodge the pre-existing, already-documented `kraken2bracken_kraken2_db=null` gap from item
4) — same clean `Success: true`. **None of this executed against a real fastq file before
being reverted.**

**If this is ever re-added**: test with at least one sample whose read count is high
enough to exercise `initial_subsample_limit` (or, cheaper, temporarily set
`--initial_subsample_limit` below your test sample's real read count to force the cap to
trigger even on a small file) — the whole point of this step is behavior that only shows
up above that threshold, which no DSL/`-preview` check can exercise at all. The user's
current plan is to farm-test the assembly lane without this step first and add it back
later only if a real run actually needs it.

## CLASSIFYING_INDEX -> MAPPING cross-classifier: consensus for species Kraken2 missed — DSL-checked only, needs a real overlap/disagreement case to prove anything

**Not farm-verified, and can't be meaningfully checked with `-preview`** — the entire
point of this feature is behavior that only shows up when a real sample has a species
Kraken2 misses but Themisto2/Metagraph catch with real breadth. A DSL/`-preview` check only
proves the wiring resolves, not that the logic is right. **Your first real run after
pulling this should specifically include (or synthesize) a sample where Kraken2 and a
sequence-index method disagree** — e.g. a species just below `min_reads_for_taxid` or
outside the Kraken2 db, but present in the Themisto2/Metagraph reference index with decent
depth.

**What it does**: when `--call_consensus_for_new_species true`, `CLASSIFYING_INDEX`
resolves a reference record for every species a sequence-index method calls with
`breadth_pct` above `new_species_min_breadth_pct` (default 10.0) and hands them all to
`MAPPING`; `MAPPING` then drops the ones `CLASSIFYING_KRAKEN2` already found for that
sample and consensuses the rest. Off by default — the cross-classifier dependency is
`MAPPING` taking `CLASSIFYING_KRAKEN2.out.identified_species_ch` as a fifth input.

**Classifier reports, MAPPING maps — the division of labour, arrived at in two steps.**
Originally `CLASSIFYING_INDEX` did the filtering, the reference resolution AND the read
pairing itself. Both moved out, for the same reason each time: `MAPPING` is what spends a
consensus, so it should decide what gets one.

- The species-level filter moved first. The invariant "one consensus per (sample, species)"
  now holds structurally for anything reaching `MAPPING`, rather than resting on each
  classifier remembering to filter itself, and a classifier no longer needs to know what a
  *different* classifier found.
- Consensus-reference resolution and read pairing followed. `CLASSIFYING_INDEX` now emits
  species + the reference record its own map-QC already picked + the supporting counts, and
  nothing else. This also removed the wasted work the first step created, where references
  were resolved for species about to be discarded as already-known.

  **SUPERSEDED — the paragraph below described the state of one commit only.** Map-QC was
  removed in the next round, so `CLASSIFYING_INDEX` now maps nothing and its calls carry no
  breadth at all; the threshold moved after `GENERATE_CONSENSUS`. See "Map-QC removed from
  both classifiers" below for what replaced this. Kept here because the distinction it
  draws — VALIDATION vs CONSENSUS mapping — is what the removal traded away, and is the
  thing to re-read before restoring either `*_MAP_QC` subworkflow.

  `CLASSIFYING_INDEX` *did* map reads: its map-QC step (`THEMISTO_MAP_QC` /
  `METAGRAPH_MAP_QC`) runs `INDEX_REFERENCE_FASTA` → `EXTRACT_REFERENCE_SUBSET` →
  `BOWTIE_INDEX` → `BOWTIE2SAMTOOLS` → `SAMTOOLS_COVERAGE`, which is exactly where the
  `breadth_pct` it reports comes from. That mapping is load-bearing, not incidental:
  breadth is what `new_species_min_breadth_pct` thresholds on, and what separates a real
  call at ~99% breadth from index noise at 3-14%. So there are two kinds of mapping in
  play — VALIDATION mapping (bowtie2, in the classifier, to decide whether a call is real)
  and CONSENSUS mapping (`params.read_aligner` + iVar, in `MAPPING`) — and a surviving
  sequence-index species goes through both. Reusing the bowtie2 BAM for consensus instead
  would build these consensuses differently from every Kraken2-side one, making the two
  incomparable, so the second pass is deliberate.

So: `MAPPING` prefers Kraken2's species and references where both classifiers agree, and
resolves + maps only the species Kraken2 missed. `SELECT_REFERENCE_RECORD_BY_NAME` fell
out of use as a result and is marked UNUSED rather than deleted — see its header for the
two situations that would bring it back.

**Two costs of this arrangement, both deliberate:**

- `MAPPING` runs its own `INDEX_REFERENCE_FASTA` pass over the reference FASTA, so the
  `SEQIDX_<n>` ids the classifier reported are valid to grep. The map-QC step inside the
  classifier already made one such pass over the same file, so an enabled run makes two
  (verified: lane off → 0 passes, lane on but feature off → 1, feature on → 2). Reusing
  the classifier's indexed output instead would mean plumbing it through `CLASSIFYING_INDEX`
  → `MAPPING`, which is undefined whenever the lane is off; a self-contained `MAPPING` was
  judged worth one extra pass.
- `params.call_consensus_for_new_species` is checked in *both* subworkflows. Redundant on
  paper — the classifier already emits nothing when it is off — but without the check in
  `MAPPING` that `INDEX_REFERENCE_FASTA` pass would run on every default run, and every run
  would require `msweep_map_reference_fasta` to exist.

Consequence for the sequence-index report: its `new_species_consensus_n` column is now
`new_species_candidates_n`, counting what the classifier can honestly know (species called
above the breadth threshold and reported, pre-filter). The post-filter truth is in the
classification report, whose per-consensus meta now carries `discovered_by:
'sequence_index'`, `discovered_by_method` (which of the three methods called it), plus
`index_reference_record`, `index_hit_count` and `index_breadth_pct` — so what a sequence
index actually contributed, and on what evidence, is visible per consensus rather than only
as a per-sample count.

**Restructured since it was first written (see "Classifier/consensus split" below).** It
originally ran its *own* `GENERATE_CONSENSUS` call inside the sequence-index lane and
published a bare consensus — no Nextclade, no SARS-CoV-2 subtyping, no classification
report row. Now both classifiers hand over to one shared `MAPPING`, so a species only a
sequence index found gets the same treatment as a Kraken2-found taxid. The
species-identity reconciliation described below is unchanged by that move; only where the
consensus runs changed.

**SECOND CORRECTION — the dedup only checked the first candidate per sample (fixed in
`4973f35`).** Separate from, and later than, the wrong-column bug documented below. The
filter used `join(identified_species_ch, remainder: true)`, but Nextflow's `join()` pairs
matching keys **one-to-one** — it does not broadcast one right-hand element across many
left-hand ones. There are many candidates per sample and exactly one identified-species
list, so every candidate after the first got a `null` right-hand side (courtesy of
`remainder: true`) and passed through unchecked. Themisto2 called 9 species on one real
sample, so this was the normal case.

Two things follow for you:

- **It would have looked exactly like a regression of the wrong-column bug** — redundant
  consensuses for species Kraken2 already found — via a completely different mechanism.
  If you see that symptom on the farm, check both.
- **`combine(by:)` broadcasts, `join()` does not.** Worth remembering generally: this
  codebase joins many-per-sample against one-per-sample channels in several places. The
  fix also had to complete the right-hand side first (`identified_by_sample_ch`,
  defaulting to `[]`), because `combine(by:)` is an inner join and would otherwise
  silently drop candidates for samples with no Kraken2 pre-report at all — the opposite
  of the intended passthrough. Verified both operators against the real channel shapes in
  a standalone script; don't swap them back without doing the same.

**Species-identity reconciliation**: Kraken2 taxids and the mSWEEP/Metagraph RVDB-index
labels are unrelated numbering schemes with no shared numeric ID, so comparison is by
normalized (`.trim().toLowerCase()`) species-name text.

**CORRECTED — the original version of this matched on the wrong column, and did so in
every case.** The claim that used to sit here, that `ref_selected` "is the only thing on
the MAPPING side comparable to the sequence-index methods' `species_label`/`species`
fields", is false. `bin/k2r_report.py` writes *two* free-text name columns into the
pre-report, at different ranks:

| column | derived from | rank | example |
|---|---|---|---|
| `virus_name` | kraken2ref's `source_taxid` | **S** | `Betacoronavirus pandemicum` |
| `ref_selected` | the selected reference taxid | S1/S2/S3 | `Severe acute respiratory syndrome coronavirus 2` |

mSWEEP and Metagraph label the RVDB index with ICTV species binomials, i.e. `virus_name`'s
vocabulary — never `ref_selected`'s. Matching on `ref_selected` alone therefore never
matched *anything*: every species MAPPING had already found looked "new", so the feature
called a redundant second consensus for it instead of a new-species one. Farm-proven on
`nf_runs/newspecies_after/`'s sample, where mSWEEP calls *Betacoronavirus pandemicum* at
29.79% breadth and MAPPING had already found it with 7297 reads under taxid 2697049.

`virus_name` is reliably species rank, not incidentally so: kraken2ref decomposes
species -> below-species by construction and its decomposed JSON carries the rank code in
its own `source` field (`source = [12, 'S']`, `target = [13, 'S1']`,
`path_as_taxids = [3418604, 2697049]` for that SARS-CoV-2 row). Checked across every
`*_decomposed.json` from all runs to date: **68 of 68 reference selections have `source`
rank `S`.**

`identified_species_ch` now matches on the **union** of `virus_name` and `ref_selected`,
not `virus_name` alone. That is deliberate: a sequence-index label that happens to match a
strain-level reference Kraken2 already selected is also genuinely already covered, and
adding entries to the "already identified" set can only ever suppress a candidate, never
invent one.

## New-species consensus: it also crashed unless ALL THREE methods were enabled

Second, independent bug in the same feature, found the moment it was first run with a
method subset. `subworkflows/classifying_index.nf`'s new-species block reached directly into
`VIRAL_METAGRAPH_ALIGN.out.map_qc` / `VIRAL_METAGRAPH_QUERY.out.map_qc`, past the
`if (params.run_*)` guards. A subworkflow that was never invoked has no `.out`, so
`--call_consensus_for_new_species true --run_msweep true` (no Metagraph) aborted in 19
seconds with:

```
Access to 'VIRAL_METAGRAPH_ALIGN.out' is undefined since the workflow
'VIRAL_METAGRAPH_ALIGN' has not been invoked before accessing the output attribute
```

The block's own comment asserted the opposite — "Each method's channel is already
`Channel.empty()` when that method is off, so nothing extra to gate here" — which
conflated the per-method `*_counts_ch` variables (which *are* set to `Channel.empty()` in
each `else`) with the subworkflow output accessor, which is not a channel at all until the
subworkflow runs. **This is exactly why the feature's original `-preview` verification
missed it: that check enabled all three methods at once**, the one combination where the
bug is invisible. Fixed by capturing each method's `map_qc` into an
`msweep_map_qc_ch`/`metagraph_align_map_qc_ch`/`metagraph_query_map_qc_ch` variable inside
its own `if/else` (matching how the file already handles the counts channels) and having
the new-species block consume those. Verified by `-preview` on both mSWEEP-only and
metagraph-query-only (`nf_runs/preview_fix/`, `nf_runs/preview_fix_mgonly/`).

**Lesson worth generalizing: when a feature is gated by N independent method flags, a
single all-flags-on `-preview` is not coverage.** Check at least one single-method subset
too — it is cheap and it is where `.out`-access bugs live. If `metagraph_align_names_dmp` ever gets a real
`names.dmp` (currently the `assets/NO_NAMES_DMP` placeholder), Metagraph's labels could
start being taxid-derived, which wouldn't change what MAPPING emits but could change
whether Metagraph's own species text still matches mSWEEP's/Kraken2's — worth a re-check
if that param is ever populated for real.

**Synchronization is per-sample, not batch-wide** — this was an explicit ask (confirmed
with the user: "just the current sample finish is okay"). `workflows/SORT_READS_BY_REF.nf`
now also emits `raw_sample_pre_report_ch` (the per-sample pre-report *file*, before it
gets exploded via `.splitCsv()` into the existing `sample_pre_report_ch`). Reading that
file directly in `subworkflows/mapping.nf` (`identified_species_ch`, plain
`.readLines()`, same style as every `count_*()` helper elsewhere in this codebase) gives
a per-sample "already identified" set with no `.groupTuple()`/`.collect()` — those
operators can't emit a group until their *whole* upstream channel closes, which would
mean waiting for Kraken2/k2r to finish for every sample in the run before any sample's
new-species consensus could start. `.unique()` (used later to de-dupe candidate new
species) doesn't have this problem — it streams, emitting each non-duplicate immediately.

**New files**:
- `bin/select_reference_record_by_name.py` + `modules/select_reference_record_by_name.nf`
  (`SELECT_REFERENCE_RECORD_BY_NAME`): a trimmed sibling of
  `bin/select_reference_records.py`/`SELECT_REFERENCE_RECORDS` (mSWEEP's own low-abundance
  validation, `modules/reference_subset.nf`) — same "longest sequence for this label"
  rule, minus the abundance-threshold gating (the caller already knows exactly which
  species it wants). Reuses `INDEX_REFERENCE_FASTA`/`EXTRACT_REFERENCE_SUBSET` from
  `modules/reference_subset.nf` unchanged, always against `msweep_ref_groups`/
  `msweep_map_reference_fasta` regardless of whether `run_msweep` is on — those params
  always have real defaults.

**Deliberately no cross-sample de-duplication**: a species found "new" in many samples
gets its reference extracted once *per sample*, not once per run — keeps everything
scoped per-sample (matching the no-batch-wait decision above) at the cost of some
redundant `seqkit`/extraction work if the same new species turns up across many samples.
Revisit only if that redundancy proves to actually matter on a real multi-sample run.

**Publish/report**: reuses `publish_consensus_files` (own `include ... as
publish_new_species_consensus_files` alias) — lands at
`${outdir}/${sample_id}/mapping/<slugified species name>/`, alongside real-Kraken2-taxid
results. `meta.taxid` for these is a filesystem-safe slug of the species name, **not a
real Kraken taxid** — confirmed acceptable with the user, since there isn't a real one
(these species were never Kraken2-sorted). A `new_species_candidates_n` count feeds into
the existing `mapping_report_prep_ch` join chain the same way every other sequence-index
count does (`GENERATE_MAPPING_REPORT`/`write_lane_summary.py` are schema-free, no changes
needed there).

**Known, deliberately out-of-scope wrinkle**: `VIRAL_MSWEEP.nf`'s own map_qc validation
is driven by mSWEEP's probabilistic abundance estimate (`SELECT_REFERENCE_RECORDS` keys
off `MSWEEP.out.abundances`), unlike the Metagraph methods, which call species directly
from raw pseudoalignment hit counts (`CALL_METAGRAPH_SPECIES`). The eventual intent
(flagged by the user, not yet built) is to restructure this the same way — Themisto2
pseudoalignment-hit-driven species calling, mSWEEP demoted to an optional add-on
estimate — but that needs a new hit-count caller script (parsing Themisto2's raw
pseudoalignment format, mirroring `bin/call_metagraph_species.py`) that doesn't exist
anywhere: checked viral-lens's own history, `rvi-viral-metagenomics-pipeline`'s `main`,
and its `feature_mGEMS`/`msweep_map_sourmash_ref` branches (and their `rvi_toolbox`
submodule pins) — none of them do this; `feature_mGEMS` replaces map_qc with mGEMS
read-binning instead (a different, unrelated approach), and `msweep_map_sourmash_ref`
only changes *which* reference record gets picked (sourmash ANI vs. longest sequence),
still triggered by mSWEEP's abundance output. **Decided out of scope for this feature** —
this cross-lane consensus work deliberately uses today's mSWEEP-abundance-driven
`map_qc` as its candidate source. If/when that Themisto2 restructuring happens, re-check
this feature's mSWEEP-sourced candidates too.

**Verification so far**: `nextflow run main.nf -preview --do_sequence_index true
--run_msweep true --run_metagraph_align true --run_metagraph_query true
--call_consensus_for_new_species true --manifest
tests/test_data/test_manifests/test_input_manifest.csv --db_path
tests/test_data/test_kraken_databases/minimal` — `Success: true`; confirmed via
`-with-dag` (the live progress display truncates too aggressively to show every process
name) that `CLASSIFYING_INDEX:SELECT_REFERENCE_RECORD_BY_NAME`,
`CLASSIFYING_INDEX:INDEX_REFERENCE_FASTA`, `CLASSIFYING_INDEX:EXTRACT_REFERENCE_SUBSET` and
the `MAPPING:GENERATE_CONSENSUS:*` process chain (one invocation, shared with the Kraken2
classifier since the split below) all appear correctly in the
resolved DAG, alongside `MAPPING`'s own separate copy of the same `GENERATE_CONSENSUS`
processes. **None of this has executed against real data or a real species
disagreement** — the whole point of this feature is unverifiable without one.

## Themisto2 is now the sequence-index lane's DEFAULT method (ported, farm-run)

Ported from `eu1/rvi_toolbox.git`'s `feature_msweep_map` branch (NB: `feature_msweep_map`,
not `feature-msweep`) as viral-lens-owned files, same as every other lane here. This is
the restructuring an earlier round flagged as "the eventual intent (flagged by the user,
not yet built)": species called directly from Themisto2 pseudoalignment read-hit counts,
with mSWEEP demoted to an optional add-on.

**That earlier round's claim that the hit-count caller "doesn't exist anywhere" was wrong.**
It checked `feature_mGEMS` and `msweep_map_sourmash_ref` but not `feature_msweep_map`,
which has both `bin/call_themisto_species.py` and the whole surrounding subworkflow. This
was a port, not a build.

New files: `workflows/VIRAL_THEMISTO_MSWEEP.nf` (from `subworkflows/themisto2-msweep.nf`),
`workflows/THEMISTO_MAP_QC.nf` (from `subworkflows/themisto_map_qc.nf`),
`modules/themisto_species_call.nf`, `modules/themisto_coverage.nf`,
`bin/call_themisto_species.py`, `bin/aggregate_themisto_coverage.py`. Differences from
upstream are confined to paths (`results_dir`->`outdir`, publish under
`<outdir>/<sample>/sequenceindex/`, `bin/` not `rvi_toolbox/bin/`, viral-lens's
`workflows/`+`modules/` include layout) -- each ported file's header says so. Everything
else (`themisto2.nf`, `msweep.nf`, `cleanup.nf`, `reference_subset.nf`, `bowtie.nf`,
`samtools_coverage.nf`) was already in viral-lens and is reused unchanged.
`workflows/VIRAL_MSWEEP.nf` was **deleted** -- fully superseded.

Gotcha that cost a round of DAG checks: upstream's `themisto2-msweep.nf` starts with a
blank line before its shebang, and Groovy rejects a shebang on line 2 ("Shebang comment
should appear at the first line"). Strip leading blank lines when porting by `git show`.

### Method selection

`run_themisto` (**default true**) is the only method flag defaulting true.

| want | flags |
|---|---|
| Themisto only (default) | `--do_sequence_index true` |
| Themisto + mSWEEP estimate | `+ --run_msweep true` |
| Metagraph only | `+ --run_themisto false --run_metagraph_align/query true` |
| both | `+ --run_metagraph_align/query true` |

**`run_msweep` changed meaning.** It used to be the Themisto2/mSWEEP method's on/off
switch; it is now an add-on *inside* the Themisto arm ("also run mSWEEP's probabilistic
abundance estimate + its `MSWEEP_MAP_QC` validation"). Species calling no longer needs
mSWEEP at all. Benign for existing commands: `--run_msweep true` still yields the same
`msweep/`+`msweep_map/` output, and now additionally yields `themisto_hits/`+`themisto_map/`.

The two prefix-parameterised report helpers were renamed
`count_metagraph_species_hits`->`count_species_hits` and
`count_metagraph_map_qc`->`count_map_qc_breadth` (with `empty_*` to match): all three
read-hit methods emit a byte-identical `species_hits`/`map_qc` schema, so the metagraph
names were misleading. New report columns: `themisto_n_species_considered/called`,
`themisto_mapqc_n_species/max_breadth_pct`.

### Verification

Gating verified by `-with-dag` on all four combinations (the truncating progress display
is useless for this -- processes appear in the DAG file as `v<n>([NAME])` and subworkflow
scopes as `subgraph NAME`): default has `THEMISTO_PSEUDOALIGN`/`CALL_THEMISTO_SPECIES`/
`THEMISTO_MAP_QC` and no mSWEEP or Metagraph; `--run_msweep` adds `MSWEEP`+`MSWEEP_MAP_QC`;
`--run_themisto false` removes the arm entirely; both-on has both.

**Real farm run** (`nf_runs/themisto_real/`, `--do_sequence_index true`, single sample):
`Success: true`, 5m03s. MAPPING lane unchanged (same 8 taxid rows). Report:
`themisto_n_species_considered=20, called=10, mapqc_n_species=10,
mapqc_max_breadth_pct=99.9632`, every mSWEEP/Metagraph column at its `EMPTY_*` default.

**Genuine improvement over mSWEEP's reference choice:** on the same sample mSWEEP's
`MSWEEP_MAP_QC` gave SARS-CoV-2 **29.79%** breadth, because `SELECT_REFERENCE_RECORDS`
picks the *longest* sequence for the label. Themisto's hit-driven pick (the *most-hit*
record) gives **99.96%** breadth at 43.9x depth. That difference is the point of the
restructuring.

### Reference-index label quality limits the calls -- READ BEFORE TRUSTING SPECIES NAMES

The port is mechanically correct (SEQIDX->label positional alignment re-verified by hand
against `rvdb_clustered_virome_species_labels.txt`: line 464646 == "Betacoronavirus
pandemicum", 741451 == "Alphainfluenzavirus influenzae", and one line either side differs,
so there is no off-by-one). The *species names* are only as good as the index's labels,
and on this sample two real defects showed up:

- **7500 of 1321608 label lines (0.57%) are the literal string `NA`** -- now FILTERED,
  see below. Upstream aggregates them into one pseudo-species "NA", which on this sample
  drew 11021 hits and **89.79% breadth**, the second-highest call in the run. Its actual
  reference record (`SEQIDX_1320074`) is
  `PZ169745.1 Severe acute respiratory syndrome coronavirus 2 isolate .../2022`: real
  SARS-CoV-2 signal wearing a useless label, not a spurious hit -- which is exactly why a
  call named "NA" is worse than no call: it carries no information, and it is not a name
  MAPPING can ever match, so it read as a brand-new species.
- **At least one record is mislabelled in RVDB/GenBank itself.** `SEQIDX_741451`, labelled
  `Alphainfluenzavirus influenzae`, is
  `OZ387637.1 Influenza A virus (A/Michigan/45/2015(H1N1)) ... chromosome: MN908947.3` --
  29903 bp, and MN908947.3 *is the SARS-CoV-2 reference accession*. So the "influenza" call
  (6103 hits, 21.05% breadth) is SARS-CoV-2 reads hitting a SARS-CoV-2 genome deposited
  under an influenza name. **Themisto did not actually detect this sample's influenza** --
  do not report that as a win; Kraken2/MAPPING is what finds the real H1N1, as 7 segments.

The remaining 7 calls (Bat coronavirus, Sarbecovirus sp., RacCS203/224/264, RpYN06,
Horseshoe bat sarbecovirus) sit at 3-14% breadth: sarbecovirus cross-mapping, correctly
separated from the 99.96% true call by the breadth validation. That is `THEMISTO_MAP_QC`
doing its job.

### `NA` labels are now skipped (deliberate deviation from upstream)

`bin/call_themisto_species.py` maps placeholder labels (`UNUSABLE_LABELS`, currently just
`"NA"`) to `None`, so they are skipped exactly like a blank line and never reach
`species_hits`, `record_ids` or `index_label_map` -- and therefore never get map-QC'd,
reported, or offered to `--call_consensus_for_new_species`. This is the one behavioural
difference from the upstream branch; `modules/themisto_species_call.nf` carries a pointer
to it so it is not silently lost on a future re-port.

Verified by re-running the script directly against the real pseudoalignment files from
`nf_runs/themisto_real/work/`: the `NA` row disappears and **every other species' hit
count is unchanged** (`Betacoronavirus pandemicum` 14127, `Alphainfluenzavirus influenzae`
6103, ...), i.e. 19 considered / 9 called instead of 20 / 10 -- exactly one fewer of each.
The per-read dedup is unaffected, as it should be: a read hitting both an `NA` record and
a real species still counts once for the real species. Confirmed end to end on the farm
(`nf_runs/real_default/`, `nf_runs/real_msweep/`): zero `NA` rows in the published
`species_hits.tsv`, report reading 19 considered / 9 called,
`themisto_mapqc_max_breadth_pct` still 99.9632.

**The filter does NOT reach mSWEEP's own output.** `MSWEEP` is a separate binary that
reads `msweep_ref_groups` directly, so `<sample>_mSWEEP_abundances.txt` still carries an
`NA` group (0.000173 on this sample). That is below `msweep_map_min_abundance` (0.001), so
`count_msweep_abundances()` does not count it and `msweep_n_groups` is 1 -- but a sample
where the unlabelled records draw real signal could surface `NA` as an mSWEEP group. Only
the Themisto read-hit path is filtered.

**Remaining consequence for `--call_consensus_for_new_species` -- still unsettled.** With
`NA` gone, two of Themisto's calls still clear the 10.0 default
`new_species_min_breadth_pct` without matching anything MAPPING found, and would each get
a new-species consensus: `Bat coronavirus` (14.08% breadth) and `Betacoronavirus sp.
RpYN06` (11.06%). Both are sarbecovirus cross-mapping from the sample's real SARS-CoV-2,
not new species. Options, none yet chosen: raise `new_species_min_breadth_pct` above the
cross-mapping band (which sat at 3-14% here), or gate on mean_depth/reads_mapped rather
than breadth alone (these two had 0.46x and 0.40x mean depth against 43.9x for the true
call, so depth separates them far more cleanly than breadth does). **Don't enable both
features together until this is decided.**

Also still unproven: Themisto's hit calls have not been compared against mSWEEP's on the
same sample with `--run_msweep true` (which would produce both tables side by side), and
no multi-sample Themisto run has been done.

## mSWEEP is abundance-estimation ONLY -- its mapping/breadth validation was removed

`MSWEEP_MAP_QC` is gone from the flow and `workflows/MSWEEP_MAP_QC.nf` is **deleted**.
When `--run_msweep true`, `VIRAL_THEMISTO_MSWEEP` now runs `MSWEEP` and emits
`abundances` only.

**Why:** it answered the same question as `THEMISTO_MAP_QC` -- "does this species' call
hold up when reads are mapped against its reference?" -- and answered it worse.
`SELECT_REFERENCE_RECORDS` chose the *longest* sequence carrying a label; Themisto chooses
the species' *most-hit* record. Measured on the same sample and the same call
(`Betacoronavirus pandemicum`): **29.79% breadth via mSWEEP's pick, 99.96% at 43.9x depth
via Themisto's**. Keeping both cost a second bowtie2 index + mapping pass per sample to
produce the inferior answer.

Removed with it: the `mapqc_n_species` / `mapqc_max_breadth_pct` report columns,
`EMPTY_MAP_QC_COUNTS`, `count_msweep_map_qc()`, and mSWEEP's contribution to
`--call_consensus_for_new_species` candidates (it has no breadth table to threshold now --
Themisto and the Metagraph methods still contribute).

**Now dead but deliberately kept, not deleted** -- restorable if abundance-driven breadth
validation is ever wanted again, and each carries a comment saying so:
`SELECT_REFERENCE_RECORDS` (`modules/reference_subset.nf`), `AGGREGATE_SPECIES_COVERAGE`
and `GENERATE_MSWEEP_MAP_SUMMARY` (`modules/samtools_coverage.nf`). Their
`withName:` entries in `nextflow.config` remain and will emit the usual harmless "no
process matching config selector" warning. `msweep_map_bowtie_threads` is unused for the
same reason but stays defined so existing command lines still validate.
`msweep_map_min_abundance` (still used by `count_msweep_abundances()`) and
`msweep_map_reference_fasta` (still used by the new-species consensus path) are both in
real use -- don't remove those. **Note `INDEX_REFERENCE_FASTA`/`EXTRACT_REFERENCE_SUBSET`
and `SAMTOOLS_COVERAGE` in those same module files are very much alive** (Themisto map-QC,
Metagraph map-QC, and the new-species path all use them) -- don't prune the files wholesale.

**Verified:** `-with-dag` on four combinations shows `--run_msweep true` resolving `MSWEEP`
with none of `SELECT_REFERENCE_RECORDS`/`AGGREGATE_SPECIES_COVERAGE`/
`GENERATE_MSWEEP_MAP_SUMMARY`. Real farm runs `nf_runs/real_default/` (3m09s) and
`nf_runs/real_msweep/` (3m06s), both `Success: true`: no `msweep_map/` directory and no
`msweep_map_summary` anywhere, `msweep/` holding only the abundances + probs files, and
the report showing `msweep_n_groups=1`, `msweep_top_group=Betacoronavirus pandemicum`,
`msweep_top_abundance=0.999273` alongside Themisto's own 9 calls.

**First side-by-side of the two methods** (`nf_runs/real_msweep/`): mSWEEP reports one
group above threshold; Themisto calls 9 species, agreeing with mSWEEP on the top call and
additionally surfacing 7 sarbecovirus relatives at 3-14% breadth / 0.09-0.46x depth, which
the breadth/depth columns separate cleanly from the 99.96% / 43.9x true call. Themisto is
the more sensitive caller; the map-QC table is what makes that sensitivity usable.

## Classifier/consensus split: one shared MAPPING for both classifiers — DSL-checked only

Requested directly by the user, and a genuine wiring change rather than a move (unlike the
earlier `subworkflows/` reorganization). `MAPPING` was doing two separable jobs: Kraken2
taxid classification, then consensus + Nextclade + subtyping + classification report. Only
the first is Kraken2-specific, but the sequence-index lane needed the second — so it had
grown its own private `GENERATE_CONSENSUS` call and published a bare consensus with none
of the downstream treatment.

Now split at the consensus boundary:

| file | workflow | does |
|---|---|---|
| `subworkflows/classifying_kraken2.nf` | `CLASSIFYING_KRAKEN2` | `SORT_READS_BY_REF` (Kraken2 + Kraken2Ref), owns `identified_species_ch` |
| `subworkflows/classifying_index.nf` | `CLASSIFYING_INDEX` | the three sequence-index methods + their own per-method report; resolves references for new species |
| `subworkflows/mapping.nf` | `MAPPING` | `GENERATE_CONSENSUS` → Nextclade → SCOV2 → `GENERATE_CLASSIFICATION_REPORT`, for **both** classifiers |

**The interface between the three — the part to preserve if you touch any of them.** The two
classifiers hand over DIFFERENT shapes, deliberately (this changed once; see "Classifier
reports, MAPPING maps" below):

- `CLASSIFYING_KRAKEN2` arrives consensus-ready, because `SORT_READS_BY_REF` resolves its
  references as part of classifying:
  - `sample_taxid_ch` — `tuple(meta, [read_1, read_2], reference_fasta)`. `meta` carries
    `id` (`"<sample_id>.<taxid>"`), `sample_id`, `taxid`, `reference_header`.
  - `sample_report_with_join_key_ch` — `[join_key, report_meta]`, `join_key` equal to the
    matching consensus's `meta.id`. `report_meta` holds the descriptive per-(sample,
    reference) fields the classification report writes out.
  - `identified_species_ch` — `[sample_id, [normalized_species, ...]]`.
- `CLASSIFYING_INDEX` hands over calls only, no reads and no extracted reference:
  - `species_calls_ch` — `[sample_id, call]`, `call` being a Map of `species_name`,
    `reference_record`, `reference_source`, `hit_count`, `method`. (This shape changed
    again when map-QC was removed — see "Map-QC removed" below. It used to carry
    `breadth_pct` and no `reference_source`.)

`MAPPING` builds the consensus-ready shape for the index side itself, then `mix()`es.

`MAPPING` takes a fifth input, `CLASSIFYING_KRAKEN2.out.identified_species_ch`, and uses
it to drop sequence-index candidates for species Kraken2 already found *before* mixing —
so the two sets are disjoint because `MAPPING` made them so, not by assumption. See "Where
the species-level check lives" above for why that enforcement sits here.

`mix()` and not `join`/`combine` for the union itself, since post-filter these are disjoint
sets of (sample, reference) pairs rather than two views of the same pair. A classifier that
isn't running passes `Channel.empty()` — `main.nf` does this explicitly in an `else` branch
rather than reaching for `CLASSIFYING_INDEX.out.*`, which is undefined when the subworkflow
was never invoked (the same trap already documented for `VIRAL_METAGRAPH_ALIGN.out`).

**A trap this move walked straight into, worth knowing before editing these filters**: a
closure's parameter count must match the tuple's width. `combine(by: 0)` leaves the joined
right-hand element on the tuple, so the `.map` *after* a `.filter` needs a trailing
throwaway param (`_identified`) — omitting it aborts the run with "Invalid method
invocation `call` with arguments ... on _closureN". Caught by running the filter against
the real channel shapes in a standalone script, not by inspection.

**What actually changes behaviorally**: a species only Themisto2/Metagraph found now gets
Nextclade, SARS-CoV-2 subtyping and a classification-report row. Subtyping in particular
works because its synthetic `report_meta` carries `ref_selected` (set to the species name),
which is what `MAPPING` branches on — so a SARS-CoV-2 infection Kraken2 missed but a
sequence index caught gets subtyped like any other.

**Two smaller things that fell out of it**, both worth knowing before editing:

- `CLASSIFYING_INDEX`'s synthetic meta now mirrors `bin/k2r_report.py`'s pre-report columns,
  because `MAPPING` uses that map as the base of its report meta regardless of which
  classifier produced it. Fields this classifier has no honest equivalent for (`virus`,
  `num_reads`, `flu_segment`, `virus_subtype`, `sample_subtype`) are deliberately left
  empty — don't "fill them in" with sequence-index figures that mean something different.
- `sample_report_with_join_key_ch`'s rows now carry `id` as a field, not only as the tuple
  key. `MAPPING`'s Nextclade input uses the row as its meta base and downstream
  per-consensus JSON keying expects `meta.id`. The old code got this from a *second*,
  near-duplicate construction of the same channel off `sample_pre_report_ch`; that
  duplication is now collapsed into one channel used for both purposes.
- `new_species_consensus_n` is now `new_species_candidates_n`, counting candidates handed
  to `MAPPING` pre-filter — the only thing the classifier can honestly measure once
  `MAPPING` owns both the dedup and the consensus. The post-filter truth is in the
  classification report, via `discovered_by: 'sequence_index'` on the meta. Renaming was
  safe because the lane reports take their columns from the union of meta keys
  (`bin/write_lane_summary.py`), so there is no fixed schema to break.

**Verification**: `-preview` across four flag combinations (baseline; `--do_sequence_index`;
`+ --call_consensus_for_new_species`; `+ --do_assembly --do_abundance`), all resolving. The
centralization itself was confirmed empirically rather than by inspection: a `-with-dag`
export with the handover enabled has **one** `GENERATE_CONSENSUS` invocation where the
previous commit had two (`grep -c initial_alignment dag.mmd`: 2 → 1). **Nothing has been
run for real** — and note that the same caveat as the cross-classifier feature itself
applies: only a sample where Kraken2 and a sequence index genuinely disagree exercises the
handover at all, so a normal run proves only that the Kraken2 side still works.

## Map-QC removed from both classifiers; the breadth gate moved after consensus

Requested directly: *"i want to remove classifiers map-QC, only rely on hit counts. for
consensus calls index-only identified species after mapping, only do if breadth_pct > 10%
by default but parametrise"*.

**What map-QC was.** Each sequence-index method used to follow its species call with a
validation mapping: `INDEX_REFERENCE_FASTA`/`EXTRACT_*_REFERENCE_SUBSET` → `BOWTIE_INDEX` →
`BOWTIE2SAMTOOLS` → `SAMTOOLS_COVERAGE` → `AGGREGATE_*_COVERAGE`, producing a per-species
`breadth_pct`. `parse_species_calls()` then thresholded on that breadth to decide which
species were worth handing to `MAPPING`. A species that survived was therefore mapped
**twice** — bowtie2 for breadth, then `params.read_aligner` for consensus.

**What it is now.** Species are called on read-hit count alone
(`themisto_align_min_hits` / `metagraph_align_min_hits`). `parse_species_calls()` takes a
third argument and joins the two files the calling step already writes:

| file | schema | role |
|---|---|---|
| `<sample>_species_hits.tsv` | `sample_id, species, hit_count, provisional_call` | which species, how many hits; `provisional_call` is Python `str(bool)` → `"True"`/`"False"` |
| `<sample>_index_label_map.tsv` | headerless `<record_id>\t<species>` | the "ideal reference" per call — written only for species clearing min-hits |

so no reference resolution was lost with map-QC: the label map already named the record.
Species are matched between the two on name, whitespace- and case-normalized. A called
species with no record is skipped with a `log.warn` — legitimate for Metagraph, which drops
a species whose best-hit record another species already claimed
(`seen_record_ids` in `call_metagraph_species.py`).

`workflows/THEMISTO_MAP_QC.nf` and `workflows/METAGRAPH_MAP_QC.nf` are **kept but
unused**, with an `UNUSED` banner naming the two situations that would bring them back;
their params in `nextflow.config` are marked the same way (`themisto_align_run_map_qc`,
`metagraph_align_run_map_qc`, `*_map_bowtie_threads`).

**The breadth gate, and the one honest compromise in it.** Nothing maps these reads before
`GENERATE_CONSENSUS` any more, so the threshold has to be applied *after* it, in
`subworkflows/mapping.nf`. Below `params.new_species_min_breadth_pct` (default 10.0) a
sequence-index-**only** species is dropped completely — no published consensus, no
Nextclade, no subtyping, no report row — and a `log.info` line says so. Kraken2-side
consensuses pass through untouched, deliberately: they are gated by Kraken2's own read-count
selection, and retro-fitting breadth onto them would silently change what the pipeline has
always reported.

**Which breadth**, and why it is not the same number map-QC produced: the QC JSON's
`percent_non_n_bases`. It is genome breadth at iVar's minimum depth
(`ivar_polish_min_depth`, 10x) over the full reference length — valid because
`samtools mpileup -aa` (`modules/run_ivar.nf`) emits every reference position including
zero-coverage ones and `ivar consensus -n N` pads them with N, so the consensus is exactly
reference-length. True depth>=1 breadth is **not available**: `bin/qc.py` buckets depth in
steps of 5 from 0 (`range(0, 101, 5)`), so "positions with any coverage" is not among the
recorded thresholds. That is arguably the better gate anyway — the index noise this exists
to reject sits at 3-14% breadth / 0.09-0.46x mean depth (measured, see the Themisto2
section above), so a 10% threshold on depth>=1 breadth would admit the top of that range
while a 10x one cannot. `consensus_breadth_pct` is recorded on the meta of **every**
consensus, both classifiers', so the report shows the number the gate was applied to.

The cost of this ordering is that a noise call's consensus is computed and then thrown
away. That is the trade for not mapping every real call twice. If discarded consensuses
ever become the expensive part, restoring `THEMISTO_MAP_QC` is the fix — which is why it
is still in the tree.

**A pre-existing bug this exposed and fixed.** `MAPPING` extracted every sequence-index
reference from `params.msweep_map_reference_fasta` via `EXTRACT_REFERENCE_RECORD`
(`seqkit grep -p SEQIDX_<n>`). But the two index families report record ids in **different
namespaces**: Themisto2 gives positional `SEQIDX_<n>` into `msweep_map_reference_fasta`,
Metagraph gives a bare taxid or an accession into `metagraph_map_reference_fasta` (a
*different* file — `C-RVDBv32.0.renamed.fasta` vs `C-RVDBv32.0.fasta`). A Metagraph-called
species' grep therefore matched nothing, the optional output was absent, and the species
vanished silently between classifier and report. Fixed by carrying `reference_source`
(`'seqidx'` / `'metagraph'`) on the call, `branch`ing on it in `MAPPING`, and adding
`EXTRACT_METAGRAPH_REFERENCE_RECORD` (`modules/metagraph_reference_subset.nf`) for the
Metagraph side. The branch has a third `unknown_ch: true` arm that `error()`s — `branch`
silently drops what no arm matched, which is exactly the failure mode being fixed.

Two related details:

- The grep pattern for a Metagraph record is built in Groovy
  (`metagraph_record_pattern()` in `mapping.nf`), not in the process's shell: it is the
  same rule as `build_record_id_pattern()` in `bin/call_metagraph_species.py` (taxid →
  anchored on the fixed `kraken:taxid|<taxid>|` prefix; accession → anchored to header
  start and required to be followed by whitespace/EOL), and a regex surviving both
  Nextflow's interpolation and the shell's quoting is not worth the risk. Metacharacters
  are escaped individually rather than wrapped in `\Q...\E` — seqkit's engine is Go's,
  which does accept `\Q...\E`, but that cannot be exercised without seqkit installed,
  and an accession's `.` silently acting as a wildcard is a worse failure than none.
- Both `EXTRACT_*_REFERENCE_RECORD` processes now `rm` a zero-byte output. `>` creates the
  file whether or not the grep matched, which defeated `optional: true` and sent an empty
  FASTA into `GENERATE_CONSENSUS`.

**Verification** (all local, no farm access this round):

- `-preview` across five flag combinations — defaults; `--do_sequence_index`;
  `+ --call_consensus_for_new_species`; `+` all three methods `+ --run_msweep`;
  `+ --do_assembly --do_abundance` — all resolving.
- `-with-dag` process counts: **zero** `*_MAP_QC` nodes; with the new-species feature on,
  exactly one `INDEX_REFERENCE_FASTA` (was two — map-QC's own pass plus `MAPPING`'s) and
  both extractors present. With the feature off, `INDEX_REFERENCE_FASTA` no longer runs at
  all. `grep -c initial_alignment dag.mmd` still 1, so the centralization held.
- `parse_species_calls()` and `metagraph_record_pattern()` exercised against fixture TSVs
  in a standalone `.nf` harness: `provisional_call == 'False'` excluded; a called species
  absent from the label map excluded with the expected `log.warn`; `SEQIDX_41`/`SEQIDX_902`
  and taxid `2697049`/accession `OZ031634.1` all producing the right record and pattern
  (`^kraken:taxid\|2697049\|`, `^OZ031634\.1(\s|$)`); missing/empty inputs → `[]`.
- The breadth filter exercised on a stand-in channel: Kraken2-side row kept regardless,
  7.4% dropped, 96.2% kept, exactly-10.0% kept (`>=`, not `>`).
- **Not run for real.** `nf-test` could not run either — no JRE this agent can reach
  (`nextflow` works, `java -version` does not) — but no existing test covers `mapping.nf`
  or the parser, and `GENERATE_CONSENSUS` itself is untouched, so its snapshot is unaffected.

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
- ~~Do the `count_genomad_summary()` / `count_vrhyme_membership()` /
  `count_checkv_quality()` / `count_vcontact3_for_sample()` helpers parse the real TSV/CSV
  files correctly?~~ **ANSWERED — yes.** They live at the bottom of
  `subworkflows/assembly.nf` now (moved out of `main.nf`), and the 3-sample
  `nf_runs/assembly_multi_post/` run produced sane non-zero counts across all four,
  including the unbinned-scaffold case. See "Verification: DONE" under "`main.nf` split
  into `subworkflows/`" above for the actual numbers.
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
counts don't collide) are at the bottom of `subworkflows/classifying_index.nf` (moved out of
`main.nf` — see "main.nf split into subworkflows/" below), each
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
