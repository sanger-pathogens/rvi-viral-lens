# rvi_integration_1 — handoff to continue on HPC

You're picking up mid-flight on branch `rvi_integration_1` of `viral-lens`. This branch is
integrating functionality from a sibling pipeline, `rvi-viral-metagenomics-pipeline`, into
`viral-lens/main.nf`. Everything below was true as of commit `d759968` (`HEAD` when this
file was written) — `git log` is the source of truth if it's since moved.

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

### 3. Build the "map reads to sequence indexes" lane — Themisto2 DONE, Metagraph outstanding

Per the route map (`docs/nf-metro/route_map.mmd`, section `seq_index_mapping`): three
methods converging on one QC Mapping step, feeding `GENERATE_MAPPING_REPORT.nf` — which
is now actually wired and running, no longer an orphan.

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
  `25.GB * task.attempt` at the user's request — **apply the same to Metagraph when you
  port it.** A broader request-vs-actual audit across every module is wanted later.
- **Sequence-to-graph alignment via Metagraph**: `VIRAL_METAGRAPH_ALIGN`
  (`subworkflows/metagraph_align.nf` on `eu1/rvi_toolbox.git`, see above) +
  `METAGRAPH_MAP_QC` (`subworkflows/metagraph_map_qc.nf`). Depends on resolving item 2
  first (this subworkflow doesn't exist on the fork viral-lens currently tracks).
- **Pseudoalign via Metagraph**: flagged in the route map as *not currently a real,
  maintained module* — `modules/metagraph.nf`'s `METAGRAPH` process used to do a plain
  k-mer query (pseudoalign) before commit `ccae756` on `eu1/rvi_toolbox.git` ("Switch
  METAGRAPH to true alignment instead of plain k-mer query") replaced it with true
  alignment. Building this as a genuine third option means reviving that mode as its own
  module (e.g. a `metagraph query` invocation), not just re-exposing something already
  there. Confirm this is still wanted before building it — it's the one piece of the route
  map that was never actually implemented anywhere.
- **QC Mapping**: reuse `MSWEEP_MAP_QC` (msweep path) / `METAGRAPH_MAP_QC` (metagraph
  paths) as-is; don't invent a new unified QC step unless there's a concrete reason to.

`GENERATE_MAPPING_REPORT` is already wired (see `main.nf`'s
`if (params.do_sequence_index && params.run_msweep) { ... }` block). When you add a
Metagraph method, feed its per-sample counts into the same `mapping_report_prep_ch`
rather than building a second report path. The two helpers that populate that meta,
`count_msweep_abundances()` / `count_msweep_map_qc()`, are at the bottom of `main.nf`
beside the assembly ones.

Note the `remainder: true` on that block's `join`: `MSWEEP_MAP_QC` emits optional
outputs, so a sample with nothing above `msweep_map_min_abundance` produces no QC table
and would otherwise vanish from the report entirely. `EMPTY_MAP_QC_COUNTS` fills it in.

### 4. Build the abundance estimation lane

Per the route map, section `abundance`: `Kraken2 + Bracken` (`KRAKEN2BRACKEN` subworkflow
— already confirmed to emit `abundance_summary`, see its subworkflow file) and
`ABUNDANCE_ESTIMATION` run in parallel; `SCRuB` sits as `KRAKEN2BRACKEN`'s downstream final
step, **not** a trunk-level preprocessing step — it consumes the whole-run Bracken summary
plus a user-supplied plate map (`params.scrub_plate_map`, mandatory when `params.run_scrub`
is true) and re-estimates relative abundance. `ABUNDANCE_ESTIMATION` doesn't go through
SCRuB. Depends on item 2 (SCRuB lives only on `eu1/rvi_toolbox.git`).

Wire `GENERATE_ABUNDANCE_REPORT` in the same pattern once this lane exists.

### 5. Widen input handling

Today's `main.nf` still only has `parse_mnf()` (single local reads manifest). The route
map's `Preprocessing` section already models three sources feeding one `MIXED_INPUT`
merge: local reads manifest, ENA download (`ENA_DOWNLOAD`), and iRODS retrieval
(`PARSE_IRODS_INPUT` → `DOWNLOAD_FROM_IRODS`, by CLI params or a lane manifest) — all in
`rvi-viral-metagenomics-pipeline/rvi_toolbox/subworkflows/mixed_input.nf` (+
`combined_input.nf`, `ena_input.nf`, `irods.nf`, `irods_manifest_parse.nf`,
`input_check.nf`). Read `mixed_input_README.md` alongside it first — it documents which
params (`--manifest_ena`, `--manifest_of_lanes`, `--studyid`/`--runid`/etc.) activate which
source. `main.nf`'s preprocessing block already has a `preprocessed_3tuple_ch` designed to
be the single downstream shape regardless of how reads arrived — feed `MIXED_INPUT`'s
output into that same shape rather than restructuring what's downstream of it.

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

- **New subworkflows/modules land in `viral-lens/workflows/`, `modules/`, `bin/` directly**,
  not inside the `rvi_toolbox` submodule, until the fork situation (item 2) is resolved.
  Module `include` paths stay relative (`../modules/x.nf` from `workflows/`) since both
  repos use the same `workflows/` + `modules/` sibling layout.
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
- Don't mark the mapping lane "done" without either resolving or explicitly flagging the
  pseudoalign-via-Metagraph gap (it doesn't exist as a real module anywhere yet, see item 3).
