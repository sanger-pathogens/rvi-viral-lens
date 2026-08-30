# Pipeline route maps with nf-metro

This directory holds `route_map.mmd` (source) and `route_map.svg` (rendered) — a
[nf-metro](https://github.com/seqeralabs/nf-metro) sketch of `main.nf` as envisioned after
the `rvi_integration_1` merge of `rvi-viral-metagenomics-pipeline` functionality into
viral-lens. Nothing here is generated from the actual Nextflow code yet — it's a planning
diagram, hand-authored to scope the integration before writing it.

## Installing nf-metro

It's a pip package, not part of either pipeline's dependencies. Isolate it with `pipx` (or
a venv) rather than installing into system Python:

```bash
pipx install nf-metro
nf-metro --version
```

Requires Python 3.10+. `pip install nf-metro` / `conda install bioconda::nf-metro` also
work if you already manage an environment.

## Authoring a map

A map is a Mermaid `graph LR` file with `%%metro` directives layered on top — plain Mermaid
comment syntax, so any Mermaid tooling that doesn't understand nf-metro just ignores them
and renders the flowchart minus the metro-line styling.

Minimum you need:

```
%%metro line: <id> | <legend label> | <hex color>

graph LR
    subgraph <section_id> [<Section title>]
        <station_id>[<Station label>]
        <station_a> -->|<line_id>,<line_id2>| <station_b>
    end
```

- **Lines** (`%%metro line:`) are colored routes through the pipeline — either genuinely
  alternative/optional paths (an aligner choice, an opt-in step) or just a readable way to
  color-code parallel branches that all run for every sample. Every edge carries a
  comma-separated list of which line(s) traverse it: `a -->|core,taxid| b`.
- **Sections** (`subgraph`) are the pipeline's logical stages, laid out left to right.
  `%%metro entry:`/`exit:` directives (`left`/`right`) hint which lines cross a section
  boundary, needed once a section has more than one inbound/outbound line.
- **File stations** (`%%metro file: id | TYPE | label`) get a file-type icon instead of a
  process box; `%%metro off_track: id, id, ...` pulls them out of the main trunk so
  inputs/outputs don't crowd the flow.
- **Fan-out/fan-in**: one line can visit several stations in parallel just by giving that
  line's id on multiple outgoing edges from the same station (see `preprocessing`'s exit
  into five parallel sections in `route_map.mmd`) — no special syntax needed.
- **Process directives** (`%%metro process: station_id | regex`) map a station to the
  Nextflow process name(s) it represents, only relevant for the live-progress feature below.

Validate before rendering — the layout engine is strict about certain topologies
(a many-way convergence of lines onto one section can trip a real bug in the router; if
`render` throws a `CurveInvariantError`, restructure the convergence rather than fighting
flags):

```bash
nf-metro validate route_map.mmd
```

## Rendering

```bash
nf-metro render route_map.mmd -o route_map.svg --format svg --embed-font --validate
```

- `--embed-font` inlines the type as base64 so the SVG looks the same on any machine.
- `--validate` re-checks the rendered geometry (not just the source) and fails loudly
  instead of silently shipping a broken layout.
- `--format html` instead produces a self-contained interactive page (pan/zoom, click a
  line in the legend to isolate it) — useful while iterating, not needed for a static
  README embed.
- Useful layout escape hatches when a map gets dense: `--diamond-style symmetric`,
  `--line-spread rails`, `--section-x-gap`/`--section-y-gap`, `--fold-threshold`.

## Live progress overlay

Once `%%metro process:` directives are in place for a map's stations (see `route_map.mmd`
for a first pass — station-to-process names are best-effort until `main.nf` actually
exists and the real process names can be checked against it), the same map can light up
stations in real time as an actual pipeline run progresses:

```bash
# serve the map, wiring up automatically and stopping when the run finishes
nf-metro serve route_map.mmd --open --shutdown-after-complete -- \
  nextflow run . -profile docker

# or point an independently-launched run's weblog at a map already being served
nf-metro serve route_map.mmd --port 8080
nextflow run . -with-weblog http://localhost:8080/events
```

Before relying on it, check the mapping actually covers the pipeline — this reports
stations whose regex matches nothing real (stale) and real processes the map can't show
(drift), and is worth wiring into CI once `main.nf` exists so a renamed process doesn't
silently stop lighting up:

```bash
nf-metro check-mapping route_map.mmd --dag with_dag_export.mmd
# or: --processes <(nextflow log <run> -f process | sort -u)
```

## What this particular map encodes

Five lines: `core` (every sample, through preprocessing and reporting), `taxid` (the
existing map-to-taxid → consensus path), and three proposed additions lifted from
`rvi-viral-metagenomics-pipeline` — `mapping`, `assembly`, `abundance`. Stations whose
label ends in "`- decision <letter>`" mark a spot the map can show the shape of but can't
resolve on its own; those are tracked separately, not duplicated here since they'll drift
out of date the moment they're decided.
