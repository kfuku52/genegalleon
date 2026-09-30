# Batched verification, audited databases and PDF rendering

The single-family API, SQLite format, type promotion rules, plot arguments and
scheduler resources remain compatible. Content reuse lasts one operation only;
no saved hash or rendering receipt is accepted as scientific completion evidence.

## Verification

`workflow_api.py capabilities` advertises
`verify_batches=shared-source-verification-v1` and `verify_batch_limit=32`.
Use `verify --batch-file PLAN.json` instead of `--family-id`, keeping the existing
`--root`, `--workspace-root`, required steps and terminal profile. A plan is:

```json
{"schema":"genegalleon-verify-batch-v1","requests":[
  {"family_id":"OG0000001","attempt":"/workspace/output/observations/ATTEMPT",
   "recorded_workspace_root":"/workspace"}
]}
```

The `verify` envelope contains `batch_schema` and ordered `results`, each with
the unchanged single-family result schema. At most 32 distinct families are
allowed. Missing or stale artifacts stay family-local unverified results;
malformed evidence or a failed shared final content fence rejects the query.
Each exact attempt, receipt, family status and workspace mapping is rechecked.
Hashes of shared inputs are read once and rehashed before returning any results.
kfauto negotiates this capability, consumes each exact result in the current
32-record collector chunk and retains scalar queries for older runtimes.

## Audited database publication

The gene-summary core uses `gene_family_database_pipeline.py` to keep audit,
private DB build and provenance construction in one process. The pipeline accepts
repeated `--audit=TOKEN`, `--database=TOKEN`, and `--record=TOKEN` options carrying
the same argv tokens as the three standalone commands. It requires one matching
workspace/store and sole database output. Standalone commands remain available.

The audit checks all manifests, required inventories and CSUBST branch identity.
The DB remains private while its record is constructed from the private bytes at
the final declared output path. Unique source contents, collection membership,
manifests and archive generation are then revalidated before publication. The
existing namespace/manifest locks, atomic create/replace modes and exact attempt
record are retained. Audit progress shows `database_build` during SQL construction
and `source_revalidation` during the final fence; audit telemetry is never a
workflow completion proof.

Uniform large TSVs spool their original inferred pandas frames to private disk
and read them back without reparsing CSV. RAM stays chunk-bounded. Any dtype
promotion retains the original second CSV pass, preserving leading zeros,
nullable Boolean values, missingness, quoted tabs and multiline fields. Spools
are private temporary files and are removed on success, error or iterator close.
Spooling stops and cached chunks are discarded once a dtype change makes reuse
impossible. The extra disk writes are the main tradeoff; mixed-type files cannot
claim this parse saving.

## PDF worker

`Rscript workflow/support/tree_plot_batch.r PLAN.json` renders 1..32 jobs in one
R process, reusing loaded packages. The plan is a JSON array:

```json
[{"id":"OG0000001","cwd":"/workspace/scratch","output":"/workspace/result.pdf",
  "args":["--stat_branch=/workspace/stat.branch.tsv",
          "--max_delta_intron_present=-0.5","--panel_widths_mm=tree:60",
          "--panel1=tree,bl_rooted,no,no,L","--show_branch_id=no",
          "--event_method=species_overlap","--species_color_table=PLACEHOLDER",
          "--pie_chart_value_transformation=identity","--long_branch_display=no"]}]
```

Supply absolute input file paths. Each job uses a fresh environment, input cache
and scratch directory. Graphics devices, options and the optional species parser
are reset between jobs; garbage collection bounds retained family data. Input
file content/signatures (including comma-containing filenames) and the exact
cached renderer source bytes are fenced
before atomic PDF publication. An individual render failure preserves its prior
output and does not contaminate later jobs.
IDs and output paths must be unique. `PLAN.json.results.json` reports ordered
per-job exit codes, including 42 for unavailable optional ggimage, with
`completion_evidence=false`. The process exits nonzero if any job fails.

A caller must still record and verify each family's provenance and exact attempt
using the workflow's normal stage functions. This renderer does not submit jobs
or change single-family array scheduling; one-family jobs retain the existing R
CLI, and their startup is not amortized by this batch interface.

## Comparable synthetic measurements

Linux arm64 Docker `local/genegalleon:dev`, Python 3.12.14; same container and
inputs, fresh processes, one warmup plus three trials for API/PDF/pipeline.
Large-file DB measurements use three trials, after prior workload warmup.
The comparison base is `3aecee4`. Times exclude fixture generation. No scheduler
workload or production completion-rate claim is involved.

| Workload | Before median | After median | Ratio |
| --- | ---: | ---: | ---: |
| 16 families, five declarations each, shared 64 MiB source | 0.971 s | 0.125 s | 7.77× |
| Eight tree PDFs, separate R processes versus one worker | 9.188 s | 2.796 s | 3.29× |
| 128 audited families, 32 MiB alignments, complete DB/record pipeline | 1.164 s | 1.019 s | 1.14× |
| One 131,072-row TSV, 24 numeric metrics, complete DB build | 1.413 s | 1.390 s | 1.02× |

PDF and DB timings were rerun after adding the renderer source fence and early
spool disposal. The many-small-file DB control changed 2.128 to 2.236 s; an earlier
run changed 2.153 to 2.129 s. Large-file DB gains also varied from 1.10× to 1.02×,
so removing the second CSV parse does not establish a reliable overall DB gain.
API parent peak RSS was about 108 MiB before/after. PDF worker peak RSS increased
about 11 MiB (229 to 240 MiB); a warmed R namespace remains resident. Uniform
large-file DB peak RSS stayed about 134 MiB. The combined audited DB
child retained about 13 MiB more (118 to 131 MiB) because both audit and DB
modules remain loaded; the parent fixture process used less memory.
Verification decisions/contracts, sorted SQLite schema/rows, recorded input
fingerprints and audit rows match. PDF drawing bytes match after removing only
CreationDate/ModDate. SQLite binary hashes can differ because parallel insertion
order varies; each generated provenance record must match its own DB bytes.

Run comparable revisions through the existing runtime wrapper, without parallel
benchmarks or tests:

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_remaining_io.py \
  --case verify-batch --output /tmp/verify.json
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_remaining_io.py \
  --case pdf-batch --output /tmp/pdf.json
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_summary_pipeline.py \
  --output /tmp/pipeline.json
bash workflow/tests/run_in_runtime.sh python workflow/tests/benchmark_orthogroup_database.py \
  --output /tmp/database.json
```

The first three accept `--support-root` for a complete saved baseline support
tree; the DB benchmark accepts `--source` for the baseline generator. Compare
logical outputs and recorded inputs, validate each record's actual output hash,
and report RSS/disk tradeoffs alongside timing.
