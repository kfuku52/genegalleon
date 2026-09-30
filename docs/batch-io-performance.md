# Batched verification, audited databases and PDF rendering

The single-family API, SQLite format, type promotion rules, plot arguments and
scheduler resources remain compatible. Content reuse lasts one operation only;
no saved hash or rendering receipt is accepted as scientific completion evidence.

## Indexed query ownership

Bulk query2family summary, storage conversion/materialization and provenance
audit now build one operation-scoped query-ID priority index. Only complete
filenames and underscore/dot boundaries are looked up; each file no longer scans
every query ID. The scalar matcher and catalog-list API remain available.
Overlapping IDs, arbitrary matcher priority, duplicates, Unicode and empty IDs
retain scalar behavior. The index snapshots its catalog and caches no filenames.

Linux arm64 Docker `local/genegalleon:dev`, Python 3.12.14, base `40b78b4`, one
warmup and three measured fresh processes: 3,003 IDs / 9,006 filename matches
took 1.929 s to 0.0102 s (190×); a complete three-directory, 3,003-row summary
took 1.502 s to 0.0619 s (24.2×). Result lists and complete TSV bytes match.
Peak RSS increased about 0.4 MiB (about 75 MiB matching / 77 MiB summary).
These are query-ownership workloads, not overall scientific-pipeline speedups.

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_query_matching.py \
  --output /tmp/query-matching.json
```

Use `--support-root` with a complete baseline support tree for comparison.

## GFF transcript statistics

Longest-transcript selection collects row positions and builds one result table,
preserving gene discovery order, annotation row order, tie handling and extra
column types. CDS phase validation indexes coordinates once per gene, retaining
every phase record at duplicate coordinates. Coordinate ordering, UTR validation,
strict/report policies, missing phases and trans-splicing behavior are unchanged.

Linux arm64 Docker `local/genegalleon:dev`, Python 3.12.14, pandas 3.0.6, base
`f6f68bf`, one warmup and three measured fresh processes: the live in-memory
selection → structure validation → statistics path for 2,000 genes took 2.607
to 0.686 s (3.80×). Selection took 0.495 to 0.0591 s; structure validation took
2.072 to 0.591 s. Median peak RSS was 129.5 to 98.8 MiB (23.7% less). The
fixture includes both strands, alternative isoforms, explicit UTRs and duplicate
CDS records. Intermediate/final tables, column types/order and warnings match.
Time excludes fixture construction and result fingerprinting; process peak RSS
includes both. This does not measure GFF download, file parsing or sequence
resolution.

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_gff_summary.py \
  --output /tmp/gff-summary.json
```

Use `--support-root` with a complete baseline support tree for comparison.

## Wide-header validation

Database and scan-schema preflight count each column name once, preserving
sorted duplicate diagnostics and rejection before DB replacement. Duplicate
checking no longer scans the full header for every column. The header reader
now follows standard TSV quoting and UTF-8 BOM handling, matching the body
reader. Quoted tabs/newlines, quotes and literal spaces in names are preserved;
duplicates remain visible and fatal instead of being renamed by pandas.

Linux arm64 Docker `local/genegalleon:dev`, Python 3.12.14, base `f6f68bf`, one
warmup and three measured fresh processes: validating 128 scan headers with
2,048 columns took 3.780 to 0.0466 s (81.1×) for raw files and 3.862 to 0.1304 s
(29.6×) for ZIP-held files. Peak RSS stayed about 105 MiB raw / 107 MiB ZIP.
The fixture includes one deliberate duplicate-column error; the complete sorted
diagnostic matches. These measure strict schema preflight, not full DB creation.

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_schema_validation.py \
  --output /tmp/schema-validation.json
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_schema_validation.py \
  --storage zip --output /tmp/schema-validation-zip.json
```

Use `--support-root` with a complete baseline support tree for comparison.

## Declaration, trait and species-tree validation

Provenance declaration keys, trait headers/selections and copy-number species
labels are counted once instead of rescanning each list for every value. Trait
selection also builds one available-name set. Sorted duplicate diagnostics,
selection order, error precedence and rejection before artifact reads remain
unchanged; file hashing and content/completion checks are retained.

Linux arm64 Docker, Python 3.12.14, pandas 3.0.6, base `3013667`, same immutable
image, one warmup and three measured fresh processes:

| Synthetic validation workload | Before median | After median | Ratio |
| --- | ---: | ---: | ---: |
| 4,096 input declarations, full contract plus cross-kind duplicate failure | 0.578 s | 0.132 s | 4.38× |
| Four-row / 4,096-column trait TSV, read and valid/duplicate selections | 0.404 s | 0.0292 s | 13.8× |
| 4,096-leaf species tree, valid and duplicate-label validation | 0.236 s | 0.0168 s | 14.1× |

Complete contract, trait TSV, selections, Newick/leaf order and diagnostics
match. Combined median process peak RSS was about 160.6 MiB before/after.
The declarations reference one small shared source, so this measures declaration
scaling rather than bulk content hashing. Fixture construction and fingerprints
are excluded from timing and included in peak RSS. These large validation
fixtures do not establish full analysis speedups or gains for small inputs.

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_declaration_validation.py \
  --output declaration-validation.json
```

Use `--support-root` with a complete baseline support tree and keep the same
`GG_CONTAINER_DOCKER_IMAGE` immutable image ID for both runs.

## Alignment-statistics I/O

Both summary readers keep at most eight pending reads per configured worker,
consume completion notifications and discard completed futures. They retain
their original worker count, named-column mapping and failure behavior. Inventory
and output tables still require memory proportional to the number of families;
only queued/completed task retention is bounded. Queued work is cancelled on
failure or iterator close, and active I/O finishes before the executor exits.

Linux arm64 Docker `local/genegalleon:dev`, Python 3.12.14, base `e56c314`, four
workers, one warmup and three measured fresh processes: 10,000 raw statistics
TSVs took 7.01 s to 6.30 s (1.11×), with median peak RSS 110.2 to 98.7 MiB
(10.5% less). A 1,000-family query2family ZIP fixture took 1.05 to 1.04 s with
about 99 MiB RSS on both versions; no speedup is established for that smaller
workload. Complete result TSV bytes match for both fixtures. Fixture generation
is excluded from time, but included in process peak RSS.

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_alignment_summary.py \
  --output /tmp/alignment-summary.json
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_alignment_summary.py \
  --reader query2family --storage zip --families 1000 --output /tmp/alignment-summary-zip.json
```

Use `--support-root` with a complete baseline support tree for comparison.

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
