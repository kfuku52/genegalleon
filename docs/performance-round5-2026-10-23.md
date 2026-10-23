# Repository performance review: data lifetime and archive reuse

This review starts from `76525a358c1bd393f53ac54abe6bb72c6ff5362a`
(0.8.388), after [round four](performance-round4-2026-10-23.md). It examines
data lifetime, redundant parsing and worker preparation across the repository,
rather than further changing numerical reductions. Three agents inspected
normal callers and existing performance reports. Baselines preceded production
edits; measurements and runtime checks were serialized.

## Candidates

| Idea | Decision / verification boundary |
| --- | --- |
| Reuse ZIP readers across alternating logical subdirectories outside snapshots | Accept: a bounded per-thread pool avoids repeatedly parsing central directories in the normal ortholog summary. |
| Release evicted compaction readers' central-directory objects | Accept: remove strong references retained by per-open cleanup callbacks while preserving bounded descriptors and complete publication checks. |
| Retain only FASTA identifiers in GFF modes that consume no sequence values | Accept: lower peak RAM in the normal whole-species GFF CLI and cold synteny cache build. |
| Read FASTA metadata without constructing discarded sequences | Accept: stream header metadata for provider planning and reused-output counting, preserving record consumption boundaries. |
| Batch ordered GRAMPA process tasks | Deferred: IPC savings versus heterogeneous tree workloads need a separate process-tree and first-error comparison. No production change or speed claim. |
| Coordinate CSUBST worker counts and per-tool thread budgets | Deferred: nested CPU budgets warrant real-tool qualification and scientific output comparison before changing scheduling. No production change or speed claim. |
| Remove a second species-store status size-stat pass | Deferred: file changes/deletions between current passes need explicit error/snapshot semantics before reusing observations. |
| Reuse command-scoped parsed worker plans/settings | Deferred: repeated JSON parsing may be reducible, but fresh worker-proof boundaries must remain explicit. |
| Reuse CSUBST scan headers between database preflight and loading | Deferred: the snapshot scope and malformed-input/mutation checks need their own complete database comparison. |

## Environment and method

The normal runtime wrapper verified the current container inputs and daily owned
upstream snapshot. Comparisons use Docker image
`local/genegalleon:performance-round3-20261023`, immutable ID
`sha256:0e46be37b7ddbdb60b331ef35017abd3f2d90820d152ce3838ec369537791752`,
runtime-input fingerprint
`8f5655e6d387393268bb292d495b28222ad68b0e072a8a622722f6daef30faa4`.
The runtime is Linux aarch64 with Python 3.12.15, pandas 3.0.6 and NumPy 1.26.4.
Docker Desktop supplies 12 CPUs and 8,214,585,344 bytes of RAM, without an
additional per-container limit.

Final comparisons use five alternating fresh-worker pairs after a discarded
warmup pair. Tables show medians and paired wall-speedup ranges. Fixtures are
prepared outside measured workers; setup/copying and verification are outside
timers. Linux `VmHWM` is captured before verification, including worker imports
and setup allocations. Parent and child measurements are distinguished where
commands spawn processes. Raw results, source/input/output hashes and logs remain
under `tmp/performance-20261023-round5/`. Synthetic warm-file-cache evidence does
not establish SIF, NAS, or whole scientific workflow performance.

## ZIP reader reuse

The non-snapshot reader previously retained one ZIP reader per thread. The
normal `query_gene_orthologs` summary alternates `stat_branch` and `cds_fasta`
members per family, repeatedly opening the same two central directories. The
snapshot pool was already optimized; this change concerns the separate normal
reader path.

The new pool retains at most eight readers per thread, evicts idle readers and
uses an uncached reader when all entries are active. Path device/inode/size/
mtime/ctime signatures are checked against the open descriptor's identity.
Index epoch changes reset every reader. Live precedence, observations, producer
locks, symlink checks and snapshot behavior remain in place. This trades up to
seven additional retained descriptors per thread and their ZIP metadata for
less repeated parsing; no general memory improvement is claimed for this path.

The benchmark runs the real ortholog summary, with fixture generation and target
imports outside its timer. Every output TSV is compared as complete bytes.

| Families / storage | Wall before → after | Paired speedup range | CPU before → after | Peak RSS before → after |
| --- | ---: | ---: | ---: | ---: |
| 8 / ZIP | 0.40589 → 0.40302 s | 0.946–1.050× | 0.40064 → 0.39578 s | 110.46 → 110.47 MiB |
| 256 / raw | 1.65876 → 1.66143 s | 0.989–1.037× | 1.65025 → 1.64989 s | 133.78 → 133.79 MiB |
| 256 / ZIP | 1.61779 → 1.39415 s | 1.13–1.20× | 1.60986 → 1.38630 s | 133.37 → 133.38 MiB |
| 1,024 / ZIP | 8.16831 → 4.48386 s | 1.77–1.83× | 8.05697 → 4.45216 s | 205.02 → 204.64 MiB |

All medium/large ZIP pairs improve. Raw and small controls cross unity; no
performance change is claimed there. A separate instrumented screen observes
512 → 2 source ZIP opens for 256 families. Final timing/RSS comparisons disable
that observer, avoiding an asymmetric instrumentation cost.

## Compaction reader lifetime

Compaction already limited open source readers to 32. However, `ExitStack`
registered a bound exit method for every opened reader, retaining every evicted
`ZipFile` and its central directory until the output chunk finished. Closing a
descriptor did not release that metadata. Cleanup now registers one callback
over the active pool; evicted readers are closed and dropped. Final cleanup uses
an exit stack over the remaining bounded set so one close failure still closes
the others. SHA256/CRC verification, compression, manifests, indexes and atomic
publication are unchanged.

Fixtures interleave family identifiers across source shards, causing the real
default 32-reader cache to cycle. The real compact operation includes reads,
compression, complete verification, index publication and source retirement.

| Source shards × members | Wall before → after | Paired speedup range | CPU before → after | Peak RSS before → after |
| --- | ---: | ---: | ---: | ---: |
| 8 × 16 = 128 | 0.07959 → 0.07475 s | 0.323–1.936× | 0.04255 → 0.04349 s | 40.84 → 40.87 MiB |
| 64 × 32 = 2,048 | 0.68239 → 0.63866 s | 1.05–1.13× | 0.60460 → 0.56758 s | 78.52 → 43.63 MiB |
| 64 × 64 = 4,096 | 1.67396 → 1.41724 s | 1.16–1.22× | 1.60019 → 1.31479 s | 188.05 → 49.09 MiB |

The largest case reduces worker peak RSS by about 74%, with all five paired
timings improving. A separate structural measurement finds 2,048 → 33 live
reader objects while the descriptor peak remains 37; the extra object is a
recently evicted local reference. Weak-reference tests also verify collection
on exit. The small control includes a retained 0.221 s candidate outlier, and
has no claimed speed or memory improvement.

Final comparisons have no tracking wrapper. The benchmark fixes only wall-clock
`time.time`/`gmtime` during compaction to make creation timestamps deterministic;
performance/process/monotonic clocks and filesystem timestamps remain real.
All actual ZIP bytes, ZipInfo compression/metadata, manifest and index/catalog
records match, and deep verification passes. The opaque random index epoch is
validated for format and replacement rather than nonce-byte equality. Across
both ZIP changes, all 84 final workers match outputs. Independent review checked
pool invalidation, active handles, cleanup and benchmark accounting.

Final ZIP data are `zip/confirmation-uninstrumented.json` and `zip/receipt.json`.

## GFF input payload lifetime

The normal cold synteny species-cache builder invokes `gff2genestat` without
sequence-consuming options. The GFF CLI previously retained every full FASTA
record throughout GFF processing. It now retains that original record list only
when `--validate-cds-length`, `--cds-resolution-dir` or `--sequence-store` needs
sequence values. The identifier-only path fully drains the same parser into the
same string-typed identifier Series. Duplicate checks still follow the complete
FASTA parse, before GFF discovery or output writes. No sequence parser,
scientific validation, structure policy or source binding is changed.

These fixtures contain 30,000 genes with 1,500nt CDS and their GFF features.
The real CLI includes target imports and main. Cold-cache measurements call
`ensure_species_gene_cache`, including its real nested GFF subprocess. Fresh
workers have no earlier children: CPU reports parent plus child CPU, while RSS
reports the main CLI process or the largest child separately. These maxima are
not added or presented as a simultaneous process-tree peak. Warm-cache setup
runs in a separate preparation process and the measured warm call starts no
GFF reader subprocess.

| Scope / transport | Wall before → after | Paired speedup range | Total CPU before → after | Main or child peak RSS before → after |
| --- | ---: | ---: | ---: | ---: |
| CLI / raw | 5.62342 → 5.48517 s | 0.993–1.064× | 5.61768 → 5.47938 s | 238.50 → 188.65 MiB |
| CLI / gzip | 5.15305 → 5.13221 s | 0.990–1.015× | 5.14826 → 5.12738 s | 238.24 → 188.93 MiB |
| Cold synteny cache / raw | 5.50661 → 5.46324 s | 0.999–1.019× | 5.49944 → 5.45633 s | child 241.22 → 191.62 MiB |
| Cold synteny cache / gzip | 5.62045 → 5.42133 s | 1.006–1.058× | 5.61261 → 5.41468 s | child 241.04 → 191.21 MiB |

The normal main/child RSS reductions are consistently 20.6–20.9%, about
49–50 MiB. The cold caller parent remains about 77.3 MiB. Adoption is based on
RAM, with no general CPU or runtime improvement claim. The 32-gene controls and
warm-cache calls have no claimed RAM improvement. Warm calls run no GFF child
and remain about 77.4 MiB; raw wall medians are 0.20455 → 0.20463 s and gzip
0.17891 → 0.17565 s, with both paired ranges crossing unity.

Across three inputs and three scopes, all 108 workers (45 measured pairs and
nine discarded warmup pairs) match complete TSV bytes, schema, typed frames and
logs. Cache receipts assert the correct actual implementation SHA before
normalizing that one expected SHA field; every other contract field and output
hash is compared. Log normalization is limited to temporary paths, datetime and
elapsed-time text. Additional differential checks match 34 complete successful
results and 38 exact errors, including unchanged previous outputs. The 19 new
tests cover all 12 consuming-flag/policy combinations and parser failure
priority; existing source/resolution and real synteny/DIAMOND consumers also
pass. An independent agent reviewed every record consumer and the measurement
accounting.

Final records are `gff/confirmation.json`, `gff/contracts-result.json` and
`gff/focused-tests.log`.

## FASTA metadata without sequence assembly

Provider discovery checks whether FASTA headers overlap GFF sequence names, and
resume discovery counts records and reads the first identifier. Both previously
assembled and whitespace-normalized complete sequence values that they discarded.
A private header scanner now serves these two metadata consumers. It yields the
previous header only after reading the next header or EOF, matching the original
record boundary and decoder/error consumption. Raw/gzip/tar opening, archive
member selection/order and closing use the same shared opener. The full-sequence
parser body is unchanged, and its public API still yields complete records.

Fixtures use seeded DNA wrapped at 80 bases: 32 small records, 16/64 MiB
chromosomes, four 16 MiB contigs with a last-contig match, and 20,000 CDS records
of 1,500nt. Raw/gzip and a 16 MiB tar.gz case are covered. Controllers prepare
files before workers; input hashes are verified before and after comparisons.
Canonical package loading and callable source-filename assertions ensure that
each worker actually executes its selected baseline or candidate.

API timers exclude imports. The planner measurement runs the real local
`plan_input_generation_tasks.main`, including fresh input hashing and JSON
publication, with imports excluded. These planner timings do not include cold
interpreter/dependency startup. Worker RSS includes imported dependencies.

| Scope / input | Wall before → after | Paired speedup range | CPU before → after | Peak RSS before → after |
| --- | ---: | ---: | ---: | ---: |
| Discovery / 64 MiB raw | 0.77698 → 0.08344 s | 8.83–9.82× | 0.77541 → 0.08234 s | 224.91 → 52.34 MiB |
| Discovery / 64 MiB gzip | 0.97441 → 0.29704 s | 3.19–3.45× | 0.97335 → 0.29608 s | 224.90 → 52.35 MiB |
| Count / 64 MiB raw | 0.76110 → 0.07775 s | 9.42–10.00× | 0.76028 → 0.07698 s | 225.88 → 51.46 MiB |
| Count / 64 MiB gzip | 0.99677 → 0.29393 s | 3.24–3.71× | 0.99633 → 0.29339 s | 225.82 → 51.46 MiB |
| Discovery / four 16 MiB contigs, gzip | 0.98583 → 0.29745 s | 3.22–3.44× | 0.98452 → 0.29643 s | 110.58 → 52.34 MiB |
| Count / 20,000 CDS, gzip | 0.36005 → 0.06256 s | 5.29–5.85× | 0.35970 → 0.06223 s | 53.06 → 51.46 MiB |
| Count / 16 MiB tar.gz | 0.29297 → 0.12158 s | 2.32–2.62× | 0.29264 → 0.12131 s | 95.65 → 51.46 MiB |
| Planner main / 64 MiB raw | 0.82269 → 0.13410 s | 6.07–6.33× | 0.81170 → 0.12426 s | 237.34 → 65.55 MiB |
| Planner main / 64 MiB gzip | 1.00513 → 0.31391 s | 3.13–3.32× | 1.00181 → 0.31124 s | 237.27 → 65.36 MiB |

The 32-record planner control is 1.933 → 2.083 ms with paired speedups
0.802–1.128×; no small-planner gain is claimed. Full-record controls retain
complete sequence fingerprints: 16 MiB raw is 0.20002 → 0.20170 s
(0.962–1.030×), and gzip 0.24850 → 0.25239 s (0.974–1.002×). Their worker peak
RSS is 110.67 → 111.47/111.52 MiB, about 0.8 MiB higher. No full-record speed or
memory gain is claimed. Large RAM savings concern wrapped chromosome inputs;
one unwrapped chromosome line still requires its line allocation.

All 324 final workers match complete counts, discovered tasks, planner JSON,
stdout/stderr or full header/sequence fingerprints. The 59 focused tests cover
raw/gzip/tar member selection and fallback, literal/empty/duplicate/Unicode
headers, strict UTF-8, missing trailers, actual CRC bit flips, record-boundary
reads and early-return closure. Existing local formatting and dry-run callers
also pass. Independent reviews checked the shared opener, source selection,
record/error consumption and the full-reader control.

Final data are `header/confirmation-api.json`, `header/confirmation-plan.json`,
`header/confirmation-records.json` and `header/evidence.json`.

## Reproduction

Save the baseline without changing branches. The complete support archive also
provides the baseline annotation package required by the metadata loader:

```bash
gg_round5_baseline=tmp/performance-20261023-round5/reproduced-baseline
mkdir -p "$gg_round5_baseline"
git archive 76525a358c1bd393f53ac54abe6bb72c6ff5362a workflow/support |
  tar -x -C "$gg_round5_baseline"
export GG_TEST_RUNTIME=docker
export GG_CONTAINER_DOCKER_IMAGE=local/genegalleon:performance-round3-20261023
```

Run sequentially, without concurrent builds, tests or benchmarks:

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_zip_readers.py \
  --before "$gg_round5_baseline/workflow/support/gene_family_output_store.py" --pairs 5 \
  --output tmp/performance-20261023-round5/zip-reproduced.json
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_gff_fasta_materialization.py \
  --before "$gg_round5_baseline/workflow/support/gff2genestat.py" \
  --profiles small-raw normal-raw normal-gzip --scopes cli cache-cold cache-warm --pairs 5 \
  --output tmp/performance-20261023-round5/gff-reproduced.json
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_fasta_metadata.py \
  --baseline "$gg_round5_baseline/workflow/support/format_species_annotation/common.py" \
  --scopes discovery,count --pairs 5 \
  --output tmp/performance-20261023-round5/header-api-reproduced.json
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_fasta_metadata.py \
  --baseline "$gg_round5_baseline/workflow/support/format_species_annotation/common.py" \
  --profiles small-raw,chr64-raw,chr64-gz --scopes plan-main --pairs 5 \
  --output tmp/performance-20261023-round5/header-plan-reproduced.json
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_fasta_metadata.py \
  --baseline "$gg_round5_baseline/workflow/support/format_species_annotation/common.py" \
  --profiles chr16-raw,chr16-gz --scopes records --pairs 5 \
  --output tmp/performance-20261023-round5/header-records-reproduced.json
```

Metadata comparisons can reuse already generated fixtures through `--fixtures`;
the original comparisons used `header/baseline-before-edit-fixtures`. With no
override, the controller creates fresh seeded inputs before measured workers.
Final measured production SHA256 values are:

| Source | SHA256 |
| --- | --- |
| `workflow/support/gene_family_output_store.py` | `c79831118c7d51417adb66f568aa0470d92d37af1293f6bbb811a5b23922c8e4` |
| `workflow/support/gff2genestat.py` | `e34c3b97098c461e3b94aade699722c40335870d84019b87665d6193decf5a6b` |
| `workflow/support/format_species_annotation/common.py` | `21e90b5df49722d4ac2aac6560f77fc07795800cf71b8fa17e62a5a99da5607e` |

## Validation

The [development guide](development-and-tests.md#choose-checks-for-a-change)
remains the check-selection reference. Commands below used the same normal
freshness policy and Docker image as the comparisons:

```bash
bash workflow/tests/run_in_runtime.sh python -m pytest -q \
  workflow/tests/test_zip_reader_pool.py workflow/tests/test_gene_family_output_store.py \
  workflow/tests/test_query_gene_orthologs.py -x -p no:cacheprovider
bash workflow/tests/run_in_runtime.sh python -m pytest -q \
  workflow/tests/test_zip_reader_pool.py -x -p no:cacheprovider
bash workflow/tests/run_in_runtime.sh python -m pytest -q \
  workflow/tests/test_gff2genestat.py workflow/tests/test_gff_source_contract.py \
  workflow/tests/test_cds_resolution.py workflow/tests/test_synteny_search_integration.py \
  workflow/tests/test_synteny_neighbor_thread_budget_runtime.py -x
bash workflow/tests/run_in_runtime.sh python -m pytest -q \
  workflow/tests/test_format_species_fasta_metadata.py workflow/tests/test_format_species_inputs.py \
  workflow/tests/test_format_species_inputs_dry_run.py -x -p no:cacheprovider
```

| Selection | Executed result |
| --- | --- |
| ZIP store/ortholog plus initial nine new cases | 261 passed in 75.90 s |
| Expanded final ZIP regression file | 15 passed in 0.15 s; six additional distinct cases |
| GFF/source/resolution/synteny, including 19 new cases | 211 passed in 18.29 s |
| FASTA metadata/local formatting/dry run, including 59 new cases | 442 passed in 46.05 s |

These are **920 distinct tests**, all passing with no selected skips. The
expanded ZIP file repeats nine cases, counted once in that total. Regression
coverage includes same-epoch replacement/in-place corruption, index refresh,
thread isolation, active nested handles and overflow, close failures, bounded
reader lifetime, source SHA/CRC failures, previous-output protection, all
sequence-consuming GFF modes and exact late parser errors.

Final host lint passed with Python 3.14.7 and Ruff 0.16.7. As in the previous
review, the ignored parser wrapper dispatches to GNU Bash in the same immutable
image with a read-only checkout and networking disabled:

```bash
env PATH="/Users/kf/miniforge3-arm64/bin:$PATH" \
  GG_LINT_BASH="$PWD/tmp/performance-20261023-round3/bash-parser" bash ./dev lint
git diff --check
```

Lint covers shell syntax, Ruff and forwarding metadata for eight entrypoints
and 25 common parameters; it is static evidence. The behavior evidence is
Docker execution. SIF, network-storage performance, a complete scientific
workflow and CSUBST/GRAMPA scheduling changes were not exercised. Scientific
settings/models, CLI and data schemas, archive layouts and upstream defaults
are preserved. VERSION remains 0.8.388 for this local-only commit.
