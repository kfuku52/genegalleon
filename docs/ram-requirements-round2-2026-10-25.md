# RAM requirements review, round 2 (2026-10-25)

Twelve proposals received separate subagent reviews, followed by independent
review and before/after measurements of five candidates. The saved baseline is
`ce0fc174a6eb4ed2f1a7c132797f5549370c4cca` (0.8.409), after the
[first RAM review](ram-requirements-review-2026-10-25.md). Published main changes
through 0.8.411 were integrated before final validation. They do not change the
five measured functions. Scientific settings, input validation, output formats,
and scheduler memory requests are preserved.

## Twelve proposals

| # | Proposal | Review and decision |
| --- | --- | --- |
| 1 | Decode trait-index rows directly from the SQL cursor | Accepted: dense queries reduce peak RSS by 13–17%. Keep the connection open through construction, restore index/dtypes, and retain source/index proof checks. |
| 2 | Reuse BH-FDR arithmetic buffers and allocate the result later | Accepted: avoid overlapping sort, rank and output arrays. Preserve stable ordering, multiply/divide order, missing values and the caller's storage. |
| 3 | Retain only IDs and lengths during constraint validation | Accepted for calls without `--output`. Consume and validate the complete FASTA/gzip input. Calls writing selected FASTA retain their validated sequence snapshot. |
| 4 | Spool multispecies promoter FASTA output | Deferred: benefit depends on multispecies inputs, adds disk I/O and complicates stdout/error preservation; each species' captured output can still be large. |
| 5 | Bound SQL input buffering for global database FDR | Accepted: SQL buffering alone reduces large-case RSS by 52–54%. Accumulate rowids and probabilities separately in typed buffers using one SQL cursor. Preserve the complete global family, exact 64-bit identities, transaction rollback and metadata. |
| 6 | Stream HGT candidate SQL queries and concatenation | Deferred: identity-only lists are small and complete family branches remain necessary. A generator around concatenation does not remove the dominant frame allocations. |
| 7 | Project selected copy-number counts while reading | Deferred: full-file numeric validation, identifier inference and parser quirks require a separate equivalence study; deleting a temporary alone may not reduce peak RSS. |
| 8 | Stream promoter N-coordinate records | Accepted: release sequence bodies after extracting coordinates/lengths, then apply the original maximum-length padding. Preserve literal `>` splitting and legacy edge cases. |
| 9 | Reduce correlation block temporaries | Deferred: the full matrix/correlation result dominates, with uncertain whole-command savings and additional special-value/order contracts. |
| 10 | Replace GFF source-contract `fetchall()` with `fetchmany(2)` | Rejected premise: the sequence store's identifier primary key already limits a valid lookup to one row. A benefit on an invalid hand-made database is not a representative optimization. |
| 11 | Stream legacy GRAMPA input lines | Deferred: the normal NWKIT path is unaffected and no representative large legacy workload was established. |
| 12 | Compact or spill scaffold-taxonomy metadata | Deferred: global seven-rank, gene, locus, host and label consistency checks require substantial state. Dense result contexts remain large; fixture allocation must also be separated from worker measurements. |

The CSV chunk readers rejected in round 1 remain rejected. This round's SQL
chunk reader does not use the affected pandas CSV parser.

## Measurement boundaries

Fixtures are synthetic and prepared outside fresh measured workers. Each
comparison alternates saved/current order, discards one warmup and reports the
median of three measured samples. Timers exclude imports, fixture preparation
and result hashing; constraint validation includes its own complete input
hashing. Linux `/proc/self/status` VmHWM includes worker imports and input
loading, and is captured before benchmark output hashing.

The freshness-checked GeneGalleon Docker image is
`sha256:340685eba10a8bf74ad03b85b1315107bdfb8dba84a37dc921b96b3e18298339`,
Linux/aarch64, Python 3.12.15, pandas 3.0.6 and NumPy 1.26.4. Docker Desktop
provided 12 CPUs and 8,214,585,344 bytes of RAM, without an additional benchmark
container limit. Our benchmarks ran serially, separately from our tests/builds.
Other work was active on the shared desktop during some samples. RSS remained
stable; timings describe local tradeoffs, not production throughput guarantees.
Raw sources, fixtures, hashes and samples are retained in
`tmp/ram-round2-20261025-pOeib1/`.

These are helper/phase measurements, not total workflow RAM requirements or SIF
qualification. Savings from the BH helper and database buffering overlap and
must not be added together.

## Results

### BH-FDR: 5,000,000 probabilities (small: 4,096)

| Profile | Peak RSS before → after (MiB) | Reduction | Wall before → after (s) | CPU before → after (s) |
| --- | ---: | ---: | ---: | ---: |
| dense | 425.55 → 286.86 | 32.6% | 0.7303 → 0.7436 | 0.7301 → 0.7434 |
| sparse | 237.59 → 209.89 | 11.7% | 0.1460 → 0.1420 | 0.1460 → 0.1419 |
| ties | 402.53 → 270.26 | 32.9% | 0.2402 → 0.2387 | 0.2401 → 0.2387 |
| small | 72.44 → 72.34 | 0.1% | 0.0004 → 0.0003 | 0.0004 → 0.0003 |

### Constraint: 128 taxa × 1,000,000 columns

| Profile | Peak RSS before → after (MiB) | Reduction | Wall before → after (s) | CPU before → after (s) |
| --- | ---: | ---: | ---: | ---: |
| Plain, validate all | 175.75 → 53.16 | 69.8% | 0.3029 → 0.3361 | 0.3029 → 0.3361 |
| Gzip, validate all | 175.91 → 53.19 | 69.8% | 0.9890 → 1.0058 | 0.9889 → 1.0058 |
| Plain, output all | 175.73 → 175.72 | 0.0% | 0.6135 → 0.6181 | 0.6133 → 0.6166 |
| Gzip, output subset | 175.76 → 175.76 | 0.0% | 1.1325 → 1.1329 | 1.1320 → 1.1324 |
| Small, plain validation | 50.87 → 50.85 | 0.0% | 0.0004 → 0.0005 | 0.0004 → 0.0005 |

The other all/subset/plain/gzip combinations also match completely. Large
validation RSS falls 69.7–69.8%; output calls stay near 176 MiB.

### Warm trait index: narrow 150,000 × 12 / wide 100,000 × 128 numeric columns

| Profile | Peak RSS before → after (MiB) | Reduction | Wall before → after (s) | CPU before → after (s) |
| --- | ---: | ---: | ---: | ---: |
| narrow | 316.81 → 263.95 | 16.7% | 1.5176 → 1.5024 | 1.3933 → 1.3711 |
| wide | 1152.47 → 1000.61 | 13.2% | 6.5882 → 6.5210 | 6.2104 → 6.1166 |
| sparse | 77.43 → 77.47 | -0.1% | 0.3469 → 0.3323 | 0.2994 → 0.3080 |
| small | 77.34 → 77.30 | 0.1% | 0.0064 → 0.0065 | 0.0042 → 0.0046 |

Sparse selects 20 IDs from the wide table; small has 64 rows and eight numeric
columns. Gene ID is an additional column. Sparse/small RSS changes are immaterial.

### Database FDR: 2,000,000 rows, SQL buffering alone (fixed old BH)

| Profile | Peak RSS before → after (MiB) | Reduction | Wall before → after (s) | CPU before → after (s) |
| --- | ---: | ---: | ---: | ---: |
| dense | 579.23 → 280.46 | 51.6% | 4.9187 → 5.1679 | 4.7012 → 4.8891 |
| textual | 558.88 → 256.77 | 54.1% | 4.7382 → 4.6510 | 4.5939 → 4.5041 |
| small | 104.26 → 106.80 | -2.4% | 0.0080 → 0.0075 | 0.0051 → 0.0048 |

### Database FDR: combined SQL buffering and BH changes

| Profile | Peak RSS before → after (MiB) | Reduction | Wall before → after (s) | CPU before → after (s) |
| --- | ---: | ---: | ---: | ---: |
| dense | 579.30 → 226.20 | 61.0% | 4.7666 → 4.6367 | 4.5777 → 4.3712 |
| missing | 535.76 → 190.62 | 64.4% | 4.4294 → 4.3376 | 4.2727 → 4.1954 |
| textual | 558.46 → 210.88 | 62.2% | 5.0360 → 5.3750 | 4.7518 → 4.7791 |
| small | 104.26 → 106.79 | -2.4% | 0.0083 → 0.0073 | 0.0044 → 0.0045 |

The missing profile has 75% NULL P values; textual mixes numeric text, floats,
integers and NULLs. Small has 128 rows. The small case adds about 2.5 MiB of
process RSS and has millisecond timings. Isolated dense buffering costs about
5% wall time (4% CPU) while halving RAM; the combined dense measurement is
slightly faster. No general throughput claim is made from these small differences.

### Promoter coordinates: 50,000 × 4,096 bases / one 64-million-base record

| Profile | Peak RSS before → after (MiB) | Reduction | Wall before → after (s) | CPU before → after (s) |
| --- | ---: | ---: | ---: | ---: |
| runs | 476.41 → 90.15 | 81.1% | 2.0404 → 1.9046 | 2.0402 → 1.9015 |
| runs-gzip | 476.44 → 90.67 | 81.0% | 2.1367 → 1.9583 | 2.1284 → 1.9520 |
| no-n | 476.42 → 86.06 | 81.9% | 0.7143 → 0.5571 | 0.7109 → 0.5537 |
| long-run | 315.66 → 264.55 | 16.2% | 0.3413 → 0.2796 | 0.3397 → 0.2783 |

Dense-coordinate control (500 × 4,096 alternating bases): 83.62 → 83.68 MiB, 0.1820 → 0.1736 s wall, 0.1776 → 0.1717 s CPU.
Its RSS is unchanged. The eight-record × 64-base controls (long-run: 1,024
bases) also match for all five profiles; RSS stays near 81 MiB and timings
are below 0.001 s. These import/AST-allocation floors are not useful speed
or memory improvement claims.

BH comparisons hash every output float as uint64 bits and check that input
storage is unchanged. Tests also cover readonly/noncontiguous arrays, nullable
inputs, signed zero, subnormals, ties and invalid probability boundaries.

Constraint comparisons match complete JSON reports and selected FASTA hashes,
including small/large, plain/gzip, all/subset and validation/output calls. Invalid
inputs must preserve existing outputs, consume gzip tails and retain error
priority. No speed improvement is claimed for the unchanged output path or
sub-millisecond controls. Plain validation trades some CPU time for lower RAM.

Trait measurements use a content-verified warm index. A measured worker fails
if the source loader runs. Every returned value, row order, duplicate, column
dtype, index dtype and index name is compared. Sparse queries and small tables
are controls. Changing the indexer source hash triggers the existing automatic
cache rebuild on the first real read; these measurements do not represent cold
construction. Decoded rows and the returned DataFrame still scale with selection.

Database comparisons hash the complete schema, metadata and every typed SQLite
value, including float bits. Fixtures include numeric text/NULLs, a shadowed
`rowid` column, negative real rowids, identities above 2**53 and near 2**63. A
fixed old BH implementation isolates SQL-buffer savings; the combined comparison
uses each revision's BH. Probability validation precedes writes, and late-invalid
chunks and injected later-update failures retain rollback coverage. The global
sort and typed buffers remain proportional to the number of candidate rows.

Promoter comparisons match complete DataFrames including dtypes and row order.
The initial streamed implementation increased RAM for a single 64-million-base
record. It was refined to release raw records/regex matches earlier and use
literal LF replacement without regex fragment allocation. The final large-record
control reduces RAM too. The parser still needs memory for the largest current
record plus coordinate strings/metadata; it is not constant-memory for arbitrary
record lengths or dense N-coordinate output.

## Reproduce

Run from the repository root, with the image above selected:

```bash
export GG_TEST_RUNTIME=docker
export GG_CONTAINER_DOCKER_IMAGE=sha256:340685eba10a8bf74ad03b85b1315107bdfb8dba84a37dc921b96b3e18298339
BASE=$(mktemp -d "$PWD/tmp/ram-round2-baseline-XXXXXX")
for file in csubst_scan_support.py trait_table_index.py species_tree_constraint.py generate_orthogroup_database.py orthogroup_statistics.py; do
  git show ce0fc174a6eb4ed2f1a7c132797f5549370c4cca:workflow/support/"$file" > "$BASE/$file"
done
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_bh_memory.py --baseline "$BASE/csubst_scan_support.py" --output "$BASE/bh.json"
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_constraint_memory.py --baseline "$BASE/species_tree_constraint.py" --scopes validate output --output "$BASE/constraint.json"
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_trait_index_memory.py --baseline-support "$BASE" --fixtures "$BASE/trait-fixtures" --output "$BASE/trait.json"
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_database_fdr_memory.py --baseline "$BASE/generate_orthogroup_database.py" --bh-source "$BASE/csubst_scan_support.py" --profiles dense textual small --output "$BASE/database-isolated.json"
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_database_fdr_memory.py --baseline "$BASE/generate_orthogroup_database.py" --output "$BASE/database-combined.json"
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_promoter_n_runs.py --baseline "$BASE/orthogroup_statistics.py" --profiles runs runs-gzip no-n long-run --records 50000 --length 4096 --long-length 64000000 --repeats 3 --output "$BASE/promoter.json"
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_promoter_n_runs.py --baseline "$BASE/orthogroup_statistics.py" --profiles alternating --records 500 --length 4096 --long-length 1024 --repeats 3 --output "$BASE/promoter-alternating.json"
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_promoter_n_runs.py --baseline "$BASE/orthogroup_statistics.py" --records 8 --length 64 --long-length 1024 --repeats 3 --output "$BASE/promoter-small.json"
```

Each new harness supports `--baseline-only`, used before its production edit.
Preserve the output directory for inspection; the commands create no curated
workflow inputs or production outputs.

## Validation

Final checks used the normal freshness-checked Docker wrapper after integrating
published main changes and bumping `VERSION` to 0.8.412. All 597 selected tests
passed, with no skips:

```bash
bash ./dev check fast workflow/tests/test_get_trait_matrix_outputs.py workflow/tests/test_gff2genestat.py workflow/tests/test_orthogroup_statistics_promoter.py workflow/tests/test_generate_orthogroup_database.py workflow/tests/test_csubst_scan_candidate_sites.py workflow/tests/test_plot_csubst_aa_change_summary.py workflow/tests/test_orthogroup_statistics_fimo.py -x
bash workflow/tests/run_in_runtime.sh python -m pytest -q workflow/tests/test_species_tree_constraint.py -x
```

The fast selection passed 551 tests in 36.61 s. Constraint checks passed 46 in
4.81 s, including the real IQ-TREE constraint integration. Initial focused runs
also passed; the final run covers the reference-lifetime refinement and its
integration with the independently published FIMO optimization.

Static checks passed separately:

```bash
GG_LINT_BASH="$PWD/tmp/ram-round2-20261025-pOeib1/lint-bash" bash ./dev lint
git diff --check
```

The macOS host has Bash 3.2. As documented in the development guide, the temporary
`GG_LINT_BASH` wrapper delegated syntax checks to Bash in the same Docker image,
with this checkout mounted at its absolute working path. Ruff and configuration
checks ran on the host; all eight entrypoints and 27 common parameters passed.
The temporary syntax-check container was stopped afterward. Static checks are
separate evidence from the executed Docker behavior tests.

Full scientific workflows, the full runtime/R lanes and SIF validation were not
run for these focused helper changes. No SIF compatibility is claimed from the
Docker results. Curated data and benchmark artifacts were preserved.
