# Workflow performance review, round four, 2026-10-23

This review starts from `f059e6644000785ad1721b73664a221b83e70a17`
(0.8.388), after the [third review](performance-round3-2026-10-23.md).
Three agents investigated PGLS correction grouping, expression-correlation
clade minima and BUSCO annotation aggregation. Baselines were measured before
production edits. Benchmarks and runtime checks ran sequentially to avoid
resource contention; accepted changes received an independent code review.

## Comparison environment

The normal freshness-checking wrapper selected GeneGalleon Docker image
`local/genegalleon:performance-round3-20261023`, immutable ID
`sha256:0e46be37b7ddbdb60b331ef35017abd3f2d90820d152ce3838ec369537791752`,
with runtime-input fingerprint
`8f5655e6d387393268bb292d495b28222ad68b0e072a8a622722f6daef30faa4`.
The environment is Linux aarch64, Python 3.12.15, pandas 3.0.6, NumPy 1.26.4
and R 4.4.3. Docker Desktop provides 12 CPUs and 8,214,585,344 bytes of RAM;
no additional container limit was applied. These Python-only workflow changes
use the mounted checkout and do not change container inputs.

Each final comparison uses five alternating fresh-worker pairs after a
discarded warmup pair. Reported times are medians; ranges are paired wall-time
speedups. Linux `VmHWM` is captured before output verification and measures
the complete worker, including fixture/setup allocations. Fixture generation
and output hashing are outside timers. Synthetic fixtures and warm file caches
do not establish whole-analysis speedups or production storage throughput.
Raw measurements, hashes, differential cases and logs remain under
`tmp/performance-20261023-round4/`. SIF was not exercised on this macOS host.

## PGLS correction grouping

Accepted. Group positions replace iteration over full diagnostic DataFrames.
Family-relative positions are mapped to the reset output index, and positional
selection reads the existing numeric p-values. Global masking and the existing
Holm/BH helper remain unchanged. Failed-row denominators, correction scopes,
diagnostic columns, input immutability and reset row order are preserved.

The fixture has 48 fits per family: two methods, three aggregations, four
responses and two predictors. The wide version uses the actual producer's
67-column shape, including its diagnostics; the narrow control has 11 columns.
Complete database
measurements call the real summary database main with 256 live or archived
families. Timers exclude imports; database timers include input reads, all table
writes and publication.

| Profile | Wall before → after | Paired speedup range | CPU before → after | Peak RSS before → after |
| --- | ---: | ---: | ---: | ---: |
| 32 families, wide | 0.02253 → 0.01599 s | 1.17–1.49× | 0.02253 → 0.01599 s | 89.87 → 89.14 MiB |
| 1,000 families, wide | 0.39555 → 0.28416 s | 1.31–1.62× | 0.39555 → 0.28413 s | 253.13 → 229.88 MiB |
| 4,000 families, wide | 1.69549 → 1.14748 s | 1.32–1.51× | 1.69539 → 1.14742 s | 697.35 → 659.90 MiB |
| 1,000 families, narrow | 0.26322 → 0.23258 s | 1.12–1.30× | 0.26320 → 0.23255 s | 134.35 → 123.87 MiB |
| 1,000 families, all failed | 0.19657 → 0.16838 s | 1.05–1.39× | 0.19658 → 0.16838 s | 228.60 → 229.88 MiB |
| 1,000 families, many response/predictor groups | 3.38258 → 2.98099 s | 1.10–1.23× | 3.38239 → 2.98075 s | 253.47 → 229.54 MiB |
| Complete database, live inputs | 1.12625 → 1.04627 s | 1.004–1.113× | 1.12930 → 1.04886 s | 238.22 → 238.18 MiB |
| Complete database, ZIP inputs | 1.12581 → 1.10785 s | 0.997–1.068× | 1.12822 → 1.11070 s | 239.76 → 239.60 MiB |

All 30 kernel comparisons improve. The ZIP database timings overlap noise;
no ZIP whole-command gain or database memory improvement is claimed. The
all-failed control uses about 1.3 MiB more peak RSS. The first prototype used
label-based family lookups and regressed on narrow inputs; positional lookup
resolved that regression before the final comparisons above.

Every kernel pair matches complete typed frames, row/index order and TSV bytes.
Database verification compares full schema/index SQL, column definitions and
every typed row, canonically ordered because concurrent readers can change
incidental insertion order. Additional differential checks cover 240 cases:
152 exact typed/output matches and 88 exact exception type/message matches,
including nullable p-values, string/Decimal inputs, unused categorical keys,
missing identities and duplicate original indices. The focused Docker tests
also execute the real R BH reference, live/ZIP readers and append behavior.

Final samples are in `pgls/confirmation-iloc.json` and `pgls/summary.json`.

## Expression-correlation clade minima

Accepted with the measured fallback overhead below. Binary clades reuse child
minima and inspect both directed cross-child blocks, retaining diagonal entries.
Clades smaller than 64 descendants and multifurcations retain the original square
reduction; the cutoff avoids expensive tiny-block reductions, rather than
changing a scientific criterion. Reuse is limited to plain C-contiguous native
float64 square arrays and ordinary valid integer indices. Other inputs preserve
the original calculation and error paths. If a combined minimum is NaN, zero or
infinite, the node uses its original reduction to preserve NaN payloads and
signed-zero selection. Both caches follow the existing branch-label overwrite
semantics.

All profiles use seeded synthetic expression values across 32 tissues. The
804-tip case reproduces the 199-node spine and deepest polytomy of the existing
`test_species_mapping_clone_handles_observed_deep_tree_without_recursion`
fixture. It does not reuse biological expression measurements. The controller
sets `OPENBLAS_NUM_THREADS=1`, `OMP_NUM_THREADS=1` and `MKL_NUM_THREADS=1` for
workers, preventing fixture correlation calculations from leaving spinning BLAS
threads in the helper timer. A stdlib-only fresh inner worker measures the real
CLI's imports, parser/main and both output tables; interpreter startup is
excluded. All installed dependencies are used, without import stubs.

| Helper profile | Wall before → after | Paired speedup range | Peak RSS before → after |
| --- | ---: | ---: | ---: |
| 512-tip pectinate | 0.18731 → 0.01495 s | 12.27–12.92× | 122.26 → 122.40 MiB |
| 1,024-tip pectinate | 1.58607 → 0.04950 s | 25.10–32.68× | 136.11 → 128.92 MiB |
| 804 tips, 199-node spine | 0.44158 → 0.01431 s | 30.20–31.52× | 127.14 → 125.42 MiB |
| 1,024-tip balanced | 14.309 → 10.562 ms | 1.31–1.41× | 135.79 → 128.95 MiB |
| 1,024-tip star | 5.049 → 5.506 ms | 0.895–0.941× | 133.43 → 132.53 MiB |
| 1,024 tips, 25% missing | 1.050 → 1.278 ms | 0.800–0.858× | 128.95 → 128.97 MiB |
| 1,024-tip pectinate, deepest row constant | 1.68876 → 1.75837 s | 0.938–0.979× | 136.55 → 136.80 MiB |

Helper CPU medians closely follow wall time. Small balanced controls with
16/32/64/128 tips and a 32-tip tree with 2,048 expression rows are also measured.
The 32-tip control costs 0.1527 → 0.1692 ms, and the extra-catalog control costs
0.1663 → 0.1838 ms. These cache/bookkeeping costs are real; no universal helper
speed or memory improvement is claimed.

| Complete CLI profile | Wall before → after | Paired speedup range | CPU before → after | Peak RSS before → after |
| --- | ---: | ---: | ---: | ---: |
| 512-tip pectinate | 0.90978 → 0.74859 s | 1.17–1.23× | 0.90419 → 0.74392 s | 161.99 → 161.32 MiB |
| 1,024-tip pectinate | 2.92366 → 1.46279 s | 1.98–2.02× | 2.91766 → 1.45794 s | 211.24 → 203.48 MiB |
| 804 tips, 199-node spine | 1.10215 → 0.67286 s | 1.58–1.68× | 1.09695 → 0.66781 s | 168.19 → 165.13 MiB |
| Deepest expression row constant | 3.12592 → 3.19627 s | 0.965–0.985× | 3.11992 → 3.19021 s | 211.25 → 211.38 MiB |

The constant deepest row generates real `corrcoef` NaNs at every ancestor. Its
preserved original reductions incur about 4.1% helper and 2.25% complete-CLI
overhead, a consistent worst-case regression in these measurements. Small,
missing-tip, extra-catalog and star CLI paired ranges cross unity, so no speed
change is claimed there. These are helper and expression-annotation CLI results,
not complete scientific workflow multipliers.

All 240 final workers match every helper key/order/scalar type/float bit pattern
or both complete CLI TSVs and stdout/stderr, including constant-expression
warnings. The 148 focused oracle cases additionally cover quiet/signaling NaN
payloads, signed zero, infinities, asymmetric blocks, unusual indices, duplicate
labels, unsupported arrays and exact warning/error behavior under NumPy warn and
raise policies. They include active large-clade reuse paths. Independent review
found and corrected preloaded imports in the first CLI harness; only the cold
import comparisons above are used for whole-command claims.

Final samples are `expression/confirmation-kernel.json`,
`expression/confirmation-cli.json`, `expression/confirmation-cli-star.json` and
`expression/confirmation-constant.json`; environment/source hashes are in
`expression/evidence.json`.

## BUSCO annotation aggregation

Accepted with the small-input tradeoff below. After the existing complete parse
and explicit string/missing-value conversion, singleton genes bypass string
aggregation. Mixed inputs aggregate repeated genes only, then restore first-seen
gene order. Fully repeated inputs retain the original aggregation. This changes
`merge_cds_annotation.load_busco`, a different path from the cross-species BUSCO
collector optimized in the previous review.

Mixed fixtures contain eight singleton genes and one two-row gene per ten
markers. Full CLI fixtures also include a CDS FASTA with four times as many
genes as BUSCO markers, UniProt metadata, two expression columns and sequence
statistics, followed by the actual joins and complete output write. Both 300nt
and 900nt CDS inputs are measured. CLI timers include imports and main; loader
timers exclude imports. These results concern annotation merging, not BUSCO
search or the complete genome annotation workflow.

| Markers / profile / CDS length | CLI wall before → after | Paired speedup range | CPU before → after | Peak RSS before → after |
| --- | ---: | ---: | ---: | ---: |
| 32 / mixed / 300nt | 0.16263 → 0.16965 s | 0.894–0.987× | 0.16100 → 0.16927 s | 74.85 → 74.85 MiB |
| 5,000 / mixed / 300nt | 0.41060 → 0.26330 s | 1.50–2.24× | 0.40997 → 0.26279 s | 88.29 → 88.39 MiB |
| 15,000 / mixed / 300nt | 0.91081 → 0.47504 s | 1.80–2.43× | 0.91017 → 0.47439 s | 111.84 → 112.66 MiB |
| 5,000 / singleton / 300nt | 0.42204 → 0.24642 s | 1.65–1.80× | 0.42144 → 0.24584 s | 88.24 → 87.93 MiB |
| 15,000 / singleton / 300nt | 0.90148 → 0.38883 s | 2.26–2.50× | 0.90081 → 0.38837 s | 112.17 → 111.47 MiB |
| 15,000 / all repeated / 300nt | 0.54297 → 0.53062 s | 1.002–1.051× | 0.54249 → 0.53025 s | 109.74 → 108.69 MiB |
| 5,000 / mixed / 900nt | 0.36777 → 0.24189 s | 1.47–1.61× | 0.36736 → 0.24155 s | 89.80 → 88.45 MiB |
| 15,000 / mixed / 900nt | 0.85317 → 0.42475 s | 1.95–2.10× | 0.85277 → 0.42438 s | 112.67 → 112.72 MiB |

All larger mixed/singleton CLI pairs improve. The small control consistently
costs about 7 ms more (4.3%); this measured regression is not dismissed as noise.
No all-repeated speedup or general memory reduction is claimed. The loader-only
900nt profiles are 0.17064 → 0.03637 s for 5,000 markers and
0.51675 → 0.09698 s for 15,000. The latter uses 85.46 → 86.04 MiB peak RSS.
The all-repeated loader's paired speedups span 0.981–1.014×.

All 192 final workers match complete output bytes and timestamp-normalized logs;
loaders also compare complete typed DataFrames, including index/column classes.
An additional 200 randomized tables and 11 edge cases match: 208 successful or
missing results and three exact exception type/message comparisons. Tests cover
literal identifiers, interleaved duplicates, Unicode/quoted tabs and newlines,
numeric missing-value spelling and late UTF-8/CSV errors that preserve an
existing output file. An independent agent reviewed the implementation, tests
and measurement harness. Final raw samples are `busco/confirm-*.json`, with
commands, runtime proof, source hashes and summary in `busco/receipt.json`.

## Reproduction

Save baseline files without changing branches, then run comparisons sequentially:

```bash
gg_round4_baseline=tmp/performance-20261023-round4/reproduced-baseline
mkdir -p "$gg_round4_baseline"
git archive f059e6644000785ad1721b73664a221b83e70a17 workflow/support |
  tar -x -C "$gg_round4_baseline"
export GG_TEST_RUNTIME=docker
export GG_CONTAINER_DOCKER_IMAGE=local/genegalleon:performance-round3-20261023
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_pgls_multiplicity.py \
  --before "$gg_round4_baseline/workflow/support/pgls_multiplicity.py" \
  --output tmp/performance-20261023-round4/pgls-reproduced.json --pairs 5
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_expression_clade_minima.py \
  --baseline "$gg_round4_baseline/workflow/support/orthogroup_statistics.py" \
  --profiles clade16,small,clade64,clade128,pectinate512,pectinate1024,deep804,balanced1024,star1024,missing1024,extra_catalog \
  --scopes kernel --pairs 5 --output tmp/performance-20261023-round4/expression-kernel-reproduced.json
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_expression_clade_minima.py \
  --baseline "$gg_round4_baseline/workflow/support/orthogroup_statistics.py" \
  --profiles small,pectinate512,pectinate1024,deep804,missing1024,extra_catalog,star1024 \
  --scopes cli --pairs 5 --output tmp/performance-20261023-round4/expression-cli-reproduced.json
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_expression_clade_minima.py \
  --baseline "$gg_round4_baseline/workflow/support/orthogroup_statistics.py" \
  --profiles constant1024 --scopes kernel,cli --pairs 5 \
  --output tmp/performance-20261023-round4/expression-constant-reproduced.json
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_busco_annotation.py \
  --baseline "$gg_round4_baseline/workflow/support/merge_cds_annotation.py" \
  --markers 15000 --profile mixed --cds-codons 300 --repeats 5 \
  --output tmp/performance-20261023-round4/busco-reproduced.json
```

For BUSCO's other controls, use marker counts 32/5,000/15,000 and profiles
`mixed`, `singleton` or `duplicated` as listed above. The default
`--cds-codons 100` produces 300nt CDS; `--cds-codons 300` produces 900nt CDS.
Use distinct output paths. Baseline/source/input/output hashes are retained in
the raw records. The final measured production SHA256 values are:

| Source | SHA256 |
| --- | --- |
| `workflow/support/pgls_multiplicity.py` | `f36a421508f10fff32bea7461bd8698fa3682d9d7920b90ff7777f58e9894ae4` |
| `workflow/support/orthogroup_statistics.py` | `bd4c5d9ab441c4e4d6a04b100d1e4ec11e224ce47f9ddc1b63c1fb880731bd67` |
| `workflow/support/merge_cds_annotation.py` | `b8bc135db8520acb4606f85739592ddcc111d8bcd090bc96d14296faf674adf6` |

## Validation

The [development guide](development-and-tests.md#choose-checks-for-a-change)
remains the check-selection reference. Focused checks used the same Docker
image and normal freshness policy as the comparisons, through
`workflow/tests/run_in_runtime.sh python -m pytest -q`:

| Selected files | Result |
| --- | --- |
| `workflow/tests/test_pgls_multiplicity.py` | 15 passed in 1.56 s; includes real R BH, live/ZIP and append |
| `workflow/tests/test_orthogroup_statistics_expression.py` and `workflow/tests/test_orthogroup_statistics_fimo.py` | 194 passed in 10.09 s; exact numeric/warning/error oracles and existing CLI/FIMO behavior |
| `workflow/tests/test_merge_cds_annotation.py` | 34 passed in 1.48 s; loaders, joins, output and failure behavior |

All **243 selected tests passed**, without selected skips. Example commands with
the environment exports above are:

```bash
bash workflow/tests/run_in_runtime.sh python -m pytest -q workflow/tests/test_pgls_multiplicity.py -x
bash workflow/tests/run_in_runtime.sh python -m pytest -q \
  workflow/tests/test_orthogroup_statistics_expression.py workflow/tests/test_orthogroup_statistics_fimo.py -x
bash workflow/tests/run_in_runtime.sh python -m pytest -q workflow/tests/test_merge_cds_annotation.py -x -p no:cacheprovider
```

Final host lint passed with Python 3.14.7 and Ruff 0.16.7. The existing ignored
GNU Bash parser wrapper uses the same immutable Docker image with networking
disabled and a read-only checkout, avoiding the macOS Bash 3.2 parser limit:

```bash
env PATH="/Users/kf/miniforge3-arm64/bin:$PATH" \
  GG_LINT_BASH="$PWD/tmp/performance-20261023-round3/bash-parser" bash ./dev lint
git diff --check
```

Lint includes shell syntax, Ruff and configuration forwarding validation for
eight entrypoints and 25 common parameters; it is static evidence. Behavior
checks and measurements above are Docker evidence. SIF and a complete scientific
workflow were not run. Scientific models/settings, public interfaces, output
schemas and upstream source defaults are preserved. VERSION remains 0.8.388;
this review is a local-only commit.
