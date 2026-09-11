# Native OU default-budget measurements

This is development tuning of computational budgets for a single family search,
not a statistical-adoption, accuracy-superiority, or speedup study. GeneGalleon
now omits the candidate-pool, refit, screening and beam-width options, leaving
those defaults to the installed NWKIT version.

## Workload and measurement

`generate.py` independently simulates a balanced 1,000-tip OU tree with 100
non-overlapping shift clades and five repeated nonbaseline optimum regimes.
One- and two-trait inputs share the first trait exactly. They have no missing
observations, process tip variance 1 and independent observation variance 0.2.
The known sampling variance supplied to inference is zero. Inference estimates
alpha, process variance and additional observation variance for every trait,
using AICc, convergence enabled, automatic cap/strategy, and no bootstrap.

`protocol.json` records the settings before tuning. `environment.json` records
the GeneGalleon Docker runtime on an Apple M2 Max. Each setting runs once in a
fresh process with one BLAS/OpenMP thread and no competing numerical tests or
benchmarks. Startup is included; there is no explicit warmup. These runs are
well separated from the one-hour decision boundary, so repeated runs are not
used to estimate a runtime distribution. The one-hour harness cutoff terminates
a run and is **not** an NWKIT wall-clock deadline.

Both settings use the same frozen NWKIT algorithm source, recorded in
`source-manifest.json` and `baseline-source.tar.gz`; only the explicit budget
arguments differ. The archive and mounted modules, and all input tree/table
hashes, were checked against their manifests. Production changes centralize and
change the defaults without changing likelihood or search algorithms.

## Results and adopted defaults

NWKIT defaults change from candidate/refit/screening budgets **24 / 48 / 2,000**
to **128 / 256 / 100,000**. Beam width remains 2; lasso iterations remain 150
and the search-memory budget remains 512 MiB. No additional proposal was needed:
the tested expanded searches finished in about 8 and 13 minutes, leaving
substantial room below the requested approximately one-hour per-family target.

| Traits | Budgets | Wall seconds | Peak MiB | Auto cap | Largest fitted / selected shifts | Refits | Quick evaluations |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | baseline | 88.15 | 219.16 | 24 | 13 / 13 | 39 | 2,000 |
| 1 | expanded | 452.54 | 330.33 | 128 | 41 / 41 | 100 | 100,000 |
| 2 | baseline | 139.61 | 258.70 | 24 | 13 / 13 | 31 | 2,000 |
| 2 | expanded | 770.29 | 333.67 | 128 | 41 / 41 | 102 | 100,000 |

The expanded runs take 5.13× and 5.52× as long as their respective baselines.
They explore further, but all four runs exhaust their screening budget. In
particular, an auto cap of 128 does **not** mean 100 or 128 shifts were fitted.

| Traits | Budgets | True / false selected branches | Missed true branches | Selected groups | Tip-mean RMSE | AICc |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 1 | baseline | 13 / 0 | 87 | 4 | 4.6012 | 3939.7395 |
| 1 | expanded | 32 / 9 | 68 | 41 | 3.5491 | 3400.8934 |
| 2 | baseline | 13 / 0 | 87 | 4 | 3.8459 | 7297.3851 |
| 2 | expanded | 35 / 6 | 65 | 41 | 2.7248 | 6537.5226 |

The truth has six total optimum regimes. Expanded results still miss many true
shifts, introduce false branches, and do not recover that regime partition;
these counts must not be presented as statistical validation of the defaults.
Among layouts shared between settings (27 for one trait, 20 for two), all
log-likelihoods agree exactly; `overlap-check.json` records this check. Selected
models intentionally differ because budgets differ. All four searches completed;
no failures or timeouts were omitted.

## Reproduction

From the repository root, extract the source into an empty temporary directory:

```sh
mkdir -p /tmp/native-ou-budget-repro
tar -xzf docs/benchmarks/native-ou-default-budgets/baseline-source.tar.gz \
  -C /tmp/native-ou-budget-repro
```

Run each command sequentially, with a new label for every run:

```sh
docker run --rm \
  -v "$PWD/docs/benchmarks:/bench" \
  -v /tmp/native-ou-budget-repro:/source:ro \
  -e PYTHONPATH=/source -e OPENBLAS_NUM_THREADS=1 \
  -e OMP_NUM_THREADS=1 -e MKL_NUM_THREADS=1 \
  -w /bench/native-ou-default-budgets \
  local/genegalleon:nwkit-ou-auto-dev \
  python measure.py repeat-baseline-1trait 1 24 48 2000
```

Replace the trailing arguments with `repeat-expanded-1trait 1 128 256 100000`,
`repeat-baseline-2traits 2 24 48 2000`, and
`repeat-expanded-2traits 2 128 256 100000`. The local image identifier and
packages are in `environment.json`; another runtime need not reproduce timings.
The measurement helper reuses the sibling
`native-ou-1000tips-100shifts-aicc/benchmark.py` process wall/CPU/RSS harness.
`python docs/benchmarks/native-ou-default-budgets/summarize.py` regenerates the
summary, checking the selected AICc formula and minimum among evaluated models.

## Interpretation limits

The automatic cap is a feasibility bound, not the largest shift count actually
fitted. All runs must report both. Finite screening and refit budgets leave
incomplete coverage; changing the candidate pool does not produce a nested set
of fitted models. More candidate exploration can also increase false positives.
Exact-clade recovery and tip-mean RMSE are descriptive for this one seeded input,
not estimates of statistical performance across a population.

Counts do not guarantee completion within one hour on other hardware, larger
or differently shaped trees, more traits, missing-data patterns or harder fits.
Calibration and support bootstraps repeat searches and are outside this timing
scope. The expanded counts deliberately leave room below one hour in the tested
workloads instead of attempting to spend the full hour.

## Integration validation

The intended NWKIT change was isolated from concurrent, unrelated worktree
edits for final validation. Its 220 Python modules were compared with the frozen
benchmark manifest: only the three default/CLI modules differ.

- GeneGalleon Docker runtime: NWKIT native tests plus `test_shift.py`: **208
  passed, 4 skipped**. The four skips require an explicitly configured legacy
  R/kfl1ou backend and are unrelated to native inference.
- GeneGalleon adapter/core integration tests: **11 passed** using that isolated
  NWKIT source, including inherited budgets, explicit adapter overrides and
  automatic-cap/shared-regime dispatch.
- GeneGalleon static shell safety/entrypoint tests: **256 passed**. Together with
  the integration tests, the final affected set has **267 passing tests**.
- Host static checks: NWKIT Ruff lint/format (479 files), mypy (220 modules) and
  repository maintainability checks passed; changed GeneGalleon/benchmark Python
  files passed Ruff and affected shell scripts passed `bash -n`.

An initial four-tip smoke input retained replicate sampling variances that
pushed a residual process-variance fit to its numerical bound. It was replaced
with complete, single-observation traits for the small dispatch fixture; the
production numerical-bound rejection was not relaxed. The affected test and
adapter suite then passed. No SIF runtime was available on this macOS host;
these are Docker runtime results, not SIF validation. The complete repository
suites and resampling-scale benchmarks were not rerun for this default-only
change.
