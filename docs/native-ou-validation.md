# Native OU implementation validation, 2026-09-10

Status: GeneGalleon uses NWKIT AICc with an automatic shift cap and convergence
search for enabled OU analysis, replacing kfl1ou at the user's request. OU analysis remains off by default. Statistical
adoption criteria remain unfulfilled; this integration decision is not evidence
that those criteria passed. See [configuration and outputs](native-ou-shifts.md).

## Automatic cap and convergence defaults, 2026-09-11

NWKIT commit `0971715` adds native `--max-shifts auto`. The requested value is
preserved in configuration/replay fingerprints; `search.shift_limit` records the
resolved cap, applicable constraints and whether budgets reduced it. Explicit
integer requests are not silently clipped. Large-tree exhaustive preflight
stops counting when its traversal bound is exceeded.

GeneGalleon now defaults to `native_ou_max_shifts="auto"`,
`native_ou_convergence=1`, `native_ou_search_strategy="auto"`, and AICc.
Its default beam cap on a 1,000-tip tree is 24 (the branch-pool budget), not an
unlimited search. Shared regimes reduce model dimension; AICc eligibility is
checked per candidate rather than imposing an unconstrained-model parameter
count on convergent layouts.

Validation for this follow-up:

- NWKIT native/legacy shift and CLI-contract suite: 264 passed, 13 skipped.
- Final limit/CLI checks, including automatic path caps and bootstrap replay:
  45 passed. Automatic and equivalent integer caps yield identical candidates,
  AICc winners and selection-support results.
- GeneGalleon adapter, shell and entrypoint checks: 265 passed, including a real
  stage run with the automatic cap and convergence defaults, shared-regime
  candidates and complete eight-file bundle publication.
- NWKIT Ruff lint/format, mypy (220 source files), and maintainability hard limits
  passed. These checks do not establish a speed or statistical-power improvement.

The runtime tests use the GeneGalleon Docker environment. An installed-package
smoke run in `local/genegalleon:nwkit-ou-auto-dev` (image
`sha256:40c8c0c545ca39445a57faac6ba8879e6e1075553426c9ee6d7a8205381b56f8`)
completed AICc/auto/convergence inference and all eight outputs without a NWKIT
source mount. All 220 installed module hashes matched the committed source.
SIF and a complete production-image rebuild remain untested. Earlier path-only timing comparisons
are not measurements of the new convergence-enabled defaults.

## Initial backend replacement validation, 2026-09-11

The initial replacement used NWKIT commit `b44a397`, which contains the native information criteria, covariance-updated
path search and optional global-null AIC gate. The initial integration selected
AICc/path, with convergence off and a ten-shift cap. The follow-up changes the
defaults to an automatic cap and convergence-capable search; OU remains opt-in.

Validation used `local/genegalleon:nwkit-ou-replacement-dev`
(image `sha256:6b1808a8e8ed078d7e7cac7edab3c01e4912c6339dfe7af0e34d1c5c17800d6b`).
This derivative of the previously validated GeneGalleon development runtime
installs the committed NWKIT source and removes the installed kfl1ou R package.
Its source manifest records the installed NWKIT commit and removes kfl1ou.
Tests use the installed NWKIT, without a NWKIT source mount.

- Relevant adapter, summary, shell, source-policy and intron tests: 289 passed.
- Entire static suite: 397 passed; both publishing workflows also passed actionlint.
- Entire fast suite: 1,784 passed (eight existing multiprocessing warnings).
- Strict runtime adapter and source-manifest tests: 10 passed.
- Final affected build/configuration checks after fixture updates: 78 passed.
- Treevis legend tests, including named NWKIT regimes with and without shifts:
  passed; real PDF rendering preserves every branch without invented optimum columns.
- Ruff lint/format, changed shell syntax and whitespace checks passed.

The adapter tests exercise AICc/path, alternative BIC and bootstrap/convergence
modes, replicate aggregation, missing coordinates, resumability, transactional
failure, numeric node names, complete bundle publication and topology-safe
branch mapping. Runtime coverage now registers the Python adapter tests instead
of the retired kfl1ou-only R integration script. Build/release triggers,
permissions, native architecture runners and validation jobs are retained.

A complete production-image rebuild, remote Actions execution and SIF execution
were not performed. Production build defaults continue to follow moving upstream
branches; the local NWKIT commit must be published before remote builds can
consume it. The earlier full NWKIT suite result was 3,986 passed, 57 skipped,
with two documented pre-existing failures (see the AICc study validation record).
No NWKIT source changed between that verification and the initial backend commit.

## Comparative studies

A subsequent [100-tip / 10-shift paired benchmark](benchmarks/native-ou-100tips-10shifts/RESULTS.md)
compares kfl1ou and native selection with observed time, peak RSS and recovery
metrics. It is a small development comparison, separate from the adoption study.
Its kfl1ou baseline was subsequently found to fail the pBIC capability probe;
the reported comparison is not evidence against the corrected pBIC backend.

The [seven-method rerun](benchmarks/native-ou-100tips-10shifts-ic/RESULTS.md) adds
NWKIT native pBIC/AIC/BIC and ordinary kfl1ou AIC. All 42 runs completed after
the corrected-pBIC probe passed. Native IC methods were faster in this setup;
accuracy varied by method and dataset. The accompanying validation records
fixed-layout score agreement and identical native candidate searches.

The [AIC candidate-path study](benchmarks/native-ou-aic-improvement/DECISION.md)
then evaluates 30 new seeds per condition plus tree/count scope checks (306 runs,
all complete). The optional `native-path` search improves weak-shift recovery
and approaches kfl1ou AIC in the primary setup. Mixed scope-check results keep
the original `auto` behavior in place. All three AIC procedures select shifts
in every global-null replicate; these results do not establish a calibrated
test for the existence of shifts.

The subsequent [global-null AIC gate pilot](benchmarks/native-ou-global-null-gate/RESULTS.md)
adds full-search plug-in calibration. In 30 null datasets with 19 inner draws,
any-shift selection drops from 30/30 to 3/30, but weak-shift recall also falls
substantially. This optional gate remains experimental: it does not establish
5% false-positive control or control false branches when real shifts exist.

The [AICc addendum](benchmarks/native-ou-aicc-comparison/RESULTS.md) adds kfl1ou
and NWKIT AICc, distinguishes beam and native-path search, and evaluates both
the original nonnull fixtures and the same thirty null fixtures. Native AICc
uses the same independent-trait correction convention as kfl1ou; fixed-alpha
cross-backend score checks pass. Existing defaults remain unchanged.

The [1,000-tip / 100-shift AICc timing study](benchmarks/native-ou-1000tips-100shifts-aicc/RESULTS.md)
compares native-path with beam search under explicit large-search budgets and
estimated covariance. Twelve sequential runs completed. The study reports
elapsed time and actual shift recovery separately; it excludes calibration,
support bootstrap and convergence search.

## Engineering evidence

Validation used the GeneGalleon Docker runtime on Linux aarch64. NWKIT changes
were overlaid on an isolated checkout of the existing NWKIT baseline to exclude
concurrent, unrelated ASR changes. No SIF runtime was available.

| Check | Result |
| --- | --- |
| Independent dense Gaussian covariance, likelihood and GLS comparisons | Passed, including missing observations, known errors and alpha limits |
| NWKIT targeted native, legacy shift and CLI convention tests | 157 passed, 4 skipped |
| Final native tests after the last fitting refactor | 84 passed |
| GeneGalleon adapter, statistics, entrypoint, shell safety and provenance tests | 321 passed |
| Adapter tests after preserving numeric internal-node names | 5 passed |
| Isolated NWKIT lint, format, types, security, dependency audit and maintainability checks | Passed |
| Wheel/sdist content and reproducibility checks | Passed |
| Installed native development image, without source mounts or R on PATH | Generated all eight artifacts |
| Both pages of the installed-image expression PDF | Rendered and visually inspected |

The intermediate full NWKIT suite had 3,851 passes, 57 skips and three failures.
The new CLI underscore-alias failure was fixed and its tests passed. The other
two failures reproduced on the unmodified baseline:
`test_pgls_raw_bootstrap_refits_automatic_gene_parameter` and
`test_archived_engine_requires_explicit_scope_and_intact_snapshot`.
Thus the full suite is not reported as clean. Combined coverage, refreshed for
changed modules, was 86%.

The local image is `local/genegalleon:native-ou-dev`. It installs the native
NWKIT snapshot and places the adapter at
`/opt/gg-native-support/detect_ou_shift_native.py`. This is a development image,
not a published release or evidence of SIF compatibility.

## Computational scope

Fixed-layout four-trait, ten-shift comparisons used the same inputs and fixed
covariance for dense GLS and tree whitening. At 1,000 tips the maximum reported
likelihood/coefficient difference was below `2e-11` for balanced trees and
`2.5e-11` for pectinate trees. These comparisons validate the numerical kernel,
not complete model selection or statistical power.

One configured search with two traits, fixed covariance and a ten-shift cap took
approximately 0.71, 2.46, 6.57 and 29.69 seconds at 32, 128, 512 and 1,000 balanced
tips, respectively. The 1,000-tip run used approximately 264 MB peak RSS.
Calibration and support resampling are excluded. These are development
measurements under concurrent machine load, not a production speedup claim.
One eight-tip comparison each with fixed and estimated covariance had zero
likelihood gap between heuristic and exhaustive search. This does not establish
general heuristic coverage or global continuous optimization.

The reproducible NWKIT entry points are `tools/benchmark_native_shift.py`,
`tools/benchmark_native_search.py` and
`tools/validate_native_shift_development.py`.

## Statistical work still required

Four small development scenarios exercised complete and missing observations,
known sampling errors, a terminal shift and a shared regime. Each scenario had
only one dataset and 19 calibration draws. The complete-data null selected a
shift, the null with known errors/missingness selected no shift, the terminal
shift was recovered, and the shared-regime example selected no shift. These
outcomes are debugging evidence; they cannot estimate false-positive rates,
power or noninferiority. Development failures exposed variance-boundary and
fixed-alpha mean-rank defects, which were corrected and rechecked.

The predeclared NWKIT `NATIVE_SHIFT_VALIDATION.md` protocol requires independent
confirmation datasets, per-stratum false-positive bounds, paired power/location/
grouping comparisons against the production method, and broader performance
and runtime validation. Those complete adoption studies have not run. The user
authorized backend replacement after reviewing the comparative results; research
metadata and the statistical limitations remain visible in the outputs.
