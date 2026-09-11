# Publication validation — 2026-09-11

CSUBST 1.16.3 is published at
[`7c5c3c2`](https://github.com/kfuku52/csubst/commit/7c5c3c28d21382d773b25257ebb9546335fecf09).
Its [complete CI run](https://github.com/kfuku52/csubst/actions/runs/34593356320)
passed, including native/fallback tests, packaging, numerical/performance parity,
sanitizers and Python 3.10–3.14 checks. The earlier sdist-support and stale
reference failures, already present before the 3Di change, were repaired in
CSUBST. Its validation record documents the scientific controls and separately
collected Linux performance baseline. The IQ-TREE NaN limitation was not fixed.

## ARM64 image

A dedicated empty BuildKit builder built the entire runtime from the current
Dockerfile and dependency definitions. This included OS/Conda/pip installation,
upstream source builds and the image's runtime checks. The first complete
build used CSUBST 1.16.1; after the CI-only CSUBST fixes were published, the same
builder refreshed the source wheel to 1.16.3 and reran runtime validation.
The final image is `local/genegalleon:3di-clean-20260911`; its exact digest and
labels are in `arm64-image.json`. The mounted workflow implementation is
GeneGalleon commit `8383ecc`, version 0.7.106.

The build resolved moving upstream branches in memory; no source commit was
written into a repository default. A temporary anonymous Docker configuration
was used for public registry pulls because this host's credential helper hung.
The user's Docker credentials were unchanged.

The final image passed both real-predictor cases with the workflow's codon
model defaults, genetic codes 1 and 2, full/trimmed coordinate differences,
gap/ambiguity handling, search/scan direct-call parity, both sites routes,
R/PDF reports, relocated bundles and offline cache reuse (`arm64-3di.log`).
It used the previously verified model resources, without installing or patching
any runtime package after the image build.

The published-source lint checks passed (`source-checks.log`), and the SIF
CI interface tests passed 36 cases (`ci-contract-tests.log`). The GitHub tests
workflow now explicitly opts into the real-predictor tests inside its AMD64 SIF
job; other users of the shared SIF action retain the default opt-out.

The complete canonical runtime suite passed: 254 runtime Python cases,
2 additional integration cases, and all 16 R validation commands. The treevis
package's `R CMD check` reported `Status: OK` (`arm64-runtime.log`). The static
lane passed 397 cases inside the same final image (`arm64-static.log`).

Commands in the final image:

```bash
python workflow/tests/run_checks.py runtime -p no:cacheprovider \
  --basetemp=/tmp/gg-final-runtime
GG_TEST_CSUBST_3DI=1 python -m pytest -q -p no:cacheprovider \
  --gg-strict-runtime workflow/tests/test_csubst_3di_runtime.py
python -m pytest -q --gg-suite static -p no:cacheprovider workflow/tests
```

`arm64-build-initial.log` records the empty-builder build, and
`arm64-build-final.log` records the final source refresh and runtime checks.
`resolved-sources.json` preserves the final build's source snapshot as evidence,
not as repository defaults. Model resources were mounted from the earlier
validation cache; the canonical runtime suite's temporary files used the
container filesystem.

Native SIF execution is unavailable on this macOS host. The publication's
GitHub tests run supplies the separate AMD64/SIF execution evidence.

## Python CI dependency correction

The first publication run, [34595419636](https://github.com/kfuku52/genegalleon/actions/runs/34595419636),
exposed a missing NWKIT dependency in the ordinary Python test environment.
The production container already included NWKIT, but the test wheel requirements
listed only CSUBST. Imports of `nwkit.file_paths` consequently failed in both
the fast and workflow-integration lanes.

Version 0.7.107 adds NWKIT's moving master branch to the test requirements and
resolves both owned dependencies once per run. Temporary wheel-build inputs and
the cache identity include both resolved commits; repository defaults still
follow moving branches. Test coverage and runtime code are unchanged.

The affected tests and CI contracts passed 232 cases in the final GeneGalleon
ARM64 image with the corrected source (`python-ci-regression.log`). Ruff and
Actionlint also passed for the changed files.

The complete test wheel set also built successfully in a separate container
from the same GeneGalleon image with build-essential installed for this check.
A fresh Python virtual environment installed the set with `--no-index`, passed
`pip check`, and imported `nwkit.file_paths` (`python-ci-wheels.log`). The
production image was unchanged. This local wheel check used ARM64; hosted CI
provides the separate AMD64 wheel and test-lane verification.


## AMD64 build correction

The first SIF job stopped before conversion: IQ-TREE compiled, but NWKIT's
external worker did not link the system zlib selected by CMake, causing
undefined `gzread`/`gzopen` symbols. NWKIT's own build code is corrected in
0.43.17; GeneGalleon continues to consume its moving source branch. Details and
the original errors are recorded in the
[NWKIT review](https://github.com/kfuku52/nwkit/tree/master/reviews/system-zlib-worker-2026-09-11).
The initial failed job is not SIF execution evidence.
