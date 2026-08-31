# Development and Tests

## Container-first validation policy

For workflow validation, integration tests, R helper checks, and
toolchain-dependent behavior, use a GeneGalleon container instead of
host-local Python, R, or command-line tools.

```bash
bash ./dev check fast
bash ./dev check runtime
bash ./dev check r
bash ./dev check full
```

`workflow/tests/run_in_runtime.sh`, used by `dev`, prefers a usable SIF on
Linux/HPC and otherwise uses `local/genegalleon:dev` through Docker. Force the
selection with `GG_TEST_RUNTIME=sif` or `GG_TEST_RUNTIME=docker`. The command
reports an error rather than treating host-local behavior as runtime evidence.
Before dispatch, it also compares the runtime's
`io.genegalleon.runtime-input` label with the current container inputs and the
repository-owned moving upstream revisions. Upstream resolution is cached for one UTC day, so
normal focused checks do not repeatedly query every repository. A mismatch
fails with the rebuild command instead of silently using an old `nwkit`,
`csubst`, or other source snapshot.

Use `GG_RUNTIME_FRESHNESS=always` to force a new upstream resolution or
`GG_RUNTIME_FRESHNESS=off` for an intentional offline check with a known older
runtime. The latter is an explicit escape hatch and is not compatibility
evidence. BUSCO and PAML remain accepted at the revisions embedded in the
runtime by default; set `GG_RUNTIME_FRESHNESS_SCOPE=all` to compare their
branch tips as well. Scheduled publishing and CI runtime-cache keys always
resolve all sources exactly.

On Linux/HPC hosts with Apptainer or Singularity, the preferred wrapper is:

```bash
bash workflow/tests/run_in_sif.sh python -m pytest -q workflow/tests/test_hgt_end_to_end.py
```

```bash
bash workflow/tests/run_in_sif.sh Rscript workflow/tests/test_treevis_main.R
```

The tree visualization helpers are an installed internal R package. Validate
its package metadata, namespace, installability, and exported API before the R
behavior test:

```bash
bash workflow/tests/run_in_sif.sh bash workflow/tests/check_treevis_package.sh
```

`workflow/tests/run_in_sif.sh` expects `./genegalleon.sif` at the repository
root. Set `GENEGALLEON_SIF=/path/to/genegalleon.sif` when using a different SIF
path. Both validation wrappers and `dev` discover versioned
Apptainer/Singularity installations under `/opt/pkg` and the legacy NIG
package path, using the same discovery helper as workflow entrypoints.

To expose an external read-only fixture or data directory at the same absolute
path inside the SIF, provide one absolute path per line:

```bash
GENEGALLEON_SIF_EXTRA_BINDS=$'/data/reference\n/shared/fixtures' \
bash workflow/tests/run_in_sif.sh python -m pytest -q workflow/tests/test_hgt_end_to_end.py
```

Relative paths, paths containing `:`, and nonexistent paths are rejected.

On macOS, build an image from the current checkout and run validation through
the same wrapper:

```bash
bash ./dev build
```

```bash
GG_TEST_RUNTIME=docker bash workflow/tests/run_in_runtime.sh \
  python -m pytest -q workflow/tests/test_hgt_end_to_end.py
```

`dev build` selects `BUILD_TARGET=development`, adding APT compilers and
headers to the distribution runtime. Direct container build wrappers and CI
default to `BUILD_TARGET=runtime`; either target can run these tests. Their
fingerprints include the target so switching profiles cannot silently reuse
the other image. Use separate tags when retaining both locally.
See [container development](container-development.md) for the independent
source stages and build-cache tradeoffs. Workflow-only edits continue to use
the mounted checkout and do not require an image rebuild.

Host-local checks are useful for quick syntax or narrow static feedback, but do
not use them as evidence of container runtime compatibility. Report validation
as SIF or Docker/container validation according to the runtime actually used.

## Install test dependencies

Python-side test dependencies are listed in:

- `workflow/tests/requirements.txt`
- `workflow/tests/requirements-smoke.txt` (minimal fail-fast CI lane)
- `workflow/tests/requirements.lock.txt` (validated exact constraints for stable
  third-party dependencies; `csubst` intentionally follows its moving `master`
  branch to match the container runtime)

Typical setup:

```bash
python -m pip install \
  --constraint workflow/tests/requirements.lock.txt \
  --requirement workflow/tests/requirements.txt
```

To keep `.pyc` files out of the repository during direct local Python runs,
set a temporary pycache prefix first:

```bash
export PYTHONPYCACHEPREFIX="${TMPDIR:-/tmp}/genegalleon_pycache"
```

GeneGalleon entrypoints and core bootstraps now do this automatically for wrapper-driven runs.
The validation wrappers instead set `PYTHONDONTWRITEBYTECODE=1`, so their
imports do not write bytecode into the mounted checkout.

CI resolves the moving `csubst` branch once for its Python wheel preparation
job. Fast and integration lanes install the resulting shared artifact with
`--no-index`; they do not rebuild the same VCS dependency independently. The
wheel cache includes the resolved source, Python/platform identity, requirements,
constraints, and preparation script. Only trusted default-branch runs save this
cache, without fallback cache keys; run artifacts expire after one day. Resolved
commits remain temporary build metadata, never repository defaults.

## Run all Python and R checks

```bash
bash ./dev check full
```

`workflow/tests/validation_manifest.json` declares the required runtime Python
files, extra runtime scenarios, integration environment, and all seven R
commands. `run_checks.py` executes this same list for `dev`, Docker image
publication, release SIFs, and the CI SIF job. `dev check full` runs all Python
tests and every declared R check; `dev check r` runs the complete R list.
`dev check runtime` also enables the real kfFractBias JCVI/LAST integrations.
Required Python checks fail the validation if they are skipped unexpectedly.
The manifest currently permits no skips, and the image includes system `unzip`
for archive interoperability checks.

To inspect the commands without starting a container:

```bash
python3 workflow/tests/run_checks.py runtime --list
```

For just the full Python suite, with the same mandatory integration policy:

```bash
bash workflow/tests/run_in_sif.sh python -m pytest -q --gg-strict-runtime workflow/tests
```

This covers a mix of:

- support-script behavior,
- output summarizers,
- manifest generation,
- trait-generation helpers,
- container build entrypoint behavior,
- path/bootstrap helpers,
- shell static-safety invariants.

Tests set `GG_INPUT_GENERATION_OUTPUT_ROOT` to a per-test temporary directory,
so CLI defaults cannot leave generated summaries or downloads in the repository
working tree.

## Run focused test subsets

Useful targeted commands:

```bash
bash ./dev check fast workflow/tests/test_shared_namespace_lock.py -k release -x
```

Arguments after the suite are forwarded to pytest, including explicit paths,
`-k`, `-x`, and `--lf`; no extra test-directory argument expands an explicitly
selected file back to the full suite. Runtime/full commands still run their
declared R checks; use a focused Python lane when only Python feedback is needed.

```bash
bash workflow/tests/run_in_sif.sh python -m pytest -q workflow/tests/test_shell_static_safety.py
```

```bash
bash workflow/tests/run_in_sif.sh python -m pytest -q workflow/tests/test_container_build_entrypoint.py
```

```bash
bash workflow/tests/run_in_sif.sh python -m pytest -q workflow/tests/test_format_species_inputs.py workflow/tests/test_build_download_manifest.py
```

```bash
bash workflow/tests/run_in_sif.sh python -m pytest -q workflow/tests/test_generate_species_trait.py workflow/tests/test_trait_input_templates.py
```

The Python suite is also divided into non-overlapping CI lanes. The same
selection works locally and in a GeneGalleon container:

```bash
bash workflow/tests/run_in_sif.sh python -m pytest -q --gg-suite fast workflow/tests
```

```bash
bash workflow/tests/run_in_sif.sh python -m pytest -q --gg-suite integration-download workflow/tests
bash workflow/tests/run_in_sif.sh python -m pytest -q --gg-suite integration-workflow workflow/tests
```

```bash
bash workflow/tests/run_in_sif.sh python -m pytest -q --gg-suite static workflow/tests
bash workflow/tests/run_in_sif.sh python -m pytest -q --gg-suite runtime workflow/tests
```

Use `--gg-suite smoke` for the minimal CI preflight. Suite membership is
defined in `workflow/tests/conftest.py`, with runtime files coming from the
shared validation manifest; pytest markers expose the same lane
metadata during full-suite collection. Strict marker and configuration checks
are enabled, so new marker names must be registered in `pyproject.toml`.
The fast and integration CI lanes use two pytest-xdist workers. `bash ./dev`
uses the same default; override it with `GG_PYTEST_WORKERS=N` or
`GG_PYTEST_WORKERS=auto`.

The runtime lane includes real integration contracts for repository-owned
upstreams. In particular, core-species selection executes the installed
`nwkit sample` command rather than a mock, and the test verifies that the
container records exact revisions for every moving source. The daily container
publisher runs this lane against each newly built architecture, so an
upstream-only change is checked without requiring a GeneGalleon commit.

CI caches a validated SIF by an exact runtime-content hash. The key includes
the platform, Singularity version, daily security epoch, all resolved upstream commits, and every
file copied into the container, but excludes GeneGalleon revision/version
metadata because the checkout is mounted at test time. A cache miss may reuse
the published daily image only when its runtime label matches exactly;
otherwise CI builds the current container before conversion. The cache is
saved only after all SIF validation checks succeed on the default branch.
The daily publisher primes that same cache from its immutable multi-architecture
image, passing the already resolved runtime identity and security epoch. It
does not resolve moving branches a second time during SIF preparation. Each
commit still validates its mounted workflow source, even on a cache hit.
Version-only and documentation-only changes do not launch the heavy SIF job.
After the SIF's embedded identity is checked against that exact input, the
validation wrapper does not re-query branch tips; a branch moving midway
through the run cannot replace the snapshot being validated. A new security
epoch or different source revision still requires a different cache entry.

## Run R-side parse checks

```bash
bash workflow/tests/run_in_sif.sh Rscript workflow/tests/test_r_scripts_parse.R
```

This is the quickest way to catch syntax regressions in bundled R helpers without running the full scientific analyses.

## Static guardrails in the test suite

The repository includes tests that intentionally enforce shell hygiene, for example:

- `set -euo pipefail` in core and entrypoint scripts,
- safe `rm -rf --` usage,
- quoted path expansions,
- avoiding fragile `for ... in $(...)` shell patterns.

When editing shell code, `workflow/tests/test_shell_static_safety.py` is often the fastest high-signal check to run first.

The editable entrypoint blocks are also the source for a machine-readable
configuration schema. Validate the forwarding registry or render current
JSON/Markdown reference data with:

```bash
bash ./dev config-check
bash ./dev config-schema json
bash ./dev config-schema markdown
```

Use `bash ./dev bump patch` (or `minor`/`major`) for an atomic Semantic Version
update; `python workflow/support/bump_version.py patch --dry-run` previews it.

`gg_genome_annotation_core.sh` and `gg_transcriptome_generation_core.sh`
support `GG_CORE_SOURCE_ONLY=1` when a
focused shell test needs to load its helper functions without starting a
workflow. Its functions and ordered stages remain in the self-contained core
file; the switch is a test boundary, not a sourced production architecture.
FASTQ recovery tests load the real transcriptome helper this way and exercise
connection interruption, process restart, changed objects, and invalid HTTP
ranges against a loopback server.

For immediate host-only syntax, Ruff, and configuration checks:

```bash
bash ./dev lint
```

This needs Python 3 and Ruff on the host and does not establish runtime
compatibility. Use the container checks above before publishing changes.

CI also runs `actionlint` against parsed GitHub Actions workflows and runs
ShellCheck at warning severity against every tracked shell script. Local
equivalents are:

```bash
actionlint
python3 workflow/tests/check_composite_actions.py
```

```bash
mapfile -d '' -t scripts < <(git ls-files -z '*.sh')
scripts+=(dev)
shellcheck -S warning -x -e SC1091,SC2034,SC2154,SC2317 "${scripts[@]}"
```

## Database memory and timing comparisons

`benchmark_orthogroup_database.py` creates deterministic many-family and
single-large-file workloads (131,072 branch rows each), runs each build in a
fresh process, and records wall time, Linux peak RSS, row counts, and a digest
of sorted SQLite contents and schema. Run both revisions in the same container,
without concurrent builds or tests, and compare the logical digests before
interpreting performance differences:

```bash
bash workflow/tests/run_in_runtime.sh python workflow/tests/benchmark_orthogroup_database.py \
  --output /tmp/database-current.json
```

Use `--source /path/to/generate_orthogroup_database.py` to measure a saved
baseline script; ensure that path is mounted in the runtime. These synthetic
workloads exercise bounded input buffering, not the peak memory of the global
BH-FDR calculation or a production workload of arbitrary size.

## Dependency-aware debug harness

`workflow/gg_all_entrypoints_debug.sh` runs all major entrypoints in a dependency-aware order and records a summary TSV.

Useful modes:

```bash
GG_ENTRYPOINT_DRY_RUN=1 bash workflow/gg_all_entrypoints_debug.sh
```

```bash
GG_ENTRYPOINT_ONLY_STEPS=gg_genome_evolution,gg_progress_summary \
bash workflow/gg_all_entrypoints_debug.sh
```

```bash
GG_ENTRYPOINT_TIMEOUT_SEC=0 bash workflow/gg_all_entrypoints_debug.sh
```

Key environment variables:

- `GG_DEBUG_LOG_DIR`: override log directory
- `GG_ENTRYPOINT_TIMEOUT_SEC`: per-step timeout in seconds
- `GG_ENTRYPOINT_DRY_RUN=1`: print commands without executing them
- `GG_ENTRYPOINT_ONLY_STEPS`: comma-separated subset of step IDs
- `GG_BENCHMARK=1`: collect timing metrics
- `GG_BENCHMARK_RAW=1`: keep raw `/usr/bin/time` outputs

Artifacts are written under:

- `workspace/output/debug_entrypoint_logs/summary.tsv`

## Building a smaller development dataset

`workflow/support/build_minimal_test_dataset.py` extracts a compact analysis-ready subset from an existing workspace.

Typical use:

```bash
python workflow/support/build_minimal_test_dataset.py --help
```

This is useful when:

- reproducing a bug with a smaller dataset,
- creating a local smoke-test workspace,
- trimming a large project down to a manageable debug case.

## Practical validation loop

A good local validation sequence after code changes is:

1. run the most relevant focused pytest file through `workflow/tests/run_in_sif.sh`,
2. run `bash workflow/tests/run_in_sif.sh Rscript workflow/tests/test_r_scripts_parse.R` if any R helper changed,
3. run `GG_ENTRYPOINT_DRY_RUN=1 bash workflow/gg_all_entrypoints_debug.sh` for wrapper-level sanity,
4. inspect `workspace/output/versions/` after a real stage run.
