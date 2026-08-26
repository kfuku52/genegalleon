# Development and Tests

## Container-first validation policy

For workflow validation, integration tests, R helper checks, and
toolchain-dependent behavior, use a GeneGalleon container instead of
host-local Python, R, or command-line tools.

```bash
bash ./dev check fast
bash ./dev check runtime
bash ./dev check r
```

`workflow/tests/run_in_runtime.sh`, used by `dev`, prefers a usable SIF on
Linux/HPC and otherwise uses `local/genegalleon:dev` through Docker. Force the
selection with `GG_TEST_RUNTIME=sif` or `GG_TEST_RUNTIME=docker`. The command
reports an error rather than treating host-local behavior as runtime evidence.

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
path. The wrapper also discovers versioned Apptainer/Singularity installations
under `/opt/pkg` and the legacy NIG package path.

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
BUILD_SIF=0 IMAGE_SOURCE=local IMAGE=local/genegalleon TAG=dev \
bash ./gg_container_build_entrypoint.sh
```

```bash
GG_TEST_RUNTIME=docker bash workflow/tests/run_in_runtime.sh \
  python -m pytest -q workflow/tests/test_hgt_end_to_end.py
```

Host-local checks are useful for quick syntax or narrow static feedback, but do
not use them as evidence of container runtime compatibility. Report validation
as SIF or Docker/container validation according to the runtime actually used.

## Install test dependencies

Python-side test dependencies are listed in:

- `workflow/tests/requirements.txt`
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

## Run the full Python test suite

```bash
bash workflow/tests/run_in_sif.sh python -m pytest -q workflow/tests
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
defined in `workflow/tests/conftest.py`; pytest markers expose the same lane
metadata during full-suite collection. Strict marker and configuration checks
are enabled, so new marker names must be registered in `pyproject.toml`.

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

`gg_genome_annotation_core.sh` also supports `GG_CORE_SOURCE_ONLY=1` when a
focused shell test needs to load its helper functions without starting a
workflow. Its functions and ordered stages remain in the self-contained core
file; the switch is a test boundary, not a sourced production architecture.

CI also runs `actionlint` against parsed GitHub Actions workflows and runs
ShellCheck at warning severity against every tracked shell script. Local
equivalents are:

```bash
actionlint
```

```bash
mapfile -d '' -t scripts < <(git ls-files -z '*.sh')
scripts+=(dev)
shellcheck -S warning -x -e SC1091,SC2034,SC2154,SC2317 "${scripts[@]}"
```

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
