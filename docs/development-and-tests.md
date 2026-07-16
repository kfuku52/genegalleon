# Development and Tests

## Container-first validation policy

For workflow validation, integration tests, R helper checks, and toolchain-dependent behavior, use the repository `genegalleon.sif` runtime instead of host-local Python, R, or command-line tools.

The preferred wrapper is:

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

`workflow/tests/run_in_sif.sh` expects `./genegalleon.sif` at the repository root. Set `GENEGALLEON_SIF=/path/to/genegalleon.sif` when using a different SIF path.

Host-local checks are useful for quick syntax or narrow static feedback, but do not use them as evidence of container runtime compatibility. If the SIF or Apptainer/Singularity is unavailable, report that limitation explicitly.

## Install test dependencies

Python-side test dependencies are listed in:

- `workflow/tests/requirements.txt`

Typical setup:

```bash
python -m pip install -r workflow/tests/requirements.txt
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

CI also runs `actionlint` against parsed GitHub Actions workflows and runs
ShellCheck at warning severity against the self-contained core drivers and the
`gg_util.sh` compatibility façade. Local equivalents are:

```bash
actionlint
```

```bash
shellcheck -S warning -x \
  -e SC1091,SC2034,SC2154,SC2317 \
  workflow/core/gg_gene_evolution_core.sh \
  workflow/core/gg_genome_evolution_core.sh \
  workflow/core/gg_transcriptome_generation_core.sh \
  workflow/support/gg_util.sh \
  workflow/support/gg_util/*.sh
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
