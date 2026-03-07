# Development and Tests

## Install test dependencies

Python-side test dependencies are listed in:

- `workflow/tests/requirements.txt`

Typical setup:

```bash
python -m pip install -r workflow/tests/requirements.txt
```

## Run the full Python test suite

```bash
python -m pytest -q workflow/tests
```

This covers a mix of:

- support-script behavior,
- output summarizers,
- manifest generation,
- trait-generation helpers,
- container build entrypoint behavior,
- path/bootstrap helpers,
- shell static-safety invariants.

## Run focused test subsets

Useful targeted commands:

```bash
python -m pytest -q workflow/tests/test_shell_static_safety.py
```

```bash
python -m pytest -q workflow/tests/test_container_build_entrypoint.py
```

```bash
python -m pytest -q workflow/tests/test_format_species_inputs.py workflow/tests/test_build_download_manifest.py
```

```bash
python -m pytest -q workflow/tests/test_generate_species_trait.py workflow/tests/test_trait_input_templates.py
```

## Run R-side parse checks

```bash
Rscript workflow/tests/test_r_scripts_parse.R
```

This is the quickest way to catch syntax regressions in bundled R helpers without running the full scientific analyses.

## Static guardrails in the test suite

The repository includes tests that intentionally enforce shell hygiene, for example:

- `set -euo pipefail` in core and entrypoint scripts,
- safe `rm -rf --` usage,
- quoted path expansions,
- avoiding fragile `for ... in $(...)` shell patterns.

When editing shell code, `workflow/tests/test_shell_static_safety.py` is often the fastest high-signal check to run first.

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

1. run the most relevant focused pytest file,
2. run `Rscript workflow/tests/test_r_scripts_parse.R` if any R helper changed,
3. run `GG_ENTRYPOINT_DRY_RUN=1 bash workflow/gg_all_entrypoints_debug.sh` for wrapper-level sanity,
4. inspect `workspace/output/versions/` after a real stage run.
