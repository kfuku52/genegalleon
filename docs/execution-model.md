# Execution Model

GeneGalleon wrappers are designed so the same `gg_*_entrypoint.sh` can be run:

- directly with `bash`,
- through SLURM,
- through UGE,
- or through PBS.

The wrapper normalizes host/scheduler/runtime details first, then streams the matching
`workflow/core/gg_*_core.sh` into the container shell.

## High-level flow

Typical execution path:

1. Launch `workflow/gg_*_entrypoint.sh`.
2. The wrapper locates `workflow/support/gg_entrypoint_bootstrap.sh`.
3. Bootstrap resolves the workflow root and default paths such as `workspace/` and `genegalleon.sif`.
4. The wrapper loads `workflow/support/gg_util.sh`, chooses `apptainer` or `singularity` (or a Docker-backed shim when requested), and normalizes scheduler variables.
5. The workspace is bind-mounted to `/workspace` and the workflow tree is bind-mounted to `/script`.
6. The wrapper changes into the workspace root and pipes the matching core script into the container shell.
7. The core script bootstraps runtime helpers from `/script/support`, resolves `workspace/input`, `workspace/output`, and `workspace/downloads`, then runs stage logic.

## Host-side bootstrap

Each entrypoint searches for `gg_entrypoint_bootstrap.sh` in this order:

1. `<submit-dir>/support/gg_entrypoint_bootstrap.sh`
2. `<submit-dir>/workflow/support/gg_entrypoint_bootstrap.sh`
3. `<entrypoint-dir>/support/gg_entrypoint_bootstrap.sh`

This makes the wrappers resilient to scheduler spool directories and direct execution from either the repo root or the `workflow/` directory.

Bootstrap then loads:

- `workflow/gg_path_defaults.sh`
- `workflow/gg_common_params.sh` for entrypoints that opt into shared defaults

Default path assumptions are:

- workflow root: auto-detected from the entrypoint location, or from the sourced bootstrap path when the scheduler rewrites the entrypoint to a spool file such as Slurm `slurm_script`
- workspace root: `workflow/../workspace`
- container image: `workflow/../genegalleon.sif`

Advanced path overrides that can be exported before launch:

- `gg_workspace_dir=/path/to/workspace`
- `gg_container_image_path=/path/to/genegalleon.sif`
- `GG_CONTAINER_RUNTIME=docker`
- `GG_CONTAINER_DOCKER_IMAGE=ghcr.io/kfuku52/genegalleon:latest`

## Scheduler normalization

Before the container starts, `gg_entrypoint_prepare_container_runtime`:

- runs a site-specific scheduler prelude from `workflow/support/gg_site_runtime.sh`,
- clears stale `SINGULARITY_*` / `APPTAINER_*` bind variables,
- optionally performs duplicate-job checks for wrappers that enable `exit_if_running`,
- detects `apptainer` or `singularity`,
- normalizes scheduler metadata so downstream scripts can rely on scheduler-neutral `GG_*` variables.

Inside the runtime, GeneGalleon consistently uses:

- `GG_TASK_CPUS`
- `GG_JOB_ID`
- `GG_ARRAY_TASK_ID`

Those are populated from SLURM/PBS/UGE values when available and then forwarded into the container as both `SINGULARITYENV_*` and `APPTAINERENV_*`.

## Container activation and bind mounts

`gg_entrypoint_activate_container_runtime` performs the second half of setup:

- resolves physical paths on non-macOS hosts,
- binds the workspace root to `/workspace`,
- binds the workflow root to `/script`,
- forwards shared `GG_COMMON_*` variables,
- forwards entrypoint-specific config variables that are registered in `workflow/support/gg_entrypoint_config_vars.sh`,
- prints a scheduler/runtime summary for debugging.

Site adapters in `workflow/support/gg_site_runtime.sh` can also add extra bind mounts or runtime flags. For example, site profiles may:

- bind scheduler spool directories,
- inject site-specific `PATH` entries,
- change the container shell command to include options such as `--contain`.

## Core script bootstrap

Inside the container, each core script loads `workflow/support/gg_core_bootstrap.sh`, which:

- resolves the support directory from `/script/support` or a local fallback,
- loads `workflow/gg_common_params.sh` when present,
- sources `workflow/support/gg_util.sh`,
- initializes the split workspace layout,
- prepares cache directories such as taxonomy DB locations under `workspace/downloads/`.

Core scripts therefore assume the following mount layout:

- `/workspace` for project data
- `/script` for the checked-out workflow code

## Where to customize behavior

For routine use, change only the top config block in each entrypoint:

```bash
### Start: Modify this block to tailor your analysis ###
...
### End: Modify this block to tailor your analysis ###
```

That block is the supported place for:

- enabling or disabling substeps,
- selecting modes,
- changing stage-specific thresholds,
- switching tools such as `blastp` vs `mmseqs2`.

`workflow/core/gg_*_core.sh` files are the implementation layer and are not the primary place for routine parameter tuning.

## What gets logged

During startup, wrappers print:

- the detected scheduler kind,
- normalized slot/task/job metadata,
- the chosen container runtime,
- bind mounts and forwarded environment counts.

On successful completion, wrappers also trigger a versions dump under:

- `workspace/output/versions/`

Those logs are often the fastest way to understand what runtime state a job actually saw.
