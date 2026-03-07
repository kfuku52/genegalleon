# Repository Layout

This page is a navigation map for the repository. When in doubt, start here and then jump to the more detailed topic pages.

## Top-level structure

- `README.md`
  - project entry page and documentation index
- `LICENSE`
  - license text
- `VERSION`
  - repository version marker
- `gg_container_build_entrypoint.sh`
  - top-level convenience wrapper for container image/SIF creation
- `docs/`
  - topic-specific user documentation
- `container/`
  - container build scripts, Dockerfile, runtime specs, and container-focused notes
- `workflow/`
  - entrypoints, core pipeline scripts, helpers, and tests
- `workspace/`
  - default project-local data root used by the wrappers
- `logo/`
  - project logo assets

## `workflow/`

This is the operational center of the repository.

Main contents:

- `workflow/gg_*_entrypoint.sh`
  - scheduler-friendly wrappers that users launch directly
- `workflow/core/gg_*_core.sh`
  - implementation scripts streamed into the container runtime
- `workflow/support/`
  - shared shell/Python/R helpers, plotting utilities, bootstrap scripts, and support tools
- `workflow/tests/`
  - pytest and R-based checks for key helpers and wrapper behavior
- `workflow/gg_common_params.sh`
  - shared defaults reused by multiple stages
- `workflow/gg_path_defaults.sh`
  - default path definitions for workspace and container image
- `workflow/gg_all_entrypoints_debug.sh`
  - dependency-aware debug runner across the main entrypoints

If you want to understand how a stage is launched, start with the entrypoint.
If you want to understand what it actually does, then read the matching core script.

## `container/`

This directory contains the container build system and related runtime checks.

Important files:

- `container/Dockerfile`
  - canonical container definition
- `container/buildx.sh`
  - Docker Buildx build helper
- `container/apptainer_from_docker.sh`
  - convert/pull a Docker image into a SIF
- `container/apptainer_local_build.sh`
  - native local SIF build helper when Docker is unavailable
- `container/gg_container_build_impl.sh`
  - implementation behind `gg_container_build_entrypoint.sh`
- `container/CAPABILITY_MATRIX.md`
  - command/runtime availability notes
- `container/scripts/`
  - install, validation, and post-build helper scripts
- `container/env/`
  - package requirement lists for the `base` runtime environment

## `docs/`

`README.md` stays intentionally short.
Detailed operational notes live under `docs/`, split by topic such as:

- execution model,
- common workflow recipes,
- container build/runtime behavior,
- input conventions,
- troubleshooting,
- development/testing.

## `workspace/`

This is the default persistent data root used by the wrappers.

High-level split:

- `workspace/input/`
  - curated inputs and manifests
- `workspace/output/`
  - stage outputs and reports
- `workspace/downloads/`
  - reusable caches, downloaded databases, and lock-managed runtime assets

See `docs/workspace-layout-and-data-model.md` for the detailed data model.

## How to navigate efficiently

Use this rough rule:

- user-facing execution behavior: `workflow/gg_*_entrypoint.sh`
- stage logic: `workflow/core/gg_*_core.sh`
- shared runtime/path behavior: `workflow/support/gg_util.sh`
- bootstrap/path defaults: `workflow/support/gg_entrypoint_bootstrap.sh`, `workflow/gg_path_defaults.sh`
- container behavior: `container/`
- usage guidance: `docs/`
