<!-- BEGIN KF AGENT POLICY: source=https://github.com/kfuku52/kf-agent-policy; version=2; sha256=9361b6d2342940b167ec77884d09c9ae337e5964b9fc94a80386fac5f1fa7e95 -->
# Common agent policy

Repository-specific instructions override these defaults.

- Work on the default branch unless the user explicitly requests another
  existing branch. Do not create or switch branches merely to commit, push,
  release, or open a pull request.
- Do not modify or recommend branch protection unless explicitly asked. If it
  blocks a requested direct push, report the blocker instead of bypassing it or
  creating a branch or pull request.
- In library metadata, use exact pins or upper bounds only for demonstrated
  incompatibility. Treat reproducibility locks separately, and prefer fixing
  and testing compatibility.
- Before adding or removing a direct dependency, confirm direct use in code,
  configuration, tests, or documentation. Validate removals in a clean
  environment.
<!-- END KF AGENT POLICY -->

# Agent / Developer Validation Policy

Use a GeneGalleon container runtime for workflow validation, integration tests, R helper checks, and toolchain-dependent behavior.

On Linux/HPC hosts with Apptainer or Singularity, prefer the repository `genegalleon.sif` runtime.

On macOS, where SIF execution is normally unavailable, use the Docker-backed GeneGalleon runtime instead. For validation of local code changes, prefer a Docker image built from the current repository rather than a stale public image:

```bash
BUILD_SIF=0 IMAGE_SOURCE=local IMAGE=local/genegalleon TAG=dev bash ./gg_container_build_entrypoint.sh
docker run --rm -i -v "$PWD:$PWD" -w "$PWD" local/genegalleon:dev python -m pytest -q workflow/tests/test_hgt_end_to_end.py
docker run --rm -i -v "$PWD:$PWD" -w "$PWD" local/genegalleon:dev Rscript workflow/tests/test_treevis_main.R
```

GeneGalleon entrypoint wrappers can also dispatch through the Docker-backed singularity shim:

```bash
GG_CONTAINER_RUNTIME=docker \
GG_CONTAINER_DOCKER_IMAGE=local/genegalleon:dev \
bash workflow/gg_progress_summary_entrypoint.sh
```

Do not treat host-local Python, R, or command-line tool behavior as authoritative for GeneGalleon runtime compatibility. Host-local checks are acceptable for quick syntax or narrow static checks, but container-backed checks are the source of truth when dependencies matter.

Preferred SIF validation entrypoint:

```bash
bash workflow/tests/run_in_sif.sh python -m pytest -q workflow/tests/test_hgt_end_to_end.py
bash workflow/tests/run_in_sif.sh Rscript workflow/tests/test_treevis_main.R
```

If `genegalleon.sif` or an Apptainer/Singularity runtime is unavailable, report that clearly and do not conclude SIF compatibility from host-local results. If Docker is available, use Docker-backed container validation and report it as Docker/container validation rather than SIF validation. If no GeneGalleon container runtime is available, report that clearly and do not conclude runtime compatibility from host-local results.

Do not add backward-compatibility workarounds for older dependency behavior.

If the root cause is in a dependency program or package, do not patch GeneGalleon to absorb it. Report it as dependency-side so the dependency can be fixed or updated instead.

# Core Workflow Architecture

Keep `workflow/core/gg_*_core.sh` as self-contained workflow implementation scripts. Do not split their functions or ordered execution stages into `workflow/core/stages/`, per-stage shell fragments, or sourced function libraries merely to reduce file size or reorganize code.

Changes to this architecture require explicit user approval. Without that approval, edit the matching core script in place and preserve its existing execution order.
