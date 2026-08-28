<!-- BEGIN KF AGENT POLICY: source=https://github.com/kfuku52/kf-agent-policy; version=6; sha256=9a1279ed25a1782ca13ed3a1bebc74d389c35e34a52f8a8ec762bb7790d05186 -->
# Common agent policy

Repository-specific instructions override these defaults.

- Before edits, inspect the worktree and preserve unrelated user changes. When
  remote state matters, update from the default branch without discarding local
  work.
- Use the default branch unless the user explicitly requests another existing
  one. Never create or switch branches solely for a commit, push, release, or
  pull request.
- Change or recommend branch protection only when explicitly asked. If it
  blocks a requested direct push, report it; never bypass it or create a branch
  or pull request.
- In library metadata, exact pins or upper bounds require demonstrated
  incompatibility. Keep reproducibility locks separate; prefer fixing and
  testing compatibility.
- Interface, option, format, filename, or schema changes must update all
  producers, consumers, tests, examples, and documentation.
- Keep each repository's top-level README concise: do not add feature-specific
  guides or extended examples. Put them in dedicated documentation or the wiki,
  linking from the README only when needed for discoverability.
- Changes confined to unpushed local commits need no backward compatibility.
- Prefer verified root-cause fixes to fallbacks or relaxed validation that only
  hide failures. Document unavoidable workarounds and their removal conditions.
- When changing GitHub Actions, preserve required coverage while minimizing
  runner use: avoid duplicate push/PR runs and high-frequency polling, cancel
  superseded CI, combine short jobs when isolation is unnecessary, and retain
  artifacts only as long as needed.
- Use platform-specific runners only when their coverage is required, gate them
  by relevant paths, schedules, releases, or manual dispatch, and never execute
  untrusted pull-request code on self-hosted runners.
- Run focused checks and, when practical, the standard suite. Directly exercise
  affected behavior or rendered artifacts; report exactly what did and did not
  run.
- For performance work, benchmark representative workloads before and after,
  verify equivalent output, and report wall time and peak memory when relevant.
- Individual local commits need no version bump. Before GitHub pushes, bump the
  version even if unrequested, using the repository's scheme or Semantic
  Versioning (`MAJOR.MINOR.PATCH`) if absent.
<!-- END KF AGENT POLICY -->

# Agent / Developer Validation Policy

## Upstream program version policy

Do not commit default versions, tags, or commit SHAs for upstream programs.
Container source defaults must follow the moving branches declared in
`container/source_branches.env`. Build wrappers may resolve those branches to
commits in memory so all platforms in one build use the same snapshot, and
callers may provide a `*_REPO_SHA` for a one-off reproduction or debugging
build, but resolved SHAs must never be copied back into repository defaults.

Cryptographic hashes used only to verify downloaded artifacts, digest-pinned
base images and GitHub Actions, and compatibility constraints with a documented
demonstrated incompatibility are outside this rule.

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
