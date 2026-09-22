<!-- BEGIN KF AGENT POLICY: source=https://github.com/kfuku52/kf-agent-policy; version=10; sha256=82e3c0eb467582a414d9a6b2feaaaf6f5c8ae330d30f2e3efbf8c303155d0e2e -->
# Common agent policy

Repository-specific instructions override these defaults.

- Follow the user's task scope within higher-priority instructions and execution
  permissions. Complete implementation through affected verification and a result
  report; a plan or investigation ends with its requested deliverable. Continue
  authorized work without repeated approval; identify actual blocking boundaries.
- Inspect the worktree and preserve unrelated changes. Refresh remote information
  when needed; do not merge, rebase, or switch branches merely to inspect it.
- Prefer the default branch when starting work without an established branch.
  Preserve an existing task branch; follow explicit user branch instructions.
  Never create or switch branches solely for a commit, push, release, or PR.
- Change or recommend branch protection only when explicitly asked. Honor explicit
  repository-specific direct-push exceptions; otherwise report a rejected push
  without bypassing protection or inventing a branch or PR.
- Unpublished implementation details may be redesigned; preserve existing public
  APIs, file formats, and saved-data compatibility unless a breaking change is
  authorized. Update affected producers, consumers, tests, examples, and docs.
- Fix verified root causes; do not hide failures with fallbacks or weaker checks.
  Document unavoidable workarounds and their removal conditions.
- Read relevant docs and run the repository's check entrypoint for the change and
  phase. Verify affected behavior; report checks run and omitted. Repeat or broaden
  successful checks only for new changes, failures, or unresolved concerns.
- For library metadata, require demonstrated incompatibility for exact pins or
  upper bounds; keep reproducibility locks separate.
- When editing READMEs, keep them concise with useful visuals inline; put extended
  guides in linked documentation.
- For GitHub push/release work, use `prepare-github-push` in `.agents/skills/`.
  Local-only commits need no version bump; GitHub pushes require one.
- For software performance work, use `benchmark-performance` in `.agents/skills/`.
  Performance claims require comparable measurements and equivalent output.
- For GitHub Actions edits, use `optimize-github-actions` in `.agents/skills/`.
  Preserve required coverage; never run untrusted PR code on self-hosted runners.
<!-- END KF AGENT POLICY -->

# Start Here

Read `README.md`, then [Development and Tests](docs/development-and-tests.md#choose-checks-for-a-change).
There is no separate CONTRIBUTING guide. For stage changes, follow the matching
`workflow/gg_*_entrypoint.sh` into `workflow/core/gg_*_core.sh`; shared helpers
live in `workflow/support/`. Use [Repository Layout](docs/repository-layout.md)
only when the owner is unclear.

Run commands from the repository root:

- `bash ./dev lint`: shell syntax, Ruff, and configuration checks on the host;
  see the development guide for prerequisites and the macOS Bash limitation.
- `bash ./dev config-check`: check entrypoint forwarding metadata (requires `python`).
- `bash ./dev check smoke`: existing small container smoke lane.
- `bash ./dev check fast workflow/tests/test_validation_runner.py -x`: example
  focused check; choose the actual file and lane from the development guide.
- `bash ./dev build`: provision a local Docker development image when needed;
  this downloads and compiles dependencies, so it is not a quick check.

There is no configured standalone type checker. The development guide is the
command and change-to-test reference; do not create a second runner or test list.
For repeated selection and result reporting, use
[validate-change](.agents/skills/validate-change/SKILL.md).

## Preserve Research Inputs and Contracts

Entrypoint editable blocks, `workflow/gg_common_params.sh`, scheduler directives,
and path overrides can contain project-specific settings. Preserve them unless
the task calls for changing them. Read [configuration](docs/configuration-and-common-parameters.md)
and [input conventions](docs/input-conventions.md) before changing forwarding,
species identifiers, CDS/protein mode, genetic codes, trait missingness, or
calibrations. Do not adjust scientific thresholds or models to make tests pass.
Preserve CLI/API behavior and output schemas, identifiers, and archive layouts;
check affected readers as well as writers.

Treat `workspace/input/` as curated data, and `workspace/output/`,
`workspace/downloads/`, `workspace/db*/`, SIFs, build caches, and benchmark/review
artifacts as data rather than incidental cleanup targets. Fixture regeneration
requires [dataset provenance](docs/test-dataset.md). Use test temporary directories
for reproductions; a normal workflow or debug harness can write persistent outputs.

Before finishing, review `git diff --check` and the final diff. Report commands,
runtime used, pass/fail/skip results, and checks blocked by missing prerequisites.
Separate static evidence from executed behavior and Docker from SIF validation.

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

Use a GeneGalleon container runtime for workflow integration, R helpers,
and toolchain-dependent validation. For those changes, read
[runtime validation](docs/agent-runtime-validation.md). Host syntax and narrow
static checks are sufficient only when runtime dependencies do not matter;
do not claim SIF compatibility from host or Docker results.

Do not absorb dependency-side defects into GeneGalleon or add workarounds for
older dependency behavior; identify and fix/update the owning dependency.

# Core Workflow Architecture

Keep `workflow/core/gg_*_core.sh` as self-contained workflow implementation scripts. Do not split their functions or ordered execution stages into `workflow/core/stages/`, per-stage shell fragments, or sourced function libraries merely to reduce file size or reorganize code.

Changes to this architecture require explicit user approval. Without that approval, edit the matching core script in place and preserve its existing execution order.
