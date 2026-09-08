<!-- BEGIN KF AGENT POLICY: source=https://github.com/kfuku52/kf-agent-policy; version=9; sha256=03e7ad2c21924fa609040d9176d1a9c3a7f0c6785f2efe97dfe03e48be13411e -->
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
- Keep top-level READMEs concise and retain useful visuals inline. Put
  feature-specific guides and extended examples in dedicated documentation or
  the wiki, linking only as needed.
- Proactively use visuals when they improve understanding.
- Changes confined to unpushed local commits need no backward compatibility.
- Prefer verified root-cause fixes to fallbacks or relaxed validation that only
  hide failures. Document unavoidable workarounds and their removal conditions.
- When changing GitHub Actions, preserve required coverage and never execute
  untrusted pull-request code on self-hosted runners.
- Run checks appropriate to the change and all repository-required checks.
  Directly verify affected behavior or artifacts; report what did and did not
  run. After success, expand or repeat checks only for new changes, failures,
  or unresolved concerns.
- Performance claims require representative before-and-after measurements and
  equivalent output.
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
