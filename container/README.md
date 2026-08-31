# Multi-Arch `GeneGalleon` Container Runtime

This directory provides a reproducible Docker-first runtime build for `GeneGalleon`
that can target both:
- `linux/amd64` (x86_64)
- `linux/arm64` (AArch64, Apple Silicon compatible runtime via Linux VM/container host)

## Why this runtime exists

From the project README and Wiki (`gg_versions`):
- `GeneGalleon` was originally assembled interactively from a miniconda3 Singularity sandbox.
- The runtime now uses a single conda `base` env, with selected tools installed from upstream GitHub at build time
  (`kfuku52/amalgkit`, `kfuku52/cdskit`, `kfuku52/csubst`, `kfuku52/nwkit`,
  `kfuku52/kfl1ou`, `kfuku52/kfFractBias`, `kfuku52/kftools`, `kfuku52/rkftools`, `kfuku52/RADTE`).
  Standard builds follow the moving branches in `source_branches.env` for every
  one of these repositories. Build wrappers resolve those branches once at the
  start of a build, in parallel, so all target architectures use the same snapshot. Explicit
  `*_REPO_SHA` variables remain available only as one-off overrides.

So this Dockerfile is designed as:
1. controlled and auditable base build,
2. architecture-aware optional package handling,
3. build-time runtime validation against required command lists,
4. explicit post-build manual steps for licensed/large assets.

The Dockerfile frontend and Micromamba base image are digest-pinned. Updating
either is an explicit maintenance change, not an implicit effect of rebuilding.
APT and conda packages are still resolved from their current repositories, so
builds are not claimed to be bit-for-bit reproducible.

See `container/CAPABILITY_MATRIX.md` for expected parity by architecture.

## Development and distribution targets

`BUILD_TARGET=runtime` is the default for build wrappers and CI publication.
It contains the scientific runtime without the APT development packages in
`apt/development.txt`. `bash ./dev build` selects `development`, which adds
those compilers, headers, and build utilities to the same runtime. Conda
packages required by R and other scientific tools remain in both targets,
including their compiler dependencies.

```bash
# Development image: local/genegalleon:dev, without SIF conversion.
bash ./dev build

# Distribution image, with a separate local tag.
BUILD_TARGET=runtime BUILD_SIF=0 IMAGE_SOURCE=local \
IMAGE=local/genegalleon TAG=runtime bash ./gg_container_build_entrypoint.sh
```

Each of the eleven moving upstream sources has an independent build stage.
With the BuildKit cache retained, changing one source revision rebuilds its
wheel, R package, or binary, then assembles and validates the runtime without
rebuilding the other sources.
The runtime receives one shared Conda environment and the small artifacts,
without source checkouts or wheel archives. Intermediate stages can stay in
the builder cache and are not shipped in Docker or SIF images. Builder disk use
can increase even when the distribution image shrinks.

`GG_BUILD_JOBS` limits native Make/R builds to two jobs by default. BuildKit
schedules independent stages separately; this is not a global concurrency
limit. Native Apptainer builds use the same artifact scripts sequentially and
remove APT development packages for `BUILD_TARGET=runtime` before validation.

See [container development](../docs/container-development.md) for cache behavior,
resource controls, and measured image sizes and rebuild times.

## Build multi-arch image with Docker Buildx

```bash
chmod +x container/buildx.sh
IMAGE=ghcr.io/<your-org>/genegalleon TAG=20260211 MODE=push ./container/buildx.sh
```

Exact source commits and checksums can be overridden at build time:

```bash
KFU52_CSUBST_REPO_SHA=<40-character-commit-sha> \
BUSCO_REPO_SHA=<40-character-commit-sha> \
PAML_REPO_SHA=<40-character-commit-sha> \
TESTNH_TARBALL_SHA256=598337183d2cec9c61cd364fab255a270062844b0ba5172913f7cf97512c43e2 \
CAFE5_TARBALL_SHA256=71871bdc74c2ffc7c1c0f4500f4742f2ff46a15cfaba78dc179d21bb1ba67ba8 \
IMAGE=ghcr.io/<your-org>/genegalleon TAG=20260211 MODE=push ./container/buildx.sh
```

Default source behavior:
- Git-sourced programs follow the moving branches recorded in `source_branches.env`
- build wrappers resolve every branch once per build in parallel and pass the resulting commits into Docker; resolution fails before building if any source is unavailable
- `/opt/pg/logs/source_revisions.tsv` records the effective source revision in both Docker and native Apptainer images
- `BUSCO` and `paml` follow their upstream `master` branches too
- `Notung`, `BioPP/testnh`, and `CAFE5` release archives are verified with SHA-256 before extraction
- GitHub/GitLab source fetches prefer release/archive downloads and fall back to `git` retry logic only when needed

Override rules:
- an explicitly supplied `*_REPO_SHA` takes precedence over the moving branch for that one build
- source SHA overrides should be full 40-character commit SHAs
- `BUSCO_MIRROR_REPO_URL` is optional and is only used as a secondary source if the primary `BUSCO_REPO_URL` fetch fails.
- if you override a repo URL to a fork, also supply a commit SHA that exists in that fork

Do not copy a resolved build commit into repository defaults. To change which
moving branch standard builds follow, edit `source_branches.env` and run the
full multi-architecture container validation suite.

`buildx.sh` runs a preflight check to ensure the conda env set used in
`workflow/core/gg_*_core.sh` is covered by env installs in `container/Dockerfile`.

Single-platform local test image:

```bash
IMAGE=local/genegalleon TAG=dev PLATFORMS=linux/arm64 MODE=load ./container/buildx.sh
```

The selected Buildx builder's internal cache is used by default. A separate
local cache export is disabled because exporting the full GeneGalleon cache can
write several gigabytes after every build. Set `USE_LOCAL_CACHE=1` only when a
portable cache directory is needed; `CACHE_DIR` defaults to `.buildx-cache`.

For `MODE=load`, `buildx.sh` fingerprints the Docker build inputs and skips
BuildKit entirely when the tagged local image already has the same fingerprint.
Set `SKIP_UNCHANGED_LOAD=0` to force a rebuild.
The same shared fingerprint is embedded by local, scheduled, and release builds
as `io.genegalleon.build-input`; OCI labels also record the repository version
and MIT license. Both fingerprints include the target, also recorded as
`io.genegalleon.build-target`, so development and runtime images cannot be
mistaken for each other. The fingerprint includes a UTC daily
`SECURITY_REFRESH_EPOCH` (default: the current `YYYY-MM-DD`). A small final
image layer refreshes Ubuntu package indexes, installs all available upgrades,
and fails if an upgrade would still remain. The epoch is also recorded as
`io.genegalleon.security-refresh-epoch` and in
`/opt/pg/logs/security_refresh_epoch.txt`.
An additional `io.genegalleon.runtime-input` label omits only repository
revision/version metadata. Validation CI keys its SIF cache by this exact
runtime fingerprint, including target, security epoch, every resolved source
revision, and the complete copied container context. Platform is a separate
part of the SIF cache key. Cached SIFs are saved only
after authoritative runtime tests succeed.

## One-command build (local/public selectable)

```bash
IMAGE_SOURCE=local IMAGE=local/genegalleon TAG=dev bash ./gg_container_build_entrypoint.sh
```

Defaults:
- `IMAGE_SOURCE=auto`
- `MODE=load`
- `BUILD_TARGET=runtime` (`development` adds APT build tools for local builds)
- `PLATFORMS`: inferred from host arch (`linux/amd64` or `linux/arm64`)
- `OUT=./genegalleon.sif`
- `IMAGE_SOURCE=local`: build from `container/Dockerfile` via Docker Buildx, or build a SIF natively with Apptainer/Singularity when Docker is unavailable
- `IMAGE_SOURCE=public`: pull `docker://IMAGE:TAG` directly with Apptainer/Singularity
- `IMAGE_SOURCE=auto`: prefer local build when Docker is available, otherwise fall back to a public image when `BUILD_SIF=1`

Useful overrides:
- `BUILD_SIF=0` to skip `.sif` conversion
- `ENGINE=singularity` to force Singularity instead of automatic runtime detection
- `NATIVE_BUILD_FAKEROOT=always` to force `--fakeroot` for native local builds on sites that support it

Use `gg_container_build_entrypoint.sh` when you want the default SIF path pinned
to `<repo-root>/genegalleon.sif`. Bare `apptainer build ...` writes to the
current working directory unless you pass an explicit output path.

Public image example:

```bash
IMAGE_SOURCE=public IMAGE=ghcr.io/kfuku52/genegalleon TAG=latest bash ./gg_container_build_entrypoint.sh
```

Native local build example on a Docker-less host:

```bash
IMAGE_SOURCE=local IMAGE=local/genegalleon TAG=dev bash ./gg_container_build_entrypoint.sh
```

If your site requires explicit rootless escalation for definition-file builds, retry with:

```bash
NATIVE_BUILD_FAKEROOT=always IMAGE_SOURCE=local bash ./gg_container_build_entrypoint.sh
```

## Run wrappers against a Docker image directly

You can keep the existing `gg_*_entrypoint.sh` interface and opt into a
Docker-backed launcher instead of requiring `genegalleon.sif`.

```bash
docker pull ghcr.io/kfuku52/genegalleon:latest
GG_CONTAINER_RUNTIME=docker \
GG_CONTAINER_DOCKER_IMAGE=ghcr.io/kfuku52/genegalleon:latest \
bash workflow/gg_gene_evolution_entrypoint.sh
```

For a local image built with `MODE=load`:

```bash
GG_CONTAINER_RUNTIME=docker \
GG_CONTAINER_DOCKER_IMAGE=local/genegalleon:dev \
bash workflow/gg_gene_evolution_entrypoint.sh
```

## CI publishing and reproducible tags

GitHub Actions now publishes GHCR images and release SIF assets:

- `.github/workflows/container-ghcr.yml`
  - schedule: daily at 04:00 JST
  - resolves every moving upstream branch, computes the expected build-input
    label for both architectures, and skips only when both match the published
    `latest` image; missing or stale published metadata triggers a rebuild
  - advances the security-refresh epoch daily, so cached system-package layers
    cannot indefinitely retain superseded Ubuntu security packages
  - tags: `YYYYMMDD-<sha7>-<source-hash12>`, `sha-<sha7>`, `latest`
- `.github/workflows/release-sif.yml`
  - trigger: Release `published`
  - tags: `<release-tag>`, `YYYYMMDD-<sha7>-<source-hash12>`, `sha-<sha7>`
  - release assets: `<repo>_<release-tag>_amd64.sif` and `.sha256`
  - oversized `SIF` files fall back to a 90-day workflow artifact

Retention policy:

- workflow `SIF` artifacts are short-lived convenience copies
- immutable GHCR tags are the long-term source of truth
- immutable tags include the repository SHA and resolved-source fingerprint and
  are never overwritten; `sha-*`, `latest`, and release tags are aliases
- historical `SIF` files should be recreated from immutable GHCR tags when needed

User-side reproducible pull example:

```bash
IMAGE_SOURCE=public IMAGE=ghcr.io/<owner>/genegalleon TAG=20260304-abcd123-8f3a2c41d905 \
bash ./gg_container_build_entrypoint.sh
```

## Convert registry or Docker-daemon image to `.sif`

```bash
chmod +x container/apptainer_from_docker.sh

# From a registry image
IMAGE=ghcr.io/<your-org>/genegalleon TAG=20260211 ./container/apptainer_from_docker.sh

# From an image already loaded in the Docker daemon
SOURCE=docker-daemon IMAGE=local/genegalleon TAG=dev ./container/apptainer_from_docker.sh
```

## Important caveats

- In the strict profile, optional manifests are intentionally empty; non-empty
  `failed_optional_*.txt` logs indicate drift or unexpected resolution issues.
- On `linux/arm64`, the current required profile excludes `Trinity` and
  `jellyfish`:
  - `Trinity` is not reliably solvable in the strict channel profile.
  - available arm64 `jellyfish` builds are Python-only and do not provide the
    required `jellyfish` CLI binary.
- Repeat-analysis pipeline steps that depended on the removed `repeat` conda env
  were dropped from `workflow/core/gg_genome_annotation_core.sh`.
- Required/optional command validation report is generated at build time:
  - `/opt/pg/logs/runtime_validation_<arch>.tsv`
- `ete4==4.4.0` is installed via `pip` in `base` because `nwkit` imports
  `ete4` modules directly for tree rooting/transfer/timetree operations.
- Database paths required by pipeline scripts must be populated manually:
  - `/usr/local/db/Pfam_LE`
  - `/usr/local/db/uniprot_sprot.pep` (and derived DIAMOND DB if needed)
  - `/usr/local/db/jaspar`
- `Notung` is downloaded at build time from the official Notung 2.9 source
  and installed as:
  - `/usr/local/bin/Notung.jar`
- `BUSCO` and `paml` are fetched from the current tips of their configured branches by default.
- `amalgkit`, `cdskit`, `csubst`, `nwkit`, `kfl1ou`, `kfFractBias`, `kftools`, `rkftools`, and
  `RADTE` install from the moving branches in `source_branches.env` by default.
- `Notung`, `BioPP/testnh`, and `CAFE5` archives are checksum-verified during build.
- The configured source is a checksum-verified upstream ZIP:
  - `NOTUNG_DOWNLOAD_PAGE=https://amberjack.compbio.cs.cmu.edu/Notung/Notung-2.9.1.5.zip`
- The corresponding default checksum is:
  - `NOTUNG_ZIP_SHA256=81cbff670ab4d2416c01eba503f81c454aa5a724b0982373dd17510113882ae6`
- `NOTUNG_DOWNLOAD_PAGE` may also point at the legacy download page if you
  want the build to resolve another `Notung-2.9.*.zip`; set its matching
  `NOTUNG_ZIP_SHA256` at the same time or verification will stop the build.
- If the official `amberjack.compbio.cs.cmu.edu` hostname has a transient DNS
  issue during build, override the fallback IP if needed:
  - `NOTUNG_DOWNLOAD_HOST_IP=128.2.205.60`

## Suggested validation per architecture

Inside built container:

```bash
if command -v micromamba >/dev/null 2>&1; then
  eval "$(micromamba shell hook --shell bash)"
elif [[ -f /opt/conda/etc/profile.d/conda.sh ]]; then
  source /opt/conda/etc/profile.d/conda.sh
fi
conda activate base
hyphy --version
iqtree --version
mapnh --help
```

Then run one pipeline smoke workflow with minimal data.

## Useful checks in repo

```bash
# Verify conda env coverage against pipeline scripts
container/scripts/check_env_coverage.sh .
```
