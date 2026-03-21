# Multi-Arch `GeneGalleon` Container Scaffold

This directory provides a reproducible Docker-first build scaffold for `GeneGalleon`
that can target both:
- `linux/amd64` (x86_64)
- `linux/arm64` (AArch64, Apple Silicon compatible runtime via Linux VM/container host)

## Why this scaffold exists

From the project README and Wiki (`gg_versions`):
- `GeneGalleon` was originally assembled interactively from a miniconda3 Singularity sandbox.
- The runtime now uses a single conda `base` env, with selected tools installed from upstream refs at build time
  (`kfuku52/amalgkit`, `kfuku52/cdskit`, `kfuku52/csubst`, `kfuku52/nwkit`,
  `kfuku52/kftools`, `kfuku52/rkftools`, `kfuku52/RADTE`).
  User-authored tools follow their configured refs by default, while explicit
  `*_REPO_SHA` pins remain available as an override when needed.

So this Dockerfile is designed as:
1. reproducible base build,
2. architecture-aware optional package handling,
3. build-time runtime validation against required command lists,
4. explicit post-build manual steps for licensed/large assets.

See `container/CAPABILITY_MATRIX.md` for expected parity by architecture.

## Build multi-arch image with Docker Buildx

```bash
chmod +x container/buildx.sh
IMAGE=ghcr.io/<your-org>/genegalleon TAG=20260211 MODE=push ./container/buildx.sh
```

Source refs, optional pins, and checksums can be overridden at build time:

```bash
KFU52_AMALGKIT_REPO_SHA= \
KFU52_AMALGKIT_REPO_REF=kfdevel \
KFU52_REPO_REF=master \
BUSCO_REPO_SHA=6278721a1916f6da310e03ec9674099028c927a4 \
PAML_REPO_SHA=8daeead6b55523f375d9ac56dcfac38373ef8a2e \
KFL1OU_REPO_REF=master \
KFL1OU_REPO_SHA= \
KFTOOLS_REPO_URL=https://github.com/kfuku52/kftools.git \
KFTOOLS_REPO_REF=master \
KFTOOLS_REPO_SHA= \
RKFTOOLS_REPO_URL=https://github.com/kfuku52/rkftools.git \
RKFTOOLS_REPO_REF=master \
RKFTOOLS_REPO_SHA= \
RADTE_REPO_URL=https://github.com/kfuku52/RADTE.git \
RADTE_REPO_REF=master \
RADTE_REPO_SHA= \
TESTNH_TARBALL_SHA256=598337183d2cec9c61cd364fab255a270062844b0ba5172913f7cf97512c43e2 \
CAFE5_TARBALL_SHA256=71871bdc74c2ffc7c1c0f4500f4742f2ff46a15cfaba78dc179d21bb1ba67ba8 \
IMAGE=ghcr.io/<your-org>/genegalleon TAG=20260211 MODE=push ./container/buildx.sh
```

Default hardening behavior:
- user-authored source installs follow their configured refs by default for `amalgkit`, `cdskit`, `csubst`, `nwkit`, `kfl1ou`, `kftools`, `rkftools`, and `RADTE`
- `BUSCO` and `paml` remain pinned by default
- `BioPP/testnh` and `CAFE5` release tarballs are verified with SHA-256 before extraction
- GitHub/GitLab source fetches prefer release/archive downloads and fall back to `git` retry logic only when needed

Override rules:
- `KFU52_REPO_REF` applies to `cdskit`, `csubst`, and `nwkit` only when the corresponding `KFU52_*_REPO_SHA` is empty.
- `KFL1OU_REPO_REF` and `RADTE_REPO_REF` default to `KFU52_REPO_REF` and are used when the corresponding `*_REPO_SHA` is empty.
- `KFTOOLS_REPO_REF`/`RKFTOOLS_REPO_REF` override only `kftools`/`rkftools` and default to `KFU52_REPO_REF` when set and the corresponding `*_REPO_SHA` is empty.
- `BUSCO_MIRROR_REPO_URL` is optional and is only used as a secondary source if the primary `BUSCO_REPO_URL` fetch fails.
- If you override a repo URL to a fork, also update the matching `*_REPO_SHA` or clear it to fall back to the ref/default branch.

For `amalgkit`, branch auto-selection can still be controlled with:
- `KFU52_AMALGKIT_AUTO_SELECT_REF=1` (default)
- `KFU52_AMALGKIT_BRANCH_CANDIDATES=master,kfdevel,devel` (default)
- `KFU52_AMALGKIT_REPO_REF=<branch>` (hard override)
This logic only takes effect when `KFU52_AMALGKIT_REPO_SHA` is empty.

`buildx.sh` runs a preflight check to ensure the conda env set used in
`workflow/core/gg_*_core.sh` is covered by env installs in `container/Dockerfile`.

Single-platform local test image:

```bash
IMAGE=local/genegalleon TAG=dev PLATFORMS=linux/arm64 MODE=load ./container/buildx.sh
```

## One-command build (local/public selectable)

```bash
IMAGE_SOURCE=local IMAGE=local/genegalleon TAG=dev bash ./gg_container_build_entrypoint.sh
```

Defaults:
- `IMAGE_SOURCE=auto`
- `MODE=load`
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
  - runs only when the previous JST day had container-related changes on the default branch
  - tags: `YYYYMMDD-<sha7>`, `sha-<sha7>`, `latest`
- `.github/workflows/release-sif.yml`
  - trigger: Release `published`
  - tags: `<release-tag>`, `YYYYMMDD-<sha7>`, `sha-<sha7>`
  - release assets: `<repo>_<release-tag>_amd64.sif` and `.sha256`
  - oversized `SIF` files fall back to a 90-day workflow artifact

Retention policy:

- workflow `SIF` artifacts are short-lived convenience copies
- immutable GHCR tags are the long-term source of truth
- historical `SIF` files should be recreated from immutable GHCR tags when needed

User-side reproducible pull example:

```bash
IMAGE_SOURCE=public IMAGE=ghcr.io/<owner>/genegalleon TAG=20260304-abcd123 \
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
- `ete4` is installed via `pip` in `base` because `nwkit` imports `ete4`
  modules directly for tree rooting/transfer/timetree operations.
- Database paths required by pipeline scripts must be populated manually:
  - `/usr/local/db/Pfam_LE`
  - `/usr/local/db/uniprot_sprot.pep` (and derived DIAMOND DB if needed)
  - `/usr/local/db/jaspar`
- `Notung` is downloaded at build time from the official Notung 2.9 source
  and installed as:
  - `/usr/local/bin/Notung.jar`
- `BUSCO` and `paml` are fetched from pinned upstream source snapshots by default.
- `amalgkit`, `cdskit`, `csubst`, `nwkit`, `kfl1ou`, `kftools`, `rkftools`, and `RADTE` follow their configured refs by default.
- `BioPP/testnh` and `CAFE5` tarballs are checksum-verified during build.
- The default source is the pinned stable ZIP:
  - `NOTUNG_DOWNLOAD_PAGE=https://amberjack.compbio.cs.cmu.edu/Notung/Notung-2.9.1.5.zip`
- `NOTUNG_DOWNLOAD_PAGE` may also point at the legacy download page if you
  want the build to resolve the latest stable `Notung-2.9.*.zip` there.
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
