# Container Build and Runtime

### One-command build (local/public selectable)

```bash
IMAGE_SOURCE=local IMAGE=local/genegalleon TAG=dev bash ./gg_container_build_entrypoint.sh
```

Defaults:
- `IMAGE_SOURCE=auto`
- `MODE=load`
- `PLATFORMS`: inferred from host arch (`linux/amd64` or `linux/arm64`)
- `BUILD_SIF`: `1` on Linux, `0` on macOS hosts where Apptainer/Singularity is typically unavailable
- `OUT=<repo-root>/genegalleon.sif`
- `IMAGE_SOURCE=local`: build from `container/Dockerfile` via Docker Buildx, or build a SIF natively with Apptainer/Singularity when Docker is unavailable and `BUILD_SIF=1`
- `IMAGE_SOURCE=public`: pull `docker://IMAGE:TAG` directly with Apptainer/Singularity
- `IMAGE_SOURCE=auto`: prefer local Docker builds when Docker Buildx is available; otherwise, when `BUILD_SIF=1`, pull a registry image directly or fall back to the published GHCR image

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

### Use CI-published GHCR images (recommended for users)

This repository now includes CI workflows that publish container images to GHCR:

- periodic publish: `.github/workflows/container-ghcr.yml`
  - schedule: daily at 04:00 JST
  - runs only when the previous JST day had container-related changes on the default branch
  - tags: `YYYYMMDD-<sha7>`, `sha-<sha7>`, `latest`
- release publish + SIF build/upload workflow: `.github/workflows/release-sif.yml`
  - tags: `<release-tag>`, `YYYYMMDD-<sha7>`, `sha-<sha7>`
  - release assets always include `.sha256`
  - `<repo>_<release-tag>_amd64.sif` is uploaded to the GitHub Release when it is under the release asset size limit; otherwise it is kept as a workflow artifact for 90 days and the Release gets a download notice instead

Operational retention policy:

- large workflow `SIF` artifacts are short-lived convenience copies and expire after 90 days
- long-term reproducibility comes from immutable GHCR tags such as `YYYYMMDD-<sha7>`
- recreate a historical `SIF` from GHCR when needed instead of storing old `SIF` artifacts indefinitely

For reproducible runs, use an immutable tag:

```bash
IMAGE_SOURCE=public IMAGE=ghcr.io/kfuku52/genegalleon TAG=20260304-abcd123 \
bash ./gg_container_build_entrypoint.sh
```

You can also use release tags:

```bash
IMAGE_SOURCE=public IMAGE=ghcr.io/kfuku52/genegalleon TAG=v1.2.3 \
bash ./gg_container_build_entrypoint.sh
```

### Build Docker image (multi-arch)

```bash
IMAGE=ghcr.io/<your-org>/genegalleon TAG=dev MODE=push ./container/buildx.sh
```

Single-platform local load example:

```bash
IMAGE=local/genegalleon TAG=dev PLATFORMS=linux/arm64 MODE=load ./container/buildx.sh
```

### Convert registry or Docker-daemon image to SIF

```bash
# From a registry image
IMAGE=ghcr.io/<your-org>/genegalleon TAG=dev ./container/apptainer_from_docker.sh

# From an image already loaded in the Docker daemon
SOURCE=docker-daemon IMAGE=local/genegalleon TAG=dev ./container/apptainer_from_docker.sh
```

Wrappers expect the SIF at repo root by default:

- `workflow/../genegalleon.sif`

### Run wrappers against a Docker image directly

This opt-in mode keeps the existing `gg_*_entrypoint.sh` interface but swaps the
container launcher for a Docker-backed shim. It is intended for local/dev hosts
that have Docker but no usable Apptainer/Singularity installation.

Published image example:

```bash
docker pull ghcr.io/kfuku52/genegalleon:latest
bash workflow/gg_gene_evolution_entrypoint.sh
```

When the default repo-root `genegalleon.sif` is missing, wrappers auto-fallback
to `ghcr.io/kfuku52/genegalleon:latest` if that image has already been pulled.
Explicit `GG_CONTAINER_RUNTIME=docker` / `GG_CONTAINER_DOCKER_IMAGE=...` still
override auto-detection.

Local image example:

```bash
IMAGE_SOURCE=local BUILD_SIF=0 IMAGE=local/genegalleon TAG=dev MODE=load \
bash ./gg_container_build_entrypoint.sh

GG_CONTAINER_RUNTIME=docker \
GG_CONTAINER_DOCKER_IMAGE=local/genegalleon:dev \
bash workflow/gg_gene_evolution_entrypoint.sh
```

Runtime notes:

- missing repo-root `genegalleon.sif` triggers Docker auto-fallback when a
  pulled Docker image is available; current fallback priority is
  `ghcr.io/kfuku52/genegalleon:latest`, then `local/genegalleon:dev`
- `GG_CONTAINER_RUNTIME=docker` forces Docker-backed wrapper mode
- `GG_CONTAINER_DOCKER_IMAGE=<image:tag>` selects the image passed to `docker run`
- `gg_container_image_path` remains the SIF path knob; do not point it at a Docker image reference.
- The wrapper still binds the workspace to `/workspace` and the workflow tree to `/script`.

Runtime profile highlights in the current container scaffold:

- single conda runtime env: `base` (`biotools`/`r` split envs are obsolete),
- `iqtree` is conda-pinned to `3.*`,
- `pigz` is included for fast compression/decompression,
- `Notung.jar` is installed during image build,
- BUSCO is installed from source (`v6.0.0`) during image build.
