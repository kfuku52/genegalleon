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
  - tags: `YYYYMMDD-<sha7>`, `sha-<sha7>`, `latest`
- release publish + SIF build/upload workflow: `.github/workflows/release-sif.yml`
  - tags: `<release-tag>`, `YYYYMMDD-<sha7>`, `sha-<sha7>`
  - release assets always include `.sha256`
  - `<repo>_<release-tag>_amd64.sif` is uploaded to the GitHub Release when it is under the release asset size limit; otherwise it is kept as a workflow artifact and the Release gets a download notice instead

For reproducible runs, use an immutable tag:

```bash
apptainer build genegalleon_20260304_abcd123.sif \
  docker://ghcr.io/kfuku52/genegalleon:20260304-abcd123
```

You can also use release tags:

```bash
apptainer build genegalleon_v1.2.3.sif \
  docker://ghcr.io/kfuku52/genegalleon:v1.2.3
```

### Build Docker image (multi-arch)

```bash
IMAGE=ghcr.io/<your-org>/genegalleon TAG=dev MODE=push ./container/buildx.sh
```

Single-platform local load example:

```bash
IMAGE=local/genegalleon TAG=dev PLATFORMS=linux/arm64 MODE=load ./container/buildx.sh
```

### Convert Docker image to SIF

```bash
IMAGE=ghcr.io/<your-org>/genegalleon TAG=dev OUT=genegalleon_20260217_amd.sif ./container/apptainer_from_docker.sh
```

Wrappers expect the SIF at repo root by default:

- `workflow/../genegalleon.sif`

Runtime profile highlights in the current container scaffold:

- single conda runtime env: `base` (`biotools`/`r` split envs are obsolete),
- `iqtree` is conda-pinned to `3.*`,
- `pigz` is included for fast compression/decompression,
- `Notung.jar` is installed during image build,
- BUSCO is installed from source (`v6.0.0`) during image build.
