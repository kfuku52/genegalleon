# Multi-Arch `GeneGalleon` Container Scaffold

This directory provides a reproducible Docker-first build scaffold for `GeneGalleon`
that can target both:
- `linux/amd64` (x86_64)
- `linux/arm64` (AArch64, Apple Silicon compatible runtime via Linux VM/container host)

## Why this scaffold exists

From the project README and Wiki (`gg_versions`):
- `GeneGalleon` was originally assembled interactively from a miniconda3 Singularity sandbox.
- The runtime now uses a single conda `base` env, with selected tools installed from GitHub at build time
  (`kfuku52/amalgkit`, `kfuku52/cdskit`, `kfuku52/csubst`, `kfuku52/nwkit`,
  `kfuku52/kftools`, `kfuku52/rkftools`, `kfuku52/RADTE`).
  `amalgkit` is selected from the newer branch commit between `master` and
  `kfdevel` by default.

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

GitHub refs can be pinned at build time:

```bash
KFU52_REPO_REF=master \
KFU52_AMALGKIT_BRANCH_CANDIDATES=master,kfdevel \
KFTOOLS_REPO_URL=https://github.com/kfuku52/kftools.git \
RKFTOOLS_REPO_URL=https://github.com/kfuku52/rkftools.git \
RADTE_REPO_URL=https://github.com/kfuku52/RADTE.git \
KFTOOLS_REPO_REF=master \
RKFTOOLS_REPO_REF=master \
RADTE_REPO_REF=master \
IMAGE=ghcr.io/<your-org>/genegalleon TAG=20260211 MODE=push ./container/buildx.sh
```

`KFU52_REPO_REF` applies to `cdskit`, `csubst`, and `nwkit`.  
`KFTOOLS_REPO_REF`/`RKFTOOLS_REPO_REF` override only `kftools`/`rkftools`
and default to `KFU52_REPO_REF` when not set.
`RADTE_REPO_REF` controls the RADTE Git checkout and defaults to repository HEAD
when not set.
For `amalgkit`, auto-selection can be controlled with:
- `KFU52_AMALGKIT_AUTO_SELECT_REF=1` (default)
- `KFU52_AMALGKIT_BRANCH_CANDIDATES=master,kfdevel` (default)
- `KFU52_AMALGKIT_REPO_REF=<branch>` (hard override)

`buildx.sh` runs a preflight check to ensure the conda env set used in
`workflow/core/gg_*_core.sh` is covered by env installs in `container/Dockerfile`.

Single-platform local test image:

```bash
IMAGE=local/genegalleon TAG=dev PLATFORMS=linux/arm64 MODE=load ./container/buildx.sh
```

## Convert Docker image to `.sif`

```bash
chmod +x container/apptainer_from_docker.sh
IMAGE=ghcr.io/<your-org>/genegalleon TAG=20260211 OUT=genegalleon_20260211_amd.sif ./container/apptainer_from_docker.sh
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
- `Notung` is downloaded at build time from the Notung 2.9 download page and
  the latest stable `Notung-2.9.*.zip` found there is installed as:
  - `/usr/local/bin/Notung.jar`
- Override the source page if needed:
  - `NOTUNG_DOWNLOAD_PAGE=https://amberjack.compbio.cs.cmu.edu/Notung/download29.html`

## Suggested validation per architecture

Inside built container:

```bash
source /home/.bashrc
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
