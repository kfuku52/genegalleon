# GeneGalleon

![GeneGalleon logo](logo/logo.png)

**Navigating the Ocean of Genomic Histories**

GeneGalleon is a container-first comparative genomics and phylogenomics workflow suite.
It provides shell-based staged pipelines for:

- transcriptome assembly and expression quantification,
- CDS/genome annotation and contamination filtering,
- orthogroup inference,
- species-tree inference and dating,
- gene-family phylogeny and trait/evolutionary annotation,
- genome-evolution downstream analyses,
- orthogroup database generation.

The repository is designed for scheduler execution (SLURM/UGE/PBS) with Singularity/Apptainer wrappers, while still allowing local runs for development and debugging.

## Quick Start

The fastest way to try GeneGalleon is to run
`gg_gene_evolution_entrypoint.sh` against the bundled test data in the default
`workspace/`.

### 1. Prepare the container image

```bash
# Apptainer + GHCR: build a local SIF from the published GHCR image
apptainer build genegalleon.sif docker://ghcr.io/kfuku52/genegalleon:latest

# Apptainer + local: build a local SIF from this repository
IMAGE_SOURCE=local IMAGE=local/genegalleon TAG=dev bash ./gg_container_build_entrypoint.sh

# Docker + GHCR: pull the published GHCR image
docker pull ghcr.io/kfuku52/genegalleon:latest

# Docker + local: build a local Docker image from this repository
IMAGE_SOURCE=local BUILD_SIF=0 IMAGE=local/genegalleon TAG=dev MODE=load bash ./gg_container_build_entrypoint.sh
```

The quick start in Step 2 expects `./genegalleon.sif`, so use one of the
Apptainer commands above for the workflow wrappers.

### 2. Run the bundled quick start

```bash
cd workflow
bash gg_gene_evolution_entrypoint.sh
```

This writes results under:

- `workspace/output/query2family`
- `workspace/downloads`

## Documentation

Detailed guides are split by topic:

- [Repository Layout](docs/repository-layout.md)
- [Execution Model](docs/execution-model.md)
- [Common Workflow Recipes](docs/common-workflow-recipes.md)
- [Container Build and Runtime](docs/container-build-and-runtime.md)
- [Workspace Layout and Data Model](docs/workspace-layout-and-data-model.md)
- [Input Conventions](docs/input-conventions.md)
- [Main Stages and What They Do](docs/main-stages-and-what-they-do.md)
- [Scheduler and Array Semantics](docs/scheduler-and-array-semantics.md)
- [Configuration and Common Parameters](docs/configuration-and-common-parameters.md)
- [Compression and FASTA Handling Policy](docs/compression-and-fasta-handling-policy.md)
- [Troubleshooting](docs/troubleshooting.md)
- [Development and Tests](docs/development-and-tests.md)

## License

This repository is distributed under the MIT License. See [LICENSE](LICENSE).
