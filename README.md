# GeneGalleon

![GeneGalleon logo](logo/logo.png)

GeneGalleon is a container-first comparative genomics and phylogenomics workflow suite.
It provides scheduler-ready staged pipelines for:

- transcriptome assembly and expression quantification,
- CDS/genome annotation and contamination filtering,
- local pairwise fractionation-bias and within-genome self-synteny retention analysis with kfFractBias,
- orthogroup inference,
- species-tree inference and dating,
- gene-family phylogeny and trait/evolutionary annotation,
- genome-evolution downstream analyses,
- orthogroup database generation.

The same entrypoints run through SLURM, UGE, PBS, or directly with Bash.

For repository development, `bash ./dev check fast` selects an authoritative
GeneGalleon SIF or Docker runtime automatically. See
[Development and Tests](docs/development-and-tests.md).

## Quick Start

Prepare either a repo-root SIF or a pulled Docker image, then run the bundled
query2family example:

```bash
# Linux/HPC with Apptainer or Singularity
IMAGE_SOURCE=public IMAGE=ghcr.io/kfuku52/genegalleon TAG=latest bash ./gg_container_build_entrypoint.sh

# Docker-only host
docker pull ghcr.io/kfuku52/genegalleon:latest

cd workflow
bash gg_gene_evolution_entrypoint.sh
```

This writes results under `workspace/output/query2family`. Docker-only,
local-build, HPC, and reproducible-tag options are covered in
[Container Build and Runtime](docs/container-build-and-runtime.md) and
[Common Workflow Recipes](docs/common-workflow-recipes.md).

## Output Examples

Query2family tree plot:

![Query2family tree plot example](docs/assets/example-plots/query2family-tree-plot.png)

Gene-family presence/absence summary:

![Gene-family presence/absence example](docs/assets/example-plots/query2family-presence-absence.png)

See [Example Plots](docs/example-plots.md) for more outputs.

## Updating with an AI agent

To update GeneGalleon without losing project settings, use:

> Update GeneGalleon to the latest `origin/main`. Preserve all project-specific
> parameter values, scheduler settings, workspace inputs, and other local
> changes. Integrate upstream changes into the current configuration instead of
> replacing customized files wholesale. If a parameter was renamed or removed,
> migrate it appropriately or ask before proceeding. Do not use destructive Git
> commands or delete workspace data. Follow `AGENTS.md`, validate the updated
> version in a GeneGalleon container, and summarize the preserved settings,
> version change, and test results.

## Documentation

Detailed guides are split by topic:

- [Repository Layout](docs/repository-layout.md)
- [Execution Model](docs/execution-model.md)
- [Common Workflow Recipes](docs/common-workflow-recipes.md)
- [Gene-Family Outputs and Progress Monitoring](docs/gene-family-outputs-and-progress-monitoring.md)
- [Species-Tree Stage ZIP Storage](docs/species-tree-stage-zip-storage.md)
- [Migrating Legacy Unzipped Workspaces to ZIP Storage (audit, conversion, and rollback)](docs/workspace-storage-management.md)
- [Example Plots](docs/example-plots.md)
- [Container Build and Runtime](docs/container-build-and-runtime.md)
- [Workspace Layout and Data Model](docs/workspace-layout-and-data-model.md)
- [Input Conventions](docs/input-conventions.md)
- [Main Stages and What They Do](docs/main-stages-and-what-they-do.md)
- [Scheduler and Array Semantics](docs/scheduler-and-array-semantics.md)
- [Site Runtime Profiles](docs/site-runtime-profiles.md)
- [SHIROKANE AGE Guide](docs/shirokane-age.md)
- [Configuration and Common Parameters](docs/configuration-and-common-parameters.md)
- [HGT Detection Research](docs/hgt-detection-research.md)
- [Compression and FASTA Handling Policy](docs/compression-and-fasta-handling-policy.md)
- [Troubleshooting](docs/troubleshooting.md)
- [Development and Tests](docs/development-and-tests.md)

## License

This repository is distributed under the MIT License. See [LICENSE](LICENSE).
