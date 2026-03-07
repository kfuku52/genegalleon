# GeneGalleon

![GeneGalleon logo](logo/logo.png)

[![Tests](https://github.com/kfuku52/genegalleon/actions/workflows/tests.yml/badge.svg)](https://github.com/kfuku52/genegalleon/actions/workflows/tests.yml)
[![GitHub Release](https://img.shields.io/github/v/release/kfuku52/genegalleon)](https://github.com/kfuku52/genegalleon/releases)
[![GitHub Last Commit](https://img.shields.io/github/last-commit/kfuku52/genegalleon)](https://github.com/kfuku52/genegalleon/commits)
[![GitHub Issues](https://img.shields.io/github/issues/kfuku52/genegalleon)](https://github.com/kfuku52/genegalleon/issues)
[![GitHub Pull Requests](https://img.shields.io/github/issues-pr/kfuku52/genegalleon)](https://github.com/kfuku52/genegalleon/pulls)
[![GitHub Stars](https://img.shields.io/github/stars/kfuku52/genegalleon?style=social)](https://github.com/kfuku52/genegalleon/stargazers)
[![GitHub Forks](https://img.shields.io/github/forks/kfuku52/genegalleon?style=social)](https://github.com/kfuku52/genegalleon/network/members)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Bash](https://img.shields.io/badge/Bash-5%2B-121011?logo=gnubash&logoColor=white)](https://www.gnu.org/software/bash/)
[![Python](https://img.shields.io/badge/Python-3.10%2B-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![R](https://img.shields.io/badge/R-4%2B-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![Container](https://img.shields.io/badge/Container-Docker%20%7C%20Apptainer-2496ED?logo=docker&logoColor=white)](https://github.com/kfuku52/genegalleon/tree/master/container)
[![Scheduler](https://img.shields.io/badge/Scheduler-SLURM%20%7C%20UGE%20%7C%20PBS-6A5ACD)](https://github.com/kfuku52/genegalleon/tree/master/workflow)
[![FASTA Policy](https://img.shields.io/badge/FASTA-.fa.gz-2EA44F)](https://github.com/kfuku52/genegalleon/blob/master/README.md#compression-and-fasta-handling-policy)

**Navigating the Ocean of Genomic Histories**

[`GeneGalleon`](https://github.com/kfuku52/genegalleon) (/dʒiːn ˈɡæliən/) is a container-first comparative genomics and phylogenomics workflow suite.
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

From the repository root:

```bash
# Recommended: build a local SIF from a published GHCR image
apptainer build genegalleon.sif docker://ghcr.io/kfuku52/genegalleon:latest

# Alternative: build a SIF locally from this repository
IMAGE_SOURCE=local IMAGE=local/genegalleon TAG=dev bash ./gg_container_build_entrypoint.sh
```

### 2. Run the bundled quick start

```bash
cd workflow
bash gg_gene_evolution_entrypoint.sh
```

This writes results under:

- `workspace/output/query2family`
- `workspace/downloads`

### 3. Try other wrappers next

After the smoke test, common next steps are:

- `gg_input_generation_entrypoint.sh`: download and format provider inputs from a manifest, and optionally build `species_trait.tsv`
- `gg_genome_evolution_entrypoint.sh`: run the combined species-tree, orthogroup, and genome-evolution workflow
- `gg_transcriptome_generation_entrypoint.sh`: transcriptome assembly and quantification
- `gg_genome_annotation_entrypoint.sh`: CDS/genome annotation and QC
- `gg_gene_database_entrypoint.sh`: build the orthogroup SQLite database
- `gg_progress_summary_entrypoint.sh`: generate summary tables and reports from completed runs

## Table of Contents

1. [Quick Start](#quick-start)
2. [Repository Layout](#repository-layout)
3. [Execution Model](#execution-model)
4. [Container Build and Runtime](#container-build-and-runtime)
5. [Workspace Layout and Data Model](#workspace-layout-and-data-model)
6. [Input Conventions](#input-conventions)
7. [Main Stages and What They Do](#main-stages-and-what-they-do)
8. [Scheduler and Array Semantics](#scheduler-and-array-semantics)
9. [Configuration and Common Parameters](#configuration-and-common-parameters)
10. [Compression and FASTA Handling Policy](#compression-and-fasta-handling-policy)
11. [Troubleshooting](#troubleshooting)
12. [Development and Tests](#development-and-tests)
13. [License](#license)

## Repository Layout

- `workflow/`
  - stage command scripts: `core/gg_*_core.sh`
  - scheduler wrappers: `gg_*_entrypoint.sh`
  - shared parameters: `gg_common_params.sh`
  - utility scripts: `workflow/support/*`
- `workflow/tests/`
  - Python and R tests for key scripts
- `container/`
  - Docker Buildx scaffold and SIF conversion helper scripts
- `workspace/`
  - project-local input/output/downloads roots
- `logo/`, `LICENSE`, `README.md`

## Execution Model

Typical execution path:

1. Launch `workflow/gg_*_entrypoint.sh`.
2. Wrapper normalizes scheduler env and starts Singularity/Apptainer shell.
3. Wrapper pipes corresponding `core/gg_*_core.sh` into the container shell.
4. Command script resolves workspace roots, validates inputs, and runs stage tasks.

Most stage behavior is controlled by the top block in each `gg_*_entrypoint.sh`:

```bash
### Start: Modify this block to tailor your analysis ###
...
### End: Modify this block to tailor your analysis ###
```

Each wrapper forwards those variables into the container environment.
`core/gg_*_core.sh` scripts consume job-supplied configuration and execute stage logic.

## Container Build and Runtime

### One-command build (local/public selectable)

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
- release publish + SIF asset upload: `.github/workflows/release-sif.yml`
  - tags: `<release-tag>`, `YYYYMMDD-<sha7>`, `sha-<sha7>`
  - release assets: `<repo>_<release-tag>_amd64.sif` and `.sha256`

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

## Workspace Layout and Data Model

### Split layout (recommended)

```text
workspace/
  input/
    species_cds/
    species_expression/
    species_gff/
    species_genome/
    species_dnaseq/
    species_trait/
    species_rnaseq/
    input_generation/
    amalgkit_metadata/
    query_sra_id/
    query_gene/
  output/
  downloads/
```

### Runtime path resolution

Runtime helpers in `workflow/support/gg_util.sh` automatically resolve roots:

- input root: `workspace_input_root()`
- output root: `workspace_output_root()`
- downloads root: `workspace_downloads_root()`

Behavior:

- layout is always split: `workspace/input`, `workspace/output`, and `workspace/downloads`.

Downloads root:

- `workspace/downloads`

## Input Conventions

### Species naming

Species IDs are inferred from filename prefixes such as `Genus_species`.
Many joins/merges across pipeline stages rely on this convention.

### FASTA file extensions

Major scripts accept:

- `.fa`, `.fas`, `.fasta`, `.fna`
- `.fa.gz`, `.fas.gz`, `.fasta.gz`, `.fna.gz`

### `workspace/input/species_cds`

- one CDS FASTA per species,
- sequence IDs should be unique,
- avoid spaces and `|` in IDs when possible,
- species prefix on IDs (`Genus_species_...`) is strongly recommended.

### `workspace/input/query_gene`

Each file is one family-level task in `mode_query2family=1`.
Accepted forms:

- amino-acid FASTA,
- in-frame CDS FASTA,
- plain gene ID list (one ID per line).

`gg_gene_evolution_core.sh` auto-detects query type:

- FASTA with DNA sequences: translated to AA query,
- FASTA with protein sequences: used directly,
- non-FASTA text: treated as gene IDs to extract CDS then translate.

### Transcriptome assembly input modes

`gg_transcriptome_generation_core.sh` supports three mutually exclusive modes:

- `mode_sraid=1`
  - input: `workspace/input/query_sra_id/GENUS_SPECIES.txt`
  - one SRA/BioProject ID per line
- `mode_fastq=1`
  - input: `workspace/input/species_rnaseq/GENUS_SPECIES/*.fastq.gz`
- `mode_metadata=1`
  - input: `workspace/input/amalgkit_metadata/GENUS_SPECIES_metadata.tsv`

For `mode_sraid=1` and `mode_fastq=1`, auto-generated amalgkit metadata is written to:

- `workspace/output/transcriptome_assembly/amalgkit_metadata/GENUS_SPECIES_metadata.tsv`

This keeps `workspace/input/amalgkit_metadata` reserved for explicit `mode_metadata=1` inputs.

### Automated provider formatting helper

Manual formatting can be replaced with `workflow/support/format_species_inputs.py`, which ports key rules from legacy raw genome formatting scripts such as:

- `data_formatting_ensemblplants.sh`
- `data_formatting_phycocosm.sh`
- `data_formatting_phytozome.sh`

Example (single provider, explicit input directory):

```bash
python workflow/support/format_species_inputs.py \
  --provider ensemblplants \
  --input-dir "/path/to/raw_genomes/20230216_EnsemblPlants/original_files" \
  --species-cds-dir workspace/output/input_generation/species_cds \
  --species-gff-dir workspace/output/input_generation/species_gff
```

Example (all providers, provider-root input directory):

```bash
python workflow/support/format_species_inputs.py \
  --provider all \
  --input-dir "/path/to/raw_genomes" \
  --species-cds-dir workspace/output/input_generation/species_cds \
  --species-gff-dir workspace/output/input_generation/species_gff
```

Notes:

- output filenames are normalized to start with `Genus_species_...`,
- formatted CDS outputs are always gzipped with `.fa.gz` extension,
- CDS IDs are prefixed with `Genus_species_...` and aggregated to one representative CDS per gene,
- common historical replacements are applied to CDS/GFF text,
- CDS are padded to codon-length multiples and transcript-level redundancies are collapsed at gene level.

Download-first workflow (manifest driven):

```bash
python workflow/support/format_species_inputs.py \
  --provider all \
  --download-manifest /path/to/download_plan.xlsx \
  --download-dir workspace/output/input_generation/tmp/input_download_cache \
  --species-cds-dir workspace/output/input_generation/species_cds \
  --species-gff-dir workspace/output/input_generation/species_gff
```

Manifest required columns:

- `provider` (required; first column in XLSX templates)
- `id` (required; second column in XLSX templates)
- `species_key` is optional
- `cds_url` and `gff_url` (or `id` to auto-resolve provider-specific URLs when supported)
  - for `provider=ncbi`, when `species_key` is omitted and `id` is given, `species_key` is inferred from NCBI species metadata (e.g. `Homo_sapiens`).
  - `provider=ncbi` accepts both `GCF_*` and `GCA_*` assembly accessions and auto-resolves NCBI assembly URLs.
  - for `provider=coge`, `id` must be CoGe `genome_id` (numeric `gid`), and CDS/GFF/Genome URLs are auto-built.
  - for `provider=cngb`, built-in inference resolves CNGB assembly IDs (`CNA...`, `cngb:...`) or linked `GCA/GCF` accessions and maps to downloadable assembly files.
  - for `provider=flybase`, `provider=wormbase`, and `provider=vectorbase`, `id` can be resolved via explicit URL columns or `GG_<PROVIDER>_*_URL_TEMPLATE`.
  - for `provider=ensembl` and `provider=ensemblplants`, `id`-only URL inference is supported via:
    - provider defaults (for example, Ensembl/EnsemblPlants index discovery),
    - or env templates: `GG_<PROVIDER>_CDS_URL_TEMPLATE`, `GG_<PROVIDER>_GFF_URL_TEMPLATE`, `GG_<PROVIDER>_GENOME_URL_TEMPLATE`,
    - and optional page template: `GG_<PROVIDER>_ID_URL_TEMPLATE`.
  - for `provider=local`, `id` can point to a local species directory (or set explicit local file paths) and formatting starts from local files.
  - `provider=phycocosm` and `provider=phytozome` are not supported in `--download-manifest`.
    Use `--input-dir` for local raw-file formatting only.

Optional columns:

- `species_key`
- `cds_filename`
- `gff_filename`
- `genome_filename`
- `cds_url_template`, `gff_url_template`, `genome_url_template`
  (`{id}`, `{species_key}`, `{provider}` placeholders)
- `local_cds_path`, `local_gff_path`, `local_genome_path`
  (for `provider=local`)

XLSX template notes:

- `download_plan.xlsx` includes drop-downs for `provider` and `id`.
- `id` drop-down values are provider-specific.
- provider drop-down order is fixed, and `local` is always listed last.
- for large provider (`ncbi`), five model-organism IDs are shown as examples (mixed `GCF_*`/`GCA_*` formats).
- for `coge` and `cngb`, IDs are also example-based by default (both show five model-organism examples).
- for `ensembl`, `ensemblplants`, `flybase`, `wormbase`, `vectorbase`, and `local`,
  IDs can be supplied from a prebuilt `id_options_snapshot.json`.
- when no snapshot is supplied, non-large providers fall back to IDs discovered from `--input-dir`.
- drop-down IDs are shown as `ID (Species name)` for non-`local` providers.
- for `provider=local`, drop-down IDs are plain local directory/path-style IDs.
- at runtime, `id` parsing uses the token before the first space (for non-`local` providers), so labels like `GCF_000001405.40 (Homo sapiens)` are accepted.
- any other value can still be typed manually.

Use `--download-only` to fetch raw files and stop before formatting.
Use `--dry-run` to preview downloads and formatting outputs without writing files.
Use `--jobs` to set download parallelism (defaults to `NSLOTS`, fallback `1`).
Resolved manifest rows (for example URL/template/local-path auto-resolution and filename filling)
are written to:
- `workspace/output/input_generation/download_plan.resolved.tsv`

Download concurrency safeguards:

- per-provider caps are applied even when `--jobs` is large.
- override per-provider cap with:
  `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_<PROVIDER>` (e.g. `..._NCBI=2`).
- NCBI E-utilities calls are throttled to documented limits by default:
  - `3 req/s` without API key,
  - `10 req/s` with `GG_NCBI_API_KEY`.
- override NCBI E-utilities rate with `GG_NCBI_EUTILS_MAX_RPS`.

Optional download authentication:

- `--auth-bearer-token-env ENV_NAME`
  - adds `Authorization: Bearer <token>` from an environment variable.
- `--http-header "Key: Value"` (repeatable)
  - add arbitrary request headers (for provider-specific requirements).

Build a manifest automatically from an existing local dataset tree:

```bash
python workflow/support/build_download_manifest.py \
  --provider all \
  --input-dir "/path/to/raw_genomes" \
  --id-options-snapshot /path/to/id_options_snapshot.json \
  --output workspace/input/input_generation/download_plan.xlsx
```

This writes `file://` URLs, so you can test the full download+format pipeline locally before replacing them with remote URLs.

Use `workspace/input/input_generation/download_plan.xlsx` as the runtime manifest input file.
Resolved/filled rows are written to `workspace/output/input_generation/download_plan.resolved.tsv`.

Latest template distribution:

- On each push to the default branch, GitHub Actions runs `build_download_manifest.py`
  and publishes the latest `download_plan.xlsx`.
- The workflow first builds `id_options_snapshot.json` (remote provider ID choices),
  then generates `download_plan.xlsx` from that snapshot.
- If remote ID fetch fails for some providers, the workflow reuses provider entries
  from the previous `id_options_snapshot.json` release asset.
- Release asset URL (rolling latest):
  `https://github.com/<owner>/<repo>/releases/download/download-plan-latest/download_plan.xlsx`
- Snapshot asset URL (rolling latest):
  `https://github.com/<owner>/<repo>/releases/download/download-plan-latest/id_options_snapshot.json`
- Workflow definition:
  `.github/workflows/download-plan.yml`

`gg_input_generation_entrypoint.sh` wrapper:

- runs `gg_input_generation_core.sh` inside the container as a single entrypoint,
- can download provider files from a manifest and format inputs in one run,
- can validate produced `species_cds` and `species_gff` consistency,
- can optionally generate `workspace/input/species_trait/species_trait.tsv`
  from configured trait databases.

Configuration:

- edit the `### Start: Modify this block ... ###` section in
  `workflow/gg_input_generation_entrypoint.sh`,
- keep `workspace/input/input_generation/download_plan.xlsx` as the runtime manifest input file.
- check resolved rows at `workspace/output/input_generation/download_plan.resolved.tsv`.

Alternative runtime overrides (without editing files) via env vars:

- `GG_INPUT_PROVIDER`, `GG_INPUT_STRICT`, `GG_INPUT_OVERWRITE`,
- `GG_INPUT_DOWNLOAD_ONLY`, `GG_INPUT_DRY_RUN`,
- `GG_INPUT_DOWNLOAD_TIMEOUT`,
- `GG_INPUT_DOWNLOAD_MANIFEST`, `GG_INPUT_INPUT_DIR`, `GG_INPUT_DOWNLOAD_DIR`,
- `GG_INPUT_SPECIES_CDS_DIR`, `GG_INPUT_SPECIES_GFF_DIR`, `GG_INPUT_SPECIES_GENOME_DIR`,
- `GG_INPUT_SPECIES_SUMMARY_OUTPUT`,
- `GG_INPUT_RESOLVED_MANIFEST_OUTPUT`,
- `GG_INPUT_AUTH_BEARER_TOKEN_ENV`, `GG_INPUT_HTTP_HEADER`,
- `GG_INPUT_SUMMARY_OUTPUT`,
- `GG_INPUT_RUN_FORMAT_INPUTS`, `GG_INPUT_RUN_VALIDATE_INPUTS`,
- trait generation:
  `GG_INPUT_TRAIT_PROFILE` (`none` or `gift_starter`; default `none`),
  `GG_INPUT_RUN_GENERATE_SPECIES_TRAIT`,
  `GG_INPUT_TRAIT_SPECIES_SOURCE` (`download_manifest` or `species_cds`; default `download_manifest`),
  `GG_INPUT_SPECIES_TRAIT_OUTPUT`,
  `GG_INPUT_TRAIT_PLAN`,
  `GG_INPUT_TRAIT_DATABASE_SOURCES`,
  `GG_INPUT_TRAIT_DATABASES`,
  `GG_INPUT_TRAIT_DOWNLOAD_DIR`,
  `GG_INPUT_TRAIT_DOWNLOAD_TIMEOUT`.
- per-provider download caps:
  `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_ENSEMBL`,
  `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_ENSEMBLPLANTS`,
  `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_PHYCOCOSM`,
  `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_PHYTOZOME`,
  `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_NCBI`,
  `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_REFSEQ`,
  `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_GENBANK`,
  `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_COGE`,
  `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_CNGB`,
  `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_FLYBASE`,
  `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_WORMBASE`,
  `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_VECTORBASE`.

Quick preset example (enable trait stage with GIFT starter config):

```bash
GG_INPUT_TRAIT_PROFILE=gift_starter bash workflow/gg_input_generation_entrypoint.sh
```
## Main Stages and What They Do

### `gg_input_generation_entrypoint.sh`

Purpose:

- optional pre-stage to automate input preparation,
- optionally download provider files from a manifest and format into
  `workspace/output/input_generation/species_cds` and
  `workspace/output/input_generation/species_gff` and
  `workspace/output/input_generation/species_genome`,
- optionally generate `workspace/input/species_trait/species_trait.tsv` using
  target species from `download_plan.xlsx` (default) or `species_cds`.

Main scripts:

- `workflow/core/gg_input_generation_core.sh`
- `workflow/support/format_species_inputs.py`
- `workflow/support/generate_species_trait.py`

Notable defaults:

- `run_validate_inputs=1` validates CDS naming rules and CDS/GFF species-set consistency,
- formatted outputs default to `workspace/output/input_generation/species_cds`, `workspace/output/input_generation/species_gff`, and `workspace/output/input_generation/species_genome`,
- trait generation is opt-in (`run_generate_species_trait=0` by default),
  and writes `workspace/input/species_trait/species_trait.tsv` when enabled,
- `trait_profile=gift_starter` is available as a quick preset:
  sets `run_generate_species_trait=1` and uses `trait_databases=gift`
  when `trait_databases` is unset or `auto`,
- each run appends a TSV record to `workspace/output/input_generation/gg_input_generation_runs.tsv`
  (configurable via `summary_output` or `GG_INPUT_SUMMARY_OUTPUT`).

Trait generation inputs:

- trait plan TSV (default: `workspace/input/input_generation/trait_plan.tsv`)
  - required columns: `database`, `source_column`, `output_trait`
  - optional columns: `value_type` (`numeric|binary|categorical`),
    `aggregation`, `positive_values` (comma-separated, for binary mapping),
    `trait_key`, `trait_key_column`
  - for `gift`, `trait_key` can be either a trait ID (`Lvl3`, e.g. `1.1.1`)
    or trait name (`Trait2`, e.g. `Woodiness_1`)
- database source map TSV (default: `workspace/input/input_generation/trait_database_sources.tsv`)
  - required column: `database`
  - optional columns: `acquisition_mode` (`bulk|species_api|gift_api`), `uri`,
    `species_column`, `delimiter`, `response_format`,
    `archive_member`, `trait_key_column`,
    `gift_version`, `gift_trait_ids`, `gift_page_size`,
    `gift_max_pages_per_trait`, `gift_bias_ref`, `gift_bias_deriv`,
    `gift_agreement_min`, `gift_versions_api`
- starter templates are bundled at:
  - `workspace/input/input_generation/trait_plan.tsv`
  - `workspace/input/input_generation/trait_database_sources.tsv`
  - bundled defaults are prewired for `austraits`, `eltontraits`, and
    `gift` (`gift_api` with official base endpoint), including starter
    trait names for `woodiness` (`Woodiness_1`, ID `1.1.1`),
    `plant_height_max_m` (`Plant_height_max`, ID `1.6.2`), and
    `seed_mass_mean_g` (`Seed_mass_mean`, ID `3.2.3`).

Trait DB retrieval policy in `generate_species_trait.py`:

- `bulk`: fetch/cache full dataset under `workspace/downloads/trait_datasets/<database>/`
  then subset target species.
  - `uri` may contain multiple comma-separated sources; all are merged.
  - ZIP sources are supported (`archive_member` to select file inside archive).
- `species_api`: query per target species and merge responses.
- `gift_api`: resolve target species `work_ID` via GIFT
  `names_matched_unique`, resolve trait tokens from `trait_key` /
  `gift_trait_ids` (ID or name via `traits_meta`), then query trait pages
  and filter to target species.

Supported DB IDs can be listed with:

```bash
python workflow/support/generate_species_trait.py --print-supported-databases
```

GIFT traits catalog can be inspected with:

```bash
python workflow/support/generate_species_trait.py --print-gift-traits --gift-trait-search wood --gift-trait-limit 20
```

Quick summary of recent gg_input_generation runs:

```bash
python workflow/support/summarize_gg_input_generation_runs.py \
  --infile workspace/output/input_generation/gg_input_generation_runs.tsv \
  --last-n 10
```

### `gg_transcriptome_generation_entrypoint.sh`

Purpose:

- metadata retrieval (amalgkit),
- FASTQ retrieval/preprocessing,
- transcriptome assembly (`Trinity` or `rnaSPAdes`),
- longest CDS extraction,
- optional contamination filtering,
- BUSCO and expression quantification summaries.

Main outputs:

- `workspace/output/transcriptome_assembly/assembled_transcripts_with_isoforms`
- `workspace/output/transcriptome_assembly/longest_cds`
- `workspace/output/transcriptome_assembly/amalgkit_quant`
- `workspace/output/transcriptome_assembly/amalgkit_merge`

Notable defaults:

- `assembly_method="rnaSPAdes"`
- `kallisto_reference="longest_cds"`
- staged FASTA outputs are written as `.fa.gz`.

### `gg_genome_annotation_entrypoint.sh`

Purpose:

- per-species CDS/genome annotation and QC,
- BUSCO (CDS/genome),
- UniProt annotation (`blastp` or `mmseqs2`),
- optional MMseqs2 taxonomy and contamination removal,
- optional genome analyses (SubPhaser, dotplot, GenomeScope).

Main outputs:

- `workspace/output/species_cds_annotation`
- `workspace/output/species_cds_busco_full`, `species_cds_busco_short`
- `workspace/output/species_genome_busco_full`, `species_genome_busco_short`

Notable defaults:

- most heavy tasks default to `0`,
- `uniprot_annotation_method="mmseqs2"` (set `blastp` to use NCBI BLASTP for UniProt annotation),
- `run_multispecies_summary=1` by default.

### `gg_genome_evolution_entrypoint.sh`

Purpose:

- unified genome-evolution entrypoint that serially runs:
  - inlined species-tree stage (formerly `gg_speciesTree_core.sh`)
  - inlined orthofinder stage (formerly `gg_orthofinder_core.sh`)
  - inlined genome-evolution stage (formerly `gg_genomeEvolution_core.sh`)
- preserves output-exists skip behavior at step level and aborts on real failures.

Main output roots:

- `workspace/output/species_tree`
- `workspace/output/orthofinder`
- `workspace/output/genome_evolution`

### Inlined Stage: Orthofinder

Purpose:

- translate CDS to species proteins,
- run OrthoFinder,
- select orthogroups for downstream analysis,
- compare orthogroup methods.

Main outputs:

- `workspace/output/orthofinder`

Temporary protein FASTA files are created under `workspace/downloads/tmp/` and removed automatically after orthogroup-related steps finish.

Notable defaults:

- `orthogroup_table="HOG"`
- `orthogroup_annotation_method="mmseqs2"` (set `blastp` to use NCBI BLASTP for representative-gene UniProt annotation)
- species-tree-aware OrthoFinder is used when species tree exists.

### Inlined Stage: Species Tree

Purpose:

- BUSCO-based single-copy extraction,
- per-gene alignments and trees,
- concatenated and ASTRAL species trees,
- IQ2MC/mcmctree dating pipeline,
- constrained-tree plotting.

Main outputs:

- `workspace/output/species_tree/species_tree_summary/undated_species_tree.nwk`
- `workspace/output/species_tree/species_tree_summary/dated_species_tree.nwk`
- BUSCO and intermediate tree/alignment directories under `species_tree/`.

Notable defaults:

- `mode_busco=1`
- `mode_orthogroup=0` (currently marked not supported)
- `busco_lineage="embryophyta_odb12"`
- `genetic_code=1`
- single-copy FASTA/alignment outputs are standardized to `.fa.gz`, and legacy
  plain `.fasta` files in species-tree intermediate directories are auto-migrated.

### `gg_gene_evolution_entrypoint.sh`

Purpose:

- per-family phylogeny inference and annotation,
- optional reconciliation/rooting/dating,
- expression and promoter analyses,
- OU/PGLS/convergence analyses,
- summary and tree plotting.

Main outputs:

- `workspace/output/query2family/*` in query2family mode,
- `workspace/output/orthogroup/*` in orthogroup mode.

Notable defaults in current snapshot:

- `mode_query2family=1`
- `run_rps_blast=1` (Pfam domain annotation is on by default),
- `run_tree_plot=1`
- `run_summary=1`
- `uniprot_annotation_method="mmseqs2"` (set `blastp` for NCBI BLASTP-based UniProt annotation),
- many advanced analyses default to `0`.

Current behavior notes:

- if `run_generax=1`, initial IQ-TREE disables UFBOOT and support is
  recomputed after GeneRax on the GeneRax topology,
- Pfam RPS-BLAST DB (`Pfam_LE`) is auto-prepared when missing, with lock-based
  synchronization for array jobs,
- gene-tree/species-tree PGLS outputs are `gene_tree_PGLS.tsv` and
  `species_tree_PGLS.tsv`.

Wrapper-specific note:

- `gg_gene_evolution_entrypoint.sh` forwards all top-block variables to the container runtime and also accepts environment-variable overrides for the same names.

### Inlined Stage: Genome Evolution

Purpose:

- BUSCO-based and orthogroup-based GRAMPA workflows,
- optional CAFE and GO enrichment analyses.

Main outputs:

- `workspace/output/genome_evolution/*`

Notable defaults:

- polyploidization-related BUSCO/GRAMPA tasks are enabled,
- `run_cafe=0`, `run_go_enrichment=0` by default,
- GO target can be specified by species name or branch ID.

### `gg_gene_database_entrypoint.sh`

Purpose:

- assemble SQLite DB from orthogroup summary tables.

Main output:

- `workspace/output/orthogroup/gg_orthogroup.db`

Required input directories:

- `workspace/output/orthogroup/stat_tree`
- `workspace/output/orthogroup/stat_branch`

### `gg_progress_summary_entrypoint.sh`

Purpose:

- summarize stage outputs into TSV reports.

Main outputs in current working directory:

- `orthogroup_summary.tsv`
- `transcriptome_assembly_summary.tsv`

Note:

- this stage runs `workflow/core/gg_progress_summary_core.sh` inside the container.

### Utility wrappers

- Versions dump:
  - collector script: `workflow/support/gg_versions.sh`
  - auto-triggered at the end of each `gg_*_entrypoint.sh` on successful completion
  - outputs are written to `workspace/output/versions/*.log`

## Scheduler and Array Semantics

Wrappers include templates for:

- SLURM
- UGE
- PBS

You can run with scheduler submission (`sbatch`, `qsub`) or direct `bash`.

Array-size rules:

- `gg_genome_evolution_entrypoint.sh`: fixed single task.
- `gg_gene_evolution_entrypoint.sh`:
  - `mode_orthogroup=1`: number of rows in `Orthogroups.GeneCount.selected.tsv` (excluding header),
  - `mode_query2family=1`: number of files in `workspace/input/query_gene`.
- `gg_genome_annotation_entrypoint.sh`: number of input species CDS files.
- `gg_transcriptome_generation_entrypoint.sh`: number of species input units for the chosen mode.

## Configuration and Common Parameters

### Per-stage configuration

Edit the top config block in each `workflow/gg_*_entrypoint.sh`:

```bash
### Start: Modify this block to tailor your analysis ###
...
### End: Modify this block to tailor your analysis ###
```

`workflow/core/gg_*_core.sh` files use `Job-supplied configuration` blocks and are not the primary place for routine parameter tuning.

### Shared common parameter file

`workflow/gg_common_params.sh` currently defines:

- `GG_COMMON_GENETIC_CODE` (default `1`)
- `GG_COMMON_BUSCO_LINEAGE` (default `embryophyta_odb12`)
- `GG_COMMON_CONTAMINATION_REMOVAL_RANK` (default `phylum`)
- `GG_COMMON_OUTGROUP_LABELS` (default `Oryza_sativa`)
- `GG_COMMON_ANNOTATION_REPRESENTATIVE_SPECIES` (default `Arabidopsis_thaliana`)
- `GG_COMMON_MCMCTREE_DIVERGENCE_TIME_CONSTRAINTS_STR` (default `Arabidopsis_thaliana,Oryza_sativa,130,-`)
- `GG_COMMON_TREEVIS_CLADE_ORTHOLOG_PREFIX` (default `Arabidopsis_thaliana_`)
- `GG_COMMON_GRAMPA_H1` (default empty)
- `GG_COMMON_TARGET_BRANCH_GO` (default `<1>`)

Current stage scripts auto-load `workflow/gg_common_params.sh` when available.
Their top config blocks use these defaults via parameter expansion, for example:

```bash
genetic_code="${genetic_code:-${GG_COMMON_GENETIC_CODE:-1}}"
busco_lineage="${busco_lineage:-${GG_COMMON_BUSCO_LINEAGE:-embryophyta_odb12}}"
outgroup_labels="${outgroup_labels:-${GG_COMMON_OUTGROUP_LABELS:-Oryza_sativa}}"
```

## Compression and FASTA Handling Policy

Current policy:

- input discovery across major stages supports both plain and gzipped FASTA,
- pipeline-tracked FASTA outputs (especially `file_*` targets in `core/gg_*_core.sh`)
  are standardized to `.fa.gz`,
- in-house scripts prefer gzipped FASTA directly; plain FASTA is generated only
  as temporary/intermediate input when a third-party tool requires uncompressed files.

Species-tree and gene-evolution intermediate directories were updated to keep
main retained FASTA/alignment artifacts gzipped with `.fa.gz` naming.

## Troubleshooting

- Stage skipped unexpectedly:
  - check `run_*` flags in corresponding `gg_*_entrypoint.sh`.
- Array task exits immediately:
  - verify `SGE_TASK_ID`/array range against actual input count.
- Optional steps not running:
  - confirm required optional inputs (expression/GFF/genome/trait files).
- Missing tree-dependent outputs in gene-family stage:
  - ensure species-tree prerequisites exist under `workspace/output/species_tree`.
- Taxonomy errors (ETE/NCBI):
  - inspect `workspace/downloads/ete_taxonomy` and rerun with clean lock state if needed.
- General environment diagnostics:

```bash
cd workflow
ls -1 ../workspace/output/versions/*.log
```

## Development and Tests

Run Python tests:

```bash
python -m pytest -q workflow/tests
```

Run R parse checks for script syntax:

```bash
Rscript workflow/tests/test_r_scripts_parse.R
```

## License

This repository is distributed under the MIT License. See `LICENSE`.
