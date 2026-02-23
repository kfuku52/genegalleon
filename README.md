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

## Table of Contents

1. [Repository Layout](#repository-layout)
2. [Execution Model](#execution-model)
3. [Container Build and Runtime](#container-build-and-runtime)
4. [Workspace Layout and Data Model](#workspace-layout-and-data-model)
5. [Input Conventions](#input-conventions)
6. [Quick Start](#quick-start)
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
    query_manifest/
    query_sra_id/
    query_gene/
    transcriptome_assembly/
  output/
  downloads/
```

### Legacy layout compatibility

Runtime helpers in `workflow/support/gg_util.sh` automatically resolve roots:

- input root: `workspace_input_root()`
- output root: `workspace_output_root()`
- downloads root: `workspace_downloads_root()`

Behavior:

- if `workspace/input` or `workspace/output` exists, split layout is used,
- otherwise single-root legacy behavior is used.

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
  - input: `workspace/input/transcriptome_assembly/amalgkit_metadata/GENUS_SPECIES.metadata.tsv`

### Automated provider formatting helper

Manual formatting can be replaced with `workflow/support/format_species_inputs.py`, which ports key rules from legacy `gfe_dataset` scripts such as:

- `data_formatting_ensemblplants.sh`
- `data_formatting_phycocosm.sh`
- `data_formatting_phytozome.sh`

Example (single provider, explicit input directory):

```bash
python workflow/support/format_species_inputs.py \
  --provider ensemblplants \
  --input-dir "/path/to/gfe_dataset/20230216_EnsemblPlants/original_files" \
  --species-cds-dir workspace/input/species_cds \
  --species-gff-dir workspace/input/species_gff
```

Example (all providers, `gfe_dataset` root):

```bash
python workflow/support/format_species_inputs.py \
  --provider all \
  --dataset-root "/path/to/gfe_dataset" \
  --species-cds-dir workspace/input/species_cds \
  --species-gff-dir workspace/input/species_gff
```

Notes:

- output filenames are normalized to start with `Genus_species_...`,
- CDS IDs are prefixed with `Genus_species_...`,
- common historical replacements are applied to CDS/GFF text,
- CDS are padded to codon-length multiples and duplicate IDs are removed.

Download-first workflow (manifest driven):

```bash
python workflow/support/format_species_inputs.py \
  --provider all \
  --download-manifest /path/to/download_manifest.tsv \
  --download-dir workspace/output/input_download_cache \
  --species-cds-dir workspace/input/species_cds \
  --species-gff-dir workspace/input/species_gff
```

Manifest required columns:

- `provider` (`ensemblplants`, `phycocosm`, `phytozome`)
- `species_key` (directory/file key used in legacy datasets)
- `cds_url`
- `gff_url`

Optional columns:

- `cds_filename`
- `gff_filename`

Use `--download-only` to fetch raw files and stop before formatting.
Use `--dry-run` to preview downloads and formatting outputs without writing files.

Optional download authentication:

- `--auth-bearer-token-env ENV_NAME`
  - adds `Authorization: Bearer <token>` from an environment variable.
- `--http-header "Key: Value"` (repeatable)
  - add arbitrary request headers (for provider-specific requirements).

Build a manifest automatically from an existing local dataset tree:

```bash
python workflow/support/build_download_manifest.py \
  --provider all \
  --dataset-root "/path/to/gfe_dataset" \
  --output workspace/input/query_manifest/download_manifest.tsv
```

This writes `file://` URLs, so you can test the full download+format pipeline locally before replacing them with remote URLs.

Use `workspace/input/query_manifest/download_manifest.tsv` as the runtime manifest file.

`gg_input_generation_entrypoint.sh` wrapper:

- runs `gg_input_generation_core.sh` inside the container as a single entrypoint,
- can build a manifest and immediately format inputs in one run,
- can validate produced `species_cds` and `species_gff` consistency.

Profile-based configuration:

- default profile file is provided at
  `workspace/input/query_manifest/inputprep_profile.sh`,
- edit values there for repeated runs,
- keep only these runtime files in `workspace/input/query_manifest`:
  - `inputprep_profile.sh`
  - `download_manifest.tsv`

Alternative runtime overrides (without editing files) via env vars:

- `GG_INPUT_PROVIDER`, `GG_INPUT_STRICT`, `GG_INPUT_OVERWRITE`,
- `GG_INPUT_DRY_RUN`, `GG_INPUT_DOWNLOAD_MANIFEST`,
- `GG_INPUT_DATASET_ROOT`, `GG_INPUT_INPUT_DIR`,
- `GG_INPUT_AUTH_BEARER_TOKEN_ENV`, `GG_INPUT_HTTP_HEADER`,
- `GG_INPUT_SUMMARY_OUTPUT`.
- `GG_INPUTPREP_CONFIG` (explicit profile path override).

Optional host dataset bind for `gg_input_generation_entrypoint.sh`:

```bash
export GG_INPUT_DATASET_HOST_PATH="/absolute/path/to/gfe_dataset"
bash workflow/gg_input_generation_entrypoint.sh
```

This binds the host path to `/external/gfe_dataset` in the container and sets `GG_DATASET_ROOT` for `gg_input_generation_core.sh`.

## Quick Start

From repository root:

```bash
cd workflow

# 1) Optional input download/format automation
bash gg_input_generation_entrypoint.sh

# 2) Unified genome-evolution pipeline
#    (runs inlined stages: speciesTree -> orthofinder -> genomeEvolution)
bash gg_genome_evolution_entrypoint.sh

# 3) Gene-family phylogeny/annotation
bash gg_gene_evolution_entrypoint.sh

# 4) Build orthogroup SQLite DB
bash gg_gene_database_entrypoint.sh

# 5) Progress summaries
bash gg_progress_summary_entrypoint.sh
```

Current default in `gg_gene_evolution_entrypoint.sh` is:

- `mode_query2family=1`
- `query_blast_method="diamond"`

So array-task cardinality is the number of files in `workspace/input/query_gene`.

## Main Stages and What They Do

### `gg_input_generation_entrypoint.sh`

Purpose:

- optional pre-stage to automate input preparation,
- build manifest tables from local `gfe_dataset` trees,
- optionally download provider files from a manifest and format into
  `workspace/input/species_cds` and `workspace/input/species_gff`.

Main scripts:

- `workflow/core/gg_input_generation_core.sh`
- `workflow/support/build_download_manifest.py`
- `workflow/support/format_species_inputs.py`

Notable defaults:

- when `run_build_manifest=1`, generated manifest is reused automatically by `run_format_inputs=1`,
- `run_validate_inputs=1` validates CDS naming rules and CDS/GFF species-set consistency,
- each run appends a TSV record to `workspace/output/inputprep/inputprep_runs.tsv` (configurable via `summary_output` or `GG_INPUT_SUMMARY_OUTPUT`).

Quick summary of recent inputPrep runs:

```bash
python workflow/support/summarize_inputprep_runs.py \
  --infile workspace/output/inputprep/inputprep_runs.tsv \
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
- DIAMOND annotation,
- optional MMseqs2 taxonomy and contamination removal,
- optional genome analyses (SubPhaser, dotplot, GenomeScope).

Main outputs:

- `workspace/output/species_cds_annotation`
- `workspace/output/species_cds_busco_full`, `species_cds_busco_short`
- `workspace/output/species_genome_busco_full`, `species_genome_busco_short`

Notable defaults:

- most heavy tasks default to `0`,
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
