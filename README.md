# GeneGalleon

![GeneGalleon logo](logo/logo.png)

[![Tests](https://github.com/kfuku52/genegalleon/actions/workflows/tests.yml/badge.svg)](https://github.com/kfuku52/genegalleon/actions/workflows/tests.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**Navigating the Ocean of Genomic Histories**

`GeneGalleon` is a container-first comparative genomics and phylogenomics workflow suite.
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
11. [Tree/Taxonomy Compatibility Notes](#treetaxonomy-compatibility-notes)
12. [Troubleshooting](#troubleshooting)
13. [Development and Tests](#development-and-tests)
14. [License](#license)

## Repository Layout

- `workflow/`
  - stage command scripts: `gg_*_cmd.sh`
  - scheduler wrappers: `gg_*_job.sh`
  - shared parameters: `gg_common_params.sh`
  - utility scripts: `workflow/script/*`
- `workflow/tests/`
  - Python and R tests for key scripts
- `container/`
  - Docker Buildx scaffold and SIF conversion helper scripts
- `workspace/`
  - project-local input/output/db roots
- `logo/`, `LICENSE`, `README.md`

## Execution Model

Typical execution path:

1. Launch `workflow/gg_*_job.sh`.
2. Wrapper normalizes scheduler env and starts Singularity/Apptainer shell.
3. Wrapper pipes corresponding `gg_*_cmd.sh` into the container shell.
4. Command script resolves workspace roots, validates inputs, and runs stage tasks.

Most stage behavior is controlled by the top block in each `gg_*_cmd.sh`:

```bash
### Start: Modify this block to tailor your analysis ###
...
### End: Modify this block to tailor your analysis ###
```

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
    species_sra_list/
    query2family_input/
    transcriptome_assembly/
  output/
  db/
```

### Legacy layout compatibility

Runtime helpers in `workflow/script/gg_util.sh` automatically resolve roots:

- input root: `workspace_input_root()`
- output root: `workspace_output_root()`
- db root: `workspace_db_root()`

Behavior:

- if `workspace/input` or `workspace/output` exists, split layout is used,
- otherwise single-root legacy behavior is used.

DB root preference:

- preferred: `workspace/db`
- fallback: `workspace/output/db` (legacy split DB location)

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

### `workspace/input/query2family_input`

Each file is one family-level task in `mode_query2family=1`.
Accepted forms:

- amino-acid FASTA,
- in-frame CDS FASTA,
- plain gene ID list (one ID per line).

`gg_geneFamilyPhylogeny_cmd.sh` auto-detects query type:

- FASTA with DNA sequences: translated to AA query,
- FASTA with protein sequences: used directly,
- non-FASTA text: treated as gene IDs to extract CDS then translate.

### Transcriptome assembly input modes

`gg_transcriptomeAssembly_cmd.sh` supports three mutually exclusive modes:

- `mode_sraid=1`
  - input: `workspace/input/species_sra_list/GENUS_SPECIES.txt`
  - one SRA/BioProject ID per line
- `mode_fastq=1`
  - input: `workspace/input/species_rnaseq/GENUS_SPECIES/*.fastq.gz`
- `mode_metadata=1`
  - input: `workspace/input/transcriptome_assembly/amalgkit_metadata/GENUS_SPECIES.metadata.tsv`

## Quick Start

From repository root:

```bash
cd workflow

# 0) Optional environment check in container
bash gg_test_job.sh

# 1) Orthogroups from species CDS
bash gg_orthofinder_job.sh

# 2) Species tree inference and dating
bash gg_speciesTree_job.sh

# 3) Gene-family phylogeny/annotation
bash gg_geneFamilyPhylogeny_job.sh

# 4) (Optional) genome-evolution analyses
bash gg_genomeEvolution_job.sh

# 5) Build orthogroup SQLite DB
bash gg_orthogroupDatabasePrep_job.sh

# 6) Progress summaries
bash gg_progressSummary_job.sh
```

Current default in `gg_geneFamilyPhylogeny_cmd.sh` is:

- `mode_query2family=1`
- `query_blast_method="diamond"`

So array-task cardinality is the number of files in `workspace/input/query2family_input`.

## Main Stages and What They Do

### `gg_transcriptomeAssembly_job.sh`

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

- `assembly_method="Trinity"`
- `kallisto_reference="longest_cds"`
- several output FASTA files are `.fa.gz`/`.fasta.gz`.

### `gg_cdsAnnotation_job.sh`

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

### `gg_orthofinder_job.sh`

Purpose:

- translate CDS to species proteins,
- run OrthoFinder and SonicParanoid,
- select orthogroups for downstream analysis,
- compare orthogroup methods.

Main outputs:

- `workspace/output/species_protein`
- `workspace/output/orthofinder`
- `workspace/output/sonicparanoid`

Notable defaults:

- `orthogroup_table="HOG"`
- species-tree-aware OrthoFinder is used when species tree exists.

### `gg_speciesTree_job.sh`

Purpose:

- BUSCO-based single-copy extraction,
- per-gene alignments and trees,
- concatenated and ASTRAL species trees,
- IQ2MC/mcmctree dating pipeline,
- constrained-tree plotting.

Main outputs:

- `workspace/output/species_tree/undated_species_tree.nwk`
- `workspace/output/species_tree/dated_species_tree.nwk`
- BUSCO and intermediate tree/alignment directories under `species_tree/`.

Notable defaults:

- `mode_busco=1`
- `mode_orthogroup=0` (currently marked not supported)
- `busco_lineage="embryophyta_odb12"`
- `genetic_code=1`

### `gg_geneFamilyPhylogeny_job.sh`

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
- `run_tree_plot=1`
- `run_summary=1`
- many advanced analyses default to `0`.

Wrapper-specific note:

- `gg_geneFamilyPhylogeny_job.sh` sets conservative local defaults for some heavy steps via env vars (for example GeneRax and Pfam-related toggles). Review wrapper exports if you want pure command-script defaults.

### `gg_genomeEvolution_job.sh`

Purpose:

- BUSCO-based and orthogroup-based GRAMPA workflows,
- optional CAFE and GO enrichment analyses.

Main outputs:

- `workspace/output/genome_evolution/*`

Notable defaults:

- polyploidization-related BUSCO/GRAMPA tasks are enabled,
- `run_cafe=0`, `run_go_enrichment=0` by default,
- GO target can be specified by species name or branch ID.

### `gg_orthogroupDatabasePrep_job.sh`

Purpose:

- assemble SQLite DB from orthogroup summary tables.

Main output:

- `workspace/output/orthogroup/gg_orthogroup.db`

Required input directories:

- `workspace/output/orthogroup/stat.tree`
- `workspace/output/orthogroup/stat.branch`

### `gg_progressSummary_job.sh`

Purpose:

- summarize stage outputs into TSV reports.

Main outputs in current working directory:

- `orthogroup_summary.tsv`
- `transcriptome_assembly_summary.tsv`

Note:

- this stage is job-wrapper-only (`gg_progressSummary_cmd.sh` does not exist).

### Utility wrappers

- `gg_test_job.sh`: validates key packages and executables in container environments.
- `gg_versions_job.sh`: dumps installed versions and DB/asset paths.

## Scheduler and Array Semantics

Wrappers include templates for:

- SLURM
- UGE
- PBS

You can run with scheduler submission (`sbatch`, `qsub`) or direct `bash`.

Array-size rules:

- `gg_speciesTree_job.sh`: fixed single task.
- `gg_geneFamilyPhylogeny_job.sh`:
  - `mode_orthogroup=1`: number of rows in `Orthogroups.GeneCount.selected.tsv` (excluding header),
  - `mode_query2family=1`: number of files in `workspace/input/query2family_input`.
- `gg_cdsAnnotation_job.sh`: number of input species CDS files.
- `gg_transcriptomeAssembly_job.sh`: number of species input units for the chosen mode.

## Configuration and Common Parameters

### Per-stage configuration

Edit only each script's top config block:

```bash
### Start: Modify this block to tailor your analysis ###
...
### End: Modify this block to tailor your analysis ###
```

### Shared common parameter file

`workflow/gg_common_params.sh` currently defines:

- `GG_COMMON_GENETIC_CODE` (default `1`)
- `GG_COMMON_BUSCO_LINEAGE` (default `embryophyta_odb12`)

If you want project-wide shared values, source this file explicitly in stage scripts:

```bash
source "${dir_myscript}/../gg_common_params.sh"
genetic_code="${GG_COMMON_GENETIC_CODE}"
busco_lineage="${GG_COMMON_BUSCO_LINEAGE}"
```

## Compression and FASTA Handling Policy

Current behavior is intentionally mixed but compatible:

- input discovery across major stages supports both plain and gzipped FASTA,
- some outputs are gzipped by design (for example transcriptome longest-CDS and contamination-removal outputs),
- some downstream outputs remain plain `.fasta` in current scripts.

Tree plotting helpers were updated to be robust with compressed alignments:

- candidate order: `.fa.gz` -> `.fasta` -> `.fa`,
- plain FASTA is generated temporarily only when required by downstream tools.

## Tree/Taxonomy Compatibility Notes

### Branch label naming

Pipeline branch labeling in outputs/scripts is standardized to `branch_id`.

If your environment uses an older `rkftools::table2phylo()` expecting `numerical_label`, tree plotting can fail with:

- `Missing required columns in table2phylo(): numerical_label`

Recommended action:

- update `rkftools` to a `branch_id`-compatible implementation.

### ETE4 requirement

Tree/taxonomy helper scripts import `ete4` directly.

ETE taxonomy DB setup in `gg_util.sh` uses ETE4 APIs and lock-protected initialization to avoid concurrent DB corruption.

## Troubleshooting

- Stage skipped unexpectedly:
  - check `run_*` flags in corresponding `gg_*_cmd.sh`.
- Array task exits immediately:
  - verify `SGE_TASK_ID`/array range against actual input count.
- Optional steps not running:
  - confirm required optional inputs (expression/GFF/genome/trait files).
- Missing tree-dependent outputs in gene-family stage:
  - ensure species-tree prerequisites exist under `workspace/output/species_tree`.
- Taxonomy errors (ETE/NCBI):
  - inspect `workspace/db/ete_taxonomy` and rerun with clean lock state if needed.
- General environment diagnostics:

```bash
cd workflow
bash gg_test_job.sh
bash gg_versions_job.sh
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
