# Workspace Layout and Data Model

GeneGalleon uses a split workspace layout and treats the workspace as the persistent project state.

## Default workspace root

By default, the host-side workspace is:

- `workflow/../workspace`

At runtime, that directory is bind-mounted into the container as:

- `/workspace`

You can relocate the workspace for a run by exporting:

```bash
gg_workspace_dir=/path/to/workspace bash workflow/gg_genome_evolution_entrypoint.sh
```

## Split layout

```text
workspace/
  input/
    species_cds/
    species_protein/
    species_genetic_code/
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

The layout is always split. Legacy mixed-layout assumptions are no longer part of the current runtime helpers.

## Directory roles

### `workspace/input/`

This is the user-managed side of the workspace.

Typical contents:

- `species_cds/`: CDS FASTA inputs for comparative/genome-evolution workflows
- `species_protein/`: optional per-species protein FASTA inputs for `gg_genome_evolution` protein mode
- `species_genetic_code/`: optional `species_genetic_code.tsv` table for per-species CDS translation overrides
- `species_gff/` and `species_genome/`: optional annotation/genome companions
- `species_expression/`: expression matrices used by gene-family downstream analyses
- `species_trait/`: trait tables such as `species_trait.tsv`
- `species_rnaseq/`, `query_sra_id/`, `amalgkit_metadata/`: transcriptome-generation inputs
- `query_gene/`: per-family query files for query2family runs
- `input_generation/`: manifests and trait-plan inputs for input-generation runs

### `workspace/output/`

This is the main results area. Major stage roots include:

- `workspace/output/input_generation`
- `workspace/output/transcriptome_assembly`
- `workspace/output/species_cds_annotation`
- `workspace/output/species_tree`
- `workspace/output/orthofinder`
- `workspace/output/genome_evolution`
- `workspace/output/query2family`
- `workspace/output/orthogroup`
- `workspace/output/csubst_site`
- `workspace/output/versions`

### `workspace/downloads/`

This is the reusable runtime cache area.

Typical contents include:

- downloaded databases,
- provider download caches,
- temporary working directories,
- lock files for array-safe downloads,
- `ete_taxonomy/` for ETE taxonomy data,
- `pfam/` for Pfam/RPS-BLAST assets,
- `trait_datasets/` for cached trait database downloads.

## Runtime path helpers

Runtime helpers in `workflow/support/gg_util.sh` define canonical paths:

- input root: `workspace_input_root()`
- output root: `workspace_output_root()`
- downloads root: `workspace_downloads_root()`
- taxonomy root: `workspace_taxonomy_root()`
- Pfam root: `workspace_pfam_root()`
- taxonomy DB file: `workspace_taxonomy_dbfile()`

These helpers are what core scripts use internally; hard-coding alternate layouts is likely to break stage interoperability.

## Ownership and lifecycle

As a rule of thumb:

- treat `workspace/input/` as human-curated inputs,
- treat `workspace/output/` as reproducible stage results,
- treat `workspace/downloads/` as disposable-but-reusable cache/state.

That means:

- editing files under `workspace/output/` by hand is usually a bad idea,
- cleaning selected caches under `workspace/downloads/` is reasonable when repairing stale downloads,
- copying final curated inputs from `workspace/output/input_generation/` into long-term project storage is often useful.

## Summary and report files outside `output/`

Some wrappers intentionally write reports at the workspace root rather than under `output/`.

Current examples:

- `workspace/orthogroup_summary.tsv`
- `workspace/transcriptome_assembly_summary.tsv`

Those are produced by `gg_progress_summary_entrypoint.sh`.

## Common relocation scenario

If you want to keep code and data separate:

```bash
gg_workspace_dir=/data/projects/my_run/workspace \
gg_container_image_path=/data/containers/genegalleon.sif \
bash workflow/gg_genome_evolution_entrypoint.sh
```

In that setup:

- the checked-out repository still provides `/script`,
- the external workspace becomes `/workspace`,
- outputs and caches stay outside the git tree.
