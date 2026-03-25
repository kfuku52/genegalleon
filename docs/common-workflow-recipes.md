# Common Workflow Recipes

This page collects practical stage combinations that are not obvious from the per-stage reference pages alone.

## 1. Run the bundled smoke test

Use this when you want to confirm the container and default workspace are working.

```bash
IMAGE_SOURCE=public IMAGE=ghcr.io/kfuku52/genegalleon TAG=latest bash ./gg_container_build_entrypoint.sh
cd workflow
bash gg_gene_evolution_entrypoint.sh
```

Default behavior:

- `mode_gene_evolution=query2family`
- queries are read from `workspace/input/query_gene`
- outputs are written under `workspace/output/query2family`

## 2. Prepare species inputs from a manifest

Use `gg_input_generation_entrypoint.sh` when you want one wrapper to:

- resolve provider downloads,
- fetch raw files,
- format species CDS/GFF/genome inputs,
- optionally build `species_trait.tsv`.

Typical run:

```bash
GG_INPUT_DOWNLOAD_MANIFEST="$PWD/workspace/input/input_generation/download_plan.xlsx" \
bash workflow/gg_input_generation_entrypoint.sh
```

If you also want the bundled GIFT starter trait workflow:

```bash
GG_INPUT_TRAIT_PROFILE=gift_starter \
bash workflow/gg_input_generation_entrypoint.sh
```

Useful outputs to inspect afterward:

- `workspace/output/input_generation/download_plan.resolved.tsv`
- `workspace/output/input_generation/gg_input_generation_runs.tsv`
- `workspace/output/input_generation/gg_input_generation_species.tsv`
- `workspace/output/input_generation/annotation_summary/annotation_summary.tsv`

Wrapper note:

- `provider=refseq` and `provider=genbank` are accepted by the wrapper as aliases of `ncbi`.

If `input_generation_mode` is unclear, treat it as follows:

- `single`: the normal mode; one run does everything.
- `array_prepare`: generate `task_plan.json` before launching array workers.
- `array_worker`: one scheduler task handles one species indexed by `GG_ARRAY_TASK_ID`.
- `array_finalize`: merge worker outputs and run shared checks and summaries once.

Typical array sequence:

```bash
GG_INPUT_INPUT_GENERATION_MODE=array_prepare \
bash workflow/gg_input_generation_entrypoint.sh
```

```bash
GG_INPUT_INPUT_GENERATION_MODE=array_worker \
sbatch workflow/gg_input_generation_entrypoint.sh
```

```bash
GG_INPUT_INPUT_GENERATION_MODE=array_finalize \
bash workflow/gg_input_generation_entrypoint.sh
```

## 3. Start from transcriptome-oriented inputs

Use `gg_transcriptome_generation_entrypoint.sh` when your primary starting point is RNA-seq.

The top config block defaults to `mode_transcriptome_assembly="auto"` when exactly one input layout is present. Set the mode explicitly when multiple layouts exist:

- `mode_transcriptome_assembly="sraid"` for `workspace/input/query_sra_id`
- `mode_transcriptome_assembly="fastq"` for `workspace/input/species_rnaseq`
- `mode_transcriptome_assembly="metadata"` for `workspace/input/amalgkit_metadata`

Then run:

```bash
cd workflow
bash gg_transcriptome_generation_entrypoint.sh
```

Main outputs:

- `workspace/output/transcriptome_assembly/assembled_transcripts_with_isoforms`
- `workspace/output/transcriptome_assembly/longest_cds`

## 4. Build the species tree, orthogroups, and genome-evolution outputs

Use `gg_genome_evolution_entrypoint.sh` when you already have multispecies CDS inputs and want the core comparative scaffold.
It also supports `input_sequence_mode=protein` when you want the species-tree and orthogroup stages to run from protein sequences only.
In protein mode, place per-species protein FASTA files under `workspace/input/species_protein`, or let GeneGalleon translate `workspace/input/species_cds` with the optional per-species override table `workspace/input/species_genetic_code/species_genetic_code.tsv`.
When `input_sequence_mode=cds`, GeneGalleon ignores `workspace/input/species_protein` and always builds temporary proteins from `workspace/input/species_cds`.
Protein mode automatically disables DNA-tree, IQ2MC/MCMCtree, and BUSCO-based genome-evolution steps that still require CDS inputs.
If you also want OMArk quality assessment, set `run_species_omark=1`; OMArk
uses the same effective protein set as the unified workflow, namely
`workspace/input/species_protein` in protein mode when present, otherwise
temporary proteins translated from `workspace/input/species_cds`.

```bash
cd workflow
bash gg_genome_evolution_entrypoint.sh
```

This is the usual upstream prerequisite for:

- orthogroup-mode gene-family analyses,
- orthogroup database generation,
- convergence analyses that depend on orthogroup outputs.

Main output roots:

- `workspace/output/species_tree`
- `workspace/output/orthofinder`
- `workspace/output/genome_evolution`

Mixed-code example:

```bash
cat > workspace/input/species_genetic_code/species_genetic_code.tsv <<'EOF'
species	genetic_code
Tetrahymena_thermophila	6
EOF
```

Then set `input_sequence_mode="protein"` in `workflow/gg_genome_evolution_entrypoint.sh` and run the wrapper normally.
Species missing from `species_genetic_code.tsv` still use the configured default `genetic_code`.

## 5. Run gene-family analyses in query2family mode

This is the default mode of `gg_gene_evolution_entrypoint.sh`.

Preparation:

- place one query file per family under `workspace/input/query_gene`
- the file basename becomes the family/task ID under `workspace/output/query2family`

Run:

```bash
cd workflow
bash gg_gene_evolution_entrypoint.sh
```

Use this mode when you want to start from hand-picked genes or families rather than all selected orthogroups.
For concrete query-file examples and the resulting per-family output layout, see
[Gene-Family Outputs and Progress Monitoring](gene-family-outputs-and-progress-monitoring.md).

## 6. Run gene-family analyses in orthogroup mode

This mode depends on prior orthogroup selection output from `gg_genome_evolution_entrypoint.sh`.

Edit the top block of `workflow/gg_gene_evolution_entrypoint.sh` so that:

```bash
mode_gene_evolution="orthogroup"
```

Then run:

```bash
cd workflow
bash gg_gene_evolution_entrypoint.sh
```

Important note:

- unlike input generation, this wrapper does not expose a dedicated host-side `GG_*` override map for these mode toggles, so changing the top block is the supported route.

## 7. Build the orthogroup SQLite database

After orthogroup-mode downstream outputs exist:

```bash
cd workflow
bash gg_gene_database_entrypoint.sh
```

Required inputs:

- `workspace/output/orthogroup/stat_tree`
- `workspace/output/orthogroup/stat_branch`

Main output:

- `workspace/output/orthogroup/gg_orthogroup.db`

## 8. Run site-level convergence analysis

After orthogroup outputs and the database are available:

```bash
cd workflow
bash gg_convergent_sites_entrypoint.sh
```

Default prerequisites:

- `workspace/output/orthogroup`
- `workspace/output/orthofinder`
- `workspace/input/species_trait/species_trait.tsv`

Main output root:

- `workspace/output/csubst_site`

## 9. Generate summary TSVs

To create summary tables for completed transcriptome or orthogroup runs:

```bash
cd workflow
bash gg_progress_summary_entrypoint.sh
```

Outputs are written to the workspace root:

- `workspace/orthogroup_summary.tsv`
- `workspace/transcriptome_assembly_summary.tsv`

Current scope:

- orthogroup summary is generated only when `workspace/output/orthogroup`
  exists and the selected gene-count table plus AMAS directories are present,
- transcriptome summary is generated when
  `workspace/output/transcriptome_assembly` exists,
- query2family runs do not currently have a dedicated auto-summary TSV; for
  manual completion-audit commands, see
  [Gene-Family Outputs and Progress Monitoring](gene-family-outputs-and-progress-monitoring.md).

## 10. Full end-to-end order

Not every project needs every step, but the common broad order is:

1. build or pull the container image
2. prepare or collect inputs under `workspace/input/`
3. optionally run `gg_input_generation_entrypoint.sh`
4. optionally run `gg_transcriptome_generation_entrypoint.sh`
5. optionally run `gg_genome_annotation_entrypoint.sh`
6. run `gg_genome_evolution_entrypoint.sh`
7. run `gg_gene_evolution_entrypoint.sh` in the mode you need
8. optionally run `gg_gene_database_entrypoint.sh`
9. optionally run `gg_convergent_sites_entrypoint.sh`
10. optionally run `gg_progress_summary_entrypoint.sh`

## 11. Dry-run and debug the full wrapper chain

For a dependency-aware wrapper sanity check:

```bash
GG_ENTRYPOINT_DRY_RUN=1 bash workflow/gg_all_entrypoints_debug.sh
```

To run only selected steps:

```bash
GG_ENTRYPOINT_ONLY_STEPS=gg_genome_evolution,gg_progress_summary \
bash workflow/gg_all_entrypoints_debug.sh
```

The debug harness writes a run summary to:

- `workspace/output/debug_entrypoint_logs/summary.tsv`
