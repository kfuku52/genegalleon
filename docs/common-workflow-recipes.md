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

For no-login GBIF distribution traits:

```bash
GG_INPUT_TRAIT_PROFILE=gbif_distribution \
bash workflow/gg_input_generation_entrypoint.sh
```

Useful outputs to inspect afterward:

- `workspace/output/input_generation/download_plan.resolved.tsv`
- `workspace/output/input_generation/gg_input_generation_runs.tsv`
- `workspace/output/input_generation/gg_input_generation_species.tsv`
- `workspace/output/input_generation/annotation_summary/annotation_summary.tsv`

When both CDS and GFF are available for a species, also inspect the adjacent
`*.fa.gz.gff-grouping.tsv` audit and the `cds_gff_*` columns in
`gg_input_generation_species.tsv` before promoting formatted inputs into a
project workspace.

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

Use `gg_genome_evolution_entrypoint.sh` when you already have multispecies CDS
inputs and want the core comparative scaffold.
CDS FASTA sequence IDs must already follow the `GENUS_SPECIES_GENEID`
convention documented in [Input Conventions](input-conventions.md).
It also supports `input_sequence_mode=protein` when you want the species-tree
and orthogroup stages to run from protein sequences only.
Most CDS-first runs do not need `workspace/input/species_protein`; GeneGalleon
can translate `workspace/input/species_cds` into temporary proteins when
needed. In protein mode, place per-species protein FASTA files under
`workspace/input/species_protein` only when curated/native proteins should be
used directly, or let GeneGalleon translate `workspace/input/species_cds` with
the optional per-species override table
`workspace/input/species_genetic_code/species_genetic_code.tsv`.
Providing correctly translated `species_protein` files can also include
lineages with different genetic codes, but codon-sequence-based analyses are
not available from protein-only inputs.
When `input_sequence_mode=cds`, GeneGalleon ignores
`workspace/input/species_protein` and always builds temporary proteins from
`workspace/input/species_cds`.
Protein mode automatically disables DNA-tree, IQ2MC/MCMCtree, and BUSCO-based
genome-evolution steps that still require CDS inputs.
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

## 5. Test phenotype associations with orthogroup copy number

After the species tree, orthogroup gene-count table, and
`workspace/input/species_trait/species_trait.tsv` are available, enable the
orthogroup copy-number trait PGLS stage from the genome-evolution wrapper:

```bash
cd workflow
GG_GENOME_EVOLUTION_RUN_CAFE=0 \
GG_GENOME_EVOLUTION_RUN_ORTHOGROUP_COPY_NUMBER_TRAIT_PGLS=1 \
GG_GENOME_EVOLUTION_ORTHOGROUP_COPY_NUMBER_TRAIT="all" \
bash gg_genome_evolution_entrypoint.sh
```

This writes the shared copy-number matrix to:

- `workspace/output/genome_evolution/orthogroup_copy_number/orthogroup_copy_number.tsv`

Trait-PGLS outputs are written under:

- `workspace/output/genome_evolution/orthogroup_copy_number/trait_pgls/`

Use `orthogroup_copy_number_trait="trait_a,trait_b"` to test selected trait
columns, and use `orthogroup_copy_number_trait_family_ids` or
`orthogroup_copy_number_trait_family_file` to restrict the orthogroups tested.
This stage tests species-level orthogroup copy numbers against species traits
with species-tree PGLS. It does not use CAFE5 result files such as
`Gamma_change.tab`; `run_cafe=1` is only needed when you also want the separate
CAFE family-size evolution analysis.

## 6. Run gene-family analyses in query2family mode

This is the default mode of `gg_gene_evolution_entrypoint.sh`.

Preparation:

- place one query file per family under `workspace/input/query_gene`
- the file basename becomes the family/task ID under `workspace/output/query2family`
- keep unrelated gene families in separate files; if query genes from different
  families are combined in one file, their homologs and phylogenetic trees can
  be artificially merged into a single unnatural family tree.

Run:

```bash
cd workflow
bash gg_gene_evolution_entrypoint.sh
```

Use this mode when you want to start from hand-picked genes or families rather than all selected orthogroups.
For concrete query-file examples and the resulting per-family output layout, see
[Gene-Family Outputs and Progress Monitoring](gene-family-outputs-and-progress-monitoring.md).

## 7. Run gene-family analyses in orthogroup mode

This mode depends on prior orthogroup selection output from `gg_genome_evolution_entrypoint.sh`.

Set the mode persistently in the entrypoint's top block, or override it for one
run:

```bash
cd workflow
GG_GENE_EVOLUTION_MODE_GENE_EVOLUTION=orthogroup \
bash gg_gene_evolution_entrypoint.sh
```

Important note:

- scoped environment variables are intended for one-off runs; edit the top
  config block when the mode is part of the persistent project configuration.

## 8. Build the orthogroup SQLite database

After orthogroup-mode downstream outputs exist:

```bash
cd workflow
GG_GENE_SUMMARY_GENE_FAMILY_SOURCE=orthogroup \
GG_GENE_SUMMARY_RUN_GENE_FAMILY_DATABASE_BUILD=1 \
GG_GENE_SUMMARY_RUN_CSUBST_SCAN_AA_CHANGE_SUMMARY=1 \
bash gg_gene_summary_entrypoint.sh
```

Required inputs:

- `workspace/output/orthogroup/stat_tree`
- `workspace/output/orthogroup/stat_branch`

Main output:

- `workspace/output/orthogroup/gg_orthogroup.db`
- `workspace/output/gene_summary/orthogroup/orthogroup_csubst_aa_change_min_support_2_summary.tsv`
- `workspace/output/gene_summary/orthogroup/orthogroup_csubst_aa_change_min_support_*.pdf`
- `workspace/output/gene_summary/orthogroup/orthogroup_csubst_aa_change_min_support_manifest.tsv`

When
`workspace/output/orthofinder/Orthogroups_filtered/Orthogroups.GeneCount.annotated.tsv`
is available, every orthogroup `min_support_*_summary.tsv` also includes the
five representative `besthit_0.05` through `besthit_0.95` annotation columns.

For query2family outputs, switch the gene-family source:

```bash
cd workflow
GG_GENE_SUMMARY_GENE_FAMILY_SOURCE=query2family \
GG_GENE_SUMMARY_RUN_GENE_FAMILY_DATABASE_BUILD=1 \
GG_GENE_SUMMARY_RUN_CSUBST_SCAN_AA_CHANGE_SUMMARY=1 \
bash gg_gene_summary_entrypoint.sh
```

## 9. Run site-level convergence analysis

After orthogroup outputs and the database are available:

```bash
cd workflow
GG_GENE_SUMMARY_GENE_FAMILY_SOURCE=orthogroup \
GG_GENE_SUMMARY_RUN_CSUBST_SITE_CONVERGENCE_SUMMARY=1 \
bash gg_gene_summary_entrypoint.sh
```

Default prerequisites:

- `workspace/output/orthogroup`
- `workspace/output/orthofinder`
- `workspace/input/species_trait/species_trait.tsv`

Main output root:

- `workspace/output/csubst_site`

The same post-analysis can be run for either gene-family mode by setting
`gene_family_source` and `run_csubst_site_convergence_summary=1`.

## 10. Generate gene-family summary figures

To combine query2family outputs into a species-tree presence/absence summary:

```bash
cd workflow
GG_GENE_SUMMARY_GENE_FAMILY_SOURCE=query2family \
bash gg_gene_summary_entrypoint.sh
```

Main output root:

- `workspace/output/gene_summary/query2family`

For orthogroups, switch the mode:

```bash
cd workflow
GG_GENE_SUMMARY_GENE_FAMILY_SOURCE=orthogroup \
bash gg_gene_summary_entrypoint.sh
```

The PDF/SVG figure combines the species tree with the gene-family
detected/undetected matrix. The full TSV matrices are written for all detected
families, while `.plot.*` files record the subset used for the figure. By
default, query2family plots all query files and orthogroup plots the first
100 orthogroups. Use `presence_absence_max_families=0` or `all` to remove the cap,
or set `presence_absence_family_ids=OG0000001,OG0000042` /
`presence_absence_family_file=selected_ogs.txt` to plot an explicit subset in the
requested order.

If the selected species tree is dated and
`workspace/output/species_tree/mcmctree_main/mcmctree_95CI.nwk` exists, the
figure adds node-age 95% HPD bars. Numeric branch-support labels are
transferred from available species-tree support files when their splits match
the plotted tree. When BUSCO full tables or a BUSCO summary table are
available, the figure also adds right-side per-species BUSCO stacked bars;
full tables include Fragmented counts. Override the auto-selected
files with `presence_absence_species_tree`, `presence_absence_species_tree_ci`,
`presence_absence_species_tree_support`, or `presence_absence_busco_table` when needed.
See [Example Plots](example-plots.md) for a compact generated figure.

## 11. Generate summary TSVs

To create summary tables for completed transcriptome, orthogroup, or
query2family runs:

```bash
cd workflow
bash gg_progress_summary_entrypoint.sh
```

Outputs are written to the workspace root:

- `workspace/orthogroup_summary.tsv`
- `workspace/query2family_summary.tsv`
- `workspace/transcriptome_assembly_summary.tsv`

Current scope:

- orthogroup summary is generated when `workspace/output/orthogroup`
  exists and the selected gene-count table is present; AMAS inputs are
  optional and may be live or ZIP-backed,
- query2family summary is generated when `workspace/output/query2family`
  and `workspace/input/query_gene` exist,
- transcriptome summary is generated when
  `workspace/output/transcriptome_assembly` exists.

## 12. Full end-to-end order

Not every project needs every step, but the common broad order is:

1. build or pull the container image
2. prepare or collect inputs under `workspace/input/`
3. optionally run `gg_input_generation_entrypoint.sh`
4. optionally run `gg_transcriptome_generation_entrypoint.sh`
5. optionally run `gg_genome_annotation_entrypoint.sh`
6. run `gg_genome_evolution_entrypoint.sh`
7. optionally enable `run_orthogroup_copy_number_trait_pgls=1` in
   `gg_genome_evolution_entrypoint.sh` for phenotype/copy-number association
   tests
8. run `gg_gene_evolution_entrypoint.sh` in the mode you need
9. optionally run `gg_gene_summary_entrypoint.sh`
10. optionally run `gg_progress_summary_entrypoint.sh`

## 13. Dry-run and debug the full wrapper chain

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
