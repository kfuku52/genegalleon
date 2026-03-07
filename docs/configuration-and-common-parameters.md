# Configuration and Common Parameters

## Routine configuration

For normal use, edit the top config block in each `workflow/gg_*_entrypoint.sh`:

```bash
### Start: Modify this block to tailor your analysis ###
...
### End: Modify this block to tailor your analysis ###
```

That is the supported place for stage toggles such as:

- `run_*` flags,
- mode switches,
- thresholds,
- tool selections,
- output-control flags such as temporary-directory cleanup.

`workflow/core/gg_*_core.sh` files contain the implementation and should usually not be edited for per-run parameter changes.

## How values are forwarded

Entry-point variables are not forwarded implicitly by name. GeneGalleon uses an explicit registry in:

- `workflow/support/gg_entrypoint_config_vars.sh`

Only variables listed there are exported into the container runtime by `forward_config_vars_to_container_env`.

This has two practical consequences:

- adding a new config variable to an entrypoint usually requires adding it to the registry,
- variables that are not in the registry remain host-local unless they are forwarded separately on purpose.

## Configuration precedence

The effective value seen by a core script usually comes from this order:

1. entrypoint config block values that are forwarded into the container,
2. stage-specific override helpers such as `GG_INPUT_*` in `gg_input_generation_entrypoint.sh`,
3. shared `GG_COMMON_*` variables,
4. hard-coded fallback values inside core scripts.

Important note:

- `gg_input_generation_entrypoint.sh` has a dedicated host-side override layer using `GG_INPUT_*`.
- other entrypoints do not currently expose a generic `GG_*` override map for all top-block variables.
- for those wrappers, routine changes should be made by editing the top block or maintaining a local wrapper copy.

## Shared common parameter file

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

These are intended for values that recur across multiple stages.

Typical examples:

- one BUSCO lineage reused by transcriptome, annotation, and genome-evolution runs,
- one genetic code reused by annotation and gene-family stages,
- one preferred outgroup reused by species-tree-aware analyses.

## How `GG_COMMON_*` is applied

Shared defaults are loaded in two places:

- host-side bootstrap for entrypoints that opt into `gg_common_params.sh`,
- core-side bootstrap via `gg_source_common_params_from_core`.

They are also forwarded into the container runtime during `set_singularityenv`.

That means you can apply a one-off shared override from the shell, for example:

```bash
GG_COMMON_BUSCO_LINEAGE=metazoa_odb12 \
GG_COMMON_GENETIC_CODE=1 \
bash workflow/gg_genome_evolution_entrypoint.sh
```

Core scripts typically consume these values with parameter expansion such as:

```bash
genetic_code="${genetic_code:-${GG_COMMON_GENETIC_CODE:-1}}"
busco_lineage="${busco_lineage:-${GG_COMMON_BUSCO_LINEAGE:-embryophyta_odb12}}"
outgroup_labels="${outgroup_labels:-${GG_COMMON_OUTGROUP_LABELS:-Oryza_sativa}}"
```

## Entry-point override patterns

### Dedicated override map: input generation

`gg_input_generation_entrypoint.sh` explicitly supports host-side overrides such as:

- `GG_INPUT_PROVIDER`
- `GG_INPUT_DOWNLOAD_MANIFEST`
- `GG_INPUT_INPUT_DIR`
- `GG_INPUT_RUN_GENERATE_SPECIES_TRAIT`
- `GG_INPUT_TRAIT_PROFILE`
- `GG_INPUT_SPECIES_CDS_DIR`
- `GG_INPUT_SPECIES_GFF_DIR`
- `GG_INPUT_SPECIES_GENOME_DIR`
- `GG_INPUT_SUMMARY_OUTPUT`

Example:

```bash
GG_INPUT_DOWNLOAD_MANIFEST="$PWD/workspace/input/input_generation/download_plan.xlsx" \
GG_INPUT_TRAIT_PROFILE=gift_starter \
bash workflow/gg_input_generation_entrypoint.sh
```

### Shared runtime paths

Advanced path knobs are handled outside the per-entrypoint config registry:

- `gg_workspace_dir`
- `gg_container_image_path`

These are resolved before the container starts and are useful when:

- you want to keep a workspace outside the repository,
- the SIF lives in a different location,
- multiple workspaces share the same checked-out code.

## Recommended configuration strategy

Use the following split in practice:

- stage-local flags in the entrypoint block,
- cross-stage defaults in `workflow/gg_common_params.sh`,
- one-off input-generation automation via `GG_INPUT_*`,
- path relocation via `gg_workspace_dir` / `gg_container_image_path`.

That keeps routine runs reproducible without forcing edits in the core implementation files.
