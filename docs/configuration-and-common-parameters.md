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
- `GG_COMMON_BUSCO_LINEAGE` (default `auto`)
- `GG_COMMON_REFERENCE_SPECIES` (default `auto`)

These are intended for values that recur across multiple stages.

For BUSCO, `GG_COMMON_BUSCO_LINEAGE=auto` resolves a dataset from species names.
For single-species stages, GeneGalleon picks the deepest BUSCO dataset mapped to that species.
For multi-species BUSCO stages, it picks the deepest BUSCO dataset shared across the dataset's species.
In `gg_genome_evolution`, the multi-species BUSCO run and BUSCO summary are shared between the
species-tree branch and the BUSCO-based genome-evolution branch. Those shared stages are controlled
by `run_species_busco` and `run_build_species_busco_summary`; the genome-evolution BUSCO steps reuse
their outputs rather than starting a second BUSCO run.
When BUSCO publishes multiple `odbN` generations, auto-resolution now uses the latest generation
for which placement mappings are available across archaea, bacteria, and eukaryota.
The first auto-resolved run may need network access to initialize the ETE taxonomy DB and download
BUSCO placement mapping files; explicit values such as `embryophyta_odb13` still bypass that logic.

For contamination removal, `contamination_removal_rank` is now configured locally in
`workflow/gg_genome_annotation_entrypoint.sh` and `workflow/gg_transcriptome_generation_entrypoint.sh`.
GeneGalleon treats `domain` as the canonical user-facing value and normalizes tool-specific
synonyms automatically (for example, `remove_contaminated_sequences.py` receives `superkingdom`).
When the sample species name is unknown but you still know the host clade, set the local
`contamination_removal_target_taxon` parameter to an NCBI-recognized taxon name such as
`Eukaryota`; the contamination-removal step will use that lineage anchor instead of the
directory or filename-derived species label.

For annotation-driven stages, `GG_COMMON_REFERENCE_SPECIES=auto` prefers model species detected
in the relevant dataset and falls back to the first available species when none of the preferred
models are present. The current priority list keeps only `Arabidopsis_thaliana` and `Oryza_sativa`
on the plant side, then checks standard cross-clade model species such as human, mouse, zebrafish,
fly, nematode, yeasts, and `Escherichia_coli`. Tree-visualization ortholog prefixes derive from
that species name downstream, so the shared common variable no longer includes a trailing underscore.

`species_tree_rooting`, `grampa_h1`, and `target_branch_go` are no longer shared `GG_COMMON_*` values.
They are now configured directly in `workflow/gg_genome_evolution_entrypoint.sh`
as genome-evolution-local parameters. `species_tree_rooting` defaults to
`taxonomy` there and accepts forms such as:

- `taxonomy`
- `taxonomy,ncbi`
- `taxonomy,ncbi,opentree,timetree`
- `outgroup,Oryza_sativa`
- `outgroup,Oryza_sativa,Amborella_trichopoda`
- `midpoint`
- `mad`
- `mv`

For backward compatibility, a bare species-label list such as
`Oryza_sativa,Amborella_trichopoda` is still interpreted as
`outgroup,Oryza_sativa,Amborella_trichopoda`, but the explicit `outgroup,...`
form is preferred for new configs.

When `grampa_h1` or `target_branch_go` are left empty, GeneGalleon skips only the
GRAMPA-related steps or the GO-enrichment step, respectively.

For duplicate-aware BUSCO genome-evolution steps, the canonical config names are
the `run_busco_dupaware_*` flags exposed in
`workflow/gg_genome_evolution_entrypoint.sh`, for example:

- `run_busco_dupaware_extract_fasta`
- `run_busco_dupaware_iqtree_dna`
- `run_busco_dupaware_notung_root_pep`
- `run_busco_dupaware_grampa_dna`

All duplicate-aware BUSCO substeps default to `0`. `run_orthogroup_grampa`
defaults to `1`, but it is still auto-disabled unless
rooted orthogroup trees are present and `grampa_h1` is non-empty.

Typical examples:

- one auto-resolved or explicit BUSCO lineage reused by transcriptome, annotation, and genome-evolution runs,
- one genetic code reused by annotation and gene-family stages,
- one annotation species reused by GO-enrichment and tree-visualization steps.

## Mixed genetic code datasets

`GG_COMMON_GENETIC_CODE` and the local `genetic_code` parameter still act as a single default code.
That default is used when:

- a stage only accepts one code for the whole run,
- `gg_genome_evolution` is translating `species_cds` and no per-species override is present,
- a species is missing from `workspace/input/species_genetic_code/species_genetic_code.tsv`.

`gg_genome_evolution_entrypoint.sh` now adds a separate switch:

- `input_sequence_mode="cds"`: normal CDS-first behavior
- `input_sequence_mode="protein"`: run species-tree and orthogroup stages from protein inputs

For mixed-code projects, the intended setup is:

1. keep `GG_COMMON_GENETIC_CODE` or local `genetic_code` as the fallback default,
2. optionally provide `workspace/input/species_genetic_code/species_genetic_code.tsv`,
3. run `gg_genome_evolution_entrypoint.sh` with `input_sequence_mode="protein"`.

In that mode, GeneGalleon prefers `workspace/input/species_protein` when present.
Most projects should leave `species_protein` absent unless curated or native
protein FASTA files should be used directly. If `species_protein` is absent, it
translates `workspace/input/species_cds` to temporary proteins, applying
per-species overrides from `species_genetic_code.tsv` first and the global
default code second.
Providing correctly translated `species_protein` files is another way to include
lineages with different genetic codes, because GeneGalleon does not translate
CDS in that path. The trade-off is that codon-sequence-based analyses are not
available from protein-only inputs.
DNA-tree and dating steps that still require CDS-only assumptions are disabled
automatically in protein mode.

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
busco_lineage="${busco_lineage:-${GG_COMMON_BUSCO_LINEAGE:-auto}}"
annotation_species="${annotation_species:-${GG_COMMON_REFERENCE_SPECIES:-auto}}"
```

For `gg_genome_evolution`, remember that this fallback is only part of the final translation rule.
The effective CDS-to-protein code priority there is:

1. `workspace/input/species_genetic_code/species_genetic_code.tsv` for matching species
2. local `genetic_code=...` in `workflow/gg_genome_evolution_entrypoint.sh`
3. `GG_COMMON_GENETIC_CODE`
4. core fallback `1`

## Entry-point override patterns

### Dedicated override map: input generation

`gg_input_generation_entrypoint.sh` explicitly supports host-side overrides such as:

- `GG_INPUT_PROVIDER`
- `GG_INPUT_DOWNLOAD_MANIFEST`
- `GG_INPUT_INPUT_DIR`
- `GG_INPUT_RUN_MULTISPECIES_SUMMARY`
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
- `GG_CONTAINER_RUNTIME`
- `GG_CONTAINER_DOCKER_IMAGE`

These are resolved before the container starts and are useful when:

- you want to keep a workspace outside the repository,
- the SIF lives in a different location,
- multiple workspaces share the same checked-out code,
- or you want to run wrappers directly against a Docker image instead of a SIF.

Docker-backed wrapper mode can still be enabled explicitly:

```bash
GG_CONTAINER_RUNTIME=docker \
GG_CONTAINER_DOCKER_IMAGE=ghcr.io/kfuku52/genegalleon:latest \
bash workflow/gg_gene_evolution_entrypoint.sh
```

When the default repo-root `genegalleon.sif` is missing, wrappers also
auto-fallback to a pulled Docker image if available. Current fallback priority
is `ghcr.io/kfuku52/genegalleon:latest`, then `local/genegalleon:dev`.

When `GG_CONTAINER_RUNTIME=docker` is set, keep `gg_container_image_path`
reserved for SIF-based runs; use `GG_CONTAINER_DOCKER_IMAGE` for the Docker image reference.

## Recommended configuration strategy

Use the following split in practice:

- stage-local flags in the entrypoint block,
- cross-stage defaults in `workflow/gg_common_params.sh`,
- species-tree rooting in `workflow/gg_genome_evolution_entrypoint.sh` via `species_tree_rooting`,
- one-off input-generation automation via `GG_INPUT_*`,
- path relocation via `gg_workspace_dir` / `gg_container_image_path`,
- direct Docker wrapper runs via `GG_CONTAINER_RUNTIME=docker` plus `GG_CONTAINER_DOCKER_IMAGE`.

That keeps routine runs reproducible without forcing edits in the core implementation files.
