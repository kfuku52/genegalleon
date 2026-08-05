# Configuration and Common Parameters

## Routine configuration

For a persistent project configuration, edit the top config block in each
`workflow/gg_*_entrypoint.sh`:

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

For one-off runs, use the entrypoint-scoped environment overrides described
below. `workflow/core/gg_*_core.sh` files contain the implementation and should
not be edited for routine parameter changes.

## How values are forwarded

Entry-point variables are not forwarded implicitly by name. GeneGalleon uses an explicit registry in:

- `workflow/support/gg_entrypoint_config_vars.sh`

Only variables listed there are eligible for scoped environment overrides and
exported into the container runtime by `forward_config_vars_to_container_env`.

This has two practical consequences:

- adding a new config variable to an entrypoint usually requires adding it to the registry,
- variables that are not in the registry remain host-local unless they are forwarded separately on purpose.

## Configuration precedence

The effective value seen by a core script usually follows this precedence,
from highest to lowest:

1. an entrypoint-scoped environment override,
2. the entrypoint config block, including any `GG_COMMON_*` value referenced by
   that block,
3. a hard-coded fallback inside the core script.

All main entrypoints support scoped overrides. The prefix identifies the
entrypoint, and the registered variable name is converted to uppercase:

- `GG_INPUT_`
- `GG_TRANSCRIPTOME_`
- `GG_GENOME_ANNOTATION_`
- `GG_GENOME_EVOLUTION_`
- `GG_GENE_EVOLUTION_`
- `GG_GENE_SUMMARY_`
- `GG_PROGRESS_SUMMARY_`

For example, `run_cafe` in `gg_genome_evolution_entrypoint.sh` becomes
`GG_GENOME_EVOLUTION_RUN_CAFE`, while `mode_gene_evolution` becomes
`GG_GENE_EVOLUTION_MODE_GENE_EVOLUTION`. Empty override values are supported
when a parameter intentionally needs to be cleared.

## Shared common parameter file

`workflow/gg_common_params.sh` currently defines:

- `GG_COMMON_GENETIC_CODE` (default `1`)
- `GG_COMMON_BUSCO_LINEAGE` (default `auto`)
- `GG_COMMON_REFERENCE_SPECIES` (default `auto`)
- `GG_COMMON_INPUT_SEQUENCE_MODE` (default `cds`)
- `GG_COMMON_CSUBST_NONSYN_RECODE` (default `no`)
- `GG_COMMON_SPECIES_LABEL_PARSER` (default `taxonomic`)
- `GG_COMMON_SPECIES_LABEL_REGEX` (default empty)
- `GG_COMMON_SPECIES_LABEL_MAP_TSV` (default empty)

These are intended for values that recur across multiple stages.

For BUSCO, `GG_COMMON_BUSCO_LINEAGE=auto` resolves a dataset from species names.
For single-species stages, GeneGalleon picks the deepest BUSCO dataset mapped to that species.
For multi-species BUSCO stages, it picks the deepest BUSCO dataset shared across the dataset's species.
In `gg_genome_evolution`, the multi-species BUSCO run and BUSCO summary are shared between the
species-tree branch and the BUSCO-based genome-evolution branch. Those shared stages are controlled
by `run_species_busco` and `run_build_species_busco_summary`; the genome-evolution BUSCO steps reuse
their outputs rather than starting a second BUSCO run.
The same per-species BUSCO short summaries are also used by the two-round
OrthoFinder core selector. By default, core species candidates must satisfy
`busco_complete_pct:ge:80` and `num_seq:le:100000`; these rules are configured
with `orthofinder_core_filters` in `workflow/gg_genome_evolution_entrypoint.sh`.
When BUSCO completeness values are unavailable because the BUSCO stage was
intentionally disabled, GeneGalleon keeps the size filter and logs a fallback
message instead of failing solely on missing BUSCO metadata.
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

Every main entrypoint applies its scoped overrides after loading the editable
config block. Examples include:

```bash
GG_GENOME_EVOLUTION_RUN_CAFE=1 \
GG_GENOME_EVOLUTION_RUN_ORTHOGROUP_COPY_NUMBER_TRAIT_PGLS=1 \
bash workflow/gg_genome_evolution_entrypoint.sh
```

```bash
GG_GENE_EVOLUTION_MODE_GENE_EVOLUTION=orthogroup \
bash workflow/gg_gene_evolution_entrypoint.sh
```

Input generation uses the shorter `GG_INPUT_` prefix. Common overrides include:

- `GG_INPUT_PROVIDER`
- `GG_INPUT_DOWNLOAD_MANIFEST`
- `GG_INPUT_INPUT_DIR`
- `GG_INPUT_RUN_MULTISPECIES_SUMMARY`
- `GG_INPUT_RUN_GENERATE_SPECIES_TRAIT`
- `GG_INPUT_TRAIT_PROFILE`
- `GG_INPUT_GENE_GROUPING_MODE`
- `GG_INPUT_GFF_REPAIR_MODE`
- `GG_INPUT_SPECIES_CDS_DIR`
- `GG_INPUT_SPECIES_GFF_DIR`
- `GG_INPUT_SPECIES_GENOME_DIR`
- `GG_INPUT_SUMMARY_OUTPUT`

Example (no-login GBIF distribution traits):

```bash
GG_INPUT_DOWNLOAD_MANIFEST="$PWD/workspace/input/input_generation/download_plan.xlsx" \
GG_INPUT_TRAIT_PROFILE=gbif_distribution \
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
- one-off stage-local changes via the matching entrypoint-scoped prefix,
- cross-stage defaults in `workflow/gg_common_params.sh`,
- species-tree rooting in `workflow/gg_genome_evolution_entrypoint.sh` via `species_tree_rooting`,
- path relocation via `gg_workspace_dir` / `gg_container_image_path`,
- direct Docker wrapper runs via `GG_CONTAINER_RUNTIME=docker` plus `GG_CONTAINER_DOCKER_IMAGE`.

That keeps routine runs reproducible without forcing edits in the core implementation files.

## CSUBST nonsynonymous-state recoding

`workflow/gg_gene_evolution_entrypoint.sh` exposes `csubst_nonsyn_recode` for
`csubst search --nonsyn_recode`. `workflow/gg_gene_summary_entrypoint.sh`
exposes `csubst_site_nonsyn_recode` for `csubst sites --nonsyn_recode` when
`run_csubst_site_convergence_summary=1`. The shared default is:

```bash
GG_COMMON_CSUBST_NONSYN_RECODE="no"
```

Set it to one of `no`, `3di20`, `dayhoff6`, `sr6`, `kgb6`, `sr4`,
`dayhoff9`, `dayhoff12`, `dayhoff15`, `dayhoff18`, `srchisq6`, or
`kgbauto6` in `workflow/gg_common_params.sh` before running gene-family
analyses. Stage-specific `csubst_nonsyn_recode` or
`csubst_site_nonsyn_recode` overrides still take priority
when supplied for a single entrypoint run.

## CSUBST binary foreground resolution

Set the following opt-in option when a `species_trait.tsv` column uses only
`0` and `1`, but disconnected foreground clades should be treated as distinct
CSUBST lineages:

```bash
csubst_resolve_binary_foreground="yes"
```

The default is `no`. When enabled, GeneGalleon uses the pruned species tree to
find maximal all-foreground clades and assigns them deterministic positive
lineage IDs before both `csubst search` and `csubst scan`. Columns that already
contain a value other than `0` or `1` are treated as manually numbered and are
left unchanged. Missing or unlisted species are treated as background when
clade boundaries are resolved. Summary PDF trait colors use every nonzero
lineage ID as foreground, including IDs of `2` or greater.

## CSUBST scan

`workflow/gg_gene_evolution_entrypoint.sh` also exposes `run_csubst_scan` for
`csubst scan`, a direct recurrent amino-acid/state-change scan that is separate
from the branch-combination search used by `run_csubst`.

Key options:

```bash
run_iqtree_anc=1
run_csubst_scan=1
csubst_scan_unit_mode="clade"
csubst_scan_match="any2spe"
csubst_scan_min_support="2"
csubst_scan_pvalue_calibration="full_scan"
csubst_scan_n_permutations="1000"
csubst_scan_site_plot="yes"
```

`csubst_scan_min_support` preserves the spelling passed to `csubst`: `"1"`
means one foreground unit, fractional values such as `"0.5"` are proportions,
and `"1.0"` means 100% of foreground units.

`csubst_scan_unit_mode` defines the independent foreground support units. The
default `clade` mode automatically splits disconnected foreground regions and
uses every branch within each foreground clade. `stem` uses the same automatic
split but scans only each clade's stem branch. `lineage` uses each positive
foreground lineage ID as one unit and follows `csubst_fg_stem_only`; use
distinct positive IDs in `workspace/input/species_trait/species_trait.tsv` when
independent lineages must remain separate in this mode.

Per-family outputs are written under `csubst_scan/`, `csubst_scan_units/`,
`csubst_scan_foreground_branch/`, `csubst_scan_plot/`, and `csubst_scan_log/`.
Database preparation imports `csubst_scan/*.tsv` into the `aa_change` table and
`csubst_scan_units/*.tsv` into `aa_change_unit`. The `aa_change` table is a
candidate state-change table, not a one-row-per-site table; GeneGalleon adds
global BH-FDR columns such as `q_rate_enrichment_global` after aggregating all
candidate substitutions into the SQLite database. In `gg_gene_summary`, set
`run_csubst_scan_aa_change_summary=1` to write ranked candidate TSV output plus
`*_csubst_aa_change_evidence_density.pdf`,
`*_csubst_aa_change_substitution_spectrum.pdf`, and
`*_csubst_aa_change_foreground_unit_support_matrix.pdf`; combine it with
`run_gene_family_database_build=1` to refresh the database and plots in one run.

GeneGalleon does not require a fixed CSUBST scan column count. Columns are
matched by header name, the union of reported columns is retained, and optional
columns absent from an individual per-family TSV are stored as `NULL` in the
SQLite database. The semantic baseline for this release requires
`site_rate_categorized` and the three empirical q-value columns in
`csubst_scan.tsv`, plus `fg_clade_branch_ids` in
`csubst_scan_units.tsv`; any number of additional columns is accepted without a
GeneGalleon change. Legacy outputs missing those baseline markers are rejected
before the SQLite database is overwritten. When a csubst upgrade changes scan
statistic semantics rather than only adding report columns, regenerate all scan
tables, plots, logs, and foreground-branch outputs together so that the result
set remains internally consistent.

For current scan outputs, `site_rate` is the posterior-mean rate and the former
IQ-TREE categorical rate is retained as `site_rate_categorized`.
`site_rate_quantile` is calculated from posterior-mean rates with tied values
sharing the same average rank. Target and other event/rate quantities cover the
full foreground clades, and q-weighted exposure uses the parent-codon posterior
and codon Q values. Analytical and empirical P/q values therefore should not be
compared directly with values from legacy scan outputs. For 3Di analyses,
`scan_rate_exposure` records the effective `state_aware` model even when the
requested default was `q_weighted`. These provenance and new empirical q-value
columns are retained in the database and the ranked AA-change summary TSV.
