# Configuration and Common Parameters

## Routine configuration

For compatible species-tree calibration auditing, reviewed calibration inputs and
opt-in prior-only/sensitivity runs, see [species-tree calibrations](species-tree-calibrations.md).

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

HGT summaries automatically return category-1 trait results with
`run_hgt_trait_focus=1`. `hgt_summary_focus_event_tsv=auto` and
`hgt_summary_focus_event_gene_tsv=auto` select native event/context tables;
explicit paths preserve a previously filtered project cohort and enriched gene
annotations as inputs. Focused results additionally apply the configured
[shared query-Pfam pair filter](gene-structure-tree-plot.md#shared-pfam-filter-for-focused-hgt).
Pfam filtering precedes trait selection, so multiple traits share the same
passing cohort and Pfam-stage counts.
For non-Arthropoda donors into Insecta, set
`hgt_summary_focus_direction_filter=non_arthropoda_to_insecta`. This additional
species-branch filter follows the Pfam pair filter and precedes every trait's
category-1 selection. `hgt_summary_focus_species_taxonomy=auto` reads the saved
`output/species_taxonomy/species_taxonomy.tsv`; an explicit path selects another
existing host-species taxonomy table. The general default `any` preserves other
projects' transfer directions. See the root `direction_event_audit.tsv` and
`direction_species_branches.tsv` for excluded, unknown and mixed branches.
See [focused HGT outputs](host-scaffold-taxonomy.md#category-1-focused-results)
for the three summary PDFs, two-page gene-tree/context PDFs, trait eligibility,
internal-branch context and event-counting rules.
`hgt_summary_focus_filter_audit_tsv` optionally supplies an existing project
direction/UFBoot audit for upstream filtering-flow counts.

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

The editable entrypoint blocks are the source of generated configuration
metadata. Run `bash ./dev config-check` to detect forwarding-registry drift,
or `bash ./dev config-schema markdown` to render a current reference table.

## Shared common parameter file

`workflow/gg_common_params.sh` currently defines:

- `GG_COMMON_GENETIC_CODE` (default `1`)
- `GG_COMMON_BUSCO_LINEAGE` (default `eukaryota_odb12`; set `auto` to infer a lineage from species names)
- `GG_COMMON_REFERENCE_SPECIES` (default `auto`)
- `GG_COMMON_INPUT_SEQUENCE_MODE` (default `cds`)
- `GG_COMMON_REPRESENTATIVE_INPUTS` (default empty; verified CDS/protein/GFF/genome manifest from [gene-model refinement](gene-model-refinement.md))
- `GG_COMMON_CSUBST_NONSYN_RECODE` (default `no`)
- `GG_COMMON_SPECIES_LABEL_PARSER` (default `taxonomic`)
- `GG_COMMON_SPECIES_LABEL_REGEX` (default empty)
- `GG_COMMON_SPECIES_LABEL_MAP_TSV` (default empty)
- `GG_COMMON_GENE_FAMILY_OUTPUT_STORAGE` (default `zip`; `zip`, `files`, or the `raw` alias for `files`; see [large-array ZIP collection](gene-family-array-archive-queue.md))
- `GG_COMMON_GENE_FAMILY_ZIP_MIN_BATCH_FILES` (default `100`; deprecated compatibility setting, no longer used by array-task cleanup)
- `GG_COMMON_GENE_FAMILY_ZIP_COMPRESSION` (default `adaptive`; `adaptive`, `deflate`, or `store`)
- `GG_COMMON_GENE_FAMILY_ZIP_COMPRESSION_LEVEL` (default `6`; `0` through `9`)
- `GG_COMMON_GENE_FAMILY_ZIP_WORKERS` (default `1`; `1` through `4`)
- `GG_COMMON_GENE_FAMILY_LARGE_ZIP_WARNING_BYTES` (default `21474836480`, 20 GiB; `0` disables conversion-report warnings)
- `GG_COMMON_GENE_FAMILY_FINAL_ZIP_MAX_BYTES` (default `0`, a final ZIP may have any size; a positive value retains named part ZIPs above that logical size)
- `GG_COMMON_GENE_FAMILY_TMP_RETENTION_DAYS` (default `7`; `0` disables the age limit)
- `GG_COMMON_GENE_FAMILY_TMP_MAX_DIRS` (default `100`; `0` disables the failed-directory count limit)
- `GG_COMMON_GENE_FAMILY_TMP_MAX_BYTES` (default `107374182400`, 100 GiB; `0` disables the byte limit)
- `GG_COMMON_GENE_FAMILY_TMP_MAX_FILES` (default `100000`; `0` disables the file-count limit)
- `GG_COMMON_SPECIES_TREE_OUTPUT_STORAGE` (default `zip`; `zip`, `files`, or the `raw` alias for `files`)
- `GG_COMMON_SPECIES_TREE_ZIP_COMPRESSION` (default `adaptive`; `adaptive`, `deflate`, or `store`)
- `GG_COMMON_SPECIES_TREE_ZIP_COMPRESSION_LEVEL` (default `6`; `0` through `9`)

These are intended for values that recur across multiple stages.

The gene-family storage setting applies to query2family and orthogroup
artifacts only. ZIP mode leaves eligible artifacts in standard ZIP shards
while downstream summaries and selected consumers use a logical live-plus-ZIP
view. `files` and its `raw` alias keep the previous one-file-per-artifact behavior.
Failed-task work directories remain below each gene-family output root rather
than using a node-wide system temporary directory. In ZIP mode, cleanup runs
only for families whose family lock is idle and removes task directories
older than the configured retention period. It also enforces directory,
aggregate-byte, and file-count limits while preferentially retaining the
newest failed-task directories, so a burst of failures is bounded without
blocking cleanup for unrelated active families.

The species-tree storage setting applies only to the five high-file-count
single-copy stage directories documented in
`docs/species-tree-stage-zip-storage.md`. Small summary, concatenated-tree,
ASTRAL, and MCMCTree directories remain directly visible for cross-workflow
consumers.

For BUSCO, the conservative default is `GG_COMMON_BUSCO_LINEAGE=eukaryota_odb12`.
Setting `GG_COMMON_BUSCO_LINEAGE=auto` explicitly resolves a dataset from species names.
For single-species stages, GeneGalleon picks the deepest BUSCO dataset mapped to that species.
For multi-species BUSCO stages, it picks the deepest BUSCO dataset shared across the dataset's species.
Placement mappings are resolved from BUSCO's standard `file_versions.tsv`,
using the latest integer ODB version available for all three domains. Dated
archives are checked against the manifest before the existing mapping-ready
stamp is published. A directory URL is not used as a listing, and acquisition
errors remain visible when no local mapping can be reused.
In `gg_genome_evolution`, the multi-species BUSCO run and BUSCO summary are shared between the
species-tree branch and the BUSCO-based genome-evolution branch. Those shared stages are controlled
by `run_species_busco` and `run_build_species_busco_summary`; the genome-evolution BUSCO steps reuse
their outputs rather than starting a second BUSCO run.
The same per-species BUSCO short summaries are also used by the two-round
OrthoFinder core selector. By default, core species candidates must satisfy
`busco_complete_pct:ge:80` and `num_seq:le:100000`; these rules are configured
with `orthofinder_core_filters` in `workflow/gg_genome_evolution_entrypoint.sh`.
`orthofinder_binary` selects a complete native executable. Its path and executable
content participate in inference provenance. For an external source runtime,
set `orthofinder_source_manifest` to its immutable source inventory; changing that
inventory invalidates cached inference under the configured stale-artifact policy.
An explicitly configured missing inventory stops before inference.
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
that species name downstream, and query2family reference-gene orthology plots
use every family-tree tip assigned to the resolved species. The shared common
variable therefore does not include a trailing underscore.

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
native MUL-tree steps or the GO-enrichment step, respectively. The legacy
GRAMPA setting names now select `nwkit mul-reconcile`; see
[replacement and output compatibility](grampa-replacement.md).

For duplicate-aware BUSCO genome-evolution steps, the canonical config names are
the `run_busco_dupaware_*` flags exposed in
`workflow/gg_genome_evolution_entrypoint.sh`, for example:

- `run_busco_dupaware_extract_fasta`
- `run_busco_dupaware_iqtree_dna`
- `run_busco_dupaware_reconciliation_root_pep`
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

MCMCtree's derived standalone NHX/Newick sidecars follow NWKit's default
`--rooting-token no` and `--rooting-nhx no` output policy. This omits declaration
tokens without changing the root, topology, branch lengths, names or support.
The original FigTree/NEXUS declarations and dated tree remain unchanged during
sidecar recovery. Explicit NWKit `--rooting-token yes` emits the rooting token;
the ON/OFF outputs must have the same interpreted tree semantics.
Conversions that cannot retain the interpreted rooting state remain rejected.
WGD classification output follows the same default-OFF token policy; its
writer/reader round trip must preserve rooted topology, lengths and node
annotations for both leading-token and NHX-rooted inputs.

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

Gene-tree rooting keeps MAD as the default. The selectable
`tree_rooting_method` values are `mad`, `reconciliation`, `midpoint`,
and `md` (`md` maps to NWKIT's `mv` method). `reconciliation` uses NWKIT's
duplication/loss-assisted rooting with the pruned species tree and the configured
species-label parser, regular expression, or mapping TSV. It does not invoke
NOTUNG. When GeneRax is enabled, the reconciliation-rooted input is preserved
with GeneRax's `--enforce-gene-tree-root`; other rooting modes retain the
existing GeneRax MAD-rooting behavior. For a one-off run:

```bash
GG_GENE_EVOLUTION_TREE_ROOTING_METHOD=reconciliation \
bash workflow/gg_gene_evolution_entrypoint.sh
```

When GeneRax is enabled, post-GeneRax UFBoot is calculated from an
**unconstrained** IQ-TREE bootstrap search. GeneGalleon then counts those
replicate-tree splits on the GeneRax target topology; it does not use the fully
resolved GeneRax tree as an IQ-TREE topology constraint. The resulting
percentages are stored in `stat_branch` as `support_generax_ufboot`.
Identical sequences are retained in the replicate trees so their tip set
remains identical to the GeneRax target.
`treevis_support_value="auto"` prefers that column, falls back to
`support_unrooted`, and suppresses node labels when neither contains values.
An explicit column name or `no` can still be configured instead of `auto`.

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

Example (no-login GBIF observation traits):

```bash
GG_INPUT_DOWNLOAD_MANIFEST="$PWD/workspace/input/input_generation/download_plan.xlsx" \
GG_INPUT_TRAIT_PROFILE=gbif_distribution \
bash workflow/gg_input_generation_entrypoint.sh
```

These are summaries of retained observations, not estimates of the true species
range. Keep the generated metadata and quality sidecars. See
[GBIF observation traits](gbif-observation-traits.md) for local downloads, filters,
explicit analysis selection, sensitivity replay and migration from older columns.


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

## HGT summary and transfer-tree visualization

Genome annotation's `run_scaffold_taxonomy=1` reuses available raw CDS taxonomy
and GFF info to generate host-scaffold composition; it does not enable either
upstream step. HGT evaluation consumes this context automatically. See
[host-scaffold taxonomy](host-scaffold-taxonomy.md) for rank-specific fractions,
candidate-free background, missing-data semantics, and counting units.

`workflow/gg_gene_summary_entrypoint.sh` exposes the HGT summary settings below:

- `hgt_summary_species_tree` (default `auto`): Newick species tree used to map
  directed HGT links; `auto` searches the standard workspace species-tree
  outputs.
- `hgt_summary_species_trait` (default `auto`): display numeric/binary tip traits
  from `workspace/input/species_trait/species_trait.tsv` when present. An explicit
  path selects another table; `none` disables the panel. Shared schema/metadata
  contracts apply. Missing observations remain NA; no ancestral states are inferred.
- `hgt_summary_transfer_tree_max_edges` (default `200`): number of mapped
  donor-to-recipient edges drawn in `plots/hgt_transfer_tree.pdf`, selected by
  alternating event-count and tree-distance rankings. Existing reverse
  directions are then added as separate arrows, so the directional count can
  exceed this initial limit; `0` draws
  all mapped edges. `plots/hgt_transfer_edges.tsv` always retains all parsed
  directed pairs.
- `hgt_summary_transfer_arrow_alpha` (default `0.55`): transfer-arrow opacity
  from `0` (transparent) to `1` (opaque), used by the summary and trait-focused
  figures. Translucent arrows make overlapping transfers easier to inspect.
- `hgt_summary_focus_require_shared_pfam` (default `1`): category-1 focused
  tables and figures require at least one retained event-linked donor/recipient
  gene pair sharing a query Pfam accession, with each gene individually passing
  the existing class scaffold-background profile. Set `0` for the earlier
  scaffold-only cohort. The input and global candidate tables are preserved.
- `hgt_summary_focus_min_shared_pfam_coverage` (default `0.5`): inclusive
  fraction required on **both proteins of the same event-gene pair**. Union
  overlapping saved query-domain intervals before dividing by that protein's
  length in amino acids. This measures shared-domain query coverage, not
  pairwise alignment coverage. Set `0` for any shared query Pfam. Review flags
  for repeat/generic binding domains, differing domain sets, >2-fold query
  length differences and proteins <100 aa do not exclude events.
- `hgt_summary_focus_allow_both_no_pfam` (default `0`): when the Pfam filter is
  enabled, set `1` to also allow a pair whose two genes both have explicit
  searched/no-hit records. Missing records and one-sided no-hits still fail.
  This explicit opt-in bypasses unavailable domain coverage; coverage stays
  unmeasured. No additional E-value threshold is imposed on saved Pfam hits.
  Saved query RPS-BLAST records are used; no additional search is run. Per-trait
  event/pair/gene audits retain exclusion reasons and source hashes, including
  when plotting is disabled.
- `hgt_summary_focus_direction_filter` (default `any`): set
  `non_arthropoda_to_insecta` to require modeled donor species branches entirely
  outside Arthropoda and recipients entirely within Insecta, after Pfam and
  before category-1 traits. Mixed, unknown and unmapped branches are withheld.
- `hgt_summary_focus_species_taxonomy` (default `auto`): existing host-species
  taxonomy for the direction filter, normally
  `output/species_taxonomy/species_taxonomy.tsv`. An enabled filter requires
  this table even without plotting.

The transfer plot counts branch-level GeneRax `Y@donor@recipient` records and
sets each arrow's constant shaft width in proportion to that direction's event
count. Reciprocal directions use separate curves. Arrowheads point to the
recipient/target. This is a count visualization, not a probability or HGT
confidence score.

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

For `3di20`, see [full-CDS inputs, model resources and report coordinates](csubst-3di.md).

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

`workflow/gg_gene_evolution_entrypoint.sh` exposes `run_csubst_scan` for
recurrent amino-acid state changes. GeneGalleon runs **analytical P only**:
`--scan_pvalue_calibration none --scan_n_permutations 0` is fixed in the core.
No candidate-fixed, full-scan maxT or parametric-bootstrap P values are computed.
The previous calibration/permutation configuration variables are no longer
forwarded. The separate arity-based `run_csubst` workflow is unchanged.

```bash
run_csubst_scan=1
csubst_scan_unit_mode="clade"
csubst_scan_match="any2spe"
csubst_scan_min_support="2"
csubst_scan_site_plot="yes"
```

`csubst_scan_min_support` preserves CSUBST's spelling: `"1"` means one unit,
`"0.5"` is a proportion, and `"1.0"` means 100%. `lineage`, `stem` and `clade`
unit modes, event threshold, control scope and nonsynonymous recoding remain
configurable. Event mass, exposure and branch-length scale use CSUBST's defaults;
the former `csubst_scan_rate_event_mode`, `csubst_scan_rate_exposure` and
`csubst_scan_rate_length` settings are no longer forwarded. The scan consumes the existing
IQ-TREE ancestral-reconstruction archive and foreground table.

### Analytical P and global BH-FDR

Database preparation imports scan candidates into `aa_change` and support
units into `aa_change_unit`. A candidate row is a state-change hypothesis,
not necessarily a unique site or orthogroup. It retains CSUBST's
`score_rate_enrichment`, `p_rate_enrichment_asymptotic`, trait × match diagnostic
q values and inference metadata, then calculates:

- `q_rate_enrichment_asymptotic_global`: Benjamini–Hochberg correction of
  `p_rate_enrichment_asymptotic` over **all finite candidate P values in this
  database**, pooling orthogroups, traits and requested match classes.
- `aa_change_fdr_metadata`: the correction method, source/destination columns,
  family definition, candidate count, finite test count and undefined count.

For ordered P values, BH uses `min(1, min_{j>=i}(m * p[j] / j))`, where `m`
is the number of finite tested candidate rows. It does not use the number of
orthogroups as the denominator. Undefined P values remain undefined; malformed,
infinite or out-of-range inputs stop the database build rather than being
clipped. Empty scan inputs retain a table schema and zero-test metadata.

These are **nominal BH-FDR estimates**. Their error-rate interpretation depends
on the validity of the analytical P model, candidate selection and dependence
assumptions. CSUBST's fractional posterior event masses and same-data candidate
selection are not repaired by BH. The correction covers the imported candidate
set, not unreported hypotheses, additional databases or settings explored later.

### Summary and candidate reports

In `gg_gene_summary_entrypoint.sh`, set `run_gene_family_database_build=1` and
`run_csubst_scan_aa_change_summary=1` to rebuild the database and create the
`*_csubst_aa_change_min_support_2_summary.tsv` and its support, substitution
spectrum and P/FDR distribution PDFs. Ranking uses global BH-FDR, breaking ties
with the rate score. The legacy `min_support_2` table retains all imported
candidates; `*_all_candidates_summary.tsv` also preserves that complete table.
Post-hoc support bounds can be changed by regenerating summaries from the
existing database, without rerunning scan. Candidates discarded during scan
discovery cannot be recovered this way.

When present in CSUBST scan output, `lineage_total`, `support_lineage_count`,
`support_lineage_fraction` and `support_lineage_ids` are retained in summaries,
candidate TSVs and manifests, and shown in candidate reports. These count
support grouped by the nonzero IDs in the existing foreground table. Paralogs
or species assigned to one phenotypic origin can share an ID; the grouped
count then counts that ID once across disconnected clade units. The total
includes only IDs with analyzable candidate branches, and the fraction divides
support by that total. This grouping does not infer phenotypic origins or
change the original P values or full-candidate global FDR. Older scan tables
without these columns remain readable when the lineage support condition is
disabled. See
[CSUBST's grouped-support definition](https://github.com/kfuku52/csubst/blob/master/docs/SCAN_INFERENCE.md#support-grouped-by-foreground-lineage-id).

The five `besthit_*` annotations are joined by orthogroup from the annotated
Orthogroups gene-count table when available. They propagate to the filtered
support views. Query2family summaries do not require these annotations.

For each integer support threshold from 3 to the observed maximum, the summary
also writes `*_min_support_<N>_summary.tsv` and its probability plot, plus
`*_min_support_manifest.tsv`. Each view calculates
`q_rate_enrichment_asymptotic_support_filtered` from the analytical P values
**after filtering on unit support** (lineage bound 0). The original P and
`q_rate_enrichment_asymptotic_global` remain unchanged as diagnostic columns.

To require both unit and species/foreground-lineage support in summary views:

```bash
csubst_scan_summary_min_unit_support=6
csubst_scan_summary_min_lineage_support=4
```

The defaults are 2 and 0; 0 disables the corresponding condition. Both bounds
are inclusive and combined with AND. `support_unit_count` counts the units
selected by scan's `unit_mode`: gene-tree clades in `clade` mode, already grouped
foreground IDs in `lineage` mode. The species/lineage bound always uses
`support_lineage_count`, not a new inference of phenotypic origins. The configured
bounds add `*_min_unit_support_<G>_min_lineage_support_<S>_summary.tsv` and three
PDFs; the unit threshold increases from G to the surviving maximum at the fixed
lineage bound S. A matching `_manifest.tsv` lists these views and both bounds.
Each (G, S) view calculates `q_rate_enrichment_asymptotic_support_filtered`
by applying BH to all finite analytical P values surviving **both** bounds,
pooling orthogroups, traits and match classes. It ranks by this filtered q value.
The manifest and candidate rows record the bounds, total candidate count, finite
test count and undefined count. Missing P values remain undefined and do not
increase the finite-test denominator. The correction family changes with the
bounds; these q values do not correct for trying multiple bound settings.
With G=0, one view is written without a unit condition.
Existing unit-only tables and the complete candidate table remain available.
An enabled lineage condition requires a valid count in every imported candidate;
missing columns or counts report the source and affected rows before replacing
summary outputs. Rebuild the database from scans containing those counts when
older and newer results have been mixed.

Candidate-site reports are opt-in:

```bash
run_csubst_scan_candidate_sites=1
csubst_scan_candidate_sites_min_support=5
csubst_scan_candidate_sites_min_lineage_support=0
csubst_scan_candidate_sites_probability_column="q_rate_enrichment_asymptotic_support_filtered"
csubst_scan_candidate_sites_probability_threshold="0.05"
csubst_scan_candidate_sites_max_candidates=0
csubst_scan_candidate_sites_pdb="none"
```

The default recomputes support-filtered analytical BH-FDR for the candidate
report's own bounds, then selects q <= 0.05. A preexisting filtered q column in
the source is overwritten using the source P values and these actual bounds.
Analytical P, full-candidate global q or CSUBST's within-trait × match analytical
q column may be selected explicitly to retain their original selection behavior.
Empirical and bootstrap columns are not accepted by this helper. For an explicitly
selected source q column, missing q never falls back to P.
Thresholds are visited from observed maximum down to the
configured minimum; a zero candidate cap retains all qualifying rows.
Candidate bounds are independent of the summary-view settings and use the same
inclusive AND rule. Set `csubst_scan_candidate_sites_min_lineage_support=4` to
require four foreground lineage IDs. Set the existing unit bound to 0 to create
one ZIP without a unit condition; bounds 0 or 1 require the complete candidate
table emitted by the updated summary helper. Older unit-only summary series
remain usable with bounds >= 2 and no lineage requirement; enabling a lineage
bound requires complete lineage counts in the candidate source tables. All
source tables are validated before existing archives or manifests are replaced.
For support-filtered BH, every threshold uses the complete candidate table;
if it is unavailable, the broadest legacy table covering the requested unit
bound supplies the family. Stale narrower views cannot remove tests from BH.
BH precedes the probability cutoff, candidate cap and skips for missing report
inputs. All surviving finite tests enter the family, including nonsignificant
rows and rows that cannot be packaged. The candidate cap is applied after
support and probability selection; it never changes the BH denominator.

Each threshold produces a ZIP such as
`<source>_csubst_aa_change_candidate_sites_min_support_<N>_q_rate_enrichment_asymptotic_support_filtered_le_0.05.zip`.
It contains manifests, candidate source rows, raw `csubst sites` outputs,
focused tree/site plots and combined reports. Shared candidate analyses are
cached once across thresholds; source summaries, selection parameters and
required input signatures govern archive reuse. Missing report inputs are
recorded as skipped candidates. With a lineage bound enabled, ZIPs and run
manifests also carry `_min_lineage_support_<S>` in their names. Both bounds are
recorded in package metadata, candidate TSVs/manifests and PDF annotations;
archive reuse verifies them even when the surviving candidates are identical.
Package metadata and run manifests record the BH policy and family counts,
including empty selections. Filtered-q archive reuse also verifies that policy
and family, the selected probability column/threshold and candidate q values.
Older archives remain reusable for explicitly selected source P/q
columns with the lineage condition disabled; their absent lineage bound is 0.
`pdb="besthit"` enables optional structure
searching. The report can be computationally expensive even without scan
permutations.

### Outputs and migration

Per-family outputs are written under `csubst_scan/`, `csubst_scan_units/`,
`csubst_scan_foreground_branch/`, `csubst_scan_plot/` and `csubst_scan_log/`.
`csubst_scan_audit/<family>_csubst_scan_audit.zip` preserves the complete scan
output, JSON inference definitions/disabled-calibration record and copies of
its input alignment, topology and foreground table. It is verified and
atomically published before scratch cleanup, including for empty scans.
Source execution paths may need relocation when reproducing a saved run.

Current score/analytical-P/analytical-q/testability/inference marker columns are
required; optional columns are unioned by header name. Unsupported legacy scan
schemas are rejected before replacing a database. Regenerate scan outputs,
then the database and summaries together. New provenance contracts invalidate
old statistical outputs.

Remove old `csubst_scan_pvalue_calibration`, `csubst_scan_n_permutations` and
`csubst_scan_permutation_seed` settings. Candidate settings formerly named
`*_q_column` / `*_q_threshold` now use `*_probability_column` /
`*_probability_threshold`; the helper CLI and manifests follow those names.
The legacy `q_rate_enrichment_global` is replaced by the explicit
`q_rate_enrichment_asymptotic_global` name.

See [integration validation](reviews/2026-09-10-csubst-scan-integration.md) for
execution coverage and the remaining analytical-model limitations.

For retired NOTUNG switches, candidate outputs, and reconciliation statistics,
see [NOTUNG replacement](notung-replacement.md).

## CDSKIT localization model

Both genome annotation and gene evolution default to `cdskit_localize_model=latest`.
On the first download, GeneGalleon queries the official CDSKIT GitHub releases and
selects the most recently published non-draft, non-prerelease `localize-*` release
with one `cdskit-localize-*.pt` checkpoint. It verifies the release asset's SHA-256
and saves the checkpoint and `selection.json` under
`workspace/downloads/cdskit_models/genegalleon-latest/` (or under
`$CDSKIT_MODEL_DIR/genegalleon-latest/` when overridden).
Subsequent runs reuse that selection without checking for newer releases, even
if the container or CDSKIT version changes. The cache lock serializes concurrent
first downloads. Failed downloads do not save a selection. A missing or corrupt
selected checkpoint is an error, not a trigger to select a different model.

The first automatic selection requires network access. With
`cdskit_localize_no_model_download=1` or `CDSKIT_OFFLINE=1`, a saved selection is
required. Explicit CDSKIT aliases or model paths remain supported and bypass
this automatic selection. Existing alias-specific caches are kept separately;
they do not pin the first `latest` selection. Model-specific backbone weights
(such as ESM2) remain managed by CDSKIT and may also need an initial download.
A newly published checkpoint must be supported by the installed CDSKIT runtime;
GeneGalleon does not silently substitute an older model on load failure.

## Genetic code in CSUBST sites

New IQ-TREE ancestral-state bundles contain `csubst.input.json` with the integer
`genetic_code` used by ASR. Ordinary bundles use schema
`genegalleon-csubst-input-v1`; full-CDS 3Di bundles use
`genegalleon-csubst-input-v2` and additionally record their structural inputs.
Both convergent-sites reports and scan-candidate sites resolve the code for each
family and pass it as `csubst sites --genetic_code`. For older bundles, the
resolver reads explicit genetic-code statements or `--seqtype CODONn`/`-st CODONn`
from `csubst.log` and `csubst.iqtree`. Conflicting or absent evidence is an error;
regenerate the ASR bundle rather than assuming code 1. Existing complete reports
are retained under the existing artifact reuse rules.
