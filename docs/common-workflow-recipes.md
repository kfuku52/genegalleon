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

## 7. Run matched RSC and species-tree PGLS

RSC tests whether species traits predict gene-expression contrasts while
reconciling paralogous gene-tree events with the species tree. Prepare
`workspace/input/species_expression`,
`workspace/input/species_trait/species_trait.tsv`, and the upstream species
tree, then enable expression generation, matched tree pruning, dating, and the
unified expression-trait stage:

```bash
cd workflow
GG_GENE_EVOLUTION_RUN_GENERATE_EXPRESSION_MATRIX=1 \
GG_GENE_EVOLUTION_RUN_TREE_PRUNING=1 \
GG_GENE_EVOLUTION_RUN_TREE_DATING=1 \
GG_GENE_EVOLUTION_RUN_EXPRESSION_TRAIT_PGLS=1 \
GG_GENE_EVOLUTION_PGLS_METHODS=rsc,species-nwkit,species-rphylopars \
bash gg_gene_evolution_entrypoint.sh
```

`pgls_methods` accepts `rsc`, `species-nwkit`, `species-rphylopars`, any
comma-separated subset, or `all`. `rsc` is the default. All selected methods
use one prepared input and the same direction, `expression ~ species trait`.
The estimands nevertheless differ: RSC models reconciled gene-lineage changes,
whereas both species-tree methods model one expression value per species after
within-biological-sample paralog aggregation. They are complementary analyses,
not interchangeable likelihood implementations.

The default response set and predictor set are both `all`; predictors are fit
separately by default. A parameterized covariance-model example is:

```bash
GG_GENE_EVOLUTION_RUN_EXPRESSION_TRAIT_PGLS=1 \
GG_GENE_EVOLUTION_PGLS_METHODS=all \
GG_GENE_EVOLUTION_RSC_RESPONSES=root,leaf \
GG_GENE_EVOLUTION_RSC_PREDICTORS=habitat,body_size \
GG_GENE_EVOLUTION_RSC_GENE_EVOLUTION_MODEL=lambda \
GG_GENE_EVOLUTION_RSC_GENE_EVOLUTION_PARAMETER=auto \
GG_GENE_EVOLUTION_RSC_SPECIES_EVOLUTION_MODEL=lambda \
GG_GENE_EVOLUTION_RSC_SPECIES_EVOLUTION_PARAMETER=auto \
bash gg_gene_evolution_entrypoint.sh
```

For species-tree PGLS, `species_expression_aggregation` selects `sum`, `mean`,
`max`, or `all`. Aggregation occurs separately within each biological sample,
before replicate uncertainty is estimated, and on the linear expression scale;
the result is transformed back according to `exp_value_type`. Thus paralogs
are never misclassified as biological replicates. The default
`species_paralog_missing=error` also prevents a partly measured paralog set
from silently changing the species-level estimand. A tied `max` uses the
largest standard error among the tied paralogs, avoiding row-order-dependent
uncertainty. For each gene family, the species comparators use the induced
species subtree containing its reconciled gene tips, rather than requiring the
family to occur in every species of the global tree.

Known-SE analyses can also set `species_paralog_sampling_covariance` to a TSV
of `response`, `gene_name_1`, `gene_name_2`, and `sampling_covariance` (plus
optional `tree_id`). NWKIT and Rphylopars then receive the correctly propagated
species-level diagonal uncertainty for `sum` and `mean`; invalid non-positive-
semidefinite combinations are rejected. Raw paired replicates need no such
file because cross-paralog co-variation is preserved in the per-sample
aggregate before replicate variance is estimated.

`auto` estimates the shape parameter for lambda, OU, kappa, delta, EB, and
ACDC models. Brownian and independent models have no shape parameter, so their
default `auto` setting means that no shape-parameter argument is supplied.
With the default original branch-length modes, both dated trees must exist;
use `rsc_gene_branch_length=unit` or `rsc_species_branch_length=unit` only as an
explicit alternative.

GeneGalleon's RSC responses are continuous expression values, so
`rsc_inference` accepts `wald` or `parametric-bootstrap`. Likelihood-ratio and
profile-likelihood options used by other NWKIT response families are not valid
for this integration. Bootstrap inference requires at least two replicates.
GeneGalleon invokes the current `nwkit regress` interface for RSC. The
`rsc_event_weighting` values are `event` (the default, equal total weight per
species-tree event) and `contrast` (equal weight per gene contrast).
`rsc_model=cluster-hc1` selects the earlier cluster-robust estimator; it cannot
propagate response/predictor sampling uncertainty or
automatically estimate a parameterized gene-evolution transform; use
`hierarchical` or `replicate-reml` for those combinations.

NWKIT normally protects large jobs from accidental dense-matrix allocation.
Set `rsc_allow_large_dense=yes` only when the normal sparse/structured route is
unavailable and the host has enough memory; this asks NWKIT to attempt the
allocation rather than guaranteeing that it will fit.

The RSC result is under `rsc_regression/`, and the species comparators are under
`pgls_species_nwkit/` and `pgls_species_rphylopars/`. The long-form
`pgls_comparison/` table places their coefficient rows together while retaining
the method and aggregation. It explicitly labels species-tip estimates as not
directly comparable to the RSC event-level estimand. `pgls_method_status/` and
`pgls_method_audit/` report requested, non-estimable, and successful methods
without silently changing models.

Rphylopars receives the same species means and available diagonal sampling
variances as NWKIT. It does not support categorical predictors in this
comparator, separate response/predictor evolutionary transformations, or full
cross-species sampling covariance. Such fits are recorded as `not_estimable`;
`rphylopars_sampling_covariance=diagonalize` is an explicit opt-in to discard
off-diagonal sampling covariance. Rphylopars 0.3.10 also fails on a singular
`phenocov_list`, so a comparison that mixes positive sampling variances with
exact zero-variance trait/species values is reported as `not_estimable` instead
of applying an artificial variance floor. Coefficients can then be compared, but
log-likelihood, AIC, BIC, optimizer reporting, and parameter counting remain
engine-specific. The adapter implements Wald inference only. If
`rsc_inference=parametric-bootstrap` is requested, the Rphylopars comparator is
reported as `not_estimable`; it is never relabeled as a successful bootstrap
or silently replaced by Wald inference.

The unified stage records the resolved NWKIT revision (and the Rphylopars
package version when requested) in its artifact manifest. Updating a moving
upstream branch therefore invalidates cached RSC/species-comparator outputs
without pinning that branch in GeneGalleon configuration.

The RSC family-level outcome is under
`rsc_status/`, and the full reconciliation/contrast/replicate audit is spread
across the other `rsc_*` directories. The summary step also writes compact
`rsc_*` and bounded `pgls_species_*` columns into `stat_tree/`. See [Input Conventions](input-conventions.md)
for biological/technical replicate and categorical-predictor metadata.
If a family has no expression rows, it receives a complete header-only RSC
bundle and `not_estimable` status instead of aborting the array. A statistical
`ValueError` confined to one predictor analysis is likewise recorded in its
NWKIT audit and the remaining predictor analyses continue; non-statistical
failures remain fatal.

## 8. Run gene-family analyses in orthogroup mode

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

## 9. Build the orthogroup SQLite database

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

To also package analytical-global-q candidates for focused site review, enable
the opt-in candidate-site stage. The summary stage may be enabled in the same
run or its existing min-support TSVs may be reused:

```bash
cd workflow
GG_GENE_SUMMARY_GENE_FAMILY_SOURCE=orthogroup \
GG_GENE_SUMMARY_RUN_CSUBST_SCAN_AA_CHANGE_SUMMARY=1 \
GG_GENE_SUMMARY_RUN_CSUBST_SCAN_CANDIDATE_SITES=1 \
bash gg_gene_summary_entrypoint.sh
```

With the defaults, this writes one ZIP per threshold from the observed maximum
down to 5 in `workspace/output/gene_summary/orthogroup`. Each threshold uses
its own recalculated `q_rate_enrichment_global <= 0.05` selection and includes
focused tree/site reports for the retained rows.

For query2family outputs, switch the gene-family source:

```bash
cd workflow
GG_GENE_SUMMARY_GENE_FAMILY_SOURCE=query2family \
GG_GENE_SUMMARY_RUN_GENE_FAMILY_DATABASE_BUILD=1 \
GG_GENE_SUMMARY_RUN_CSUBST_SCAN_AA_CHANGE_SUMMARY=1 \
bash gg_gene_summary_entrypoint.sh
```

## 10. Run site-level convergence analysis

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

## 11. Generate gene-family summary figures

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

The standard PDF/SVG figure combines the species tree with the gene-family
detected/undetected matrix. In query2family mode, a second
`query2family_reference_gene_orthologs.pdf`/`.svg` figure expands each family
into every gene from the selected reference species. The reference species is
resolved from `GG_COMMON_REFERENCE_SPECIES`; `auto` prefers a supported model
species present in the plotted dataset. Dark cells are orthologs specific to
one reference gene, while a light bar spanning adjacent reference-gene cells
marks a copy inferred to predate their duplication; the printed number is the
copy count. A non-contiguous assignment is shown in orange and identified in
the legend rather than being presented as an ordinary pre-duplication span.
Local-synteny evidence is calculated separately for every
candidate/reference-gene pair. A cell can contain multiple copies, so its
summary uses the highest number of distinct shared flanking-gene similarity
groups while the exact per-copy values remain in `.synteny.tsv`. One shared
anchor is single-anchor evidence, while two or more are classified as supported
in the evidence table. The legend reports the evaluated neighborhood radius as
the number of upstream and downstream gene models inspected on each side.
Missing local-synteny support does not override the reconciled-tree orthology
assignment. Gene-tree UFBoot evidence uses the non-root speciation-MRCA branch
that underlies each candidate/reference assignment. A multi-copy cell uses the
lowest per-copy value and is numeric only when every non-reference copy has an
evaluable branch; missing `support_generax_ufboot` values (or the generic
`support_unrooted` fallback) and root MRCAs make the summary unavailable. Exact
pairwise values and unavailability reasons remain in `.ufboot.tsv`, together
with the support-column source. GeneRax values come from
an unconstrained IQ-TREE UFBoot search whose split frequencies are subsequently
mapped onto the GeneRax topology. This is branch-topology support, not a confidence score for the
reconciliation event classification itself. A
compact reference-gene tree above the matrix is aligned to the
same columns. Duplication nodes in this compact tree remain family-colored
filled circles. On the reconciled species tree, each mapped branch instead gets
a small family-colored bar plot: events from one family on the same branch are
collapsed into one bar, and bar height is proportional to their duplication
count. Bars for different families are spread along that branch. The color key
and compact trees use the input query filename. When the plotted tree is dated,
internal reconciliation labels are transferred by matching clades from the
corresponding undated mapping tree, so mapped events are not lost when its node
labels differ. The
column labels are gene IDs verified against each family's `cds_fasta` output;
the species prefix is omitted for display while the exact CDS FASTA ID remains
in the tables. Multiple families remain separate column blocks, with a family
title above each reference-gene tree and a vertical separator between matrices. The
corresponding `.columns.tsv`, `.glyphs.tsv`, `.tree.tsv`, `.synteny.tsv`, and
`.ufboot.tsv` files record the
reference-gene order, plotted inference, compact tree topology, family-local
duplication index, reconciled species-node mapping, and candidate/reference
local-synteny and gene-tree branch-support evidence. The full TSV matrices are written for all detected families,
while `.plot.*` files record the subset used for the figures. By
default, query2family plots all query files and orthogroup plots the first
100 orthogroups. Use `presence_absence_max_families=0` or `all` to remove the cap,
or set `presence_absence_family_ids=OG0000001,OG0000042` /
`presence_absence_family_file=selected_ogs.txt` to plot an explicit subset in the
requested order.

Detected cells reserve narrow evidence bands along their upper and lower edges
by default. The upper band shows the highest per-copy local-synteny anchor count,
and the lower band shows the lowest per-copy decisive-branch Gene tree UFBoot
support. This leaves the center of a multi-column pre-duplication glyph
continuous and keeps its copy number unobstructed. Both numeric scales map zero
to the same darkest viridis color and their respective maxima to the same
lightest color. A gray slash marks unavailable evidence. Reference-self cells
show white edge bands with a short horizontal state mark; their central
presence/copy-number glyph is unchanged.
Set `presence_absence_evidence_layout=rail`
for the earlier right-edge rail, `glyph` for the diamond/circle overlay, or
`off` to suppress both evidence encodings and their legends. The legacy `rail`
layout also marks the reference gene itself with a white dash. In the
optional `glyph` layout, a viridis diamond shows local-synteny anchors and an
Inferno circle shows Gene tree UFBoot; reference-self and unavailable UFBoot
cells have no circle.

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

## 12. Generate summary TSVs

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

## 13. Full end-to-end order

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

## 14. Dry-run and debug the full wrapper chain

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
