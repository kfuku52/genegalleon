# Main Stages and What They Do

### `gg_input_generation_entrypoint.sh`

Purpose:

- optional pre-stage to automate input preparation,
- optionally download provider files from a manifest and format into
  `workspace/output/input_generation/species_cds` and
  `workspace/output/input_generation/species_gff` and
  `workspace/output/input_generation/species_genome`,
- optionally generate `workspace/input/species_trait/species_trait.tsv` using
  target species from an auto-discovered `download_plan` manifest or `species_cds`.

Main scripts:

- `workflow/core/gg_input_generation_core.sh`
- `workflow/support/format_species_inputs.py`
- `workflow/support/generate_species_trait.py`

Notable defaults:

- when `GG_INPUT_DOWNLOAD_MANIFEST` is unset, exactly one non-empty CSV, TSV,
  or XLSX `download_plan*` file with `provider` and `id` columns is selected;
  multiple valid candidates stop the run and require an explicit manifest path,
- `run_validate_inputs=1` validates CDS naming rules, including the required
  `GENUS_SPECIES_GENEID` sequence-ID pattern, CDS/GFF species-set consistency,
  and CDS-to-GFF mapping compatibility, with species-level mapping checks
  parallelized via `validate_cds_gff_mapping.py --nthreads`,
- `run_species_busco=1` runs BUSCO on formatted CDS inputs by default,
- `run_multispecies_summary=1` generates BUSCO plots and `annotation_summary.tsv` under `workspace/output/input_generation/annotation_summary/`,
- formatted outputs default to `workspace/output/input_generation/species_cds`, `workspace/output/input_generation/species_gff`, and `workspace/output/input_generation/species_genome`,
- when a species has both provided CDS and GFF inputs, CDS isoforms are grouped
  through the GFF parent/alias graph before the longest CDS per resolved gene is
  selected; `gene_grouping_mode=rescue_overlap` can rescue compatible fragmented
  models, while `strict` preserves the provider's model boundaries,
- GFF-backed grouping writes adjacent `*.gff-grouping.json` and
  `*.gff-grouping.tsv` audits and records mapped, unmapped, ambiguous, and
  coordinate-rescued counts in `gg_input_generation_species.tsv`,
- per-run metadata defaults to `workspace/output/input_generation/gg_input_generation_runs.tsv`,
  `workspace/output/input_generation/gg_input_generation_species.tsv`, and
  `workspace/output/input_generation/download_plan.resolved.tsv`,
- trait generation is opt-in (`run_generate_species_trait=0` by default),
  and writes `workspace/input/species_trait/species_trait.tsv` when enabled,
- `trait_profile=gift_starter` is available as a quick preset:
  sets `run_generate_species_trait=1` and uses `trait_databases=gift`
  when `trait_databases` is unset or `auto`,
- `trait_profile=gbif_distribution` is available as a no-login GBIF occurrence
  preset: sets `run_generate_species_trait=1`, uses `trait_databases=gbif`
  when `trait_databases` is unset or `auto`, and writes built-in distribution
  traits such as northern latitude limit, breadth, and occupied grid area,
- each run appends a TSV record to `workspace/output/input_generation/gg_input_generation_runs.tsv`
  (configurable via `summary_output` or `GG_INPUT_SUMMARY_OUTPUT`).

Wrapper-specific notes:

- the wrapper accepts dedicated host-side overrides such as `GG_INPUT_DOWNLOAD_MANIFEST`,
  `GG_INPUT_PROVIDER`, `GG_INPUT_TRAIT_PROFILE`, and `GG_INPUT_SPECIES_*`,
- `provider=refseq` and `provider=genbank` are accepted by the wrapper as aliases of `ncbi`.

`input_generation_mode` semantics:

- `single`: one process runs formatting, validation, per-species BUSCO, optional trait generation, and the final multi-species BUSCO summary.
- `array_prepare`: prepares downloads if needed, discovers species tasks, and writes `workspace/output/input_generation/tmp/task_plan.json`; this is the setup step before scheduler array workers run.
- `array_worker`: each array task reads one row from `task_plan.json` using `GG_ARRAY_TASK_ID`, formats one species, and optionally runs BUSCO for that species; outputs are written as shard files under `workspace/output/input_generation/tmp/`.
- `array_finalize`: a single follow-up run merges shard outputs, validates the merged species set, checks BUSCO counts, then runs the shared post-processing steps such as trait generation and `run_multispecies_summary`.

Practical rule:

- use `single` unless you explicitly want scheduler array parallelism across species.

Trait generation inputs:

- trait plan TSV (default: `workspace/input/input_generation/trait_plan.tsv`)
  - required columns: `database`, `source_column`, `output_trait`
  - optional columns: `value_type` (`numeric|binary|categorical`),
    `aggregation`, `positive_values` (comma-separated, for binary mapping),
    `trait_key`, `trait_key_column`
  - for `gift`, `trait_key` can be either a trait ID (`Lvl3`, e.g. `1.1.1`)
    or trait name (`Trait2`, e.g. `Woodiness_1`)
- database source map TSV (default: `workspace/input/input_generation/trait_database_sources.tsv`)
  - required column: `database`
  - optional columns: `acquisition_mode` (`bulk|species_api|gift_api`), `uri`,
    `species_column`, `delimiter`, `response_format`,
    `archive_member`, `trait_key_column`,
    `gift_version`, `gift_trait_ids`, `gift_page_size`,
    `gift_max_pages_per_trait`, `gift_bias_ref`, `gift_bias_deriv`,
    `gift_agreement_min`, `gift_versions_api`
- starter templates are bundled at:
  - `workspace/input/input_generation/trait_plan.tsv`
  - `workspace/input/input_generation/trait_database_sources.tsv`
  - bundled defaults are prewired for `austraits`, `eltontraits`, and
    `gift` (`gift_api` with official base endpoint), including starter
    trait names for `woodiness` (`Woodiness_1`, ID `1.1.1`),
    `plant_height_max_m` (`Plant_height_max`, ID `1.6.2`), and
    `seed_mass_mean_g` (`Seed_mass_mean`, ID `3.2.3`).

Trait DB retrieval policy in `generate_species_trait.py`:

- `bulk`: fetch/cache full dataset under `workspace/downloads/trait_datasets/<database>/`
  then subset target species.
  - `uri` may contain multiple comma-separated sources; all are merged.
  - ZIP sources are supported (`archive_member` to select file inside archive).
- `species_api`: query per target species and merge responses.
- `gift_api`: resolve target species `work_ID` via GIFT
  `names_matched_unique`, resolve trait tokens from `trait_key` /
  `gift_trait_ids` (ID or name via `traits_meta`), then query trait pages
  and filter to target species.

Supported DB IDs can be listed with:

```bash
python workflow/support/generate_species_trait.py --print-supported-databases
```

GIFT traits catalog can be inspected with:

```bash
python workflow/support/generate_species_trait.py --print-gift-traits --gift-trait-search wood --gift-trait-limit 20
```

Quick summary of recent gg_input_generation runs:

```bash
python workflow/support/summarize_gg_input_generation_runs.py \
  --infile workspace/output/input_generation/gg_input_generation_runs.tsv \
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

- `assembly_method="auto"` (selects an available supported assembler; set `rnaSPAdes` or `Trinity` to force one)
- `kallisto_reference="longest_cds"`
- staged FASTA outputs are written as `.fa.gz`.

### `gg_genome_annotation_entrypoint.sh`

Purpose:

- per-species CDS/genome annotation and QC,
- BUSCO (CDS/genome),
- UniProt annotation (`blastp` or `mmseqs2`),
- optional `cdskit localize` targeting-peptide and peroxisome localization prediction,
- optional MMseqs2 taxonomy and contamination removal,
- optional genome analyses (SubPhaser, dotplot, GenomeScope).

Main outputs:

- `workspace/output/species_cds_annotation`
- `workspace/output/species_cds_cdskit_localize`
- `workspace/output/species_cds_busco_full`, `species_cds_busco_short`
- `workspace/output/species_genome_busco_full`, `species_genome_busco_short`

Notable defaults:

- most heavy tasks default to `0`,
- `uniprot_annotation_method="mmseqs2"` (set `blastp` to use NCBI BLASTP for UniProt annotation),
- `cdskit_localize_organism_group="auto"` infers plant/non-plant mode from `busco_lineage` where possible,
- `run_multispecies_summary=1` by default.

### `gg_fractionation_bias_entrypoint.sh`

Runs fully local pairwise fractionation-bias and within-genome self-synteny
retention analyses with
[`kfFractBias`](https://github.com/kfuku52/kfFractBias). Each scheduler array
task selects one data row from
`workspace/input/fractionation_bias_pairs.tsv`, resolves that row's target and
query files from `workspace/input/species_cds` and
`workspace/input/species_gff`. Rows with `mode=compare` run `kffractbias
compare`; rows with `mode=self` run `kffractbias selfcompare`. The analysis
uses local JCVI MCscan, LAST (or BLAST), and QUOTA-ALIGN; it does not contact
CoGe or another analysis service at runtime.

Required TSV columns:

- `analysis_id`: unique output-safe identifier,
- `target_species` and `query_species`: GeneGalleon species labels matching
  the corresponding CDS and GFF filenames,
- `quota`: explicit target:query syntenic depth such as `1:2`.

`mode` is optional and defaults to `compare`. For `mode=self`, set the same
species in `target_species` and `query_species` and use a symmetric quota such
as `1:1` or `2:2`. Its common value is passed to `selfcompare --depth`. The
target GFF override and `target_seqids` are used for the shared genome;
`query_feature`, `query_attribute`, and `query_seqids` must be empty.

Optional columns and defaults are:
`mode=compare`, `window_size=100`, `step_size=1`, `denominator=all`, empty
`target_seqids`, empty `query_seqids`, empty `exclude_seqid_regex`,
`cscore=0.7`, `aligner=last`, empty target/query GFF feature and attribute
overrides, `minimum_mapping_fraction=0.5`, `self_hit_percent=98`, and
`diagonal_bound=300`. The last two values apply only to self mode and control
JCVI's near-self percent-identity filtering and minimum intrachromosomal
gene-rank distance from the self diagonal. Feature and attribute overrides must
be supplied as a pair. A complete template is available at
`workspace/input/fractionation_bias_pairs.example.tsv`.

For a row whose `analysis_id` is `sorghum_maize`, outputs are written below
`workspace/output/kffractbias/sorghum_maize/`:

- `sorghum_maize.genes.tsv`: one row per analyzed target gene and query
  sequence combination,
- `sorghum_maize.windows.tsv`: sliding-window retention statistics,
- `sorghum_maize.summary.json`: parameters, counts, hashes, and metadata,
- `sorghum_maize.plot.pdf` and `sorghum_maize.plot.png`,
- `sorghum_maize.synteny.zip`: the local JCVI working data needed to inspect
  or reproduce the retained synteny result.

Self rows publish the same six-file bundle below
`workspace/output/genome_evolution/self_fractionation_bias/<analysis_id>/`.
Their summaries label the result as `self_synteny_retention`: it is a
within-genome retention-asymmetry measurement conditional on extant annotated
genes, not an outgroup-based estimate of ancestral losses.

Set each scheduler's array range to the number of data rows in the TSV. For a
single local row, the default is sufficient:

```bash
cp workspace/input/fractionation_bias_pairs.example.tsv \
  workspace/input/fractionation_bias_pairs.tsv
bash workflow/gg_fractionation_bias_entrypoint.sh
```

`artifact_stale_policy` applies in the same way as other GeneGalleon stages;
the six published outputs are replaced as one recoverable transaction.

### `gg_genome_evolution_entrypoint.sh`

Purpose:

- unified genome-evolution entrypoint that serially runs:
  - inlined species-tree stage (formerly `gg_speciesTree_core.sh`)
  - inlined orthofinder stage (formerly `gg_orthofinder_core.sh`)
  - inlined genome-evolution stage (formerly `gg_genomeEvolution_core.sh`)

The optional `run_self_fractionation_bias=1` stage incorporates and aggregates
completed `mode=self` array results. Heavy LAST/JCVI jobs remain scheduler-array tasks in
`gg_fractionation_bias_entrypoint.sh`; the unified genome-evolution job only
validates and summarizes them. The default configuration table is
`workspace/input/fractionation_bias_pairs.tsv` and can be overridden with
`self_fractionation_bias_table`.

The aggregation writes:

- `workspace/output/genome_evolution/self_fractionation_bias_summary/self_fractionation_bias_summary.tsv`,
- matching PDF and PNG summary plots.

The headline metric takes the best interchromosomal retention track for each
target window, then averages those window values. Intrachromosomal retention is
reported separately so a residual same-sequence track cannot inflate the main
within-genome cross-sequence summary.
- preserves output-exists skip behavior at step level and aborts on real failures.
- accepts either CDS-first input (`input_sequence_mode=cds`) or protein-only
  input (`input_sequence_mode=protein`); protein-only input can include
  correctly translated proteins from lineages with different genetic codes, but
  codon-sequence-based analyses are unavailable.

Main output roots:

- `workspace/output/species_tree`
- `workspace/output/orthofinder`
- `workspace/output/genome_evolution`

Additional optional QC outputs:

- `workspace/output/genome_evolution/omamer_search`
- `workspace/output/genome_evolution/omark`
- `workspace/output/genome_evolution/omark_summary_table/omark_summary.tsv`

### Inlined Stage: Orthofinder

Purpose:

- translate CDS to species proteins, or reuse `workspace/input/species_protein` in protein mode,
- run OrthoFinder,
- select orthogroups for downstream analysis,
- compare orthogroup methods,
- plot single-copy ortholog decay as line plots with SD across random species subsamples.

Main outputs:

- `workspace/output/orthofinder`
- `workspace/output/orthofinder/single_copy_ortholog_decay/single_copy_ortholog_decay_plot.pdf`
- `workspace/output/orthofinder/single_copy_ortholog_decay/single_copy_ortholog_decay_summary.tsv`

Temporary protein FASTA files are created under `workspace/downloads/tmp/` and removed automatically after orthogroup-related steps finish.
When `workspace/input/species_genetic_code/species_genetic_code.tsv` is present, it overrides the global `genetic_code` on a per-species basis during CDS-to-protein translation; species missing from the table still use the default `genetic_code`.

Notable defaults:

- `orthogroup_table="HOG"`
- `orthogroup_annotation_method="mmseqs2"` (set `blastp` to use NCBI BLASTP for representative-gene UniProt annotation)
- `run_single_copy_ortholog_decay_plot=1`
- `orthogroup_decay_replicates=1000`
- species-tree-aware OrthoFinder is used when species tree exists.
- before OrthoFinder starts, protein-input species and species-tree leaves are
  parsed with the same configured species-label parser, regex, and mapping.
  Exact parsed labels are matched first; remaining labels may match by the
  parser's taxonomy query only when that mapping is one-to-one. True omissions,
  duplicate parsed labels, and ambiguous genus-level matches stop before
  OrthoFinder.
- when the species count exceeds `max_orthofinder_core_species`, the core
  species set is selected with size and BUSCO filters by default:
  `orthofinder_core_filters="busco_complete_pct:ge:80,num_seq:le:100000"`.
  If a species tree is available, GeneGalleon calls `nwkit sample` with
  `orthofinder_core_method="max-pd"` to keep phylogenetic diversity in the
  OrthoFinder core run; otherwise it falls back to ranked selection.
- OrthoFinder is reported as successful only after its whole-dataset output is
  validated. GeneGalleon accepts the current complete `Orthogroups.tsv` layout,
  falls back to the root-level `Phylogenetic_Hierarchical_Orthogroups/N0.tsv`
  when required by OrthoFinder 3.1+ output differences, and rejects clade-level
  `N*.tsv` files as substitutes for the root table.

### Inlined Stage: Species Tree

Purpose:

- BUSCO-based single-copy extraction,
- optional OMArk proteome quality assessment using the effective protein set
  (provided `species_protein` inputs in protein mode, otherwise temporary
  proteins translated from `species_cds`),
- per-gene alignments and trees,
- concatenated and ASTRAL species trees,
- IQ2MC/mcmctree dating pipeline,
- constrained-tree plotting.

Main outputs:

- `workspace/output/species_tree/species_tree_summary/undated_species_tree.nwk`
- `workspace/output/species_tree/species_tree_summary/dated_species_tree.nwk`
- BUSCO and intermediate tree/alignment directories under `species_tree/`.

Notable defaults:

- species-tree substeps are enabled in the unified wrapper by default, including BUSCO extraction, ASTRAL tree inference, and the IQ2MC/mcmctree dating steps,
- the same multi-species BUSCO run also feeds the BUSCO-based genome-evolution branch; BUSCO-based genome-evolution steps reuse the BUSCO outputs driven by `run_species_busco` and `run_build_species_busco_summary`,
- OMArk is optional (`run_species_omark=0` by default) and runs after
  OrthoFinder so it can reuse the current effective protein inputs,
- shared defaults such as `busco_lineage` and `genetic_code` are loaded from `workflow/gg_common_params.sh`,
- species-tree rooting is configured locally in `workflow/gg_genome_evolution_entrypoint.sh` via `species_tree_rooting`,
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
- optional localization tables under `workspace/output/query2family/cdskit_localize/`
  or `workspace/output/orthogroup/cdskit_localize/`.
- optional `csubst scan` outputs under `csubst_scan/`, `csubst_scan_units/`,
  `csubst_scan_foreground_branch/`, `csubst_scan_plot/`, and `csubst_scan_log/`.

Notable defaults:

- `mode_gene_evolution="query2family"`
- `run_rps_blast=1` (Pfam domain annotation is on by default),
- `run_tree_plot=1`
- `run_summary=1`
- `uniprot_annotation_method="mmseqs2"` (set `blastp` for NCBI BLASTP-based UniProt annotation),
- `cdskit_localize_organism_group="auto"` infers plant/non-plant mode from `busco_lineage` where possible,
- many advanced analyses default to `0`.

Current behavior notes:

- if `run_generax=1`, initial IQ-TREE disables UFBOOT and support is
  recomputed after GeneRax on the GeneRax topology,
- Pfam RPS-BLAST DB (`Pfam_LE`) is auto-prepared when missing, with lock-based
  synchronization for array jobs,
- gene-tree PGLS remains a separate legacy analysis; matched expression-trait
  results are written to `rsc_pgls`, `pgls_species_nwkit`,
  `pgls_species_rphylopars`, and `pgls_comparison`,
- `run_expression_trait_pgls=1` runs the methods selected by `pgls_methods`
  (`rsc`, `species-nwkit`, `species-rphylopars`, or `all`) with gene expression
  as the response and species traits as predictors; repeated paralog contrasts
  mapping to one species-tree event are handled by event weighting and the
  hierarchical model instead of being treated as independent species
  observations,
- species-tree methods aggregate paralogs within each biological sample before
  replicate modeling. `species_expression_aggregation` controls sum, mean,
  max, or all, and `species_paralog_missing=error` prevents incomplete copy
  coverage from being silently accepted,
- RSC uses a dated gene tree and dated species tree with the default
  `original` branch-length modes. Set the corresponding mode to `unit` only
  when the dated tree is intentionally unavailable. `rsc_event_source="auto"`
  uses preserved GeneRax NHX speciation/duplication annotations when present
  and otherwise performs LCA reconciliation,
- RSC records families with missing/empty expression and statistically
  non-estimable predictor fits with explicit reasons rather than aborting the
  family array. Its inspectable tables are written under the
  `rsc_*` output directories. Species-method statuses, aggregation audits, and
  side-by-side coefficients are written under `pgls_*`; compact `rsc_*` and
  `pgls_species_*` fields are copied into each
  family's `stat_tree` row. This summary has a fixed-width status/count schema
  plus the best converged `rsc_best_*` row; all fitted rows remain in
  `rsc_pgls`,
- the compact `stat_tree` best-p fields are multiplicity-aware: raw minima are
  retained under `*_p_value_raw`, while unsuffixed p-values use Holm correction
  over all usable associations for the family and method; BH-adjusted values
  and the tested-association count/scope are reported alongside them,
- `run_csubst_scan=1` uses existing CSUBST ancestral-reconstruction inputs and
  writes candidate amino-acid/state changes; it is independent from
  `run_csubst`, which runs branch-combination search.
- `csubst scan` columns are matched by name rather than position or count;
  additional columns and per-file optional-column differences are retained,
  with unavailable values stored as `NULL`. Legacy outputs missing the current
  semantic marker columns are rejected during database preparation.

Wrapper-specific note:

- `gg_gene_evolution_entrypoint.sh` forwards registered top-block variables to
  the container runtime,
- persistent changes belong in the top config block,
- one-off changes can use the `GG_GENE_EVOLUTION_` prefix, for example
  `GG_GENE_EVOLUTION_MODE_GENE_EVOLUTION=orthogroup`.

### Inlined Stage: Genome Evolution

Purpose:

- BUSCO-based and orthogroup-based GRAMPA workflows,
- optional CAFE, orthogroup copy-number trait PGLS, and GO enrichment analyses.

Main outputs:

- `workspace/output/genome_evolution/*`

Notable defaults:

- duplicate-aware BUSCO genome-evolution substeps default to `0`,
- `run_orthogroup_grampa=1`, but GRAMPA is auto-disabled unless rooted
  orthogroup trees exist and `grampa_h1` is set,
- `run_cafe=0`, `run_orthogroup_copy_number_trait_pgls=0`, and
  `run_go_enrichment=0` by default,
- `orthogroup_copy_number/orthogroup_copy_number.tsv` is a shared
  species-by-orthogroup copy-number matrix used by CAFE and by trait PGLS,
- `run_orthogroup_copy_number_trait_pgls=1` tests associations between
  orthogroup copy numbers and `workspace/input/species_trait/species_trait.tsv`
  with species-tree PGLS; it does not use CAFE5 result files such as
  `Gamma_change.tab`,
- orthogroup copy-number trait-PGLS outputs are written under
  `workspace/output/genome_evolution/orthogroup_copy_number/trait_pgls/`,
- `grampa_h1` and `target_branch_go` default to empty strings; leaving them empty skips GRAMPA or GO enrichment only,
- GO target can be specified by species name or branch ID.

### `gg_gene_summary_entrypoint.sh`

Purpose:

- run mode-aware cross-family summaries for `query2family` or `orthogroup` outputs,
- generate species-tree presence/absence and copy-number matrices for `query2family` or `orthogroup`,
- optionally build the selected source's `gg_orthogroup.db`,
- optionally run HGT and convergent-site summaries against the selected source.

Main Gene-Family Source:

- `gene_family_source="query2family"` or `gene_family_source="orthogroup"`

Main presence/absence outputs:

- `workspace/output/gene_summary/query2family/query2family_presence_absence.tsv`
- `workspace/output/gene_summary/query2family/query2family_copy_number.tsv`
- `workspace/output/gene_summary/query2family/query2family_presence_absence.pdf`
- `workspace/output/gene_summary/*/*_presence_absence.plot_selection.tsv`
- `workspace/output/gene_summary/*/*_presence_absence.plot.long.tsv`

Notable defaults:

- `run_presence_absence_summary=1`
- query-gene names in the presence/absence plot come from
  `workspace/input/query_gene` file basenames
- `presence_absence_max_families=auto` plots all query2family queries but only the first 100 orthogroups by default
- `presence_absence_family_ids` and `presence_absence_family_file` select an explicit plotted subset for either mode
- `presence_absence_species_tree=auto` prefers query2family-pruned dated species
  trees when available
- `presence_absence_species_tree_ci=auto` adds MCMCtree 95% HPD node-age bars for
  dated species trees when `mcmctree_95CI.nwk` is available
- `presence_absence_species_tree_support=auto` transfers numeric branch-support
  labels from matching species-tree support files
- `presence_absence_busco_table=auto` adds per-species BUSCO stacked bars when
  BUSCO full tables or `workspace/output/species_tree/busco_summary_table/busco_summary.tsv`
  are available; full tables include Fragmented counts
- `presence_absence_plot_width=7.2`; the plotter caps figure width at 7.2 inches
- `run_gene_family_database_build=0`, `run_csubst_scan_aa_change_summary=0`, `run_csubst_scan_candidate_sites=0`, `run_hgt_candidate_summary=0`, `run_hgt_summary_plots=0`, and `run_csubst_site_convergence_summary=0`
- database, CSUBST scan AA-change summary, HGT, and CSUBST site convergence flags are valid for both sources and use the selected source's gene-family output directory,
- `run_gene_family_database_build=1` assembles the selected source's SQLite DB from
  logical `stat_tree` and `stat_branch` inputs, whether they are live files or
  members of GeneGalleon ZIP shards; structural columns are required, while
  analysis-specific columns that are absent from individual families are stored
  as SQL `NULL`; before the DB is replaced, content-based provenance manifests
  and CSUBST/`stat_branch` descendant-tip branch identities are audited across
  the selected source, with results written to
  `gene_summary/<source>/<source>_artifact_provenance_audit.tsv`; legacy outputs
  without manifests are reported as `legacy_untracked` and remain usable,
- overwrite builds use a temporary SQLite file and replace the published database
  only after every input has been loaded and indexed successfully,
- when present, `csubst_scan/` is imported as DB table `aa_change`, and
  `csubst_scan_units/` is imported as `aa_change_unit`; `aa_change` receives
  global BH-FDR columns after all candidate substitutions are loaded,
- `run_csubst_scan_aa_change_summary=1` writes
  `*_csubst_aa_change_min_support_2_summary.tsv` ranked by the global FDR
  columns and, when `aa_change` candidates are available,
  CSUBST scan plots for foreground-support-binned significance rates and the
  substitution spectrum, plus a four-panel comparison of analytical,
  candidate-level empirical, and full-scan empirical maxT P values with their
  corresponding global BH-FDR q values; all primary filenames include
  `min_support_2`, matching the upstream default; the same directory also
  receives the integer support-threshold series from 3 to the observed maximum,
  with filtered TSVs, recalculated BH q values, P/q-distribution PDFs, and a
  manifest; orthogroup TSVs include the five representative `besthit_*`
  annotation columns from `Orthogroups.GeneCount.annotated.tsv` when available,
  and those columns are retained throughout the threshold series; use it with
  `run_gene_family_database_build=1` to refresh the DB and plots in one run,
- `run_csubst_scan_candidate_sites=1` selects candidates independently from each
  recalculated min-support summary (default: analytical global
  `q_rate_enrichment_global <= 0.05`, thresholds from the observed maximum down
  to 5), analyzes each unique candidate once, and writes a self-contained ZIP
  per threshold directly under `gene_summary/<source>`; every candidate entry
  contains its annotated row, raw `csubst sites` outputs, an exact-site tree
  PDF, and a combined report; PDB search is disabled by default,
- `run_csubst_site_convergence_summary=1` runs site-level convergence screening with
  `csubst_site_wrapper.py`, combines orthogroup results with species traits,
  and writes convergence outputs under the selected source's `csubst_site`
  summary directory.

### `gg_progress_summary_entrypoint.sh`

Purpose:

- summarize stage outputs into TSV reports.

Main outputs in the workspace root (`workspace/`):

- `workspace/orthogroup_summary.tsv`
- `workspace/query2family_summary.tsv`
- `workspace/transcriptome_assembly_summary.tsv`

Note:

- this stage runs `workflow/core/gg_progress_summary_core.sh` inside the container.
- orthogroup summary generation is skipped when the selected gene-count table is absent; AMAS inputs are optional and may be live or ZIP-backed.
- query2family summary generation is skipped when `workspace/output/query2family`
  or `workspace/input/query_gene` is absent.

### Utility wrappers

- Gene-family ZIP manager:
  - `workflow/gg_gene_family_archive.sh`
  - lists, verifies, deletes, undeletes, restores, and explicitly archives logical query2family/orthogroup artifacts
- Species-tree stage ZIP manager:
  - `workflow/gg_species_tree_archive.sh`
  - converts, verifies, and reports the high-file-count `species_tree/single_copy_*` stage directories
- Versions dump:
  - collector script: `workflow/support/gg_versions.sh`
  - auto-triggered at the end of each `gg_*_entrypoint.sh` on successful completion
  - outputs are written to `workspace/output/versions/*.log`
- Dependency-aware debug runner:
  - `workflow/gg_all_entrypoints_debug.sh`
  - can dry-run or benchmark the main entrypoint chain
  - writes `workspace/output/debug_entrypoint_logs/summary.tsv`
- Minimal dataset builder:
  - `workflow/support/build_minimal_test_dataset.py`
  - extracts a smaller development dataset from an existing workspace
