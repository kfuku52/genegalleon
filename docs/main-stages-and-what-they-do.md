# Main Stages and What They Do

### `gg_input_generation_entrypoint.sh`

Purpose:

- optional pre-stage to automate input preparation,
- optionally download provider files from a manifest and format into
  `workspace/output/input_generation/species_cds` and
  `workspace/output/input_generation/species_gff` and
  `workspace/output/input_generation/species_genome`,
- optionally generate `workspace/input/species_trait/species_trait.tsv` using
  target species from `download_plan.xlsx` (default) or `species_cds`.

Main scripts:

- `workflow/core/gg_input_generation_core.sh`
- `workflow/support/format_species_inputs.py`
- `workflow/support/generate_species_trait.py`

Notable defaults:

- `run_validate_inputs=1` validates CDS naming rules, CDS/GFF species-set consistency, and CDS-to-GFF mapping compatibility, with species-level mapping checks parallelized via `validate_cds_gff_mapping.py --nthreads`,
- formatted outputs default to `workspace/output/input_generation/species_cds`, `workspace/output/input_generation/species_gff`, and `workspace/output/input_generation/species_genome`,
- per-run metadata defaults to `workspace/output/input_generation/gg_input_generation_runs.tsv`,
  `workspace/output/input_generation/gg_input_generation_species.tsv`, and
  `workspace/output/input_generation/download_plan.resolved.tsv`,
- trait generation is opt-in (`run_generate_species_trait=0` by default),
  and writes `workspace/input/species_trait/species_trait.tsv` when enabled,
- `trait_profile=gift_starter` is available as a quick preset:
  sets `run_generate_species_trait=1` and uses `trait_databases=gift`
  when `trait_databases` is unset or `auto`,
- each run appends a TSV record to `workspace/output/input_generation/gg_input_generation_runs.tsv`
  (configurable via `summary_output` or `GG_INPUT_SUMMARY_OUTPUT`).

Wrapper-specific notes:

- the wrapper accepts dedicated host-side overrides such as `GG_INPUT_DOWNLOAD_MANIFEST`,
  `GG_INPUT_PROVIDER`, `GG_INPUT_TRAIT_PROFILE`, and `GG_INPUT_SPECIES_*`,
- `provider=refseq` and `provider=genbank` are accepted by the wrapper as aliases of `ncbi`.

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

- `assembly_method="rnaSPAdes"`
- `kallisto_reference="longest_cds"`
- staged FASTA outputs are written as `.fa.gz`.

### `gg_genome_annotation_entrypoint.sh`

Purpose:

- per-species CDS/genome annotation and QC,
- BUSCO (CDS/genome),
- UniProt annotation (`blastp` or `mmseqs2`),
- optional MMseqs2 taxonomy and contamination removal,
- optional genome analyses (SubPhaser, dotplot, GenomeScope).

Main outputs:

- `workspace/output/species_cds_annotation`
- `workspace/output/species_cds_busco_full`, `species_cds_busco_short`
- `workspace/output/species_genome_busco_full`, `species_genome_busco_short`

Notable defaults:

- most heavy tasks default to `0`,
- `uniprot_annotation_method="mmseqs2"` (set `blastp` to use NCBI BLASTP for UniProt annotation),
- `run_multispecies_summary=1` by default.

### `gg_genome_evolution_entrypoint.sh`

Purpose:

- unified genome-evolution entrypoint that serially runs:
  - inlined species-tree stage (formerly `gg_speciesTree_core.sh`)
  - inlined orthofinder stage (formerly `gg_orthofinder_core.sh`)
  - inlined genome-evolution stage (formerly `gg_genomeEvolution_core.sh`)
- preserves output-exists skip behavior at step level and aborts on real failures.
- accepts either CDS-first input (`input_sequence_mode=cds`) or protein-only input (`input_sequence_mode=protein`).

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
- compare orthogroup methods.

Main outputs:

- `workspace/output/orthofinder`

Temporary protein FASTA files are created under `workspace/downloads/tmp/` and removed automatically after orthogroup-related steps finish.
When `workspace/input/species_genetic_code/species_genetic_code.tsv` is present, it overrides the global `genetic_code` on a per-species basis during CDS-to-protein translation; species missing from the table still use the default `genetic_code`.

Notable defaults:

- `orthogroup_table="HOG"`
- `orthogroup_annotation_method="mmseqs2"` (set `blastp` to use NCBI BLASTP for representative-gene UniProt annotation)
- species-tree-aware OrthoFinder is used when species tree exists.

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
- the same multi-species BUSCO run also feeds the BUSCO-based genome-evolution branch; BUSCO-based genome-evolution steps reuse the BUSCO outputs driven by `run_species_busco` and `run_species_get_busco_summary`,
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

Notable defaults:

- `mode_gene_evolution="query2family"`
- `run_rps_blast=1` (Pfam domain annotation is on by default),
- `run_tree_plot=1`
- `run_summary=1`
- `uniprot_annotation_method="mmseqs2"` (set `blastp` for NCBI BLASTP-based UniProt annotation),
- many advanced analyses default to `0`.

Current behavior notes:

- if `run_generax=1`, initial IQ-TREE disables UFBOOT and support is
  recomputed after GeneRax on the GeneRax topology,
- Pfam RPS-BLAST DB (`Pfam_LE`) is auto-prepared when missing, with lock-based
  synchronization for array jobs,
- gene-tree/species-tree PGLS outputs are `gene_tree_PGLS.tsv` and
  `species_tree_PGLS.tsv`.

Wrapper-specific note:

- `gg_gene_evolution_entrypoint.sh` forwards all top-block variables to the container runtime,
- routine changes should be made by editing the top config block,
- unlike `gg_input_generation_entrypoint.sh`, this wrapper does not expose a separate host-side `GG_*` override map for all top-block variables.

### Inlined Stage: Genome Evolution

Purpose:

- BUSCO-based and orthogroup-based GRAMPA workflows,
- optional CAFE and GO enrichment analyses.

Main outputs:

- `workspace/output/genome_evolution/*`

Notable defaults:

- duplicate-aware BUSCO genome-evolution substeps default to `0`,
- `run_orthogroup_grampa=1`, but GRAMPA is auto-disabled unless rooted
  orthogroup trees exist and `grampa_h1` is set,
- `run_cafe=0`, `run_go_enrichment=0` by default,
- `grampa_h1` and `target_branch_go` default to empty strings; leaving them empty skips GRAMPA or GO enrichment only,
- GO target can be specified by species name or branch ID.

### `gg_gene_database_entrypoint.sh`

Purpose:

- assemble SQLite DB from orthogroup summary tables.

Main output:

- `workspace/output/orthogroup/gg_orthogroup.db`

Required input directories:

- `workspace/output/orthogroup/stat_tree`
- `workspace/output/orthogroup/stat_branch`

Notable defaults:

- `run_database_prep=1`
- the stage skips database generation if the required `stat_tree/` or `stat_branch/` directories are missing.

### `gg_progress_summary_entrypoint.sh`

Purpose:

- summarize stage outputs into TSV reports.

Main outputs in the workspace root (`workspace/`):

- `workspace/orthogroup_summary.tsv`
- `workspace/transcriptome_assembly_summary.tsv`

Note:

- this stage runs `workflow/core/gg_progress_summary_core.sh` inside the container.
- orthogroup summary generation is skipped when the selected gene-count table or AMAS directories are absent.

### `gg_convergent_sites_entrypoint.sh`

Purpose:

- run site-level convergence screening with `csubst_site_wrapper.py`,
- combine orthogroup results with species traits,
- generate convergence summary tables and branch-level reports.

Main output root:

- `workspace/output/csubst_site`

Notable defaults:

- `arity_range="2-10"`
- `trait="all"`
- `skip_lower_order="yes"`
- `file_trait`, `dir_orthogroup`, `dir_orthofinder`, and `dir_out` auto-resolve to workspace defaults unless explicitly overridden.

### Utility wrappers

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
