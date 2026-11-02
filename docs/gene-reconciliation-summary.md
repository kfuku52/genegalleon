# Reconciliation database and downstream evidence

`gg_gene_summary` can import existing GeneRax recPhyloXML into the selected
gene-family database, then analyse that database without changing it. It does
not rerun GeneRax, CAFE, RNA quantification, expression fitting or MS searches.

Enable `run_gene_family_database_build=1` and
`gene_family_database_reconciliation=1`. Set `reconciliation_species_tree` to
the complete rooted Newick tree used for reconciliation. XML is read through
the normal raw/ZIP store under `generax_xml/<family>_generax.xml`; it is never
expanded into thousands of temporary files. Existing provenance audits,
selected-catalogue scope and publication fences still apply. Rebuild rather
than append when importing reconciliation.

The optional `reconciliation_method_manifest` JSON preserves method metadata:

```json
{
  "gene_model_contract_sha256": "<SHA-256 of the frozen model-set contract>",
  "default": {"program": "GeneRax", "MAD": true},
  "families": {
    "OG0000001": {
      "core_commit": "<actual commit>",
      "binary_sha256": "<actual binary SHA>",
      "adoption_decision_sha256": "<actual decision SHA>",
      "upstream_convention_validated": false
    }
  }
}
```

Use recorded parameters, commits, qualification and adoption references from
the actual run. Missing metadata stays unspecified; the importer does not
invent a default model or treat a project adoption as upstream validation.
The model-set contract must bind the exact gene-model/protein inputs shared
by reconciliation, expression and proteomics, not just a project name.

The added tables are `reconciliation_meta`, `reconciliation_species_branch`,
`reconciliation_family`, `reconciliation_clade`, `reconciliation_event` and
`reconciliation_gene`. They preserve XML SHA, event order (including multiple
events on one gene-tree clade), parent relationships, exact extant gene IDs,
species and family membership. They are separate from `tree` and `branch`
statistics produced by NWKIT's LCA reconciliation. Missing XML is explicitly
`missing_xml`, rather than zero loss events. Invalid XML, unknown events,
duplicate genes or inconsistent species trees fail the private build.
When present, each family's original `artifact_provenance/<family>.generax.json`
is preserved with its SHA, parameters, input/output identities and diagnostics.
Its declared XML SHA must match. A missing legacy producer manifest stays SQL
`NULL`; separately supplied method/adoption metadata does not manufacture one.

Branch identity uses descendant species sets. Internal node renumbering and
sibling order do not change the mapping. The complete rooted topology and
species set must agree; no pruning or numeric-label guessing is performed.
XML `loss` locations describe extinction of reconstructed lineages. They do
not identify a molecular deletion, an extant lost-gene ID or its exact timing.
Extant descendants are reported as observations, not inferred ancestral copy
numbers or proof of uninterrupted retention.

## Per-gene product, GO and localisation database

Set `run_gene_annotation_database=1` and `gene_annotation_manifest` to a
SHA-bound manifest of existing `gg_genome_annotation` TSVs. This stage runs after
the model-bound reconciliation DB is available and before its downstream summaries.
Enable `gene_family_database_reconciliation=1` when building that base DB. It copies every
existing table into `gg_orthogroup.annotated.db`, adds the gene annotations, and
uses that annotated DB for subsequent stages in the same invocation. The base
`gg_orthogroup.db` remains unchanged. GeneRax, annotation searches and predictors
are not repeated. This also supports a completed base DB from an earlier run.

```json
{
  "schema": "gene_annotation_database_v1",
  "gene_model_contract_sha256": "<same model contract as the base DB>",
  "family_catalog_sha256": "<same family catalogue as the base DB>",
  "source_database_sha256": "<optional exact base DB SHA>",
  "require_complete_genes": true,
  "sources": [
    {
      "species": "Species_name",
      "path": "annotations/Species_name_annotation.tsv",
      "sha256": "<exact TSV SHA>",
      "product_column": "sprot_recname"
    }
  ]
}
```

Paths resolve against the manifest. By default the sources must cover the exact
species set and all reconciled genes. Gene IDs are joined literally with species;
no prefix stripping, guessed aliases or family-name substitution is performed.
Duplicate genes or inconsistent family assignments fail before publication.
Explicit `require_complete_genes=false` permits partial annotations and records
the number of unannotated reconciliation genes. It does not fill missing values.

- `gene_annotation` retains each native TSV column as queryable text, including
  Swiss-Prot hits, GO strings, targeting signals, localisation scores, model SHA,
  score kind, calibration and decision state. Columns absent from a species's
  input are SQL NULL; a present but empty source cell stays empty. The key is
  `(species,gene_id)`. The `annotation_source_id` joins its source metadata.
- `orthogroup` contains only membership in the DB's selected catalogue;
  `source_orthogroup` also preserves original assignments outside that catalogue.
  Annotation-only genes remain available, without inventing a GeneRax leaf.
- `gene_annotation_product` provides a consistent `product` view, using each
  source's declared product column (default `sprot_recname`).
- `gene_annotation_go` preserves ordered GO ID, aspect, name and evidence tuples.
  Default columns are `go_ids`, `go_aspects`, `go_terms`, `go_evidence`, separated
  by semicolons. An optional `go_columns` object maps `ids`, `aspects`, `terms`,
  `evidence` to other source columns. Inconsistent tuple lengths fail rather than
  silently pair unrelated names or evidence codes.
- `gene_annotation_source` and `gene_annotation_meta` record exact source paths,
  hashes, original column schemas, model/family contracts and coverage audit.

Prediction scores retain their original semantics. For example, an uncalibrated
CDSKIT `sigmoid_score` is not relabelled a calibrated probability. Product names
describe extant genes, not identified extinct copies. The reconciliation summary
uses embedded product annotations for its exact reconciled members; an external
product TSV that contradicts them is rejected. GO/localisation queries can use
the same annotated DB directly, without reopening the original TSVs.

The public helper is also available for an admitted standalone workflow:

```bash
python workflow/support/gene_annotation_database.py \
  --source-db workspace/output/orthogroup/gg_orthogroup.db \
  --manifest gene_annotations.json \
  --output-db workspace/output/orthogroup/gg_orthogroup.annotated.db \
  --scratch /path/to/private/scratch
```

Import uses bounded batches and one database output, with no per-gene files.
Private scratch avoids random database writes on shared storage. The published
base DB is held under its normal shared build guard, and the derivative uses its
own exclusive guard; unrelated databases remain available. Changed inputs,
active journals or invalid annotations preserve an existing annotated output.

## Optional downstream stage

The optional downstream stage remains disabled in debug mode unless its
configuration is supplied; debug mode does not invent an evidence cohort.

Set `run_reconciliation_summary=1` and `reconciliation_summary_config` to a
JSON file. Paths below resolve relative to that file. Every external input
requires its exact SHA-256. Example structure:

```json
{
  "schema": "reconciliation_summary_v1",
  "family_catalog_sha256": "<reconciliation_meta.family_catalog_sha256>",
  "gene_model_contract_sha256": "<same contract as the database>",
  "require_complete_reconciliation": true,
  "expression_p_value_semantics": "one_sided_trap_enrichment_positive_logfc_v1",
  "expression_fdr": 0.05,
  "expression_min_logfc": 1.0,
  "proteome_species": "Species_name",
  "inputs": {
    "branches": {"path": "branches.json", "sha256": "<SHA>"},
    "cafe": {"path": "cafe_comparison.tsv", "sha256": "<SHA>"},
    "expression_transcript": {"path": "transcript_fits.tsv", "sha256": "<SHA>"},
    "expression_cds": {"path": "cds_fits.tsv", "sha256": "<SHA>"},
    "proteome_members": {"path": "searched_species_genes.tsv", "sha256": "<SHA>"},
    "proteome": {"path": "unique_protein_support.tsv", "sha256": "<SHA>"},
    "shared_peptides": {"path": "shared_family_psms.tsv", "sha256": "<SHA>"},
    "gene_annotations": {"path": "gene_products.tsv", "sha256": "<SHA>"}
  }
}
```

Only `branches` is mandatory; omit unused evidence inputs and associated
parameters. Functional evidence requires matching nonempty model-set contracts.
`require_complete_reconciliation=true` prevents downstream publication if any
selected family lacks XML. Its default is false, with missing values retained.

- Branch JSON is a list of objects with unique `branch_id` strings and exact
  `descendants` species arrays. Each target must be a clade in the DB tree.
  Taxa, branch scope and labels belong in project configuration, not the code.
- CAFE TSV columns: `branch_id`, `orthogroup`, integer `copy_change`, and
  `significant_decrease` (`0` or `1`). Include every target branch/family pair,
  with significance derived from the original declared CAFE test/cutoff. An
  optional `status` column defaults to `analysed`. For families excluded from
  the original CAFE analysis, explicitly use `not_analysed` and leave both
  numeric fields blank. Every branch/family pair is still required; an omitted
  row is an error. Unknown results remain SQL NULL, and the output reports
  analysed, unanalysed and mutually comparable family counts for each branch.
  Unanalysed families cannot become `generax_extinction_only`, and they remain
  eligible for GeneRax-based expression analysis. A
  significant decrease requires negative copy change. Missing analysis is not
  silently converted to a nonsignificant result. Net CAFE copy change and
  GeneRax lineage extinction measure different quantities; overlap is reported
  descriptively, not as validation of one model by the other.
- Expression TSV columns: `orthogroup`, `status`, `logFC`, `p_value`. Include
  the full family catalogue, including untested rows and their original status.
  `status=success` requires finite logFC and a directional P value in [0,1].
  Existing fits must use the same models, contrasts and biological cohort.
  For each target branch, BH is recomputed among successfully tested families
  with a GeneRax extinction event at that branch. No correction across branches
  is added. Enrichment requires `q < expression_fdr` and
  `logFC > expression_min_logfc`. Raw fits are reused; quantification and fitting
  are not repeated. Evidence is family-level, not expression of an extinct copy.
- Proteome membership TSV: `gene_id`, `orthogroup`, covering *every* gene of
  `proteome_species` in the DB. Unique protein support additionally requires
  `unique_unmodified_peptides_IL_equivalent`, `best_peptide_q` and
  `best_protein_q`; rows must have at least one unique peptide and both q values
  at most 0.01. This is an admitted search-result summary, not a raw PSM filter.
  Identified proteins explicitly unassigned to a selected family (blank
  `orthogroup` and no DB membership) are excluded from family joins and counted
  in the evidence audit; contradictory or unknown nonempty assignments fail.
- Shared family PSM TSV: `orthogroup`, `peptide`, `MS_file`, `spectrum`,
  `proteins`, `peptide_q`, `protein_q`. `proteins` lists two or more unique gene
  IDs separated by semicolons, all assigned to that family in the searched
  species. Each file/spectrum is unique and both q values must be at most 0.01.
  Shared peptides support a family, never an individual gene.
  A family without a searched-species gene is `not_assessable_no_gene`.
  `not_detected` is not proof of absence or lack of secretion.
- Optional annotation TSV: unique `gene_id` and `product`, using DB gene IDs.
  Products describe extant family members, not identified extinct genes.

The stage opens the source DB read-only and publishes four files under
`gene_summary/<source>/reconciliation/`: `reconciliation_summary.db`,
`branch_family.tsv`, `branch_summary.tsv` and `manifest.json`. File count is
independent of family or branch count. No per-family result directories are
created. Previous results survive validation failure, and changed source inputs
abort publication. The manifest records source/evidence hashes and interpretation
limits; normal stage provenance also binds all nested configuration inputs.

For a separate admitted workflow, the same helper is available as:

```bash
python workflow/support/reconciliation_summary.py \
  --dbpath workspace/output/orthogroup/gg_orthogroup.db \
  --config project_reconciliation_summary.json \
  --output-dir workspace/output/gene_summary/orthogroup/reconciliation
```
