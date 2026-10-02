# Gene-family presence/absence plots

Run the `gg_gene_summary` stage to reuse existing reconciled family trees, CDS
FASTA files and optional local-synteny tables. No sequence search or phylogeny
reconstruction is required. See [stage outputs and defaults](main-stages-and-what-they-do.md#gg_gene_summary_entrypointsh).

The default remains `presence_absence_ortholog_basis=reference_species`, with
all original queries retained when the query-gene basis is selected. Species
and gene identifiers in scientific tables are never changed by display labels.

## Prefer nearer query-source species

Set these parameters in `workflow/gg_gene_summary_entrypoint.sh`, or use their
`GG_GENE_SUMMARY_` scoped environment overrides:

```bash
presence_absence_ortholog_basis="query_gene"
presence_absence_query_selection="closest"
presence_absence_target_species="Target_species"
presence_absence_query_metadata="/path/to/query_metadata.tsv"
```

Source identity must be explicit. A query FASTA header can supply it:

```fasta
>query_id | Human-readable name | species=Source species
MERNLLS
```

For gene lists, or source species absent from the species tree, provide a TSV:

```tsv
family_id	query_id	source_species	tree_species
family_A	query_id	Source_relative	Source_species
```

`tree_species` is optional: it explicitly maps an unsequenced source relative
to a species-tree tip for ranking. The original `source_species` is preserved.
Conflicting source metadata and unknown query IDs within a selected family are
errors. Unknown/unmapped sources are retained; gene-ID prefixes and best-hit
anchor species are **not** used to infer query origin.

Nearness is ranked by shared ancestry: `distance` is the number of edges from
the target tip back to its MRCA with the source in the supplied species tree,
not total source/target path length, gene-tree distance or an inferred divergence
time. Extra sampling below the source clade does not make it appear more distant.
Only a strictly nearer retained query can replace a farther query. Both must
map to the same tip or have a reconciled speciation MRCA, and the replacement
must cover every exact ortholog gene represented by the farther query.
Equal-distance queries and distinct paralog lineages remain. By default this
coverage check includes all species in the family; optionally restrict it with
`presence_absence_selection_species="Species_one,Species_two"`. The union of
gene IDs in those species is checked before and after selection.

The `.selection.tsv` records every original query, source, distance, decision
and replacement reason. `.query_map.tsv` retains the selected records and their
anchors, with additive `source_species` and `hog_ids` columns. The original
query files, family trees and saved artifacts are never rewritten.

## Combine saved query families and OG/HOG anchors

Set `presence_absence_family_manifest` to an ordered TSV, and select
`presence_absence_ortholog_basis="query_gene"`. Its rows define the complete
query-anchor plot, independently of the ordinary family-summary subset.
`presence_absence_species_tree` can specify the plotting/ranking tree explicitly.

Required manifest columns are `family_id`, `source_dir`, `source_family_id`.
Each row must have **either** `query_file` **or** `anchor_species`:

```tsv
family_id	source_dir	source_family_id	query_file	anchor_species	hog_table	hog_ids
enzyme_queries	../output/query2family	enzyme	../input/query_gene/enzyme
contracted_groups	../output/orthogroup	OG0000123		Species_one	N0.tsv	HOG0000100;HOG0000101
```

Paths are relative to the manifest's directory, unless absolute. Saved stores
can be live files or GeneGalleon's existing ZIP-backed artifacts. Plot block
IDs (`family_id`) must be unique; `source_family_id` is the original artifact
ID. A missing tree, CDS FASTA, HOG or membership join is an error, not absence.

`anchor_species` selects **all** saved tree tips of that species. Optional
`hog_table` and semicolon-separated `hog_ids` restrict anchors to those HOGs.
The table is an OrthoFinder hierarchical-orthogroup TSV with `HOG`, `OG`,
optional `Gene Tree Parent Clade`, and species columns containing comma-separated
gene IDs. Exact CDS IDs are preferred; unprefixed gene IDs are accepted only
when their species column resolves to exactly one saved tree tip. Every selected
HOG must belong to `source_family_id`, and its members must exist in both tree
and CDS FASTA. Set `presence_absence_query_label="label"` to show HOG membership
in generated anchor labels.

HOG membership is annotation, not an orthology or monophyly assumption. The
**full original OG tree** defines S/D orthology and topology support, including
when two HOGs are not reciprocally monophyletic. Prefer one combined manifest
row for multiple HOGs in the same OG. `.overlap.tsv` identifies unique genes
represented in multiple plot blocks, and the producer warns against summing
those blocks as independent counts. The query-anchor `.long.tsv` describes full
source-family tip counts, not HOG-restricted or selected-query copy counts.

## Display settings

```bash
presence_absence_plot_width="auto"  # Default remains 7.2 inches; explicit wider values are honored.
presence_absence_query_label="label"  # Default: id.
presence_absence_focus_species="Species_one,Species_two"
presence_absence_legend_columns="auto"  # auto, 1, 2, or 3.
presence_absence_label_map="/path/to/display_labels.tsv"
```

Auto width reserves at least 0.14 inches per column plus tree/label space on
one page. Rotated labels use font-metric bounds with safety padding. Explicit
widths are not silently capped. The ortholog legend adapts to available width;
focus rows are outlined without changing counts or evidence.

The display map has `kind`, `id`, `label` columns. `kind=family` uses the plot
block/family ID; `kind=column` uses the exact anchor CDS FASTA ID:

```tsv
kind	id	label
family	contracted_groups	Contracted candidate family
column	Species_one_gene_123	HOG0000100: gene 123
```

The plotting helper also accepts these as `--width=auto`, `--label_map=PATH`,
`--focus_species=Species_one,Species_two`, and `--legend_columns=auto`. Empty
`--out_pdf=` / `--out_svg=` explicitly disable that format; omitted output
options also leave that format disabled.

## Interpret evidence conservatively

Undetected orthologs are not proof of genomic deletion. Copy-number contraction
does not require complete absence, and multiple anchor relationships are not
independent loss events. Local-synteny support requires the existing minimum
of two distinct shared neighbor groups; one shared group is separate evidence.
Unavailable neighborhoods/topology support remain unavailable rather than zero.
Self-anchor comparisons are separate from evaluated pair evidence. None of the
selection or display settings changes these scientific thresholds.
