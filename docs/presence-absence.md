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

## Interpret duplication bars

Colored bars count distinct original species-overlap **D nodes with displayed
genes in both child subtrees**, grouped by family and mapped species-tree branch.
Displayed genes include strict orthologs and additional ortholog candidates in
the matrix. Each gene ID contributes once to subtree membership even when it
appears in several query columns or glyphs, and each qualifying D node contributes
one event. D nodes with displayed genes in only one child are excluded.

The original full-tree S/D calls remain unchanged; the method does not rerun
species overlap on a pruned tree or turn candidate-associated D nodes into S.
Internal-node mapping still uses the saved GeneRax assignment, or species
coverage when that assignment is unavailable. Query/anchor/HOG selection and
the candidate cutoff can change which genes are displayed and therefore which
D nodes qualify. Standalone plots also restrict membership to the species rows
actually displayed. The legend says **Bar height = displayed-gene duplication
count**. Its numbered bars are height references; an upper legend value need
not occur in the data.

The `.tree.tsv` preserves all original D rows and adds
`displayed_child1_gene_ids` and `displayed_child2_gene_ids` for auditability.
`displayed_gene_ids` records the complete displayed gene set once per family,
on the compact anchor-tree root. The plot checks that this set agrees with its
glyph table. Legacy tree tables without these three columns remain readable
with a warning and the **full-family duplication count** legend; regenerate the
summary tables to apply the displayed-gene scope.

Counts are inferred events, not present-day gene copy numbers or confirmed
biological duplications. Gene-tree uncertainty and the underlying gene models
still matter; no bootstrap filter is applied to these counts.

## Flag candidates across low-confidence duplication nodes

Set the threshold in the editable configuration block of
`workflow/gg_gene_summary_entrypoint.sh`:

```bash
presence_absence_dup_conf_score_threshold="0.05"
```

For a one-off run showing every original query, use scoped overrides from the
repository root after configuring the workspace/input/output paths:

```bash
GG_GENE_SUMMARY_RUN_PRESENCE_ABSENCE_SUMMARY=1 \
GG_GENE_SUMMARY_PRESENCE_ABSENCE_ORTHOLOG_BASIS=query_gene \
GG_GENE_SUMMARY_PRESENCE_ABSENCE_QUERY_SELECTION=all \
GG_GENE_SUMMARY_PRESENCE_ABSENCE_PLOT_WIDTH=auto \
GG_GENE_SUMMARY_PRESENCE_ABSENCE_DUP_CONF_SCORE_THRESHOLD=0.05 \
bash workflow/gg_gene_summary_entrypoint.sh
```

Use `0.1` or `0.2` in place of `0.05` to compare cutoffs. Set a distinct
`summary_output_dir` (or `GG_GENE_SUMMARY_SUMMARY_OUTPUT_DIR`) for each run to
retain all versions; the same output directory otherwise replaces the summary
files. No rerun of sequence search, tree inference or reconciliation is needed.
The query-gene PDF/SVG and `.dup_conf.tsv` files use the
`query2family_query_gene_orthologs` prefix in that directory.

The default `0` preserves strict orthology calls. A positive threshold adds
**orange additional ortholog candidate glyphs** when a gene and an anchor from
different species have a D MRCA whose duplication confidence score is less
than or equal to the threshold. Existing blue orthologs and the original D
events remain intact. Both reference-species and query-gene bases support this
setting, including saved OG/HOG sources. Query selection still uses strict
orthology; candidates cannot authorize replacing a query or merging anchors.

The legend labels these as **Additional ortholog candidate**, with
`duplication confidence score <= threshold` on a separate line. The saved
`weak_duplication` relation identifier is retained for compatibility.

The score is the Jaccard overlap of species under the two children, recomputed
on the full saved gene tree. It is not a bootstrap value or a probability that
the duplication is real. Same-species paralogs are never added by this setting.
Gapped anchor sets are split into separate glyphs so intervening columns are
not painted. A copy spanning several anchors is counted once in its glyph;
column totals are not independent copies. Family-level presence and sequence
counts remain unchanged.

The new `.dup_conf.tsv` companion records every added gene/anchor pair, its
original D MRCA, shared/union species counts, score, threshold and raw branch
UFBoot when available. Orthology UFBoot remains unavailable for these D pairs,
with reason `weak_duplication`; branch support does not establish orthology.
Local synteny is evaluated separately. The collector accepts
`--dup_conf_score_threshold 0.05 --out_dup_conf candidates.dup_conf.tsv`, and
requires `--out_dup_conf` whenever the cutoff is positive. This prevents saving
additional candidate glyphs without their audit table. Similarly,
standalone candidate plots require both `--dup_conf_score_threshold=0.05` and
`--ortholog_dup_conf_table=candidates.dup_conf.tsv`. The summary stage forwards
these automatically. The plot checks the pair identities, scores and cutoff
against the glyphs, and checks supplied synteny/UFBoot identities and candidate
classifications. Mixed results are errors. Added candidates must have a saved
CDS sequence.
Every additional candidate glyph must have one D MRCA and a positive species
overlap. When a tree table is supplied, the plot checks that MRCA against its
original D nodes, and with displayed-descendant provenance also checks that the
gene/anchor pair lies in opposite child subtrees, even without a UFBoot table.
Ortholog plots reserve enough vertical space for copy-number text, including
strict-only and stacked candidate lanes. Automatic height can therefore
increase, especially with evidence bands.
An explicit height below the required minimum is rejected with that minimum.

## Interpret evidence conservatively

Undetected orthologs are not proof of genomic deletion. Copy-number contraction
does not require complete absence, and multiple anchor relationships are not
independent loss events. Local-synteny support requires the existing minimum
of two distinct shared neighbor groups; one shared group is separate evidence.
Unavailable neighborhoods/topology support remain unavailable rather than zero.
Self-anchor comparisons are separate from evaluated pair evidence. None of the
selection or display settings changes these scientific thresholds.
