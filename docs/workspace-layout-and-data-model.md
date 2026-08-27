# Workspace Layout and Data Model

GeneGalleon uses a split workspace layout and treats the workspace as the persistent project state.

## Default workspace root

By default, the host-side workspace is:

- `workflow/../workspace`

At runtime, that directory is bind-mounted into the container as:

- `/workspace`

You can relocate the workspace for a run by exporting:

```bash
gg_workspace_dir=/path/to/workspace bash workflow/gg_genome_evolution_entrypoint.sh
```

## Split layout

```text
workspace/
  input/
    species_cds/
    species_protein/
    species_genetic_code/
    species_expression/
    species_gff/
    species_genome/
    species_dnaseq/
    species_trait/
    species_rnaseq/
    input_generation/
    amalgkit_metadata/
    query_sra_id/
    query_gene/
  output/
  downloads/
```

The layout is always split. Legacy mixed-layout assumptions are no longer part of the current runtime helpers.

## Directory roles

### `workspace/input/`

This is the user-managed side of the workspace.

Typical contents:

- `species_cds/`: CDS FASTA inputs for comparative/genome-evolution workflows;
  sequence IDs must follow `GENUS_SPECIES_GENEID` as described in
  [Input Conventions](input-conventions.md)
- `species_protein/`: optional per-species protein FASTA inputs for `gg_genome_evolution` protein mode
- `species_genetic_code/`: optional `species_genetic_code.tsv` table for per-species CDS translation overrides
- `species_gff/` and `species_genome/`: optional annotation/genome companions
- `species_expression/`: expression matrices used by gene-family downstream analyses
- `species_trait/`: trait tables such as `species_trait.tsv`
- `fractionation_bias_pairs.tsv`: scheduler-array manifest for pairwise
  `compare` and within-genome `self` kfFractBias analyses
- `species_rnaseq/`, `query_sra_id/`, `amalgkit_metadata/`: transcriptome-generation inputs
- `query_gene/`: per-family query files for query2family runs
- `input_generation/`: manifests and trait-plan inputs for input-generation runs

### `workspace/output/`

This is the main results area. Major stage roots include:

- `workspace/output/input_generation`
- `workspace/output/transcriptome_assembly`
- `workspace/output/species_cds_annotation`
- `workspace/output/species_cds_cdskit_localize`
- `workspace/output/species_tree`
- `workspace/output/orthofinder`
- `workspace/output/genome_evolution`
- `workspace/output/query2family`
- `workspace/output/orthogroup`
- `workspace/output/csubst_site`
- `workspace/output/{query2family,orthogroup}/csubst_scan`
- `workspace/output/{query2family,orthogroup}/csubst_scan_units`
- `workspace/output/{query2family,orthogroup}/csubst_scan_plot`
- `workspace/output/versions`

`workspace/output/input_generation/gg_input_generation_species.tsv` is the rolling species summary table
from `gg_input_generation_entrypoint.sh`. When taxonomy metadata is available, it includes per-species
`taxid` plus nuclear, mitochondrial, and plastid genetic-code columns. It also records the GFF repair
mode/status, audit path, repaired gene/reference counts, and unresolved ambiguity/collision counts.
When provided CDS and GFF files are grouped together, it records the grouping
mode/source, the adjacent `*.gff-grouping.tsv` audit path, mapped/unmapped/
ambiguous CDS counts, and coordinate-rescue counts.
Formatted source GFF files may have neighboring `*.gff.gz.repair.json` audit files; raw downloaded GFF
files remain unchanged.

When localization prediction is enabled, per-species CDS tables are written to
`workspace/output/species_cds_cdskit_localize/`, and per-family tables are
written under `workspace/output/query2family/cdskit_localize/` or
`workspace/output/orthogroup/cdskit_localize/` depending on the active mode.

### `workspace/downloads/`

This is the reusable runtime cache area.

Typical contents include:

- downloaded databases,
- provider download caches,
- temporary working directories,
- lock files for array-safe downloads,
- `ete_taxonomy/` for ETE taxonomy data,
- `pfam/` for Pfam/RPS-BLAST assets,
- `trait_datasets/` for cached trait database downloads.

## Runtime path helpers

Runtime helpers in `workflow/support/gg_util.sh` define canonical paths:

- input root: `workspace_input_root()`
- output root: `workspace_output_root()`
- downloads root: `workspace_downloads_root()`
- taxonomy root: `workspace_taxonomy_root()`
- Pfam root: `workspace_pfam_root()`
- taxonomy DB file: `workspace_taxonomy_dbfile()`

These helpers are what core scripts use internally; hard-coding alternate layouts is likely to break stage interoperability.

## Ownership and lifecycle

As a rule of thumb:

- treat `workspace/input/` as human-curated inputs,
- treat `workspace/output/` as reproducible stage results,
- treat `workspace/downloads/` as disposable-but-reusable cache/state.

For the newer protein-mode workflow, this distinction matters:

- `workspace/input/species_protein/` is curated input that you provide on purpose,
- `workspace/downloads/tmp/species_protein/` is temporary derived state created from `species_cds`,
- the temporary directory is regenerated when GeneGalleon detects input-sequence-mode or genetic-code changes.

That means:

- editing files under `workspace/output/` by hand is usually a bad idea,
- cleaning selected caches under `workspace/downloads/` is reasonable when repairing stale downloads,
- copying final curated inputs from `workspace/output/input_generation/` into long-term project storage is often useful.

## Summary and report files outside `output/`

Some wrappers intentionally write reports at the workspace root rather than under `output/`.

Current examples:

- `workspace/output/gene_summary/query2family/query2family_presence_absence.tsv`
- `workspace/output/gene_summary/query2family/query2family_presence_absence.pdf`
- `workspace/output/gene_summary/query2family/query2family_presence_absence.svg`
- `workspace/output/gene_summary/query2family/query2family_reference_gene_orthologs.pdf`
- `workspace/output/gene_summary/query2family/query2family_reference_gene_orthologs.svg`
- `workspace/output/gene_summary/query2family/query2family_reference_gene_orthologs.columns.tsv`
- `workspace/output/gene_summary/query2family/query2family_reference_gene_orthologs.glyphs.tsv`
- `workspace/output/gene_summary/query2family/query2family_reference_gene_orthologs.tree.tsv`
- `workspace/output/gene_summary/query2family/query2family_reference_gene_orthologs.synteny.tsv`
- `workspace/output/gene_summary/query2family/query2family_reference_gene_orthologs.ufboot.tsv`
- `workspace/output/gene_summary/query2family/query2family_query_gene_orthologs.pdf`
- `workspace/output/gene_summary/query2family/query2family_query_gene_orthologs.svg`
- `workspace/output/gene_summary/query2family/query2family_query_gene_orthologs.columns.tsv`
- `workspace/output/gene_summary/query2family/query2family_query_gene_orthologs.glyphs.tsv`
- `workspace/output/gene_summary/query2family/query2family_query_gene_orthologs.tree.tsv`
- `workspace/output/gene_summary/query2family/query2family_query_gene_orthologs.synteny.tsv`
- `workspace/output/gene_summary/query2family/query2family_query_gene_orthologs.ufboot.tsv`
- `workspace/output/gene_summary/query2family/query2family_query_gene_orthologs.query_map.tsv`
- `workspace/orthogroup_summary.tsv`
- `workspace/query2family_summary.tsv`
- `workspace/transcriptome_assembly_summary.tsv`

The `workspace/output/gene_summary/*` reports are produced by
`gg_gene_summary_entrypoint.sh`; gene-family figures combine a species tree,
branch support/dated-tree uncertainty when available, a detected/undetected
gene-family matrix, and per-species BUSCO stacked bars when BUSCO full tables
or a BUSCO summary table are available. For large orthogroup sets, full TSV
matrices are retained while `.plot.*` files record the default or user-selected
subset used for the PDF/SVG. Query2family summaries additionally expand each
family into every gene from `GG_COMMON_REFERENCE_SPECIES`; a light multi-column
glyph represents a copy inferred to predate the duplication of those reference
genes, and every detected glyph
contains its copy count. The compact reference-gene tree above those columns is
derived from the reconciled family tree. Duplication nodes use family-colored
circles in the compact query trees and matching family-colored bars in the
species-tree mapping.
When `presence_absence_ortholog_basis=query_gene` or `both`, the parallel
`query2family_query_gene_orthologs.*` artifacts instead use unique mapped query
anchors as columns. `query_map.tsv` preserves all original query records when
several resolve to the same tip.
The species tree receives every reconciled duplication from the complete gene
tree, not only duplication nodes in the reference-gene subtree, using
`spnode_generax` when available and `spnode_coverage` otherwise. The
`.tree.tsv` records the target as `mapped_species_node`; `in_reference_tree`
distinguishes nodes used in the compact reference-gene tree from mapping-only
duplication rows. Events from one family mapped to the same species-tree branch
are collapsed into one bar whose height is proportional to their count, while
bars for different families are distributed along the branch. Both the
color key and reference-gene tree titles use the input query filename directly.
For dated display trees, reconciled internal-node labels are transferred from
the matching undated species tree by descendant clade before the duplication
rows are mapped.
Column labels use gene IDs verified against `cds_fasta`, and multiple families
are isolated into titled blocks with separate query-tree topologies.
By default, a thin upper viridis band records the largest per-copy
local-synteny anchor count separately in each reference-gene column. A thin
lower viridis band runs across the complete ortholog glyph and records the Gene
tree UFBoot value for its single non-root, orthology-defining speciation-MRCA
branch. All candidate/reference rows expanded from one glyph must agree on that
branch, availability state, and support value; producers and plotters reject
inconsistent tables. Reference-self glyphs show white upper and lower bands
with a short horizontal state mark. The optional
`glyph` layout instead uses an upper-right viridis
diamond and a lower-right Inferno circle. Pair-level values and explicit
unavailability reasons are stored in `.synteny.tsv` and `.ufboot.tsv`,
respectively; the latter also records whether `support_generax_ufboot` or the
generic `support_unrooted` fallback was used. `stat_branch` stores the explicit
post-GeneRax percentage in
`support_generax_ufboot`, mapped by unrooted split onto the rooted GeneRax
branch-table topology.
The workspace-root summary TSVs are produced by
`gg_progress_summary_entrypoint.sh`.

## Common relocation scenario

If you want to keep code and data separate:

```bash
gg_workspace_dir=/data/projects/my_run/workspace \
gg_container_image_path=/data/containers/genegalleon.sif \
bash workflow/gg_genome_evolution_entrypoint.sh
```

In that setup:

- the checked-out repository still provides `/script`,
- the external workspace becomes `/workspace`,
- outputs and caches stay outside the git tree.
