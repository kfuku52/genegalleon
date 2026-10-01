# Dated-tree presentation

The dated-species-tree workflow adds geological-period backgrounds using the
ICS International Chronostratigraphic Chart [2026/06](https://stratigraphy.org/ICSchart/ChronostratChart2026-06.pdf).
Boundaries are plotted on the same absolute Ma axis as the saved tree, without
changing ages or intervals. Precambrian is shown as one broad interval. The
bundled table has a light presentation palette; numerical boundaries are source
data and do not represent new dating constraints.
Full period names appear vertically above their corresponding background bands
in the tree panel, without a separate geological legend. Neogene uses an ochre
background and Quaternary uses blue so the adjacent bands remain distinguishable.
On deep-time trees, crowded names spread within the header with short leaders
to their original bands, preserving the time scale.

The standalone helper also makes a dated-tree/BUSCO composite:

```bash
python workflow/support/plot_dated_tree.py \
  --infile /data/species_tree/mcmctree_95CI.nwk \
  --outfile /data/dated_species_tree_busco.svg \
  --busco-summary /data/annotation_summary/annotation_summary.tsv \
  --geological-background period \
  --layout-report /data/dated_species_tree_busco.layout.json
```

Run this inside the GeneGalleon container. The legacy two-argument CLI remains
available; `--geological-background none` selects the presentation renderer with
an unshaded background. PDF and SVG remain vector outputs, and PNG is supported
for previews. The default font is Helvetica 8 pt. The tree's Ma axis and the
BUSCO gene-count axis share the bottom edge of the panels, while BUSCO percentages
remain at the top. Geological names sit 4 pt above the tree panel. Age intervals
are opaque and drawn behind tree branches. Species rows use 9 pt spacing with the default
font, and the figure height is computed from the rows and header/footer text.
`--row-spacing` sets spacing in points; `--figure-height` overrides the computed
height. Use `--figure-width`, `--font-family` and `--font-size` for publication
requirements. The default width is 4.8 inches; widths down to 3.6 inches
are supported when the labels and panels fit. The renderer measures species
labels to reserve their physical width and splits a crowded legend into two rows.

BUSCO lineage identity is read from `busco_cds_lineage` in the summary, or from
BUSCO's own full/short-result headers. Canonical neighbouring
`species_cds_busco_short`/`species_cds_busco_full` directories are discovered
automatically; use `--busco-results DIRECTORY` for a copied legacy summary.
The axis reads `Number of BUSCO genes` followed by `(embryophyta_odb12)` on a
second line when that dataset is recorded for every plotted species. Unknown
identity is left unlabelled rather than inferred from the gene count. Mixed
datasets or unequal totals are rejected because a shared completeness axis
would be misleading. `annotation_summary.r` now preserves per-species lineage
in `busco_cds_lineage`/`busco_genome_lineage` and includes it in its BUSCO plots.

Highest posterior density intervals are spelled out in the legend, including
their probability level. Equal-tailed or unspecified credible intervals retain
their corresponding descriptions. Age labels are optional: `--node-ages all`,
`--node-ages root`, or `--age-clades FILE` (TSV column `descendant_species` with
comma-separated exact clades). Text is positioned around its own node; distant
labels have thin leader lines. A crowded figure fails with a request for more
space or fewer labels rather than publishing overlapping text.

`--tip-order FILE` takes a TSV column `species_id` in top-to-bottom order,
which must be compatible with the topology. `--tip-annotations FILE` takes
`species_id`, `colour` and `font_weight`. These options change presentation
only. The layout report records ages, intervals, geological boundaries, BUSCO
counts, dataset identity and tip-row coordinates. The plot and requested report
are staged together, with rollback on a handled publication error.
