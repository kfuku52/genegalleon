# Pairwise genome synteny

`gg_genome_evolution` can infer gene-based synteny between two species with
JCVI/MCscan and DIAMOND, then draw a chromosome ribbon plot and a dotplot.
The ribbon plot serves the chromosome-comparison use case of a GENESPACE
riparian plot; it does not construct GENESPACE synteny-constrained orthogroups
or pan-gene sets. No species tree, BUSCO database or OrthoFinder run is required.

## Inputs and first run

Provide one matching annotation release per species in `workspace/input/species_gff`
and protein FASTA files in `species_protein`, or CDS FASTA files in `species_cds`.
Use the [species naming conventions](input-conventions.md). `auto` prefers a
protein source for each species and translates CDS where no protein source exists.
CDS translation uses `genetic_code` and the existing per-species
`species_genetic_code/species_genetic_code.tsv` overrides. Non-triplet CDS,
internal stops, duplicate FASTA IDs, ambiguous ID aliases and conflicting
annotation coordinates are errors.

The planning CLI's `--representative-map` and `--representative-inputs` paths
are optional. Omission uses ordinary inputs; an explicit path must still name
a valid file and retains the existing hash-bound bundle checks.

Create `workspace/input/synteny_pairs.tsv` with these required columns:

```tsv
analysis_id	target_species	query_species
triphyophyllum_ancistrocladus	Triphyophyllum_peltatum	Ancistrocladus_abbreviatus
```

The same row is supplied in `workspace/input/synteny_pairs.example.tsv`.
`analysis_id` must be unique. Both species must differ. Each row is one comparison;
rows run sequentially inside this single-task workflow using its allocated CPUs.

```bash
GG_GENOME_EVOLUTION_GENOME_EVOLUTION_MODE=synteny \
  bash workflow/gg_genome_evolution_entrypoint.sh
```

`genome_evolution_mode=synteny` enables this stage and exits before species-tree
and orthogroup setup. Existing outputs from those stages are preserved. In
ordinary `all` mode, set `run_pairwise_synteny=1` to add the stage to the pipeline;
its default is `0`. Scheduler resource directives are unchanged.

Optional pair-table columns are `target_fasta`, `target_gff`, `query_fasta`,
`query_gff`, `target_feature`, `target_attribute`, `query_feature`,
`query_attribute`, `target_seqids`, `query_seqids`, `target_cds`, `query_cds`,
`target_genome`, `query_genome`, `target_sizes`, and `query_sizes`.
CDS columns supply optional paths for dS coloring. Paths are absolute or
workspace-relative. Multiple same-species source files require explicit paths;
explicit FASTA paths also require `synteny_sequence_mode=protein` or `cds`.
GFF feature/attribute overrides must be supplied together. Seqids are exact,
comma-separated chromosome names in the desired ribbon display order.
Unlisted chromosomes remain in the analysis and can appear in the dotplot if
they meet its physical-length filter.

The dotplot includes only scaffold/chromosome lengths **at least 1,000,000 bp**
by default (`synteny_dotplot_min_length=1000000`). Supply genome FASTA paths
in `*_genome`, or chromosome-size tables in `*_sizes` (two whitespace-separated
columns: exact seqid and positive length; `.fai` files also work). Do not supply
both for the same side. Otherwise a unique matching FASTA in `species_genome`
is used, or GFF `##sequence-region seqid 1 length` declarations. GFFs must
declare the whole sequence span for this use; partial regions need assembly
FASTA/size inputs. Missing lengths are errors, never estimated from gene counts
or the last annotated gene; they are checked before dS alignment starts.
Set the minimum to `0` explicitly to disable this
display-only filter. The analysis and dS estimates still cover all anchors.

`synteny_dotplot_sort=homoeolog` (default) places supported 2×2 chromosome
groups next to each other on both axes. Groups must have unique anchors in all
four chromosome comparisons; a deterministic greedy search prioritizes their
harmonic-mean support, with ribbon order breaking ties. Weak cross-connections
remain visible. This is a display heuristic, not a biological homoeology/WGD
test or a globally optimal grouping. `karyotype` follows the ribbon order,
appending unlisted eligible chromosomes in BED order; `none` preserves BED
order. Dotplot modes do not reverse chromosomes or change within-chromosome gene order.
`dotplot_order.json` records the groups, lengths, order and excluded anchor
counts; `dotplot_display.tsv` records each chromosome's length and selection.

The default `synteny_karyotype_sort=both_length` reorders both ribbon tracks
to reduce the sum of **ribbon width × connection length**. `target_length`
fixes the query order and reorders the target; `query_length` fixes the target
and reorders the query. Width is the mean of the two endpoint widths in JCVI
gene-rank display coordinates; length
is the straight line between their centers, measured on the 20×8-inch canvas.
All displayed blocks contribute, including secondary connections. Chromosome
gene counts, JCVI gaps, and block positions within chromosomes are included,
so this is not simply a strongest-partner ordering. The objective does not
minimize crossing count or the arc length of the decorative Bézier curves.

For either fixed-track mode with up to 20 moving chromosomes, exact subset
dynamic programming finds the global minimum of this fixed-track objective,
to floating-point precision.
The start position depends only on the widths of the preceding subset, allowing
the solver to avoid enumerating all permutations. Exact cost ties prefer input
order. For more than 20 chromosomes, bounded adjacent-swap improvement avoids
exponential resource use; this is explicitly recorded as **not globally optimal**.
Use an explicit main-chromosome seqid list when exact optimization is needed.

For `both_length`, small problems enumerate one track and solve the other
exactly (at most 2,000,000 combined work units). Larger problems, including
18×18 chromosomes, use six deterministic multistart alternating searches;
each track update is exact when it has at most 20 chromosomes. Joint global
optimality is **not guaranteed** for this bounded search. The reported result
cannot be worse than the initial one-track optimum. Inspect the recorded
objective and solver status rather than interpreting the order as ancestry.

`target` and `query` retain the earlier dominant-partner heuristic: unique lifted
anchor-pair counts choose the strongest fixed partner, then mean partner gene
rank and input order break ties; chromosomes with no displayed partner go last.
`none` preserves the supplied lists or natural order. In these legacy modes the lists
select the displayed chromosomes; all ribbons, chromosome orientations, gene
order and analysis are unchanged. The dotplot can follow this ordering using
`synteny_dotplot_sort=karyotype`. `karyotype_order.json`
records the input/final orders, solver, objective values and optimality status
(or anchor support for the dominant-partner heuristic).

For Triphyophyllum (target) versus Ancistrocladus (query), list `scaffold1`
through `scaffold18` explicitly in the pair table's `query_seqids` and use:

```bash
GG_GENOME_EVOLUTION_GENOME_EVOLUTION_MODE=synteny \
GG_GENOME_EVOLUTION_SYNTENY_KARYOTYPE_SORT=both_length \
  bash workflow/gg_genome_evolution_entrypoint.sh
```

FASTA IDs match GFF attributes exactly or after removal of the exact species
prefix. Ambiguous matches fail. The annotation mapper bundled in kfFractBias
resolves known gene loci; the longest available protein/CDS per locus is selected,
with identifier ordering breaking length ties. Unknown locus relationships are
reported in the summary instead of being guessed. All original FASTA IDs,
selected representatives, normalized JCVI IDs, excluded isoforms and unmapped
records remain in the ID mapping tables. The default required mapping fraction
is `1`; inspect source-release and identifier mismatches before changing it.

## Analysis and outputs

Default JCVI filters are `synteny_cscore=0.7`, `synteny_min_anchors=4`, and
`synteny_search_distance=20`. There is no quota filter: many-to-many blocks
remain available, including duplications. The result alone does not establish
whole-genome duplication or ancestral chromosome structure. Tandem-duplicate
filtering follows JCVI's algorithm. Both seed and lifted anchors are retained.
Karyotypes include a compact ribbon-symbol legend using the saved analysis
parameters and execution log, not current/default settings. It reports the
DIAMOND E-value threshold, seed C-score, tandem distance, minimum unique genes
in each genome, maximum gene-rank chaining gap and liftover distance. The
chaining distance is not a fixed sliding window or a limit on total block size.
Liftover recruits original protein hits by Manhattan gene-rank distance, without
reapplying the seed C-score filter. Ribbons have no dS filter, and only blocks
whose endpoints are on displayed chromosomes are drawn. `karyotype_style.json`
records the legend text, criteria and source hashes. Standalone rendering can
supply `--analysis PATH`; without provenance it does not invent threshold values.

Outputs are under `workspace/output/genome_evolution/synteny/`:

| Directory | Contents |
| --- | --- |
| `analysis/<analysis_id>/` | BED and protein inputs, ID maps, seed/lifted anchors, `anchors.tsv`, `blocks.tsv`, `summary.json`, commands and logs |
| `ds/<analysis_id>/` | Optional CDSKIT `ds.tsv`, audited codon alignments (`aligned_pairs.tsv.gz`) and method/source summary |
| `plots/<analysis_id>/` | `karyotype` and `dotplot` in PDF/SVG/PNG, `display.tsv`, reproducible layout/seqids and copied plotting inputs |

Both plots use **gene rank**, rather than physical chromosome length. BED/block
coordinates use 0-based, half-open intervals. Ribbons are colored by the target
chromosome; dotplot colors distinguish block orientation by default. Chromosome labels and
their original orientation are preserved by default; natural ordering puts Chr2 before
Chr10. Display exclusions are recorded, and the dotplot retains every anchor
between its eligible chromosomes without downsampling. The 1 Mbp filter does
not change ribbon selection; fragmented assemblies may still need an explicit
ribbon seqid list for readability. Full plotting inputs are retained alongside
separate filtered/reordered `dotplot.*.bed` and `dotplot.anchors` files.

Dotplot PDFs have a total page width of 3.6 inches, including labels and margins,
and a square physical plot area (not equal x/y data units). Numeric major ticks
point outwards from the bottom and left edges; their labels stay outside the
plot area. X-axis chromosome labels are vertical; all text is 8 pt Helvetica,
with italic species names and upright `(gene rank)` suffixes. The `dS` token in
the dotplot title and colourbar label is also italic; the surrounding text
remains upright.
The layout is fitted without scaling the text. Both plots use black 8 pt
Helvetica text and restrained soft palettes. Karyotypes italicise only species names.
SVGs retain editable Helvetica text; raster previews use the available system
font substitute if Helvetica is unavailable. Karyotype PDFs have a 7.2-inch page
width and 0.06-inch outer padding. Original chromosome labels are vertical to
avoid crowding. Only the opt-in signed layout mode below reverses chromosome display directions.
Each track positions its species label 4 pt beyond the outer edge of its own
longest rendered chromosome/scaffold label. The measurement uses the actual
font metrics and rotated label bounds, so long labels on one track do not
create unused space on the other. Scale bars follow the species labels,
with the gene-count text above each bar.
All karyotype chromosomes have rectangular outlines and fills, including short
chromosomes; their ends are never rounded according to chromosome length.
Karyotypes use a **shared gene-count scale** by default
(`synteny_karyotype_scale=shared`): one gene occupies the same physical width
on both tracks, and a single scale bar serves both species. The largest track
fills the available width; the shorter one is left-aligned and is not stretched
to match it. Native chromosome gaps are preserved, so the scale measures
chromosome gene counts, not inter-chromosome whitespace.
Set `GG_GENOME_EVOLUTION_SYNTENY_KARYOTYPE_SCALE=independent` for the previous
per-species normalization with one scale bar per track. These are not Mbp
scales: both plots still use gene-rank coordinates.
The ribbon-length objective uses the normalized 20:8 design geometry; fitting
and cropping the page do not change the order or block endpoints. Its endpoint
widths and positions use the selected shared/independent scaling mode.

Chromosome ribbon colors are separate by default
(`synteny_karyotype_color=chromosome`). Set
`GG_GENOME_EVOLUTION_SYNTENY_KARYOTYPE_COLOR=homoeolog` to give each supported
2x2 group a shared color on both tracks and their ribbons. Grouping is based on
displayed chromosomes, independently of the dotplot's physical-length filter,
and does not reorder either track. This is opt-in display grouping, not a
biological homoeology call. `karyotype_colors.json` records the full color map;
`karyotype_style.json` records the scaling mode, scale unit/value, number of
scale bars, track gene counts, track ratios and per-track species-label placement.
Changing the scaling mode
invalidates only plots, not the analysis or dS estimates.

The target track is above the query track by default
(`synteny_karyotype_track_order=target-query`). Set
`GG_GENOME_EVOLUTION_SYNTENY_KARYOTYPE_TRACK_ORDER=query-target` to reverse the
display. This moves only the vertical track coordinates; BEDs, gene ranks,
chromosome order, colors, ribbon endpoints and analysis roles remain unchanged.
It invalidates only plots. The standalone `pairwise_synteny_karyotype.py` helper
accepts `--track-order query-target`; `pairwise_synteny.py plan` accepts
`--karyotype-track-order query-target`. The style report records top-to-bottom
`track_order` and `display_species_order`; per-track counts, ratios and label
metadata stay in target/query analysis order (`track_metadata_order`). Shared
scales place their single bar below the bottom displayed track.
Saved layout BED/block paths are resolved relative to the layout directory,
so `--directory` can be used from any working directory. Duplicate chromosome
or BED gene IDs and selected chromosomes absent from the BED are rejected.
Ribbon endpoints must belong to one chromosome per track, with a positive block
score and a `+` or `-` orientation.
The standalone CLI stages the plot and `karyotype_style.json` together and
restores both if installation fails; input layout, BEDs and blocks are protected
from output replacement.

The summary records source hashes, selected/unmapped gene counts, syntenic gene
fractions, block/anchor counts, genetic codes, algorithm parameters, tool versions
and annotation-mapper identity. A run with no blocks fails before publication.
Analysis, optional dS estimation and plots have separate standard artifact-provenance contracts and
recoverable bundle publication. A stage lock prevents concurrent runs from
mixing analysis and plots. Completed analysis is retained if a plot fails.

## Order and orientation, including three or more species

Set `synteny_karyotype_sort=order_orientation` to optimize **both the order and
display direction** of a pairwise karyotype. This is a different objective from
`both_length`: it gives equal weight to three normalized terms, ribbon width
times straight midpoint distance, weighted midpoint crossings, and the width
of ribbons whose endpoint directions disagree. Width is the mean endpoint span
in native JCVI gene-rank coordinates. Zero-width blocks remain present with zero
weight. Crossings exclude tied midpoints and are not intersections of Bézier
areas. All selected blocks contribute; no synteny or dS thresholds change.
Pairwise optimization uses the same 20:8 design geometry as the legacy solver;
multi-track optimization uses its actual rendered geometry.

For a single diagram with multiple species, provide an optional JSON file
through `synteny_karyotype_panels_file`. Every named comparison must already be
listed in `synteny_pairs.tsv`, for example:

```tsv
analysis_id	target_species	query_species
anc_dio	Ancistrocladus_abbreviatus	Dioncophyllum_thollonii
dio_tri	Dioncophyllum_thollonii	Triphyophyllum_peltatum
```

```json
{
  "panels": [{
    "panel_id": "anc_dio_tri",
    "tracks": [
      {"species": "Ancistrocladus_abbreviatus"},
      {"species": "Dioncophyllum_thollonii"},
      {"species": "Triphyophyllum_peltatum"}
    ],
    "comparisons": ["anc_dio", "dio_tri"],
    "mode": "order_orientation",
    "seed": 1,
    "restarts": 4,
    "max_rounds": 8
  }]
}
```

The `tracks` list fixes the top-to-bottom species order. A middle track has
**one shared order and orientation**, scored against every incident comparison.
Repeated species must have identical BEDs across the selected comparisons;
different annotation releases or coordinate systems are rejected. Add `seqids`
to each track to select readable chromosome/scaffold subsets; omission includes
every annotated sequence. Selection affects the diagram only.

Each track can also specify:

| Option | Meaning |
| --- | --- |
| `groups` | Disjoint chromosome lists that must remain adjacent, for example `[["scaffold3", "scaffold16"], ["scaffold11", "scaffold1"]`. Groups must be contiguous in the supplied `seqids` order; group order and member order may change together. No groups or homoeology are inferred by the optimizer. |
| `fixed_order` | Boolean, default `false`; retain this track's chromosome order. Directions may still change unless locked separately. |
| `orientation` | Initial display directions, for example `{"scaffold3":"-"}`. Unspecified sequences start forward. |
| `fixed_orientation` | Locked display directions, for example `{"scaffold3":"+"}`. Keys must be selected sequences; conflicting explicit initial directions are errors. |
| `colors` | Optional seqid-to-matplotlib-color map; identities remain the original seqids. |
| `labels` | Optional seqid-to-display-label map, for example `{"scaffold3":"3"}`. Labels must be unique within a track. BED IDs, saved seqids and optimization identities remain unchanged. |

Panel `mode=order` optimizes order while retaining initial directions;
`mode=none` retains the supplied layout. The `order_orientation` mode is the
panel default. No mode changes the species order, BEDs, sequence data, gene IDs,
alignments, or inferred block orientation. Reversed chromosome labels have
`(rev)` and the reproducible JCVI `seqids` file uses its trailing `-` display
control. Pairwise dotplots keep their original directions even when following
the reordered chromosome list.

At most 20,000 permitted joint order/direction states are explored exhaustively.
Larger problems use a seeded, bounded multistart search with weighted-center
order proposals, group-preserving swaps and direction changes. `restarts` is
bounded to 1–16 and `max_rounds` to 1–100. Crossings are evaluated using linear
storage rather than an all-block-pairs matrix. The returned objective cannot
be worse than the constrained starting layout to floating-point tolerance;
equal scores prefer fewer direction changes and the input order. The report
states `globally_optimal=false` for bounded search and for `mode=none` rather
than calling an unsearched layout optimal. This is a visualization criterion,
not evidence of ancestral orientation or structural rearrangement.

```bash
GG_GENOME_EVOLUTION_GENOME_EVOLUTION_MODE=synteny \
GG_GENOME_EVOLUTION_SYNTENY_KARYOTYPE_PANELS_FILE=input/synteny_panels.json \
GG_GENOME_EVOLUTION_SYNTENY_PLOT_ONLY=1 \
artifact_stale_policy=rebuild \
  bash workflow/gg_genome_evolution_entrypoint.sh
```

Outputs are `plots/panels/<panel_id>/karyotype.{pdf,svg,png}` with 7.2-inch width,
8 pt Helvetica and one shared gene-count scale by default. Plots use
8 pt gaps between chromosomes so adjacent rotated labels remain separated.
Select fewer chromosomes if their labels and gaps cannot fit the fixed width.
The bundle includes
unchanged source BEDs/blocks, relative-path `manifest.json`, native `layout` and
`seqids`, colors, and `karyotype_layout.json`. The report records input/final
orders and directions, constraints, geometry, seed, solver, before/after score
components, optimality status, source hashes and displayed/excluded block
counts. Pairwise signed layouts use the existing `karyotype_order.json`.
Panel settings invalidate plots only, preserving completed analysis and dS.
The pairwise `analysis_id=panels` is reserved only when panels are enabled.

For saved inputs outside a workspace, the same backend accepts a single
manifest containing `tracks` (each with `species`, `bed` and optional settings
above), `edges` (`{"tracks":[0,1],"blocks":"anc_dio.simple"}`), and optimizer
options. Paths resolve beside the manifest. A new destination is required:

```bash
python workflow/support/synteny_multitrack.py \
  --manifest saved_layout.json --output new_layout --formats pdf,svg,png
```

The generated bundle's `manifest.json` can be supplied to this command to
reproduce its layout without protein search, synteny detection or dS fitting.

## Dotplot colored by dS

```bash
GG_GENOME_EVOLUTION_GENOME_EVOLUTION_MODE=synteny \
GG_GENOME_EVOLUTION_SYNTENY_DOTPLOT_COLOR=ds \
GG_GENOME_EVOLUTION_SYNTENY_DS_COLOR_MAX=2 \
  bash workflow/gg_genome_evolution_entrypoint.sh
```

This requires a runtime containing CDSKIT's `dnds` command and MAFFT; an older
runtime fails explicitly. There is no alternate estimator fallback. CDS sources
must translate exactly to the selected synteny proteins. CDS-mode inputs are
reused; protein-mode inputs use `species_cds` or explicit `target_cds`/`query_cds`.
Both species must use the same genetic code for this pairwise estimator.
ID-map headers and original IDs must be unique: identical protein translations
cannot justify swapping distinct synonymous CDS variants. Malformed dS tables,
unknown result statuses, and mismatched report/codon-semantics versions are errors,
not missing biological evidence.

Each unique lifted-anchor gene pair is aligned by MAFFT in protein space and
back-translated with CDSKIT. Native batched YN00 uses equal path weights
(`weighting=0`), pair-specific F3x4 frequencies and kappa. Missing/gapped codons
are removed jointly. Alignment and estimator batch workers use the allocated
CPUs; BLAS/OpenMP threads within the dS process are limited to one per worker
to avoid nested parallelism. Unestimable or saturated distances are `NA` and
appear in gray, never as zero. A gray square and the missing-pair count sit
directly beside the colorbar at the same vertical centerline. All anchors on
eligible chromosomes remain visible,
including those above the color limit; the upper limit clips colors only.
The continuous dS palette runs from muted medium blue to dark blue (higher dS
is darker), with no pale yellow. Every colour in the scale has at least 3:1
contrast against the white background, keeping small points visible.
`dotplot_ds.json` records counts.
Missing chromosome lengths or no anchors passing the physical-length filter
stop the run before alignment/dS estimation; completed synteny analysis is retained.
These pairwise distances are not a codeml likelihood analysis or a WGD test.

## Combined pairwise dS dotplots

To compare all three species pairs in one figure, list the three comparisons in
the ordinary pair table and set `synteny_dotplot_panels_file` to a JSON file.
The workflow performs the same native analysis and CDSKIT dS stages for each
pair, then combines their gene-pair coordinates with **one shared dS colour
scale and missing-value key**. It does not rerun alignment to combine plots.
The [three-species pair table](examples/synteny-three-species.tsv) and
[panel configuration](examples/synteny-dotplot-panels.json) are complete examples.

```json
{
  "panels": [{
    "panel_id": "anc_dio_tri",
    "comparisons": ["anc_dio", "anc_tri", "dio_tri"],
    "tracks": [
      {"species": "Ancistrocladus_abbreviatus"},
      {"species": "Dioncophyllum_thollonii"},
      {"species": "Triphyophyllum_peltatum"}
    ],
    "columns": 2
  }]
}
```

```bash
GG_GENOME_EVOLUTION_GENOME_EVOLUTION_MODE=synteny \
GG_GENOME_EVOLUTION_SYNTENY_PAIRS_FILE=input/synteny_pairs.tsv \
GG_GENOME_EVOLUTION_SYNTENY_DOTPLOT_COLOR=ds \
GG_GENOME_EVOLUTION_SYNTENY_DOTPLOT_PANELS_FILE=input/synteny_dotplot_panels.json \
  bash workflow/gg_genome_evolution_entrypoint.sh
```

The `comparisons` order sets the subpanel order. `tracks` lists each participating
species once and fixes the visual axis roles: the earlier species is on x and
the later species is on y. This can swap axes relative to the analysis roles;
gene IDs, anchor pairs and dS estimates remain unchanged. Each species uses the
same sequence order and direction throughout the figure. By default its first
comparison supplies the eligible sequence order. Optional per-track `seqids`,
`orientation` and `labels` have the same meanings as the multi-track layout:

```json
{"species": "Ancistrocladus_abbreviatus", "seqids": ["scaffold16", "scaffold3"],
 "orientation": {"scaffold16": "-"},
 "labels": {"scaffold16": "16", "scaffold3": "3"}}
```

Explicit sequences must pass the ordinary dotplot physical-length filter.
Use the saved karyotype layout's `display_order`, `display_orientation` and
`display_labels` to populate these settings when matching an existing figure.
Reversals are marked `(rev)` and affect only this combined display. Native
scaffold boundaries remain separate; no chromosome or homoeolog assignments
are inferred. Repeated species must have identical full BED/CDS inputs and
eligible sequence sets across comparisons; inconsistent releases are errors.

`columns` defaults to two (one for a single comparison) and accepts 1–8.
Each column is 3.6 inches wide, retaining 8 pt Helvetica text. The common key
uses a spare grid cell when available; otherwise it sits below the panels.
Long/crowded sequence labels that cannot fit at 8 pt require explicit sequence
selection or shorter display labels. Counts report **unique** displayed anchor
pairs, including missing distances and values above `synteny_ds_color_max`.
All retained pairs are drawn without downsampling; the colour limit never
filters anchors. The ordinary per-pair dotplots retain their existing layout.

Outputs are `plots/dotplot_panels/<panel_id>/dotplot.{pdf,svg,png}` in the selected
formats, copied plotting inputs, `manifest.json`, and `dotplot_panels.json`.
The report records sequence layouts, axis swaps, pair/NA/clipped-colour counts,
the common scale and source hashes. `analysis_id=dotplot_panels` is reserved
only when this option is enabled. Changing panel layout or the colour limit
invalidates plots only: use the existing `synteny_plot_only=1` and
`artifact_stale_policy=rebuild` to reuse current analysis and dS results.

## Redraw and resume

An unchanged rerun reuses both analysis and plots. A changed display order or
output format invalidates only the plots. Sorting-setting changes also affect
only plots, as do separate assembly FASTA/size inputs and dotplot filter/order
settings. Changing the annotation GFF still invalidates analysis.
Use this after changing `*_seqids`, `synteny_karyotype_sort`, or
`synteny_karyotype_scale`:

```bash
GG_GENOME_EVOLUTION_GENOME_EVOLUTION_MODE=synteny \
GG_GENOME_EVOLUTION_SYNTENY_PLOT_ONLY=1 \
GG_GENOME_EVOLUTION_SYNTENY_PLOT_FORMATS=pdf,svg,png \
artifact_stale_policy=rebuild \
  bash workflow/gg_genome_evolution_entrypoint.sh
```

Plot-only requires a complete, current analysis and, for dS coloring, a current
dS stage. It never starts alignment or dS estimation. A changed dS color limit
invalidates only plots; changing CDS or estimator source invalidates dS and plots,
without rerunning an otherwise current protein-based synteny analysis.
Changed FASTA/GFF inputs or
analysis settings cannot silently use the previous analysis. The normal
`artifact_stale_policy=stop` stops on drift; an explicitly selected `rebuild`
regenerates the affected stage. Plot-only does not rebuild analysis. Review
`summary.json`, both figures and `display.tsv` before interpreting real data.

For plot-only, stale analysis is rejected even with `artifact_stale_policy=reuse`.
Input hashes are rechecked before publication to detect changes during a run.
This guard includes the pair table without adding its byte hash to every phase's
cache key, so a subsequent display-only table edit still reuses scientific analysis.
