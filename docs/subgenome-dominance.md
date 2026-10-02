# Subgenome retention and expression contrasts

An optional `gg_genome_evolution` stage consumes completed synteny/phasing and
expression results. Set `genome_evolution_mode="subgenome"` for an independent
run, or `run_subgenome_dominance=1` in `all` mode. It does not assign subgenomes
from gene number or expression, and does not rerun quantification or synteny.

`subgenome_manifest` defaults to `workspace/input/subgenome_analyses.tsv`.
Its TSV columns are `analysis_id`, `species`, `mapping_file`, and optional
`retention_file`, `homoeolog_file`, `expression_file`, `samples_file`,
`mapping_validated` and `retention_validated` (both 0 by default). Input paths are absolute or relative to the
manifest. Analysis IDs must be unique safe path components.
`expression_unit` defaults to `TPM`; `FPKM` is also accepted for within-sample
copy ratios. Other units, including counts and log transforms, are rejected.
Optional `reference` and `pair_set` columns describe the input dataset for
comparison-figure filters. They are metadata, not instructions to reconstruct
callability or change homoeolog eligibility.

Input tables:

| File | Required TSV columns | Contract |
| --- | --- | --- |
| mapping | gene_id, group_id, subgenome, assignment_scope, assignment_basis, evidence | One row per gene. Scope is local or global. Basis is synteny, phylogeny, kmer or curated, supported independently of retention/expression. At least two labels per group. |
| retention | group_id, block_id, locus_id, subgenome, callable, retained | One row per ancestral/reference opportunity per expected subgenome, including explicit uncallable rows. Boolean values are 0/1. All sides of a locus share one nonoverlapping block. |
| homoeolog | pair_id, block_id, gene_id | One independently inferred homoeolog gene per subgenome per set. Genes cannot occur in multiple sets. Tandem copies and ambiguous many-to-many assignments must be resolved upstream. |
| expression | gene_id, sample columns | Raw, finite, nonnegative TPM or FPKM for the mapped representative transcript. Do not use OG sums, log-transformed values, counts without effective-length correction, or inconsistent annotations. |
| samples | column, tissue, biological_id | Columns must match the expression table. Technical replicates sharing a tissue/biological_id are averaged. Tissue-specific biological IDs are required. |

Native kfFractBias `*.genes.tsv` is an intermediate, not a callable-opportunity
table: zero outside the expected collinear interval is not evidence of loss.
Define matched ancestral intervals independently; exclude assembly gaps and
unresolved annotation/alignability. Compare the same reference loci on both
sides, using outgroup and window-size sensitivity analyses. Self-synteny
retention alone is conditional on surviving genes and cannot establish loss.
Set `retention_validated=1` only after those callability and annotation checks.
Otherwise the result is explicitly exploratory syntelog detection, and a zero
must not be described as proven ancestral gene loss.

Each local group is analysed separately. Labels A/B in distinct local groups
must not be pooled or interpreted as parental subgenomes. `global` requires
independent phasing across groups; it adds `ALL_GROUPS` contrasts using pooled
loci with group-specific resampling blocks. Ancestral locus IDs must be unique
across groups for global retention. Local assignments leave genome-wide
dominance untested. More than two labels are
supported by pairwise contrasts. Display ordering of a synteny plot is not a
validated assignment.

For retention the effect is the paired difference in retained fractions.
For expression it is the mean log2(TPM_A/TPM_B), after averaging technical
replicates and then biological-sample log ratios. Both copies must be positive
in every biological sample of a tissue; exclusions are recorded in
`expression_coverage.tsv`. No pseudocount is added. The separate
`expression_detection_difference` metric compares the paired fraction of
biological samples with positive expression. Missing gene IDs remain unassessable; genuine
numeric zeros contribute to detection, not a log ratio. Analyse
mapping ambiguity separately, using unique diagnostic reads or simulation
before setting `mapping_validated=1`.

Statistics resample **nonoverlapping independently defined blocks**, not
overlapping sliding windows. At least three blocks are required for a 95%
cluster bootstrap interval and a two-sided block sign-flip test.
Integer block sums use exact subset-sum dynamic programming when the integer
sum range, reduced by its greatest common divisor, fits
`subgenome_exact_max_states` (default 262144). Other cases with up to 16 blocks
use exhaustive enumeration; larger tests use exact meet-in-the-middle when
each half fits that state ceiling. Remaining cases use
`subgenome_permutation_replicates` Monte Carlo draws (default 100000), with the
plus-one p-value correction. The null model is the same in all four methods.
Monte Carlo results include a 95% Wilson interval for the sampled null tail
probability, distinct from the biological effect CI. This interval does not
propagate Monte Carlo uncertainty through the multiple-testing correction.
Benjamini–Hochberg adjustment is within each metric/analysis across groups,
contrasts and tissues; each row records that family's number of estimable tests.
These intervals describe loci conditional on the sampled tissues, annotation
and mapping; they do not measure between-population biological uncertainty.
Block independence/exchangeability is a scientific assumption to check upstream.
Non-significance is not evidence that dominance is absent.

`subgenome_bootstrap_replicates` (default 2000) controls only the CI.
`subgenome_seed` (default 1) is expanded by SHA256 into separate bootstrap and
permutation streams keyed by analysis ID, metric, group, both copy labels and
tissue. Canonical block order makes results independent of input row order.
Removing retention, adding another contrast or changing bootstrap replication
does not consume another comparison's permutation stream. Changing analysis
IDs intentionally changes its streams; q values can change when the actual
test family changes. Individual stream seeds and inference version 2 are
recorded. Version 2 preserves the effect and null definitions but can change
CIs and p/q values relative to the former shared-stream, 2000-draw version.
Recompute the whole metric family when upgrading; retain prior published runs
as separate versioned artifacts.

Outputs under `workspace/output/genome_evolution/subgenome_dominance/<analysis_id>`:
`statistics.tsv`, `expression_pairs.tsv`, `expression_coverage.tsv`,
`summary.json`, and PNG/SVG/PDF contrasts when estimable. The stage retains input
SHA256s, seed, NumPy version, input/implementation provenance, and uses the
existing transaction lock, stale-policy and atomic publication mechanisms.
Missing retention or expression data remain explicit `not_estimable` results.

`contrasts_absolute.png` and `.svg` also show A/B-invariant magnitudes: the
absolute retention/detection difference or absolute **mean** log2 expression
ratio, rather than a mean of absolute per-gene ratios. Confidence intervals are
images of the signed 95% intervals under the absolute-value transform; an
interval spanning zero has a lower bound of zero. These plots discard direction
and do not show whether retention and expression favour the same copy.
The signed statistics and their P/q values are retained unchanged by plotting.

## Comparison figures and plot configuration

The run root also contains `comparison.{png,svg,pdf}` and
`comparison_absolute.{png,svg,pdf}`. Columns follow manifest order by default;
each metric has common x limits across columns. Multiple analyses of one species
receive separate columns labelled with the analysis ID. Points represent a
**group-specific subgenome contrast**, not a window. Multiple comparisons in a
group are vertically offset, with uniform circles; labels within local groups
do not imply a shared parental identity across groups. The points TSV preserves
copy labels, sample/block counts, signed statistics and displayed intervals.

Set `subgenome_plot_config` to a JSON file. Its display filters apply to the
combined comparison figure; individual analysis plots retain every estimable
comparison of the selected metrics. Filters never recompute p/q values or
reduce the original Benjamini–Hochberg family.

```json
{
  "font_family": "DejaVu Sans",
  "font_size": 8,
  "width_pt": 650,
  "height_pt": 320,
  "formats": ["png", "svg", "pdf"],
  "point_colour": "black",
  "error_bar_colour": "black",
  "significance_colour": "#D55E00",
  "reference": "Beta",
  "pair_set": "all_pairs",
  "tissue": "leaf",
  "metrics": ["retention_difference", "expression_log2_ratio"],
  "species_order": ["Ancistrocladus_abbreviatus", "Triphyophyllum_peltatum", "Nepenthes_gracilis"],
  "contrast_anchors": {"Nepenthes_gracilis": "D"}
}
```

The species names, Beta reference, leaf and subset above are an example,
not global defaults. `species_order` must name every selected species once;
`analysis_ids` optionally selects particular analyses. `reference`/`pair_set`
match analysis metadata; `tissue` filters expression rows. `contrast_anchors`
selects comparisons containing that label for each named species and does not
infer which copy is dominant. `group_labels` maps species to group/label objects.
`x_limits` and `absolute_x_limits` separately map metrics to `[lower, upper]`
in displayed units: percentage points for fraction differences and log2 units
for expression. Explicit limits must contain zero and all displayed CIs.

Style options also include `point_size`, `error_bar_width`, `capsize`, `dpi`,
`show_significance` and `panel_labels` (both default true), and
`significance_threshold` (default 0.05). `width_pt`/`height_pt` size the combined
figure; individual plots use automatic dimensions unless
`individual_width_pt`/`individual_height_pt` are specified.
The star legend spells out Benjamini–Hochberg. All text uses the configured
family and point size, including the star and italic species titles. Physical
figure dimensions are in points (72 pt per inch); PDF export preserves them
without tight-crop scaling. Downstream resizing changes the effective text size.
The renderer checks actual fonts, text sizes and canvas bounds, and records font
and output hashes in `*_provenance.json`. Too-small canvases fail with a request
to increase dimensions rather than silently shrinking text.

For Helvetica, set `font_family` to `Helvetica` and provide locally licensed
regular, bold and italic font files through `font_files`, an array of paths
relative to the config file or absolute paths. These files must be accessible
inside the container. No font binaries are distributed with GeneGalleon.
Missing families fail rather than silently falling back. PDF embeds the resolved
fonts; SVG retains editable text and needs those fonts on the viewing system.
Config and supplied font files are part of the stage's hashed input contract.

Existing statistics can be rendered without rerunning inference:

```bash
python workflow/support/subgenome_dominance.py report \
  --manifest comparison_manifest.tsv --plot-config plot.json --output comparison_figures
```

The comparison manifest requires `analysis_id`, `species`, `statistics_file`
and accepts context columns such as `reference`, `pair_set`, `expression_unit`,
`assignment_scope` and `retention_status`. Paths are relative to that manifest.
Both signed and absolute figures, points tables and provenance are exported;
input file hashes are checked again after rendering. Version 1 statistics tables
remain readable. The report command does not recalculate statistical inference.

Reference: [Saul et al. (2023), Nature Plants](https://www.nature.com/articles/s41477-023-01562-2).
Its retention and expression evidence should be assessed separately; a gene
copy with higher expression is not automatically the dominant subgenome.
