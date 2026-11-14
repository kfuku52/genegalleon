# Native WGD Candidates and Duplication Origins

This experimental analysis infers species-branch candidates from gene counts;
known WGD events are not required inputs. It combines native NWKIT DL+genome
multiplication inference, kfFractBias unquota self-synteny, CDSKIT YN00 dS,
and NWKIT multiple-outgroup Ks correction. DupGen_finder, Whale, and ksrates
are not runtime dependencies.

## Inputs and Execution

Use an updated GeneGalleon runtime containing the new `nwkit wgd-count`,
`nwkit ksrate`, and kfFractBias raw self-evidence commands. Local development
validation can mount the owned source repositories read-only using
`GENEGALLEON_DOCKER_EXTRA_BINDS` and set `PYTHONPATH`; this does not update an
existing image or establish SIF compatibility.

In the editable block of `workflow/gg_genome_evolution_entrypoint.sh`, set
`genome_evolution_mode="wgd"`, then run that entrypoint in the usual runtime.
This standalone mode leaves existing species-tree and OrthoFinder outputs
untouched. Alternatively set `run_wgd_ssd=1` in `all` mode.

Default inputs are:

- `output/species_tree/species_tree_summary/dated_species_tree.nwk`, otherwise
  the existing `undated_species_tree.nwk`. The tree must be rooted, have unique
  tip names, and have nonnegative branch lengths. Count-model rates use these
  branch units; an undated input does not produce absolute WGD ages.
  Rootedness is read automatically: `[&U]` and root NHX `nwkit_rooted=no` or
  `unknown` are rejected, even for a binary top-level node. An unmarked binary
  root retains NWKIT's legacy rooted interpretation; a root polytomy needs an
  explicit `[&R]` or root NHX `nwkit_rooted=yes` declaration. These declarations
  identify the supplied root; the stage does not reroot or arbitrarily orient a tree.
- `output/orthofinder/Orthogroups/Orthogroups.GeneCount.tsv` and matching
  `Orthogroups.tsv`. Override with `wgd_counts_file` and `wgd_members_file`.
  Species columns must match the full tree. Count-dependent prior filtering
  is not modeled: do not use copy-number-selected families for confirmatory
  inference. Counts and memberships must represent one copy per gene locus,
  not alternative transcript isoforms. A missing count is not zero.
- Optional `input/wgd_genomes.tsv`. Without this file, statistically
  supported candidates can be reported from gene-family counts alone.
  Individual gene duplication origins remain unresolved.

The genome manifest contains unique species-tree tips, with these columns:

```text
species                 mode      fasta   gff   cds   feature   attribute
Arabidopsis_thaliana     cds
Oryza_sativa            protein
```

Use actual tabs. Only `species` is mandatory; `mode` defaults to the configured
CDS/protein mode. Empty FASTA/GFF fields use the existing exact-species input
discovery in `input/species_cds`, `species_protein`, and `species_gff`.
Explicit paths are relative to the workspace, or absolute. An explicit `cds`
path can supply matching CDS for protein inputs. No CDS means no Ks for that
species, not a fabricated zero. The species genetic-code table and configured
default code apply, as in other sequence stages. Contrasts between different
genetic codes are not estimated in this initial integration.

GFF identifiers must map unambiguously to FASTA IDs. The longest annotated
locus isoform is used; discarded isoforms and unmapped gene-tree tips do not
receive confident origin labels. CDS translations must exactly match the
selected proteins. Gene membership IDs must be the canonical species-prefixed
IDs used by the workflow and belong to their declared species column. Duplicate
count families, malformed counts, and inconsistent `Total` columns are errors.
Before fitting counts, supplied annotations are checked for multiple membership
IDs from the same locus, even when they were split across families or lack a
supplied sequence. Correct representative-locus membership/count inputs are
required; synteny's isoform selection does not silently rewrite family counts.
Without a genome annotation this contract remains the input provider's responsibility.

## Candidate Inference

Set `wgd_gene_trees_file` to a rooted GeneGalleon family-tree ZIP (members
ending in `_root.nwk`) to select the alternative duplication-score mode.
This mode uses canonical species-prefixed tips and scores internal non-root
branches without fitting the count model. A family is eligible when its tree
includes every daughter lineage and at least one species from the immediate
sister lineage. It supports the branch when a binary duplication spans that
exact clade, with every daughter lineage represented in both copy groups.
Each family votes at most once per branch. `wgd_min_duplication_score`
(default `0.01`) is the minimum supporting / eligible fraction. This is
an exploratory score, not a P-value, and the result depends on family
selection, taxon coverage, and gene-tree rooting. It does not classify
terminal events. At least one descendant species in every daughter
lineage must have observable synteny for a combined assessment.
With insufficient synteny opportunity, passing branches
are labeled `duplication-only-candidate`; with matching synteny, they are
`duplication-and-synteny-candidate`. Neither label assigns a gene origin.
This mode writes `duplication_candidates.tsv` instead of the count-model
outputs and cannot be used with the native per-gene tree likelihood.
Optionally provide a GeneRax UFBoot tree ZIP as `wgd_cross_tree_support_file`,
with or without cross-species synteny. `duplication_tree_support.tsv` counts
each voting family once per branch, separating absent trees, different leaf
sets, and missing exact split support. The last three columns count families
with at least one matching duplication whose own split and both copy splits
all have UFBoot support at or above 50, 70, or 95. These are descriptive
levels, not WGD significance thresholds or a test of the root position.

Set `wgd_duplication_taxon_holdout=1` to write
`duplication_taxon_holdout.tsv`. For each original branch and excluded species,
it reports eligible and voting families and gains and losses relative to the
full scan. A branch is unidentifiable if the excluded species was its only
representative in a daughter group or the immediate sister group. If no
families are eligible, the score is `NA`, not zero. The audit prunes the
archived topology in memory; it does not re-infer gene trees, change calls, or
provide a P-value.

Set `wgd_duplication_coverage_audit=1` to write
`duplication_coverage.tsv`. It splits eligible families by whether all
descendant species of the branch are present in their archived tree and
reports a complete-case score. For partially covered families, a separate
score counts duplications whose species set equals the branch's observed
descendants, with both copy groups still covering every daughter lineage.
The immediate sister lineage must still be represented. This sensitivity
check can recover votes lost to missing taxa, but presence-dependent
selection, gene-tree error, and sparse sampling can also change these
scores. It does not alter default calls or provide a P-value.

Set `wgd_duplication_root_audit=1` to write
`duplication_root_orientations.tsv`. For each family, examine which binary
duplications could appear by changing only the root direction on its
existing topology, without rearranging splits or resolving polytomies.
Report the number of extra eligible families and how many already vote on
a different branch under the archived root. The any-orientation score is a
descriptive upper bound on votes on one branch, not a preferred rooting,
a root confidence interval, or a replacement for the archived-root score.

`wgd_duplication_family_panels` optionally partitions the gene-tree archive
into 2 to 64 disjoint panels by a seeded hash of the family ID. This makes
the assignment independent of archive order and keeps all eligible and
supporting votes from a family together. The per-panel and summary tables
record sample sizes, branch scores, and how often the configured cutoff is
met. An untestable panel is not a negative. Panel agreement is not a
P-value, particularly for linked families, and does not change WGD calls.

### Sparse Cross-Species Synteny

In duplication-score mode, set `wgd_cross_synteny_file` to a TSV with
columns `species_a`, `species_b`, `bed_a`, `bed_b`, and `anchors`.
Paths are absolute or workspace-relative. Provide each genome's self
comparison (same species on both sides) and one or more direct cross
comparisons spanning the daughter lineages at each internal split.
Species and BED gene IDs must match the input species and gene trees.
BED gene IDs must be globally unique, and each archive member must name
a family in the membership table and contain only its genes (subsets allowed).
For this mode the membership column supplies the species identity, including
when one species name is a prefix of another.

`duplication_synteny_family_overlap.tsv` reports how many distinct families
contribute both a duplication vote and a matched cross-synteny path. It also
shows the score with those families removed from the numerator and denominator.
The subtraction uses all matched paths even in the strict-phase sensitivity
run. Exclusion is conditioned on the evidence and is not an independent
test. The default scores and calls are unchanged.
Use one BED per species and original `###`-delimited JCVI anchor files
(plain or gzip), not lifted anchors or preassigned branch scores. This
stage reads existing comparisons; it does not launch genome alignments.
The input tree assigns each cross pair to its MRCA; a root comparison
is validated but not scored because there is no root stem.

In 100-gene windows of the first species, an observable cross window
has at least five distinct genes on both sides of one direct block.
With both self comparisons present, a pair needs at least
`wgd_min_cross_windows` (default 3) observable windows. If exact copy
support falls short, the same window floor must also be met by direct
anchors connecting genes in the same archived family tree. Otherwise the
pair is `insufficient_archive_anchors`, not observed negative evidence.
An empty self comparison is observed, not missing. The
correspondence search requires five distinct genes in each copy slot,
within the same original blocks. It retains one-copy paths separately,
without inferring a loss, and tests the exact linked genes against the
two daughter copies of a duplication in the input gene tree.

A pair supports a candidate when it has at least `wgd_min_cross_windows`
concordant windows and at least as many distinct matching families.
For combined assessment or support, every daughter lineage must have
an assessed or supporting pair, respectively. With this manifest,
cross-block correspondence determines branch synteny status instead of
the self-synteny/Ks rule; optional genome inputs can still supply the
separate Ks outputs. If correspondence can be assessed but falls short,
the branch stays unresolved. Only branches without sufficient
observability can use the duplication-only rule. The thresholds
are exploratory and are not calibrated WGD probabilities.

With `wgd_cross_orientation_audit=1`, `cross_synteny_orientation.tsv` repeats
the 100-gene window calculation in the second species using the same anchors
and the same family-tree pass. It reports differences in observability and
in per-pair support, with separate window and family counts. It is a
sensitivity audit: the supplied orientation still determines primary calls,
and two orientations do not count as independent evidence.
`cross_synteny_orientation_quartets.tsv` additionally counts exact
copy pairings and families recovered in both orientations, with separate
window counts for each reference. When unique-phase sensitivity is enabled,
it also reports pairings that are unique in both orientations.

Optionally set `wgd_cross_require_both_orientations=1` alongside the unique
phase and orientation options. A pair must meet the minimum family count
and the minimum window count in **each** reference using the same
uniquely concordant copy pairings. Cross-block dS and quartet selection use
only the common set in this mode. Observability and gene-only fallback
are unchanged; a positive signal in one orientation without corroboration
in the other does not become a duplication-only call. The tree additionally
shows `R`, the sum across species pairs of the smaller of their two shared
unique window counts. This is a sensitivity filter, not independent evidence.
The region audit adds a `both_orientations` mode. Regions sharing a window
in either reference are joined; after deleting one region the support rule
still enforces the window floor in both references. The window totals
shown in the region tables remain in the supplied reference orientation.
Existing region and replication modes retain their original scopes.
This option also writes `cross_synteny_anchor_nonreuse.tsv`. It counts
reference windows that contain five gene-disjoint direct pairs within one
original block, both before and after requiring archive-family agreement.
This diagnostic does not require five disjoint four-copy paths or change
the primary definition of an observable window.

`cross_synteny_copy_phase.tsv` separates four-gene sets for the same species
pair and matched duplication node into concordant-only, discordant-only,
and ambiguous paths. Ambiguous means the original blocks admit both ways
of pairing the same four genes across species. Families and windows can
contain multiple quartets and can appear in more than one category. On the
species tree, `W` counts matched paired windows and `U` counts windows with
at least one concordant-only quartet. This audit does not change the primary
candidate decision or turn phase agreement into an age estimate.

For a stricter sensitivity run, set `wgd_cross_require_unique_phase=1`.
Every supporting pair must meet the configured window and family floors
using concordant-only four-copy paths. Ambiguous paths and one-copy paths
do not count toward these floors. Observable pairs that fail the stricter
rule remain assessed, so they cannot trigger a duplication-only fallback.
The pair, window, and orientation audits still show the underlying paths
and baseline support. The downstream cross-block copy-pair dS analysis
also uses only concordant-only paths in this mode. The region holdout tables
add a `unique_phase` mode to report dependence on any one genomic region.

If a focal species has at least two assessed comparisons on the same branch,
`cross_synteny_copy_replication.tsv` shows whether their paired paths use the
same exact focal gene pairs in the same archived families. Strict mode also
reports the overlap using only concordant-only paths. Shared focal genes or
family trees make these checks dependent; they are not independent P-values.

To expand an existing manifest without all-vs-all comparisons, run the sparse
planner in the GeneGalleon runtime:

```bash
python workflow/support/wgd_synteny_plan.py \
  --workspace workspace --species-tree /path/to/species_tree.nwk \
  --cross-synteny /path/to/wgd_cross_synteny.tsv \
  --pair-evidence /path/to/wgd_ssd/cross_synteny_pairs.tsv \
  --pairs-out /path/to/next_pairs.tsv --audit-out /path/to/branch_coverage.tsv
```

The planner visits non-root splits from descendants toward the root, using
species with self-synteny in the input manifest. It proposes new pairs only
for daughter lineages that are not covered, preferring pairs that cover more
missing daughters, then shorter tree distance. With `--pair-evidence`, only
assessed comparisons count toward coverage; unassessed pairs are not repeated.
With `--pair-evidence`, an optional `--replication-budget N` proposes up to
`N` additional near-to-root comparisons on branches with one assessed pair.
It changes only one member of that pair, so the shared species' exact gene
pairs can be compared. It does not replace comparisons for missing coverage
or change any existing candidate calls.
Without it, coverage means only that a pair was supplied, not that it was
observable. The audit distinguishes missing self-synteny from exhausted
cross-pair choices. When the plan is nonempty, `next_pairs.tsv` uses the
`synteny_pairs_file` columns for the separate pairwise-synteny stage. The
planner does not create anchors or modify the existing WGD manifest. After
generating and validating the original anchors, add those comparisons to
the manifest before rerunning WGD. An empty plan leaves just the header.

### Relative Age From Cross-Block dS

With a cross-synteny manifest, optionally set `wgd_cross_ks_cds_file` to a
TSV with columns `species` and `cds` (absolute or workspace-relative FASTA
paths, optionally gzip-compressed). An empty `cds` cell uses the usual
`input/species_cds` discovery. Use the same representative gene IDs in the CDS.
A species with CDS but no BED can supply ortholog
distances for divergence calibration; it does not supply synteny evidence.

For combined candidates, GeneGalleon measures dS for exact within-species
copy pairs from paired cross-synteny windows that match the same duplication
in a gene tree. Partial paths are not aged. `wgd_max_pairs` caps the number
of distinct copy pairs selected per branch and focal species; the selection
spans windows, and an identical request is estimated only once. Ortholog
calibration uses archived members of single-copy families, up to
`wgd_max_ks_families` per comparison. For each nearby boundary, choose
one nearest-outgroup trio by the smallest shared-family count among its
three comparisons (larger is better), before measuring dS. Ties are
broken by input species IDs. It does not run an all-versus-all sequence search.
Cross-code pairs are not estimated. Low-information alignments (fewer than
60 retained codons, less than 70% coverage of the longer CDS, or fewer than
20 synonymous sites) and non-estimable dS values are audited and excluded.

`wgd_cross_ks_method` defaults to `yn00`. Use `gy94` or
`gy94-regularized` for a separate estimator sensitivity run across all
cross-block, ortholog and optional terminal/paralog distances in that run.
The procedure does not substitute YN00 estimates for failed likelihood
fits or treat GY94 parameter boundaries as measured dS. Pairs with
identical coding sequences have zero estimated distance but no fitted
synonymous-site count under GY94; they are recorded as `identical` and
excluded from the quality-filtered summaries. Alignment reuse is independent
of the subsequent estimator. This setting does not change self-synteny
distances outside the cross-block module.

On the species tree, inverted primary boundaries remain visible as dashed
intervals and crosses. When matched-family estimates are available for such a
row, open diamonds and dotted intervals show the alternative split estimates.
These are diagnostics and do not change the primary event-age decision.

Set `wgd_max_ks_outgroup_families` above zero to check the boundary ordering
with a common external species for both splits. Up to two outgroups outside
both splits are chosen by phylogenetic proximity and archived-family
coverage, before looking at dS. The same seeded family sample is measured for
both outgroups, and only families with all required distances contribute.
At the corresponding focal lane on the species tree, a separate signed
axis shows the median of the per-family older-minus-younger corrected dS,
its simultaneous order-statistic interval (when estimable), and the family
count. A negative gap reverses the expected order; an interval crossing zero
does not resolve it. The diagnostic assumes independent family observations
and is not used to call or date a WGD. The summaries and input pair results
are in `cross_ks_common_outgroups.tsv` and `cross_ks_common_outgroup_pairs.tsv`.
The default zero skips the extra alignments.

To rerun with different WGD settings, `wgd_cross_ks_reuse_summary` accepts
one or more newline-separated `summary.json` paths in separate, completed
run directories with their original `cross_ks_aligned_pairs.code*.tsv`
files. The CLI accepts repeated `--cross-ks-reuse-summary` arguments. No source
bundle may be the destination. GeneGalleon checks the source hashes, MAFFT and backtranslation
identity, genetic code, ordered gene IDs, and the exact current CDS for each
reused row. Identical alignments can be shared by bundles; conflicting
verified alignments for the same pair and current CDS stop the run instead
of choosing one by input order. New or changed pairs are aligned normally. All dS estimates and
downstream tests are recomputed. Summaries without alignment provenance or
unsupported MAFFT layouts cannot be used for reuse.

The NWKIT `ksrate` model puts the younger and older species divergence boundaries
on each focal species' dS scale, using conservative simultaneous
pair-median intervals. Copy-pair observations sharing either genomic
window or a family are grouped before taking the median across components.
Its 95% order-statistic interval needs at least six components (or
`wgd_min_cross_windows`, when larger); with fewer, only a point is shown.
This interval assumes independent identically distributed component medians;
grouping reduces dependence but does not establish that assumption. A focal
supports placement only when its whole interval lies between the two divergence intervals. Branch-level
support requires at least one supporting focal in every daughter lineage
and no focal with a clearly older or younger interval. Missing or overlapping
boundaries do not establish placement, but a usable boundary can still flag a
copy interval entirely outside the candidate. Divergent copy ages are flagged,
not used to reject or move a WGD automatically: dS dates copy divergence,
which can precede allopolyploidization. The diagnostic is not an independent
P-value and does not change the original branch or gene-node classification.

Set `wgd_ks_component_holdout=1` to display the minimum and maximum copy
median after excluding any one linked component. This range is shown just
below the focal dS point when the extrema differ, and written to
`cross_ks_component_holdout.tsv` along with the observation and component
counts. For zero or one component there is no holdout range. It measures
influence, not confidence, and does not change the primary interval.

With `wgd_max_ks_quartets` above zero, `cross_ks_quartets.tsv` adds a
within-family check for combined candidates with paired paths and same-code
CDS in both species. The limit is per branch and species pair; each selected
window/family contributes one quartet, sampled before looking at dS. The six
dS values compare the three possible pairs of opposing paths. A positive
`paired_margin` means the two synteny-matched copy paths have the smallest
sum; `other_sums_difference` records deviation from an additive tree.
Missing or saturated pairs yield an incomplete quartet. This is a check
of local copy correspondence, not a clock correction, age estimate, or test
of the branch or number of WGD events. The default zero skips additional
alignments; pairs already requested for the age analysis are reused.

Set `wgd_max_terminal_ks_families` above zero to inspect terminal copy splits
in species with supplied self anchors and CDS. The analysis samples at most
one same-family self-anchor pair per family, in rounds across original anchor
blocks with the configured seed. A candidate quartet requires a binary
focal-only ancestor for that pair, a same-code relative from the immediate
species sister group at the adjacent gene node, and an external copy outside
both groups in both trees. It records which selections fail and whether the
six-pair dS comparison is consistent with the focal copies forming a clade.
Triangles on the terminal rows of the species tree show the median copy
dS across original-block medians, the 95% order-statistic interval when
estimable, and the numbers of measured families and blocks. Blocks need not
be independent, and quartet agreement does not distinguish terminal SSD
from WGD or establish multiplicity. This option adds diagnostics, not calls.
See `terminal_ks_quartets.tsv` and `terminal_ks_summary.tsv`.

Without a gene-tree ZIP, the default is the calibrated gene-count mode.

NWKIT models linear duplication/loss, positive-geometric ancestral counts,
terminal/internal background rate groups, and four fixed-shape gamma family
rate categories. At an event, each ancestral copy contributes its original
copy plus binomially retained extra copies. `wgd_multiplicity=2` tests doubling;
3 is a separate triplication alternative, not an estimated ploidy.

Optionally supply `wgd_count_detection_tsv` with `leaf_name` and
`detection_probability` columns, covering every tip exactly once with a
known probability in `(0, 1]`. It models independently missed copies;
it does not correct extra annotation copies or turn a missing value into
zero. Do not treat BUSCO completeness as a per-copy probability without
independent justification. Supply `wgd_count_rate_groups_tsv` with `branch_id`
and `regime` columns to replace the terminal/internal background groups.
Include every non-root branch of the input tree once. Specify groups and
probabilities before inspecting candidate support; rerunning with a different
table is a sensitivity analysis with its own null fit and bootstrap.

The scan searches every positive-length non-root branch at fractions
0.25, 0.5, and 0.75 by default. Set `wgd_count_event_fractions` to a specific
comma-separated grid measured from each branch's parent, for example `0.5`
for a midpoint-only scan. `wgd_count_family_gamma_shape` changes the fixed
shape, or `none` disables family-rate heterogeneity. Changing either is a
different model/search and requires its own bootstrap and checkpoint.
It fits one event at a time. State doubling must meet a
per-family likelihood tolerance; optimization or convergence failures stop
the stage. A separately fitted branch-specific duplication/loss burst provides
an AIC diagnostic against coordinated genome multiplication.

`wgd_count_bootstrap=199` refits null simulations and calibrates the maximum
statistic over the searched branch/time grid. These are conditional plug-in
bootstrap P-values, not posterior probabilities or guaranteed error control
against arbitrary heterogeneous SSD histories. Zero bootstrap replicates are
an exploratory scan and cannot support a WGD origin. Inspect Monte Carlo
error, nuisance parameter boundaries, and the competing branch-burst fit.
Use `wgd_count_bootstrap_workers` with a NWKIT runtime that supports
`wgd-count --bootstrap-workers` to run more than one null refit at a time.
Use `wgd_count_scan_workers` with `wgd-count --scan-workers` to parallelize
the preceding observed candidate fits. The two phases run in sequence and
each defaults to one worker. Do not exceed the allocated CPUs or memory
for concurrent fits. With a runtime supporting `wgd-count --bootstrap-checkpoint`,
set `wgd_count_bootstrap_checkpoint` to a file in an existing directory outside
curated inputs, task scratch, and managed WGD outputs. Reusing the same path
with unchanged input, settings, and NWKIT runtime refits only unfinished
null draws. An incompatible checkpoint is rejected rather than overwritten.
Set `wgd_alpha` to change the search-wide conditional P-value cutoff
(default `0.05`). Both combined and count-only support require explicit
non-boundary diagnostics for the event,
background, and branch-burst fits. Missing diagnostics or a numerical bound in
any competing fit leave the event unresolved; an unbounded event fit alone
does not establish an adequate comparison.

Each branch records its synteny opportunity separately from its support.
A successful zero-block scan is assessed, not missing. For a nonempty
scan, a species is observable on a branch when it has a usable anchor Ks
estimate and a resolvable branch interval; qualifying positive evidence
also establishes observability. If fewer than two descendant species
(one for a terminal branch) are observable, a calibrated count-supported
branch is labeled `count-only-candidate`. When enough species are
observable but support falls short, the event remains `unresolved`.
Count-only candidates do not assign individual genes a WGD origin.

Only families observed in every root-child clade are retained; the likelihood
and simulations condition on this rule separately for each missingness mask.
This does not correct all OrthoFinder clustering, annotation, or other
selection errors. Detection probabilities and detailed branch-rate regimes
are available in the standalone NWKIT model; BUSCO percentages are not
automatically treated as known per-gene detection probabilities.

Counts, single-copy comparison families, and self-anchor pairs have explicit
uniform-sampling limits (`wgd_max_count_families`, `wgd_max_ks_families`,
`wgd_max_pairs`) and reproducible selection audits. Pairwise MAFFT protein
alignment and CDSKIT backalignment precede YN00 dS. Saturation, invalid sites,
unestimated pairs, and missing comparisons remain distinct. Calibration excludes
families with known multiple copies anywhere in the full membership/count tables;
each compared species must have an observed count of one and one selected CDS.
Unknown counts are not assumed single-copy. These family contrasts can still
contain hidden paralogy; they are not asserted to be gene-tree-validated orthologs.
Sampling lowers the observed branch-matched coverage; it is not extrapolated
to unsampled pairs. A sample of `N` pairs can cover at most `2*N` loci. For
large annotations, increase `wgd_max_pairs` if this upper bound cannot meet
the configured coverage rule; do not interpret a sampling-limited unresolved
result as evidence against WGD.

Unquota self-block arms, depth, and coverage are observed extant annotation
metrics, not ploidy or ancestral loss estimates. The initial combined evidence
rule requires calibrated count support, preference over the branch-burst AIC
alternative, and branch-matched anchor coverage of at least 0.2 in at least
three locus-disjoint block components per species. Blocks sharing anchor loci
are joined into one component so redundant block IDs cannot inflate replication.
Only positioned, different-locus, nonoverlapping, non-tandem anchors with a valid
positive Ks and supported branch placement count. Coverage counts unique anchor
loci over all annotated loci, including loci
without selected sequences; position ranks also include intervening annotated
loci, and coordinate bounds cover the full annotated locus, not just its selected
transcript. Missing annotated-locus totals cannot be replaced by sequence counts.
Internal candidates require two descendant species;
terminal candidates require one. These are configurable experimental evidence
rules, not probability thresholds. Local or redundant segmental blocks alone
are insufficient.

Ks placement must lie strictly between the primary `ci_lower`/`ci_upper`
intervals of neighboring focal-lineage-corrected speciation boundaries.
`ks_boundaries.tsv` records `ci_method`; the evidence summary reports all methods
present. Optional bootstrap diagnostic intervals never replace the primary bounds.
GeneGalleon explicitly requests `pair-median-bonferroni`: simultaneous
noninterpolated order-statistic intervals for all used pair population medians,
propagated through focal/sister/outgroup correction and the median across trios.
This assumes independent identically distributed families within each pair,
not independence between pairs. Small samples can have unbounded intervals and
therefore unresolved placement. Bootstrap percentile intervals remain diagnostics;
standalone NWKIT's default method is unchanged. Negative,
inconsistent, nonmonotone, missing, and boundary-overlapping corrections stay
unresolved. No external outgroup means the root boundary cannot be resolved;
this limits placement on the oldest branches. A pair's dS point estimate has
no substitution-model confidence interval in this initial implementation.

## Gene-Node Classification

After completing genome inference, set `run_wgd_ssd_classification=1` in the
gene-evolution entrypoint. `wgd_evidence_dir` defaults to
`output/genome_evolution/wgd_ssd`. The stage reconciles the chosen rooted gene
tree against the **full** matching species tree, not the family-pruned tree.
Both trees must satisfy the rootedness contract above before the classification
output directory is created or reconciliation runs.
NWKIT reconciliation additionally requires strictly bifurcating trees; a declared
root polytomy accepted for count inference is not arbitrarily resolved for classification.
GeneRax NHX events are respected when present; otherwise NWKIT LCA events are
used. Transfer/speciation events are not relabeled as duplication origins.
Relative `wgd_evidence_dir` paths are resolved from the workspace.

- `WGD-supported`: a cross-child anchor pair matches the reconciled species
  branch, and that event has combined count/synteny/Ks evidence support.
- `SSD-supported`: a fully mapped terminal two-tip duplication has direct
  different-locus, nonoverlapping tandem adjacency on the matching terminal
  species branch and no conflicting collinearity evidence.
- `unresolved`: everything else, including insufficient coordinates, same-locus
  annotation ambiguity, overlapping loci, parser/coordinate species conflicts,
  distant/dispersed copies without origin evidence,
  event/age mismatch, or tandem/collinearity conflict.

Proximal/distant/interchromosomal features are reported, but none is converted
to SSD simply because a synteny anchor is absent. Labels are conditional on one
input gene tree; input support is retained, but these labels are not bootstrap
stability estimates or posterior probabilities. They do not identify auto-
versus allopolyploid origin, donor lineages, or rediploidization histories.
The separate [native MUL-tree search](grampa-replacement.md) ranks
auto-/allopolyploid topology hypotheses by D+L parsimony; it does not add a
donor-network likelihood or posterior to these WGD-origin assignments.
Tandem adjacency supports an SSD hypothesis but cannot exclude WGD-derived copies
brought together by later rearrangement. Missing collinearity is not proof that
copies were never WGD-derived. Broad segmental duplication can also resemble WGD
evidence; this experimental combined rule is not a validated discriminator for
all such histories.

### Optional MUL Node Diagnostics

Set `wgd_mul_diagnostics=1` alongside `run_wgd_ssd_classification=1`, with an
explicit `wgd_mul_h1` selector in the **full, original-label** species tree.
`wgd_mul_h2` optionally restricts second parents; empty searches all admissible
placements. This runs a separate per-family `nwkit mul-reconcile --score-model
dl` search on the same unmodified rooted gene/species inputs. It does not join
the genome-stage GRAMPA bundle, change WGD/SSD rules, or infer DL+ILS node origins.

Every globally tied best candidate and every optimal assignment is retained in
`mul/nodes.tsv`, with `mul/scores.tsv` and `mul/results.json`. Rooted topology
IDs and descendant-tip clade IDs verify the join independently of child order
and branch lengths. Leaf species, candidate metadata, assignment completeness,
score totals and uniqueness must agree; mismatches fail the stage. Original
input/evidence hashes remain checked, and the entire family ZIP publishes
atomically, preserving previous results on failure.

`node_diagnostics.tsv` has one row per internal gene node, including
non-duplication events with classification `NA`. It retains the original
reconciliation event, classification/reason, and adds `mul_mapping_status`
(`consistent` or `ambiguous`), `mul_dl_duplication` (`all`, `none`, `ambiguous`),
candidate/mapping counts and all distinct mapped occurrence clades/cost states.
Consistency means that the mapped tip occurrences, species and local D+L costs
agree across all co-optimal assignments, not biological certainty or stability
under tree/model uncertainty. Counts are enumerations, not posterior weights.
Even `mul_dl_duplication=none` does not overwrite `SSD-supported` evidence.

The classified NHX additionally stores `mul_mapping_status`,
`mul_dl_duplication`, `mul_gene_node`, `mul_best_hypotheses` and
`mul_optimal_mappings`. `duplication_origins_mul.pdf` uses the original origin
colors and a circle for consistent MUL mappings or a triangle for ambiguous
ones. Its `N<number>` labels join to the run-local `gene_node`; clade IDs are
the child-order-independent identifiers. Existing origin TSV columns and the
original PDF are unchanged. Stage-owned MUL attributes are cleared on rerun,
including when diagnostics are disabled.

The default is off. `wgd_mul_max_candidates`, `wgd_mul_max_state_pairs` and
`wgd_mul_max_maps` default to NWKIT's 10,000 candidates, 10,000,000 state pairs
per gene/candidate and 100,000 optimal assignments per gene/candidate.
Exceeding a cap fails without filtering candidates or truncating assignments.
Settings support the normal `GG_GENE_EVOLUTION_` overrides and participate in
artifact invalidation. Disabling the option rebuilds the family bundle under
the existing stale-artifact policy; cached diagnostics are not treated as new
disabled-mode results.

Set `wgd_native_tree_likelihood=1` to additionally run `nwkit wgd-tree` using
the count-fit parameters and each single-event candidate. This optional model
uses a rooted binary species-colored gene topology, not observed gene branch
lengths or sequence likelihood. It integrates duplication/loss, sampling,
the positive-geometric root prior, family-rate categories and WGD retention.
Branch likelihood flow is an analytic positive polynomial, evaluated in log
space with double/extended-arithmetic agreement checks. These numerical
differences are diagnostics, not certified roundoff bounds.
Unsupported transfers, unmatched tips, nonbinary trees, triplication candidates, excessive topology
size, or numerical failures stop this optional stage explicitly.

`native_conditional_wgd_probability` is the latent origin probability for the
reconciled event/node pair **given the fixed topology, supplied parameters,
root prior, sampling and one candidate doubling**. It is not the probability
that WGD occurred, a parameter/tree posterior, or an independent significance
test. It does not change the evidence-rule classification. The complete
per-candidate table is retained even when its candidate differs from the
LCA-mapped branch. Multiple-event histories are not jointly inferred.

## Outputs and Validation

Genome results are in `output/genome_evolution/wgd_ssd/`:

- `count_candidates.tsv`, `count_model.json`, `counts.tsv`, and
  `count_family_selection.tsv`: model fits, conditional calibration, and inputs.
- `self_synteny/<species>/`: unquota raw/lifted anchors, block/depth tables,
  summary, selected-ID mapping, and tool logs. Valid zero-block outputs are
  distinct from failed alignment/scan commands.
- `pair_selection.tsv`, `ks_family_selection.tsv`, `aligned_pairs.code*.tsv`, `pair_ds.code*.tsv`, and
  `family_ks_contrasts.tsv`: sequence-selection and dS audit.
- `ks_boundaries.tsv`, `ks_trios.tsv`, `ks_model.json` when contrasts exist.
- `anchor_evidence.tsv`, `gene_positions.tsv`, `wgd_events.tsv`,
  `wgd_candidates.pdf`, and `summary.json`: integrated evidence and provenance.
- `wgd_species_tree.pdf` and `wgd_species_tree.png`: the full input species
  tree, with native branch IDs and scores or P-values for supported candidates.
  Teal branches indicate combined support; with cross synteny, dash-dot teal
  indicates a candidate whose support requires one-copy paths. `W` counts
  matched paired windows (not independent events or gene pairs). Dashed orange
  indicates duplication-only or count-only candidates. With a cross
  manifest, dotted purple marks indicate unobserved synteny and a
  duplication score below the cutoff. Labels use only the supplied
  tree and inference output. Gray does not mean WGD was ruled out.
- `cross_synteny_pairs.tsv` and `cross_synteny_windows.tsv` (with a cross
  manifest): per-pair coverage and exact matched-family evidence.
- `cross_synteny_archive_coverage.tsv`: counts direct anchor pairs and
  observable windows before and after requiring both ends of each link to
  be in the same archived gene tree. This diagnostic does not change the
  candidate rule or require any particular duplication topology.
- `cross_synteny_region_summary.tsv` and `cross_synteny_regions.tsv`: merge
  tree-matched windows that share a source block or a window on either copy.
  Report paired-only and all paths separately, with the effect of omitting
  each region on the synteny rule. These are sensitivity diagnostics; a
  region is not necessarily biologically independent, and the original
  candidate decision is unchanged.
- `cross_synteny_nested_copies.tsv`: compare paired hits on nested branches
  that pass the duplication-score cutoff. In species with paired paths on
  both branches, report overlaps of unordered copy-window pairs, exact
  gene pairs, and families separately. Overlap does not by itself establish
  a common duplication; missing comparisons do not imply independence.
- Optionally set `wgd_cross_tree_support_file` to GeneGalleon's GeneRax
  UFBoot ZIP (members named `*_generax.ufboot.nwk`) to write
  `cross_synteny_tree_support.tsv`. This compares the duplication and both
  copy splits for each paired hit by exact unrooted bipartition on the same
  gene set. Missing or unmatched values stay empty; no threshold is applied
  to the candidate decision. This does not test the root position.
- With cross-block CDS inputs, `cross_ks_paralogs.tsv`,
  `cross_ks_orthologs.tsv`, and `cross_ks_focals.tsv` audit genes, selection,
  and focal placement. `cross_ks_boundaries.tsv` and `cross_ks_trios.tsv`
  report corrected divergence when contrasts exist. The species-tree plot
  marks age support and conflicts without changing candidate colors.
  Inconclusive dS estimates have no branch marker. Dotted guides connect each annotated branch to its focal
  species' dS tracks in the same PNG and PDF. Dots show copy-pair medians;
  colored ticks show rate-corrected species divergences, and horizontal
  segments show their available intervals. Each focal label gives the number
  of copy pairs with usable dS (`n`) and the number of grouped components (`c`).
  Nonmonotone split estimates remain visible with crosses and dashed
  intervals, but are not treated as usable branch boundaries.
  The dS tracks use a separate axis
  and do not change the tree's branch units. `wgd_cross_ks.pdf` and
  `wgd_cross_ks.png` retain the standalone age diagnostics.
- `cross_ks_matched_boundaries.tsv` re-estimates the gap between the two
  boundary values on the same complete ortholog families, resampling each
  family jointly. It diagnoses boundary ordering, not WGD age or a rate
  difference between retained duplicates and single-copy genes. It does not
  replace the primary boundary estimates or the branch call.
- `cross_ks_matched_trios.tsv` reports the median of individual signed
  three-distance corrections in complete families, with simultaneous
  order-statistic intervals and counts of negative corrections. It targets
  the complete-family distribution, which may differ from the primary
  pair-median estimate if missingness or substitution rates vary by family.
  It does not replace the primary intervals or the figure's placement rule.

Each family has an additional `wgd_ssd/<family>_wgd_ssd.zip` archive containing
`duplication_origins.tsv`, `pair_evidence.tsv`, full-tree `reconciliation.tsv`,
`classified_gene_tree.nhx`, `duplication_origins.pdf`, and `summary.json`
under its `<family>/` directory. This file participates in normal gene-family
archive collection; the inner table schemas are the same in raw and ZIP storage.
Older uncollected `wgd_ssd/<family>/` directories are preserved. Rebuild their
stale classification stage to produce the managed family archive.
Input NHX annotations are retained, except that this stage's own origin and
conditional-probability properties are cleared and recomputed. The family PDF is topology-spaced, not a
branch-length or age plot. The candidate PDF displays up to 40 highest count
statistics; all candidates remain in the TSV.
Optional native results are `native_tree_origins.tsv` and
`native_tree_likelihood.json`. Evidence hashes reject edited/stale tables before
classification, and hashes are checked again before successful completion.
`gene_positions.tsv` coordinates are zero-based, half-open full-locus bounds.
Existing table columns are retained; `synteny_block_definition` records the
component-counting rule. Normal artifact
stop/rebuild/reuse policy applies.

Numerical transition/pruning references, seeded count simulations, negative
classification cases, and Docker integration tests provide separate kinds of
evidence. Full empirical false-positive-rate benchmarking and real-genome
validation are not supplied by a passing integration fixture. Fixed-topology
likelihood has independent small-case/count-likelihood checks, but is not
claimed equivalent to Whale. Tree ensembles, RR-family reference analysis,
and donor-network likelihood are not implemented by this stage.

Standard container builds verify native WGD dependency imports and commands,
complete annotation locus ordering, and raw self-block coverage exports. Missing
dependency capabilities fail the build before a workflow is started.
