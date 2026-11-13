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

### Sparse Cross-Species Synteny

In duplication-score mode, set `wgd_cross_synteny_file` to a TSV with
columns `species_a`, `species_b`, `bed_a`, `bed_b`, and `anchors`.
Paths are absolute or workspace-relative. Provide each genome's self
comparison (same species on both sides) and one or more direct cross
comparisons spanning the daughter lineages at each internal split.
Species and BED gene IDs must match the input species and gene trees.
Use one BED per species and original `###`-delimited JCVI anchor files
(plain or gzip), not lifted anchors or preassigned branch scores. This
stage reads existing comparisons; it does not launch genome alignments.
The input tree assigns each cross pair to its MRCA; a root comparison
is validated but not scored because there is no root stem.

In 100-gene windows of the first species, an observable cross window
has at least five distinct genes on both sides of one direct block.
With both self comparisons present, a pair is assessed when it has
at least `wgd_min_cross_windows` (default 3) observable windows.
An empty self comparison counts as assessed, not missing. The
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

### Relative Age From Cross-Block dS

With a cross-synteny manifest, optionally set `wgd_cross_ks_cds_file` to a
TSV with columns `species` and `cds` (absolute or workspace-relative FASTA
paths, optionally gzip-compressed). An empty `cds` cell uses the usual
`input/species_cds` discovery. In this mode archive members must name a
family from the membership table and contain only genes in that family
(subsets are accepted). Use the same representative gene IDs in the CDS.
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

The NWKIT `ksrate` model puts the younger and older species divergence boundaries
on each focal species' dS scale, using conservative simultaneous
pair-median intervals. Copy-pair observations sharing either genomic
window or a family are grouped before bootstrapping their component medians.
At least three components (or `wgd_min_cross_windows`, when larger) and
40 bootstrap draws are required. A focal supports placement only when
its whole interval lies between the two divergence intervals. Branch-level
support requires at least one supporting focal in every daughter lineage
and no focal with a clearly older or younger interval. Missing or overlapping
boundaries do not count against a candidate. Divergent copy ages are flagged,
not used to reject or move a WGD automatically: dS dates copy divergence,
which can precede allopolyploidization. The diagnostic is not an independent
P-value and does not change the original branch or gene-node classification.

Without a gene-tree ZIP, the default is the calibrated gene-count mode.

NWKIT models linear duplication/loss, positive-geometric ancestral counts,
terminal/internal background rate groups, and four fixed-shape gamma family
rate categories. At an event, each ancestral copy contributes its original
copy plus binomially retained extra copies. `wgd_multiplicity=2` tests doubling;
3 is a separate triplication alternative, not an estimated ploidy.

The scan searches every positive-length non-root branch at fractions
0.25, 0.5, and 0.75. It fits one event at a time. State doubling must meet a
per-family likelihood tolerance; optimization or convergence failures stop
the stage. A separately fitted branch-specific duplication/loss burst provides
an AIC diagnostic against coordinated genome multiplication.

`wgd_count_bootstrap=199` refits null simulations and calibrates the maximum
statistic over the searched branch/time grid. These are conditional plug-in
bootstrap P-values, not posterior probabilities or guaranteed error control
against arbitrary heterogeneous SSD histories. Zero bootstrap replicates are
an exploratory scan and cannot support a WGD origin. Inspect Monte Carlo
error, nuisance parameter boundaries, and the competing branch-burst fit.
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
  Solid teal branches indicate combined support; dashed orange branches
  indicate duplication-only or count-only candidates. With a cross
  manifest, dotted purple marks indicate unobserved synteny and a
  duplication score below the cutoff. Labels use only the supplied
  tree and inference output. Gray does not mean WGD was ruled out.
- `cross_synteny_pairs.tsv` and `cross_synteny_windows.tsv` (with a cross
  manifest): per-pair coverage and exact matched-family evidence.
- With cross-block CDS inputs, `cross_ks_paralogs.tsv`,
  `cross_ks_orthologs.tsv`, and `cross_ks_focals.tsv` audit genes, selection,
  and focal placement. `cross_ks_boundaries.tsv` and `cross_ks_trios.tsv`
  report corrected divergence when contrasts exist. The species-tree plot
  marks age support and conflicts without changing candidate colors.
  Inconclusive dS estimates have no branch marker. Dotted guides connect each annotated branch to its focal
  species' dS tracks in the same PNG and PDF. Dots show copy-pair medians;
  colored ticks show rate-corrected species divergences, and horizontal
  segments show their available intervals. Each focal label gives the number
  of copy pairs with usable dS (`n`), not the number of bootstrap components.
  Nonmonotone split estimates remain visible with crosses and dashed
  intervals, but are not treated as usable branch boundaries.
  The dS tracks use a separate axis
  and do not change the tree's branch units. `wgd_cross_ks.pdf` and
  `wgd_cross_ks.png` retain the standalone age diagnostics.

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
