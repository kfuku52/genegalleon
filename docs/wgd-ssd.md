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
- Optional `input/wgd_genomes.tsv`. Without this file, only count candidates
  are reported, and origins remain unresolved.

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
Combined WGD support requires explicit non-boundary diagnostics for the event,
background, and branch-burst fits. Missing diagnostics or a numerical bound in
any competing fit leave the event unresolved; an unbounded event fit alone
does not establish an adequate comparison.

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
