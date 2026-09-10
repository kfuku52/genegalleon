# Species-tree calibrations and sensitivity analysis

Existing projects keep their automatic/manual calibration settings and dated-tree
paths. Updating the workflow does not enable extra MCMC runs or silently change
their calibration policy. The workflow adds an inventory of observed calibration
labels; it does not retrospectively certify their source evidence or convergence.

## Interpretation

TimeTree's intervals describe variation across studies. Its
[FAQ](https://timetree.temple.edu/faqs) distinguishes an empirical-rule interval
for at least five studies from a min–max range for fewer studies. The precise
meaning of an API record must be checked against its source studies and metadata.
Do not label every `precomputed_ci_low/high` pair a posterior 95% credible interval.

NWKIT converts the selected bounds into MCMCtree `B` labels. The current automatic
workflow uses soft bounds, with the NWKIT default tail probabilities of 0.025 on
each side, and excludes clades smaller than 20% of the tree. These operations do
not establish independence between studies, fossils, sequences or calibrated nodes.
The interval endpoints do not by themselves justify the shape or tail probabilities
of a calibration density. Secondary calibrations remain supported.

Distinguish the specified calibration density, the effective joint/marginal prior
after tree-ordering constraints, and the posterior conditional on the sequence
likelihood and model. The [PAML documentation](https://github.com/abacus-gene/paml/wiki/MCMCtree)
recommends sampling without sequence data before interpreting posterior dates.
Agreement with TimeTree after using TimeTree as a prior is not independent validation.

## Existing projects

The defaults remain:

```bash
timetree_constraint=1
mcmctree_divergence_time_constraints_str=""
mcmctree_calibration_manifest=""
run_mcmctree_calibration_diagnostics=0
```

The empty manifest preserves the existing artifact contract, automatic API call,
manual-bound alternative, public age units and output names. No new required
input is added to legacy projects. Normal artifact freshness rules still apply.

Whenever `constrained_tree/constrained.nwk` is present, an additive inventory is
written under `constrained_tree/calibration_audit/<tree-sha256>/`. It includes the
observed tree, `candidates.tsv` and `status.json`. Snapshots are immutable; existing
ones are reused. `calibration_audit/current_status.json` identifies the current
observed tree and its matching review record, if one exists. A review record for
a different topology or set of calibration distributions is not applied to it.

Every inventory TSV is an **unreviewed template**, including when the observed
tree was originally generated from a reviewed manifest. Labels alone cannot
recover the TimeTree query date, API response, study count, underlying studies,
shared fossils or shared sequence data. Unknown evidence stays unknown. This
inventory does not refetch TimeTree or implement a second TimeTree parser.
Automatic source-response provenance belongs in NWKIT's retrieval interface;
until available, record those details from the primary sources during review.

## Select reviewed calibrations

1. Copy `candidates.tsv` into a project input file. Never edit the snapshot.
2. Preserve `topology_sha256`, `node_id` and `descendant_tips`. The latter is a
   sorted JSON list of exact tip labels. IDs are based on descendant sets, not
   traversal order or potentially unstable node numbers. A changed rooted
   topology requires a new inventory. Rows for currently uncalibrated internal
   nodes are also included so an independently justified calibration can be added.
3. Set `decision` to `accept`, `exclude` or `unreviewed`. Only `accept` rows are used.
4. Complete the evidence fields below and explicitly specify the distribution.

| Field | Accepted values or required evidence |
|---|---|
| `calibration` | Explicit `B(lower,upper,pL,pU)`, `L(lower,offset,scale,pL)` or `U(upper,pU)`; nonzero bounded interval and valid tails |
| `time_unit` | `Ma`; public age units are preserved |
| `source_type` | `primary`, `secondary`, `timetree` |
| `source_references` | Study DOI/ID and primary evidence references; add query date, study count and response identity when available |
| `interval_kind` | `empirical_rule`, `min_max`, `study_posterior`, `fossil_bounds`, `user_defined` |
| `dependency_groups` | Nonempty JSON array, e.g. `["study-A","fossil-X"]`; membership may overlap |
| `sequence_overlap` | `none`, `shared`, `partial`, `unknown`, including overlap with the target analysis |
| `calibration_overlap` | Same categories, including shared primary calibrations and secondary-calibration ancestry |
| `node_basis` | Target MRCA/crown/stem assignment, taxon substitution, fossil placement and supporting references |
| `distribution_basis` | Why these endpoints, distribution shape and tails are appropriate; a citation to a CI alone is insufficient |
| `reviewer`, `review_note` | Reviewer and explanation of the decision, unresolved dependencies and limitations |

An unknown overlap can be acknowledged with a justified review note. It is not
treated as evidence of independence. Use a common uncertainty group when several
calibrations may share an unresolved source. A different DOI does not establish
independence; overlapping taxa alone do not establish duplicated sequence data.

Point calibrations and other distribution families are intentionally unavailable
in this new reviewed interface. The legacy interface remains unchanged. New
families require coordinated NWKIT parsing, GeneGalleon unit conversion, examples
and runtime tests, rather than silent approximation.

For example, set these variables in the genome-evolution entrypoint config block:

```bash
mcmctree_calibration_manifest="/absolute/path/to/project/input/calibrations.reviewed.tsv"
artifact_stale_policy="stop"
```

The manifest overrides automatic/manual selection, without requiring a change to
`timetree_constraint`. A malformed manifest is rejected before replacing the tree.
`reviewed_calibrations.json` stores the manifest digest and review evidence.
Its status is `reviewed_inputs`, **not** scientific validation or MCMC convergence.

The default stale policy stops when an existing analysis needs rebuilding. To
perform that rebuild explicitly, use `artifact_stale_policy="rebuild"`. Before
replacing an existing constrained tree, the workflow copies the previous
`constrained_tree`, `mcmctree_parameter_estimation`, `mcmctree_main` directories and
species-tree provenance records into `species_tree/calibration_history/`.
Account for the extra storage, particularly Hessians and retained MCMC files.
Dependent outputs follow their existing input-digest freshness checks.
`artifact_stale_policy="reuse"` is rejected with reviewed input or enabled
diagnostics because it could combine a new requested calibration with stale dates.

## Run additional diagnostics

```bash
run_mcmctree_calibration_diagnostics=1
mcmctree_calibration_diagnostic_chains=4
mcmctree_calibration_diagnostic_seed=1729
```

The existing `mcmc_burnin`, `mcmc_sampfreq`, `mcmc_nsample` settings, clock and
birth–death/rate priors flow through the IQ2MC control file. Runs are sequential
to avoid multiplying memory demand. Each target and chain has a distinct recorded
positive seed. These experiments run in a separate directory and never replace
the main dated tree or promote a favorable sensitivity result into the main analysis.

The scenarios are the supplied baseline, each calibration omitted in turn, and
each evidence group omitted when a matching reviewed manifest is supplied.
Each scenario runs `usedata=0` and `usedata=2` with the same remaining priors.
Omitting the last calibration is recorded as `skipped_no_remaining_calibration`:
the generic IQ2MC `RootAge` must not silently become the sole absolute-age anchor.
The retained `RootAge` control is recorded for all other scenarios and must also
be considered when interpreting root-calibration removal.

Results are under `species_tree/mcmctree_calibration_diagnostics/<run-sha256>/`:

- Frozen inputs and a contract describing inputs, adapters, executable and versions.
- Every chain's control, tree, raw `mcmc.txt`, standard output/error and timing.
- `diagnostic_samples.tsv`: regularly spaced retained draws. PAML also prints
  iteration 1 outside that grid when `sampfreq > 1`; it stays in the raw file but
  is not included in autocorrelation diagnostics. Incomplete generation grids fail.
- `summary.tsv`: per-node/parameter mean, median, **95% equal-tail** interval,
  rank-normalized/folded split R-hat, bulk/tail ESS, MCSE and numerical status.
- Per-target node correlations and mass outside the specified calibration bounds.
- `comparisons.tsv`: posterior-versus-prior and posterior-versus-baseline median
  differences in Ma and interval-width ratios. These are descriptive sensitivities.
- `status.json`: run status, skipped scenarios, source limitations and output hashes.

Numerical screening uses R-hat ≤1.01, bulk/tail ESS ≥400, and MCSE(mean)/SD ≤0.05
for every sampled parameter except the known constant prior likelihood, if printed.
Nonfinite diagnostics or constant age/rate chains do not pass. Diagnostics use
the R [posterior package](https://mc-stan.org/posterior/reference/rhat.html), included
in updated GeneGalleon containers. `diagnostics_pass` is a numerical screen only;
it does not certify calibration evidence, model fit or scientific robustness.
`not_converged` preserves all outputs and makes their status explicit. Execution
errors record `failed` and stop the requested diagnostic stage. No failure is
silently replaced by a smaller calibration set or point estimate.

Unchanged completed experiments are reused only after verifying their output
hashes. Failed/incomplete runs are preserved and not overwritten; inspect them
and use a new output directory or an explicitly changed seed/control for a new
run. Main-analysis samples are not retroactively recovered by this feature.

The standalone runner also accepts separately justified alternative inputs:

```bash
python workflow/support/mcmctree_calibration_experiments.py \
  --control /absolute/path/to/normalized-public-unit.ctl \
  --tree /absolute/path/to/iq2mc.rooted.nwk \
  --alignment /absolute/path/to/iq2mc.dummy.phy \
  --hessian /absolute/path/to/iq2mc.mcmctree.hessian \
  --manifest /absolute/path/to/calibrations.reviewed.tsv \
  --outdir /absolute/path/to/independent-experiments \
  --time-scale 100 --chains 4 --seed 1729
```

Here `--time-scale` is the divisor from public Ma to the internal unit of the
control's rate priors and the original IQ2MC run; `100` is only an example.
The control's RootAge and supplied calibration ages must be in public units, as
in GeneGalleon's exported IQ2MC artifacts. Never change the divisor independently
of those priors. Reuse a Hessian only with the exact same fixed topology, tip order,
alignment partitions and fitted substitution model. Regenerate IQ2MC inputs when
those change. Group exclusion does not change topology or sequence data.

## Scientific acceptance and cost

Compare effective priors with the specified distributions before interpreting
posteriors. Node correlations in prior-only output reflect the implemented prior,
not a reconstruction of unmodeled correlations between source studies. Do not
choose calibrations because they yield the preferred posterior dates.

Predefine focal nodes and a scientifically meaningful age tolerance. Compare
direct-evidence, audited secondary, source-group exclusion and defensible
alternative distribution scenarios. A lack of significant differences does not
prove robustness. Interval-width ratios have no automatic biological pass threshold.

The runtime tests exercise real IQ2MC/MCMCtree, scale conversion, source preservation,
cache integrity and diagnostic failure detection. They do not establish general
95% coverage. A separate simulation study must generate source studies as well as
target data, vary shared sequences/fossils and node dependence, and measure truth
bias, RMSE, interval width, coverage and failures across independent datasets.
Start with 20–50 pilot replicates per condition; about 456 independent replicates
give roughly ±2 percentage points at 95% coverage, before allowance for design effects.

With K calibrations and G evidence groups there are at most `1 + K + G` scenarios,
minus explicitly skipped ones. At C chains, each scenario adds `2*C` MCMC runs.
Use pilot elapsed time and ESS/second rather than extrapolating only from tip count.
Each chain uses approximately `burnin + sampfreq*nsample` iterations. Frozen Hessian
and alignment inputs are copied once per experiment bundle, not once per chain.
Propagating accepted alternative dated trees to RADTE, CAFE, OU and PGLS remains a
separate downstream analysis; this feature does not propagate their uncertainty.

## Implementation validation (2026-09-10)

Validation used the current checkout mounted into the local GeneGalleon Docker
image `local/genegalleon:calibration-review-dev`, based on
`local/genegalleon:standard-iqtree-dev` with R `posterior` installed. This is Docker
overlay validation, not a full rebuild of the production image or a SIF result.

Overlapping focused runs passed:

- 264 tests: shell safety, existing MCMCtree time scaling and calibration runtime tests.
- 128 tests: calibration runtime, entrypoint configuration, artifact provenance
  and existing genome-evolution protein-mode workflows.
- 34 tests: updated calibration runtime tests (18 cases, including real
  IQ2MC/MCMCtree in two public/internal unit scales) and container source policy.
- All 27 R support scripts parsed; changed Python files passed Ruff and changed
  shell files passed Bash syntax checks.

The genome-evolution entrypoint dry-run passed. The all-entrypoint dry-run could
not complete its orthogroup gene-evolution step because this checkout lacks
`Orthogroups.GeneCount.selected.tsv`; no input was fabricated to bypass that
prerequisite. No complete real-project analysis, general interval-coverage study,
production image rebuild or SIF execution was performed.
