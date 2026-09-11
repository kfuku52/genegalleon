# NWKIT OU shift analysis

`run_native_ou=1` in `gg_gene_evolution_entrypoint.sh` enables NWKIT's native
multivariate OU shift analysis. OU analysis remains disabled unless requested.
NWKIT replaces the former kfl1ou workflow and container dependency. The default
selection is **AICc with convergence enabled**, with an automatic shift cap.
This workflow choice does not establish statistical superiority: AICc selects a
predictive model and does not control false shift detections. See the
[validation record](native-ou-validation.md) for measured accuracy and limits.

The model shares shift locations across traits. The default uses separate OU
covariance parameters per trait; `native_ou_trait_covariance="full"` estimates
evolutionary covariance across traits, shared across regimes. It supports
sampling errors, additional estimated observation variance, missing coordinates
and fixed/stationary roots. Trees must be rooted, binary,
ultrametric and have positive branches; invalid trees are rejected.

## Configuration

```bash
run_native_ou=1
native_ou_criterion="AICc"
native_ou_max_shifts="auto"
native_ou_convergence=1
native_ou_calibration_replicates=199
native_ou_calibration_level="0.05"
native_ou_bootstrap=0
native_ou_seed=1
native_ou_bootstrap_seed=2
native_ou_root_model="OUfixedRoot"
native_ou_trait_covariance="diagonal"
native_ou_alpha_model="trait-specific"
native_ou_estimate_measurement_error="yes"
native_ou_search_strategy="auto"
native_ou_replicate_separator="_"
treevis_branch_color="ou_native_regime"
```

`auto` searches the full location/regime space when it fits NWKIT's
5,000-candidate exhaustive budget, otherwise using beam/local search with
shared regimes and nested returns. `native_ou_max_shifts="auto"` starts at the
structural limit of N−2 for N tips. For beam search it is capped further by
NWKIT's candidate pool and refit budget minus one. GeneGalleon omits candidate,
refit, screening and beam-width options so the installed NWKIT version supplies
their defaults. These budgets are recorded in the NWKIT model configuration.
An explicit `native_ou_max_shifts` integer still sets the requested shift cap. AICc chooses the final
model among evaluated candidates; the automatic cap is not an estimate of the
true number of shifts. The model JSON records the requested and resolved limits
and whether a computational budget reduced the cap.

`native-path` remains available with `native_ou_convergence=0` and criterion
`AIC` or `AICc`. It generates joint configurations along group-lasso paths,
refits retained layouts without shrinkage, and selects the lowest criterion.
Its automatic cap uses N−2 and the refit budget, without the branch-pool limit.
It does not search shared regimes. Both heuristic strategies are approximate;
finite budgets and optimizer iterations can limit coverage.

The convergence-capable strategies also support `BIC`, `pBIC`, and `bootstrap`.
`pBIC` is not available for full covariance or shared-alpha joint fits.
The last uses the experimental sequential plug-in bootstrap with
`native_ou_calibration_replicates` draws per test; it is not a proven uniform
error guarantee. Positive `native_ou_bootstrap` repeats the entire selected
procedure for stability frequencies, including inner calibration when
`native_ou_criterion="bootstrap"`.

### Evolutionary covariance across traits

Opt in with a current NWKIT runtime containing the joint-covariance implementation:

```bash
native_ou_trait_covariance="full"
native_ou_alpha_model="trait-specific"
native_ou_criterion="bootstrap"
```

`shared` estimates one alpha; `trait-specific` estimates one per trait. Full
covariance adds p(p−1)/2 cross-trait parameters. The model JSON exports
`joint_covariance.process_tip_covariance` and `diffusion_covariance` in original
trait/time units. Joint likelihood is reported once, not as additive per-trait
likelihoods. These choices are included in artifact fingerprints and resume checks.

Keep `native_ou_alpha_model="trait-specific"` as the default: it allows traits
to have different evolutionary time scales. Use `native_ou_alpha_model="shared"`
when a common time scale is scientifically justified, or as a sensitivity
analysis. Alpha sharing and evolutionary covariance are independent choices;
full covariance does not require shared alpha. For a sensitivity comparison,
keep the data, covariance structure, observation-error treatment, root model,
selection criterion and search settings the same, and use separate output
directories. Compare selected shifts and their stability, not just runtime.
The small error-and-missingness benchmark compared AIC/BIC searches; it does
not establish superiority under GeneGalleon's default AICc criterion.

Complete data without known or estimated observation error permit an exact,
fast covariance profile for shared alpha. GeneGalleon's default replicate
sampling errors and additional observation-error estimation use the general
joint fit. Current NWKIT selects dense GLS with analytic gradients for small
joint inputs and tree-based computation for larger inputs, supporting both
alpha models. Do not remove those errors merely to enable the separable
profile. Full covariance does not automatically estimate cross-trait sampling
errors; replicate-mean variances remain diagonal plug-in estimates.

The full-covariance AICc option uses observed tip vectors as its sample-size
convention and is a heuristic score, not a calibrated detection test. Bootstrap
regenerates correlated data and repeats candidate screening and covariance
fitting. Its finite-sample calibration still requires validation. See the
[100-tip experiment and implementation proposal](benchmarks/ou-trait-covariance-100tips/README.md);
the original simplified-model timing does not predict this workflow's runtime.

NWKIT also provides `shift-simulate` to generate trait TSVs and known shift/covariance
truth from explicit parameters or a completed native model. This is an
unconditional simulation, distinct from posterior ancestral-state draws.

Former `native_ou_candidate_pool`, `native_ou_refit_budget`,
`native_ou_screening_budget` and `native_ou_beam_width` entrypoint settings are
removed. For a standalone reproduction, the Python adapter still accepts explicit
`--candidate-pool`, `--refit-budget`, `--screening-budget` and `--beam-width`
overrides; omitted options inherit NWKIT defaults. Computational budgets are
not wall-clock deadlines: larger trees, more traits and resampling add cost.
See the [budget measurements](benchmarks/native-ou-default-budgets/README.md)
for the paired 1,000-tip runs used to tune NWKIT defaults.

## Migration from kfl1ou

Replace `run_l1ou=1` with `run_native_ou=1` and remove the old `l1ou_*`,
`large_tree_num_gene`, and `large_tree_max_nshift` settings. Use
`native_ou_max_shifts="auto"` or an explicit integer cap instead. Old RData fits
and `l1ou_*` tables are not reused. New results use `ou_native/<family>_ou_native.*`; downstream summaries and
tree coloring read the new model. Set `treevis_branch_color="ou_native_regime"`.
Old output files can remain on disk but are no longer read by this workflow.
The old derived optimum-expression/tau columns are not synthesized from
unidentifiable OU optima; use the exported effects and regime-parameter tables.

Expression columns such as `root_1,root_2,leaf_1,leaf_2` become `root,leaf`.
Only the final separator-delimited suffix is removed; use an empty separator
to treat every input column as a distinct trait. Means ignore missing replicates;
the sampling variance of a mean is unbiased sample variance divided by the
observed replicate count. No observations produce a missing coordinate. One
observation has no estimable sampling variance: the supplied known component is
zero and the audit explicitly marks it unavailable, while the default additional
trait-level observation variance is estimated. These sampling variances are
plug-in estimates treated as fixed during calibration, not known population
variances. Replicate dependence or batch effects need a separate model.

## Outputs and reruns

The stage publishes a transactional bundle under `ou_native/<family>_ou_native.*`:

| Suffix | Contents |
| --- | --- |
| `model.json` | NWKIT parameters, likelihoods, calibration, search coverage, support, input/source fingerprints and replicate audit |
| `regime-map.tsv` | Complete NWKIT branch-ID/regime map, including root |
| `effects.tsv` | Trait/shift mean effects; optimum effects only where supported |
| `regimes.tsv` | Trait/regime parameters; unsupported optima are NA |
| `tips.tsv` | Observations, fitted means, residuals and sampling SE by tip/trait |
| `replicates.tsv` | Observed replicate counts, means and sampling variance status |
| `branch-summary.tsv` | NWKIT branch IDs, shared regimes, shift indicators and stability frequencies |
| `pdf` | One observed/fitted expression tree view per trait |

NWKIT IDs are level-order indices; GeneGalleon uses clade ranks. The orthogroup
summary validates the complete dated-tree topology, names and branch lengths
before translating IDs into `stat_branch` columns `ou_native_regime`,
`ou_native_is_shift`, `ou_native_selection_frequency`, and `ou_native_research_only`.
Tree summaries include `ou_native_num_shift` and `ou_native_num_regime`.
Branch colors can use `ou_native_regime`. No unsupported optima are converted to
zero or used as expression pie charts. The dedicated PDF uses fitted tip means.

Stage provenance includes input expression and tree, adapter source, NWKIT source
identity and all exposed search settings. Stale stage artifacts follow the normal
`artifact_stale_policy` rules. The standalone helper also accepts
`--resume-model previous.model.json`: reuse requires matching raw replicates,
adapter, normalized inputs, configuration and numerical implementation. This
reuses completed analyses; interrupted calibration checkpoints are not supported.
Failures preserve the previous artifact bundle.

The helper can be invoked directly inside a GeneGalleon container:

```bash
python /script/support/detect_ou_shift_native.py \
  --tree dated.nwk --traits expression.tsv --output-prefix family_ou \
  --criterion AICc --search-strategy auto --max-shifts auto --convergence
```

An optional `--regime-map` fits a supplied map without discovery. The standalone
helper accepts fixed `--alpha` and `--process-tip-variance` for reproducible
numerical comparisons. The NWKIT [native guide](https://github.com/kfuku52/nwkit/blob/master/NATIVE_SHIFT.md)
specifies the model and boundary conventions; the [adoption protocol](https://github.com/kfuku52/nwkit/blob/master/NATIVE_SHIFT_VALIDATION.md)
defines the remaining statistical validation. The research-status metadata is
preserved even though GeneGalleon now uses NWKIT as its OU backend. Docker
validation does not establish SIF compatibility.
