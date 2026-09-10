# NWKIT OU shift analysis

`run_native_ou=1` in `gg_gene_evolution_entrypoint.sh` enables NWKIT's native
multivariate OU shift analysis. OU analysis remains disabled unless requested.
NWKIT replaces the former kfl1ou workflow and container dependency. The default
selection is **AICc with native-path**, with at most ten shifts per family.
This workflow choice does not establish statistical superiority: AICc selects a
predictive model and does not control false shift detections. See the
[validation record](native-ou-validation.md) for measured accuracy and limits.

The model shares shift locations across traits and allows separate OU covariance
parameters per trait. It supports sampling errors, additional estimated
observation variance, missing coordinates and fixed/stationary roots. It does
not estimate full cross-trait covariance. Trees must be rooted, binary,
ultrametric and have positive branches; invalid trees are rejected.

## Configuration

```bash
run_native_ou=1
native_ou_criterion="AICc"
native_ou_max_shifts=10
native_ou_convergence=0
native_ou_calibration_replicates=199
native_ou_calibration_level="0.05"
native_ou_bootstrap=0
native_ou_seed=1
native_ou_bootstrap_seed=2
native_ou_root_model="OUfixedRoot"
native_ou_estimate_measurement_error="yes"
native_ou_search_strategy="native-path"
native_ou_candidate_pool=24
native_ou_refit_budget=48
native_ou_screening_budget=2000
native_ou_beam_width=2
native_ou_replicate_separator="_"
treevis_branch_color="ou_native_regime"
```

`native-path` generates joint shift configurations along group-lasso regularization
paths, refits retained layouts without shrinkage, and selects the smallest AICc.
It supports `AIC` or `AICc` with convergence disabled. Search is approximate:
refit/screening budgets and finite optimizer iterations can limit coverage.
The shift cap applies at every tree size and does not guarantee that many true
shifts will be recovered.

For shared regimes and nested returns, set `native_ou_convergence=1` and choose
`auto`, `exhaustive`, or `lasso`. `auto` enumerates small spaces within NWKIT's
5,000-candidate budget, otherwise using beam/local search. These strategies also
support `BIC`, `pBIC`, and `bootstrap`. The last uses the experimental sequential
plug-in bootstrap with `native_ou_calibration_replicates` draws per test; it is
not a proven uniform error guarantee. Positive `native_ou_bootstrap` repeats the
entire selected procedure for stability frequencies, including inner calibration
when `native_ou_criterion="bootstrap"`.

## Migration from kfl1ou

Replace `run_l1ou=1` with `run_native_ou=1` and remove the old `l1ou_*`,
`large_tree_num_gene`, and `large_tree_max_nshift` settings. Configure the single
`native_ou_max_shifts` cap instead. Old RData fits and `l1ou_*` tables are not
reused. New results use `ou_native/<family>_ou_native.*`; downstream summaries and
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
  --criterion AICc --search-strategy native-path --max-shifts 10
```

An optional `--regime-map` fits a supplied map without discovery. The standalone
helper accepts fixed `--alpha` and `--process-tip-variance` for reproducible
numerical comparisons. The NWKIT [native guide](https://github.com/kfuku52/nwkit/blob/master/NATIVE_SHIFT.md)
specifies the model and boundary conventions; the [adoption protocol](https://github.com/kfuku52/nwkit/blob/master/NATIVE_SHIFT_VALIDATION.md)
defines the remaining statistical validation. The research-status metadata is
preserved even though GeneGalleon now uses NWKIT as its OU backend. Docker
validation does not establish SIF compatibility.
