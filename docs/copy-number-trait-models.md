# Copy-number associations and joint predictor selection

Genome evolution has two independently enabled stages. Both use the shared
orthogroup copy-number matrix and species traits; neither runs by default.

## Prespecified orthogroup-by-trait tests

```bash
run_orthogroup_copy_number_trait_pgls=1
orthogroup_copy_number_trait="height,present,other_copy_count"
orthogroup_copy_number_trait_response_families="height=gaussian,present=binomial,other_copy_count=negative-binomial"
```

The historical `trait_pgls` stage name and result paths are retained. The model
is `trait ~ log1p(orthogroup_copy_number) + intercept`, one orthogroup at a time, using
the dated species-tree Brownian covariance. The new mapping selects the
**response** distribution by trait column:

| Family | Trait | Estimation |
|---|---|---|
| `gaussian` | Continuous numeric | NWKIT PGLS, REML |
| `binomial` | Numeric 0/1 | NWKIT logit phylogenetic GLMM, Laplace ML; reference 0 |
| `poisson` | Non-negative integer | NWKIT log-link phylogenetic GLMM, Laplace ML |
| `negative-binomial` | Non-negative integer with dispersion | NWKIT log-link NB2 phylogenetic GLMM, Laplace ML |

Unspecified traits retain Gaussian behavior. Numeric coding never silently
selects a distribution. Invalid binary/count values and duplicate or unknown
mapping keys fail. Both association and selection stages use `log1p(copy_number)` (natural log)
before fitting; selection then standardizes within each training fold. The raw
copy-number matrix remains unchanged. Coefficients, PCC and OLS slopes refer
to the transformed predictor, not a one-copy increase. Output
`predictor_transform=log1p` records this convention. The transform compresses
large counts but does not eliminate extrapolation or guarantee better prediction.
Coefficient regularization is explicitly disabled for these ordinary Wald
coefficient tests. Sparse binary data may have unavailable inference; inspect
`status`, `inference_status` and optimizer diagnostics rather than interpreting
missing P-values as non-significance. Successful finite tests receive the
existing global and per-trait BH adjustments; the significant table uses the
global adjusted P-value cutoff. This is not a variable-selection stage.

The result adds `response_family`, `link_function` and `coefficient_penalty`.
`covariance_estimator` is Gaussian REML or Laplace ML as appropriate. The old
PCC and OLS slope remain descriptive summaries; they do not define the GLMM
test. Non-Gaussian fields such as R-squared follow NWKIT's definitions or may be
unavailable. Distributions, mapping and NWKIT source identity are part of the
artifact provenance contract.

## Exploratory joint selection

```bash
run_orthogroup_copy_number_trait_selection=1
orthogroup_copy_number_trait="height,present,other_copy_count"
orthogroup_copy_number_trait_response_families="height=gaussian,present=binomial,other_copy_count=negative-binomial"
orthogroup_copy_number_trait_selection_folds="input/species_trait/phylogenetic_folds.tsv"
orthogroup_copy_number_trait_selection_strengths="1,0.1,0.01"
orthogroup_copy_number_trait_selection_l1_ratios="1,0.5"
orthogroup_copy_number_trait_selection_prediction="conditional"
```

The folds file is workspace-relative or absolute and has `leaf_name` and `fold`
columns. Supply at least three phylogenetic groups. Each nested training split
must have at least four tips and a variable trait (both binary classes for
binomial). Grouping does not imply phylogenetic independence. Choose biologically
meaningful clades rather than grouping on the observed trait.

This stage calls `nwkit regress-select` and fits all selected orthogroups jointly.
It supports lasso (`l1_ratio=1`) and elastic net (`0<l1_ratio<1`), with an
unpenalized intercept. Scaling, exclusion of constant predictors and fitting
occur inside each training fold. Nested group CV tunes penalties internally,
then evaluates untouched outer groups. The full-data model is tuned separately.
The generic NWKIT command also supports prespecified unpenalized covariates;
the GeneGalleon adapter currently supplies only orthogroup predictors.

Both stages share `orthogroup_copy_number_trait_family_ids`,
`orthogroup_copy_number_trait_family_file` and
`orthogroup_copy_number_trait_max_families`. These are user-defined scope
restrictions, not a P-value filter. Do not feed in candidates selected using the
same trait data and then describe nested-CV performance as free of that screening
bias. The matrix-preparation stage runs when selection alone is enabled.

Missing response species are excluded per trait and counted in the manifest.
Both stages preserve the original tree root, shared ancestral branches and
double-precision branch lengths when excluding species, so subsetting does not
change the retained species' Brownian covariance. Selection input tables reject
empty/duplicate headers and rows with inconsistent field counts.
Predictor counts must be finite non-negative integers. All tree species must
have a fold assignment. The adapter uses a predictor-name file to support
large OG sets without command-line length limits. Current NWKIT selection uses
a dense covariance and accepts at most 500 analyzed species per trait; predictor
count can exceed species count. No job is launched by merely changing these
default-off settings.

`conditional` prediction includes the training random-effect mode projected
through the tree covariance; `fixed` uses only fixed effects. Both apply the
inverse link to a plug-in predictor rather than integrating latent uncertainty.
They evaluate different prediction tasks. Count scores are squared error and
binary scores are log loss. Outer-fold baseline predictions use the same tree
without the orthogroup predictors, allowing predictive improvement to be assessed.
The mean nested-CV loss is descriptive predictive
performance, not a coefficient hypothesis test.

Outputs are under
`workspace/output/genome_evolution/orthogroup_copy_number/trait_selection/`:

- `manifest.tsv`: traits, response families, analyzed/missing species counts
  and current result prefixes (`trait_0001`, etc.). This manifest is authoritative
  for the current run; files from older prefixes are not implicitly current.
- `predictors.tsv`: internal predictor names mapped to original orthogroup IDs.
- Per trait: `.coefficients.tsv`, `.path.tsv`, `.cv.tsv`, `.predictions.tsv`,
  `.stability.tsv` and `.metadata.json`. Coefficient/path/stability tables also
  include `Orthogroup` for direct interpretation.

Every trait must complete before the adapter publishes the bundle. Handled
fitting or write failures preserve previous outputs. The whole output directory
is covered by provenance so deleted or changed members invalidate the cache.
Unrelated and historical files are preserved.

No post-selection P-values, adjusted P-values or confidence intervals are
produced. Outer-fold selection frequency is descriptive, not formal FDR control.
In particular, refitting selected OGs and applying BH does not repair selection
bias. Selection-aware inference is reserved for a separately validated extension;
existing fixed-model bootstrap tests do not account for this selection process.

When joint selection is the only enabled workflow step (all other `run_*`
settings are zero), the entrypoint accepts existing dated-tree, selected
copy-number and trait inputs without requiring CDS/protein FASTA files. Mixed
runs retain their normal sequence-input validation.

## Confidence intervals and calibration

Result tables preserve `confidence_interval_lower`, `confidence_interval_upper`,
and `confidence_level=0.95` from NWKIT. Intervals are pointwise on the model's
coefficient/link scale, not simultaneous intervals adjusted by global BH.
Significance continues to use the global adjusted P value. See
[regression calibration](regression-calibration.md) for small-sample evidence and limitations.
