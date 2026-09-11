# NWKIT evolutionary covariance and simulation proposal

Status: the core model, joint screening, exact shared-alpha fast path, vector
pruning, correlated bootstrap, and `nwkit shift-simulate` are implemented locally
in NWKIT as of 2026-09-11. GeneGalleon forwards the opt-in covariance/alpha options;
workflow defaults remain unchanged. See the [integration guide](../../native-ou-shifts.md)
and NWKIT's `SHIFT_COVARIANCE.md` and `reviews/trait-covariance-2026-09-11/` for
the implemented contract and implementation-based validation. Shrinkage, low-rank
estimators and the broader statistical adoption study below remain future work.

The following retains the original design proposal based on local NWKIT
inspection and the controlled 100-tip prototype in this directory.

## Model contract

Add full evolutionary diffusion covariance to native shift fitting. Keep the
distinction between evolutionary covariance and measurement covariance explicit.
Initially share the evolutionary covariance across regimes; regime-specific
matrices introduce a separate model-selection problem and substantially more
parameters. Shared alpha and trait-specific positive alpha should be explicit
model choices, not selected silently to obtain a faster fit.

For diagonal attraction A=diag(alpha_i), branch length t, diffusion D and regime
optimum theta, use transition F=diag(exp(-alpha_i*t)), intercept (I-F)theta,
and innovation Q_ij=D_ij*(1-exp(-(alpha_i+alpha_j)*t))/(alpha_i+alpha_j).
Use stable expm1 evaluation and the defined limit at zero. Parameterize positive
definite D by Cholesky factors; do not estimate an arbitrary stationary covariance
and assume it implies valid diffusion for unequal alpha. For stationary roots,
C_ij=D_ij/(alpha_i+alpha_j); fixed roots have a different covariance contract.
Zero-process and alpha-boundary models need explicit support/diagnostics rather
than added numerical jitter or silently clipped estimates.

Do not initially estimate a fully coupled attraction matrix. Full diffusion
covariance already answers the requested question; off-diagonal attraction adds
different biological assumptions and an additional identifiability burden.

Observation model: retain known sampling variances and estimated additional
diagonal observation variances. Known correlated sampling covariance can be an
explicit input. Estimating an unrestricted measurement covariance simultaneously
with evolutionary covariance requires an identifiability assessment and suitable
replicate data; do not conflate these matrices.

## Reuse in NWKIT

- `nwkit/vector_gaussian.py`: vector affine-Gaussian pruning and conditioning;
  supports joint error covariance, with O(nodes * traits^2) storage.
- `nwkit/vector_ou_fit.py`: OU transition construction, parameter optimization,
  and incomplete-data machinery. Its stationary, single-optimum fitting interface
  is not directly a fixed-root, multiple-regime shift fitter; reuse kernels, not
  unmodified fit semantics.
- `nwkit/vector_simulation.py`: joint process simulation and observation errors
  while retaining missing masks.
- `nwkit/shift_native_screen.py`: group-lasso candidate proposal, currently based
  on per-trait covariance whitening; extend to joint trait whitening.
- `nwkit/shift_native_bootstrap.py`: full-search replay and seed handling;
  simulation must draw correlated evolutionary innovations.
- Extend fixed-layout likelihood, search scoring, output schema, provenance,
  resume validation, bootstrap, tests and documentation together. The existing
  output currently declares diagonal covariance explicitly.

## Simulation interface

Expose a documented simulation API and a CLI (name to be finalized) accepting
either explicit generating parameters or a fitted model JSON. Export observed
trait TSV, latent tip/node values, true branch/regime mapping, both covariance
matrices, root treatment, complete generating parameters and independent seeds.
Support no shift, multiple/nested shifts, shared regimes, missingness, known
sampling error and replicate generation. Differentiate unconditional draws from
posterior conditional draws.

Share process definitions with inference for consistency, but validate them
against an independent dense covariance/likelihood oracle on small trees and
independent branch-recursion simulations. Self-consistency alone does not prove
correctness.

## Improving speed without removing covariance

1. Shared-alpha, complete-data, no-observation-error exact fast path: separable
   tip/trait covariance, tree whitening, profiled GLS means and closed-form ML
   covariance. Cache tree geometry; reuse factorizations only when their inputs
   are unchanged. This is the algebra exploited by the preceding benchmark.
2. General path for trait-specific alpha, noise or missingness: vector tree
   pruning instead of a dense (n*p)-square covariance. At fixed parameters and
   bounded tree degree, small dense p-square operations give approximately
   O(n*p^3) time and O(n*p^2) storage per likelihood evaluation. Optimization and
   regime-mean solves add costs; these are not end-to-end runtime predictions.
3. Jointly whiten evolutionary residuals across traits when proposing shifts.
   Do not screen solely with the diagonal model: the benchmark's contrast shifts
   can disappear before a full-covariance refit is attempted. Retain a controlled
   union of proposal paths, then refit all finalists with their covariance free.
4. Warm-start nearby layouts, use multi-start safeguards, cache valid geometry,
   and update mean-design factorizations after local changes. Check finalists
   against exhaustive small-tree references to quantify missed candidates.
5. Parallelize independent bootstrap draws with reproducible per-draw seeds and
   controlled BLAS threads. Every draw regenerates candidates and refits covariance;
   reusing the observed-data selected branches would change calibration.

## Improving statistical performance

Evaluate shrinkage toward diagonal covariance and a low-rank-plus-diagonal model
as explicit additional estimators. They can reduce noisy covariance estimation
with many traits, but improvement in shift detection is an empirical question.
Ten traits already have 55 full covariance parameters; twenty have 210.

Estimate covariance from phylogenetically adjusted, fitted-mean residuals, not
raw tip correlations. Fit mean shifts and covariance together so that genuine
shifts are not permanently absorbed into the null covariance. Tune shrinkage or
rank inside the procedure replayed by bootstrap; do not attach ordinary ML AICc
parameter counting unchanged to a penalized estimator. Avoid blindly discarding
low-variance principal components: they carry the strongest signal in the
opposed-shift experiment.

References: [mvMORPH pruning documentation](https://search.r-project.org/CRAN/refmans/mvMORPH/html/pruning.html),
[multivariate penalized GLS documentation](https://jclavel.r-universe.dev/mvMORPH/doc/manual.html),
[Ledoit–Wolf covariance shrinkage](https://ledoit.net/ole1a.pdf).
These support algorithmic directions, not performance guarantees for NWKIT.

## Validation sequence

1. Fixed layouts: dense likelihood agreement, diagonal/one-trait nesting where
   model assumptions match, trait permutation and unit-rescaling invariance,
   fixed/stationary roots, missing values and known observation covariance.
2. Generative tests: empirical means/covariances, seed reproducibility,
   export/import consistency and independent-oracle checks.
3. Exhaustive small-tree shift searches: likelihood and selected-layout agreement;
   verify inference parameter counts and diagnostics separately.
4. 100-tip paired experiments with operational fitted-null bootstrap, not the
   previous generating-null oracle calibration. Vary correlation strength and
   sign/block structure, signal direction/size, trait count, tree shape,
   measurement error, missingness and multiple shifts. Include unequal-alpha
   generating conditions when assessing the shared-alpha model.
5. Report true-branch recall, false discovery proportion, null familywise error,
   branch-localization error, covariance-estimation error, failures, wall time and
   peak memory. Compare exact and accelerated versions with equivalent outputs;
   quantify any screening approximation separately.
6. Integrate into GeneGalleon only after the owning NWKIT implementation and
   runtime validation exist. Previous prototype speed ratios do not predict
   production runtime. Preserve Docker versus SIF validation distinctions.

Recommended order: correlated simulation and fixed-layout full-covariance
reference; exact shared-alpha fast path and general observation-error path;
joint candidate search and fitted-null calibration; optional regularization;
GeneGalleon integration.
