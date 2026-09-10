# NWKIT tree and regression integration

GeneGalleon uses NWKIT for dated-tree drawing, MAD/midpoint rooting and root
comparison, and species-tree regression. NOTUNG still generates candidate roots.

## Dated species trees

The existing `run_plot_mcmctreer` switch and public PDF filename are retained.
`workflow/support/plot_dated_tree.py` draws with `nwkit draw`, preferring the
public-unit `mcmctree_95CI.nhx` when available, otherwise the dated Newick tree.
Supplied credible intervals are drawn; absent intervals are not invented.
Branch lengths and the scale bar are in Ma. Layout follows NWKIT, including a
branch-length scale bar instead of the former R plot's time axis.

## Gene-tree rooting

`workflow/support/species_tree_guided_gene_tree_rooting.py` calls `nwkit root`
for MAD and midpoint roots. Selection keeps this priority:

1. MAD if its root edge matches a NOTUNG candidate.
2. Midpoint if its root edge matches a NOTUNG candidate.
3. The first naturally ordered NOTUNG candidate.
4. MAD when there are no NOTUNG candidates.

Compatibility compares the root edge independently of the position along that
edge. Candidate tip sets must match. NOTUNG inference is unchanged. Beside the
existing `.root.txt` log, `.root.tsv` and `.root.pdf` record NWKIT's MAD/midpoint
comparison; the log records the selected root and NOTUNG compatibility. The
comparison figure covers MAD and midpoint, not every NOTUNG candidate.

## Regression

Expression-trait `pgls_methods` accepts `rsc`, `species-nwkit`, their
comma-separated combination, or `all`. `all` runs these two methods. The retired
`species-rphylopars` method, `rphylopars_sampling_covariance` option, and
`pgls_species_rphylopars` output member are no longer produced or accepted.
The unused legacy gene-tree R regression helper is also removed. See
[the expression-trait recipe](common-workflow-recipes.md) for supported models,
replicate uncertainty, categorical predictors, and species aggregation.

Orthogroup copy-number/trait regression defaults to `trait_value ~ log1p(copy_number)`
with NWKIT Brownian restricted-maximum-likelihood (REML) GLS and Wald inference.
Copy numbers are exact predictors. For one complete observation per species,
the former Rphylopars fit automatically disabled phenotype error and used REML;
the new fit preserves its coefficient and standard-error test. Using ML here
would underestimate the residual variance relative to the former test. R remains
responsible for preparation and summary graphics. Likelihood and engine-specific
summary fields retain NWKIT's definitions and should not be compared directly
with historical Rphylopars likelihood summaries.

The table retains `R2`, `pval`, `logLik`, `PCC`, `OLS_slope`, adjusted p-values,
and fit status. Engine-specific `R2adj`, `sigma`, `Fstat`, `AIC`, and `BIC`
columns are removed. Native coefficient, standard error, statistic, degrees of
freedom, evolutionary rate, covariance/model/inference metadata, and optimizer
and small-sample diagnostics are included. `fit_mode` is `nwkit_brownian_reml`.
Use the full result table and fit status when interpreting screened results.

Binary, Poisson and negative-binomial response families can now be selected
per trait; non-Gaussian models use NWKIT Laplace ML with the appropriate link.
The optional joint-selection stage calls `nwkit regress-select` and emits
exploratory coefficients, nested-CV predictions and selection frequencies.
See [copy-number models](copy-number-trait-models.md) for configuration and the
explicit distinction between ordinary tests and selection outputs.

The optional [BUSCO quality diagnostic stage](copy-number-quality-diagnostics.md)
adds separate `log1p`-count PGLS comparisons with BUSCO as a covariate or a
cohort restriction. It preserves the original root/shared Brownian history
when dropping species and uses full-precision branch lengths. This is a
sensitivity analysis, not a copy-number measurement-error model.

## Failure handling and cache invalidation

The dated-tree PDF and its summary copy are published together. Rooting computes
all family results in temporary directories before replacing the report/tree
directories. Copy-number regression publishes its three TSVs and two figures as
one bundle; species-tree PGLS likewise stages all ten outputs. Ordinary write or
plot failures preserve the previous bundle. Output paths that alias inputs or
other output members are rejected. A missing NOTUNG candidate directory,
duplicate candidate tip labels, mismatched tip sets, or an unresolved candidate
root is an error; these conditions do not silently select MAD.

Dated plotting, rooting, and copy-number regression include the NWKIT package
version, installed Python-source digest, and available build revision in their
artifact contracts. Rooting also records the adapter source. Dependency or
adapter changes are therefore subject to the configured stale-artifact policy.
The significant-results table (including a header-only result) and SVG summary
are required copy-number bundle members, so deleting either invalidates the
cached stage.

## Validation

A review after migration reproduced and fixed an unintended ML setting in the
GeneGalleon adapter (the previous implementation used REML),
stale dependency caches, acceptance of malformed NOTUNG candidate inputs,
input/output path collisions, and partial replacement after late failures.

For complete one-row-per-species data, direct comparisons with the previous
Rphylopars call agree on coefficient, standard error, p-value and R² for 4-, 8-,
and 12-species examples. The original ML migration gave p=0.08870047 instead of
0.15733516 in the four-species example; REML restores p=0.15733516. A permanent
regression test independently computes Brownian GLS using a covariance matrix,
including the residual degrees of freedom. Native REML likelihood values differ
from Rphylopars' likelihood summary and remain explicitly engine-specific.

The full strict Python suite passed in GeneGalleon's Docker runtime, followed
by a real NOTUNG-to-NWKIT integration check. The complete declared R validation
suite passed, including package checks and parsing all 23 R helper scripts.
For the final R run, workflow sources were copied into an isolated validation
snapshot and its treevis package was installed there to avoid mixing concurrent
workspace edits with an older installed package. Dated-tree and root-comparison
PDFs were also inspected visually. NWKIT drawing tests passed 259 cases with two
unrelated skips (CairoSVG/native Cairo support and a case-sensitive filesystem).

The runtime was `local/genegalleon:nwkit-migration-dev`, derived from
`local/genegalleon:intron-asr-dev` with the local NWKIT `draw_render.py` fix for
clipped root credible intervals. That dependency fix remains in the sibling
NWKIT checkout and needs publishing there before a fresh build from the moving
upstream branch includes it. No revision was pinned in GeneGalleon defaults.
SIF execution and a fresh full container build were not run.

Shell syntax, changed-file Python lint, and Actions/composite YAML and expression
checks passed. Host actionlint with external ShellCheck did not finish within
the timeout; the YAML/expression check was rerun with that external integration
disabled. Repository-wide Python lint has one unrelated existing import-order
finding in `test_synteny_cutoff_metadata.py`.
