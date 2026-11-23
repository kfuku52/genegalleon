# Species dating: parallel MCMC and advisory diagnostics

Saved dated trees can be redrawn with spelled-out interval legends, geological
period backgrounds and dataset-labelled BUSCO panels. See
[Dated-tree presentation](dated-tree-presentation.md).

The container builds IQ-TREE 3 and PAML's MCMCtree from their respective upstream
branches. GeneGalleon uses IQ-TREE's IQ2MC mode to estimate branch lengths,
the gradient and Hessian, then adapts the generated control for upstream MCMCtree.
The execution control passes the Hessian path in `usedata = 2` and makes the
IQ2MC conditional birth–death prior explicit. An explicit `C` or `M` in the
input control is preserved. The original IQ2MC files remain reusable together,
including the matching rooted tree and dummy alignment. Updating MCMCtree alone
does not require regenerating the IQ2MC inputs.

Species dating calls `nwkit mcmctree-run` for four independently initialized
chains. `mcmc_parallel_jobs` sets the worker budget for sampling and diagnostics
and is capped by the allocated CPUs. Sampling uses up to four workers, with one
thread per child. Additional workers can compute R diagnostics in parallel.
The default budget is four. `mcmc_seed=1729` supplies the first seed. The burn-in,
sampling interval, and minimum saved samples use
`mcmc_burnin`, `mcmc_sampfreq`, and `mcmc_nsample`. The default minimum is
2,020,000 generations per chain.

NWKIT diagnoses all four complete chains with rank-normalized R-hat < 1.01,
bulk ESS >= 400 and tail ESS >= 400. It checks the initial run, extends when
needed, and checks again using all evenly spaced draws. Intermediate checks
compute R-hat and ESS. The final diagnostics also report MCSE. Set
`mcmc_max_generations` to bound the total generations per chain including
burn-in. With `auto`, the cap is the burn-in plus 64 times the initial sampling
generations within PAML’s limits. Near the stopping criteria, NWKIT checks more
frequently. Set `mcmc_burnin_percentages="2.5,5,10,20,40"` to try excluding
the same initial window from all chains when the full traces fall short. The
smallest passing cutoff is checked again after adding draws, and a full trace
that passes has priority. The `sample_window` in the tree’s convergence JSON
records the original burn-in, effective cutoff and generations used. See the
[NWKIT guide](https://github.com/kfuku52/nwkit/blob/master/docs/guides/MCMCTREE_RUN.md).

NWKIT resumes each chain from its own final checkpoint and pools the draws to
generate a FigTree tree with 95% HPD intervals. GeneGalleon converts that tree to
the original analysis units. At the generation limit, an unmet convergence
criterion produces a provisional tree with an advisory status. A failed chain,
diagnostic, or pooled summary preserves the previous published tree.

## Evidence and restart

`output/species_tree/mcmctree_main/runs/<run-id>/` retains inputs in internal
units, executed controls, samples, logs, checkpoints and diagnostics. The run
manifest records input and binary hashes. Rebuilding with the same inputs and
settings reuses completed chain segments. An interrupted segment starts again
from its preceding verified checkpoint. Changing the inputs or sampling
settings creates a separate run directory.

`mcmctree_main/convergence.json` records the advisory state and the run location.
After conversion, `dated_species_tree.nwk.convergence.json` travels beside the
canonical dated tree. Tree hashes prevent unrelated/replaced trees from inheriting
a `passed` result. Missing or mismatched evidence yields `legacy_unverified` and
warns; diagnostic status does not block old trees. Existing artifact stale-input
policies still apply when an implementation or input changes. The genome and gene workflows warn when consuming
provisional dated trees, including runs that skip the MCMC producer. Raw public
`iq2mc.mcmctree.out` remains a public-unit tree summary for existing consumers;
the unmodified PAML output is retained under `runs/`. GeneGalleon publishes the
public summary and convergence record together after conversion and validation.
Cached tree conversions refresh their advisory metadata when the verified source
tree is unchanged. SIGTERM and SIGINT stop active sampler processes and retain
completed segments for restart.

The separate, opt-in calibration experiments documented in
[species-tree-calibrations.md](species-tree-calibrations.md) remain independent of
this main-run diagnostic workflow and do not replace its public tree.

## Runtime and interpretation

The container R environment includes `r-posterior`. Rebuild the runtime to install
it. The runner requires this package to compute convergence diagnostics.
Runtime tests exercise the NWKIT command on real IQ2MC inputs and check the
public time scale and advisory status.

The four chains increase MCMC CPU time roughly fourfold at unchanged lengths;
wall time depends on concurrent resources and memory pressure. No speedup is
promised without representative measurements. Estimate storage from actual sample
row width and parameter count; retain raw samples rather than thinning them merely
to make diagnostics pass.

Dates and intervals remain conditional on topology, calibration, clock and other
model assumptions. Prior-only and calibration-sensitivity analyses remain needed.
Downstream use of a summary tree does not propagate species-date uncertainty.

References: [PAML MCMCtree documentation](https://github.com/abacus-gene/paml/wiki/MCMCtree),
[rank-normalized R-hat and ESS](https://avehtari.github.io/rhat_ess/rhat_ess.html),
[`posterior` R-hat](https://mc-stan.org/posterior/reference/rhat.html).
