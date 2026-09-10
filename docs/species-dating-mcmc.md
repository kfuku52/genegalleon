# Species dating: parallel MCMC and advisory diagnostics

Species dating runs four independently initialized MCMCtree chains, with separate
burn-in and distinct explicit seeds. `mcmc_parallel_jobs=4` runs all four
concurrently when at least four CPUs are allocated. The workflow caps concurrent
chains at `GG_TASK_CPUS`; lower this setting further to reduce memory use.
Each child uses one OpenMP thread. `mcmc_seed=1729` supplies the first seed; the
remaining seeds are 1730, 1731, and 1732. Changing the seed changes the run identity.
The existing burn-in, sampling-frequency and sample-count options apply per chain.

Diagnostics are **warnings, not downstream gates**. All four complete chains are
checked using R `posterior`: rank-normalized split/folded R-hat < 1.01, bulk ESS
>= 400 and tail ESS >= 400 for every sampled parameter except the iteration
index. Mean and 2.5%/97.5% quantile MCSE are also written. Missing dependencies,
undefined diagnostics (including constant parameters), and unequal chain lengths
are reported as `diagnostic_failed`. Diagnostic thresholds do not prove model
adequacy or convergence. Tail ESS uses 5%/95% quantiles; it does not measure HPD
endpoint precision. The quantile MCSE columns are in the original internal units.

Diagnostics are retained in separate `diagnostics-*` attempts, including the R
session information; `convergence.json` identifies the current diagnostic directory.
Even with fewer than four complete chains, available chains are diagnosed when at
least two remain, but the overall status stays `incomplete`.

The MCMCtree sample file already excludes burn-in; GeneGalleon does not discard it
a second time. All valid recorded rows are retained, including any initial saved
state emitted by PAML. Files must be finite, have increasing iteration indices,
and reach the configured final iteration before a chain is marked complete.

Complete chains are pooled and summarized with PAML's `print=-1` mode. The published
FigTree tree and its HPD intervals therefore come from the pooled sample, not an
average of separately computed interval endpoints. If some chains fail, complete
chains still produce a provisional summary with state `incomplete`. If none
complete, or the pooled tree cannot be generated/validated, no new usable result
can be published. Numerical or file failures are not converted into valid trees.

## Evidence and restart

`output/species_tree/mcmctree_main/runs/<run-id>/` permanently retains the exact
internal-unit inputs, executed controls, independent raw sample files, full PAML
output, stdout/stderr, checkpoints, seeds, exit/failure state, and diagnostic
outputs. `manifest.json` hashes the inputs, executable and diagnostic adapters
and records the public/internal time factor. Each chain has separate attempt
directories. The ordinary temporary-directory cleanup does not remove this store.
Rebuilds do not erase earlier runs. Archive/delete evidence only deliberately.

Rerunning identical inputs reuses completed chains after verifying their stored
hashes (including stdout/stderr). Corrupt or incomplete restart receipts cannot
certify a chain. Incomplete or corrupted chains are rerun from their original seed into a
new attempt, retaining the previous attempt. This is a chain-level restart, not
binary checkpoint continuation. A resumed segment is never counted as a new
independent chain. PAML binary checkpoint continuation requires separate validation
of the installed dependency's RNG and sample-boundary behavior; checkpoints are
retained as evidence but are not automatically spliced by GeneGalleon.

`mcmctree_main/convergence.json` records the advisory state and the run location.
After conversion, `dated_species_tree.nwk.convergence.json` travels beside the
canonical dated tree. Tree hashes prevent unrelated/replaced trees from inheriting
a `passed` result. Missing or mismatched evidence yields `legacy_unverified` and
warns; diagnostic status does not block old trees. Existing artifact stale-input
policies still apply when an implementation or input changes. The genome and gene workflows warn when consuming
provisional dated trees, including runs that skip the MCMC producer. Raw public
`iq2mc.mcmctree.out` remains a public-unit tree summary for existing consumers;
the true unmodified PAML output is retained under `runs/`. Failed summary creation
is recorded as `failed`, rather than leaving a successful convergence label. Public
tree files are replaced as a bundle only after successful summary conversion and
validation; a failed attempt preserves the previous public files. Cached tree
conversions refresh their advisory metadata when the verified source tree is unchanged.
SIGTERM/SIGINT stop the child process groups and retain interruption evidence.

The separate, opt-in calibration experiments documented in
[species-tree-calibrations.md](species-tree-calibrations.md) remain independent of
this main-run diagnostic workflow and do not replace its public tree.

## Runtime and interpretation

The container R environment includes `r-posterior`. Rebuild the runtime to install
it; absence is an advisory diagnostic failure, not evidence of convergence.
Tests in `test_mcmctree_chains.py` exercise real PAML, warning-only publication,
chain reuse, corruption handling and diagnostic fixtures in a GeneGalleon runtime.

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
