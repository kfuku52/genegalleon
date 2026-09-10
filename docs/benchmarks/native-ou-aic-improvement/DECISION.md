# Decision: expose the native AIC path explicitly

The new candidate path materially improves the primary 100-tip/10-shift weak
condition and reaches kfl1ou-like point estimates there. It is retained as
`--selection native --criterion AIC --search-strategy native-path`.
The default `auto` continues to use exhaustive search when feasible, otherwise
the existing beam search. No information-criterion penalty or alpha bound changed.

## Independent evidence

All 270 primary and 36 scope-check runs completed successfully.

- Weak shifts: F1 rises from 0.217 to 0.434 (paired difference 0.217,
  95% seed-bootstrap interval [0.172, 0.263]); mean RMSE falls from 1.119 to
  0.774 (difference -0.346, interval [-0.408, -0.285]). kfl1ou gives F1 0.409
  and RMSE 0.752. The path/kfl1ou RMSE ratio is 1.029, interval [0.994, 1.065].
- Strong shifts: F1 rises from 0.583 to 0.638 (difference 0.054,
  interval [0.014, 0.095]); kfl1ou gives 0.644. Mean RMSE is 0.675 for the path,
  0.842 for the original and 0.758 for kfl1ou, but paired RMSE differences
  remain uncertain. The path-minus-kfl1ou F1 interval is [-0.028, 0.017].
- Nonnull median time is 2.11 seconds for the path, 4.41–4.55 for the original,
  and 8.46–9.81 for kfl1ou. Native peak RSS medians are about 112 MiB.
  These are different configured procedures, not equivalent-output kernels.

The prespecified kfl1ou point-estimate margins are met in both primary nonnull
conditions. Their uncertainty intervals do not establish those margins uniformly;
this is not a formal equivalence or general-superiority claim.

## Why auto was not changed

The smaller scope check is mixed. On balanced trees with five strong shifts,
F1 decreases from 0.747 to 0.667 and RMSE increases from 0.570 to 0.619. On
random trees with five strong shifts, RMSE increases from 0.638 to 0.812 even
though F1 improves slightly. There are only three seeds per condition, but this
provides insufficient evidence to replace the general default.

All three AIC procedures select shifts in 30/30 global-null datasets. Mean
selected counts are 9.93 (original), 10.00 (path) and 7.57 (kfl1ou). The path's
null mean RMSE is also slightly worse: 0.717 versus 0.707 and 0.678. Ordinary
AIC selection here must not be interpreted as calibrated evidence that shifts
exist. This experiment does not validate the separate bootstrap procedure.

## Integration after measurement

The frozen candidate snapshot contains the measured search implementation.
The benchmark adapter calls `sparse_native_search` directly. After reviewing
independent results, only the CLI/runner default dispatch was changed back to
beam; explicit `native-path` retains the measured implementation. The path,
whitening, likelihood, AIC and heuristic modules remain byte-identical to the
measured snapshot. `final-integration.patch` records this dispatch-only change;
`final-source-manifest.json` fingerprints the delivered source. Focused routing
and search tests were rerun after this decision. No independent data were used
to retune the search algorithm.
