# AICc added to the OU shift comparison

Same 100-tip/10-shift nonnull fixtures as the original seven-method comparison; three seeds per effect. Precision, recall, F1 and RMSE are averages across datasets. Higher F1 and lower RMSE are preferable. Native beam and path searches are distinct methods.

| Method | Weak F1 | Weak RMSE | Strong F1 | Strong RMSE |
|---|---:|---:|---:|---:|
| NWKIT beam bootstrap (B=19) | 0.000 | 1.275 | 0.412 | 1.145 |
| NWKIT beam pBIC | 0.100 | 1.207 | 0.378 | 1.504 |
| NWKIT beam BIC | 0.000 | 1.261 | 0.439 | 0.920 |
| NWKIT beam AIC | 0.142 | 1.054 | 0.433 | 0.946 |
| NWKIT beam AICc | 0.159 | 1.061 | 0.439 | 0.935 |
| NWKIT path AIC | 0.367 | 0.653 | 0.433 | 0.595 |
| NWKIT path AICc | 0.367 | 0.653 | 0.433 | 0.595 |
| kfl1ou pBIC | 0.096 | 1.356 | 0.367 | 1.448 |
| kfl1ou BIC | 0.000 | 1.261 | 0.390 | 0.723 |
| kfl1ou AIC | 0.296 | 0.700 | 0.421 | 0.746 |
| kfl1ou AICc | 0.296 | 0.679 | 0.421 | 0.746 |

## No-shift cases

Same thirty no-shift fixtures for every row. The global-null gate row uses 19 draws from the preceding study. No null data were used to tune the AICc formula or search.

| Method | Any false shift | Mean false branches | Mean RMSE |
|---|---:|---:|---:|
| NWKIT beam AIC | 30/30 | 9.93 | 0.707 |
| NWKIT beam AICc | 30/30 | 7.50 | 0.657 |
| NWKIT path AIC | 30/30 | 10.00 | 0.717 |
| NWKIT path AICc | 30/30 | 6.47 | 0.622 |
| kfl1ou AIC | 30/30 | 7.57 | 0.678 |
| kfl1ou AICc | 30/30 | 5.33 | 0.618 |
| NWKIT path AIC + gate (B=19) | 3/30 | 1.00 | 0.180 |

AICc reduces the average number of false branches and null prediction error, but all three AICc procedures still select at least one shift in 30/30 null datasets (95% Wilson interval for the rate: 88.65–100%). Thus it does not solve the false-positive problem in this setting. The gate reduces null detections, but its substantial nonnull power cost was measured on different nonnull fixtures; see the [gate study](../native-ou-global-null-gate/RESULTS.md).

## Full nonnull recovery table

| Effect | Method | Precision | Recall | F1 | Mean false branches | Mean K | RMSE |
|---|---|---:|---:|---:|---:|---:|---:|
| 2 | NWKIT beam bootstrap (B=19) | 0.000 | 0.000 | 0.000 | 0.667 | 0.667 | 1.275 |
| 2 | NWKIT beam pBIC | 0.100 | 0.100 | 0.100 | 6.000 | 7.000 | 1.207 |
| 2 | NWKIT beam BIC | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 1.261 |
| 2 | NWKIT beam AIC | 0.153 | 0.133 | 0.142 | 7.667 | 9.000 | 1.054 |
| 2 | NWKIT beam AICc | 0.198 | 0.133 | 0.159 | 6.000 | 7.333 | 1.061 |
| 2 | NWKIT path AIC | 0.367 | 0.367 | 0.367 | 6.333 | 10.000 | 0.653 |
| 2 | NWKIT path AICc | 0.367 | 0.367 | 0.367 | 6.333 | 10.000 | 0.653 |
| 2 | kfl1ou pBIC | 0.178 | 0.067 | 0.096 | 3.667 | 4.333 | 1.356 |
| 2 | kfl1ou BIC | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 | 1.261 |
| 2 | kfl1ou AIC | 0.333 | 0.267 | 0.296 | 5.667 | 8.333 | 0.700 |
| 2 | kfl1ou AICc | 0.333 | 0.267 | 0.296 | 4.667 | 7.333 | 0.679 |
| 6 | NWKIT beam bootstrap (B=19) | 0.429 | 0.400 | 0.412 | 2.667 | 6.667 | 1.145 |
| 6 | NWKIT beam pBIC | 0.392 | 0.367 | 0.378 | 5.333 | 9.000 | 1.504 |
| 6 | NWKIT beam BIC | 0.444 | 0.433 | 0.439 | 4.000 | 8.333 | 0.920 |
| 6 | NWKIT beam AIC | 0.433 | 0.433 | 0.433 | 5.667 | 10.000 | 0.946 |
| 6 | NWKIT beam AICc | 0.444 | 0.433 | 0.439 | 5.333 | 9.667 | 0.935 |
| 6 | NWKIT path AIC | 0.433 | 0.433 | 0.433 | 5.667 | 10.000 | 0.595 |
| 6 | NWKIT path AICc | 0.433 | 0.433 | 0.433 | 5.667 | 10.000 | 0.595 |
| 6 | kfl1ou pBIC | 0.367 | 0.367 | 0.367 | 6.333 | 10.000 | 1.448 |
| 6 | kfl1ou BIC | 0.417 | 0.367 | 0.390 | 4.000 | 7.667 | 0.723 |
| 6 | kfl1ou AIC | 0.444 | 0.400 | 0.421 | 4.333 | 8.333 | 0.746 |
| 6 | kfl1ou AICc | 0.444 | 0.400 | 0.421 | 4.333 | 8.333 | 0.746 |

All 120 newly measured runs completed. The six repeated kfl1ou AIC fits reproduce the original predictions, likelihoods and recovery metrics. The seven original procedures retain their frozen accuracy results; the AICc/path additions use the same inputs and runtime image. Runtime measurements come from separate batches, some under concurrent validation load; no speed ranking is inferred here. Individual wall/CPU/RSS records remain in the raw runs.

The nonnull comparison has only three seeds per condition, the cap equals the true shift count, and these fixtures were used during native-path development. Null rates use thirty seeds; see the Wilson intervals in `summary.json`. These results do not establish general superiority, a calibrated significance test, or production adoption. See [protocol](README.md), [raw metrics](metrics.csv), and [score validation](validation/native-kfl-score-agreement.json).
