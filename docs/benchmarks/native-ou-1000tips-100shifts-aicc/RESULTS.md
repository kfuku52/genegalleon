# 1,000-tip / 100-shift AICc timing comparison

One trait; balanced tree; unknown alpha and process variance refitted; max 100 shifts; three datasets per effect. Cache warmups excluded. The two procedures can select different models.

| Effect | Method | Complete | Wall median (s) | Wall range (s) | CPU median (s) | Peak RSS median (MiB) |
|---|---|---:|---:|---:|---:|---:|
| 2 | path AICc | 3/3 | 61.29 | 46.89–61.68 | 61.19 | 246.8 |
| 2 | beam AICc | 3/3 | 340.84 | 322.05–346.28 | 340.77 | 252.6 |
| 6 | path AICc | 3/3 | 66.10 | 63.53–67.13 | 66.08 | 248.8 |
| 6 | beam AICc | 3/3 | 341.45 | 316.58–343.96 | 341.34 | 252.2 |

In these workloads, path used about one fifth of the beam search wall time: ratios of beam/path medians were 5.56 and 5.17 at effect scales 2 and 6, respectively. The procedures produced different models; this is a configured-procedure timing comparison.

Medians are among completed runs; interpret them with completion/timeout counts. No censored run is discarded. Ratios are reported only when both methods completed all three cases.

| Effect | Method | Mean selected K | Mean true shifts recovered | Precision | Recall | F1 | Mean RMSE | Mean refits | Largest fitted K (mean) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 | path AICc | 100.000 | 33.333 | 0.333 | 0.333 | 0.333 | 0.950 | 56.333 | 100.000 |
| 2 | beam AICc | 79.667 | 9.667 | 0.125 | 0.097 | 0.108 | 1.322 | 211.667 | 100.000 |
| 6 | path AICc | 100.000 | 51.333 | 0.513 | 0.513 | 0.513 | 0.663 | 52.333 | 100.000 |
| 6 | beam AICc | 97.000 | 41.667 | 0.429 | 0.417 | 0.423 | 1.797 | 216.000 | 100.000 |

6 completed path runs contain candidate-path points that did not meet convergence within the configured iteration budget. All retained final covariance-fit modes passed the search completion checks; the candidate generator remains approximate.

These are configured search costs, not times to guarantee recovery of all 100 true shifts. The cap equals truth. Candidate paths may hit their 150-iteration limit; final selected coefficients are unpenalized fits. Per-run candidate/convergence/budget metadata is preserved. Shared regimes, sampling error, multiple traits, null calibration and support bootstrap are excluded.

See [protocol](README.md), [raw metrics](metrics.csv), [summary](summary.json), and [input validation](input-qa.json). This is GeneGalleon Docker validation on a shared host, with no concurrent assistant benchmark or test suite. No SIF claim is made.
