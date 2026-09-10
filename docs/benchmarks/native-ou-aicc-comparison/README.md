# AICc addendum: kfl1ou and NWKIT native shift search

This adds AICc to the original seven-method comparison and includes the newer
NWKIT native-path search. See [comparison tables](RESULTS.md).

## Prespecified measurements

- Same six nonnull fixtures from `../native-ou-100tips-10shifts-ic/data`: 100 tips,
  ten planted shifts, effects 2/6, three paired seeds per effect. These six fixtures were used during native-path
  development; they are not held-out tests.
- Five fresh procedures per fixture: NWKIT beam AICc, NWKIT path AIC, NWKIT path
  AICc, kfl1ou AICc, and kfl1ou AIC replay. The other seven-method accuracy values
  remain frozen. The replay checks the previous kfl1ou AIC outputs.
- Thirty null fixtures from `../native-ou-aic-improvement/data/effect0-*`:
  NWKIT beam AICc, NWKIT path AICc and kfl1ou AICc. Ordinary AIC null outcomes
  are taken from that previous study with exactly the same inputs.
- 120 new runs, sequential, with cyclic order rotation, one BLAS/OpenMP thread,
  same GeneGalleon Docker image, and a 1,800-second timeout per run. Whole-suite
  validation may run concurrently; timing is not a controlled comparison against
  the original frozen batch. The comparison table emphasizes recovery and false
  detections rather than mixed-batch runtime ratios.

No truth, input, search cap, alpha bound or penalty was tuned from these results.
All models estimate alpha and process variance. Native cap 10, pool 24, refit
budget 48, screening budget 2000, beam width 2, lasso iterations 150; kfl1ou
uses the same adapter and default search as the prior comparison. The native
second path updates covariance using the requested AIC/AICc winner; its candidate
set can therefore differ by criterion. AICc is not applied retrospectively to
only the AIC winner.

## AICc definition and implementation check

NWKIT now supports AICc for fixed layouts, exhaustive/beam search, and native-path
search without convergence (beam/exhaustive retain their convergence support).
It adds `2p(p+1)/(n-p-1)` to AIC, using kfl1ou's independent-trait convention:
shared locations count once in p; n sums observed scalar coordinates across
traits. Explicitly fixed covariance parameters are excluded. Nonpositive
n-p-1 produces an ineligible candidate with an explicit status. The global-null
bootstrap gate remains AIC-only. No default was changed.

Before timing, the existing corrected-pBIC probe passed all 22 checks and AIC/
AICc fixed-alpha scores agreed between NWKIT and kfl1ou in eight combinations
of root treatment and shared/distinct regimes, within 2e-6. See `validation/`.
Unit tests independently cover missing observations, shared parameter counting,
fixed versus estimated covariance, finite-sample boundaries, selection and CLI
support/resume routing.

## Reproduction and interpretation

Run `validate.py`, then `run.py`, inside the recorded image with NWKIT mounted at
`/src`, this repository at `/gg`, and `PYTHONPATH=/src`. All BLAS/OpenMP thread
variables are set to one. `run.py` refuses to overwrite its results file. Use
an empty output copy with the same sibling fixture directories for a rerun.
`source.tar.gz` and `source-manifest.json` preserve the NWKIT implementation;
`implementation.patch` is relative to the previous global-null-gate source archive.
Per-run input hashes, scores, selected clades, predictions, CPU/RSS and failures
are retained in raw JSON. Run `summarize.py` to regenerate both tables.

Precision/recall/F1 use exact descendant-tip clades, with empty precision zero;
RMSE compares predicted tip means against independently simulated true means.
Accuracy values average over datasets; FPR is the proportion of null datasets
with any selected shift. Three nonnull seeds per effect and thirty null datasets
are development evidence, not proof of population reliability or adoption.
The ten-shift cap equals truth in the nonnull comparison. No measurement error,
missingness, convergence or multivariate performance is assessed here. Docker
validation does not establish SIF compatibility.

After summarizing, `verify.py` independently checks every fresh score against
`-2 logLik + 2(2K+3)` plus the AICc correction when requested (all benchmark
fits have one fully observed trait, estimated alpha/variance and distinct regimes).
It also rechecks input hashes, selected-clade metrics, RMSE and source archives.
