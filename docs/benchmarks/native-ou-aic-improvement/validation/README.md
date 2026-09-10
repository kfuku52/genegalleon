# Validation evidence

- Cross-backend diagnosis: all 18 native likelihoods at R-fitted alpha agree
  within 2e-5. The same fixed-layout likelihood and AIC implementation is retained.
- Baseline replay: all six original selected clade lists, log likelihoods and
  fitted tip means reproduced at absolute tolerance 1e-10 before new timing.
- Final focused suite: 68 passed, covering the new path, CLI routing, support
  replay, old beam search, quick profiles and AIC/BIC/pBIC scoring. This follows
  the final covariance-seed validation, unused-pool budget guard and opt-in
  dispatch decision.
- Full NWKIT suite: 3,964 passed, 57 skipped, two failed. Both are the same
  previously established baseline failures: the PGLS Pandas string-column test
  and the archived null-contract replay audit. See `full-pytest.log` and the
  [previous baseline evidence](../../native-ou-100tips-10shifts-ic/validation/README.md).
  The full suite ran before the last small input/budget guard and default-dispatch changes; the final
  focused suite covers those changes. The full suite is not clean.
- Ruff lint/format and mypy passed. Maintainability hard limits passed with
  pre-existing warnings. No new runtime dependency was added.

The test runtime is the recorded GeneGalleon Docker image with current NWKIT
mounted via PYTHONPATH. Hypothesis was installed only in the ephemeral test
container for full-suite collection; the measurement image was not changed.
The selected candidate and baseline are separately frozen before independent
measurement. This is Docker validation; SIF was not available or tested.

`delivered-routing-replay.json` additionally verifies that delivered auto and
explicit native-path reproduce their frozen baseline/candidate clades, likelihoods
and tip means on one weak and one strong independent dataset (tolerance 1e-10).
