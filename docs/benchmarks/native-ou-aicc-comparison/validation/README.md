# Validation

- Fixed-alpha AIC/AICc scores agree between NWKIT and kfl1ou over eight
  root/shared-regime combinations. Maximum AICc score error is 1.42e-14;
  acceptance tolerance is 2e-6. The corrected-pBIC preflight also passed 22 checks.
- Focused native IC/path/CLI tests: 63 passed before the final added CLI dispatch
  and invalid-null cases. The full suite below includes those final cases.
- Final full NWKIT suite: **3,986 passed, 57 skipped, two failed** in 462.20 s.
  Both failures match the established baseline: PGLS raw bootstrap string-dtype
  assignment and the archived null-contract seeded replay. See `full-tests.log`
  and [previous validation](../../native-ou-global-null-gate/VALIDATION.md).
- Ruff lint passed, Ruff formatting passed, mypy passed all 219 modules;
  maintainability hard limits passed with previously present growth warnings.
- Tests and cross-backend checks used GeneGalleon Docker. Hypothesis was installed
  only into `/tmp/testdeps` in the ephemeral full-suite test container. No
  production dependency or existing default was changed. SIF was not tested.

The source archive is verified against its manifest. Per-run score/metric/input
checks are recorded in the study's `independent-qa.json`. The fresh kfl1ou AIC
replay uses the unchanged adapter and verifies old selected predictions, likelihoods
and accuracy at numerical tolerance, without mixing criteria.
