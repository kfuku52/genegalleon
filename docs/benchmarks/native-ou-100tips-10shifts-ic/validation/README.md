# Validation evidence

The fixed-alpha agreement check compares native and corrected kfl1ou likelihoods
and all three information criteria in 12 root/layout/criterion combinations.
All errors are below 2e-6. The independent optimum-coordinate pBIC probe also
passes all 22 required checks before measurement.

- NWKIT final criterion/CLI tests: 29 passed, including fixed layouts, selection,
  support replay without calibration, changed-criterion resume rejection,
  within-complexity pBIC ranking, zero-alpha eligibility and independent dense
  information determinants. Related search/fit/model tests also passed.
- kfl1ou related pBIC, convergence and multivariate tests: 176 passed, 2 skipped
  because optional genlasso was unavailable. Final AIC tests: 9 passed, including
  public fitting, selection and full-covariance sensitivity/convergence.
- NWKIT full check: 3,940 passed, 57 skipped, **2 failed**. Both failures also
  reproduce on the unmodified HEAD, in the same runtime: the PGLS bootstrap
  Pandas string-column assignment test and the archived null-contract replay
  audit. The full suite is not clean; these failures are independent of this
  implementation. Full and baseline logs are retained here.
- Ruff lint/format, mypy, dependency checks, Bandit and dependency vulnerability
  audit passed. Separate coverage reporting gives 86%, and maintainability
  checks pass their hard limits with existing warnings.

The final frozen image additionally passes the installed-package pBIC and
cross-backend score probes. Installed native source hashes match the source
manifest; kfl1ou installed-library hashes are recorded in the pBIC attestation.
This is GeneGalleon Docker validation. SIF execution was not available.

All six post-run native searches had identical candidate layouts and likelihoods
across the four selection rules. Default-bootstrap branches, tip means and
search counts were also identical to the original runs; new timing measurements
were collected independently. See `native-search-equivalence.json`.
