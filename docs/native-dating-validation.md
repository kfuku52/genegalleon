# Native dating validation (2026-09-10)

The default native GY94/F3x4 + G4 dating path remains exploratory. Its nominal
95% profile intervals do not consistently attain 95% coverage in independent
sequence simulations. Species ages are fixed throughout this study.

## Changes

GeneGalleon preserves the actual estimator, sequence model, uncertainty status,
nominal level and diagnostics in orthogroup statistics. The interpretation field
identifies experimental conditional estimates. NWKIT plots state this scope and
that substitution parameters are fitted. Missing and calibration-limited
intervals remain explicit.

An independent internal-duplication family exposed a numerical initialization
failure in NWKIT: a branch-only starting estimate collapsed a duration to its
constraint boundary. The exact sequence optimizer then reported incompatible
constraints despite a feasible chronology. The fix uses the existing interior
chronology as the sequence starting point when that boundary condition occurs.
It does not alter the objective, constraints, minimum duration or rate-SD fit.
The diagnostic is `sequence_initial_ages_reset_from_duration_boundary`.

## Evidence

The primary pre-fix study contains 400 independent AliSim families, each with
one prespecified duplication. A further 250 stress families test longer
alignments, larger families, combined loss/copy-rate shifts, incorrect mappings,
and low rate variation. The NWKIT evidence directory
`examples/radte/default-profile-validation/` contains frozen protocols, per-family
results, source hashes and an independently audited summary and coverage plot.

| Primary condition | Trials | Points | Returned intervals | Truth covered / returned |
|---|---:|---:|---:|---:|
| Root duplication | 200 | 200 | 189 | 176/189 (93.1%) |
| High rate variation | 100 | 100 | 98 | 85/98 (86.7%) |
| Internal duplication | 100 | 99 | 99 | 95/99 (96.0%) |

All 650 input/result records passed an independent audit. In stress tests,
longer alignments covered 41/49 returned intervals (83.7%); larger families
covered 48/48 with two profile timeouts; loss plus copy-rate shift covered
49/50. Wrong mapping covered 41/43 returned intervals, all calibration-limited,
with seven strict-clock limits. Low rate variation covered 35/37 returned
intervals, with 13 strict-clock limits (35/50 correct returns overall).

Failures and missing intervals also remain in all-trial denominators. These
conditional fractions must not be interpreted as the chance of returning a
correct interval on any attempted family. Calibration-limited intervals are
included among returned intervals and separately flagged.

The original failed family remains a failure in the pre-fix study. Its retained
regression fixture fails before the initialization fix and succeeds afterward;
independent restarts agree on the unchanged objective. Its repaired age is 9.341
with interval [7.388, 10], against truth 7.5. This one repair is not a coverage
validation. Twenty additional held-out internal families on the fixed code
returned 20 intervals, 19 covering truth; none triggered the initialization reset.
They are a separate regression check, not part of the pre-fix coverage study.

A paired diagnostic supplies the true rate SD to 50 existing root families:
coverage changes from 43/48 returned intervals to 47/50. Auto can also select a
different estimator, so this oracle comparison does not identify a unique cause
or provide a deployable calibration fix.

## Scope and verification

Simulations supply true topology, use LCA reconciliation and externally generated
codon sequences, and refit the substitution model. They do not validate topology
inference, ILS, empirical data or uncertainty in species ages. No interval
inflation or method selection was fitted to these outcomes.

Validation uses Docker image `local/genegalleon:standard-iqtree-dev` with local
NWKIT sources mounted on `PYTHONPATH`. It does not certify the installed image
snapshot, a new release or SIF execution. NWKIT owns the numerical implementation
and detailed reproducibility record; GeneGalleon owns workflow/report integration.

Verification completed: 247 GeneGalleon summary/static checks, 119 NWKIT
numerical/profile regression checks, and nine GeneGalleon native/IQ-TREE
integration scenarios. The repaired example plot was also rendered and inspected.

## Final commit review

The independent commit snapshot passed 354 NWKIT dating, plot, CLI and example
tests and 256 GeneGalleon dating integration, summary and shell safety tests.
Five additional smoke tests passed. NWKIT Ruff/format, mypy, security/static
checks and complexity limits passed. Wheel/source distributions passed content
and reproducibility checks. These are scoped local checks, not the full hosted
CI matrix; SIF execution and cross-platform runtime validation were not run.

The review also fixed a validation-runner failure: nonfinite returned ages now
remain failed observations, and an invalid profile result preserves a completed
point-only fit. The original 650-family results are unchanged. The initialization
boundary check was extracted into a helper to meet the complexity limit without
changing its decision or the fitted objective.
