# Worktree integration audit (2026-09-11)

Reviewed all 14 registered GeneGalleon worktrees against main `e202f76`.
Nine directories existed; five of those contained uncommitted changes. Five
registrations pointed to missing directories. Missing working files could not be inspected.
Existing worktrees and recovery branches were preserved.

## Changes retained

- Native OU covariance/alpha configuration, adapter forwarding, provenance,
  resume checks, documentation and the historical 100-tip covariance prototype.
  Trait-specific alpha remains the default; shared alpha is an explicit
  scientific assumption or sensitivity comparison.
- The adapter test now verifies the current small-input dense GLS engine for
  both alpha models, including resume and configuration-mismatch rejection.
- The original scientific-validity review, its GO reproduction script and its
  saved validation log. A historical-status note distinguishes the original
  findings from the current implementation; the old log was not regenerated.

## Changes already incorporated or superseded

| Worktree or recovery commit | Disposition |
|---|---|
| `3751/genegalleon` (GO) | Initial CAFE branch-specific implementation was committed in `8973ae5`, replaced by unmodified CAFE report screening in `c004244`, then hardened in `526df6f`. Do not restore the old likelihood-ratio implementation. |
| `3a0e/genegalleon` (dating) | Documentation and reporting landed in `6102877`; runtime tests were subsequently updated in `d67adb8` and the reconciliation migration. Current tests handle both positive variance and the strict-clock boundary; the old fixture-specific expectation should not replace them. |
| GBIF review `checkout` and `merged` | Main contains the finalized observation-trait implementation (`7706fc8`) and later fixes. A three-way comparison of the merged worktree produced only conflicts in tests referencing the removed kfl1ou stage/dependency; no clean residual changes. The checkout is an earlier version of that work. |
| `a162657` (scheduler bootstrap) | Current entrypoints already search submission, working and script directories, with additional managed-runtime handling. |
| `d83f5ec` (MCMCtree cache) | Current workflow validates cached trees with NWKIT, rejects invalid reuse, and retains explicit rebuild handling. Do not restore the older serialization helper. |
| `4611e87` (RSC output name) | Current producers and consumers already use `rsc_regression` and `.regression.tsv`. |
| `17fcd1b` (symlinked FASTA directory) | Current annotation discovery uses `find -H` and has runtime/static tests for a symlinked search root. The older helper replacement is unnecessary. |

Other existing detached worktrees were clean and their HEADs were ancestors of
main. No worktree cleanup, branch changes, or remote push was performed.
Generated `Rplots.pdf` and `workspace/.gg_cache/` are not source changes.

## Validation

GeneGalleon Docker image `c8c69075b5d5`, with committed NWKIT `0caec89`
mounted through `PYTHONPATH`, one BLAS/OpenMP thread:

- `test_shell_static_safety.py` and `test_gg_util_paths.py`: 343 tests passed.
- Final `test_native_ou_shift.py`: 13 tests passed. The initial run exposed the
  obsolete pruning-engine assertion; this was repaired before the final run.
- Shell syntax and `git diff --check` passed.

The historical benchmark was not rerun for this integration. This was scoped
Docker validation, not the full CI matrix or SIF validation.
