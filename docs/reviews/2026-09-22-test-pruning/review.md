# Test pruning review — 2026-09-22

The review covered the 177 Python test modules, 15 R test scripts, shared
fixtures/lane selection and validation commands. Historical test logs and
scientific benchmark evidence under `docs/` are records, not active test suites.
No production behavior or validation skip policy is changed by this cleanup.

The decision is the plausible defect a test detects, weighed against its
maintenance cost, process startup and dependence on source spelling. Neither
case counts nor coverage percentages are retention targets.

## Removed and consolidated

| Area | Decision and remaining evidence |
| --- | --- |
| Shell source snapshots | Remove exact quoted/unquoted command copies. The existing delivery ShellCheck scans every tracked shell script; actual path, publication and FASTQ tests remain. Keep project-specific safety checks that ShellCheck does not supply. |
| Shell helper implementation checks | Remove duplicate checks of copy/move, private pycache, environment forwarding, version-report reuse, lock and semaphore internals. `test_gg_util_paths.py`, `test_gg_shared_lock_shell.py`, `test_shared_lock.py` and namespace-lock tests execute failures, contention and publication recovery. |
| Genome evolution wiring snapshots | Remove duplicates of protein input selection, genetic codes, legacy defaults, required species trees, storage modes and core sampling. `test_genome_evolution_protein_mode.py` executes these workflows. Keep wiring assertions for stages without an equivalent execution check. |
| Native analysis wiring | Remove copied OU/RSC/root-extraction command strings; their runtime tests execute the core stage and check output bundles, tree routing and failed publication. |
| Module organization and style | Remove facade-import reconstruction, module line-count limits, wildcard-import bans, declaration/comment ordering and error-message quoting snapshots. These do not define usable behavior. The explicit self-contained-core architecture policy check remains. |
| Plot/summary source snapshots | Delete `test_r_plot_static.py` and `test_transcriptome_summary_static.py`. Literal warning text, palette strings and source fragments do not establish numerical filtering or rendered output correctness. Real R plot tests and multifamily PDF/SVG tests remain. There is no claim that they replace every removed transcriptome-summary semantic check: those checks never executed the algorithm. |
| Copied numerical implementation | Delete `test_aligned_taxid_fallback_logic_matches_expected`: it only ran copied NumPy expressions and could not detect a production change. Production lineage-resolution tests remain. |
| CLI help | Restrict discovery to 105 actual entrypoints rather than executing all 137 Python files. Remove fake kftools imports for two commands with dedicated tree/statistics tests; they are not fast-lane CLI smoke targets. Retained smoke tests require real usage output and no unexpected files. |
| BLAST coverage CLI | Merge the `--ncpu 2` case into the existing output-column/value test; serial/parallel equivalence remains separately tested. |
| FASTA discovery | Combine suffix acceptance, hidden/non-file exclusion and deterministic order into one fixture deliberately created in reverse order. |
| R interval merging | Fold the simple reversed interval into the wider bridging regression, retaining mixed strands and adding exact normalized start/end assertions. Preserve separate disjoint/transitive cases. |

The Python reduction is 134 collected cases (100 removed functions, two fewer
file-discovery cases and 32 fewer help cases). One R scenario is consolidated.
This is an accounting of changes, not a performance measurement.

## Retention decisions across the rest of the suite

- **Archive/provenance/workflow API:** retain same-size/same-mtime corruption,
  interrupted transactions, inode exhaustion, symbolic/hard-link aliasing,
  lock ownership, concurrent readers/writers, relocation and optional-output
  tests. These protect existing user results; superficially similar error cases
  enter different publication or recovery paths.
- **Download/providers/FASTQ:** retain loopback HTTP and provider payload tests,
  ETag/range changes, resumed gzip members, partial reads, redirects, credentials,
  rate limits and failed-cache publication. Different providers have different
  schemas. Controlled remote responses exercise local decisions rather than
  merely asserting that a mock returns its fixture.
- **GFF/CDS/species identity:** retain isoform selection, strand/phase handling,
  trans-splicing, ambiguous aliases, qualified species labels and source binding.
  Collapsing these cases could silently assign sequences to the wrong gene or
  species.
- **Scientific calculations:** retain analytical references and real upstream
  integration for reconciliation, dating, covariance, copy-number selection,
  missingness and multiplicity. Empty/invalid/unknown values are biologically
  different from zero. Do not remove a numerical oracle just because it uses
  the same mathematical identity as the implementation.
- **R and rendered figures:** retain data-to-label/color mapping, interval units,
  clipping, fixed physical widths, legends and actual driver output tests. The
  layout and normalization code belongs to GeneGalleon, not to ggplot2 alone.
- **Build/runtime/scheduler/CI:** retain source policy, artifact identity,
  architecture/loader compatibility, source-resolution failure, resource
  budgets, signal propagation and scheduler directives. Runtime integration
  cannot by itself guarantee these deployment contracts.
- **Templates, summaries and public CLI:** retain checked-in input validity,
  schema/column contracts, identifier boundaries, missing-data behavior and
  actual output checks. Simple-looking assertions on these externally consumed
  formats are valuable compatibility checks.

A per-function deletion ledger is in [removed-tests.tsv](removed-tests.tsv).

## Validation

Validation used a freshly built ARM64 GeneGalleon development Docker image,
`local/genegalleon:test-pruning`, with the normal runtime freshness gate enabled.
SIF compatibility was not tested on this macOS host.

- The full Python invocation begun before pruning collected the original suite:
  2,820 passed; two real 3Di predictor cases failed because the host UID could not
  write `/.cache`. These cases were retained and rechecked using writable
  `CSUBST_CACHE_DIR` and `HF_HOME` locations: both passed.
- Consolidated BLAST/file-discovery tests and CLI smoke checks: 132 passed before
  dropping the two kftools-only smoke targets; the final CLI selection was
  rechecked separately: all 105 passed.
- Static lane: 292 passed; the final affected static module was rechecked after
  its last redundant assertion/helper removal: 162 passed.
- All 16 R validation commands passed, including the consolidated interval case.
- ShellCheck over all tracked shell scripts, actionlint, composite-action lint,
  Ruff, entrypoint configuration schema verification and container Bash syntax
  checks passed. The host's Bash 3.2 cannot run the Bash 4+ syntax lane; its
  constituent checks were run with the appropriate tools instead.

The shared worktree also contained an independent repository-audit task's
changes. Those production and fixture changes are outside this test-pruning
commit; its validation results should not be read as a clean-checkout benchmark.
