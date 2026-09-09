# Unmodified official IQ-TREE 3 validation

NWKIT requires the ordinary `iqtree3` executable (major version 3 or later) for
the initial model fit. Repeated evaluations use either standard IQ2MC exports
or an optional external worker linked to the unmodified official library.
IQ-TREE itself needs no patch or custom CLI option. NWKIT's adapter owns the
worker protocol and keeps the fitted model, alignment and topology resident.
The Python distribution contains only NWKIT's adapter source, with no IQ-TREE
source, library or executable.

The Docker overlay was built from the unmodified official
[IQ-TREE 3 repository](https://github.com/iqtree/iqtree3), including its recursive
submodules. The runtime reports version 3.1.4. Source revisions and binary
checksums are recorded under `/opt/pg/logs` in `local/genegalleon:iqtree3-dev`.
The default source follows the moving branch in `container/source_branches.env`.

Initial CLI validation in that Docker image:

- 186 NWKIT dating tests passed, covering native models, IQ-TREE models and
  derivatives, input validation, clock inference and uncertainty calculations.
- One additional tiny-codon-branch regression passed against independent
  matrix-exponential pruning, at the standard IQ2MC text-export tolerances.
- Seven GeneGalleon workflow integration tests passed (257.51 seconds),
  including default Gamma and FreeRate profile intervals, fixed species ages,
  and published tree/table/manifest/PDF contracts.
- A two-evaluation benchmark smoke test completed; this is an execution check,
  not a performance comparison or throughput claim.
- The overlay's installed executable resolves to `/opt/pg/iqtree3/iqtree3`;
  its help has no custom `--likelihood-session` interface.

GeneGalleon container source-policy and shell-static tests: 250 passed.
Changed NWKIT Python files passed Ruff checks. Shell syntax and whitespace
checks passed. SIF execution was not available on this macOS host; these results
do not establish SIF compatibility.

Historical custom-session timing reports do not apply to this implementation.

## Optional external library validation

The `local/genegalleon:iqtree-library-dev` overlay builds both the ordinary CLI
and `BUILD_LIB=ON` library from the same unmodified official source snapshot.
It then builds NWKIT's adapter and installs `nwkit-iqtree-worker` separately
under `/usr/local/bin`. Library, adapter and worker identities are recorded in
`/opt/pg/logs/iqtree3_library_worker.json`; the tested library reports 3.1.4.

The focused worker, discovery and IQ-TREE likelihood suite passed 69 tests.
These checks cover numerical agreement with the CLI, repeated evaluations in
one process, model/alignment reuse, independent bootstrap workers, tiny branch
lengths, and failure handling. Two additional GeneGalleon library integration
tests passed, including Gamma and FreeRate profile intervals, fixed species
ages, output artifacts and worker identity in cache provenance.

The final combined GeneGalleon run passed 259 tests (361.62 seconds): all nine
dating integration cases, container source-policy checks and shell-static
safety checks. This includes both explicit CLI and external library execution.

An isolated NWKIT source snapshot containing these changes passed its complete
test suite: 3,005 passed, 5 skipped (572.98 seconds), with 85% total coverage.
The snapshot excludes unrelated concurrent development in the shared checkout.
Ruff, formatting, Mypy, dependency consistency, Bandit, dependency auditing
and the repository's complexity limits also passed.
Wheel/sdist contents and reproducibility checks passed. Installing the generated
wheel in a separate directory preserved its MIT metadata and successfully
discovered the independently installed library worker.

These are Linux ARM64 Docker results. Native macOS, Windows and SIF execution
have not been validated. No new performance comparison is claimed.
See [setup and numerical limits](gene-tree-dating.md) and the NWKIT
[external library build guide](https://github.com/kfuku52/nwkit/wiki/IQ-TREE-library).
