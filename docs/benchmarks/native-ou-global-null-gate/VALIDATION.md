# Validation evidence

- Targeted native bootstrap/path/IC/CLI and CLI contract/wiki suite: **165 passed**.
  Tests cover plus-one p-values, ties, deterministic full-search replay,
  acceptance/rejection, draw resolution, failure propagation, CLI exports,
  fixed-layout/backend rejection, resume configuration changes, and nested support
  repeating the gate with distinct seeds. See `targeted-tests.log`.
- Full NWKIT suite on the final implementation: **3,975 passed, 57 skipped,
  2 failed** in 445.57 seconds. The two failures match the previously established
  baseline: `test_pgls_raw_bootstrap_refits_automatic_gene_parameter` (Pandas string
  dtype assignment) and `test_archived_engine_requires_explicit_scope_and_intact_snapshot`
  (seeded archived replay disagreement). See `full-tests.log` and
  [previous baseline evidence](../native-ou-aic-improvement/validation/README.md).
  The full suite is not clean; these failures were not relaxed or hidden.
- Host static checks: Ruff lint passed; Ruff format checked 477 files; mypy passed
  for all 219 modules; maintainability hard limits passed (existing growth warnings).
- Container measurements: 30 fixed-alpha diagnoses, 48 B=19 paired comparisons,
  3 B=199 representative comparisons, and all 1,509 inner searches completed.
- All 51 pairs' clade metrics and RMSEs independently verified against archived
  AIC predictions and freshly refitted null models. Bootstrap p-values and gate
  decisions recomputed from saved statistics. Source archive hashes verified.

Container commands mounted current NWKIT as `/src` with `PYTHONPATH=/src` and
used the GeneGalleon image in `image-id.txt`. The full suite initially could not
collect without Hypothesis; it was then installed into `/tmp/testdeps` only in
an ephemeral test container, using `PYTHONPATH=/src:/tmp/testdeps`. No production
runtime dependency was added. This is Docker validation; SIF was not available
or tested.

Targeted test command inside the container:

```sh
python -m pytest -q tests/test_shift_native_bootstrap.py tests/test_shift_native_cli.py \
  tests/test_shift_native_path.py tests/test_shift_native_ic.py \
  tests/test_cli.py tests/test_cli_contracts.py tests/test_wiki_examples.py
```

Full suite command: `python -m pytest -q` after ephemeral test dependency setup.
The independent numeric audit command is `python verify_metrics.py` inside the
same container with the benchmark directory structure from the protocol.
