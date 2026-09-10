# Analytical-only CSUBST scan integration — 2026-09-10

GeneGalleon now runs CSUBST scan with calibration disabled and zero
permutations. It computes global Benjamini–Hochberg FDR from
`p_rate_enrichment_asymptotic`, pooling all imported candidate rows across
orthogroups, traits and match classes. No empirical/maxT/bootstrap P value is
computed by the workflow. This supersedes the earlier resampling integration.

The result column is `q_rate_enrichment_asymptotic_global`.
`aa_change_fdr_metadata` records the method, columns, correction scope, finite
test count and undefined count. Missing values remain missing and invalid
probabilities cause failure. Empty scans retain their schemas and zero-test
metadata. Support views preserve the original P and global FDR; they do not
recompute BH. Candidate reports default to the global FDR column.

Current CSUBST inference metadata and an atomic per-family audit ZIP are
retained. The ZIP includes input copies and the disabled-calibration JSON
record. Old scan/database/summary contracts are invalidated by a new
analytical-BH provenance parameter. Legacy schemas are rejected before database
replacement. Migration is described in
[configuration and common parameters](../configuration-and-common-parameters.md#csubst-scan).

## Validation

Runtime: `local/genegalleon:csubst-scan-inference-dev`, a local GeneGalleon
Docker image with CSUBST `75b68b7` compiled from its clean checkout (package
version 1.15.0). Python 3.12.14 and IQ-TREE 3.1.4, Linux ARM64. These identifiers
record validation inputs; they are not upstream defaults.

The focused suite covers:

- Current CSUBST analytical output through DB import, BH and candidate selection.
- The actual core command, with precomputed IQ-TREE files: calibration is `none`,
  requested replicates are zero, analytical P values exist and all empirical
  P columns are undefined.
- Empty scans and zero-test metadata, global BH across OG/trait/match boundaries,
  ties, missing values, invalid probabilities and correction denominator.
- Unchanged FDR under support filtering, report ZIP/cache behavior and shell
  configuration/provenance contracts.

```bash
docker run --rm -v "$PWD:/work:ro" -w /work \
  -e PYTHONDONTWRITEBYTECODE=1 local/genegalleon:csubst-scan-inference-dev \
  python -m pytest -q -p no:cacheprovider \
  workflow/tests/test_generate_orthogroup_database.py \
  workflow/tests/test_plot_csubst_aa_change_summary.py \
  workflow/tests/test_csubst_scan_candidate_sites.py \
  workflow/tests/test_csubst_scan_runtime_integration.py \
  workflow/tests/test_gg_util_paths.py \
  workflow/tests/test_shell_static_safety.py
```

**445 passed, no skips (36.02 seconds)** on an isolated snapshot of the staged
commit contents, excluding unrelated worktree changes. Host shell syntax and
Ruff checks also passed. No whole-genome
end-to-end analysis, full repository suite or SIF execution was performed.
Docker validation does not establish SIF compatibility.

## Analytical P limitations

BH is implemented as requested, but its output is a nominal FDR estimate whose
validity depends on the input P values and selection/dependence assumptions.
Fractional posterior event mass is not an independent integer Poisson sample;
small-sample asymptotic tails and candidate selection using the same data need
additional modeling or calibration in CSUBST.

These problems are addressable, but not by changing the BH formula or rounding
posterior mass. A dependency-side solution needs a justified observation/
latent-event model, treatment of ancestral-state uncertainty and exposure,
and either a prespecified hypothesis family or selection-aware inference.
An analytically evaluated, valid model can support very small P values without
per-OG permutation P calculations. Independent simulation is still needed to
check its approximation and complete decision procedure; a passing execution
suite is not that scientific validation.
