# Workflow helper efficiency

The follow-up audit addresses twelve repeated-work and allocation costs. It
preserves biological models, resampling streams, solver search limits, literal
identifiers and output schemas.

| Area | Change |
| --- | --- |
| Promoters | One streaming SeqKit BED extraction per species; restore request order and preserve duplicate contigs, gzip and strand/edge handling. |
| Cluster bootstrap | Bound the temporary draw matrix while retaining exact RNG state and replicate sums. |
| Synteny geometry | Group genes by scaffold once. |
| Joint synteny ordering | Reuse geometry, ID lookups and parsed blocks within the invocation; retain every track solve and tie rule. |
| Expression tables | Share typed, content-bound SQLite row indexes across family queries. |
| Annotation counts | Skip unused columns while parsing; retain original NA/comment/quote semantics. |
| Gene-tree rooting | Reuse exact rooted trees from NWKIT's comparison pass instead of repeating MAD/midpoint; singleton-bearing inputs retain native tree shape. |
| K-mer distances | Store only the used upper triangle, retaining pair order and calculations. |
| HGT context | Index scaffold start/end orders and overlapping intervals, including nested loci. |
| Resolved CDS traits | Query selected rows through the same table index, retaining input and selected-sequence hashes. |
| Dated-tree tip positions | Build a name-to-rank lookup once. |
| BUSCO summary | Collect species tables and bind once in file order. |

The shared table index is optional for direct helpers and enabled by family
workflows at `workspace/.gg_cache/trait_tables`. Each query freshly hashes the
source and database, then fences their identities through retrieval. A shared-filesystem namespace lock
serializes publication/retrieval for each content key; a partial or corrupt
publication is rebuilt. Interrupted lock owners require the same explicit
reconciliation as other shared-workspace namespace locks; ownership is never
stolen based on age. Parser source, indexer source, pandas version and string
inference settings bind the key. SQLite stores literal IDs, source row indices,
column types and duplicate rows without executing serialized Python objects.
Wide rows are compressed independently to reduce index I/O while retaining
random access.
Cold construction adds latency and disk use; reuse is intended for many families.
For internal singleton nodes, rooting still computes native trees separately
because NWKIT normalizes those nodes in comparison output; those inputs retain
the original tree shape and branch placement.

## Reproduce

Save the baseline support directory from commit `3ce51d28`, then run:

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_workflow_efficiency.py \
  --baseline-support /path/to/baseline/workflow/support \
  --repeats 3 --output /path/to/comparison.json
```

The runner alternates baseline/current in fresh processes, discards one warmup
per variant, compares output hashes and records all trials, wall time and Linux
peak RSS. Fixture creation and imports are excluded; cold table-index creation,
fresh hashes and retrieval are included. Rooting includes its function-local
imports and comparison rendering. PDF metadata dates are normalized only for
that equivalence comparison. The k-mer compiler is excluded from timing and its
RSS is separated from the measured kernel.

Interpret each workload at its measured boundary: R cases cover preprocessing,
not plotting; dated-tree timing covers tip positioning; the bootstrap null uses
100 draws to isolate bootstrap allocations. Empty k-mer sketches isolate pair
storage and output, with real nonempty sketches covered by runtime regressions.
The resolved-CDS workload covers 200 families and a verified 40,000-row table;
its external-source map is empty to isolate table retrieval. Real FASTA/GFF
source verification remains required. The expression workload covers 20 families
and 80 traits. Short runs can cost more because they pay cold construction.
Both variants use the same runtime and source fingerprints are recorded.
These are helper measurements, not end-to-end analysis speedups. Docker results
do not establish SIF compatibility.

## Measured results

Linux ARM64 Docker, Python 3.12.15, NumPy 1.26.4, pandas 3.0.6, NWKIT 0.43.53,
SeqKit 2.14.0, R 4.4.3 and GCC 13.3.0. Values are medians of three measured
trials; the [complete record](benchmarks/workflow-efficiency-2026-10-15.json)
contains all samples, output hashes and source fingerprints.

| Workload | Baseline seconds | Updated seconds | Ratio |
| --- | ---: | ---: | ---: |
| Promoters (60 genes, 6.4 MB genome) | 3.141 | 0.031 | 100.7× |
| Bootstrap (10,000 blocks, 2,000 draws) | 0.289 | 0.125 | 2.3× |
| Synteny geometry (40,000 genes, 512 scaffolds) | 0.571 | 0.013 | 44.4× |
| Joint ordering (20,000 genes/track, 121 solves) | 3.401 | 0.112 | 30.5× |
| Expression (20 families, 80 traits) | 7.117 | 3.272 | 2.2× |
| Resolved CDS (200 families) | 11.505 | 3.846 | 3.0× |
| HGT coordinates (12,000 genes, 520 events) | 6.633 | 0.324 | 20.5× |
| Tip positioning (12,000 tips) | 0.561 | 0.005 | 113.7× |
| K-mer pair storage/output (1,600 species) | 1.019 | 0.957 | overlapping ranges |
| Root selection and comparison (4 tips) | 4.479 | 2.185 | 2.0× |
| R annotation counts (3 × 30,000 rows, 80 columns) | 13.568 | 8.779 | 1.5× |
| R BUSCO preprocessing (400 × 1,000 records) | 5.793 | 1.261 | 4.6× |

The k-mer wall-time ranges overlap; no time improvement is claimed. Its
measured kernel peak RSS falls from **100.9 to 51.9 MiB**, with identical bytes.
At the 10,000-species limit, the result array itself changes from 4,000,000,000
to 1,999,800,000 bytes; this is an allocation calculation, not a measured full run.

Bootstrap peak RSS falls from **349.9 to 60.3 MiB**. R annotation preprocessing
falls from 205.8 to 147.6 MiB, and BUSCO preprocessing from 174.6 to 129.6 MiB.
Other workloads show similar worker RSS, except rooting: drawing now runs in
the parent, whose RSS changes from 109.5 to 188.3 MiB. The baseline rendering
child already peaked at 187.2 MiB. These are per-process peaks; they cannot be
added to infer whole-job memory, and no rooting memory improvement is claimed.

Cold construction remains a tradeoff: the first resolved-CDS family takes
about 0.290 seconds versus 0.060 without the index; reuse over 200 families
reduces total time. For expression, the first family costs about 1.948 seconds versus 0.350;
reuse over 20 families offsets that build cost. An unrelated four-second website
job affected the original expression pass; the entire case was repeated in a
quiet window. The record retains the excluded original samples and reason.

## Validation

Implementation validation used the Docker image recorded above, with its built-in NWKIT 0.43.53;
no source override was used in the final checks. `GG_RUNTIME_FRESHNESS=always`
matched the image to current inputs before execution.

- `bash ./dev check runtime-python` on promoter, pairwise synteny, NWKIT tree workflows, BUSCO guide-tree and workflow CDS validation files: **277 passed**.
- `bash ./dev check fast` on the affected promoter, trait, CDS, GFF, subgenome, HGT-coordinate and reconciliation-rooting files: **238 passed**.
- Additional table-index/CDS/namespace-lock selection: **76 passed**; covers concurrent builders, corruption, source replacement, compression aliases and literal values/types.
- Outgoing-commit checks for container-build selection, artifact provenance/audit, rescue stores and raw validation: **310 passed, 1 skipped**. The skip is the optional zstd codec because `zstandard` is absent; gzip/auto paths passed.
- `bash ./dev lint`: passed, using the documented `GG_LINT_BASH` Docker wrapper for Bash 4+; ShellCheck on tracked shell scripts and 162 shell-static tests also passed.
- Runtime R checks: 30 source parses, actual BUSCO summary/annotation-count regressions and dated-tree interval/NHX regressions passed.
- Follow-up CI collection check: all 11 promoter runtime cases are required; validation-runner/development-tooling checks: **65 passed**.
- Follow-up legacy-input routing fixture: supply the core's workspace variable
  and check cache arguments, including spaces in paths.
  `GG_RUNTIME_FRESHNESS=always bash ./dev check runtime-python workflow/tests/test_gene_model_refinement_downstream.py -x`:
  **30 passed** in rebuilt Docker image
  `sha256:a28e1760f1cf2cca40798943bc8e48d594b160147c00baf0441820ada07e41a7`,
  runtime input `f0f256153dba58726e3cb60f15f81002adead91a9c458d8fc6da75bad40d8d31`.
  The old image was rejected as stale; it was used only to reproduce the
  isolated fixture failure before rebuilding.
- `git diff --check`: passed.

NWKIT's owning change (`88d3cba6`, version 0.43.53) passed its `tools/check.py
release` gate on macOS ARM64/Python 3.14: **5,167 passed, 127 skipped**, 86%
combined line/branch coverage, security checks, and reproducible wheel/sdist
checks. Optional studies/backends are separate from the rooting coverage.
Local Docker evidence does not establish SIF behavior; published-runtime/SIF
qualification is handled by the repository's delivery workflows.
