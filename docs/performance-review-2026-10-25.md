# Twelve execution-speed candidates, 2026-10-25

Twelve separate agents evaluated one proposal each. Promising proposals received
an independent review, complete output comparisons and repeated measurements.
Scientific models, thresholds, genetic codes, identifiers, public schemas and
workflow stage order are unchanged. Curated workspace inputs were not used or
modified. An isolated worktree kept unrelated concurrent work out of this change.

## Decisions

| # | Proposal | Decision |
| --- | --- | --- |
| 1 | Build CDS annotation membership maps without exploded DataFrames | Accept. Ordered dictionary insertion retains the first assignment, including missing groups. `setdefault` avoids a second hash lookup. |
| 2 | Index orthogroup member selection | Reject after inspection: normal paths already use indexed maps, sets and `isin`; no repeated full membership scan remained. |
| 3 | Aggregate FIMO annotation columns from one group-position index | Accept. Avoid six sets of per-group Series callbacks; retain first-seen sequence and within-sequence hit order. |
| 4 | Cache normalized donor weights during representative selection | Reject. Five-pair normal and longer-alignment cases regress about 11% and 10%; gains in wider/ineligible cases do not establish a general benefit. |
| 5 | Remove repeated synteny-neighbor scans / accelerate coordinate bounds | Reject. Neighbor lookup is already indexed. A separate NumPy bounds prototype has lower medians but paired ratios cross unity in every profile. |
| 6 | Recognize canonical structural GFF attributes before detailed normalization | Accept. A full match covers only fields that cannot enter metadata repair. Other attributes and source/feature type errors retain the original path. |
| 7 | Avoid repeated transcript ORF input normalization | Reject at screening. The normal actual-CDSKIT call is unchanged in time; no material benefit was established to pursue. The published adapter is unchanged. |
| 8 | Reuse more provenance/input-proof reads | Reject after inspection. Command-local reuse already exists; remaining full rereads detect mutation at the final verification boundary. |
| 9 | Group refinement-review selections once by species | Accept. Complete verified collection avoids one full selection scan per species. Counts retain duplicates; the gallery still uses the last selection. |
| 10 | Share CSUBST scan headers between preflight and database loading | Defer. Direct/live inputs lack a universal immutable-content boundary; unconditional reuse changes mutation/schema error behavior. |
| 11 | Group fixed-track genes once for chromosome ordering | Accept for explicit `target`/`query` modes. Per-chromosome sorting and natural-ID ties remain unchanged. Default `both_length` bypasses this change. |
| 12 | Start indexed FASTA extraction from requested identifiers | Accept for at most 1,024 distinct normalized requests. Larger sets retain the original planner choice; all request sizes and output records remain supported. |

Rejected production prototypes remain outside the published workflow. Their
screening/confirmation data are retained under `tmp/performance-20261025/` in
the review worktree. Proposals 2, 8 and 10 have inspection evidence only.

## Measurement conditions

- Source baseline: `610adefe2d432ca74085496d8db82a9c51f62b50` (0.8.408).
- Freshness-checked GeneGalleon Docker image
  `local/genegalleon:performance-20261025`, immutable ID
  `sha256:4ed6901c185a713a4e73e1c8dd1a835fcf88af3fba88b7e25abede7d64c0b899`;
  runtime-input identity
  `77ee9e8c1f7e8f45ac9863db5757afc0488079b9f9fea199f71cbf84ead121dc`.
- Linux aarch64, Python 3.12.15, pandas 3.0.6, NumPy 1.26.4. Docker Desktop
  provides twelve CPUs and about 7.65 GiB RAM. Benchmark workers and descendants
  use CPU 11; BLAS worker counts are one where applicable.
- One discarded warmup pair, then five alternating fresh-process pairs. Every
  trial is retained. Fixtures, imports and output comparisons are outside the
  operation timers. Timers include the actual affected reads/writes and proof
  checks as described by each benchmark. Interpreter startup is excluded.
- Benchmarks and workflow tests in this review run sequentially after the build.
  Other chats' Docker jobs were observed on this shared host. CPU affinity does
  not provide host-wide isolation; short overlapping controls receive no speed
  claim. Large claims require consistently improving pairs.
- Linux `VmHWM` measures the complete worker before output verification. It can
  include imports, fixture construction and input-fingerprint allocations; in
  particular, chromosome-order fixtures dominate their worker peak. These are
  not whole-workflow memory measurements.

Synthetic inputs and warm caches do not establish production NAS throughput,
whole-analysis speedups or SIF compatibility. Review fixtures use the same storage
for each paired comparison. Membership, GFF normalization and refinement
collection use native container temporary files; FIMO, chromosome ordering and
indexed FASTA fixtures use the mounted checkout.

## Confirmed affected operations

Values are medians. Ranges describe paired before/after wall-time ratios.

| Operation / synthetic input | Before | After | Paired speedup |
| --- | ---: | ---: | ---: |
| Membership loader, 50,000 single-copy rows | 0.03343 s | 0.02605 s | 1.23–1.31× |
| Membership loader, 50,000 rows × four copies | 0.08442 s | 0.06153 s | 1.26–1.45× |
| Membership loader, 1,000 rows × 200 members | 0.06992 s | 0.04870 s | 1.31–1.53× |
| Merge parser/main, 50,000 CDS IDs, OG annotation only | 0.07430 s | 0.06719 s | 1.04–1.38× |
| FIMO branch including loading, 512 sequences × 16 hits | 0.07228 s | 0.05521 s | 1.24–2.10× |
| FIMO, 512 singleton hits | 0.03902 s | 0.03017 s | 1.03–1.90× |
| GFF normalization, 20,000 genes × eight rows, raw | 0.30191 s | 0.17018 s | 1.73–1.80× |
| Same GFF, gzip transport | 0.87921 s | 0.77226 s | 1.02–1.31× |
| GFF with gene metadata, 10,000 genes × eight rows | 0.19791 s | 0.13001 s | 1.38–1.61× |
| Verified review collector, 128 species × 500 loci | 0.48846 s | 0.20570 s | 1.99–2.57× |
| Chromosome ordering, 24 chromosomes × 2,000 genes | 0.39533 s | 0.19052 s | 1.88–2.49× |
| Chromosome ordering, 512 contigs × 128 genes, target | 5.81607 s | 0.23071 s | 19.83–26.65× |
| Same ordering, query mode | 3.87859 s | 0.17625 s | 20.34–23.96× |
| Indexed FASTA, 100,000 records, 32 requested IDs | 0.09885 s | 0.01111 s | 7.85–10.42× |
| Same store, 512 requested IDs | 0.11563 s | 0.03332 s | 3.01–3.94× |
| Same store, 1,024 requested IDs | 0.17885 s | 0.05505 s | 1.93–5.39× |
| 2,048 case-paired requests, 1,024 normalized keys | 0.13349 s | 0.04955 s | 2.64–2.82× |

The [measurement record](benchmarks/workflow-speed-2026-10-25.json) retains every
measured pair's wall time, CPU time and worker peak RSS for accepted and measured
rejected proposals, together with final source hashes and FASTA query plans.

Small, sparse/wide membership and small/dense FIMO controls include overlapping
pairs; no stable improvement is claimed there. GFF's all-annotated control is
0.30359→0.29469 s with paired ratios 0.85–1.08×, so no difference is claimed.
The small four-species collector is 0.00436→0.00776 s, while the 32-species case
is 0.06638→0.04794 s with one slower pair. The accepted collector claim is
therefore scoped to the larger cohort.

The giant-membership fixture reduces worker peak RSS from 106.80 to 99.77 MiB.
FIMO's 512-sequence case increases worker peak RSS by about 0.3 MiB. The review
index retains one extra reference per selection plus list/dictionary overhead;
it trades bounded additional memory for avoiding repeated scans. No general
memory improvement is claimed.

The initial unbounded FASTA join forced indexed lookups even for large request
sets; its 12,500-request median regressed 1.159→1.318 s. The final implementation
counts the already populated TEMP request table, after case normalization and
deduplication, and limits the fixed lookup order to 1,024 keys. This is a query-cost
choice, not a biological limit. Actual `EXPLAIN QUERY PLAN` output matches the
baseline for 1,025, 4,096, 12,500 and 100,000 keys. Larger cases do not receive a
speed claim: the 1,025-key boundary is 0.13881→0.14996 s (ratios 0.87–1.03×),
12,500 keys are 0.47301→0.53354 s (0.82–1.20×), and full extraction is
4.40439→4.22268 s (0.89–1.39×). The extra count adds work, and shared-host I/O
remains variable. Variant queries retain their separate full species scan and
also have no reliable speedup here. FASTA bytes, missing-request codes/logs,
case collisions, lexical identifiers and source/database hashes match.

## Reproduce

Save the baseline without changing branches, then use a matching runtime:

```bash
gg_speed_baseline=tmp/speed-20261025-baseline
mkdir -p "$gg_speed_baseline"
git archive 610adefe2d432ca74085496d8db82a9c51f62b50 workflow/support |
  tar -x -C "$gg_speed_baseline"
export GG_TEST_RUNTIME=docker
export GG_CONTAINER_DOCKER_IMAGE=local/genegalleon:performance-20261025
```

Run comparisons sequentially, without builds or tests. `taskset -c 11` reproduces
this host's affinity; choose an available CPU on another machine. Freshness
checks remain enabled. Rebuilding the runtime resolves the declared moving
branches; the immutable image above identifies this review's environment.

```bash
bash workflow/tests/run_in_runtime.sh taskset -c 11 python workflow/benchmarks/benchmark_orthogroup_annotation_map.py \
  --baseline "$gg_speed_baseline/workflow/support/merge_cds_annotation.py" --repeats 5
bash workflow/tests/run_in_runtime.sh taskset -c 11 python workflow/benchmarks/benchmark_fimo_annotation.py \
  --baseline "$gg_speed_baseline/workflow/support/orthogroup_statistics.py" --repeats 5
bash workflow/tests/run_in_runtime.sh taskset -c 11 python workflow/benchmarks/benchmark_gff_structural_attributes.py \
  --baseline "$gg_speed_baseline/workflow/support/gff_attribute_syntax.py" --repeats 5
bash workflow/tests/run_in_runtime.sh taskset -c 11 python workflow/benchmarks/benchmark_refinement_review.py \
  --baseline "$gg_speed_baseline/workflow/support/plot_gene_model_refinement.py" --repeats 5
bash workflow/tests/run_in_runtime.sh taskset -c 11 python workflow/benchmarks/benchmark_karyotype_order.py \
  --baseline "$gg_speed_baseline/workflow/support/pairwise_synteny.py" --repeats 5
bash workflow/tests/run_in_runtime.sh taskset -c 11 python workflow/benchmarks/benchmark_fasta_query_plan.py \
  --baseline "$gg_speed_baseline/workflow/support/fasta_sequence_store.py" --records 100000 --repeats 5
```

## Validation and delivery

Before delivery, the worktree fast-forwarded to
`1b108ff2f0dce9c623c35bb831fd8fea6fc9c710` (0.8.410), preserving the concurrent
published work. None of the six measured implementation files changed upstream;
their published hashes still match the measured code. The patch version is
0.8.411. The complete final selection below passed **1,182 tests in 123.99 s,
with no skips**, after integration:

```bash
bash workflow/tests/run_in_runtime.sh python -m pytest -q \
  workflow/tests/test_merge_cds_annotation.py \
  workflow/tests/test_orthogroup_statistics_fimo.py \
  workflow/tests/test_format_species_inputs.py \
  workflow/tests/test_paired_species_curation.py \
  workflow/tests/test_gff2genestat.py \
  workflow/tests/test_pairwise_synteny.py \
  workflow/tests/test_plot_gene_model_refinement.py \
  workflow/tests/test_generate_orthogroup_database.py \
  workflow/tests/test_performance_helpers.py \
  workflow/tests/test_gene_evolution_query_mode.py -x
```

The earlier split selections also passed (1,070 and 112 tests). After the lint
requirement for explicit `zip(strict=True)` was applied, the final membership
implementation received a fresh five-pair comparison and its focused tests passed
again. The permanent FASTA benchmark also passed a bounded smoke comparison of
small, 1,025-request and normalized-duplicate profiles against production code.

`GG_LINT_BASH="$PWD/tmp/performance-20261025/bash-parser" bash ./dev lint`
passed: shell syntax, Ruff 0.16.7, eight entrypoints and 27 common parameters.
`git diff --check` and final source/diff review also passed. The host's Bash 3.2
cannot parse all container scripts, so `GG_LINT_BASH` selects an untracked wrapper
delegating Bash syntax checks to the same immutable Docker image with the
checkout mounted read-only. Ruff and the entrypoint configuration checker run
on the host. These are static checks, separate from executed Docker behavior.

The development guide remains the check-selection authority; this section is
an execution record. SIF is unavailable on this macOS host. R helpers and
container inputs are unchanged, so the substantial R/full scientific lanes were
not run. No SIF or whole scientific workflow compatibility claim is made here.
