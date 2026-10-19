# Repository efficiency review and reference annotation measurements

The 2026-10-15 review screened Python support code for repeated tree traversals,
file reads inside loops, whole-table retention and queued parallel work. It also
inspected core/entrypoint orchestration, R/treevis loops, archive/provenance
boundaries, container build stages, test routing and the six parsed GitHub Actions
workflows. Baseline source was GeneGalleon `759b644c`.

| Surface reviewed | Finding and disposition |
| --- | --- |
| WGD anchor placement | A fresh `CladeIndex` and species-tip lookup were built for every Ks anchor. Reuse one command-local `KsBranchPlacer` for an unchanged species tree. |
| WGD duplication classification | Each left-child tip traversed all right-child tips before rejecting different species. Group positioned right-child tips by species once per node and emit the original pair order. |
| Swiss-Prot annotation in genome annotation and gene evolution | All reference descriptions and metadata were retained for each invocation. Retain descriptions for selected accessions and parse metadata in 16,384-row chunks, keeping matching rows. |
| Input staging, refinement, rescue and verification | Recent changes already batch source proofs and preparation. Retain command/publication boundaries and input-change detection; see [batch I/O measurements](batch-io-performance.md#runtime-and-genome-index-verification-boundaries). |
| Summaries, FASTA/GFF and plotting | Existing indexed/batched paths cover many previous hotspots; see [batch I/O measurements](batch-io-performance.md). The two-pass soft-mask reporting helper has no production call site, so it was not changed. |
| Archive storage and database generation | Preserve transaction, verification and schema contracts; previous bounded I/O work is documented in [remaining I/O performance](remaining-io-performance.md). |
| Container, validation and Actions | Existing source-artifact/runtime caches, wheel sharing, shards, superseded-run cancellation and runtime-change gating are retained. The parsed workflows use hosted runners. No CI coverage or publishing behavior was changed. |

This source review and the bounded measurements below do not profile every
scientific program or establish that all remaining inefficiencies are removed.

## Behavior preserved

`branch_for_ks(focal, ks, tree, bounds)` remains available for individual calls.
The batch placer reads the current bounds on every call, including uncertainty,
monotonicity and missing-data checks. A new batch owns its own index; no persistent
cache is introduced. Callers must keep the indexed tree unchanged within a batch.

Pair selection retains left-tip order, then matching right-tip order, and excludes
unmapped tips. It still emits every eligible same-species pair. If every tip is
from one species, the required output itself remains quadratic.

Swiss-Prot descriptions retain the last duplicate accession; metadata retains the
first. The existing unfiltered Python reader calls remain supported. Both plain
and gzip references preserve literal identifiers such as `NA` and `001`, quoted
tabs/newlines, blank fields and the output column schema. Complete reference files
are still parsed; this change reduces retained objects rather than disk-read
volume. Metadata memory grows with one chunk plus matching rows, so a dense match
set provides less memory saving.

## Comparable measurements

The benchmark uses fresh child processes, one warmup per implementation and three
measured trials in alternating baseline/current order, without concurrent tests
or builds. Fixtures and imports are outside the timed region. The Ks timer
includes construction of the reusable placer, and the Swiss-Prot timer executes
the complete annotation `main()` through TSV writing. Peak RSS is for the whole
worker process, including imports and fixture preparation. Output hashing happens
after timing. Every warmup and measured output hash matched across implementations.

Environment: Docker Linux/aarch64, Python 3.12.15, pandas 3.0.6, NWKIT 0.43.52;
the runtime matched the daily owned upstream snapshot and current container
inputs. The immutable image ID, source hashes, individual trials and exact output
hashes are in [the measurement record](benchmarks/reference-annotation-2026-10-15.json).

| Synthetic workload | Before median | After median | Ratio | Before / after process peak RSS |
| --- | ---: | ---: | ---: | ---: |
| 4,000 Ks placements on a balanced 128-species tree | 0.9932 s | 0.00846 s | 117.5× | 109.9 / 109.6 MiB |
| Eight cross-child selections, 512 tips per child, 128 species, four copies per species | 1.7438 s | 0.00770 s | 226.3× | 118.4 / 118.5 MiB |
| Complete annotation: 200,000 references, 1,000 queries, 146.2 MB reference input | 2.3407 s | 1.4408 s | 1.62× | 582.5 / 145.6 MiB |

Index construction drops from 4,000 calls to one. Pair selection preserves all
16,320 output pairs across the eight selections. Swiss-Prot annotation reduces
median process peak RSS by 75.0% and elapsed time by 38.4% on this sparse workload.
The first two rows measure the affected operations rather than complete WGD
analysis, reconciliation, alignment or scientific model fitting.

A dense control uses 20,000 references and 20,000 queries. Before/after elapsed
time ranges overlap (0.523–0.629 s and 0.474–0.832 s); medians are 0.546 and
0.654 s, with 133.9 and 128.9 MiB peak RSS. No speedup is established for that
case. Chunk filtering has overhead when most references are required. These
warm-cache synthetic workloads do not establish whole-workflow, NAS, GPU or SIF
performance.

## Reproduce

Save the baseline support tree without switching branches, then run the benchmark
in the same GeneGalleon runtime as the current implementation:

```bash
baseline="$(mktemp -d /tmp/gg-reference-before.XXXXXX)"
git archive 759b644c workflow/support | tar -x -C "${baseline}"
GENEGALLEON_DOCKER_EXTRA_BINDS="${baseline}" \
  bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_reference_annotation.py \
  --baseline-support "${baseline}/workflow/support" --output tmp/reference-comparison.json
```

For the dense control, use the same command with `--cases uniprot --references
20000 --queries 20000` and a different output path. Defaults are 128 species,
4,000 anchors, four copies, eight pair selections, 200,000 references, 1,000
queries and three measured repetitions. `--baseline-only` allows recording an
unchanged implementation before editing. Fixtures use temporary directories.

## Validation

The following focused checks use the Docker runtime; none establishes SIF
compatibility:

```bash
bash workflow/tests/run_in_runtime.sh python -m pytest -q \
  workflow/tests/test_reformat_uniprot_diamond.py workflow/tests/test_wgd_evidence.py \
  workflow/tests/test_wgd_ssd_contracts.py -x --durations=10
bash ./dev check fast workflow/tests/test_merge_cds_annotation.py \
  workflow/tests/test_build_uniprot_metadata_from_dat.py -x
bash workflow/tests/run_in_runtime.sh python -m pytest -q \
  workflow/tests/test_wgd_ssd_runtime.py workflow/tests/test_wgd_ssd_gene_stage.py -x --durations=10
```

The commands passed 191, 28 and 7 tests respectively (226 total, no skips).
The final command exercises real count/synteny/Ks tools, classification and the
gene-evolution cache/archive integration; it completed in 198.22 seconds.

Host `ruff check workflow container/scripts` and
`python3 workflow/support/entrypoint_config_schema.py --check` passed.
`bash ./dev lint` was blocked by macOS Bash 3.2; all 75 tracked shell scripts and
`dev` passed syntax checks separately with container Bash. `git diff --check`,
benchmark source/output hashes and report links also passed. No full scientific
analysis or SIF run is part of these measurements.
