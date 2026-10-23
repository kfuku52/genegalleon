# Workflow helper performance review, 2026-10-23

Ten candidates were assigned to ten separate agents; seven were accepted and
three rejected. Measurements and runtime tests were serialized within this
review. Unrelated host/container jobs were active, so short control differences
are treated as noise when paired results
overlap or change direction. No scientific thresholds, models, identifiers,
output schemas, archive layouts or entrypoint settings were changed.

## Environment and comparison method

- Baseline: `a9f19346771e540a4a669a6470981dff91ec4fdc` (checkout version 0.8.380).
- GeneGalleon Docker image: `local/genegalleon:performance-0.8.374`, immutable ID
  `sha256:ffdb1ab64997f0746620ed6c2e1ffde4b2ec9808316d949ebc3a9292be08f0c0`.
  The runtime wrapper verified current container inputs and the daily owned
  upstream snapshot. The default `dev` image was stale and was not used.
- Linux aarch64, Python 3.12.15, pandas 3.0.6, NumPy 1.26.4; Docker Desktop with
  12 virtual CPUs and 8,214,585,344 bytes of available memory. No additional
  per-container CPU or memory limit was applied.
- Fresh worker processes, discarded warmups and alternating before/after order.
  Fixtures and output verification are outside operation timers. Peak RSS is a
  worker high-water mark and can include imports, input creation and output
  verification; each benchmark documents its boundary. It is not whole-workflow
  memory. Large on-disk fixtures are written with bounded buffers. The initial
  completion-log RSS run inherited a large parent high-water mark and was
  discarded; its reported replacement uses Linux `VmHWM`.
- Comparisons require identical complete outputs. DataFrame comparisons include
  types and row/column order. Rooting comparisons restore exact bipartitions;
  complete adapter checks compare Newick, TSV and PDF, normalizing only temporary
  paths in captured logs and PDF creation/modification dates.
- These synthetic, warm-cache results establish neither production NAS
  throughput, whole-analysis speedups nor SIF compatibility.
  Genome fixtures use repeated ACGTN/IUPAC patterns; their gzip compression
  ratios do not represent a production assembly.

## All ten candidates

| # | Candidate | Decision and evidence |
| --- | --- | --- |
| 1 | Gene-count formatting: one object column at a time, integer counts and one result construction | Accept. Avoids repeated Arrow scalar access while retaining bounded handling of giant cells. Complete count/CLI TSVs match; 15 additional boundary comparisons passed. |
| 2 | Reconciliation root candidates: integer bipartitions | Accept. Shared unique tip bits replace descendant-name sets, complements and sorting. All candidate validation, rooting priority and publication fences remain. |
| 3 | Sign-flip resampling: smaller sign population | Accept. An int8 `[-1, 1]` population replaces int64 signs; multiplication and reductions remain float64. Results and consumed random state match. |
| 4 | UniProt descriptions: select keys before formatting | Accept. Unselected descriptions skip formatting; the complete input is still decoded, including broken tails. Dense and unfiltered calls retain their behavior. |
| 5 | BLAST coverage: singleton groups | Accept. A single HSP avoids small temporary arrays and reductions. Strict string joining, NaN lengths, reversed coordinates and integer boundary behavior remain. |
| 6 | HGT context: family-local source hashes | Accept. Each already-read immutable XML/stat text is hashed once per family. Fresh calls still observe changed sources; no persistent cache was added. |
| 7 | GFF summaries: reuse metadata and contig sets | Reject. Five calls on 2,000 genes × 3 blocks changed 0.3508→0.3349 s, but 300 genes × 100 blocks changed 1.0076→1.0764 s. RSS changed less than 0.4%; 179 exact comparisons passed. No consistent benefit. |
| 8 | Synteny summaries: remove a copy and replace merge with map | Reject. Three calls on 80,000 rows changed 0.561→0.575 s; adding 48 columns changed 3.384→3.302 s. RSS changed less than 0.5%; 18 exact comparisons passed. No material benefit. |
| 9 | OrthoFinder completion: stream the log | Reject. Completion checking occupied only 0.0483% of a 1,000-OG archive conversion. A 64 KiB log became 59% slower; a 64 MiB stress log used 208.62→24.01 MiB but became 62% slower. All sizes were synthetic assumptions, not observed production log sizes. 47 contract cases passed. |
| 10 | Genome length input: compiled IUPAC validation | Accept. Replaces per-line uppercase/set construction with an explicit ASCII character class, retaining the ASCII guard. 1,754 typed result/error comparisons passed. |

Candidates 7–9 left production code unchanged. Their private experiments remain
under `tmp/performance-20261023/`; weak screening results were not expanded into
large integration runs or permanent benchmark fixtures.

## Accepted measurements

Values are medians after discarded warmups. Gene counts, topology validation and
resampling received a second five-pair confirmation after the initial screening.

| Case | Before | After | Evidence/limit |
| --- | ---: | ---: | --- |
| Gene counts, 50,000 rows × 64 species, dense / sparse | 0.6205 / 0.5764 s | 0.4422 / 0.4268 s | Five pairs; every dense pair improves 17.5–32.2%. |
| Gene counts, 10,000 rows, wide (256 species) / giant cells (32 species) | 0.5841 / 0.0989 s | 0.3717 / 0.0543 s | Five pairs; complete TSV bytes match. |
| Gene counts, complete hog2og CLI, 20,000 × 64 | 0.7667 s | 0.6718 s | Five pairs; every paired call improves. |
| Gene counts, complete sp2og CLI control | 1.0455 s | 1.1090 s | Paired changes span −11.4% to +8.8%, paired median −0.1%; no consistent whole-CLI change claimed. |
| Root validation, 512 balanced tips × 128 candidates | 1.5674 s; 174.6 MiB | 0.05999 s; 109.8 MiB | Five confirmation pairs; 26.1×. |
| Root validation, 512 pectinate tips × 128 candidates | 4.6368 s; 181.1 MiB | 0.06167 s; 110.0 MiB | Five confirmation pairs; 75.2×. |
| Sign flips, 17 blocks × 100,000 draws | 0.00994 s | 0.01008 s | Short control; no difference claimed. |
| Sign flips, 1,001 blocks × 100,000 draws | 0.71875 s | 0.55727 s | Five confirmation pairs; complete result/state hashes match. |
| Sign flips, 10,001 blocks × 100,000 draws | 6.7850 s; 55.34 MiB | 5.4207 s; 48.73 MiB | Five confirmation pairs; 20.1% less time and 11.9% less peak RSS. |
| Complete resampling inference, 1,001 blocks | 0.6504 s; 55.76 MiB | 0.5604 s; 49.07 MiB | Five pairs; full result and final random-state hashes match. |
| Unchanged bootstrap control, 10,001 blocks | 0.10012 s | 0.10024 s | No performance difference claimed. |
| UniProt complete CLI, 300,000 records, raw / gzip | 0.692 / 0.711 s | 0.473 / 0.494 s | Three pairs; all TSV bytes match. Dense/unfiltered control gains are not claimed. |
| BLAST singleton / mixed, 60,000 rows | 2.041 / 1.457 s | 0.957 / 0.894 s | Five pairs, every paired run improves. Whole-worker RSS increases about 2.0% / 1.2%. |
| BLAST multi-HSP control, 300,000 rows | 1.658 s | 1.619 s | Five fixed-hash-seed pairs span −19.5% to +13.5%; no performance difference claimed. |
| HGT context, 512 events, small / 1 MB XML | 85.14 / 339.99 ms | 38.67 / 96.43 ms | Five pairs; complete event/link JSON bytes match. |
| HGT context, four events / 828 KB XML | 35.58 ms | 32.70 ms | Five pairs; complete output bytes match. |
| HGT context, one event | 0.337 ms | 0.349 ms | Five pairs with batched complete calls; overlapping distributions, no difference claimed. |
| Genome validation, 32 MiB, raw 64 / 80 columns | 0.3220 / 0.2784 s | 0.1485 / 0.1271 s | Five pairs; 2.17 / 2.19×. |
| Genome validation, 32 MiB, gzip 64 / 80 columns | 0.3609 / 0.3079 s | 0.1835 / 0.1564 s | Five pairs; 1.97×. |
| Genome validation, 32 MiB, one raw line | 0.1292 s; 118.6 MiB | 0.0541 s; 86.0 MiB | Five pairs; worker `VmHWM`. Mixed IUPAC and late-invalid inputs also match. |

Rooting's full 128-tip adapter also produced identical output bundles in three
measured pairs (4.413→4.075 s); native rooting/rendering costs and noise prevent
interpreting the candidate-validation speedup as a whole-workflow multiplier.
BLAST also passed 300 deterministic random DataFrame/TSV comparisons. UniProt
passed 70 reader comparisons and complete raw/gzip CLI comparisons. HGT context
passed 12 complete output/error comparisons, including later source changes.

## Reproduce accepted comparisons

From the repository root, save the baseline without switching branches:

```bash
gg_review_baseline=tmp/performance-20261023/baseline
mkdir -p "$gg_review_baseline" tmp/performance-20261023/{c01,c02,c03,c04,c05,c06,c10}
git archive a9f19346771e540a4a669a6470981dff91ec4fdc workflow/support |
  tar -x -C "$gg_review_baseline"
export GG_TEST_RUNTIME=docker
export GG_CONTAINER_DOCKER_IMAGE=local/genegalleon:performance-0.8.374
```

Run these sequentially, without concurrent tests or builds. Baseline files are
mounted with the checkout, and the runtime wrapper retains freshness checking.

```bash
bash workflow/tests/run_in_runtime.sh env PYTHONHASHSEED=0 python workflow/benchmarks/benchmark_gene_counts.py \
  --baseline "$gg_review_baseline/workflow/support/orthogroup_table_formatter.py" \
  --repeats 5 --output tmp/performance-20261023/c01/confirmation.json
bash workflow/tests/run_in_runtime.sh env PYTHONHASHSEED=0 python workflow/benchmarks/benchmark_root_candidate_topology.py \
  --baseline "$gg_review_baseline/workflow/support/species_tree_guided_gene_tree_rooting.py" \
  --cases topology --repeats 5 --output tmp/performance-20261023/c02/confirmation.json
bash workflow/tests/run_in_runtime.sh env PYTHONHASHSEED=0 python workflow/benchmarks/benchmark_subgenome_resampling.py \
  --baseline "$gg_review_baseline/workflow/support/subgenome_dominance.py" \
  --repeats 5 --output tmp/performance-20261023/c03/confirmation.json
bash workflow/tests/run_in_runtime.sh env PYTHONHASHSEED=0 python workflow/benchmarks/benchmark_uniprot_descriptions.py \
  --baseline "$gg_review_baseline/workflow/support/reformat_uniprot_diamond.py" \
  --profiles raw-cli gz-cli gz-dense gz-unfiltered --trials 3 \
  --output tmp/performance-20261023/c04/confirmation.json
bash workflow/tests/run_in_runtime.sh env PYTHONHASHSEED=0 python workflow/benchmarks/benchmark_blast_coverage.py \
  --baseline "$gg_review_baseline/workflow/support/annotate_blast_coverage.py" \
  --scale 3 --repeats 5 --output tmp/performance-20261023/c05/confirmation.json
bash workflow/tests/run_in_runtime.sh env PYTHONHASHSEED=0 python workflow/benchmarks/benchmark_hgt_context_hashes.py \
  --baseline-support "$gg_review_baseline/workflow/support" --candidate-support workflow/support \
  --repeats 5 --output tmp/performance-20261023/c06/confirmation.json
bash workflow/tests/run_in_runtime.sh env PYTHONHASHSEED=0 python workflow/benchmarks/benchmark_dotplot_genome_validation.py \
  --baseline-support "$gg_review_baseline/workflow/support" --candidate-support workflow/support \
  --bases 33554432 --trials 5 --output tmp/performance-20261023/c10/confirmation.json
```

## Validation

Each selection below used `bash workflow/tests/run_in_runtime.sh python -m pytest
-q SELECTION -x` with the Docker environment above. These are executed checks for
this review; the [development guide](development-and-tests.md#choose-checks-for-a-change)
remains the change-to-check reference.

| Selection | Result |
| --- | --- |
| `workflow/tests/test_orthogroup_table_formatter.py` | 11 passed |
| `workflow/tests/test_nwkit_tree_workflows_runtime.py -k 'root_selection or root_candidate_masks or root_adapter or root_comparison_reuse or root_reuse'` | 17 passed, 60 deselected |
| `workflow/tests/test_subgenome_dominance.py` | 71 passed, including final RandomState and strided-input cases |
| `workflow/tests/test_reformat_uniprot_diamond.py workflow/tests/test_focus_hgt_context_verification.py` | 42 passed |
| `workflow/tests/test_annotate_blast_coverage.py workflow/tests/test_extract_gene_id_from_blast_table.py workflow/tests/test_annotate_stat_branch_query_markers.py` | 26 passed |
| `workflow/tests/test_summarize_hgt_transfer_context.py` | 23 passed |
| `workflow/tests/test_pairwise_synteny.py` | 318 passed |

The accepted changes' pytest selections total **508 passed**, with no selected
test skipped. Rooting tests cover all priority choices, invalid
late candidates, singleton shape and input replacement before publication.
Resampling tests compare complete random draws, output probabilities, unchanged
inputs and final generator state for five NumPy bit generators; final additional
cases cover RandomState with integer/float32/float64 strided inputs.

`bash ./dev lint` passed with host Python 3.14.7 and Ruff 0.16.7 on PATH. Because
the host only supplied Bash 3.2, the documented `GG_LINT_BASH` override selected
an untracked `tmp/performance-20261023/bash-parser` wrapper delegating to the
same immutable Docker image's Bash 5.2.21 with the checkout mounted read-only.
Git, Ruff and the configuration checker ran on the host. The check validated all
75 tracked shell scripts plus `dev`, eight entrypoints and 25 common parameters.
`git diff --check` and final source/diff review passed. No version bump was
needed for this local-only commit.

SIF validation is unavailable on this macOS host; all executed workflow evidence
in this review is Docker evidence. R helpers, container inputs and entrypoints
were unchanged, so their substantial scientific integration lanes were not run.
