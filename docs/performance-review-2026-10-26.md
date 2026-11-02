# Twelve execution-speed candidates, 2026-10-26

Three subagents evaluated four proposals each. Promising changes received
repeated before/after measurements, complete output comparisons and review by
another agent. The parent reviewed the final code and rejected changes whose
benefit did not justify their complexity. Seven proposals were accepted and five
were rejected. Scientific models, thresholds,
identifiers, public formats and workflow stage order remain unchanged.

## Decisions

| # | Proposal | Decision |
| --- | --- | --- |
| 1 | Further vectorize RSC expression wide-to-long preparation | Reject after inspection: tuple iteration, column positions and observation templates already implement this optimization. |
| 2 | Reuse RSC usable-result masks and association rows within status aggregation | Accept. Compute selections once, preserving model-count/error order and empty-input behavior. |
| 3 | Remove per-value finite callbacks for bounded numeric P-values | Reject after measurement. About 2% incremental status improvement is noisy; the complete stat-tree read does not improve. Additional dtype branches are unjustified. |
| 4 | Select OMArk's top detected clade without sorting detected rows | Reject after measurement. Ordinary eight-clade controls regress in CPU; wall ratios cross unity and RSS is unchanged. |
| 5 | Retain only sole Complete BUSCO IDs while parsing full tables | Accept. Discard repeated-row bodies while validating every row/status, every marker and sole IDs in the original order. |
| 6 | Use set membership when selecting BUSCO marker records | Accept. Preserve the ordered marker list for JSON/sketches and preserve source record order. |
| 7 | Reuse binary-foreground tree traversal plans across traits | Accept. A lazy plan belongs to one resolve call; foreground flags remain independent and later calls observe tree changes. |
| 8 | Reduce repeated intron-ASR site/node aggregation | Reject after inspection: NWKIT owns inference/tree transfer/drawing and observed counts already use vectorized mapping. |
| 9 | Avoid repeated chromosome-color membership scans | Reject after inspection: distinct chromosome sets and dictionary membership already exist. |
| 10 | Count karyotype BED genes once per chromosome | Accept. One scan counts only selected chromosome IDs, bounding the Counter by the selection; exact IDs and selected order remain intact. |
| 11 | Encode scaffold plot size groups as categories | Accept for preparation. Preserve every plot coordinate, label and RNG call; complete drawing has no established speed gain. |
| 12 | Remove per-record generators and duplicate length calculation in FASTQ profiling | Accept. Preserve complete decoded hashes, validation, EOF behavior and the mutation fence. |

## Environment and measurement boundaries

The source baseline is `c6164ebe8709e15f552ccee7dddababd3052d781`.
Concurrent policy/test-only commits were integrated before delivery; none changed
the measured support modules. Existing unrelated GenomeScope work was excluded.
All fixtures are synthetic, outside curated workspace inputs and production
outputs. Raw fixtures, prototypes and initial records remain under
`tmp/performance-20261026-ykv8mdhu/`.

The old image failed the normal freshness check. `GG_BUILD_JOBS=4 bash ./dev build`
rebuilt the development runtime; its required tool checks and `pip check` passed.
Normal daily owned-source freshness then passed for immutable Docker image
`sha256:875a6b5926b47eabed44cd8de5f6820b416a378da2945e26edaf379317be438a`,
runtime-input identity
`10c8654f12e6d5fc8024d78f621ba57764f0a3512f4986d4f170bad61ab0c62a`.
Measurements used Linux/aarch64, Python 3.12.15, pandas 3.0.6 and NumPy 1.26.4.
The [measurement record](benchmarks/workflow-speed-2026-10-26.json) retains
paired samples, complete-output checks and final production source hashes.

Baseline-only measurements preceded production edits. Paired workers alternated
baseline/current order, discarded one warmup pair and retained three screening
pairs or five confirmation pairs. This review serialized its benchmarks, builds
and tests. Other desktop work was present, so noisy controls receive no speed
claim. Timers exclude imports, fixture construction and output verification;
complete render/reader cases include the named file I/O. Status/trait fixtures
load before timing. Peak RSS includes imports and loaded inputs, captured before
result hashing. Preliminary status baseline-only RSS had an inherited parent
allocation floor and is excluded from memory claims.

These are helper/phase measurements with warm caches. They do not establish
whole-workflow speedups, production NAS throughput or SIF compatibility.
Improvements from different phases must not be added together.

## Accepted measurements

Confirmation medians use five pairs, or three pairs for rows marked †.
Ratios are paired baseline/current wall times.

| Operation / synthetic workload | Wall before → after (s) | CPU before → after (s) | Peak RSS before → after (MiB) | Paired wall ratio |
| --- | ---: | ---: | ---: | ---: |
| RSC status, 200,000 result rows | 0.74796 → 0.34917 | 0.74787 → 0.34916 | 269.02 → 245.57 | 1.53–2.69× |
| BUSCO table, 10,000 singleton markers | 0.01230 → 0.00872 | 0.01187 → 0.00834 | 54.79 → 51.96 | 1.28–2.04× |
| BUSCO table, 27,500 rows with wide metadata | 0.04098 → 0.02128 | 0.03993 → 0.02092 | 69.82 → 51.91 | 1.63–2.70× |
| Marker projection, 64 species × 5,000 markers / 200 selected | 0.42443 → 0.00935 | 0.42441 → 0.00935 | 60.14 → 60.14 | 32.2–48.7× |
| Binary foregrounds, balanced 1,024 tips × 64 traits | 1.02337 → 0.23143 | 1.02332 → 0.23139 | 47.69 → 47.67 | 3.74–8.56× |
| Karyotype counts, 200,000 genes / 32 selected chromosomes † | 0.61584 → 0.05518 | 0.61581 → 0.05517 | 185.51 → 185.43 | 10.84–12.36× |
| Karyotype counts, 200,000 genes / 100,000 contigs / two selected † | 0.00748 → 0.00484 | 0.00748 → 0.00484 | 191.73 → 191.80 | 1.47–1.59× |
| Complete karyotype SVG/report, two 200,000-gene / 32-chromosome tracks † | 8.17824 → 7.58669 | 8.16217 → 7.57233 | 635.55 → 635.70 | 1.08–1.12× |
| Scaffold preparation and three jitter scatters, 1,000,000 rows | 0.76244 → 0.51514 | 0.71295 → 0.47467 | 269.71 → 271.11 | 1.06–1.57× |
| FASTQ scan, 500,000 reads × 150 bases, plain | 0.97072 → 0.87593 | 0.91021 → 0.79655 | 19.88 → 19.82 | 1.06–1.21× |
| Same FASTQ, gzip † | 1.00387 → 0.82539 | 1.00291 → 0.82474 | 19.99 → 20.13 | 1.10–1.57× |

The BUSCO wide profile reduces peak RSS 25.6%; all five confirmation pairs improve
in wall and CPU time. The initial three-pair screen had one wall-time reversal.
Single-trait and tiny/nonbinary controls show no stable speed gain. Tree-plan
RSS is essentially unchanged. The marker set and plan add O(selected markers)
and O(tree nodes) state; BUSCO retention is O(marker count) plus one split row.

RSC comparisons match complete DataFrames, including dtypes, index, row order
and input immutability. The source matches the isolated mask-only measurement.
The rejected finite-callback prototype gave 0.31851→0.31197 s against the mask-only
baseline, with ratios 0.95–1.15×; complete stat-tree reading gave
0.56140→0.58571 s. OMArk's 500-species/eight-clade descending control increased
CPU 0.05498→0.06252 s, with every paired CPU ratio below one. Neither prototype
is included in the delivered implementation.

The initial karyotype Counter retained unselected chromosomes too. Independent
review prompted a selected-only Counter and a 100,000-contig control; retained
count state is now O(selected chromosomes). Missing selected chromosomes still
fail before publication. Complete final-source SVG/report comparisons are
included in the table. BED parsing and drawing are outside count-phase timers.
The tiny 100-gene/two-chromosome count control regresses from 11.75 to 32.58
microseconds; no general small-input speedup is claimed.

Scaffold preparation times the real `run` prefix plus its three jittered scatter
helpers. The fourth scatter and complete rendering are covered separately by
full SVG/PDF byte comparisons. Complete 20,000-row rendering has paired ratios
0.98–1.11×, so no full-plot speedup is claimed. Preparation adds about 1.4 MiB of
peak RSS. Tiny scaffold/FASTQ and unchanged FASTA controls show no stable gain.
FASTQ comparisons include every returned hash/statistic; invalid late records,
partial EOF, truncated gzip trailers and mutation fences retain their behavior.

## Reproduce

Run serially from the repository root with the image above available:

```bash
export GG_TEST_RUNTIME=docker
export GG_CONTAINER_DOCKER_IMAGE=sha256:875a6b5926b47eabed44cd8de5f6820b416a378da2945e26edaf379317be438a
mkdir -p tmp
GG_PERF_BASE=$(mktemp -d "$PWD/tmp/performance-baseline-XXXXXX")
git archive c6164ebe8709e15f552ccee7dddababd3052d781 workflow/support | tar -x -C "$GG_PERF_BASE"

bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_summary_candidates.py \
  --baseline "$GG_PERF_BASE/workflow/support" --mode status --rows 200000 --repeats 5 \
  --fixtures "$GG_PERF_BASE/status" --output "$GG_PERF_BASE/status.json"

bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_guide_foregrounds.py \
  --baseline "$GG_PERF_BASE/workflow/support/busco_guide_tree.py" \
  --candidate workflow/support/busco_guide_tree.py --repeats 5 \
  --profiles rows-singletons rows-wide markers-default \
  --work-dir "$GG_PERF_BASE/busco" --output "$GG_PERF_BASE/busco.json"

bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_guide_foregrounds.py \
  --baseline "$GG_PERF_BASE/workflow/support/resolve_binary_trait_foregrounds.py" \
  --candidate workflow/support/resolve_binary_trait_foregrounds.py --repeats 5 \
  --profiles traits-many traits-single \
  --work-dir "$GG_PERF_BASE/traits" --output "$GG_PERF_BASE/traits.json"

bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_karyotype_input_profile.py \
  --prepare --fixtures "$GG_PERF_BASE/plot-inputs"

bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_karyotype_input_profile.py \
  --baseline "$GG_PERF_BASE/workflow/support/pairwise_synteny_karyotype.py" \
  --source workflow/support/pairwise_synteny_karyotype.py --case count-large --repeats 3 \
  --fixtures "$GG_PERF_BASE/plot-inputs" --output "$GG_PERF_BASE/count.json"

bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_karyotype_input_profile.py \
  --baseline "$GG_PERF_BASE/workflow/support/scaffold_size_histogram.py" \
  --source workflow/support/scaffold_size_histogram.py --case scatter-large --repeats 5 \
  --fixtures "$GG_PERF_BASE/plot-inputs" --output "$GG_PERF_BASE/scatter.json"

bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_karyotype_input_profile.py \
  --baseline "$GG_PERF_BASE/workflow/support/resource_profile.py" \
  --source workflow/support/resource_profile.py --case resource-fastq-large --repeats 5 \
  --fixtures "$GG_PERF_BASE/plot-inputs" --output "$GG_PERF_BASE/fastq.json"
```

Use `count-sparse`, `count-small`, `render-large` and `render-small` for the other
karyotype cases. Use `scatter-small`, `scaffold-render-large` and
`scaffold-render-small` with the scaffold module, and `resource-fastq-gzip`,
`resource-fastq-small` and `resource-fasta-control` with the resource module.
These use three measured pairs. Complete render comparisons fix the SVG hash
salt and source date and compare all bytes plus the complete karyotype report.

BUSCO projection executes the actual `build` assignment extracted from its AST,
including set construction; this timing excludes sketching/tree construction.
Complete runtime tests exercise those native integrations separately. The row
parser compares ordered IDs and the full marker universe, projection compares
ordered species/records, and trait comparisons include every changed table cell.

## Validation

The final production sources passed **136 Docker tests with zero skips**:

```bash
bash workflow/tests/run_in_runtime.sh python -m pytest -q \
  workflow/tests/test_busco_guide_tree.py workflow/tests/test_busco_guide_tree_runtime.py \
  workflow/tests/test_resolve_binary_trait_foregrounds.py -x

bash workflow/tests/run_in_runtime.sh python -m pytest -q \
  workflow/tests/test_reconciled_speciation_contrast.py workflow/tests/test_summarize_omark.py -x

bash workflow/tests/run_in_runtime.sh python -m pytest -q \
  workflow/tests/test_resource_planning.py workflow/tests/test_scaffold_size_histogram.py \
  workflow/tests/test_pairwise_synteny_runtime.py \
  -k 'karyotype or profile or scaffold_plot or metrics or success_publishes or core_parallel' -x
```

The groups passed 58 tests in 27.81 s, 31 in 3.29 s, and 47 in 5.09 s;
the final selection deselected 63 unrelated tests. The 58 include real
k-mer/RapidNJ/BUSCO integration. Thirty-four differential summary edge/error
comparisons also passed, including missing/nullable values, complex values,
duplicate indexes, empty inputs and error priority. The final benchmark sources
match the delivered production hashes.

A supplementary broad three-file run passed 110 tests in 352.75 s without skips,
with three existing dual-codon warnings. The selected-only Counter landed during
that run, so those 110 are not counted as final-source evidence. The final 47
tests were started afterward. No production source changed after its final tests.

`bash ./dev lint` passed, including all shell syntax, host Ruff and configuration
checks for eight entrypoints and 27 common parameters. The macOS host has Bash
3.2; a temporary `GG_LINT_BASH` wrapper delegated syntax parsing to Bash in the
same Docker image, as allowed by the development guide. Ruff/config ran on the
host; the parser container was stopped afterward. `git diff --check` passed.
These static results are separate from executed Docker behavior.

Publication uses main base `ed33c16c992b3fb6a3cfc801eee93ffab9013d82` and patch
version 0.8.422. Existing unrelated tracked edits were compared byte-for-byte
against the initial saved diff and preserved. Full scientific workflows, the
complete runtime/R lanes and SIF validation were not run for these focused
helper changes.
