# Workflow performance follow-up, 2026-10-23

Four further ideas were evaluated after the [first review](performance-review-2026-10-23.md).
Three were accepted; the running-maximum proposal was rejected. Three candidates
were delegated to separate agents, and scaffold attachment received an additional
agent review. Production changes preserve scientific calculations and output
schemas. A discovered index-grouping defect was fixed before extending bulk
assignment to branch annotations.

## Comparison conditions

- Baseline: `a6d7d7f30c8a407de8533d0270d6d8e92a0332a0`, including the first review's changes.
- Freshness-checked GeneGalleon Docker image `local/genegalleon:performance-0.8.374`,
  immutable ID `sha256:ffdb1ab64997f0746620ed6c2e1ffde4b2ec9808316d949ebc3a9292be08f0c0`.
  Linux aarch64, Python 3.12.15, pandas 3.0.6; Docker Desktop, 12 vCPUs and
  8,214,585,344 bytes of memory, with no additional container limit.
- Baselines were measured before implementation. Final comparisons discarded a
  warmup pair and used five alternating pairs of fresh workers. Small homoeolog
  profiles batch complete calls. Imports, fixture creation and result verification
  are outside operation timers. Peak RSS is Linux worker `VmHWM` through the
  operation, including imports and input construction, before output hashing.
- Review measurements and runtime tests were serialized. Unrelated machine work
  remained possible: even unchanged controls varied appreciably. Medians alone
  are not treated as proof when paired changes cross zero.
- These are synthetic, warm-cache helper measurements. They do not establish
  whole-workflow speed, production storage throughput or SIF compatibility.

## Decisions and measurements

### 1. Enumerate supported homoeolog tiles

Accepted. `homoeolog_order` builds chromosome partners in one pass and uses a
query-to-target-position index to visit only target pairs sharing at least two
query partners. Exact rational harmonic scores, tie order, anchor deduplication,
complete input validation and greedy grouping remain unchanged. Position-based
enumeration also preserves repeated selected chromosomes.

| Synthetic input | Wall median before → after | Paired speedup range | Peak RSS before → after |
| --- | ---: | ---: | ---: |
| 500 scaffolds per track, sparse tiles | 47.21 → 10.14 ms | 2.87–5.79× | 24.64 → 24.92 MiB |
| 500 scaffolds, diagonal/unplaced | 34.86 → 2.43 ms | 11.98–20.01× | 23.55 → 23.82 MiB |
| 4,000 scaffolds per track, sparse tiles | 3.6323 → 0.07690 s | 37.73–51.98× | 44.83 → 45.91 MiB |
| Dense, two targets × 500 queries | 1.3831 → 1.3991 s | 0.85–1.21× | 107.22 → 107.46 MiB |

CPU changes track wall time. Normal 24-chromosome, dense eight-chromosome and
500-target/two-query controls had overlapping paired distributions; no speed
change is claimed for them. Sparse enumeration trades about 1.1 MiB at 4,000
scaffolds for the measured reduction in work. All complete typed-result hashes
match. An additional 2,702 differential comparisons cover random graphs, ties,
duplicates, alternate sequence containers and exact exceptions.

### 2. Construct promoter N ranges directly

Accepted. `get_N_coordinate` uses contiguous `[Nn]+` matches instead of retaining
every N position and joining adjacent positions in Python. A substring guard
skips the regex when neither N nor n occurs. FASTA parsing, gzip, inclusive
one-based coordinates, terminal runs and existing outside-scaffold ranges remain.

Ordinary profiles contain 2,000 promoters of 2,000 bases, wrapped at 80 columns.

| Synthetic input | Wall median before → after | CPU median before → after | Peak RSS before → after |
| --- | ---: | ---: | ---: |
| Run-rich raw FASTA | 307.22 → 46.14 ms | 307.12 → 46.15 ms | 85.90 → 84.85 MiB |
| Same sequences, gzip | 323.49 → 53.32 ms | 323.44 → 53.32 ms | 85.92 → 84.91 MiB |
| N-free control | 47.83 → 16.67 ms | 47.83 → 16.67 ms | 85.84 → 84.90 MiB |
| Alternating N/non-N | 488.46 → 433.53 ms | 488.44 → 433.52 ms | 91.26 → 89.54 MiB |
| One 2-million-base all-N stress record | 539.22 → 17.45 ms | 539.21 → 17.44 ms | 153.24 → 82.93 MiB |

Run-rich, N-free and stress inputs improved in every pair. Alternating input
was slower in one pair; its smaller median difference is variable. The stress
record is not an ordinary promoter. Normal-profile memory differences are small.

The unguarded proposal made N-free input 48.2% slower in the first full comparison;
the guard addresses that measured regression. The final candidate was remeasured
across all profiles. Output schemas, dtypes, row order and complete TSV hashes
match. Differential verification includes all 5,461 short sequences over ANnC
of lengths zero through six, raw/gzip, wrapped input, empty records and headers.
Permanent tests also execute the join into the branch-annotation CLI output.

### 3. Attach scaffold context using required columns and bulk assignment

Accepted. Branch mapping now iterates the two required gene columns instead of
constructing records containing every context column. Branch metrics are assigned
by column through the existing bulk helper. Candidate-locus exclusions, recipient
selection, unique-scaffold pooling and fraction calculations remain unchanged.

| Synthetic input | Wall median before → after | CPU median before → after | Peak RSS before → after |
| --- | ---: | ---: | ---: |
| 64 taxonomy genes/species, 33 branches | 0.2890 → 0.2474 s | 0.2868 → 0.2454 s | 89.96 → 89.91 MiB |
| 2,048 taxonomy genes/species, 33 branches | 1.3980 → 1.2489 s | 1.3939 → 1.2454 s | 101.36 → 101.43 MiB |
| 512 taxonomy genes/species, 528 branches | 2.9441 → 0.8260 s | 2.9410 → 0.8234 s | 102.99 → 103.96 MiB |

Each input has two species. The many-branch case has 514 candidate gene rows;
every pair improves by 2.25–5.05×, with a 3.56× ratio of medians and about 1 MiB
additional peak memory. Smaller controls cross unity (0.74–1.45× and
0.75–1.64×); no reliable speed change is claimed for those controls. Complete
DataFrame schemas, dtypes, indices, row/column order and TSV bytes match.

Review found that the existing gene bulk helper merged different missing index
labels with `factorize`. A regression reproduced incorrect propagation between
`None`, NaN, NA and NaT. Grouping now uses the index's own unique-label lookup,
matching scalar `.at` updates. This intentionally fixes the existing mixed-missing
gene-index defect as well as making the new branch path correct. Scalar-update
oracles cover genes and branches, partial updates, categorical, nullable integer,
datetime, duplicate and unique tuple indices, and unchanged caller inputs.
Thirty-two additional baseline comparisons cover branch results/errors and
empty, partial and missing-orthogroup cases, including mixed missing labels.

### 4. Keep only running window maxima

Rejected. Keeping one interchromosomal maximum per window reduced peak RSS on
480,000 rows from 46,168 to 26,740 KiB (42.1%). Sixty exact controls passed,
including signed-zero ties, aggregation input order and errors after valid rows.
However, a five-call-per-worker confirmation increased many-query CPU time in
every pair, by 6.4–46.9% (median 21.7%). A subsequent three-pair speed-only check
within fresh workers remained noisy, as did the unchanged intrachromosomal path.
The roughly 19 MiB saving did not justify an unresolved CPU cost. Production code
and existing tests for this candidate were restored exactly; experiments remain
under `tmp/performance-20261023-round2/self-fractionation/`.

## Reproduce accepted comparisons

From the repository root, save the baseline without changing branches:

```bash
gg_followup_baseline=tmp/performance-20261023-round2/baseline
mkdir -p "$gg_followup_baseline"
git archive a6d7d7f30c8a407de8533d0270d6d8e92a0332a0 workflow/support |
  tar -x -C "$gg_followup_baseline"
export GG_TEST_RUNTIME=docker
export GG_CONTAINER_DOCKER_IMAGE=local/genegalleon:performance-0.8.374
```

Run sequentially, without concurrent tests/builds:

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_homoeolog_order.py \
  --baseline-support "$gg_followup_baseline/workflow/support" --candidate-support workflow/support \
  --trials 5 --output tmp/performance-20261023-round2/homoeolog-reproduced.json
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_promoter_n_runs.py \
  --baseline "$gg_followup_baseline/workflow/support/orthogroup_statistics.py" \
  --output tmp/performance-20261023-round2/promoter-reproduced.json
bash workflow/tests/run_in_runtime.sh env PYTHONHASHSEED=0 python workflow/benchmarks/benchmark_scaffold_context.py \
  --baseline-support "$gg_followup_baseline/workflow/support" --genes 512 --branch-repeats 16 \
  --repeats 5 --output tmp/performance-20261023-round2/scaffold-many-reproduced.json
```

For scaffold controls, replace `--genes 512 --branch-repeats 16` with
`--genes 64 --branch-repeats 1` or `--genes 2048 --branch-repeats 1` and use
separate output paths. Its original single-implementation CLI and JSON fields
remain available. The promoter benchmark extracts the actual nested functions
from each source AST; it includes file reading and DataFrame construction.
Homoeolog timing includes complete anchor parsing. Raw samples, source/input
hashes and verification logs are retained in `tmp/performance-20261023-round2/`.

## Validation

The following selections use the freshness-checked Docker runtime above:

```bash
bash workflow/tests/run_in_runtime.sh python -m pytest -q \
  workflow/tests/test_orthogroup_statistics_promoter.py workflow/tests/test_orthogroup_statistics_fimo.py \
  workflow/tests/test_orthogroup_statistics_support.py workflow/tests/test_shell_static_safety.py -x -p no:cacheprovider
bash workflow/tests/run_in_runtime.sh python -m pytest -q \
  workflow/tests/test_pairwise_synteny.py workflow/tests/test_scaffold_taxonomy.py \
  workflow/tests/test_score_hgt_candidates.py -x
```

The first selection passed **284 tests in 23.21 s**; the second passed **404 in
59.09 s**: **688 passed**, with no selected skips. The mixed-missing regression
first failed against the old grouping, then passed after its repair. Production
helper SHA256 values match the sources used for final measurements.

`bash ./dev lint` also passed, including Ruff and entrypoint configuration checks
(eight entrypoints, 25 common parameters). The host has Bash 3.2, so the actual
invocation used the documented `GG_LINT_BASH` override with a Docker GNU Bash
5.2.21 parser and host Python 3.14.7/Ruff 0.16.7:

```bash
env PATH="/Users/kf/miniforge3-arm64/bin:$PATH" \
  GG_TEST_RUNTIME=docker GG_CONTAINER_DOCKER_IMAGE=local/genegalleon:performance-0.8.374 \
  GG_LINT_BASH="$PWD/tmp/performance-20261023/bash-parser" bash ./dev lint
git diff --check
```

The local parser wrapper is retained with the review artifacts and dispatches
`bash "$@"` in the immutable image above with the checkout mounted read-only.
Final diff and whitespace checks passed. Lint is static evidence; the pytest
results are executed Docker evidence. SIF and a complete scientific workflow
were not run. Container inputs were unchanged, so no image rebuild was needed.
