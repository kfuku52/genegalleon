# Workflow performance review, round three, 2026-10-23

This review starts from `74d4bfd71a53bab71a32a5fe8c8936ebeaa3aa31`
(0.8.388), after the [previous follow-up](performance-followup-2026-10-23.md).
Three agents investigated six further ideas. Measurements and runtime checks
were serialized; production changes preserve calculations, validation, ordering,
identifiers, output schemas and publication behavior.

## Comparison environment

Final comparisons and validation use the freshness-checked GeneGalleon Docker
image `local/genegalleon:performance-round3-20261023`, immutable ID
`sha256:0e46be37b7ddbdb60b331ef35017abd3f2d90820d152ce3838ec369537791752`:
Linux aarch64, Python 3.12.15, pandas 3.0.6, NumPy 1.26.4, R 4.4.3,
ggplot2 4.0.3, ggtree 4.0.4 and ape 5.8.1. Docker Desktop provides 12 CPUs and
8,214,585,344 bytes of memory; no extra container limit was applied.

Baselines were established before editing production code. The R package change
is a container input, so the normal local build created the image above with
`BUILD_SIF=0 IMAGE_SOURCE=local IMAGE=local/genegalleon
TAG=performance-round3-20261023 BUILD_TARGET=development bash
./gg_container_build_entrypoint.sh`. The build initially missed dependency-layer
caches; this is not a build-speed measurement. `GG_RUNTIME_FRESHNESS=always`
refreshed the previous daily owned-source snapshot and verified the new image.
Subsequent checks retained the normal freshness policy. All final baseline/
candidate pairs, including a repeat of HGT after the rebuild, use this same image.

Synthetic fixtures and warm file caches do not establish production storage
throughput, whole-analysis speedups or SIF compatibility. Raw samples, complete
source/input/output hashes and logs are retained under
`tmp/performance-20261023-round3/` and `tmp/performance-round3-20261023/`.

## Candidates

| Idea | Decision and scope |
| --- | --- |
| Index HGT migration gene records by orthogroup | Accept: remove one full gene-dictionary scan per orthogroup summary. This accelerates migration of previously written HGT tables, not fresh HGT scoring or database lookup. |
| Use the built-in first nonmissing BUSCO metadata reduction | Accept: retain independent first URL/description selection across ordered species without per-group Python callbacks. |
| Cache valid-amino-acid masks per aligned sequence | Accept: avoid repeated alphabet membership work while retaining pairwise means, missing-tip behavior and CDS translation. |
| Reuse child-clade expression-correlation minima | Deferred: promising asymptotic reduction, but requires separate verification of NaN, signed-zero and reduction-order semantics. No production change or speed claim. |
| Filter coverage before parsing legacy BLAST evalues | Deferred: the normal producer already provides `min_evalue`, so the proposed gain mainly concerns a legacy input path. No production change or speed claim. |
| Normalize FASTA headers during collection | Deferred: the reader already streams input; a material whole-command gain is unproven. No production change or speed claim. |

## HGT taxonomy migration

The existing `(orthogroup, gene_id)` dictionary still retains the first duplicate
record. Its values are indexed once by orthogroup, in insertion order. Branch
candidate order, rank aggregation, existing values, row counts and temporary-file
publication remain unchanged.

These synthetic fixtures contain ten interleaved genes per orthogroup and
prepopulated taxonomy fields, with taxonomy database lookup disabled. Timings
include the real parser/main, target imports, all table I/O and replacements;
interpreter/harness startup, fixture generation and verification are excluded.
Five alternating fresh-worker
pairs follow a discarded warmup pair. Linux `VmHWM` is read before output hashing.

| Orthogroups / genes | Wall median before → after | Paired speedup range | CPU median before → after | Peak RSS before → after |
| --- | ---: | ---: | ---: | ---: |
| 32 / 320 | 0.5623 → 0.5541 s | 0.99–1.03× | 0.5569 → 0.5497 s | 140.98 → 140.96 MiB |
| 2,000 / 20,000 | 4.3049 → 2.0591 s | 1.89–2.18× | 4.2999 → 2.0546 s | 302.79 → 303.13 MiB |
| 4,000 / 40,000 | 13.7664 → 3.6955 s | 3.41–4.02× | 13.7592 → 3.6903 s | 467.25 → 468.25 MiB |

Every medium/large pair improved. The small control overlaps noise and has no
claimed speed change. The index costs about 0.3–1.0 MiB in these larger cases.
All three output TSVs and stdout/stderr match byte for byte in every measured
worker. Regression tests cover interleaved
groups, duplicate keys, the same gene in different groups, blank/missing groups,
existing fields, real resolver logic with a fake NCBI backend, and failure cleanup.

The final raw comparison is
`tmp/performance-20261023-round3/hgt-taxonomy-final.json`; the earlier screen and
confirmation in the previous image are retained separately.

## Sequence similarity

Validity masks are computed once per aligned amino-acid sequence. Each pair
still uses the same intersected positions and `mean(a[valid] == b[valid])`.
FASTA validation, CDS translation, genetic codes, missing-tip handling and plot
construction are unchanged. A second agent reviewed the change and benchmark.

Five alternating fresh-R-process pairs follow a discarded warmup pair, across
nine profile/mode combinations. The kernel runs the actual producer with only
the final panel constructor replaced by data capture. Panel timings include
FASTA reading, similarity calculation and plot construction, excluding package
loading, fixture generation, output hashing and rendering. Batched controls are
reported per call below. CPU medians closely follow wall time.

| Tips × amino-acid sites / scope | Wall median before → after | Paired speedup range | Peak RSS before → after |
| --- | ---: | ---: | ---: |
| 7 × 300, panel | 75.67 → 77.00 ms | 0.98–1.02× | 216.23 → 216.29 MiB |
| 64 × 300, panel | 120.0 → 116.0 ms | 1.02–1.10× | 229.79 → 238.55 MiB |
| 256 × 1,000, kernel | 1.372 → 0.981 s | 1.36–1.41× | 229.64 → 230.68 MiB |
| 256 × 1,000, panel | 1.466 → 1.061 s | 1.32–1.43× | 236.33 → 246.43 MiB |
| 256 × 1,000, ambiguous residues, panel | 1.294 → 0.859 s | 1.49–1.54× | 238.28 → 240.02 MiB |
| 512 × 300, kernel | 2.729 → 2.255 s | 1.17–1.23× | 253.11 → 253.84 MiB |

Every large paired trial improves. The small panel control crosses unity and has
no claimed speed change. This trades additional mask storage for CPU time: peak
RSS rises about 4.3% in the dense 256-tip panel, and up to 5.4% in the separate
64-tip kernel profile. No memory saving or rendering-time improvement is claimed.
The complete raw comparison is retained in
`tmp/performance-round3-20261023/treevis-comparison.json`.

All 90 measured workers have matching complete output fingerprints. Independent
tests compare full typed panel data, order, NA values, plot layers, PNG bytes and
grob dimensions. They cover all standard amino acids, ambiguous residues, gaps,
stop characters, Unicode, FASTA order and missing tips, plus seven exact error
contracts. The original CDS/genetic-code and real FIMO integrations also pass.

## BUSCO metadata

The collector uses the built-in first nonmissing reduction after the existing
`"-"` → `pd.NA` conversion, then restores all-missing values to `"-"`.
URL and description are still chosen independently from ordered species.
The existing first-row metadata choice within each species, duplicate sequence
joins, identifier handling and `first_non_missing` callable remain unchanged.
A second agent independently reviewed the production diff, tests and benchmark.

The same new Docker image and five alternating fresh-worker pairs are used.
The complete CLI timer includes imports and `main`, excluding interpreter startup,
fixture generation and output verification. The metadata kernel is extracted
from the actual `main` AST; its inputs are prepared using the normal reader.

| Input markers × species / workers | CLI wall median before → after | Paired speedup range | CPU median before → after | Peak RSS before → after |
| --- | ---: | ---: | ---: | ---: |
| 32 × 2 / 1 | 0.1833 → 0.1752 s | 0.98–1.11× | 0.1818 → 0.1733 s | 74.85 → 75.03 MiB |
| 2,000 × 10 / 1 | 0.4923 → 0.3340 s | 1.45–1.54× | 0.4908 → 0.3327 s | 81.14 → 81.21 MiB |
| 10,000 × 50 / 1 | 4.6887 → 3.9529 s | 1.15–1.23× | 4.6867 → 3.9509 s | 222.66 → 222.31 MiB |
| 2,000 × 10 / 2 | 0.5446 → 0.3644 s | 1.35–1.63× | 0.5452 → 0.3563 s | 81.80 → 82.10 MiB |

All medium/large CLI pairs improve. The small CLI control crosses unity, and RSS
differences are small; no small-input speed or memory improvement is claimed.
The metadata-only reductions are larger: 0.17693 → 0.006498 s for 2,000 × 10
and 0.97938 → 0.14083 s for 10,000 × 50. These are not whole-command multipliers.

Complete TSV bytes, schemas and dtypes match in every benchmark. Serial stdout
matches exactly; two-worker progress lines retain their variable completion
order and are compared as multisets, including occurrence counts. An additional
500 typed metadata-frame comparisons and 50 complete CLI/error comparisons match.
Tests cover complementary metadata, first-row duplicate suppression, all-empty
species, literal NA-like values, leading-zero identifiers, quoted Unicode/tab/
newline fields and invalid UTF-8 before output replacement. Raw results are in
`tmp/performance-20261023-round3/busco/`.

## Reproduce accepted comparisons

Save the baseline without switching branches, then use the matching runtime:

```bash
gg_round3_baseline=tmp/performance-20261023-round3/reproduced-baseline
mkdir -p "$gg_round3_baseline"
git archive 74d4bfd71a53bab71a32a5fe8c8936ebeaa3aa31 workflow/support |
  tar -x -C "$gg_round3_baseline"
export GG_TEST_RUNTIME=docker
export GG_CONTAINER_DOCKER_IMAGE=local/genegalleon:performance-round3-20261023
```

Run sequentially, without concurrent tests or builds:

```bash
bash workflow/tests/run_in_runtime.sh env PYTHONHASHSEED=0 python workflow/benchmarks/benchmark_hgt_taxonomy_migration.py \
  --baseline "$gg_round3_baseline/workflow/support/add_hgt_taxonomy_columns.py" \
  --profiles small screen large --repeats 5 \
  --output tmp/performance-20261023-round3/hgt-reproduced.json
bash workflow/tests/run_in_runtime.sh Rscript --vanilla workflow/benchmarks/benchmark_treevis_sequence_similarity.R \
  --baseline-source "$gg_round3_baseline/workflow/support/treevis/R/05_pairwise_similarity.R" \
  --candidate-source workflow/support/treevis/R/05_pairwise_similarity.R \
  --profiles small,dense64,dense256,ambiguous256,dense512 --repeats 5 \
  --output tmp/performance-20261023-round3/treevis-reproduced.json
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_busco_metadata.py \
  --baseline "$gg_round3_baseline/workflow/support/collect_common_BUSCO_genes.py" \
  --markers 2000 --species 10 --ncpu 1 --repeats 5 \
  --output tmp/performance-20261023-round3/busco-reproduced.json
```

For BUSCO controls, use `--markers 32 --species 2`, or
`--markers 10000 --species 50`, with distinct output paths. Its parallel comparison
uses `--markers 2000 --species 10 --ncpu 2 --cases cli`.

## Validation

The [development guide](development-and-tests.md#choose-checks-for-a-change)
remains the check-selection reference. The final Docker commands were:

```bash
bash workflow/tests/run_in_runtime.sh python -m pytest -q \
  workflow/tests/test_add_hgt_taxonomy_columns.py workflow/tests/test_score_hgt_candidates.py \
  workflow/tests/test_hgt_taxonomy_missingness.py workflow/tests/test_busco_copy_selection.py -x
bash workflow/tests/run_in_runtime.sh python -m pytest -q \
  workflow/tests/test_collect_common_BUSCO_genes.py -x -p no:cacheprovider
```

The first passed **71 tests in 7.28 s**; the second passed **10 in 2.85 s**:
**81 passed, no selected skips**. R validation parsed all three changed/new R
files, ran `bash workflow/tests/check_treevis_package.sh` with **Status: OK**,
installed the working package into a temporary library and asserted that it
loaded the candidate function. The full
`Rscript --vanilla workflow/tests/test_treevis_sequence_cis_similarity.R`
then passed for the installed candidate and again with the saved baseline source.
The existing empty-colour fixture warnings were identical in both runs.
The exact executed R sequence is retained in
`tmp/performance-round3-20261023/validate-treevis.sh`, with its output in
`treevis-validation.log`.

Host lint also passed using Python 3.14.7, Ruff 0.16.7 and the documented
`GG_LINT_BASH` override. Its ignored parser wrapper dispatches to GNU Bash in
the same immutable image with a read-only checkout and networking disabled:

```bash
env PATH="/Users/kf/miniforge3-arm64/bin:$PATH" \
  GG_LINT_BASH="$PWD/tmp/performance-20261023-round3/bash-parser" bash ./dev lint
git diff --check
```

Lint includes Ruff, shell syntax and configuration checks (eight entrypoints,
25 common parameters); it is static evidence. The behavior evidence above is
Docker evidence. SIF and a complete scientific workflow were not run for these
changes. No scientific settings, upstream source defaults or version number
were changed. This is a local-only commit.
