# Twelve execution-speed candidates, round two, 2026-10-25

Twelve separate agents evaluated another set of execution-speed proposals after
the first review. Six improvements were selected after repeated measurements,
complete output comparisons and compatibility review.
Scientific models, thresholds, identifiers, public schemas and workflow stage
order remain outside this optimization scope. Curated workspace inputs are not
benchmark fixtures.

## Decisions

| # | Proposal | Decision / evidence |
| --- | --- | --- |
| 1 | Aggregate repeated BUSCO annotation rows without per-gene Series callbacks | Accept. Five-pair confirmation improves repeated-hit inputs and complete merge main; singleton-only input has no stable gain. Typed tables, TSVs and parser errors match. |
| 2 | Cache node-label strings and ancestor strings when constructing tree relations | Accept. Descendant order still comes from the original traversal. Balanced and ladder trees improve in five-pair confirmation; tiny and star controls have no stable gain. |
| 3 | Group FASTA identifiers once by species when planning GFF work | Accept for large multispecies planning. Preserve the species parser, FASTA/file order, duplicate rejection and exact per-file lists. Complete CLI comparisons match but have no stable whole-command gain. |
| 4 | Reuse lazy species aliases in presence/absence matrix construction | Reject at screening. Complete raw/ZIP calls have no meaningful gain. |
| 5 | Batch promoter requests and discover each genome once | Reject after inspection. Requests are already grouped with one SeqKit call per species; hoisting discovery would shift a mutable-filesystem observation boundary. |
| 6 | Avoid allocating row dictionaries for unselected WGD evidence rows | Accept. Sparse, dense and wide inputs improve in five-pair confirmation. All source rows still receive global validation, including rows outside the requested family. |
| 7 | Use the native indexed genome's reference membership lookup | Accept for fragmented assemblies. The native `pysam.FastaFile` lookup avoids a linear reference-name scan. Other genome implementations retain their existing protocol; complete catalog outputs match. |
| 8 | Accelerate trait-table type reconstruction or query planning | Reject after inspection. The indexed wanted-table CROSS JOIN already exists; source/index hashes and mutation fences remain necessary. |
| 9 | Stream database TSV chunks or batch inserts further | Reject after inspection. Small-file parsing, uniform-type spooling, executemany, deferred indexes and bounded workers already exist. Promotion reparsing preserves lexical identifiers such as `001` and nullable Boolean values. |
| 10 | Accelerate genome sequence digests with a byte parser | Reject after inspection. Existing 64 KiB batches already use native upper/translate/hash operations. Stream/gzip reads and stat fences remain necessary; byte handling can change Unicode uppercase semantics, including `ß` and `ſ`. |
| 11 | Reuse duplicate candidate taxon normalization once per branch | Reject at screening. The retained complete context-attachment pair regresses; no stable benefit is established. |
| 12 | Decode GFF exceptions once and skip source-overlap parsing when its key is absent | Accept for block interpretation. Ordinary and metadata-rich helper profiles improve in five-pair confirmation; complete GFF main has no stable speed gain. |

Private prototypes, fixtures and raw records remain under
`tmp/performance-round2-20261025/` in the review worktree. Inspection-only
rejections carry no measured speed claim.

The [measurement record](benchmarks/workflow-speed-round2-2026-10-25.json)
retains the paired measurements, CPU times, worker peak RSS and source hashes.

## Measurement conditions

- Source baseline: `0ede811852306f29e65a6b03b6431a3d9cfb00a1` (0.8.411), the
  published result of the first review. Support sources are saved without
  switching branches. Before delivery, concurrent RAM changes were integrated
  from `c51152a52b7fbd87b4af0ac115978ff03f165fe3` (0.8.412). Final tree timings
  below compare that revision with the combined file; the other five measured
  production modules were unchanged upstream. Initial tree measurements are
  retained in the record alongside the integration comparisons.
- Freshness-checked GeneGalleon Docker image
  `local/genegalleon:performance-20261025`, immutable ID
  `sha256:4ed6901c185a713a4e73e1c8dd1a835fcf88af3fba88b7e25abede7d64c0b899`;
  runtime-input identity
  `77ee9e8c1f7e8f45ac9863db5757afc0488079b9f9fea199f71cbf84ead121dc`.
- Linux aarch64, Linuxkit 7.0.12, glibc 2.39, Python 3.12.15, pandas 3.0.6
  and pysam 0.23.3.
  Docker Desktop provides twelve CPUs and about 7.65 GiB RAM. Workers and their
  descendants use CPU 11; numerical-library worker counts are one where relevant.
- Screening uses one measured pair after one discarded warmup pair. Confirmation
  uses one discarded warmup pair and five alternating fresh-process pairs.
  Edge/error corpora run separately and are not repeated speed estimates.
- Fixtures, imports and complete output verification are outside operation
  timers unless the harness explicitly says otherwise. Wall time, process CPU
  time and Linux worker peak RSS are retained. RSS includes imports and retained
  inputs; it is not whole-workflow memory use.
- This review serializes its benchmarks and tests. Other Docker jobs were
  observed on the shared host. CPU affinity does not isolate shared storage or
  host-wide memory pressure, so noisy controls receive no speed claim.

Workloads are synthetic and use warm caches. Results do not establish NAS
throughput, whole-analysis speedups or SIF compatibility. GFF planning profiles
execute the real `main` prefix through task/source assignment, then return the
entire task list at an AST cutoff before `frames = []`; GFF processing and output
publication are excluded. Complete CLI profiles execute the unmodified `main`
and compare complete TSVs, parsed frame types, logs and errors.

## Confirmed operations and controls

Values are five-pair medians; ranges are paired baseline/candidate wall-time
ratios. Complete orthogroup CLI timings include interpreter startup, imports,
all statistics and both TSV writes; helper timings exclude startup and imports.

| Operation / synthetic input | Before | After | Paired speedup |
| --- | ---: | ---: | ---: |
| BUSCO loader, 50,000 genes with mixed hit counts | 0.20019 s | 0.11630 s | 1.66–1.77× |
| BUSCO loader, 10,000 genes × two hits | 0.40524 s | 0.08395 s | 4.71–4.99× |
| BUSCO loader, 10,000 genes × eight hits | 0.55688 s | 0.18402 s | 2.94–3.07× |
| Merge main using the two-hit BUSCO input | 0.43245 s | 0.10281 s | 4.00–4.21× |
| Tree relations, balanced 2,048-leaf tree | 0.02492 s | 0.01391 s | 1.72–1.83× |
| Tree relations, ladder 512-leaf tree | 0.10095 s | 0.02902 s | 3.21–3.61× |
| Tree relations, ladder 1,536-leaf tree | 0.89808 s | 0.25350 s | 3.47–3.66× |
| Complete orthogroup CLI, balanced 512-leaf tree | 0.60428 s | 0.59971 s | 0.94–1.08× |
| Complete orthogroup CLI, ladder 512-leaf tree | 0.88210 s | 0.78913 s | 1.10–1.15× |
| Complete orthogroup CLI, ladder 1,536-leaf tree | 2.50394 s | 1.80091 s | 1.23–1.49× |
| GFF task planning, 128 species × 2,000 IDs | 1.79815 s | 0.73284 s | 2.41–2.53× |
| GFF task planning, 512 species × 1,000 IDs | 10.32837 s | 1.49947 s | 6.61–7.18× |
| WGD sparse, 200,000 positions / 400,000 anchors | 0.94072 s | 0.47021 s | 1.90–2.03× |
| WGD dense, 20,000 positions / 40,000 anchors | 0.12990 s | 0.10603 s | 1.18–1.42× |
| WGD wide, 100,000 positions / 200,000 anchors | 1.26543 s | 0.82098 s | 1.52–1.55× |
| CDS reconstruction, 10,000 contigs / 20,000 models | 0.80877 s | 0.17880 s | 4.15–5.02× |
| Complete catalog, 10,000 contigs / 3,000 models | 1.14847 s | 1.06641 s | 1.01–1.14× |
| GFF block interpretation, 20,000 genes × eight rows | 0.29818 s | 0.21562 s | 1.33–1.51× |
| Same GFF geometry with twenty metadata fields per row | 1.04055 s | 0.63728 s | 1.47–1.67× |

BUSCO singleton-only input is 0.07883→0.07909 s (0.96–1.02×), so no gain is
claimed there. Tree controls include small, single-tip, duplicate-name and star
topologies. The initial star median is 0.00603→0.00716 s (0.53–1.06×), while
the integration comparison is 0.00599→0.00538 s (1.02–1.15×); no general star
improvement is claimed. Complete balanced-tree CLI pairs also overlap unity.
Full GFF multispecies CLI is 0.46164→0.42070 s with
ratios 0.94–1.13×; single-species planning is 0.27103→0.27994 s (0.67–1.00×).
These do not establish a stable whole-command or single-species improvement.

The native-genome small helper consistently regresses from 0.000406 to
0.000430 s (0.93–0.98×), about 24 microseconds; chromosomal and exon-rich helper
controls have overlapping pairs. Complete chromosomal catalog construction is
0.93289→0.93041 s (0.97–1.00×), so its gain is scoped to fragmented references.
GFF block parsing also improves trans-splicing and source-overlap helper profiles,
but complete GFF main is 0.70459→0.68998 s (0.92–1.04×); no whole-command
improvement is claimed for that proposal.

GFF planning worker peak RSS drops from 134,396 to 118,924 KiB for 128 species
and 179,336 to 144,976 KiB for 512 species. These include retained IDs/task lists.
The BUSCO two-hit loader adds about 0.64 MiB and complete merge main about
0.76 MiB of worker peak RSS. The final 1,536-leaf ladder adds about 0.55 MiB. These
scoped tradeoffs do not establish a general memory improvement.
WGD worker peak RSS is effectively unchanged in sparse inputs (77,972→77,984
KiB), rises by 200 KiB in the dense case and drops by 144 KiB in the wide case.
Fragmented CDS reconstruction is 71,028→71,032 KiB, while complete fragmented
catalog construction is 117,976→117,824 KiB. GFF block-helper peaks differ by
12–16 KiB in ordinary/metadata-rich controls. No general memory gain is claimed
for these proposals.

The six final production source hashes match the measured implementations.
Comparisons isolate each changed module with shared dependencies held equal
within each pair; the improvements are not additive whole-workflow estimates.
Permanent harnesses also passed bounded comparisons against the final code,
including complete merge/GFF/catalog calls and the WGD error corpus.

## Screening evidence

Screening retains one pair and cannot establish a stable improvement.

| Candidate / profile | Before | After | Status |
| --- | ---: | ---: | --- |
| Presence matrix, 64 species × 512 families, raw | 0.24869 s | 0.25890 s | No meaningful gain. |
| Presence matrix, 256 species × 512 families | 0.59707 s | 0.59502 s | No meaningful gain. |
| Presence matrix, 128 species × 512 families, dense | 0.44497 s | 0.43353 s | No meaningful gain. |
| Presence matrix, 64 species × 512 families, ZIP | 0.88945 s | 0.91154 s | No meaningful gain. |
| Scaffold context, two species × 512 taxonomy genes / 528 branches | 0.82071 s | 0.88844 s | Reject; no stable benefit established. |

## Reproduce

Save the fixed baseline and use a matching runtime.

```bash
gg_speed_baseline=tmp/speed-round2-20261025-baseline
gg_speed_tree_baseline=tmp/speed-round2-20261025-tree-baseline
mkdir -p "$gg_speed_baseline" "$gg_speed_tree_baseline"
git archive 0ede811852306f29e65a6b03b6431a3d9cfb00a1 workflow/support |
  tar -x -C "$gg_speed_baseline"
git archive c51152a52b7fbd87b4af0ac115978ff03f165fe3 workflow/support |
  tar -x -C "$gg_speed_tree_baseline"
export GG_TEST_RUNTIME=docker
export GG_CONTAINER_DOCKER_IMAGE=sha256:4ed6901c185a713a4e73e1c8dd1a835fcf88af3fba88b7e25abede7d64c0b899
```

Run sequentially, without builds/tests. Choose an available CPU instead of CPU
11 on another host. Freshness stays enabled; rebuilds follow the declared moving
upstream branches.

```bash
bash workflow/tests/run_in_runtime.sh taskset -c 11 python workflow/benchmarks/benchmark_busco_annotation_groups.py \
  --baseline "$gg_speed_baseline/workflow/support/merge_cds_annotation.py" --repeats 5 \
  --output tmp/round2-busco.json
bash workflow/tests/run_in_runtime.sh taskset -c 11 python workflow/benchmarks/benchmark_tree_relations.py \
  --baseline "$gg_speed_tree_baseline/workflow/support/orthogroup_statistics.py" --repeats 5 \
  --output tmp/round2-tree-relations.json
bash workflow/tests/run_in_runtime.sh taskset -c 11 python workflow/benchmarks/benchmark_gff_species_tasks.py \
  --baseline "$gg_speed_baseline/workflow/support/gff2genestat.py" --repeats 5 \
  --output tmp/round2-gff-tasks.json
bash workflow/tests/run_in_runtime.sh taskset -c 11 python workflow/benchmarks/benchmark_wgd_evidence_rows.py \
  --baseline "$gg_speed_baseline/workflow/support/wgd_ssd.py" --repeats 5 \
  --output tmp/round2-wgd-evidence.json
bash workflow/tests/run_in_runtime.sh taskset -c 11 python workflow/benchmarks/benchmark_catalog_references.py \
  --baseline "$gg_speed_baseline/workflow/support/gene_model_catalog.py" --repeats 5 \
  --output tmp/round2-catalog-references.json
bash workflow/tests/run_in_runtime.sh taskset -c 11 python workflow/benchmarks/benchmark_gff_block_attributes.py \
  --baseline "$gg_speed_baseline/workflow/support/gff_feature_structure.py" --repeats 5 \
  --output tmp/round2-gff-block-attributes.json
```

## Validation and delivery

The delivery target is `main`, version **0.8.413**. Its fixed pre-push base is
`c51152a52b7fbd87b4af0ac115978ff03f165fe3` (0.8.412). Upstream RAM changes
and their documentation were preserved by fast-forwarding the review checkout
and reapplying its saved changes. All six final source hashes match measured code;
the measurement record also identifies the final reproduction harnesses.

Before integration, the full focused Docker selection below passed **1,672 tests
in 479.13 s, with no skips**. Its 19 Biopython warnings concern explicit tests of
dual-coding genetic codes. The development guide remains the check-selection
authority; this is an execution record.

```bash
GG_TEST_RUNTIME=docker \
GG_CONTAINER_DOCKER_IMAGE=sha256:4ed6901c185a713a4e73e1c8dd1a835fcf88af3fba88b7e25abede7d64c0b899 \
bash workflow/tests/run_in_runtime.sh python -m pytest -q \
  workflow/tests/test_merge_cds_annotation.py \
  workflow/tests/test_orthogroup_statistics_fimo.py \
  workflow/tests/test_orthogroup_statistics_support.py \
  workflow/tests/test_orthogroup_statistics_relations.py \
  workflow/tests/test_gff2genestat.py \
  workflow/tests/test_gff_translation_exceptions.py \
  workflow/tests/test_format_species_inputs.py \
  workflow/tests/test_paired_species_curation.py \
  workflow/tests/test_pairwise_synteny.py \
  workflow/tests/test_gene_model_catalog.py \
  workflow/tests/test_gene_model_catalog_phase_conventions.py \
  workflow/tests/test_gene_model_refinement.py \
  workflow/tests/test_cloned_cds_models.py \
  workflow/tests/test_wgd_ssd_contracts.py \
  workflow/tests/test_wgd_ssd_gene_stage.py \
  workflow/tests/test_generate_orthogroup_database.py \
  workflow/tests/test_gene_evolution_query_mode.py \
  workflow/tests/test_gff_repair.py -x
```

After integration, the overlapping tree implementation and upstream promoter
implementation were checked for preservation, and the tree comparisons were
repeated against the integrated baseline. The combined code passed **503 tests
in 53.84 s, with no skips**, covering the changed shared readers and their callers:

```bash
bash workflow/tests/run_in_runtime.sh python -m pytest -q \
  workflow/tests/test_orthogroup_statistics_fimo.py \
  workflow/tests/test_orthogroup_statistics_support.py \
  workflow/tests/test_orthogroup_statistics_relations.py \
  workflow/tests/test_orthogroup_statistics_promoter.py \
  workflow/tests/test_gff2genestat.py \
  workflow/tests/test_get_trait_matrix_outputs.py \
  workflow/tests/test_generate_orthogroup_database.py \
  workflow/tests/test_gene_evolution_query_mode.py \
  workflow/tests/test_gff_translation_exceptions.py -x
```

`GG_LINT_BASH="$PWD/tmp/performance-20261025/bash-parser" bash ./dev lint`
passed before and after integration: Bash syntax, Ruff, eight entrypoints and 27
common parameters. The wrapper
delegates syntax parsing to Bash in the same immutable Docker image; Ruff and
configuration checks run on the host. `git diff --check` and the staged diff
check passed. These are static evidence, separate from Docker behavior.

Benchmark preflights exposed two harness defects: unnamed internal CLI fixture
nodes triggered the baseline's existing name guard, and the promoted WGD worker
did not receive its required baseline argument. Both harnesses were corrected
and their comparisons passed; no production validation was weakened.

SIF is unavailable on this macOS host. R helpers and container inputs are
unchanged, so the substantial R/full scientific lanes were not run. Docker
results do not establish SIF compatibility or whole-analysis speedups.
