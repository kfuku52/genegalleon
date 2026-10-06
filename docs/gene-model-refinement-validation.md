# Gene-model refinement validation and measurements

Measured on 2026-10-05. This covers all-isoform admission, conserved
representative selection, existing-locus prediction, and selected downstream
views. It is implementation evidence for an experimental annotation aid;
biological precision/recall and scientific thresholds are not calibrated.

## Additional correctness audit

The implementation received a second audit of scientific eligibility, source
identity, publication dependencies and actual downstream tool inputs. The
following defects were reproduced and corrected; the final measurements below
distinguish that revision from the earlier implementation.

| Boundary | Reproduced problem | Correction |
| --- | --- | --- |
| Representative selection | Unusable donors influenced scores; a long ineligible candidate changed the quality normalization and adoption margin; longest could select forbidden predicted paths; a sole eligible repair lacked a baseline margin. | Exclude nonvoting paths from conservation scores and eligibility-based quality normalization; evaluate the source baseline on the same quality scale; enforce prediction eligibility in either policy; measure the actual source-baseline difference. |
| Component capacity and reproducibility | Invalid correspondence could exceed the component cap; edge order changed the last floating-point digit of the objective. | Determine copy ambiguity globally, use voting edges for components, and accumulate weights in canonical order. |
| Rescue followed by refinement | Refinement used the original sources and lost accepted missing models; combined arrays could freeze the plan before rescue finalized. | Verify augmented publications, retain their source inputs and genetic codes, and schedule refinement preparation after rescue finalization. |
| Partial CDS and genetic codes | Raw partial CDS was translated from an incompatible frame downstream; a selected code override could be ignored. | Publish admitted CDS/GFF with only known partial bases clipped, assert exact protein agreement, pass effective codes, and reject mixed codes for a single codon model. |
| Publication integrity | QC missed absent stages and foreign/obsolete keys; copied manifest metadata and derived SQLite/annotation content could escape verification. | Reconstruct the required dependency graph; bind exact manifest bytes and every consumed publication receipt/content before and after building. |
| Annotation identity and ownership | Ancestor axes could disagree with CDS; implicit noncoding owners were missed; substring-based prediction Parents and mixed GTF/GFF output corrupted graphs. | Verify axes, aggregate implicit ownership, emit exact IDs and normalized derived graphs, and archive original decoded bytes exactly. |
| Parallel selected analyses | Another bundle could replace a shared FASTA, synteny or search index between build and read. | Separate indices, CDS resolution, traits and their receipts by manifest SHA and sequence view; preserve legacy cache paths. |
| Query/reference translation | Protein-mode gene-ID queries and DIAMOND references globally retranslated raw partial/mixed-code CDS. | Use admitted protein inventories for selected queries and reference databases, including exclusion membership. |
| Selected TBLASTN | Nucleotide databases retained raw negative/partial paths and used one global genetic code for a mixed-code bundle. | Index admitted analysis CDS separately, apply each species' code, and use the total admitted nucleotide database size in each search and its provenance. |
| External pairwise inputs | A hash-consistent external manifest could associate an unrelated admitted protein with a selected CDS. | Compare the admitted protein with the exact phase-aware, species-code translation before preparing a pairwise input. |
| Fractionation reader | Qualified CDS FASTA IDs did not match original selected GFF gene/transcript IDs. | Build a verified canonical coordinate view from the exact representative map; exercise the real reader and official validation CLI. |
| Context-dependent genetic codes | Dual-coding stop paths could bypass native `translation_uncertain` admission and be labeled intact. | Withhold observed in-frame dual codons, including terminal ones, from comparative/protein/coding admission and phase inference; preserve DNA/source annotations and the existing nondual translation convention. |

Biological score thresholds were not weakened to pass these regressions.

## Environment and reproducibility

Docker image `local/genegalleon:cds-admission-20261004`, immutable image ID
`sha256:8a6d69ec10defff91117d59591eeaabdbd76b41218525a7c038de815af54f87e`;
Linux aarch64, Python 3.12.14, Biopython 1.88, pysam 0.23.3. Host: Apple M2 Max,
12 logical CPUs, 64 GiB RAM; Docker VM: 12 CPUs, approximately 7.65 GiB RAM.
The runtime wrapper confirmed freshness against the
current container inputs and daily owned upstream snapshot. No SIF was executed.

From the repository root, use:

```bash
export GG_TEST_RUNTIME=docker
export GG_CONTAINER_DOCKER_IMAGE=local/genegalleon:cds-admission-20261004
bash workflow/tests/run_in_runtime.sh python -m pytest -q \
  workflow/tests/test_gene_model_catalog.py \
  workflow/tests/test_gene_model_selection.py \
  workflow/tests/test_gene_model_store.py \
  workflow/tests/test_gene_model_refinement.py -x
bash workflow/tests/run_in_runtime.sh python -m pytest -q \
  workflow/tests/test_pairwise_synteny.py \
  workflow/tests/test_representative_selection.py \
  workflow/tests/test_synteny_cutoff_metadata.py \
  workflow/tests/test_gene_model_refinement_downstream.py \
  workflow/tests/test_gene_evolution_query_mode.py \
  workflow/tests/test_gene_model_refinement_runtime.py -x
bash workflow/tests/run_in_runtime.sh env KFFRACTBIAS_RUN_INTEGRATION=1 \
  python -m pytest -q workflow/tests/test_fractionation_bias_integration.py -x
bash ./dev check static
bash ./dev check smoke
```

Audit validation used the same immutable Docker runtime. Counts below describe
separate, overlapping runs and must not be added.

| Check | Result | Wall time / revision |
| --- | --- | --- |
| Catalog, selector, SQLite store and refinement (four complete files) | 217 passed | 71.31 s; final eligibility normalization and dual-codon guard |
| Real refinement/rescue predictor, selected downstream and canonical fractionation reader (three files) | 50 passed | 99.04 s; after final catalog guard; before the subsequent TBLASTN fix |
| Wider selected readers, CDS resolution, synteny, array/tooling, caches and gene-ID query mode (13 files) | 479 passed | 149.92 s; before the final dual-codon guard and TBLASTN fix |
| Final pairwise/maps/synteny/selected queries and real refinement runtime (six complete files) | 136 passed | 79.19 s; after the TBLASTN and external protein-consistency fixes |
| Existing input-generation single/parallel-array/resume suite | 35 passed | 314.24 s |
| Existing genome-evolution protein-mode suite | 80 passed | 608.23 s |
| Real fractionation-bias integration, including external scratch and self mode | 5 passed; no skips | 31.43 s; `KFFRACTBIAS_RUN_INTEGRATION=1` inside the container |
| Static lane | 294 passed | 3.29 s; final source |
| Smoke lane | 5 passed; 35 deliberately deselected | 7.27 s |

The four-file primary check retains 19 expected Biopython warnings for
context-dependent genetic-code fixtures; their paths are explicitly withheld.
The broader legacy formatter GFF/GTF/attribute selection also passed 95 cases
(208 deliberately deselected) during initial implementation, before this audit.

The real predictor tests verify exact restoration of held-out exons on both
strands, addition of a complete coding RNA path, homology-only versus RNA-supported
adoption, actual JCVI synteny reuse for an existing anchor-excluded locus, and
preservation of genuine loss/internal-stop/assembly-gap/pseudogene negatives.
Other cases cover source mismatch exclusions and snapshots, GTF graph export,
copy ambiguity, noncoding ownership, conflicting predictions, source/RNA/tool
mutation, corrupt/stale receipts, strict transcript maps, cache provenance,
selected CDS/protein/GFF agreement and scheduler dependencies. Small exact
selection oracles and star/weighted graphs check the bounded heuristic and
its degree-independent adoption margins.

Host Ruff and configuration checks passed; Docker Bash syntax passed for 48
workflow/dev files. The combined `dev lint` entrypoint could not complete because
the host has Bash 3.2 and the image lacks Ruff. Its Bash/Ruff/configuration
components were executed separately. The full unrelated runtime/R/download
suite and a real Slurm cluster were not executed; scheduler submissions were
verified using isolated fixtures. No curated inputs were regenerated.

## Existing annotation admission

The [final-source read-only admission audit](assets/benchmarks/gene-model-refinement/admission-post-audit-final.json)
records all source hashes before/after and mismatch examples. Reproduce it with
an output filename in a new temporary/output directory:

```bash
bash workflow/tests/run_in_runtime.sh python \
  docs/assets/benchmarks/gene-model-refinement/audit_catalog.py \
  --input workspace/input --output /tmp/gg-admission-new.json \
  --runtime-image sha256:8a6d69ec10defff91117d59591eeaabdbd76b41218525a7c038de815af54f87e
```

The seven curated species contain 2,160 source CDS records and loci, all uniquely
mapped to source representatives. The catalog preserves **2,685 coding
candidates**: **2,572 usable** and **113 withheld**. It identifies 465 source-bound
complete ORFs with uniquely inferred missing phase and 814 exact preceding-CDS
matches with the established terminal-stop-to-`NNN` convention.

There are 87 source/genomic mismatch records: 85 have equal lengths with source
masking, and two have length/coordinate disagreement. Nine Cephalotus masking
cases conceal genuine internal genomic stops; Dionaea cases include IUPAC
ambiguity and the two unresolved length differences. They remain withheld.
All **21 source files** have unchanged hashes. This admission audit does not
measure repair accuracy or isolated performance.

## Selector performance

The table in this section records the correctness-audit revision before the
subsequent optimization below; its implementation hashes remain in each report.

Generated workloads contain nine species, four isoforms per locus, 300-aa core
proteins and a sparse degree-four correspondence graph. One warmup precedes
three measured fresh-process repetitions, except the explicitly labeled single
capacity observations. Output equivalence includes choices, scores and audit
results, excluding performance counters. The timed region is the selection
call; interpreter startup and generated-input construction are excluded.
Process RSS includes imports and the loaded workload. The complete-pipeline
measurements below include startup and input verification.

| Workload / implementation | Wall time | Peak RSS | Evidence |
| --- | --- | --- | --- |
| 32 families / uncached reference | median 35.852 s; 35.713–38.235 s | 50.45 MB | [post-audit paired measurement](assets/benchmarks/gene-model-refinement/selector-cache-post-audit-32.json) |
| 32 families / cached selector | median 1.650 s; 1.639–1.690 s | 50.52 MB | Same input/output hashes as reference; 21.73× |
| 32 families / SQLite components | median 1.583 s; 1.583–1.589 s | 48.44 MB | Same scientific output as both reference and cached selector |
| 256 families, 2,304 loci / cached in-memory | median 21.328 s; 19.972–21.352 s | 89.99 MB | [post-audit paired store measurement](assets/benchmarks/gene-model-refinement/selector-store-post-audit-256.json) |
| 256 families / SQLite components | median 14.935 s; 12.070–19.368 s | 58.77 MB | Same input and scientific output as cached in-memory; three fresh-process repeats each |
| 1,000 families, 9,000 loci / in-memory | 88.498 s, one observation | 226.84 MB | [capacity observation](assets/benchmarks/gene-model-refinement/selector-memory-1000.json) |
| 1,000 families / SQLite components | median 51.691 s; 50.409–83.467 s | 85.04 MB | [repeated store measurement](assets/benchmarks/gene-model-refinement/selector-store-1000.json); same scientific output |

MB is decimal. Store measurements use the fresh worker's `/proc/self/status`
VmHWM, excluding the driver's input preparation. Final 256-family SQLite
preparation took 0.345 s, with a 20.74 MB database. The historical 1,000-family
preparation took 1.28 s with an 80.96 MB database. Sequences were loaded nine
loci at a time; the component pair cache peaked at 85 entries. Global identifier,
edge and result metadata still grows with input size.

The final 256-family SQLite timing and historical 1,000-family timing have
substantial variation. The paired 256-family result establishes equivalent
output and lower peak RSS (89.99 to 58.77 MB); it is not a precise runtime
speedup estimate. The 1,000-family observations support capacity and memory
feasibility rather than a precise speedup over a single in-memory observation.
The cache comparison is an equivalent-work/output optimization. Longest-CDS
selection has different behavior and is not an equivalent performance baseline.
The final paired measurement was repeated after all correctness fixes,
including eligibility-based quality normalization. All nine measured
reference/cached/SQLite runs retained the regular graph's earlier exact
scientific hash; implementation hashes and runtime freshness evidence are in
the report. Earlier measurements remain historical evidence rather than
before/after timing baselines for this revision. The 1,000-family capacity
observations precede the final normalization; they quantify storage and
component loading, with no precise final-version runtime extrapolation.

```bash
bash workflow/tests/run_in_runtime.sh python workflow/tests/benchmark_gene_model_selection.py \
  --families 32 --species 9 --length 300 --implementations reference cached store \
  --repeats 3 --warmups 1
bash workflow/tests/run_in_runtime.sh python workflow/tests/benchmark_gene_model_selection.py \
  --families 256 --species 9 --length 300 --implementations cached store \
  --repeats 3 --warmups 1
```

## Subsequent optimization and equivalent-output checks

Profiling the audited selector identified repeated eligibility checks and
whole-graph reconstruction for each small component. Eligibility is now
compiled once per invocation; trial states and voting weights cover only their
component. A second profile identified pairwise scoring overhead. Each score
cache now reuses one lazily constructed aligner, projects only requested splice
boundaries, and computes an identical-protein result directly after applying
the existing alignment-cell cap. Coordinate updates reuse the scores already
computed for that node's unchanged neighborhood. Thresholds, gap scores,
floating-point accumulation order, synchronous adoption and capacity limits
are retained.

The direct pre-optimization native alignment/full-residue-mapping oracle was
compared using exact JSON: 1,922 exhaustive short-sequence/cap combinations,
500 randomized directional comparisons, partial/minus/internal-stop cases,
terminal splice phases, and all 128 ASCII characters. Weighted irregular
graphs, excluded donors, candidate/component caps, LRU eviction and invocation
boundaries also passed. The optimized nine-file Docker check completed with
**444 passed, no skips, 19 expected dual-codon warnings, in 161.57 s**:

```bash
bash workflow/tests/run_in_runtime.sh python -m pytest -q \
  workflow/tests/test_gene_model_catalog.py \
  workflow/tests/test_gene_model_selection.py \
  workflow/tests/test_gene_model_store.py \
  workflow/tests/test_gene_model_refinement.py \
  workflow/tests/test_performance_helpers.py \
  workflow/tests/test_representative_selection.py \
  workflow/tests/test_gene_model_refinement_runtime.py \
  workflow/tests/test_gene_model_refinement_downstream.py \
  workflow/tests/test_pairwise_synteny.py -x
```

These checks execute real miniprot/JCVI and DIAMOND/TBLASTN inputs and warm
database reuse. Static checks passed separately: `dev check static` **294 passed
in 2.88 s** after the final provenance-label correction, Ruff over `workflow` and `container/scripts`, configuration checks
for eight entrypoints/24 parameters, and Docker Bash syntax for 75 tracked
shell/dev files. Host `dev lint` remains blocked by Bash 3.2; its components
were executed separately. No SIF was used.

Input verification now uses the existing `digest_paths` boundary to freshly
hash every unique publication/role target once. The full receipt, copied TSV
binding and external role paths are retained. Parsing is fenced before and
after the batch, and changes to previously hashed files are rejected. No
validated content hash is reused across invocations or stage boundaries.

Selected TBLASTN nucleotide totals are recorded when FASTA blobs are built,
using the same decompressed UTF-8 byte semantics as the previous full scan.
The private aggregate is checked against database/source identity and record
counts. Its reader verifies before/after identities without writing the index.
Existing schema-2 indices receive one locked atomic rebuild; open readers keep
their original inode. Public manifest bytes and extracted FASTA remain exact.
Malformed/missing metadata, database/manifest/source mutation and portable
copies are covered by executed regressions.
The gene-evolution provenance parameter now declares the actual FASTA-store
schema 2; its old literal schema 1 label was stale. Manifest bytes are unchanged.
Correcting this label invalidates older query/primary artifact receipts once;
their reusable sequence indices are retained. The focused native DIAMOND/TBLASTN
check, including warm reuse and schema provenance, passed **2 tests in 7.26 s**
after this correction.

All optimization measurements used the same immutable runtime, one warmup
and three fresh-process repetitions per implementation. Our benchmark/test
jobs ran serially. Filesystem caches were warm; no cache flush was performed.
Before/after scientific inputs, outputs, and implementation hashes are recorded
in the linked reports.

| Operation / workload | Before median (range) | Final median (range) | Peak RSS before → final | Evidence |
| --- | --- | --- | --- | --- |
| In-memory selector, 256 families × 9 species | 19.375 s (18.231–19.747) | 10.938 s (10.825–11.308), 1.77× | 90.09 → 90.03 MB | [before](assets/benchmarks/gene-model-refinement/selector-optimization-before-256.json), [final](assets/benchmarks/gene-model-refinement/selector-optimization-final-256.json) |
| SQLite selector, same 2,304 loci | 11.939 s (11.849–18.375) | 11.169 s (11.067–11.199), 1.07× | 58.78 → 58.30 MB | Same selector reports |
| Published-input verification, 3 × 64 MiB genome FASTA | 0.595 s (0.593–0.603) | 0.443 s (0.432–0.503), 1.34× | 61.20 → 61.29 MB | [before](assets/benchmarks/gene-model-refinement/representative-inputs-optimization-before.json), [after](assets/benchmarks/gene-model-refinement/representative-inputs-optimization-after.json) |
| FASTA byte-total retrieval, 54 MB / 18,000 records | 185.16 ms (182.98–195.08) | 2.21 ms (1.13–4.38), approximately 80× | 25.66 → 24.86 MB | [before](assets/benchmarks/gene-model-refinement/fasta-totals-optimization-before.json), [final](assets/benchmarks/gene-model-refinement/fasta-totals-optimization-final.json) |

Selector equivalence covers every choice, score and audit field, with exact
scientific SHA-256 `7a6cbee4005b48c4b4716091e18a0d7b37009259f29373351de3dbad8ebff0f8`.
Pair requests fell from 632,832 to 503,808; unique score evaluations remained
21,760. SQLite still loads at most nine loci and 85 cached pairs for this
workload. Its modest speed improvement is much smaller than the in-memory
improvement; the earlier before run also contained a slow outlier.

The [first component-only measurement](assets/benchmarks/gene-model-refinement/selector-optimization-after-256.json)
gave 12.509 s in memory and 11.934 s in SQLite. Only the in-memory path improved
at that stage. [Before profiling](assets/benchmarks/gene-model-refinement/selector-optimization-profile-before.txt)
and [component-only profiling](assets/benchmarks/gene-model-refinement/selector-optimization-profile-component.txt)
showed approximately 86.6 million → 17.0 million calls and then native alignment
as the dominant remaining cost. Profiler times are diagnostic and are excluded
from the timing baselines.

The verification workload is a real off/longest publication of the three-species,
seven-locus truth fixture. Each genome has an appended unannotated contig with
80-base lines; the original genome prefix, CDS and GFF remain byte-identical.
Its timing includes actual `verify-inputs --field layout` imports and all content
checks. Exact layout, manifest and scientific bundle hashes remain unchanged;
fresh worker wall including interpreter/instrumentation startup is 0.655 →
0.499 s. Fixture creation/publication is excluded. The initially invalid padding
fixture failed faidx and was preserved; no measurement from it is used.

FASTA total timing excludes imports. In the final paired run the unchanged
legacy scan took 230.36 ms (220.16–231.18); the aggregate took 2.21 ms, roughly
two orders of magnitude faster. Millisecond wall latency varies, so this is
an approximate scale rather than an exact predicted speed ratio. Fresh worker
wall including startup was 278.61 → 54.12 ms in that pair. Total, source identity,
public manifest and sampled FASTA hashes match the original before workload.
The [first after report](assets/benchmarks/gene-model-refinement/fasta-totals-optimization-after.json)
records the one-time atomic migration at 3.38 s; later `ensure` was unchanged.
The old cold build was 2.89 s. Migration needs temporary build space and occurs
once for each pre-aggregate index; its cost is excluded from warm retrieval.

Reproduction uses retained benchmark sources outside curated `workspace/input`:

```bash
bash workflow/tests/run_in_runtime.sh python workflow/tests/benchmark_gene_model_selection.py \
  --families 256 --species 9 --length 300 --implementations cached store \
  --repeats 3 --warmups 1 \
  --expected-output-sha256 7a6cbee4005b48c4b4716091e18a0d7b37009259f29373351de3dbad8ebff0f8
bash workflow/tests/run_in_runtime.sh python workflow/tests/benchmark_representative_inputs.py \
  --genome-mib 64 --repeats 3 --warmups 1 \
  --work-dir workspace/output/gene_model_refinement_performance/representative-inputs-v2 \
  --expected-report docs/assets/benchmarks/gene-model-refinement/representative-inputs-optimization-before.json \
  --output workspace/output/gene_model_refinement_performance/verification-replay.json \
  --runtime-image-label local/genegalleon:cds-admission-20261004 \
  --runtime-image-id sha256:8a6d69ec10defff91117d59591eeaabdbd76b41218525a7c038de815af54f87e
bash workflow/tests/run_in_runtime.sh python workflow/tests/benchmark_fasta_sequence_totals.py \
  --work-dir workspace/output/gene_model_refinement_performance/fasta-totals \
  --implementations legacy aggregate \
  --expected-result docs/assets/benchmarks/gene-model-refinement/fasta-totals-optimization-before.json
```

For a fresh checkout, choose new empty work directories and omit expected
reports on the preparation run; capture clean JSON reports before changing
the implementation. Retained indices include private absolute paths and are
validated through `ensure` when moved.

## Complete pipeline and restart behavior

The original table below is historical. The
[optimized complete-pipeline report](assets/benchmarks/gene-model-refinement/complete-pipeline-optimization-final.json)
repeats the same fixture/styles/policies three times: conserved/conservative
batched cold 1.381 s (1.360–1.386), batched resume 0.462 s (0.456–0.664),
staged cold 3.242 s (3.207–3.557), and staged resume 2.242 s (2.237–2.409).
Every scientific output hash matches the historical report, including both
repaired representatives, the accepted exon-skipped path and all four negatives.

Because the first staged observation was slower than historical measurements,
the exact old main/selector/FASTA helpers were restored in an isolated tree and
SHA-256-gated against the before reports. Current and old frozen trees then ran
in alternating order (before/after, after/before, before/after) in one Docker
runtime. The [paired summary and six raw reports](assets/benchmarks/gene-model-refinement/complete-pipeline-paired-summary.json)
retain exact scientific outputs and show substantial launch/predictor variation:

| Conserved/conservative paired run | Before median (range) | After median (range) |
| --- | --- | --- |
| Batched cold | 1.591 s (1.358–1.936) | 1.369 s (1.354–1.487) |
| Batched resume | 0.541 s (0.483–0.586) | 0.455 s (0.452–0.490) |
| Staged cold | 3.740 s (3.001–4.103) | 2.991 s (2.943–4.135) |
| Staged resume | 2.527 s (2.064–2.564) | 2.266 s (2.068–2.728) |

These overlapping ranges do not support a precise whole-pipeline speedup or
regression claim. The operation-level comparisons above are the quantified
optimization evidence. Native prediction and repeated interpreter startup
remain important for small fixtures. Maximum individual-process RSS remained
approximately 68.06 MiB in the paired conserved cold runs. Baseline snapshots
are retained under `workspace/output/gene_model_refinement_performance/baseline-recheck`;
the active implementation was never replaced during this audit.

The [post-audit complete pipeline report](assets/benchmarks/gene-model-refinement/complete-pipeline-post-audit.json)
records three repeats for each execution style/cache state/policy. The generated
truth fixture has three species × seven
loci, 280-aa proteins, two held-out exon models, a reviewed whole coding RNA path,
and four intact/biological-negative loci. miniprot 0.18-r281 runs at two CPUs.
This fixture uses explicit frozen correspondence; it does not time discovery of
new synteny comparisons. Actual synteny reuse is covered by the runtime tests.

| Conserved/conservative execution | Wall time median (range) | CPU median | Maximum individual process RSS |
| --- | --- | --- | --- |
| One `run` invocation, new output | 1.341 s (1.324–1.353) | 1.374 s | 68.05 MiB |
| One `run` invocation, completed-output resume | 0.464 s (0.438–0.472) | 0.405 s | 57.72 MiB |
| Separate stage CLI calls, new output | 2.903 s (2.865–2.904) | 2.786 s | 68.05 MiB |
| Separate stage CLI calls, completed-output resume | 1.996 s (1.995–2.002) | 1.799 s | 57.71 MiB |

Cold means a new workflow output directory; OS caches were not flushed. Wall
includes fresh interpreter startup, plan/catalog/correspondence/selection/
prediction/finalization and source/receipt verification. Docker startup and
fixture generation are excluded. CPU includes completed predictor descendants;
RSS is the largest individual CLI/child process, not simultaneous aggregate
memory. Separate stages are invoked sequentially, without Slurm submission or
queue time. For this small fixture, interpreter and validation overhead matter.

Separate-stage cold medians: plan 0.339 s, catalog 0.519 s, correspondence
0.411 s, selection 0.402 s, prediction 0.730 s, finalization 0.492 s. Stage medians
are descriptive and need not sum to the median total.

Every same-policy repeat/style/restart has identical scientific output hashes
and unchanged source bytes. Conserved/conservative exactly restores **2/2**
held-out exon models and adopts both repaired representatives, accepts the exact
**1/1** missing exon-skipped path, and changes **0/4** intact/negative loci.
The added path remains unselected in this ambiguous fixture; acceptance does not
force representative adoption. Longest/off takes a batched cold median 1.212 s
and makes no repairs/additions; it is a different scientific policy, not an
output-equivalent speed baseline.

```bash
bash workflow/tests/run_in_runtime.sh python workflow/tests/benchmark_gene_model_refinement.py \
  --repeats 3 --cpus 2 \
  --output workspace/output/gene_model_refinement_performance/complete-pipeline-replay.json
```

## Delivery integration (0.8.148)

Before delivery, the checkout was fast-forwarded to `origin/main` at `f296ae7`.
The upstream coding-only GFF reader and alignment-ancestor closure were retained
alongside exact representative-map selection. Runtime registration retains all
upstream checks. Refinement imports the native synteny reader through rescue,
so its contract tests now run in the container runtime lane. Its CLI help check
also runs there, with real dependencies and an empty temporary working directory;
it is excluded from the fast environment that does not install kfFractBias.

All checks below used the same fresh immutable Docker image documented above.
These are separate runs; their overlapping counts must not be added.

| Delivery check | Result | Wall time |
| --- | --- | --- |
| Refinement, selected readers, arrays, synteny and upstream integration (17 files) | 761 passed; no skips; 19 expected dual-codon warnings | 225.26 s |
| Static lane | 294 passed | 3.77 s |
| Smoke lane | 5 passed; 35 deliberately deselected | 6.39 s |
| Fast routing and support CLI checks | 143 passed | 17.67 s |
| Strict runtime refinement CLI help | 1 passed; 106 deliberately deselected | 0.99 s |

```bash
bash workflow/tests/run_in_runtime.sh python -m pytest -q \
  workflow/tests/test_gene_model_catalog.py \
  workflow/tests/test_gene_model_selection.py \
  workflow/tests/test_gene_model_store.py \
  workflow/tests/test_gene_model_refinement.py \
  workflow/tests/test_performance_helpers.py \
  workflow/tests/test_representative_selection.py \
  workflow/tests/test_gene_model_refinement_runtime.py \
  workflow/tests/test_gene_model_refinement_downstream.py \
  workflow/tests/test_pairwise_synteny.py \
  workflow/tests/test_fractionation_selected_inputs.py \
  workflow/tests/test_cds_resolution.py \
  workflow/tests/test_gff2genestat.py \
  workflow/tests/test_input_generation_array_scripts.py \
  workflow/tests/test_synteny_cutoff_metadata.py \
  workflow/tests/test_synteny_search_integration.py \
  workflow/tests/test_input_recovery_regressions.py \
  workflow/tests/test_rescue_model_evidence.py -x
bash ./dev check fast workflow/tests/test_validation_runner.py \
  workflow/tests/test_support_script_help_smoke.py -x
bash workflow/tests/run_in_runtime.sh python -m pytest -q \
  --gg-strict-runtime --gg-suite=runtime \
  workflow/tests/test_gene_model_refinement.py \
  workflow/tests/test_gene_model_refinement_runtime.py \
  -k 'cli_help or longest_only_mode' -x
```

Ruff over `workflow`/`container/scripts` and the saved admission-audit script,
configuration checks (eight entrypoints, 24 common parameters), Docker Bash
syntax (75 tracked shell/dev files), staged whitespace checks, and local
documentation/evidence paths passed. Host `dev lint` is blocked by Bash 3.2;
its components were checked separately. SIF, the full unrelated runtime/R/
download suite, and real Slurm submission remain unverified in this delivery.
Performance reports retain their measured source hashes and historical
revisions; they were not regenerated after the upstream integration.

## Native real-input measurements (0.8.169, 2026-10-06)

These additional measurements use the fresh SIF runtime on audrey1 with QNAP
inputs. They are separate from the historical Docker measurements above.
The baseline modules were saved from 0.8.168 before editing. Input, module,
binary and output hashes accompany the native measurement records.

For the 1,187,329,241-base, 131-contig Ancistrocladus assembly, three runs of
the previous QNAP copy/index path took 34.01, 22.82 and 22.47 seconds.
Verified local staging took 17.40 seconds cold; three warm uses took 3.63,
3.62 and 3.62 seconds. Warm medians differ by about 6.3-fold. Every run produced
the same decompressed FASTA SHA-256. The optimized timings include complete
source and cached FASTA/index verification. This measures genome staging,
not whole-workflow speed. A local cache needs storage for decompressed genomes;
final self-contained publication copies are still required.

The prior rescue publication contains 4,960 loci and 3,661 distinct proteins.
Swiss-Prot baseline search/publication took 76.28 seconds cold and 25.38 seconds
warm. An initial split-cache implementation regressed to about 74–77 seconds
warm: 27,882 small SQLite calls incurred shared-disk locking overhead.
Batched accession/protein lookups removed that overhead. Three warm runs took
22.44, 18.81 and 23.43 seconds; primary support classifications agreed exactly
with the baseline. Changing the minimum aligned length from 50 to 55 took
27.05 seconds and searched zero proteins, reusing all 3,661 raw alignments.
This demonstrates threshold-independent search reuse. It does not establish
hundreds-of-species end-to-end speed or a cold-cache speed improvement.

Known-annotation holdouts used official NCBI reference genome/CDS/GFF bundles
for C. elegans WBcel235 (GCF_000002985.6) and S. cerevisiae R64
(GCF_000146045.2). Download MD5s and frozen input SHA-256s were verified.
Twenty deterministic single/multiple-exon-stratified loci per organism were
chosen from exact intact CDS/genome matches with 100–600-residue proteins.
The default identity 0.5 and coverage 0.95 gates restored 18/20 worm and 14/20
yeast coding paths exactly. One additional worm path was admitted inexactly:
miniprot retained a short annotated intron as coding sequence. Both organisms
had zero additions in 20 deliberate assembly-gap controls each.

An explicit benchmark profile with identity 0.98 rejected that inexact worm
path while retaining the same 18 and 14 exact restorations and zero gap-control
additions. This exact-source calibration is not a recommended threshold for
divergent interspecies donors and does not change global defaults. The misses
include very short initial coding exons; invalid start/stop paths remain
withheld. These tests evaluate local prediction/admission against annotation,
not independent whole-cohort discovery precision, functional-gene validity,
or independent interspecies/RNA evidence. Short-protein Swiss-Prot support
thresholds also remain unchanged and their failed length/coverage reasons
are reported separately.

Reproduce the holdout component with
`workflow/support/benchmark_gene_model_holdout.py --inputs FILE --output DIR`;
optional `--species-profiles FILE` freezes explicitly chosen target thresholds.
Run through `workflow/tests/run_in_runtime.sh` with the intended runtime.

## General limits

These workloads do not establish whole-genome or hundreds-of-species runtime,
annotation accuracy, or SIF compatibility. Local windows and prediction CPU
limits bound predictor work; nomination rate, genome I/O, family size and protein
length affect cost. Complete effective genome/source snapshots also require
storage proportional to inputs. Fresh downstream invocations verify complete
publication hashes, including genome/source snapshots; repeated verification
can dominate I/O for very large bundles. The optimization above separately
measures a 192 MiB synthetic genome payload as well as the small complete
pipeline; it does not establish multi-gigabase or hundreds-of-species runtime.
Components larger than 5,000 loci abstain and
per-locus candidate/alignment/cache limits prevent unbounded searches. A bounded
heuristic may miss a global optimum. Donor sampling is not clade balanced;
true lineage-specific isoforms and one-to-many correspondence need biological
review. The release consumes reviewed whole coding RNA paths, not raw BAM files.
Translation uses the existing linear codon-table convention; it does not model
complete-CDS alternative initiators (TTG/GTG retain L/V). Context-dependent dual
stops are explicitly withheld rather than interpreted as complete-CDS stops.
Selected TBLASTN searches retain the declared E-value cutoff and use the total
admitted nucleotide length via `-dbsize`, avoiding artificially small database
sizes. Separate species searches still have BLAST's length-adjustment and
subject-count heuristics and a 50,000-subject limit per search. They do not
promise bit-identical E-values to a concatenated single-code database. External
pairwise manifests declare their own admission membership; sequence agreement
is checked, but these reviewed inputs do not acquire the producer's complete
scientific admission audit merely by supplying hashes.
