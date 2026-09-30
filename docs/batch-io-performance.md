# Batched verification, audited databases and PDF rendering

The single-family API, SQLite format, type promotion rules, plot arguments and
scheduler resources remain compatible. Content reuse lasts one operation only;
no saved hash or rendering receipt is accepted as scientific completion evidence.

## Indexed query ownership

Bulk query2family summary, storage conversion/materialization and provenance
audit now build one operation-scoped query-ID priority index. Only complete
filenames and underscore/dot boundaries are looked up; each file no longer scans
every query ID. The scalar matcher and catalog-list API remain available.
Overlapping IDs, arbitrary matcher priority, duplicates, Unicode and empty IDs
retain scalar behavior. The index snapshots its catalog and caches no filenames.

Linux arm64 Docker `local/genegalleon:dev`, Python 3.12.14, base `40b78b4`, one
warmup and three measured fresh processes: 3,003 IDs / 9,006 filename matches
took 1.929 s to 0.0102 s (190×); a complete three-directory, 3,003-row summary
took 1.502 s to 0.0619 s (24.2×). Result lists and complete TSV bytes match.
Peak RSS increased about 0.4 MiB (about 75 MiB matching / 77 MiB summary).
These are query-ownership workloads, not overall scientific-pipeline speedups.

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_query_matching.py \
  --output /tmp/query-matching.json
```

Use `--support-root` with a complete baseline support tree for comparison.

## Family presence/absence tables

Long-table construction reads each unique matrix column once and reuses species
display names within one operation. Family-major order, requested species order,
missing values and column types remain unchanged. Ambiguous axes and missing
labels retain the original scalar lookup behavior and first error.

Linux arm64 Docker, Python 3.12.14, pandas 3.0.6, base `a41fa93`, same immutable
image, one warmup and three measured fresh processes: 2,048 families × 64 species
produced 131,072 rows in 2.642 to 0.405 s (6.52×). Every seventh family lacks its
branch output and retains missing copy/presence values. Complete table values,
columns/types/order and source matrices match. Median process peak RSS fell
from 228.8 to 212.7 MiB (7.1% less). Timing excludes fixture construction and
fingerprinting; peak RSS includes both. This measures long-table construction,
not branch-file reading, plotting or the complete summary stage.

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_presence_absence.py \
  --output presence-absence.json
```

Use `--support-root` with a complete baseline support tree and keep the same
`GG_CONTAINER_DOCKER_IMAGE` immutable image ID for both runs.

## Alignment-derived MAPNH frequencies

For sequences of at least 128 codons, IQ-TREE-to-MAPNH conversion counts repeated
codons before accumulating the three nucleotide-position frequencies. Shorter
sequences retain the original loop. RNA/case normalization, whole-codon rejection
for ambiguity/gaps, incomplete terminal codons and leaf matching remain unchanged;
this dispatch affects runtime only.

Linux arm64 Docker, Python 3.12.14, base `46b00fc`, same immutable image, one
warmup and three measured fresh processes: reading and computing F3X4 frequencies
for a full 1,024-gene alignment and its two half-size subroots, with 512 codons
per gene, took 0.469 to 0.169 s (2.77×). The fixture has biased nucleotide
frequencies, all 64 valid codons, ambiguity/gaps, Unicode invalid bases,
RNA/lowercase sequences, incomplete tails and normalized leaf labels. Frequencies,
derived theta values and source bytes match. A 32-gene/32-codon control remained
about 1.11 ms on both versions. Process peak RSS stayed about 143 MiB (137.7 MiB
for the control), with no memory saving established. Timing includes three FASTA
reads, matching and frequency calculation, excluding fixture construction,
theta conversion and fingerprints; it does not measure IQ-TREE inference or
the full MAPNH conversion.

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_alignment_frequencies.py \
  --output alignment-frequencies.json
```

Use `--support-root` with a complete baseline support tree and keep the same
`GG_CONTAINER_DOCKER_IMAGE` immutable image ID for both runs. Add
`--genes 32 --codons 32` for the short-sequence control.

## CSUBST candidate identity and input-state annotation

Candidate input-state annotation reads compact tuples of its three needed
identifier columns, retaining the existing row behavior for non-string IDs or
ambiguous columns. Input signatures, missing-input text, analysis-key hashes and
cache names remain exact; verification and cache-completion rules are unchanged.

Linux arm64 Docker, Python 3.12.14, pandas 3.0.6, base `545bcf0`, same immutable
image, one warmup and three measured fresh processes: 8,192 candidate rows across
256 families, 32 additional probability columns, 1,184 candidates with missing
inputs and nondefault source-row indexes:

| Input-state annotation | Before median | After median | Ratio |
| --- | ---: | ---: | ---: |
| 8,192 candidates | 0.101 s | 0.0160 s | 6.32× |
| 32 candidates / 32 families | 0.936 ms | 0.928 ms | No material gain |

Candidate-ID assignment separately avoids per-row Series construction and named
cell lookups when columns are unique/flat and orthogroup IDs are strings. It uses
the same mixed array values as the original loop, preserving scalar precision,
nullable values, identity JSON, analysis/cache hashes and duplicate rejection.
Other inputs retain the original row coercion and error order. Compared with
`cb9c3f2` using the same environment/process sampling and full identity fields:

| Candidate-ID assignment | Before median | After median | Ratio |
| --- | ---: | ---: | ---: |
| 8,192 candidates / 256 families | 0.278 s | 0.143 s | 1.94× |
| 32 candidates / 32 families | 1.79 ms | 1.24 ms | 1.44× |

Every output value, type, column and index/order matches; the input table remains
unchanged. Large-case peak RSS stayed about 139 MiB, without an established memory
saving (about 143–144 MiB with ID assignment). Timing includes copying and
annotation but excludes fixture construction and fingerprinting; peak RSS includes
both. This does not measure input discovery,
artifact hashing, scientific scans or site analysis; the first table excludes ID
assignment.

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_candidate_input_state.py \
  --output candidate-input-state.json
```

Use `--support-root` with a complete baseline support tree and keep the same
`GG_CONTAINER_DOCKER_IMAGE` immutable image ID for both runs. Add
`--candidates 32 --families 32` for the small control.
Add `--assign-ids` to measure ID assignment separately and fingerprint its complete
output through subsequent input-state annotation.

## GFF transcript statistics

Longest-transcript selection collects row positions and builds one result table,
preserving gene discovery order, annotation row order, tie handling and extra
column types. CDS phase validation indexes coordinates once per gene, retaining
every phase record at duplicate coordinates. Coordinate ordering, UTR validation,
strict/report policies, missing phases and trans-splicing behavior are unchanged.

Linux arm64 Docker `local/genegalleon:dev`, Python 3.12.14, pandas 3.0.6, base
`f6f68bf`, one warmup and three measured fresh processes: the live in-memory
selection → structure validation → statistics path for 2,000 genes took 2.607
to 0.686 s (3.80×). Selection took 0.495 to 0.0591 s; structure validation took
2.072 to 0.591 s. Median peak RSS was 129.5 to 98.8 MiB (23.7% less). The
fixture includes both strands, alternative isoforms, explicit UTRs and duplicate
CDS records. Intermediate/final tables, column types/order and warnings match.
Time excludes fixture construction and result fingerprinting; process peak RSS
includes both. This does not measure GFF download, file parsing or sequence
resolution.

A further pass (base `04dc6dd`, same workload/runtime/repetitions) iterates
coordinate columns directly rather than constructing a small coordinate frame
and tuple iterator for each gene. Missing/duplicate-column behavior and coordinate
casts remain unchanged. Complete selection → structure → summary time fell from
0.686 to 0.433 s (1.58×); structure validation alone fell from 0.592 to 0.336 s.
Selected/annotated/final table fingerprints and warnings match. Median process
peak RSS was 98.9 to 99.3 MiB, with no material memory saving established.

Gene-local phase validation also avoids casts for exact Python string values
and columns already using the native integer dtype. Other types retain the
original pandas conversion, and each gene is validated before advancing to the
next. A further comparison (base `b3272f5`, same image/runtime/repetitions) gives:

| Complete in-memory GFF workload | Before median | After median | Ratio |
| --- | ---: | ---: | ---: |
| 2,000 genes, 3 exons each | 0.425 s | 0.321 s | 1.32× |
| 64 genes, 128 exons each | 0.118 s | 0.121 s | 0.97× |

Structure-validation medians were 0.327 to 0.226 s for the first case and about
0.0376 s on both versions for the second; no speedup is established for the
many-exon control. All intermediate/final tables, types/order and warnings match.
Process peak RSS remained about 99.1 / 96.9 MiB for the two workloads.
Add `--genes 64 --exons 128` to the command below to reproduce the control;
the existing three-exon default is unchanged.

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_gff_summary.py \
  --output /tmp/gff-summary.json
```

Use `--support-root` with a complete baseline support tree for comparison.

## CDS coordinate-compatibility reporting

Sequence/coordinate compatibility reporting updates scalar cells directly,
retaining collection/MultiIndex selection behavior. Exact length and terminal-N
checks, mismatch reasons, cleared coordinates/structure fields, value types,
gene order and partial updates before a later error remain unchanged.
With unique columns and string gene IDs, iteration reads compact tuples of the
two needed columns; legacy inputs retain their original row/error behavior.

Linux arm64 Docker, Python 3.12.14, pandas 3.0.6, base `efebae1`, same immutable
image, one warmup and three measured fresh processes:

| Complete compatibility reporting | Before median | After median | Ratio |
| --- | ---: | ---: | ---: |
| 4,096 genes, 2,925 compatible / 1,171 mismatched | 1.056 s | 0.233 s | 4.53× |
| 4,096 genes, all compatible | 0.246 s | 0.0758 s | 3.25× |
| 32 genes, mixed lengths | 10.4 ms | 3.27 ms | 3.18× |

The subsequent tuple iteration change, compared separately with the scalar-cell
implementation in `d886351`, used the same image and process sampling:

| Further row-iteration improvement | Before median | After median | Ratio |
| --- | ---: | ---: | ---: |
| 4,096 genes, 2,925 compatible / 1,171 mismatched | 0.244 s | 0.165 s | 1.48× |
| 4,096 genes, all compatible | 0.0743 s | 0.0224 s | 3.32× |
| 32 genes, mixed lengths | 3.19 ms | 2.61 ms | 1.22× |

Complete tables, columns/types, indexes/order and input sequence records match.
Process peak RSS stayed about 90 MiB (86 MiB for the small control), with no
material memory saving established. Timing excludes fixture construction and
fingerprinting; peak RSS includes both. This does not measure GFF parsing,
transcript selection, sequence resolution or the entire annotation stage.

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_structure_compatibility.py \
  --output structure-compatibility.json
```

Use `--support-root` with a complete baseline support tree and keep the same
`GG_CONTAINER_DOCKER_IMAGE` immutable image ID for both runs. Add
`--scenario compatible` for the all-compatible case or `--genes 32` for the
small mixed-length control.

## CDS validation and source selection

CDS admission checks ambiguous bases once over the exact existing internal
codon span. It still excludes the first incomplete codon and terminal codon;
empty stop/dual-coding codon sets avoid unnecessary membership scans. Genetic
codes, padding, phase constraints, stop counts, error precedence and selection
policy are unchanged. Source and selected-sequence hashes remain identical.

Linux arm64 Docker, Python 3.12.14, pandas 3.0.6, base `a5bb48b`, same immutable
image, one warmup and three measured fresh processes: 1,024 genes with 500 body
codons (normally 1,506 bases), including whitespace/lowercase, internal stops,
ambiguity, missing supplied sequences, extensions, partial CDS and disagreeing
sources:

| CDS workload | Before median | After median | Ratio |
| --- | ---: | ---: | ---: |
| Unconstrained admission, including CDSKit padding | 0.281 s | 0.214 s | 1.31× |
| Admission with explicit GFF phases | 0.127 s | 0.0590 s | 2.15× |
| Complete supplied/genomic two-source selection | 0.420 s | 0.301 s | 1.40× |

Every acceptance/rejection, selected source, reason, padded sequence and record
hash matches. Median process peak RSS stayed about 99 MiB. Timing excludes
fixture construction and JSON fingerprints; peak RSS includes both. These
measurements exclude GFF/genome parsing and artifact verification/publication.

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_cds_evaluation.py \
  --output cds-evaluation.json
```

Use `--support-root` with a complete baseline support tree and keep the same
`GG_CONTAINER_DOCKER_IMAGE` immutable image ID for both runs.

## Intron-ASR result validation

Native node-name/class/parent validation iterates compact tuples of the needed
columns instead of creating a pandas Series per node. Input row order and parent
index lookup remain unchanged, including inferred CSV indexes and rejection of
ambiguous duplicate indexes. Native-ID to clade-rank translation, probability,
observed-count and imputation checks, output columns/types and error precedence
are retained.

Linux arm64 Docker, Python 3.12.14, pandas 3.0.6, base `3a4be10`, same immutable
image, one warmup and three measured fresh processes: a 2,048-tip balanced tree
and 4,095 native rows in reverse order, with missing leaf observations, internal
probabilities, literal identifier values and additional columns. Complete table
read, tree parse, validation and ID translation took 0.0821 to 0.0293 s (2.80×).
Every output value, column type, row/index order and source byte matches. Median
process peak RSS stayed about 151.5 MiB. Timing excludes fixture construction,
imports and fingerprints; this does not measure ancestral-state inference or
the full branch-statistics stage.

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_asr_intron_loading.py \
  --output asr-intron-loading.json
```

Use `--support-root` with a complete baseline support tree and keep the same
`GG_CONTAINER_DOCKER_IMAGE` immutable image ID for both runs.

## Wide-header validation

Database and scan-schema preflight count each column name once, preserving
sorted duplicate diagnostics and rejection before DB replacement. Duplicate
checking no longer scans the full header for every column. The header reader
now follows standard TSV quoting and UTF-8 BOM handling, matching the body
reader. Quoted tabs/newlines, quotes and literal spaces in names are preserved;
duplicates remain visible and fatal instead of being renamed by pandas.

Linux arm64 Docker `local/genegalleon:dev`, Python 3.12.14, base `f6f68bf`, one
warmup and three measured fresh processes: validating 128 scan headers with
2,048 columns took 3.780 to 0.0466 s (81.1×) for raw files and 3.862 to 0.1304 s
(29.6×) for ZIP-held files. Peak RSS stayed about 105 MiB raw / 107 MiB ZIP.
The fixture includes one deliberate duplicate-column error; the complete sorted
diagnostic matches. These measure strict schema preflight, not full DB creation.

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_schema_validation.py \
  --output /tmp/schema-validation.json
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_schema_validation.py \
  --storage zip --output /tmp/schema-validation-zip.json
```

Use `--support-root` with a complete baseline support tree for comparison.

## Declaration, trait and species-tree validation

Provenance declaration keys, trait headers/selections and copy-number species
labels are counted once instead of rescanning each list for every value. Trait
selection also builds one available-name set. Sorted duplicate diagnostics,
selection order, error precedence and rejection before artifact reads remain
unchanged; file hashing and content/completion checks are retained.

Linux arm64 Docker, Python 3.12.14, pandas 3.0.6, base `3013667`, same immutable
image, one warmup and three measured fresh processes:

| Synthetic validation workload | Before median | After median | Ratio |
| --- | ---: | ---: | ---: |
| 4,096 input declarations, full contract plus cross-kind duplicate failure | 0.578 s | 0.132 s | 4.38× |
| Four-row / 4,096-column trait TSV, read and valid/duplicate selections | 0.404 s | 0.0292 s | 13.8× |
| 4,096-leaf species tree, valid and duplicate-label validation | 0.236 s | 0.0168 s | 14.1× |

Complete contract, trait TSV, selections, Newick/leaf order and diagnostics
match. Combined median process peak RSS was about 160.6 MiB before/after.
The declarations reference one small shared source, so this measures declaration
scaling rather than bulk content hashing. Fixture construction and fingerprints
are excluded from timing and included in peak RSS. These large validation
fixtures do not establish full analysis speedups or gains for small inputs.

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_declaration_validation.py \
  --output declaration-validation.json
```

Use `--support-root` with a complete baseline support tree and keep the same
`GG_CONTAINER_DOCKER_IMAGE` immutable image ID for both runs.

## Expression replicate preparation

Raw expression preparation compiles observation metadata and source-column
positions once, then reads tuple rows. Missingness checks also avoid creating
one pandas Series per gene. Gene/observation order, literal column names,
pairing, technical IDs, batches, skipped-response reasons and numeric validation
remain unchanged. Known-SE and unreplicated output formats remain compatible.

Linux arm64 Docker, Python 3.12.14, pandas 3.0.6, base `49a966b`, same immutable
image, one warmup and three measured fresh processes: 3,000 genes, eight retained
responses and three replicates, plus constant/missing-leaf control responses:

| Preparation workload | Before median | After median | Ratio | Median peak RSS before/after |
| --- | ---: | ---: | ---: | ---: |
| Automatic independent replicates, 64,166 output rows | 0.429 s | 0.0904 s | 4.75× | 99.5 / 99.2 MiB |
| Explicit paired replicates with technical IDs/batches, 8,823 rows | 0.381 s | 0.0655 s | 5.82× | 81.5 / 81.4 MiB |

Complete TSV contents, column types/order and metadata match. Fixtures include
awkward literal headers, reversed sample metadata, partial missingness and
entirely missing observations. Timing includes input reads and expression
validation/conversion; it excludes fixture construction, result fingerprinting,
predictor preparation, output publication and model fitting. Process peak RSS
includes construction and fingerprinting. These are preparation measurements,
not whole-analysis speedups.

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_expression_formatting.py \
  --case unpaired --output expression-unpaired.json
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_expression_formatting.py \
  --case paired --output expression-paired.json
```

Use `--support-root` with a complete baseline support tree and keep the same
`GG_CONTAINER_DOCKER_IMAGE` immutable image ID for both runs.

## Query branch markers

Direct query matching indexes normalized exact IDs and the gene-list priority
once. Each tip checks suffixes at the existing underscore, hyphen and dot
boundaries. Exact sources retain priority, gene-list sources retain catalog
order, and FASTA IDs still require exact matching. Branch annotation reads the
node/tip columns directly rather than constructing Series for every branch.
Best-hit groups retain the existing global query order without sorting again
per node. BLAST filtering, best-hit selection and marker columns are unchanged.

Linux arm64 Docker, Python 3.12.14, pandas 3.0.6, base `5f8cf15`, same immutable
image, one warmup and three measured fresh processes: 4,096 tips plus four edge
controls, 512 base query IDs plus overlapping suffixes and exact FASTA IDs,
internal branches and competing BLAST hits. Direct source matching took 1.704
to 0.0283 s (60.3×); complete 4,613-row branch annotation, including input/output
and BLAST processing, took 1.818 to 0.0506 s (35.9×). Complete direct-source
maps, output TSV bytes and reported marker count match. Median peak RSS was
about 78.8 / 78.6 MiB. Fixture construction/fingerprinting are excluded from
timing and included in peak RSS. Gains for smaller query catalogs will differ;
this does not measure the entire gene-evolution stage or plotting.

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_query_markers.py \
  --output query-markers.json
```

Use `--support-root` with a complete baseline support tree and keep the same
`GG_CONTAINER_DOCKER_IMAGE` immutable image ID for both runs.

## Seeded single-copy ortholog decay

Repeated decay calculations use cumulative Boolean intersections instead of
integer counters for all-present/all-single-copy metrics. When permutations
reuse more columns than a complete table pass, species-major presence and
single-copy masks are computed once. Small partial runs avoid that full-table
work. Species permutations, seed handling, requested count order, all/selected
metrics and summary/plot formats remain unchanged; input arrays are untouched.
Input tables require finite, non-negative integer counts within int64 range;
fractional values are rejected before conversion rather than silently truncated.

Linux arm64 Docker, Python 3.12.14, pandas 3.0.6, base `dff3b95`, same immutable
image, one warmup and three measured fresh processes: 8,192 orthogroups, 64
species, 128 permutations, five requested subset sizes; selected mode also
includes 4,096 selected orthogroups. Calculation plus summary took 0.133 to
0.0216 s (6.13×) in all mode and 0.163 to 0.0265 s (6.17×) in selected mode.
Every seeded replicate value and complete summary TSV/type/order matches.
Median process peak RSS was 78.6 to 79.9 MiB (about 1.3 MiB more). Cached masks
use two Boolean bytes per all-table cell plus one per selected-table cell;
cumulative state uses less space than the previous integer counters.
Fixture construction/fingerprinting are excluded from timing and included in
peak RSS. These measurements exclude file parsing and figure rendering/export.

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_ortholog_decay.py \
  --output ortholog-decay.json
```

Use `--support-root` with a complete baseline support tree and keep the same
`GG_CONTAINER_DOCKER_IMAGE` immutable image ID for both runs.

## Unrooted branch-support mapping

Canonical split keys use the same smaller-side and lexicographic tie rules.
A cardinality bound avoids constructing a provably larger complement, and only
the selected side is sorted when sizes differ. Leaf uniqueness, topology/leaf
set agreement, complete support coverage, root-edge consistency and support
range checks are retained. Missing leaf names are rejected rather than converted
into the literal string `"None"`; a real leaf named `"None"` remains valid.

Linux arm64 Docker, Python 3.12.14, base `989e44b`, same immutable image, one
warmup and three measured fresh processes. Two trees have identical topology
with reversed child order, explicit support and deterministic branch IDs:

| Complete mapping workload | Before median | After median | Ratio |
| --- | ---: | ---: | ---: |
| Balanced tree, 2,048 tips | 0.931 s | 0.0145 s | 64.3× |
| Comb tree, 1,024 tips | 0.208 s | 0.0681 s | 3.05× |

Every mapped branch value and diagnostic matches, and input trees are unchanged.
Median process peak RSS stayed about 148.5 MiB (balanced) / 167.7 MiB (comb).

A further pass (base `70b0185`, same workloads/runtime/repetitions) uses integer
bit masks over one sorted, operation-local tip catalog. It avoids descendant
name sets and sorted tuple keys during normal mapping. Equal-size complements
still select the lexicographically smaller side; errors decode the original
tuple previews, and the existing tuple-key helper APIs remain available.
All mapped values, diagnostics and input trees again match:

| Complete mapping workload | Before median | After median | Ratio | Process peak RSS |
| --- | ---: | ---: | ---: | ---: |
| Balanced tree, 2,048 tips | 0.0145 s | 0.0103 s | 1.41× | 148.4 → 147.8 MiB |
| Comb tree, 1,024 tips | 0.0668 s | 0.00497 s | 13.4× | 167.7 → 145.3 MiB |

The comb workload uses 13.4% less total process peak memory. Masks still have
quadratic worst-case bit storage in tree size; this is a compact representation,
not a linear-memory guarantee. Timing includes both trees' validation, split construction, support
checks and final mapping; it excludes fixture construction/import, fingerprints,
tree inference, support estimation and the rest of branch-statistics generation.

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_support_mapping.py \
  --output support-mapping-balanced.json
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_support_mapping.py \
  --shape comb --tips 1024 --output support-mapping-comb.json
```

Use `--support-root` with a complete baseline support tree and keep the same
`GG_CONTAINER_DOCKER_IMAGE` immutable image ID for both runs.

## HGT candidate branch summaries

Candidate summaries index the first retained leaf taxon once per branch instead
of filtering the leaf table for every candidate gene. Candidate order and
duplicates, missing/absent taxa, first-leaf selection, evidence counts, lineage
resolution and representative annotations remain unchanged.
Gene aggregation computes group counts once, converts branch IDs once, uses
tuple first rows and reuses taxonomy column names. Existing sort/tie policies,
ID formatting, categorical group order and optional-column defaults are retained.

Linux arm64 Docker, Python 3.12.14, pandas 3.0.6, base `26a21cc`, same immutable
image, one warmup and three measured fresh processes: two overlapping candidate
branches over 2,048 genes, reversed candidate order, an absent gene, mixed
expression/intron/synteny evidence and contamination annotations. Branch and
raw-gene summarization took 0.778 to 0.143 s (5.44×). Complete raw/final branch,
gene and orthogroup tables, types and order match; median peak RSS stayed about
204.5 MiB. Further gene/orthogroup aggregation optimization (base `79dd19f`,
same workload/runtime/repetitions) took 0.912 to 0.137 s (6.65×), with complete
tables/types/order again matching. Corresponding median process peak RSS was
204.4 to 191.8 MiB (6.2% less). Timing excludes fixture construction, final fingerprints, SQLite reads,
scaffold context and plotting. The local taxonomy resolver has no database, so
this does not measure NCBI lookup cost or the entire HGT stage.

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_hgt_summary.py \
  --output hgt-summary.json
```

Use `--support-root` with a complete baseline support tree and keep the same
`GG_CONTAINER_DOCKER_IMAGE` immutable image ID for both runs.

## HGT host-scaffold classification

Isoform disagreement is determined with one grouped distinct-label count rather
than constructing a small table for each locus/rank. The same scaffold,
count-unit, locus and rank boundaries are retained, including agreement at one
rank and disagreement at another. Input filtering, taxonomy validation, locus
counting, missing-group behavior and scaffold composition are unchanged.

Linux arm64 Docker, Python 3.12.14, pandas 3.0.6, base `5e07d5d`, same immutable
image, one warmup and three measured fresh processes: 2,048 input genes over
32 scaffolds and seven ranks, including shared isoforms, conflicts, missing
taxonomy, absent coordinates and trans-splicing. Complete classification and
aggregation took 0.441 to 0.0758 s (5.82×). All 14,161 gene/rank rows and 224
scaffold/rank rows, types, columns, order and index match; inputs remain unchanged.
Median process peak RSS was 84.3 to 82.7 MiB. The local rank resolver returns
fixed lineages; this excludes NCBI database lookup, file parsing, fixture
construction and fingerprints, and is not a full HGT-stage measurement.

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_scaffold_construction.py \
  --output scaffold-construction.json
```

Use `--support-root` with a complete baseline support tree and keep the same
`GG_CONTAINER_DOCKER_IMAGE` immutable image ID for both runs.

## HGT host-scaffold context

Gene context is assigned by column rather than by individual cell. Shared
scaffold totals are referenced without copying each metric dictionary per gene.
Per-species input streaming, strict taxonomy-table validation, rank missingness,
union-of-candidate-loci background exclusion and recipient-only branch pooling
are unchanged. Output column types, row/index order and scalar duplicate-index
update behavior are retained.

Linux arm64 Docker, Python 3.12.14, pandas 3.0.6, base `1ce15e2`, same immutable
image, one warmup and three measured fresh processes: two species with 2,048
taxonomy genes each, seven ranks, 64 total scaffolds, 2,050 candidate gene rows
and 33 candidate branches. Complete attachment, including file reads, strict
validation, full/background composition and branch pooling, took 1.520 to
0.782 s (1.94×). Complete branch/gene TSV contents, types and index/order match.
Median process peak RSS was 102.7 to 103.0 MiB (about 0.3 MiB more). The fixture
includes shared isoform loci, reused IDs across species, missing mappings and
unresolved recipients. Timing excludes fixture construction and fingerprinting;
peak RSS includes both. No taxonomy database lookup, candidate discovery or
full HGT-stage speedup is measured.

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_scaffold_context.py \
  --output scaffold-context.json
```

Use `--support-root` with a complete baseline support tree and keep the same
`GG_CONTAINER_DOCKER_IMAGE` immutable image ID for both runs.

## Alignment-statistics I/O

Both summary readers keep at most eight pending reads per configured worker,
consume completion notifications and discard completed futures. They retain
their original worker count, named-column mapping and failure behavior. Inventory
and output tables still require memory proportional to the number of families;
only queued/completed task retention is bounded. Queued work is cancelled on
failure or iterator close, and active I/O finishes before the executor exits.

Linux arm64 Docker `local/genegalleon:dev`, Python 3.12.14, base `e56c314`, four
workers, one warmup and three measured fresh processes: 10,000 raw statistics
TSVs took 7.01 s to 6.30 s (1.11×), with median peak RSS 110.2 to 98.7 MiB
(10.5% less). A 1,000-family query2family ZIP fixture took 1.05 to 1.04 s with
about 99 MiB RSS on both versions; no speedup is established for that smaller
workload. Complete result TSV bytes match for both fixtures. Fixture generation
is excluded from time, but included in process peak RSS.

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_alignment_summary.py \
  --output /tmp/alignment-summary.json
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_alignment_summary.py \
  --reader query2family --storage zip --families 1000 --output /tmp/alignment-summary-zip.json
```

Use `--support-root` with a complete baseline support tree for comparison.

## Verification

`workflow_api.py capabilities` advertises
`verify_batches=shared-source-verification-v1` and `verify_batch_limit=32`.
Use `verify --batch-file PLAN.json` instead of `--family-id`, keeping the existing
`--root`, `--workspace-root`, required steps and terminal profile. A plan is:

```json
{"schema":"genegalleon-verify-batch-v1","requests":[
  {"family_id":"OG0000001","attempt":"/workspace/output/observations/ATTEMPT",
   "recorded_workspace_root":"/workspace"}
]}
```

The `verify` envelope contains `batch_schema` and ordered `results`, each with
the unchanged single-family result schema. At most 32 distinct families are
allowed. Missing or stale artifacts stay family-local unverified results;
malformed evidence or a failed shared final content fence rejects the query.
Each exact attempt, receipt, family status and workspace mapping is rechecked.
Hashes of shared inputs are read once and rehashed before returning any results.
kfauto negotiates this capability, consumes each exact result in the current
32-record collector chunk and retains scalar queries for older runtimes.

## Audited database publication

The gene-summary core uses `gene_family_database_pipeline.py` to keep audit,
private DB build and provenance construction in one process. The pipeline accepts
repeated `--audit=TOKEN`, `--database=TOKEN`, and `--record=TOKEN` options carrying
the same argv tokens as the three standalone commands. It requires one matching
workspace/store and sole database output. Standalone commands remain available.

The audit checks all manifests, required inventories and CSUBST branch identity.
The DB remains private while its record is constructed from the private bytes at
the final declared output path. Unique source contents, collection membership,
manifests and archive generation are then revalidated before publication. The
existing namespace/manifest locks, atomic create/replace modes and exact attempt
record are retained. Audit progress shows `database_build` during SQL construction
and `source_revalidation` during the final fence; audit telemetry is never a
workflow completion proof.

Uniform large TSVs spool their original inferred pandas frames to private disk
and read them back without reparsing CSV. RAM stays chunk-bounded. Any dtype
promotion retains the original second CSV pass, preserving leading zeros,
nullable Boolean values, missingness, quoted tabs and multiline fields. Spools
are private temporary files and are removed on success, error or iterator close.
Spooling stops and cached chunks are discarded once a dtype change makes reuse
impossible. The extra disk writes are the main tradeoff; mixed-type files cannot
claim this parse saving.

## PDF worker

`Rscript workflow/support/tree_plot_batch.r PLAN.json` renders 1..32 jobs in one
R process, reusing loaded packages. The plan is a JSON array:

```json
[{"id":"OG0000001","cwd":"/workspace/scratch","output":"/workspace/result.pdf",
  "args":["--stat_branch=/workspace/stat.branch.tsv",
          "--max_delta_intron_present=-0.5","--panel_widths_mm=tree:60",
          "--panel1=tree,bl_rooted,no,no,L","--show_branch_id=no",
          "--event_method=species_overlap","--species_color_table=PLACEHOLDER",
          "--pie_chart_value_transformation=identity","--long_branch_display=no"]}]
```

Supply absolute input file paths. Each job uses a fresh environment, input cache
and scratch directory. Graphics devices, options and the optional species parser
are reset between jobs; garbage collection bounds retained family data. Input
file content/signatures (including comma-containing filenames) and the exact
cached renderer source bytes are fenced
before atomic PDF publication. An individual render failure preserves its prior
output and does not contaminate later jobs.
IDs and output paths must be unique. `PLAN.json.results.json` reports ordered
per-job exit codes, including 42 for unavailable optional ggimage, with
`completion_evidence=false`. The process exits nonzero if any job fails.

A caller must still record and verify each family's provenance and exact attempt
using the workflow's normal stage functions. This renderer does not submit jobs
or change single-family array scheduling; one-family jobs retain the existing R
CLI, and their startup is not amortized by this batch interface.

## Comparable synthetic measurements

Linux arm64 Docker `local/genegalleon:dev`, Python 3.12.14; same container and
inputs, fresh processes, one warmup plus three trials for API/PDF/pipeline.
Large-file DB measurements use three trials, after prior workload warmup.
The comparison base is `3aecee4`. Times exclude fixture generation. No scheduler
workload or production completion-rate claim is involved.

| Workload | Before median | After median | Ratio |
| --- | ---: | ---: | ---: |
| 16 families, five declarations each, shared 64 MiB source | 0.971 s | 0.125 s | 7.77× |
| Eight tree PDFs, separate R processes versus one worker | 9.188 s | 2.796 s | 3.29× |
| 128 audited families, 32 MiB alignments, complete DB/record pipeline | 1.164 s | 1.019 s | 1.14× |
| One 131,072-row TSV, 24 numeric metrics, complete DB build | 1.413 s | 1.390 s | 1.02× |

PDF and DB timings were rerun after adding the renderer source fence and early
spool disposal. The many-small-file DB control changed 2.128 to 2.236 s; an earlier
run changed 2.153 to 2.129 s. Large-file DB gains also varied from 1.10× to 1.02×,
so removing the second CSV parse does not establish a reliable overall DB gain.
API parent peak RSS was about 108 MiB before/after. PDF worker peak RSS increased
about 11 MiB (229 to 240 MiB); a warmed R namespace remains resident. Uniform
large-file DB peak RSS stayed about 134 MiB. The combined audited DB
child retained about 13 MiB more (118 to 131 MiB) because both audit and DB
modules remain loaded; the parent fixture process used less memory.
Verification decisions/contracts, sorted SQLite schema/rows, recorded input
fingerprints and audit rows match. PDF drawing bytes match after removing only
CreationDate/ModDate. SQLite binary hashes can differ because parallel insertion
order varies; each generated provenance record must match its own DB bytes.

Run comparable revisions through the existing runtime wrapper, without parallel
benchmarks or tests:

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_remaining_io.py \
  --case verify-batch --output /tmp/verify.json
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_remaining_io.py \
  --case pdf-batch --output /tmp/pdf.json
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_summary_pipeline.py \
  --output /tmp/pipeline.json
bash workflow/tests/run_in_runtime.sh python workflow/tests/benchmark_orthogroup_database.py \
  --output /tmp/database.json
```

The first three accept `--support-root` for a complete saved baseline support
tree; the DB benchmark accepts `--source` for the baseline generator. Compare
logical outputs and recorded inputs, validate each record's actual output hash,
and report RSS/disk tradeoffs alongside timing.
