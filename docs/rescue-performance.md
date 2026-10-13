# Rescue performance

The optimization keeps each synteny interval independent, caps alignment threads
at the worker CPU budget and preserves candidate order. Whole-genome query
deduplication uses exact protein sequence equality; results are expanded back
to every original candidate before validation. The comparison cache binds
actual BED/protein inputs, comparison options and the relevant toolchain.
See [gene-model rescue](gene-model-rescue.md) for controls and cache guarantees.

## BUSCO resource bounds

The refinement comparison admits at most `--jobs` species at once. The first
observed failure stops further admission and waits for admitted species to
finish naturally. Successful results retain the original species order.

The BUSCO compatibility helper sets the child's `MMSEQS_NUM_THREADS` to the
requested `--cpu`/`-c` value, or one when it is absent, while retaining a positive
lower inherited limit. MetaEuk subcommands that omit `--threads` can otherwise
use the host's online CPU count. The helper preserves the parent environment,
BUSCO arguments and exit status. These changes do not alter model selection or
BUSCO scoring thresholds; no whole-workflow speedup is asserted.

The ordinary comparison contract includes the evaluator and wrapper hashes.
Upgrading these files requires a new report directory. Changing only `--jobs`
keeps the existing contract; completed leaves still require their exact source,
summary and full-table hashes before reuse.

## Baseline evidence

The 2026-10-04/05 run used GeneGalleon 0.8.124, x86_64 SIF, miniprot
0.18-r281, 8 CPUs and 64 GB per model worker, with four species concurrently.
The 154 synteny jobs used four CPUs each, eight jobs concurrently, and completed
in 46 minutes. In the model stage:

| Species | Independent intervals | Interval elapsed time | Fallback queries | Distinct protein sequences |
|---|---:|---:|---:|---:|
| Ancistrocladus abbreviatus | 41,149 | about 64 min | 332,914 | 148,189 |
| Triphyophyllum peltatum | 42,920 | about 61 min | 369,423 | 140,630 |
| Vitis vinifera | 28,579 | about 41 min | 218,036 | 107,665 |

Interval elapsed times above are spans between the first and last miniprot
log writes, including inter-process and file overhead. Twenty evenly spaced
interval logs per species had median miniprot wall times of 0.080–0.090 s and
CPU/wall ratios near one, despite eight requested threads. Ancistrocladus
fallback used 1,842 s; its BUSCO used 813 s. These are diagnostic observations,
not a before/after speedup measurement.

The highest observed miniprot index RSS among 21 completed species was 24.836 GB
in Dionaea. That excludes the Python parent and BUSCO, so it does not justify a
smaller whole-worker memory request. Six workers at the existing 64 GB request
are a site-specific concurrency choice, not a portable default.

## Measured search improvements

On the same x86_64 SIF and eight-CPU budget, one warmup and three measured
repetitions on Ancistrocladus evidence produced the following median wall times:

| Operation | Legacy | Optimized | Observed ratio |
|---|---:|---:|---:|
| 256 evenly spaced intervals | 19.372 s | 3.395 s | 5.71× |
| First 4,096 whole-genome queries (2,626 unique sequences) | 23.594 s | 16.895 s | 1.40× |

Every restored GFF byte matched: combined interval SHA256
`cd67363141a2cd8343d9066399f3600ff4a673d49968514fef88582ee3237008` and genome
SHA256 `ccd6ad3cbc1b53bc1e939f943a0f245a8e28824a97cf64c50ac03f827cd7991a`.
The SIF SHA256 was
`2685031e1bea7b425a27fdfa5d03f611c9eaedffe7d60a34a79c2d9a0bb82389`;
miniprot was 0.18-r281. Search samples were interleaved after the warmup and the
shared genome index was built once. Query deduplication/expansion are included
in optimized mapping times. Genome mapping samples varied (legacy 23.52–25.19 s,
optimized 16.20–20.57 s), so the ratio describes this bounded workload and is
not a predicted full-workflow speedup. The Python parent's maximum RSS was
about 132 MiB, excluding miniprot and the shared index build.

The full Ancistrocladus equivalence check also passed: all 41,149 interval GFFs
and the expanded whole-genome evidence for 332,914 original queries (148,189
unique proteins) matched the original 0.8.124 producer receipts. The combined
interval SHA256 was
`3757ba950e54dbf5b3e4cf2544459956e98f3803f88b7b19cd6b378f3579dbd1`;
the full genome GFF SHA256 was
`d7867cb9640187dedc132cd2ddf76c7fe0faeab9f3ce0fe09a15508f5e417136`.
Frozen source files and receipts were rechecked unchanged afterward. This was
one optimized run without a legacy rerun or warmup, so its timings are diagnostic
and do not establish an additional controlled speedup ratio.

## Reproducible bounded benchmark

Use the repository runtime wrapper with completed real rescue evidence. It
checks sampled files against the producer's hashed receipt and rechecks inputs
afterward, including the fallback query FASTA. The receipt must belong to the
selected plan/species. Recorded interval numbers must be contiguous and every
recorded interval must still contain its region, queries and original GFF;
missing directories cannot silently reduce the coverage of a full check.
Every trial compares SHA256 of all restored GFF bytes, including
query names, order, model IDs and unmapped evidence. A mismatch stops the run.
The default covers 256 evenly spaced real intervals and the first 4,096 fallback
queries, one warmup plus three measured repetitions for each method.

```bash
GENEGALLEON_SIF_EXTRA_BINDS=/data/rescue_run \
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_rescue_search.py \
  --evidence /data/rescue_run/gene_model_rescue \
  --species Ancistrocladus_abbreviatus --cpus 8 \
  --output workspace/output/rescue_benchmark
```

The output path must be new. It retains commands, logs, GFF files, the copied
genome/index and `result.json`, including source/runtime identities, exact sample
hashes, timings and medians. Fixture extraction and one shared index build are
outside the search timings. Deduplication timings include writing unique queries
and expanding GFF evidence. RSS fields distinguish the Python parent from the
alignment subprocess; they are not total concurrent-worker memory estimates.
The bounded results do not establish full-workflow speed or the sensitivity of
new gene models. Resource overrides and cache reuse/invalidation are also covered
by the existing array and real-tool rescue tests.

`--check-existing` performs one optimized run over every original interval and
fallback query, comparing each interval GFF and the complete expanded genome
GFF with the original receipt hashes. It skips the legacy run and warmups. Use
this for a full real-data equivalence check; its single timings are diagnostic,
not a controlled speedup benchmark.

When the producer did not run an interval or fallback search, that phase is
listed in `result.json` under `skipped`, with no timings or index build.
Partially recorded or missing evidence fails the check. Output equivalence and
input-change checks remain active under Python's `-O` option.

Verified prediction reuse also avoids generating unused local-search inputs.
Only uncached windows are fetched and written; a completely cached local search
does not create `regions.fa` or `queries.fa`. Genome indexing and its mandatory
input QC still run. A predictor-facing genome alias is created only when an
actual genome-wide or GeMoMa search needs it.

A bounded benchmark with 2,048 candidates sharing 64 windows used one warmup
and three measured trials. With 99% cached candidates, local FASTA output fell
from 67.5 MB to 0.659 MB, and median writing plus eight-thread hash verification
fell from 0.754 s to 0.011 s. Complete reuse wrote no local-search FASTA bytes;
uncached inputs retained identical bytes. These measurements cover input
preparation, not alignment or total workflow runtime.

Final CDS/GFF export likewise streams `models.json`, retaining accepted models
only, while preserving producer-receipt and source verification. On a 4.76 GB
Ancistrocladus result, every trial retained the same 379 models with identical
canonical JSON SHA256. One warmup plus three alternating measured trials per
method gave a median 42.7 s / 12.1 GiB peak RSS for full-array loading versus
25.5 s / 38 MiB for streaming. These figures measure JSON loading/filtering,
not the whole export, genome validation or alignment.

Atomic JSON output streams the same sorted, indented document through a 1 MiB
buffer. Temporary files are removed after encoding, writing or publication
failures, and immutable plans retain their existing comparison semantics.
This avoids keeping an additional whole-document string and encoder chunk list
in memory. During the 23-species run, the Simmondsia worker reached 176.6 GiB
virtual-memory peak while materialising a 32.2 GB JSON result before writing it.

A separate normal-SIF benchmark used 2,048 actual Simmondsia records, one warmup
and three trials per method. All output bytes were identical. Median peak RSS
was 66.0 MiB for full-document encoding and 38.8 MiB for buffered streaming;
the loaded input itself used about 39 MiB. Median wall times were 0.238 s and
0.247 s, so this benchmark demonstrates lower additional memory, not a timing
speedup. The full-run peak is an observation, not a controlled comparison.

## Streaming rescue audit tables

The rescue producer writes `quality_flags.tsv` and `audit.tsv` from iterators
instead of materialising another row list for each table. Model order, unresolved
candidate order, columns and CSV quoting remain unchanged. The validated model
collection still resides in memory; this change reduces table-writer overhead.

A normal x86_64 SIF benchmark repeated 2,048 frozen Simmondsia records 128 times
(262,144 model rows). One warmup and three alternating measured trials compared
the exact production writer statements. Both complete output files were byte
identical, including an unresolved fixture containing tabs, quotes and newlines.
Median process peak RSS fell from 170.9 MiB to 56.8 MiB. The loaded sample is
included in RSS, while input loading, AST extraction and process startup are
outside writer wall time. Legacy wall times ranged from 5.55 to 7.21 seconds
and streaming times from 4.90 to 5.78 seconds on the shared server. These
measurements support a memory reduction for this serialization workload,
without establishing whole-producer, cold-cache or 500-species performance.

Use immutable source and sample directories and the sample's full SHA256:

```bash
GENEGALLEON_SIF_EXTRA_BINDS=$'/data/frozen_source\n/data/model_sample' \
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_rescue_tables.py \
  --source /data/frozen_source --candidate "$PWD" \
  --models /data/model_sample/models.json --sample-sha256 "$SAMPLE_SHA256" \
  --repeat 128 --output workspace/output/rescue_table_benchmark
```

The output directory must be new. The benchmark retains implementation hashes,
sample identity, commands, trial RSS and timing measurements, output hashes and
row-count validation. A changed sample or different output bytes fails the run.

## Compact model-store envelope encoding

The compact writer validates and encodes each final candidate envelope once,
then reuses the same immutable bytes for its accepted overlay. Explicit partial
matching retains its distinct ordinal-free encoding. Direct helper validation,
record bounds (including the JSONL newline), collision checks, output order,
compression and verified reader contracts remain unchanged.

A bounded normal x86_64 SIF benchmark used one CPU, 3,000 synthetic loci
(6,000 donor-specific models), 750 accepted models and 150 revisions. Each
condition had one warmup per method followed by six alternating trials, three
per method. All published member bytes matched, including SQLite databases,
compressed shards and manifests; an independent baseline reader restored every
model, accepted path, partial record and revision.

| Codec and partial mode | Previous median writer time | Single-pass median writer time | Observed reduction |
|---|---:|---:|---:|
| gzip, automatic | 1.182 s | 0.934 s | 20.9% |
| zstd, automatic | 1.018 s | 0.891 s | 12.5% |
| gzip, explicit reordered/repeated | 1.554 s | 1.263 s | 18.8% |
| zstd, explicit reordered/repeated | 1.229 s | 1.100 s | 10.6% |

Writer timing includes compression, SQLite, checksums, fsync and publication;
input decoding, startup and readback validation are excluded. This is a storage
component measurement on a shared server, with only three measured trials per
method and condition. The zstd automatic ranges overlapped (previous
0.949–1.047 s; single-pass 0.891–1.199 s). All process high-water memory readings
were 153.3 MiB and already reached before the timed writer, so this measurement
does not resolve writer memory differences or establish a memory reduction.
These results do not predict whole-producer or 500-species runtime.
