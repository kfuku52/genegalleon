# Rescue performance

The optimization keeps each synteny interval independent, caps alignment threads
at the worker CPU budget and preserves candidate order. Whole-genome query
deduplication uses exact protein sequence equality; results are expanded back
to every original candidate before validation. The comparison cache binds
actual BED/protein inputs, comparison options and the relevant toolchain.
See [gene-model rescue](gene-model-rescue.md) for controls and cache guarantees.

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

## Reproducible bounded benchmark

Use the repository runtime wrapper with completed real rescue evidence. It
checks sampled files against the producer's hashed receipt and rechecks inputs
afterward. Every trial compares SHA256 of all restored GFF bytes, including
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
