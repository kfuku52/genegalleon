# Input-aware resources

GeneGalleon describes the scientific input and divides the allocated CPU/memory
budget among its tools. kfauto selects a scheduler request before a new JobSpec
is sealed. This feature never enlarges a running allocation.

## Profile the intended work

Run the standard-library profiler once before planning, with the exact inputs
and scientific settings of the next execution. Do not glob all historical
workspace artifacts. For an array, describe one representative member and pass
the array size separately; heterogeneous members require separate profiles.

```bash
python workflow/support/resource_profile.py \
  --workflow gg_genome_evolution \
  --fasta workspace/input/species_protein/species_a.fa \
  --fasta workspace/input/species_protein/species_b.fa \
  --settings analysis-settings.json \
  > workspace/input/resource-profile.json
```

Repeated `--fasta`/`--fastq` arguments support gzip. FASTQ currently requires
four-line records. `--species` overrides the default file count (or sequence
count for `gg_gene_evolution`). Record the unfinished execution scope using
repeated `--stage` arguments, using the exact names printed by `gg_step_start`.
Without a nonempty explicit scope, timing is inspection-only. Automatic learning
requires the started stages to match that scope exactly and none of those stages
to be skipped. No-work and partial retries cannot replace comparable history.
The profiler records counts, total residues,
largest sequence, input content hashes, and settings identity. It does not
decide which scientific stages may be skipped. Refresh the profile when inputs,
settings, or remaining work change. kfauto checks file size/mtime without
rescanning sequence contents on each observation; this is an acceleration hint,
not a replacement for GeneGalleon's artifact provenance validation.

## Genome-evolution parallelism

The allocation is still `GG_TASK_CPUS` and `GG_MEM_TOOL_GB`; no process is granted
additional scheduler cores by changing an environment variable.

* `orthofinder_algorithm_threads=auto`: max(1, allocated cores / 8), additionally
  capped by tool memory divided by `orthofinder_memory_gb_per_thread` (default 4).
  Search parallelism continues to use the full allocated CPU budget.
* `genome_parallel_jobs=auto`: per-gene concurrency is bounded by both allocated
  cores and tool memory divided by `genome_parallel_memory_gb_per_job` (default 2).
  The per-worker IQ-TREE/Notung memory budget uses this concurrency.
* Explicit values remain available and cannot exceed the allocation or the
  corresponding estimated memory cap. BUSCO retains its existing independent
  concurrency/memory controls.

These memory estimates are initial policy values, not guaranteed peak usage.
Use larger per-worker estimates for large alignments. The normal registered
environment overrides expose the new settings under `GG_GENOME_EVOLUTION_`.

## Execution evidence

The standard `exec` container bridge records each core attempt under
`workspace/output/resource_metrics/attempt-<id>.json`. It runs the real command
with inherited stdin/stdout/stderr and preserves failure and signal status.
Non-core helpers and legacy shell-only site adapters keep their existing path.
`GG_RESOURCE_METRICS=0` disables the wrapper for a controlled comparison.

Records include wall time, child CPU time, allocated CPU/total memory, exit code,
and the OS-reported maximum **per-process** RSS. This RSS is not simultaneous
aggregate memory across parallel workers and must not be used to reduce an
array's total memory budget. Stage boundaries record elapsed time and reaped
child CPU deltas from the owning shell; overlapping/background stages are not
separately attributed, and stage peak RSS is not available.

For comparable history, supply `GG_RESOURCE_PROFILE` (a container-visible path),
`GG_RESOURCE_RUNTIME_ID` (the exact runtime identity), and
`GG_RESOURCE_SERVER_ID`. kfauto's adaptive scalar adapter supplies these.
Successful classified attempts with verified execution scope atomically publish a latest-per-CPU index
under `resource_metrics/history/`, keyed by workflow/input/runtime/server/stage.
This enables the next resource selection without scanning all historical jobs.
Records include `expected_stages`, `skipped_stages`, and `learning_eligible`.
Older records without this eligibility evidence are excluded from learning.
Without the identity fields, timing remains useful for inspection but is not eligible for
automatic resource learning. Completed output equality and actual 8/16/32-core
scaling still require representative scientific benchmarks; no speedup is
asserted merely from a larger allocation or a synthetic test.

OrthoFinder's separate search/analysis parallelism is described in its
[official documentation](https://github.com/davidemms/OrthoFinder/blob/master/README.md#parallelising-orthofinder-algorithm).
