# Input-generation performance

Native array preparation freezes the plan, verifies local inputs and shared
resources, and checks the donor's sealed plan/settings. Each worker imports only
its own species after acquiring its target task lock. Existing prepares that
already imported checkpoints remain usable. `input_generation_stage_resume.py
import` without `--task-index` still imports the whole cohort; `check-source`
checks the donor contract without copying species outputs.

Imports use the workflow's shared-filesystem namespace locks: a shared donor
phase lock excludes prepare/finalize, and a shared donor task lock excludes
writers of that species. Different species can import concurrently. The target
has the same phase/task protection. Native workers pass their existing target
ownership token rather than reacquiring their own exclusive lock. Never mix
these locks with an older flock-only runtime in one workspace.

A current format checkpoint is sufficient without reading unrelated completion
receipt outputs. Fresh source hashes are reused only inside the import; sources
are fully rehashed after copying. Copies hash their stream, and destination
files are independently rehashed and compared before publishing the checkpoint.
Independent source-gene ownership and mapping QC remain required for validation
reuse. Formatting and validation contracts and BUSCO lineage settings are
unchanged.

For bound local staging, the binding step uses the current invocation's full
preflight read with file-identity fences. Gzip validation still runs, and
publication independently rehashes the inputs after intervening work. No
persistent size/mtime checksum cache certifies mutable research inputs.

## Timing records

Native input generation writes advisory JSONL under
`output/input_generation/tmp/performance/<core-pid>/`. Planning, staging,
checkpoint import and genome reconstruction/indexing record wall time, inclusive
SHA-256/copy byte counters, and process/child peak RSS, with timestamp/host/job tags.
RSS is the process-lifetime
high-water mark, not an isolated phase allocation. Each process has its own file;
these records never authorize completion or reuse. A telemetry write failure
warns without changing scientific results. Other support-script invocations can
opt in with `GG_PERFORMANCE_DIR`.
BUSCO/ETE dataset preparation also records wall time, including download/lock
waits; those shell phase records explicitly omit RSS and byte counters.

Version-image timing is recorded under `output/versions/performance/`. Version
inventory caches require a fresh image identity, repository version, inventory
script hash, inspect result, runtime command/bindings and environment digest.
Environment values are never written to cache metadata. Diagnostic job IDs are
excluded; tool-resolution environment changes invalidate the cache. Every call
also writes a dated host/job record under `output/versions/run.*.log`.

The default private inventory cache is `downloads/version_inventory`. Set
`GG_VERSIONS_CACHE_DIR` on the entrypoint host to share a cache; it must be owned
by the executing user with mode 700. Different bindings or relevant environments
collect separately. Only successful inventories are cached, their content hash
is checked on reuse, and collection failures still fail the entrypoint.

Ordinary SIFs are fully SHA-256 hashed on every invocation. An already-enabled
Linux fs-verity image can use its kernel-enforced content identity in constant
time. GeneGalleon only measures it and never enables fs-verity or changes an
image. The verity identity is explicitly labelled and is distinct from a normal
file SHA-256. See the [Linux fs-verity documentation](https://docs.kernel.org/filesystems/fsverity.html).

## Comparisons

Use the same qualified container and filesystem, without concurrent tests or
builds. Both harnesses alternate baseline/current runs, discard one warmup,
report three measured trials by default, and require identical output
fingerprints. These synthetic workloads measure staging and checkpoint transfer,
not whole-project elapsed time or BUSCO runtime.

```bash
GENEGALLEON_DOCKER_EXTRA_BINDS=/path/to/baseline \
  bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_input_staging.py \
  --baseline-support /path/to/baseline/workflow/support --genome-mib 256 \
  --output tmp/staging-comparison.json

GENEGALLEON_DOCKER_EXTRA_BINDS=/path/to/baseline \
  bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_stage_resume.py \
  --baseline-support /path/to/baseline/workflow/support --genome-mib 256 \
  --output tmp/resume-comparison.json
```

The resume harness includes streaming-copy hashes in SHA-256 counters and reports
copy bytes separately. A hashed copy contributes to both metrics; do not sum
them as disjoint I/O. Fingerprints exclude lock-owner diagnostics. Genome indexes
already live in task scratch and are reused within a normalizer. Broader index
sharing should follow measurements of the new indexing phase rather than
changing independent gene-selection parsers.
