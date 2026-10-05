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
Validation imports compare every destination file with the source validation
proof, including the rewritten summary's target format proof. Aliased output
paths (including symlinks and hard links to donor or raw files) are rejected
before copying; unused directory settings do not block an import.
Summary rewriting parses a freshly hashed in-memory source, so changing and
restoring the source during parsing cannot introduce unverified paths.
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
Docker inventories use the currently resolved immutable image ID and collect
with that ID, so replacing a tag invalidates its cached inventory.

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

On 2026-10-05, comparison with v0.8.132 in the same qualified Docker runtime
(Linux arm64, Python 3.12.14) used a 1,024 MiB nominal synthetic genome, warm
filesystem caches, one warmup and three alternating measured trials:

| Operation | Baseline median | Updated median | SHA-256 bytes, baseline → updated |
| --- | ---: | ---: | ---: |
| Initial local staging | 1.469 s | 1.014 s | 3,222,271,619 → 2,148,181,556 |
| Native checkpoint import | 7.602 s | 5.323 s | 12,876,516,732 → 8,584,346,448 |

Both comparisons produced identical fingerprints. Median peak process RSS
stayed about 66–67 MiB. Hash counters in the import comparison cover the parent
process, including copy-stream hashes; the unchanged metadata subprocess also
reads inputs. These results establish approximately 31%/30% lower times for
these isolated operations, not an HPC or whole-project speedup. Other staging
phases were similar within run-to-run variation. No production speedup is
claimed for inventory caching or fs-verity; the qualified test filesystem did
not provide enabled fs-verity images.

A subsequent audit compared the hardened importer with v0.8.134 using the same
runtime and 1,024 MiB method. Output fingerprints matched, with median times of
5.127 s before and 4.990 s after; this small difference is not a new speedup
claim. Hash reads increased by only the 124-byte summary shard, from
8,584,346,448 to 8,584,346,572 bytes. Copy bytes were unchanged and peak RSS
remained about 67 MiB.
