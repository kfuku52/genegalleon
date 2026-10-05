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
receipt outputs. The worker resolves metadata and imports in one `start-worker` invocation.
Fresh hashes are shared only within separate pre-transfer and post-transfer
batches, including physical hard links; sources are fully rehashed after copying. Copies hash their stream, and destination
files are independently rehashed and compared before publishing the checkpoint.
Validation imports compare every destination file with the source validation
proof, including a hash of the exact newly serialized summary bytes. Aliased output
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
Cache initialization serializes concurrent creators and explicitly restricts
new directories when an NFS server overrides the requested creation mode.
Existing caches with unsafe permissions remain an error.
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

## Scoped worker checks and streaming genomes

`check-prepared --task-index INDEX` verifies the sealed plan/settings, shared
resources, unknown marker entries and that worker's staging receipts. Other
workers' known staging receipts are deferred to the unchanged global
prepare/finalize checks. Default `prepared` and whole-cohort imports still check
all marker files. Worker imports select the corresponding donor index.
Metadata fields use one JSON read per group and NUL-delimited shell assignments,
with no `eval`; null/missing fields and error propagation retain their behavior.
The native worker's 18 individual field reads become four grouped reads.

FASTA genome formatting reads bounded sequence chunks and writes through the
existing seqkit compressor and atomic publisher. Header normalization, organelle
exclusion, empty/duplicate IDs, Unicode whitespace, archive member selection and
80-column input wrapping are preserved. CDS selection and the GBFF path retain
their existing record interface. Seqkit still retains a chromosome internally.

A 2026-10-05 comparison against v0.8.139 (`1140bbd`) used the immutable Docker
runtime above, Linux arm64 / Python 3.12.14, warm filesystem caches, one warmup
and three alternating fresh processes. No other tests/builds ran concurrently.
All compared outputs matched; worker fingerprints cover genome/CDS/GFF,
fx2tab/BUSCO outputs and complete worker QC results. The implementation identity
field deliberately changes with code and is excluded from scientific QC parity.

| Operation | Before median | After median | Scope |
| --- | ---: | ---: | --- |
| Checkpoint import, 256 MiB nominal genome | 1.632 s | 0.838 s | Full checkpoint/output fingerprint |
| Genome formatting, 64 MiB single chromosome | 1.343 s | 0.440 s | Real seqkit, one thread; decoded FASTA identical |
| Native BUSCO-lineage-change worker, 64 MiB genome | 2.224 s | 1.668 s | Real core/validators; fake BUSCO/seqkit/Rscript |
| Worker preparation check, 560 staged species | 0.0235 s | 0.0126 s | Two staging receipts per species |
| Worker preparation check, 639 staged species | 0.0273 s | 0.0198 s | Two staging receipts per species |

The genome benchmark's median sampled **whole process-tree** peak RSS fell from
484,184 to 354,700 KiB (about 473 to 346 MiB, 27%). Sampling at 10 ms can miss
short peaks; this does not establish a safe scheduler memory request. The
parent-only drop is not reported as a whole-worker saving.

Import parent SHA counters fell from 2,146,091,980 to 1,609,570,254 bytes, with
268,260,730 copy bytes unchanged. Baseline counters omit its metadata child;
the updated in-process metadata read is included. Neither parser reads nor
whole-worker subprocess I/O are represented by these counters.

In the staging fixture, per-worker full-file hashes fell from 1,122 to 4 for
560 species and from 1,280 to 4 for 639 species. Across all workers this is
628,320 to 2,240 and 817,920 to 2,556 reads respectively, plus any shared-resource
checks in real projects. Plan/settings still require full fresh reads. Timing
on production NAS, SIF and real BUSCO remains unmeasured; array concurrency,
lineages, scientific thresholds and scheduler resources were not changed.

Reproduce the genome/scoped-check/final-QC/native-worker comparisons with a
read-only baseline repository visible through the runtime wrapper:

```bash
GENEGALLEON_DOCKER_EXTRA_BINDS=/path/to/baseline \
  bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_native_efficiency.py \
  --baseline-repo /path/to/baseline --mode genome --genome-mib 64 \
  --output tmp/genome-efficiency.json
```

Other modes are `prepared --species 560` (or 639), `final --genome-mib 64`, and
`worker --genome-mib 64`. Each output path must be new. Use the existing resume
harness above with `--genome-mib 256` for the checkpoint comparison. Fixtures
are temporary; benchmark JSON is advisory and never authorizes reuse.

The optional GeMoMa refinement also loads a donor's protein FASTA once, lazily
when there is a model to align, rather than once per model. A six-model,
two-query comparison with v0.8.140 read the FASTA six times before and once
after, with identical full model/coverage/identity results. That bounded check
uses real FASTA parsing and Bio alignment with fake Java/model validation;
no full GeMoMa runtime or memory improvement is claimed.
