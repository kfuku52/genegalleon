# Full provenance audit performance and progress

## Revalidation with uneven read latency

Run-scoped content and logical-source revalidation keep an ordered rolling
window of at most eight submissions per existing reader. A completed source
allows the next source into that window without waiting for an entire batch.
Reader counts, inventory limits, the content-before-guard barrier, every final
content read and every identity check are unchanged. The first failing input
still determines the reported error; all submitted readers join before it
escapes. No result from an earlier observation replaces a read.

A Linux aarch64 / Python 3.12.15 Docker comparison against `7775c7db` used
four CPUs, 8 GiB, 1,024 real 8 KiB sources and 1,024 logical guards. One warmup
per revision/case preceded three alternating fresh-process trials. The tail
case adds 0.5 ms per final read and 20 ms instead for each 32nd source/guard.
This models uneven storage latency, not an end-to-end HPC measurement.

| Final-read latency | Before (s) | After (s) | Median peak RSS before/after (MiB) |
| --- | ---: | ---: | ---: |
| local, no added delay | 0.12673 | 0.12998 | 22.72 / 23.39 |
| modeled uneven latency | 2.16572 | 1.70619 | 22.79 / 23.41 |

The modeled case took 21.2% less time (1.27x). Local timings overlap and do not
establish a local-storage speedup. Every complete source-set SHA256 matched;
each content source was read twice and every logical guard was checked once.
OS caches were not forcibly cleared. Both revisions used the same existing
runtime; this standard-library benchmark is separate from runtime freshness
and SIF compatibility qualification.

Run the same current runner against immutable baseline/candidate support trees
under the same CPU/memory limits, alternating revisions and discarding warmups:

```bash
python -B workflow/benchmarks/benchmark_audit_revalidation.py \
  --support-root /path/to/workflow/support --latency tails --workers 4
```

## Final metadata fence

The final cross-query namespace fence retains every fresh `lstat` signature,
including absent paths. With a usable `GG_TASK_CPUS` hint it can overlap at most
four metadata readers. Fences with fewer than 64 paths, unconfigured hints and
calls inside existing workers stay serial. The first 32 signatures are checked
serially; when their average cost is below 100 microseconds the rest also remain
serial to avoid thread overhead on warm local storage. The sample is part of the
fence: no signature is skipped, duplicated or reused from an earlier query.

Parallel work uses batches of at most eight paths per reader. Ordered results
preserve the first failing path, and all workers join before errors escape.
Regular-file replacement, removal, absent-to-present changes, symlinks, nested
workers and CPU bounds remain covered by the existing workflow API test module.
This changes only signature-read scheduling; content rehashes and archive fences
remain mandatory.

A Linux aarch64 / Python 3.12.15 Docker fixture compared `7d8864d7` with this
implementation using the same 4,000 files, four CPUs and 8 GiB. One discarded
warmup per revision and three alternating fresh-process trials per scenario
produced these medians. Every complete fresh-signature-set SHA256 matched.

| Per-signature modeled delay | Before (s) | After (s) | Peak RSS before/after (KiB) |
| --- | ---: | ---: | ---: |
| none | 0.01479 | 0.01406 | 26,496 / 26,596 |
| 1 ms | 4.9765 | 1.3925 | 26,544 / 27,492 |

OS caches were not forcibly cleared. The warm local measurements overlap and
do not establish a speedup. The latency fixture improved by 3.57x with about
1 MiB more peak RSS; this is not an end-to-end audit or measured HPC speedup.
Heterogeneous metadata latency can leave the serial mode selected after a fast
sample; correctness and the complete signature fence are unchanged.

Create a private fixture once, then alternate support roots with the same
current benchmark runner and `GG_TASK_CPUS=4`:

```bash
python -B workflow/benchmarks/benchmark_observation_fence.py \
  --workspace /tmp/gg-observation-fixture --create
GG_TASK_CPUS=4 python -B workflow/benchmarks/benchmark_observation_fence.py \
  --workspace /tmp/gg-observation-fixture --support-root /path/to/workflow/support \
  --metadata-delay-ms 1
```

## Large inventory guard lookup

The temporary SQLite source-guard store now keeps its key index separate from
large inventory identities. A `WITHOUT ROWID` table placed multi-megabyte
identities in the key B-tree, causing unrelated small guard lookups and inserts
to repeatedly traverse overflow pages. The ordinary rowid table retains unique,
non-null keys, unique reader indices, type-preserving identities, all final
reader checks, the 512 MiB disk budget and the two-million-reader bound. This is
private query storage; no persisted workflow format or scientific output changes.

A paired Linux aarch64 / Python 3.12.15 Docker fixture used 107,841 inventory
members and 1,000 smaller guards, four CPUs and 8 GiB. Against `8cf1e617`, one
discarded warmup per revision and three alternating trials in fresh processes
measured median insertion time 1.3383 / 0.03385 seconds. Median process read
characters fell from 9,498,598,116 to 32,116 and read syscalls from 2,322,002 to
2,002. These counters include OS-cached reads, not physical disk traffic.
Median peak RSS was 65,708 / 65,672 KiB. Every complete identity set SHA256 and
reader result matched; OS caches were not forcibly cleared. This isolated
lookup measurement is not an end-to-end audit or HPC speedup claim.

Run the same current benchmark runner in separate processes with old/new
`--support-root` paths, alternating revisions under the same runtime limits:

```bash
python -B workflow/benchmarks/benchmark_guard_spill.py \
  --support-root /path/to/workflow/support
```

## Scalar subdirectory preflight reads

Read-only preflight profiling found 144,073 metadata calls consuming 37.79 of
80.0 sampled seconds; file reads consumed another 19.49 seconds. The sample
was deliberately interrupted for diagnosis and is not a completion proof.
Subdirectory inputs were still serial, and collection traversal resolved each
member again before the query's existing source resolver. Scalar subdirectory
digests now overlap at most four readers within `GG_TASK_CPUS`; callers already
running in audit worker threads remain serial to avoid multiplying concurrency.
Member traversal shares the original query resolver. All content reads, final
rehashes, raw overrides, inventory/epoch guards and failure propagation remain.

A paired Linux aarch64 / Python 3.12.15 Docker benchmark used four CPUs, 512
families, two provenance steps, a 1 MiB shared input and a synthetic 1 ms
per-member I/O delay. Against `e266142db9cfffa238f5b50ce685ea10b2bb40d5`,
one discarded warmup and three alternating trials per layout gave these medians.
Every logical digest, byte count and 512 final content revalidations matched.

| Layout | Before (s) | After (s) | Peak RSS before/after (KiB) |
| --- | ---: | ---: | ---: |
| raw | 1.8063 | 1.1090 | 44,360 / 44,560 |
| ZIP | 1.8594 | 1.1282 | 45,132 / 45,468 |
| mixed | 1.7476 | 1.2938 | 45,056 / 45,308 |

These fixture measurements do not establish end-to-end HPC speedup. Each trial
starts a fresh query; the OS cache is not forcibly cleared. More overlapping
readers increase simultaneous I/O and a small amount of memory. Reproduce with
the existing fixture creator below and this operation, alternating old and new
support roots while retaining the same current benchmark runner:

```bash
python -B workflow/benchmarks/benchmark_provenance_audit.py \
  --workspace /tmp/gg-audit-fixture --support-root /path/to/support \
  --operation subdir-digest --subdir stat_branch --query-observation \
  --workers 4 --member-delay-ms 1
```

`artifact_provenance.py audit` preserves the TSV columns, sorting, scientific
hashes, stale-policy classifications and CSUBST branch checks. It hashes shared
declared sources once during classification, then rereads every unique source
and manifest before publishing the report. Persistent digests from earlier runs
cannot substitute for this audit's content reads. Final content checks are needed
because shared-filesystem attributes alone did not detect a rapid same-size
rewrite in the actual SIF regression fixture.

One producer read lock and an initial/final metadata fence replace repeated
archive-index scans. Raw/ZIP source selection, pending-index, symlink, tombstone
and generation checks remain enforced. Thread pools and their submission windows
are bounded; the default is up to four workers, capped by `GG_TASK_CPUS`.
`--workers 1` selects serial operation. Run-scoped caches hold at most 250,000
sources/guards; either capacity overflow fails instead of dropping a check.
Repeated reads reject a changed source identity instead of rebinding its digest.
The audit holds up to eight ZIP readers per thread and closes them at the fence.

## Evidence, 2026-09-30

Synthetic fixtures repeat a 24 MiB shared table across families and a 16 KiB
family table across 15 manifests. Each run uses Linux Python 3.12.14 in the same
`local/genegalleon:dev` Docker runtime, limited to 4 CPUs/4 GiB. The baseline is
GeneGalleon `63efea916bee99eeae98d2ae30ead62b374b91f4`. Separate warmups precede
three measurements per version; baseline trials use an ABBAAB sequence.
Final candidate trials were repeated three times after the capacity and
source-identity guards were tightened. Run-scoped RAM caches start empty; the
baseline retains its persistent digest
cache warmed by prior trials. The new full audit intentionally ignores those
persistent digests. OS caches are not forcibly cleared.

| Layout | Families / manifests | Baseline median | New median | Speedup |
|---|---:|---:|---:|---:|
| Raw | 1,000 / 15,000 | 45.54 s | 8.16 s | 5.58× |
| ZIP | 300 / 4,500 | 90.91 s | 2.64 s | 34.44× |
| Mixed raw/ZIP | 300 / 4,500 | 63.57 s | 2.55 s | 24.90× |

Every baseline/candidate TSV SHA256 matches within its fixture. These are
complete fixture audits; the fixture has no CSUBST branch data. Existing branch
identity regressions exercise those unchanged scientific checks separately.
These numbers are not an ETA for a running scientific job.

Seven 6,815-family / 102,225-manifest runs took 60.30–65.40 seconds with a maximum
262.65 MiB peak RSS. Every TSV SHA256 matched. Instrumented atomic progress
publication consumed at most 0.27% of wall time. In the six paired trials,
progress-on/off medians were
64.33/61.77 seconds, with overlapping ranges: this jitter does not independently
establish an end-to-end overhead bound. The 0.27% measurement covers snapshot
publication rather than total filesystem traffic or all audit bookkeeping.

Reproduce private fixtures outside any scientific workspace:

```bash
python -B workflow/benchmarks/benchmark_provenance_audit.py \
  --workspace /tmp/gg-audit-benchmark-raw --create
python -B workflow/benchmarks/benchmark_provenance_audit.py \
  --workspace /tmp/gg-audit-benchmark-zip --create --families 300 \
  --layout zip --support-root "$PWD/workflow/support"
python -B workflow/benchmarks/benchmark_provenance_audit.py \
  --workspace /tmp/gg-audit-benchmark-raw --support-root workflow/support \
  --workers 4 --label trial-1
```

Fixture creation refuses an existing destination. For baseline measurements,
extract its immutable `workflow/support` tree and pass that directory; omit
`--workers` because the old CLI has no such option. Keep the runtime, CPU/memory
limits, fixture, warmups and cache conditions comparable.

## Machine-readable observations

The default `--progress-interval 10` writes atomic owner-private
`REPORT.tsv.progress.json`; zero disables progress snapshots. Phases are
inventory, manifest hashing, legacy inventory, branch identity, source
revalidation and report publication. Counters and ETA describe only the current
phase; heartbeat time and last advancement time are separate.

`REPORT.tsv.result.json` records source-closure, manifest-inventory and report
SHA256, status counts, exit code, elapsed time and bounded cache counters.
`digest_bytes_read` counts declared-source hash stream bytes, including the final
reread; it excludes manifest JSON and physical/compressed filesystem traffic.
`progress_write_ns` measures progress publication time.

Observed jobs additionally write `progress-artifact_audit.json` and
`result-artifact_audit.json` beside their exact `run.json`.
The optional [workflow API](workflow-inspection-api.md) progress capability
returns these as attempt evidence, with strict bounds and identity validation.
Neither progress completion nor a successful audit result proves that the full
workflow or database generation finished. Existing jobs retain their immutable
source and observation coverage.

## Whole-store fingerprints, 2026-10-17

Database reuse and HGT transfer-context preparation also fingerprint the whole
logical family store before their stage begins. This previously resolved and
hashed every member serially. Whole-store fingerprints now use up to four I/O
workers, bounded by `GG_TASK_CPUS` (one when unset). Direct Python callers may
request `gene_family_store_digest(root, workers=1)` for serial execution.
Subdirectory fingerprints remain serial inside already parallel family audits.

Every logical member still participates in the same sorted SHA256 stream, byte
count and member count. Submission windows contain at most eight tasks per
worker. Workers retain read-only, digest-proof and runtime-path contexts; all
workers stop before the store closes snapshot ZIP readers. Pending archive
updates, final index/content fences, live overrides and read failures still
invalidate the complete operation. This does not reuse older verification
receipts or alter a running job's pinned source.

Private 1,000-member / 16,384,000-byte fixtures used the same Linux/amd64
GeneGalleon Docker image on macOS ARM64, Python 3.12, 4 CPUs and 16 GiB. The
baseline is `3b2ad80510a0c6aa02292c72d7927627eefa0d54`. Each scenario used a
fresh process, one discarded warmup and three measurements; no OS cache was
forcibly cleared and read-only fingerprints used no persistent digest cache.

| Layout / per-member delay | Serial median | Four-worker median |
| --- | ---: | ---: |
| Raw / none | 0.331 s | 0.337 s |
| Mixed raw/ZIP / none | 0.398 s | 0.426 s |
| Raw / synthetic 1 ms | 2.552 s | 0.546 s |
| Mixed raw/ZIP / synthetic 1 ms | 2.499 s | 0.626 s |

All 32 measurements, including warmups, returned the same SHA256, byte count
and member count. Peak parent RSS was below 56 MiB. The 4.0–4.7x speedup models latency between member reads; it is not a measured
HPC speedup or a running-job ETA. The small warm local trials had overlapping
ranges and no established benefit; their medians increased by 6–28 ms. Compare
on the target filesystem before extrapolating.

Use the existing fixture creator above, then run each revision against the same
private workspace and support tree:

```bash
python -B workflow/benchmarks/benchmark_provenance_audit.py \
  --workspace /tmp/gg-audit-benchmark-raw --support-root workflow/support \
  --operation store-digest --workers 4 --member-delay-ms 1 --label trial-1
```

Omit `--member-delay-ms` for a real filesystem measurement. The optional delay
is confined to this synthetic benchmark; production fingerprints never add it.
