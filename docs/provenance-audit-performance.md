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
`--workers 1` selects serial operation. Full audits and the audited database
pipeline retain at most 16,384 digest entries and logical-source resolutions in
RAM. Resolution uses an LRU cache so later families can reuse their own sources
after the cache fills. Once the digest inventory exceeds that RAM threshold,
all digest identities are indexed in a private temporary SQLite database;
the RAM cache keeps only recent entries. Final validation streams the complete
inventory in insertion order without copying the full entry/guard lists.
Every unique source is still read initially and rehashed before publication.

The digest index has a 1 GiB disk budget and a two-million-reader bound.
Source guards retain their separate 4,096-entry RAM threshold, 512 MiB spill
budget and two-million-reader bound. Reader callables remain in process, so
total memory still grows with the number of sources; the RAM cache limit is
not a whole-process memory cap. Indexes are deleted on success or failure and
are never reused across runs. Storage failures, inventory loss and capacity
overflow fail the whole audit instead of dropping checks or reporting families
as corrupt. Small read-only queries retain their existing finite capacity.
Repeated reads reject a changed source identity instead of rebinding its digest.
The audit holds up to eight ZIP readers per thread and closes them at the fence.

Re-registering a source guard checks its identity and retains its first reader.
That reader belongs to the original query snapshot; replacing it with a later,
short-lived fingerprint store caused per-source lock/epoch checks after that
store closed. A 128-family diagnostic reduced final-stage epoch reads from
1,564 to 28 while preserving all 1,538 content revalidations and 2,062 logical
guards. Initial/final archive epoch fences, fresh content reads and late-mutation
rejection remain required. The audited-database benchmark records `index_epochs`
alongside ZIP opens and metadata-header checks to expose this regression.

`*.tree_plot.args.json` files are renderer argument outputs, not provenance
manifests. The audit excludes this known suffix only when selecting manifests;
the files remain visible through the logical store and its inventory fence.
Their owning `*.tree_plot.json` manifests still check their bytes and presence.
Other malformed manifests, including missing `family_id` or `step`, remain
errors. This works with existing raw, ZIP and metadata-SQLite storage without
moving or regenerating saved files.

## Digest spill qualification, 2026-11-04

Measurements use the same Linux aarch64 / Python 3.12.15 Docker image
`sha256:f0632dbd422a3ee87e3c2a82c3612cf64413a06c932140e69c178a16ef73b708`,
limited to four CPUs and 8 GiB. OS caches were not forcibly cleared. The image's
freshness check fails against current runtime inputs/upstream revisions; these
are Python behavior and benchmark results in that known image, not current
upstream-runtime or SIF compatibility evidence.

The paired audited-database fixture used 512 families, four manifests and twelve
unique artifacts per family, with ZIP artifacts and SQLite-packed metadata.
Against `b8996c14df2074c48cf5b0fb7465d790c21fb235`, one discarded warmup per
revision and three alternating fresh-process trials gave median times of
16.05 / 9.28 seconds (42.2% less time) and peak RSS of 179.53 / 177.39 MiB.
Epoch reads fell from 6,201 to 57. This fixture stays below the digest spill
threshold, isolating the repeated source-guard snapshot checks. All complete
logical database, audit-report and manifest-inventory
SHA256 values matched; each trial retained all 6,146 final content reads and
8,206 source guards.

Full measurements, counters, source hashes and validation commands are retained
in [the machine-readable evidence](benchmarks/provenance-audit-spill-2026-11-04.json).

Use the same runner, runtime and CPU/memory budget for both revisions:

```bash
python -B workflow/tests/benchmark_audited_database.py \
  --baseline-support /path/to/baseline/workflow/support \
  --families 512 --samples 3 --warmups 1 --layout packed \
  --output /tmp/gg-audited-db-paired.json
```

The expanded qualification includes renderer arguments and deliberately crosses
the former 250,000-digest ceiling with independent artifacts, rather than only
increasing the number of manifests referencing shared sources:

A complete final-code run reused those immutable synthetic inputs in a fresh
temporary root, with copied lock namespaces omitted. It published the audited
database with 18,811 rows each in `tree` and `branch`, 94,055 `current` audit rows
and 18,811 `legacy_untracked` rows under the unchanged legacy policy. All
319,789 unique digests were reread before publication and all 413,860 guards
were retained. The RAM digest cache held 16,384 entries. The single qualification
took 556.37 seconds with 2,302 MiB peak RSS, 13 metadata-header checks and 73
archive-epoch reads. This is synthetic execution/scale evidence, not a real
scientific analysis, a before/after large-run speedup or a whole-process RAM cap.

```bash
python -B workflow/tests/benchmark_audited_database.py \
  --families 18811 --files-per-family 14 --plot-arguments \
  --samples 1 --warmups 0 --layout packed \
  --output /tmp/gg-audited-db-large.json
```

For an isolated memory/time comparison, run the same revalidation fixture in
alternating fresh processes. Both revisions use source-guard spill; enable
`--spill` only for the candidate's new digest spill. Compare the complete
source-set SHA256 and every initial/final read counter as well as peak RSS,
collection time (including private fixture writes) and final validation time:

For 80,000 real 8 KiB files, one discarded warmup per revision and three
alternating trials reduced median peak RSS from 183.13 to 140.33 MiB (23.4%).
Median collection-plus-validation wall time increased from 20.52 to 24.96
seconds (21.6%): SQLite spill trades additional work for lower RAM use and
larger inventories. All complete source-set SHA256 values matched; both
revisions read every file twice, checked all 80,000 guards and hashed
1,310,720,000 source bytes. This isolated result is not a whole-workflow memory
claim or an end-to-end speedup.

```bash
python -B workflow/benchmarks/benchmark_audit_revalidation.py \
  --support-root /path/to/baseline/workflow/support \
  --files 80000 --workers 4 --latency local
python -B workflow/benchmarks/benchmark_audit_revalidation.py \
  --support-root workflow/support --spill \
  --files 80000 --workers 4 --latency local
```

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
`cache_entries` counts the current RAM digest cache; `digest_spilled_entries`
counts the complete disk-indexed digest inventory (including hot entries), or
zero before spilling. `unique_digest_reads` and `digest_revalidation_reads`
continue to count initial and final content reads, respectively.
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
