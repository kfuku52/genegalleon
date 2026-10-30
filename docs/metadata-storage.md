# Offline metadata collection

`workflow/support/metadata_store.py` collects completed observation JSON,
resource metrics, artifact provenance, version reports, version-inventory cache,
output-inventory journals and native ZIP index buckets into `.gg-metadata.sqlite3`
inside each existing metadata root. Scientific annotation/statistics TSV bodies,
FASTA, trees and ZIP members keep their existing formats. Family state, archive
queues and distributed locks remain under their existing ownership protocols.

Scheduler stdout/stderr use [closed batch ZIP storage](scheduler-log-storage.md)
instead. They can exceed the metadata body bound and are written by an external
scheduler, so collect them only after its output delivery finishes.

Stop **all** workspace writers, including older runtimes, before collection or
restoration. Install the complete new reader/writer runtime first and select it
explicitly for subsequent work. No online shared writable SQLite DB is used:
the collector builds locally, verifies every original body, copies and fsyncs
an immutable DB, atomically publishes it, and only then removes unchanged source
inodes. Live raw files take precedence over archived records. Closed attempt
directories remain to preserve API identities; running attempts are retained.

```bash
python workflow/support/metadata_store.py pack \
  --root workspace/output/observations --scratch /local/scratch \
  --writers-stopped --delete-originals --operation-id exact-maintenance-id
```

Repeat only the same operation ID after an interrupted retirement. Its persisted
source identities and checksums exclude new or replaced files. Preserve that DB
until the operation is verified. A new collection requires a new operation ID.
The filesystem publication assumes cooperative exclusion of nonparticipating
writers; a flag alone cannot stop an older process or a manual writer.

Readers verify application ID `0x47474d44`, user version `1`, metadata kind,
owner/mode and standalone publication. `files` preserves zlib-compressed exact
bytes, SHA-256, mode and mtime; `records` indexes workflow, attempt/family/step,
runtime, scheduler identity and resource fields. `archive_entries` indexes ZIP
locators and member hashes without importing their scientific TSV contents.
Database identity fences invalidate observations spanning a replacement or an
in-place SQLite commit. Read-only connections and the 100-byte SQLite header
change counters detect same-size commits even when filesystem timestamps have
only second precision; lookups do not hash the complete database.

Restore before running an older reader. Restoration retains the DB and refuses
to overwrite changed raw files; native provenance tombstones and newer ZIP
members remain authoritative.

```bash
python workflow/support/metadata_store.py restore \
  --root workspace/output/observations --writers-stopped
```

For `version_inventory` only, a separate `pack --cache-limit N` call after
collection expires regenerable contexts, retaining the newest N (1–1024)
64-hex log/JSON pairs. Newly collected publications receive increasing values
in the existing `generation` field. Retention orders by that generation, then
observed mtime and path, so a new publication survives equal timestamps. An
initial legacy batch uses its observed mtimes and stable path ties; original
file metadata is preserved exactly. Unknown names remain. Once activated, the producer uses
one cache lock and reuses the persisted retention policy. Historical per-run
version reports are separate lossless records. Existing lock namespaces must
not be removed while any older writer can hold their inodes.

A local 1,000-attempt fixture (five JSONs each), Python 3.13.16 on macOS arm64,
used 6,001 directory/file entries before collection and 1,005 afterward
(83.25% fewer). Public status record digests matched. Three reads took
1.014/1.027/1.053 seconds before and 1.612/1.742/1.599 after the snapshot reader
optimization; peak child RSS changed from 63.5 to 68.4 MB. This is a local
file-count/compatibility measurement, not a Shirokane throughput claim. The
tradeoff favors fewer shared-filesystem entries; SQLite reads are not claimed
to be faster. Runtime tests cover publication/retirement races, checksum failure,
restore, mixed raw/DB reads, index update/repair, provenance generations and
cache producer/consumer behavior.


Collector publication retains the same source/body checks and namespace exclusion.
Retirement pins one immutable DB snapshot and fsyncs each changed directory once,
including successful removals before an interruption. Verified collector roots
avoid repeating ancestor lookups for every source. These changes reduce repeated
filesystem operations; they do not permit packing with active workspace writers.

For an owned engineering workspace, reproduce collector measurements with:

```bash
python workflow/benchmarks/benchmark_metadata_store.py --work-dir /owned/workspace --attempts 1000
```

An optional `--support /baseline/support` selects the comparison implementation.
A local 5,000-file fixture used one 100-attempt warmup and three fresh
1,000-attempt repetitions on Python 3.13.16/macOS arm64. Collection took
3.042/2.865/2.766 seconds before and 2.229/2.093/2.095 after (median 27% lower).
Peak RSS ranges were 39.7–40.5 MB and 40.2–41.7 MB; every logical body digest
matched. These are local collector measurements, not Shirokane throughput.

Use `--layout version-cache` to include a legacy coordination tree for every
flat cache context. Flat metadata collection preserves these opaque trees and
only scans the payload level.

Closed-observation run headers are prefetched with at most 16 readers and 32
pending tasks. The existing inventory, eligibility checks, symlink rejection,
source signatures and publication checks still determine the exact collected set.
All workers join before errors return or publication proceeds.

For an inventory-only comparison, `--inventory-only --run-read-latency-ms 1`
simulates one millisecond per run-header read. Three 500-attempt/five-JSON local
repetitions on Python 3.13.16/macOS arm64 took 0.818/0.851/0.853 seconds before
and 0.115/0.114/0.122 afterward, with identical body digests. Warm inventory
medians were 0.081/0.112 seconds; complete warm collection medians were
0.999/1.098 seconds, with peak child RSS medians 33.8/35.4 MiB. This is a
controlled latency experiment and a local overhead measurement, not a cluster
throughput prediction. On Shirokane, a read-only comparison of all 366,053
eligible observation files returned the exact same membership as the original
inventory; different cache states prevent a comparable wall-time claim there.

Small metadata bodies are also prefetched during collection and the final
prepublication source verification. At most 16 reads/results, each bounded by
8 MiB, run ahead of the current record. SQLite inserts remain ordered on the
caller thread. Both source-reader and caller failures join every worker before
cleanup. The 128-MiB native index formats retain serial reads. Database readback,
source signatures/body checks and immediate retirement checks are unchanged.

For a controlled source-read latency comparison, add
`--source-read-latency-ms 1`; `--keep-originals` isolates publication from
retirement. Three 500-attempt/five-JSON repetitions on Python 3.13.16/macOS
arm64, compared with 0.8.433, had median complete collection times of
17.191/7.692 seconds and publication-only times of 10.219/1.101 seconds.
Every body digest and structured-record row digest matched. Complete warm
collection medians were 1.023/1.267 seconds: read-ahead adds local scheduling
overhead when storage is already fast. Modified runs used 35.7–39.5 MiB peak
RSS for these small bodies; the read-ahead byte bound is larger for maximum-size
records. The latency is simulated, so these results do not predict cluster
throughput or the duration of a whole migration.

## Captured runtime sources

Completed `workspace/output/runtime` snapshots can be collected explicitly after
all workspace writers stop. This is off by default. The existing metadata store
keeps each captured source body once by SHA-256 in `runtime_bodies`, with exact
path, mode and mtime references in `files`. It fully verifies body references,
source identities and unchanged bytes before unlinking the raw captures.

```bash
python workflow/support/runtime_snapshot_storage.py workspace/output/runtime \
  --collect --writers-stopped --operation-id exact-maintenance-id
```

Only recognized `gg_*_entrypoint/RUN/` captured scripts and `VERSION` are eligible.
Unknown files, coordination trees and the immutable shared source pool remain.
Only empty run directories represented in the verified DB are pruned. Readonly
hardlinks are verified with link counts adjusted for this operation's own
unlinks; external link or content changes stop retirement.

The current entrypoint restores only its archived run, after verifying that
all current source names and bytes match the capture. Other archived runs stay
in SQLite. The existing full metadata restore command restores all captures
before use with an older runtime. An archived capture never selects a scientific
runtime or grants execution authority.

The comparable local benchmark uses 2,000 runs, 24 source versions and four
files per run, with one warmup and three repetitions in Linux aarch64/Python
3.12.15. It retained identical path, byte, mode and mtime digests while reducing
8,096 file paths and 2,002 directories to 98 files and four directories. Median
logical inventory time decreased from 0.195 to 0.0485 seconds. Original files
already shared 96 source inodes, so path-count reduction is not an inode-quota
claim; removing the 2,000 run directories is the main inode saving. The SQLite
DB occupied 6.4 MB. These fixture results do not predict cluster throughput.

```bash
python workflow/benchmarks/benchmark_runtime_snapshot_storage.py \
  --snapshots 2000 --versions 24
```
