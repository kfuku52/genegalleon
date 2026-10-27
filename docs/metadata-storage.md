# Offline metadata collection

`workflow/support/metadata_store.py` collects completed observation JSON,
resource metrics, artifact provenance, version reports, version-inventory cache,
output-inventory journals and native ZIP index buckets into `.gg-metadata.sqlite3`
inside each existing metadata root. Scientific annotation/statistics TSV bodies,
FASTA, trees and ZIP members keep their existing formats. Family state, archive
queues and distributed locks remain under their existing ownership protocols.

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
Database identity fences invalidate observations spanning a replacement.

Restore before running an older reader. Restoration retains the DB and refuses
to overwrite changed raw files; native provenance tombstones and newer ZIP
members remain authoritative.

```bash
python workflow/support/metadata_store.py restore \
  --root workspace/output/observations --writers-stopped
```

For `version_inventory` only, a separate `pack --cache-limit N` call after
collection expires regenerable contexts, retaining the newest N (1–1024)
64-hex log/JSON pairs. Unknown names remain. Once activated, the producer uses
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
