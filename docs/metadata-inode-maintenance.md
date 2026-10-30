# Metadata inode maintenance

Default performance traces use 64 append stripes under
`workspace/output/versions/performance`, protected by fixed namespace locks.
The JSONL record schema is unchanged. Offline `metadata_store.py pack`
losslessly moves legacy PID files and stripes into that root's immutable SQLite
DB. A new collection generation gives subsequent stripes new names, so a later
collection cannot overwrite historical records. `performance_metrics.py read
--directory DIR` streams archived and live records. Custom `GG_PERFORMANCE_DIR`
locations retain their historical per-process layout. Annotation/statistics TSV
files are excluded.

With all project writers stopped, the owning maintenance CLI can retire only
empty directories whose archived observation says the process exited at least
24 hours ago, or unused legacy version-cache lock markers/namespaces. Observation
completion validation stays independent: an exited run is not a successful
scientific result. Unknown, recent, nonempty or active entries remain.

```bash
python workflow/support/metadata_maintenance.py prune-observations \
  workspace/output/observations --journal /private/recovery/observations \
  --operation-id exact-offline-operation --writers-stopped
python workflow/support/metadata_maintenance.py retire-version-locks \
  workspace/downloads/version_inventory --journal /private/recovery/cache-locks \
  --operation-id exact-offline-cache-operation --writers-stopped
python workflow/support/metadata_maintenance.py verify \
  workspace/output/observations --journal /private/recovery/observations
```

The caller must independently establish offline project/queue/claim ownership;
`--writers-stopped` records that prerequisite, not proof by itself. Upgrade the
project's exact default workflow API before removing observation directories.
The new API reconstructs archived attempt directories from the verified SQLite
DB; missing unarchived paths, corrupt data and permission errors still fail.
Cache lock retirement requires the existing SQLite producer's common lock.
The DB and its common lock/metadata guard always remain.

Each private recovery journal binds original identities and DB SHA before the
first removal. Only the original operation can resume missing members. Changed
inodes, databases or recovery records stop execution. `restore ROOT --journal
JOURNAL --writers-stopped` recreates the selected empty directories/markers with
original feasible metadata, and refuses new contents. It does not modify the
DB. Restore physical observation directories before downgrading the reader;
restore archived file bodies with the existing metadata collector's restore CLI
before a program downgrade that cannot read SQLite. Keep all recovery journals.

## Repeated fixture measurement

Run in the ordinary fresh runtime:

```bash
bash workflow/tests/run_in_runtime.sh python \
  workflow/tests/benchmark_metadata_inodes.py --count 2048 --trials 3
```

Docker aarch64, Python 3.12.15: all three trials retained identical logical bodies
and verified restoration. Counts include target-root descendants and permanent
coordination files; the two external recovery journals add six retained files.

| Target, 2,048 records | Before | After |
| --- | ---: | ---: |
| Archived observation directories | 2,052 | 4 |
| Legacy version-cache controls | 6,148 | 7 |
| Legacy PID traces after collection | 2,048 | 4 |
| New live performance traces, including guards | 2,048 | 256 |
| New striped traces after collection | 256 | 196 |

Median append time was 0.041 s for historical plain appends and 0.743 s for the
same 2,048 records with stripe/collector coordination (about 0.36 ms per record).
This trades extra telemetry coordination for bounded live inodes. Median offline
collection/retirement time for the fixture was 2.036 s. These are local Docker
fixture measurements, not shared-filesystem throughput claims. Existing legacy
controls are retired only once; fixed current controls are deliberately retained.
