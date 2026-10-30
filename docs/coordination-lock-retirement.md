# Offline coordination guard retirement

`workflow/support/coordination_lock_storage.py` provides `snapshot`, `retire`
and `verify_retired` for exact idle per-family store guards and run-lease guards.
The caller must fence every old producer start and verify that no project job
is queued or running before passing `writers_stopped=True`. An empty marker,
its age, or an unavailable process is never sufficient authorization.

Retirement additionally holds the existing archive-reader, archive-maintenance
and exclusive family gates. Any primary run lease, occupied namespace, changed
identity, unknown recorded shape or pending storage conversion stops it.
Scientific files, publication guards and stable global gates stay unchanged.

Seal the physical output root and `snapshot(root)` before execution. Supply an
owned private journal directory. A verified, fsynced gzip journal records every
original empty marker and namespace identity before deletion. Keep this journal
and reuse the identical snapshot to resume a partially removed operation; never
reseal a replacement marker. `verify_retired` checks that each recorded guard is
absent. These empty synchronization objects contain no scientific results.

Existing runtimes recreate the same per-family paths when needed. This offline
operation reduces idle metadata at a maintenance boundary; subsequent work can
create guards again. It does not change concurrency, lease ownership, archive
formats, annotations or statistical TSVs.

A Docker regression benchmark uses 1,000 synthetic families, one 32-family
warmup and three identical trials. `workflow/benchmarks/benchmark_coordination_lock_storage.py`
counts the complete temporary tree, including retained bucket directories, newly
created global gates and the journal: 6,259 entries before, 269 afterward;
6,000 guard entries removed. Protected scientific fixture SHA-256 matched in
all trials. On the recorded development image, retirement wall times were
0.711, 0.828 and 0.858 seconds (median 0.828); process peak RSS was 41,404 KiB.
These are local synthetic Docker measurements, not Shirokane timing or SIF proof.
