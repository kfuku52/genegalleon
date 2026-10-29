# Closed scheduler-log storage

Scheduler stdout/stderr are not gene-family artifacts or metadata JSON. Route
new submissions to a dedicated batch directory, then pack that directory using
the existing verified ZIP64 directory store. No scientific outputs, exit status,
or scheduler execution records are rewritten. Existing scheduler log locations
are unchanged unless the submitter explicitly selects the new layout.

```text
workspace/output/scheduler_logs/batch_001/
  scheduler_output/       # scheduler writes here while the batch runs
  scheduler_output.zip   # replaces the directory after verified collection
```

Create `scheduler_output/` before submission and pass its absolute path to the
scheduler's stdout/stderr options (UGE `-o`/`-e`, or SLURM `--output`/`--error`
with job/task placeholders). Use a fresh batch directory for every submission,
including retries. Do not append a new job to an archived batch.

The collector must run **after all batch tasks and scheduler output delivery have
finished**. A worker's EXIT trap is too early: the scheduler can still write its
final output. Scheduler dependencies can arrange one finite collector per
submission batch. Verify terminal state/output delivery through the site's normal
execution controller before passing `--writers-stopped`; this acknowledgement is
not a scheduler check, and file age is not evidence that a job stopped. Failed
jobs' closed logs are eligible too; packing never marks science successful.

```bash
# Run inside the GeneGalleon runtime after the batch is closed.
python workflow/support/directory_output_store.py pack \
  --kind scheduler_logs --root /path/to/batch_001 --writers-stopped
python workflow/support/directory_output_store.py verify \
  --kind scheduler_logs --root /path/to/batch_001

# Read one log without expanding all logs. stdout is the original binary body.
python workflow/support/directory_output_store.py read \
  --kind scheduler_logs --root /path/to/batch_001 \
  --relative-path scheduler_output/gene_job.e123.7

# Restore historical physical paths when an older reader requires them.
python workflow/support/directory_output_store.py materialize \
  --kind scheduler_logs --root /path/to/batch_001 --writers-stopped
```

ZIP members preserve bytes, SHA-256, file modes and nanosecond modification times.
Large logs stream with bounded buffers instead of the metadata JSON size limit.
Publication verifies the new archive and rechecks sources before removing them;
changed sources, symlinks, insufficient space and mixed conflicting forms fail
without discarding the recovery copy. An interrupted retirement can resume
against the same archive. Ordinary ZIP tools can inspect the archive; do not
edit it in place. Existing directory-store readers understand logical log paths.

The lock covers one physical batch root. Collecting batch_001 does not acquire a
workspace-wide or project-wide fence and does not prevent batch_002 from running.
The scheduler itself does not participate in these locks; the closed-writer
precondition is essential. Shared lock namespaces are retained under their
existing ownership protocol.

For an entirely idle workspace, `gg_workspace_storage.sh` also discovers
`output/scheduler_logs/<batch>/scheduler_output` and its ZIP:

```bash
bash workflow/gg_workspace_storage.sh convert --workspace /path/to/workspace \
  --target scheduler_logs --to zip --writers-stopped --verification deep
```

Do not use that all-batch command while another discovered batch is active.
During an array, collect only the exact closed batch via `directory_output_store`.
Use bounded submission batches if the live-file peak matters. There is no
background collector and no compression in worker exit traps.

Closed execution JSON, resource metrics, version reports, artifact provenance
and family inventories use [offline metadata SQLite collection](metadata-storage.md).
Their writers must stop before collection. Keep state, queues, owner records,
claims and lock namespaces under their original protocols. ZIP/SQLite collection
reduces payload files, not every persistent directory or coordination inode.

## Reproducible file-count check

```bash
python workflow/benchmarks/benchmark_scheduler_log_storage.py --tasks 1000 --repeats 3
```

A synthetic Linux arm64/Docker fixture with Python 3.12.15 used two 4-KiB logs
per task. After one discarded warmup, three fresh-process trials packed 2,000
logs into one ZIP. Including the retained lock guard, regular files changed from
2,000 to 2 and directories from 2 to 3. Every body hash, size, mode and mtime
matched. Packing took 0.307/0.291/0.289 seconds, with peak process RSS
39.4–39.5 MB. These are local fixture costs and file-count reductions, not a
SHIROKANE throughput measurement or a prediction of GeneRax computation time.
