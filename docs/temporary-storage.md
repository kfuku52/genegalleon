# Temporary computation storage

`GG_COMMON_TMP_ROOT` selects disposable workflow scratch storage. The default is
`auto`: NIG execution nodes use their independent `/data1` filesystem; all other
environments, including Shirokane and audrey1, use workspace storage. Selection
happens on the execution host before container launch, and the requested and
selected values are logged. NIG detection requires an explicit
`GG_SITE_PROFILE=nig`, or a recognized NIG node name together with the NIG Lustre
home layout. A hostname beginning with `m` alone is insufficient. NIG `auto`
fails if `/data1` is unavailable, unwritable, or on the operating-system
filesystem; it does not silently consume workspace quota or system `/tmp`.
Set `GG_COMMON_TMP_ROOT=workspace` explicitly when workspace storage is intended.

Workspace computation paths remain below `<repository>/workspace/output`; an
explicit `gg_workspace_dir` changes that workspace root. Auxiliary tool
`TMPDIR`, `TMP`, `TEMP` and Python bytecode caches use a private
`.genegalleon-runtime-<uid>` directory in the selected storage, including
workspace mode. External task supervisors further isolate tool temporary files
by run. These auxiliary directories are caches, not scientific checkpoints.

| Value | Location |
| --- | --- |
| `auto` (default) | NIG `/data1`; otherwise workspace |
| `workspace` | Existing workflow-specific output `tmp` directories |
| `/scratch/user` | Private GeneGalleon directories inside that existing directory |
| `env` | Execution host's `TMPDIR`, resolved immediately before container launch |
| `/tmp` | Explicit opt-in; check its filesystem and capacity first |

For example, in a job script running on the compute node:

```bash
GG_COMMON_TMP_ROOT=env bash workflow/gg_gene_evolution_entrypoint.sh
# Or choose an existing writable directory explicitly:
GG_COMMON_TMP_ROOT=/scratch/user bash workflow/gg_genome_evolution_entrypoint.sh
```

`env` fails when `TMPDIR` is unset. An explicit root must be absolute, existing,
writable and searchable. Commas, colons and newlines are rejected because the
container bind syntax cannot represent them safely. Invalid roots fail without
falling back to workspace or `/tmp`. The selected root is mounted at `/gg_tmp`
with Apptainer, Singularity, and the Docker-backed runtime. Custom adapters must
use `exec`; `shell` adapters are rejected for external scratch. Use the entrypoints;
direct core execution does not set up or supervise external scratch.

Each external run uses a private directory scoped by user, workspace, workflow,
gene-evolution mode (where applicable), and execution. The host path, container path and available bytes are logged.
Child tools receive `TMPDIR`, `TMP` and `TEMP` pointing inside that run. Main
computation directories in input generation, genome annotation, transcriptome
generation, gene evolution, genome evolution and fractionation bias are routed
there too, including the derived protein staging directories used by genome
evolution. Tools that insist on writing next to their output may still write
there; this option is not a promise to relocate every temporary byte.

Input-generation task plans, array summary shards, shared download caches,
validation stamps, locks and ZIP materialization receipts remain in the
workspace. Output publication staging also remains beside the destination so
atomic rename works across filesystems. A completed scratch computation is not
a durable checkpoint: only published outputs are authoritative.

Successful runs remove their external run directory unless `delete_tmp_dir=0`.
Failures and interrupted runs retain scratch. A forwarded termination signal
remains a failure even if a child's signal handler returns exit code zero.
Gene evolution also honors `delete_preexisting_tmp_dir`: 1 removes idle scratch
for the same task number, mode and input assignment before a new run; 0 reuses
the most recent matching idle directory. Query mode also requires an unchanged
query filename inventory and sorting locale; orthogroup mode checks the family
assigned to that table row. Active runs are protected with
file locks inherited by the computation, including through the metrics monitor.
Completion cleanup takes the same scope lock as startup and retention scans,
so a finishing job cannot remove a directory another launcher is inspecting.
Termination while waiting for that lock retains scratch and returns failure.
Changing workspace, workflow, scratch root or compute node does not
recover node-local files from a different location.

Gene evolution's existing `gene_family_tmp_retention_days`,
`gene_family_tmp_max_dirs`, `gene_family_tmp_max_bytes` and
`gene_family_tmp_max_files` limits also apply to idle external runs, scoped by
workspace, workflow and gene-evolution mode at that scratch location.
Transcriptome generation has the corresponding `transcriptome_tmp_retention_days`
(default 7 days), `transcriptome_tmp_max_dirs` (100),
`transcriptome_tmp_max_bytes` (1 TiB) and `transcriptome_tmp_max_files` (200,000)
controls for failed or interrupted task directories. A per-task
`.gg_active.lock` prevents cleanup while a retry is running. Set a limit to `0`
to disable that particular bound. Cleanup occurs on invocation and completion,
not as a background service. Other workflows retain failed runs until manually
removed or removed by the site's policy. Never remove an active run.
`delete_tmp_dir=0` prevents successful-run cleanup by GeneGalleon but does not
disable retention limits or site cleanup.

Transcriptome getfastq files can be retained in a reusable cache by setting
`transcriptome_getfastq_cache_dir` (or `GG_TRANSCRIPTOME_GETFASTQ_CACHE_DIR`) to
an absolute directory outside `transcriptome_assembly`. The cache is partitioned
by species and written only after the completion manifest, metadata SHA-256,
getfastq parameters and FASTQ filesystem identity have been checked. A later job
with the same metadata and parameters reuses that cache before invoking
`amalgkit_getfastq`; `run_amalgkit_getfastq=0` fails closed if the contract is
missing or stale. With this cache enabled, `remove_amalgkit_fastq_after_completion=1`
preserves the cached FASTQs. The entrypoint binds the configured host directory
at a stable container path, so it may be outside the selected workspace. It is
intentionally outside the disposable assembly workspace and is not removed by
`delete_tmp_dir`.

The default workspace paths are:

| Workflow | Default computation directory (relative to workspace) |
| --- | --- |
| Input generation | `output/input_generation/tmp` |
| Genome annotation | `output/tmp/<task>_<species>` |
| Transcriptome generation | `output/transcriptome_assembly/tmp/<task>_<species>` |
| Gene evolution | `output/{orthogroup,query2family}/tmp/<task>_<family>` |
| Genome evolution | `output/species_tree/tmp`, `downloads/tmp/species_protein` and its staging variants |
| Fractionation bias | `output/tmp/kffractbias/<task>_<analysis>.<unique>` |

Node-local storage may disappear at job termination even when GeneGalleon keeps
it. Choose workspace storage when retained intermediate files must remain
accessible after the job; choose site-approved scratch explicitly when disposable
local computation is desired. No performance improvement is claimed without a
representative benchmark.

## Retained files and legacy workspace cleanup

Genome evolution stages derived OMAmer query FASTAs and the resolved genetic-code
table inside the job's computation scratch. Successful OMArk output validation
and provenance recording precede query deletion; `delete_tmp_dir=0` keeps query
scratch for debugging. Failed searches retain it. OMArk summaries fingerprint
the result directory, which contains no newly generated query FASTAs.

Older workspaces can contain `output/genome_evolution/omark/<species>/<species>.query.fa`
and `downloads/tmp/species_genetic_code.resolved.tsv`. Inventory these before
cleanup. Verify retained inputs, complete `.omamer` and `.sum` results, and current
provenance, and exclude running or queued consumers. Deleting a legacy query
changes the OMArk directory fingerprint: preserve the old manifest and regenerate
the summary contract from unchanged `.sum` files instead of ignoring a stale
contract or rerunning OMAmer unnecessarily.

Download `*.corrupt.*` quarantines can be removed after a good replacement is
validated and diagnosis is finished. `.part` files remain resumable downloads.
Provider `.archive_cache` ZIPs, configured getfastq caches, and `orthofinder/core`
results are retained for reuse and require an explicit project-specific decision;
ordinary scratch cleanup does not remove them. Input-generation task plans,
completion receipts, and MCMCtree chain evidence must remain available.

For an existing project's retired failed scratch, first list exact paths, owners,
file counts, bytes and consumers. Reuse the existing gene-family `cleanup-tmp`
and transcriptome retention guards where applicable. For other retired scratch,
take and checksum-verify a recovery archive, recheck ownership, path identity and
job inactivity, then remove only the inventoried paths. Do not apply wildcard
deletion to `workspace/output`, `workspace/downloads`, or another user's data.

Generate a read-only JSON inventory without launching a workflow:

```bash
python workflow/support/workspace_cleanup_inventory.py --workspace /absolute/project/workspace
```

For an older `gfe_data` layout, pass that directory with `--legacy-output-root`.
The report distinguishes scratch, quarantines and retained results/caches, records
foreign-owned entries and inaccessible paths, and never follows symlinks. It does
not establish job inactivity or authorize deletion. Nested cache/quarantine records
can overlap; do not sum their counts as independent storage. Configured external scratch
and FASTQ caches must be inventoried separately at their recorded host paths.
