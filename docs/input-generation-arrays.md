# Input generation with species arrays

`array_prepare` is one download/prepare job: it freezes a species task plan,
downloads manifest inputs with database-specific parallel queues, hashes the
local files, and prepares shared taxonomy/BUSCO resources. Existing gzip
downloads are checked through fail-closed validation receipts under the shared
staged-cache parent; a matching source identity and file inode/size/mtime
allows a retry to skip the full gzip read. A missing, malformed, or stale
receipt always falls back to full validation. Each `array_worker`
uses those staged files for formatting, validation, fx2tab, and BUSCO; it does
not fetch missing reference files. `array_finalize`
requires verified completion receipts for every planned species before publishing
the merged species summary and resolved download manifest.

For existing raw files, a manifest row can opt into `bind_local_sources=1`.
Every supplied role must then use an absolute `file://` URL without an archive
member. Prepare validates and binds the frozen source hashes directly; it does
not create another raw-file copy. Workers and subsequent prepare checks reject
changed/missing sources. CoGe GFF content checks still apply. Default rows retain
the ordinary isolated staged-copy behavior. Keep bound files visible in the
same container namespace until all workers and finalization finish.

SHA-256 checks read each unique path once within a verification phase, including
roles sharing a file and original/staged references to bound sources. Preflight,
binding, and publication remain independent full-content checks; later invocations
never trust a persistent size/mtime hash cache. Worker preflight checks both the
original plan and resolved receipt in one pass, rejecting conflicting hashes.
Hashing checks the open file and current path identity to reject modifications or
replacement during the read. Plan and receipt formats are unchanged.
At worker completion, raw sources and declared outputs share one fresh hashing
pass. Canonical paths are read once even when a raw source is also an output or
has a symlink alias; each original receipt key is retained, and alias replacement
or a change during the batch rejects publication. An annotation with no nuclear
CDS records fails formatting before genome output or BUSCO. GFF-derived CDS
also detects an empty coding-feature set before loading the genome into memory.
Repeated CoGe internal feature IDs with an identical source transcript name are
collapsed only after their complete CDS models agree on coordinates, strand,
phase and annotation attributes. Both CDS derivation and GFF repair apply this
proof; conflicting or partial models still fail.

For a new workspace/runtime, bound rows can reuse successful native staging
receipts from a previous prepare, including a partially failed prepare. Supply
all five manifest columns: `reuse_staged_plan` (absolute old plan path),
`reuse_staged_plan_sha256`, `reuse_staged_workspace` (absolute old workspace),
`reuse_staged_task_index` (one-based), and `reuse_staged_receipt_sha256` (SHA-256
of the old `task_plan.json.tasks/N.json`). The old staged version-2 plan and
receipt must agree on plan hash, index, provider, species, every role and its
source hash. Old `/workspace/...` receipt paths are mapped through the supplied
workspace; other paths must match exactly. Incomplete or inconsistent proof
is rejected, rather than silently using it.

Planning freezes those previously validated hashes without rereading the raw
files. Prepare then performs one fresh full SHA-256 read per unique source and
reuses gzip validation only after matching the sealed receipt. That read is
fenced through binding, publication and the end of staging by device, inode,
size, mtime and ctime checks. A change, replacement or symlink invalidates the
attempt, including changes after an earlier species was published. No checksum
is cached across invocations; a retry and every worker still verify content
afresh. CoGe GFF checks remain current. This reuses raw staging validation only;
formatting, CDS/GFF validation and BUSCO results are not certified by it.

BUSCO lineage downloads are extracted into a temporary directory inside the
workspace download cache and published only after BUSCO succeeds. An incomplete
download remains there for diagnosis and is not treated as a ready lineage.
ETE taxonomy database intermediates are also built in the workspace taxonomy
directory, so container-home overlay capacity does not limit preparation.

```mermaid
flowchart LR
  D[Single download/prepare job<br/>Independent database queues] --> L[Hashed local inputs]
  L --> A[Species compute array<br/>Multiple nodes]
  A --> F[Single finalize job]
```

This applies to `gg_input_generation`. RNA `amalgkit getfastq` belongs to
`gg_transcriptome_generation` and is not changed by this execution model.

## Slurm submission

Use the usual `GG_INPUT_*` overrides and shared workspace configuration. Run the
helper from the same directory as the normal entrypoint. For example:

```bash
export GG_INPUT_PROVIDER=all
export GG_INPUT_DOWNLOAD_MANIFEST=/shared/project/download_plan.tsv
python workflow/gg_input_generation_array.py \
  --task-plan /shared/project/workspace/output/input_generation/tmp/task_plan.json \
  --partition compute --max-running 8 --cpus 4 --memory 32G
```

This is a **dry run**: it prints commands without submitting or downloading.
When the plan already exists, it also displays exact species/task counts. Before
prepare has run, the count is shown as `N`, to be read from the resulting plan.
`--memory` is total memory per task. The helper uses `sbatch --wrap`, so resource
settings come from these options rather than scheduler headers in the entrypoint.
Set the workspace through the project's common configuration as usual; the plan
path alone does not change the workspace.

Download workers and compute resources are independent. Use `--prepare-cpus`
to set the download/prepare job's CPUs and internal worker count,
`--prepare-memory` for its total memory, and `--prepare-partition` for a
network-enabled partition when needed. Omitted values inherit `--cpus`,
`--memory`, and `--partition`. For example, `--prepare-cpus 8 --prepare-memory 8G`
can be combined with `--cpus 4 --memory 32G --max-running 8`. These are resource
requests, not database connection limits; choose them from measured usage.

Add `--submit` to run prepare with `sbatch --wait`, then submit the species array
and a finalizer with `afterok` on that array. Keep the helper running until those
IDs have been printed. If prepare fails, no workers are submitted. Scheduler
rejection details and any returned job ID are preserved in the error report.
Contiguous task IDs are compressed into Slurm ranges; the cluster
`MaxArraySize` still limits the largest index. If a worker
fails, its dependent finalizer will not run. The printed job IDs let you inspect
or cancel those jobs with normal scheduler commands.

```bash
python workflow/gg_input_generation_array.py \
  --task-plan /shared/project/workspace/output/input_generation/tmp/task_plan.json \
  --partition compute --max-running 8 --cpus 4 --memory 32G --retry --submit
```

`--retry` skips prepare and selects only missing or invalid completion receipts.
For a submitted retry with missing workers, the native helper refuses to run
while any `gg_input_array_worker` job owned by the caller is active or pending.
This is conservative across projects because older Slurm jobs do not record
their task-plan identity. Previewing `--retry` remains read-only.
It also works after the helper was interrupted following successful prepare.
An atomic prepare-completion marker prevents retry from bypassing failed or
incomplete shared setup; rerun without `--retry` in that case. If all workers
are complete, it submits just finalize. Old failed-array dependent finalizers may
remain pending in Slurm; cancel those by their printed IDs when superseding them.
UGE/PBS users can continue invoking the three existing modes directly, with one
worker per 1-based task index. Automated submission in this helper is Slurm only.

## Frozen inputs and restart behavior

`GG_INPUT_REQUIRE_GENOME=1` opts into requiring a nonempty formatted genome
FASTA for every selected species. The default is `0`, so projects using only
CDS/annotation remain supported. In array mode, prepare rejects a staged
species without a genome FASTA or genomic GBFF source before workers launch;
workers and finalize verify the formatted output. The setting is frozen with
the array plan. Enable it only in a fresh workspace/plan, rather than changing
settings on an active or previously prepared array.
The output check reads each required FASTA or GFF through the end, including
the gzip trailer. It rejects empty records, unsupported nucleotide symbols,
and records made only of gap/missing symbols. Required GFF features must have
nine populated fields, positive ordered coordinates, a finite score or `.`,
and valid strand and phase values. These structural checks do not establish
CDS/GFF identity or scientific annotation quality. The check does not determine
whether a genome is nuclear; nuclear-only datasets
must review assembly provenance and exclude organelle-only sources in their
manifest.

`GG_INPUT_REQUIRE_CDS=1` and `GG_INPUT_REQUIRE_GFF=1` independently require
formatted CDS FASTA and GFF with at least one feature. Both default to `0`.
CDS derived from GFF plus genome or GBFF counts; GFF derived from GBFF counts.
`gg_input_generation` already needs a usable source of CDS for its species tasks;
leaving the CDS flag off does not make it a CDS-free workflow. The flag adds
an explicit check that formatting produced a FASTA sequence for every species.
The GFF requirement rejects a species without a GFF or GBFF source during
prepare. Workers check required outputs before issuing completion receipts,
and finalize checks every species again. These settings are frozen with the
array plan. An intentionally CDS-only run also needs `GG_INPUT_RUN_VALIDATE_INPUTS=0`,
because CDS-to-GFF mapping validation requires GFF regardless of these flags.

`GG_INPUT_BUSCO_TIMEOUT_SECONDS` defaults to `0` (no timeout). A positive
value bounds each species BUSCO invocation and is frozen into the array plan.
If the limit is reached, that worker fails without a completion receipt so a
later retry can select it; the timeout does not turn missing BUSCO output into
a successful species. Worker logs keep the first 10,000 BUSCO lines and the
last 50 lines, while preserving the BUSCO exit status. This avoids unbounded
logs from a repeatedly failing predictor. Choose the limit after measuring
large representative species, not from the scheduler wall time alone.

Every selected manifest row must have an explicit, valid `species_key`, a
supported `provider`, and an `id`. Duplicate output species prefixes are rejected,
even across providers. Species keys and explicit download filenames must be
non-hidden filename components, without directory separators or control characters. Prepare embeds the selected rows; workers do not reread a
mutable manifest. Local source references are resolved before embedding. Local
raw inputs, including file URLs in manifests, are hashed during planning
(or frozen from sealed staging evidence and freshly hashed by prepare);
downloaded inputs are hashed by prepare. Resolved tasks and manifests are bound
to the prepare-completion marker. Changed raw inputs or missing staged receipts
are rejected instead of silently downloading during a worker run.

Plans and execution settings are immutable. A workspace and its custom output
directories are bound to one plan; a different plan cannot reuse their shard
namespace or write concurrently through another workspace. Output-directory
locks and ownership records live beside those directories. Repeating the same prepare preserves
completed work; changed inputs or settings require a new output workspace. Do not
edit plan/settings/receipt files or remove lock files while jobs are active.
Completion receipts are written atomically only after all enabled worker stages
succeed and include hashes of raw inputs, formatted outputs, enabled fx2tab/BUSCO
outputs, and summary/statistics shards. When CDS/GFF validation runs, its per-task
`tmp/task_stats_shards/N.mapping.json` records phase and UTR conflict counts and
is included in the completion receipt; finalization aggregates available QC
shards without requiring them from workers completed by older runtimes.
Finalize also publishes `species_mapping_qc.tsv` with one row per species;
`not_recorded` means mapping QC is absent (for example, an older completed
worker or a run without validation), not that its annotation was clean.
Finalize checks exact shard indices,
species identities, receipts, and outputs; incomplete or stale results leave the
canonical species summary intact. The original manifest is not reread during
workers or finalization; trait species are reconstructed from the frozen plan.
Trait configuration files and the resolved shared BUSCO lineage are also
checked for changes. Canonical tables are staged beside their destinations and
renamed only after the optional final shared stages succeed. Each table rename
is atomic; publication of multiple files is not a filesystem-wide transaction. Shared stages and workers also hold workspace
locks to prevent simultaneous publication or cleanup.

### Resume completed worker stages

With `overwrite=0`, workers retain independent successful formatting and
CDS/GFF-validation checkpoints under `tmp/stage_checkpoints/`. A BUSCO failure
does not invalidate those upstream stages. On retry, fresh content hashes of
the declared inputs, outputs, QC and summary/statistics shards must match the
checkpoint, along with formatting/validation parameters. A matching checkpoint
skips the formatter or validators entirely; output existence alone never
certifies success. fx2tab and BUSCO use their existing artifact contracts.
`overwrite=1` explicitly reruns the enabled stages.

For a changed BUSCO lineage, prepare a separate output workspace as usual and
set these entrypoint parameters (or their `GG_INPUT_*` environment overrides):

- `resume_from_task_plan`: the absolute frozen donor plan path;
- `resume_from_task_plan_sha256`: the SHA-256 of that exact plan;
- `resume_from_input_generation_root`: its `output/input_generation` directory.

Native prepare imports verified formatting, validation and available fx2tab
outputs for matching species and identical raw content. Formatting parameters
and required-output settings must agree. It remaps shard indices and output
paths, writes new stage checkpoints and preserves the donor plan/workspace.
BUSCO results are not imported: the new lineage is assessed normally and
produces new BUSCO provenance. The donor workspace must be inactive; its phase
lock prevents copying during workers, prepare or finalize.

Completed older workers can be imported when their native completion receipt
and format provenance certify the needed files and parameters. A failed older
worker without independent stage checkpoints cannot prove validation succeeded;
its uncertified stages run once. Missing mapping QC cannot be treated as clean
annotation. Subsequent retries retain the new checkpoints.

Array mode retains `tmp/task_plan.json`, settings, staged downloads, and receipts
for auditing/retry. Storage can be reclaimed after the run is no longer needed,
with no jobs active. Shared lock/ownership sidecars must also be preserved while
a plan remains in use. Retrying after deleting raw downloads requires a new plan;
those raw files are part of the completion evidence. A failed prepare preserves
task receipts for species whose complete source bundles were staged
successfully. Download validation or merge errors that cannot be attributed
to an individual species prevent new task receipts for that attempt; discovery
errors prevent new receipts for the affected provider. Downloaded files remain
available for retry. The next prepare retries only unresolved species and rewrites the
pending manifest; a successful prepare rerun verifies staged files without
contacting their original servers. The validation receipt directory is
plan-independent so hardlinked staging directories can reuse it. Use a fresh
workspace if previously frozen files have changed. For clusters without
compute-node internet access, select a
network-enabled `--prepare-partition`; workers use local references and the
shared taxonomy/BUSCO resources prepared there.

The low-level plan script supports on-demand manifest tasks for direct callers.
The core always requests `--stage-downloads`; its staged plans cannot fall back
to worker-side downloading. Do not switch runtimes for an active plan.

On audrey1, a launcher that sets `GG_CONTAINER_PROJECT_ROOT_BIND` to the same
absolute project directory on both sides of the bind runs Apptainer with
`--contain`. This keeps Python/BUSCO semaphores in a private `/dev/shm` while
preserving project-local absolute paths. The project root must contain all
absolute inputs and shared resources needed inside the container; existing
launches without this bind retain their previous runtime behavior.

## Shared database request limits

All input-generation jobs in the same workspace default to the shared directory
`workspace/.gg_cache/input_download_limits`. To coordinate different workspaces,
set the same `GG_INPUT_DOWNLOAD_LIMIT_DIR` in all jobs. This path must be visible
at the same absolute path inside their containers, on a shared filesystem with
atomic `mkdir` and exclusive file creation. Node-local `/tmp` is unsuitable.

```bash
export GG_INPUT_DOWNLOAD_LIMIT_DIR=/shared/project/download_limits
export GG_INPUT_MAX_CONCURRENT_DOWNLOADS_NCBI=2
export GG_INPUT_REQUEST_INTERVAL_NCBI=0.4
```

The concurrency limit covers request opening and streaming until response close.
Start intervals are in seconds, independent of the number of worker jobs or CPUs.
The default is two requests per logical database, with 0.4 seconds between starts.
NCBI API/FTP/www domains share `NCBI`; Ensembl/EnsemblGenomes share `ENSEMBL`.
CoGe, CNGB, GWH, DDBJ, Figshare, FlyBase and WormBase also have domain groups.
EnsemblGenomes' `ensemblgenomes.ebi.ac.uk` endpoints share `ENSEMBL`;
unrelated EBI services do not. Other supported providers use
their provider name for `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_<PROVIDER>` and
`GG_INPUT_REQUEST_INTERVAL_<PROVIDER>`. A recognized destination database takes
precedence over a provider hint. Unrecognized CDN destinations inherit the
originating logical database, including its 429/503 cooldown. Standard Authorization/Cookie headers are stripped on cross-origin redirects.
Provider hints are isolated per request thread; Ensembl variants and RefSeq/GenBank share their
parent database group. Direct/local requests to unknown hosts use independent
host buckets with the `DIRECT` limits. These limits cover the Python
input-provider request path; shared BUSCO/taxonomy setup uses its existing
download locks during prepare.

HTTP 429/503 responses share `Retry-After` cooldowns (seconds or HTTP dates), with
up to three retries. Missing/invalid Retry-After uses 1, 2, then 4 seconds.
`GG_INPUT_DOWNLOAD_LIMIT_WAIT` controls admission timeout, including metadata-lock contention (default 3600
seconds).
Different policies in the same bucket are rejected. To change a policy, stop all
clients before selecting a fresh shared limit directory for all of them.

The dispatcher groups transfers by destination database, using the provider only
as a fallback for unknown hosts. Different `direct` hosts have independent
queues; local copies do not consume network slots. It rotates runnable queues
and checks shared slot/cooldown readiness before assigning workers. The check is
advisory: another downloader can win a slot in between, so the transport still
performs authoritative admission on each request and redirect.

Slots, input-generation output locks, and version-report locks use atomic shared
namespace operations, including on Lustre mounts with `localflock`. There is no
age-based eviction that could admit extra work while an owner is alive. Normal
completion and ordinary command failures release ownership. Signals, signal-like
exit statuses (128 and above), or node failure leave fail-closed owner records,
because child processes may still be writing. Inspect records with
`shared_namespace_lock.py inspect PATH`, stop all
clients, verify the owning job and its children have ended, and reconcile the
specific ownership record before resuming. Never remove locks merely because a
timeout expired. Request interval/cooldown timestamps require synchronized clocks.

The admission protocol uses a `namespace-v1` subdirectory. **Stop all old
flock-based clients before migration**: old and new protocols do not coordinate,
even if configured with the same parent limit directory. Use a fresh output
workspace for the new runtime and preserve old plans/receipts for audit. Validate
normal admission, release, crash behavior, and output locking across the target
nodes before increasing compute concurrency.

See the [implementation review and validation limits](input-generation-array-review.md).

Shared worker lock acquisition retries transient reader-registration contention for up to 30 seconds; exclusive phase and duplicate-worker locks remain nonblocking. A timeout leaves existing ownership untouched.


Source CDS overlaps carrying a consistent `low-quality sequence region` note
can receive the same complete genome-to-publisher-CDS proof as unannotated
overlaps; their source quality note remains unchanged. Unproved, mixed, and
other biological exceptions retain their strict existing checks. Pseudogene
CDS length validation accepts only the formatter's exact terminal `N` padding
to the next multiple of three. Genomic feature lengths stay unchanged, and
pseudogenes do not acquire a coding phase or an inferred intron model.

CDS-only species can run validation and BUSCO with `require_cds=1` and
`require_gff=0`. Workers still validate longest-CDS selection; CDS-to-GFF
mapping QC is produced and receipt-bound only when that task supplies a GFF.
A supplied or required GFF that is missing still fails validation.
