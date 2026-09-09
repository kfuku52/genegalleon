# Input generation with species arrays

`array_prepare` freezes a species task plan and prepares shared taxonomy/BUSCO
resources. With a download manifest, each `array_worker` downloads its own
species and runs formatting, validation, fx2tab, and BUSCO. `array_finalize`
requires verified completion receipts for every planned species before publishing
the merged species summary and resolved download manifest.

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
It also works after the helper was interrupted following successful prepare.
An atomic prepare-completion marker prevents retry from bypassing failed or
incomplete shared setup; rerun without `--retry` in that case. If all workers
are complete, it submits just finalize. Old failed-array dependent finalizers may
remain pending in Slurm; cancel those by their printed IDs when superseding them.
UGE/PBS users can continue invoking the three existing modes directly, with one
worker per 1-based task index. Automated submission in this helper is Slurm only.

## Frozen inputs and restart behavior

Every selected manifest row must have an explicit, valid `species_key`, a
supported `provider`, and an `id`. Duplicate output species prefixes are rejected,
even across providers. Species keys and explicit download filenames must be
non-hidden filename components, without directory separators or control characters. Prepare embeds the selected rows; workers do not reread a
mutable manifest. Local source references are resolved before embedding. Local
raw inputs, including file URLs in manifests, are hashed during planning; worker-resolved downloads are hashed before
formatting. Changed raw inputs are rejected instead of silently reusing a plan.

Plans and execution settings are immutable. A workspace and its custom output
directories are bound to one plan; a different plan cannot reuse their shard
namespace or write concurrently through another workspace. Output-directory
locks and ownership records live beside those directories. Repeating the same prepare preserves
completed work; changed inputs or settings require a new output workspace. Do not
edit plan/settings/receipt files or remove lock files while jobs are active.
Completion receipts are written atomically only after all enabled worker stages
succeed and include hashes of raw inputs, formatted outputs, enabled fx2tab/BUSCO
outputs, and summary/statistics shards. Finalize checks exact shard indices,
species identities, receipts, and outputs; incomplete or stale results leave the
canonical species summary intact. The original manifest is not reread during
workers or finalization; trait species are reconstructed from the frozen plan.
Trait configuration files and the resolved shared BUSCO lineage are also
checked for changes. Canonical tables are staged beside their destinations and
renamed only after the optional final shared stages succeed. Each table rename
is atomic; publication of multiple files is not a filesystem-wide transaction. Shared stages and workers also hold workspace
locks to prevent simultaneous publication or cleanup.

Array mode retains `tmp/task_plan.json`, settings, worker downloads, and receipts
for auditing/retry. Storage can be reclaimed after the run is no longer needed,
with no jobs active. Shared lock/ownership sidecars must also be preserved while
a plan remains in use. Retrying after deleting raw downloads requires a new plan;
those raw files are part of the completion evidence. For clusters without compute
node internet access, download in a network-enabled job first, then run arrays
with `input_dir` over those local inputs and no download manifest.

## Shared database request limits

All input-generation jobs in the same workspace default to the shared directory
`workspace/.gg_cache/input_download_limits`. To coordinate different workspaces,
set the same `GG_INPUT_DOWNLOAD_LIMIT_DIR` in all jobs. This path must be visible
at the same absolute path inside their containers, on a shared filesystem with
cross-node `flock` support. Node-local `/tmp` is unsuitable.

```bash
export GG_INPUT_DOWNLOAD_LIMIT_DIR=/shared/project/download_limits
export GG_INPUT_MAX_CONCURRENT_DOWNLOADS_NCBI=2
export GG_INPUT_REQUEST_INTERVAL_NCBI=0.4
```

The concurrency limit covers request opening and streaming until response close.
Start intervals are in seconds, independent of the number of worker jobs or CPUs.
The default is two requests per logical database, with 0.4 seconds between starts.
NCBI API/FTP/www domains share `NCBI`; Ensembl/EnsemblGenomes share `ENSEMBL`.
CoGe, CNGB, GWH and DDBJ also have domain groups. Other supported providers use
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

Slots use kernel-held locks, released on process exit, including SIGKILL. There
is no age-based eviction that could admit extra transfers while an old owner is
alive. Lock files persist and must not be deleted while clients exist. Request
interval/cooldown timestamps require synchronized node clocks. Docker tests cover
multiple processes and SIGKILL; cross-node lock behavior on the target HPC mount
and SIF must be validated before relying on these limits there. A filesystem that
does not honor remote locks cannot be made safe by this implementation alone.

See the [implementation review and validation limits](input-generation-array-review.md).
