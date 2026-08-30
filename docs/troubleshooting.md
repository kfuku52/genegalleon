# Troubleshooting

## Start with these diagnostics

The quickest first-pass checks are:

```bash
cd workflow
ls -1 ../workspace/output/versions/*.log
```

```bash
cd workflow
test -f ../genegalleon.sif && echo "SIF found"
```

```bash
python -m pytest -q workflow/tests/test_gg_util_paths.py
```

## Common problems

### Container runtime not found

Symptom:

- the wrapper reports that neither `apptainer` nor `singularity` was found on `PATH`.

What to check:

- confirm the compute-node environment exposes the runtime; some sites do not
  make it available on login nodes,
- inspect the startup `site profile` line and set `GG_SITE_PROFILE` explicitly
  if automatic detection selected the wrong profile,
- if your site uses a module system, load the container runtime before launching,
- on NIG, the site profile searches versioned Apptainer/Singularity packages
  under `/opt/pkg` and then the legacy package path,
- if you are on macOS, build and validate with the Docker-backed GeneGalleon
  runtime instead of treating host-local results as authoritative.

See [Site Runtime Profiles](site-runtime-profiles.md) for the recognized
environments.

### `genegalleon.sif` missing

Symptom:

- the wrapper starts normally but fails when entering the container because the default image path does not exist.

What to do:

- build or pull the image with `bash ./gg_container_build_entrypoint.sh`,
- or export `gg_container_image_path=/path/to/genegalleon.sif` before launching the wrapper.
- or `docker pull ghcr.io/kfuku52/genegalleon:latest` and relaunch; wrappers now auto-fallback to the pulled image when repo-root `genegalleon.sif` is missing.
- or switch the wrapper to Docker mode explicitly with `GG_CONTAINER_RUNTIME=docker` and `GG_CONTAINER_DOCKER_IMAGE=<image:tag>`.

### Wrong workspace is being used

Symptom:

- expected inputs are "missing" even though files exist somewhere else on disk.

What to check:

- the wrapper prints the resolved `gg_workspace_dir` during startup,
- the runtime binds that resolved path to `/workspace`,
- if you use an external workspace, export `gg_workspace_dir=/path/to/workspace` explicitly.

### Array task exits immediately

Symptom:

- only task 1 runs, or array tasks fail because inputs are not found.

What to check:

- array size must match the actual input count,
- `gg_gene_evolution_entrypoint.sh` expects one task per query file or per selected orthogroup row,
- `gg_genome_annotation_entrypoint.sh` expects one task per species CDS file,
- `gg_transcriptome_generation_entrypoint.sh` expects one task per chosen input unit,
- the runtime summary prints the normalized `GG_ARRAY_TASK_ID` and `GG_TASK_CPUS`.

### `sbatch` rejects mutually exclusive memory options

Symptom:

- `sbatch` fails before the job starts with an error such as
  `--mem, --mem-per-cpu, and --mem-per-gpu are mutually exclusive`.

What to check:

- do not combine `sbatch --mem=...` with entrypoints that already request `#SBATCH --mem-per-cpu=...`,
- unset inherited `SBATCH_MEM`, `SBATCH_MEM_PER_CPU`, or `SBATCH_MEM_PER_GPU` variables before submission if your shell profile exports them,
- current entrypoints include `#SBATCH --ignore-pbs` because some Slurm builds also parse embedded `#PBS` lines unless told not to.

### Job is killed near the requested memory limit

Symptom:

- SLURM/PBS/UGE reports an out-of-memory kill even though a tool was configured with the job's requested RAM.

What to check:

- scheduler allocation is recorded as `GG_MEM_TOTAL_GB`,
- tools should receive `GG_MEM_TOOL_GB`, which leaves `GG_MEM_TOOL_RESERVE_GB` for container-side overhead and helper processes,
- increase the scheduler memory request if `GG_MEM_TOOL_GB` is still too small for the tool,
- reduce `GG_MEM_TOOL_RESERVE_GB` only after confirming the extra processes and tool overhead fit inside the job allocation.

### Stage skipped unexpectedly

Symptom:

- the wrapper finishes quickly and logs that no run was necessary.

Likely causes:

- the corresponding `run_*` flag is disabled,
- output-exists skip logic detected already-finished stage outputs,
- an upstream prerequisite directory or table is missing and the stage chose to skip rather than fail hard.

Examples:

- `gg_gene_summary_entrypoint.sh` skips database generation if logical `stat_tree` or `stat_branch` inputs are absent from both live files and ZIP storage,
- `gg_progress_summary_core.sh` skips orthogroup summary generation if the selected gene-count table is absent; AMAS inputs may be live or ZIP-backed.

GeneGalleon records content-based provenance manifests for generated artifacts,
including gene-family, genome-annotation, transcriptome, species-tree,
orthogroup, and gene-summary stages. Existing legacy outputs without a manifest
are reused rather than regenerated. When all declared inputs are still
available, GeneGalleon backfills a manifest from the current files and
parameters so later changes can be detected. Optional outputs record either a
present or absent state, so a tool that legitimately produces no result (for
example, no valid CSUBST foreground branch combination) is still complete.

`artifact_stale_policy` controls a detected input, output, or parameter mismatch:

- `stop` (default) prints the mismatched family, stage, manifest, and reason,
  then exits before modifying outputs;
- `reuse` continues with the stale output and does not rewrite its manifest;
- `rebuild` enables the affected stage and regenerates it without confirmation.

Only declared data content and output-affecting parameters participate in the
freshness decision. Tool, container, and GeneGalleon versions are retained as
diagnostics and do not cause regeneration. Raw and ZIP-backed managed output
directories use the same logical content digest, so storage conversion alone
does not cause regeneration.

When `run_gene_family_database_build=1`, inspect
`gene_summary/<source>/<source>_artifact_provenance_audit.tsv` for the exact HOG,
step, status, and reason. `legacy_untracked` is informational and does not stop
database generation. `semantic_mismatch` means the CSUBST and `stat_branch`
branch IDs do not identify the same descendant-tip clades and the affected
upstream chain must be regenerated.

Developers can inventory shell stages that rely on existence-, mtime-, count-,
or readiness-only cache guards, plus Python functions that return/continue on
an existing output, with:

```bash
python workflow/support/audit_cache_guards.py \
  --baseline workflow/tests/data/artifact_cache_guard_baseline.txt \
  --output-tsv cache_guard_audit.tsv
```

The committed baseline records known migration debt without hiding it from the
TSV. It is currently empty. CI scans both `workflow/core` and
`workflow/support` and fails when a new unprovenanced cache guard is introduced.

### Gene-family storage conversion is pending

Symptom:

- `gg_gene_evolution` or progress summary refuses to start because
  `.gg_store/storage-conversion.pending` exists. A store created by the first
  experimental ZIP implementation may use
  `.gg_archives/storage-conversion.pending` until `migrate-layout` completes.

What to do:

- stop old array jobs that may still write below the affected query2family or orthogroup root,
- inspect the saved phase and counts with `gg_gene_family_archive.sh conversion-status --root <output-root>`,
- rerun the same `gg_gene_family_archive.sh convert-storage ... --to zip|raw` command,
- or add `--resume` to require that the pending conversion exists and matches the requested direction,
- raw-to-ZIP resume cleans partial shards and repairs a pending index automatically; use the standalone `repair` command only when no conversion marker exists,
- do not remove the conversion marker manually; it records the required target direction.

### A one-off parameter override is ignored

Symptom:

- an inline environment value does not change the effective entrypoint config.

What to check:

- use the entrypoint-scoped prefix rather than the raw shell variable name,
  for example `GG_GENOME_EVOLUTION_RUN_CAFE=1` or
  `GG_GENE_EVOLUTION_MODE_GENE_EVOLUTION=orthogroup`,
- confirm that the parameter is registered in
  `workflow/support/gg_entrypoint_config_vars.sh`,
- inspect the effective config summary in the job log.

See [Configuration and Common Parameters](configuration-and-common-parameters.md)
for the complete prefix scheme.

### Taxonomy or database cache issues

Symptom:

- taxonomy-related steps hang, fail, or repeatedly rebuild caches.

Relevant locations:

- `workspace/downloads/ete_taxonomy`
- `workspace/downloads/pfam`
- `workspace/downloads/locks`

What to do:

- inspect whether another job is still using the shared lock,
- remove stale lock files only when you are sure no active job is using the cache,
- rerun after cleaning only the specific broken cache subtree rather than the whole workspace.

### Shared-filesystem locking fails or behaves inconsistently

ZIP-backed gene-family stores and their active-task/temporary-directory guards
use atomic shared-filesystem namespace locks, not `flock`. This protects active
tasks even on Lustre `localflock` mounts. The `.gg_store_locks/` tree (including
`.lock.namespace-v1` sidecars) is persistent coordination state outside the
movable store metadata and must not be deleted by cleanup tools.
Other workflows' array finalizers still require cross-node `flock` semantics.

What to check:

- verify the same namespace lock from two compute nodes on the actual workspace
  filesystem; a writer must block while any reader is active,
- confirm the Lustre, BeeGFS, or NFS mount options with the site administrator;
  workflows still using `flock` cannot use node-local locking,
- stop all users of a gene-family store before switching from a flock-only
  runtime; old and new protocols must never run concurrently,
- namespace locks are never stolen based on age or a remote PID. After a crash,
  timeout is fail-closed: verify that every job and reader using that exact
  store has stopped, retain the owner records for audit, and reconcile only the
  abandoned lock. Never remove a lock just to make a retry proceed.

### A replacement transcriptome task reaches the shared summary concurrently

Current GeneGalleon releases serialize the complete multispecies-summary
transaction across scheduler job IDs. Each writer builds in a private sibling
directory and atomically publishes the finished `annotation_summary` tree.
Seeing `getcwd() failed` together with a missing `expression.imputed.tsv`
indicates an older workflow release that allowed one replacement job to remove
another writer's current directory. Update the immutable workflow release
before retrying; repeatedly resubmitting that older release can reproduce the
same race.

### A large public FASTQ download stops at the same byte offset

The public-original fallback keeps a hidden `*.download.part` file and resumes
gzip downloads with a validated HTTPS byte range. Do not delete that partial
file between bounded retries. GeneGalleon verifies `Content-Range`, total size,
the provider checksum when available, gzip/FASTQ structure, and then atomically
publishes the completed file. A server that ignores Range is handled as a full
restart; a checksum mismatch or non-FASTQ response is never appended to the
saved partial file.

### `gg_genome_evolution` protein mode does not behave as expected

Symptoms:

- protein mode exits before OrthoFinder starts,
- DNA species-tree steps seem to disappear,
- changing `species_genetic_code.tsv` appears to have no effect.

What to check:

- `input_sequence_mode` is actually set to `protein` in `workflow/gg_genome_evolution_entrypoint.sh`,
- either `workspace/input/species_protein/` or `workspace/input/species_cds/` exists,
- if both exist, protein mode prefers `species_protein` and ignores `species_genetic_code.tsv`,
- if you want per-species translation overrides to matter, remove or relocate `species_protein` and let the run translate from `species_cds`,
- protein mode can include correctly translated proteins from lineages with
  different genetic codes,
- protein mode intentionally disables DNA-only and codon-sequence-based steps
  such as DNA IQ-TREE, IQ2MC, and MCMCtree.

Useful log messages:

- `species_genetic_code.tsv is ignored because species_protein inputs are provided`
- `Ignoring species_protein inputs in cds mode`
- `Shared protein input signature changed for ...`

Those messages indicate GeneGalleon is applying the current mode rules and invalidating stale derived proteins when the effective inputs change.

### `species_cds` validation fails

Symptoms:

- a genome annotation, genome evolution, or input-generation run exits during CDS validation,
- the log reports that sequence names are inconsistent with the species name parsed from the FASTA filename,
- the log reports duplicate sequence names or prohibited characters.

What to check:

- every `workspace/input/species_cds` FASTA filename starts with the species label, such as `GENUS_SPECIES_...`,
- every FASTA header follows the required `GENUS_SPECIES_GENEID` pattern,
- the `GENUS_SPECIES` prefix in each sequence ID matches the prefix parsed from that file's name,
- sequence IDs are unique within each FASTA,
- sequence IDs do not contain whitespace or special punctuation such as `|`.

Common invalid examples:

- `>AT1G08465.1`: missing the required species prefix,
- `>Arabidopsis_thaliana|AT1G08465.1`: contains the prohibited `|` character,
- `>Arabidopsis_thaliana_AT1G08465.1` inside `Oryza_sativa_IRGSP.fa.gz`: sequence prefix and filename prefix do not match.

### Input-generation summary lacks `taxid` or genetic-code metadata

Symptom:

- `gg_input_generation_species.tsv` exists, but taxonomy columns are blank.

Current summary file:

- `workspace/output/input_generation/gg_input_generation_species.tsv`

Older one-off exports or ad hoc copies may still be named `species_summary.tsv`,
but that is no longer the canonical filename written by the wrapper.

What to check:

- species names must match an NCBI scientific name or synonym closely enough for taxonomy lookup,
- the shared ETE taxonomy DB under `workspace/downloads/ete_taxonomy/` must be readable,
- the wrapper logs a warning and continues when taxonomy cache preparation fails, so blank metadata is not fatal by itself.

Current behavior:

- nuclear and mitochondrial genetic codes come from NCBI taxonomy metadata,
- plastid genetic code is a lineage-based best-effort default, not a direct species-specific NCBI field.

### Optional analyses never run

Symptom:

- expression, GFF, genome, or trait analyses are silently absent from downstream outputs.

What to check:

- confirm the matching optional inputs exist under `workspace/input/`,
- confirm the corresponding `run_*` flag is enabled,
- confirm any mode-specific prerequisite stage has already been completed.

### Species-tree-dependent outputs missing in gene-family runs

Symptom:

- tree-guided downstream outputs are absent in `gg_gene_evolution_entrypoint.sh`.

What to check:

- `workspace/output/species_tree` exists and contains the expected summary trees,
- the gene-family run is configured to use tree-dependent analyses,
- orthogroup mode prerequisites from `gg_genome_evolution_entrypoint.sh` are already present.

## Log locations worth checking

- `workspace/output/versions/*.log`: tool/runtime inventory after successful runs
- scheduler stdout/stderr from `#SBATCH --output` / `#$ -cwd` jobs
- `workspace/output/debug_entrypoint_logs/summary.tsv`: summary table from `workflow/gg_all_entrypoints_debug.sh`
- `workspace/output/input_generation/gg_input_generation_runs.tsv`: run history for input-generation runs

## Safe cleanup targets

If you need to clean state, prefer targeted cleanup:

- stale provider download cache under `workspace/output/input_generation/tmp/`
- stale taxonomy cache under `workspace/downloads/ete_taxonomy/`
- stale locks under `workspace/downloads/locks/`

Avoid deleting the whole workspace unless you really want to rebuild everything.
