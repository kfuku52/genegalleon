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

- confirm the cluster login environment actually exposes the runtime,
- if your site uses a module system, load the container runtime before launching,
- if you are on macOS, note that SIF creation is normally skipped by default and many runs are expected to happen on Linux hosts.

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

### Stage skipped unexpectedly

Symptom:

- the wrapper finishes quickly and logs that no run was necessary.

Likely causes:

- the corresponding `run_*` flag is disabled,
- output-exists skip logic detected already-finished stage outputs,
- an upstream prerequisite directory or table is missing and the stage chose to skip rather than fail hard.

Examples:

- `gg_gene_database_core.sh` skips database generation if `stat_tree/` or `stat_branch/` is missing,
- `gg_progress_summary_core.sh` skips orthogroup summary generation if the selected gene-count table or AMAS directories are absent.

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
