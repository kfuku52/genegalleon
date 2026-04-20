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

### `sbatch` rejects mutually exclusive memory options

Symptom:

- `sbatch` fails before the job starts with an error such as
  `--mem, --mem-per-cpu, and --mem-per-gpu are mutually exclusive`.

What to check:

- do not combine `sbatch --mem=...` with entrypoints that already request `#SBATCH --mem-per-cpu=...`,
- unset inherited `SBATCH_MEM`, `SBATCH_MEM_PER_CPU`, or `SBATCH_MEM_PER_GPU` variables before submission if your shell profile exports them,
- current entrypoints include `#SBATCH --ignore-pbs` because some Slurm builds also parse embedded `#PBS` lines unless told not to.

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
