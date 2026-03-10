# Scheduler and Array Semantics

Wrappers include scheduler headers for:

- SLURM
- UGE
- PBS

You can submit wrappers with scheduler commands such as `sbatch` or `qsub`, or invoke them directly with `bash`.

## Direct execution vs scheduler execution

Direct execution:

```bash
bash workflow/gg_genome_evolution_entrypoint.sh
```

Scheduler submission examples:

```bash
sbatch workflow/gg_gene_evolution_entrypoint.sh
```

```bash
qsub workflow/gg_gene_evolution_entrypoint.sh
```

When no scheduler metadata is present, GeneGalleon falls back to local defaults such as:

- `GG_TASK_CPUS=1`
- `GG_ARRAY_TASK_ID=1`
- `GG_JOB_ID=1`

## Scheduler normalization

Internally, GeneGalleon normalizes SLURM/PBS metadata to scheduler-neutral variables:

- `GG_TASK_CPUS`
- `GG_ARRAY_TASK_ID`
- `GG_JOB_ID`

Those normalized values are what core scripts use downstream, regardless of the original scheduler.

The startup log prints a runtime summary that shows:

- the detected scheduler kind,
- the original scheduler variables,
- the normalized values forwarded into the container.

## Fixed single-task wrappers

These wrappers are intended to run as one task:

- `gg_input_generation_entrypoint.sh`
- `gg_gene_convergence_entrypoint.sh`
- `gg_gene_database_entrypoint.sh`
- `gg_genome_evolution_entrypoint.sh`
- `gg_progress_summary_entrypoint.sh`

## Array-size rules

- `gg_gene_evolution_entrypoint.sh`:
  - `mode_gene_evolution=orthogroup`: number of rows in `workspace/output/orthofinder/Orthogroups_filtered/Orthogroups.GeneCount.selected.tsv` excluding the header
  - `mode_gene_evolution=query2family`: number of files in `workspace/input/query_gene`
- `gg_genome_annotation_entrypoint.sh`: number of input species CDS files
- `gg_transcriptome_generation_entrypoint.sh`: number of species input units for the selected mode

## Practical submission notes

- update the array range in the scheduler header so it matches the current input count,
- for local runs, wrappers effectively behave like a single-task launch unless you export scheduler-style variables yourself,
- if an array run exits immediately, confirm both the array range and the underlying input files/tables.

## Site-specific behavior

`workflow/support/gg_site_runtime.sh` can apply site-aware behavior such as:

- changing into `PBS_O_WORKDIR`,
- adding site-specific container bind mounts,
- selecting a modified container shell command.

That logic is automatic; you do not normally need to edit scheduler headers just to get the runtime bindings right.

## Debugging scheduler problems

The fastest places to look are:

- the wrapper stdout/stderr file declared in the scheduler header,
- the startup runtime summary printed by the wrapper,
- `workspace/output/versions/*.log` after successful jobs.
