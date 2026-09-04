# Scheduler and Array Semantics

Wrappers include scheduler headers for:

- SLURM
- UGE
- PBS

On current Slurm, `sbatch` can also read `#PBS` directives unless `--ignore-pbs` is set.
GeneGalleon entrypoints therefore include `#SBATCH --ignore-pbs` so the embedded PBS header
remains available for `qsub` without conflicting with Slurm resource flags.

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

The UGE directives in each entrypoint default to the SHIROKANE `ljob`
resource. Use standard `qsub` options for array ranges and resource overrides:

```bash
qsub -terse -t 1-N workflow/gg_gene_evolution_entrypoint.sh
```

See the [SHIROKANE AGE Guide](shirokane-age.md) for SIF transfer, AGE
verification, and resource selection.

When no scheduler metadata is present, GeneGalleon falls back to local defaults such as:

- `GG_TASK_CPUS=1`
- `GG_ARRAY_TASK_ID=1`
- `GG_JOB_ID=1`

## Scheduler normalization

Internally, GeneGalleon normalizes SLURM/PBS metadata to scheduler-neutral variables:

- `GG_TASK_CPUS`
- `GG_ARRAY_TASK_ID`
- `GG_JOB_ID`
- `GG_ARRAY_JOB_ID`

Those normalized values are what core scripts use downstream, regardless of the original scheduler.

`GG_JOB_ID` identifies an individual job. `GG_ARRAY_JOB_ID` identifies the
array as a whole: on Slurm it comes from `SLURM_ARRAY_JOB_ID`, while the
individual `SLURM_JOB_ID` can differ for every array element. Both values are
forwarded into the container. Shared finalizers group ready markers by the
array ID so the last completed task runs the multispecies summary once.

Slurm memory can be requested with either `--mem-per-cpu` or `--mem`.
GeneGalleon reads the corresponding `SLURM_MEM_PER_CPU` or
`SLURM_MEM_PER_NODE` value in MiB and derives `GG_MEM_TOTAL_GB` in whole GiB,
rounding down after summing a per-CPU allocation. Explicit `GG_MEM_TOTAL_GB`
or `GG_MEM_PER_CPU_GB` settings, and the legacy memory aliases, take precedence
over automatic detection. For example, `--cpus-per-task=8 --mem=8G` gives an
8 GiB total allocation, rather than multiplying a default per-CPU value by 8.

The startup log prints a runtime summary that shows:

- the detected scheduler kind,
- the original scheduler variables,
- the normalized values forwarded into the container.

## Fixed single-task wrappers

These wrappers are intended to run as one task:

- `gg_genome_evolution_entrypoint.sh`
- `gg_gene_summary_entrypoint.sh`
- `gg_progress_summary_entrypoint.sh`

`gg_input_generation_entrypoint.sh` is not fixed single-task only. It has staged modes controlled by `input_generation_mode`:

- `single`: run the whole wrapper in one process
- `array_prepare`: create `workspace/output/input_generation/tmp/task_plan.json`
- `array_worker`: run one task-plan row per scheduler array index
- `array_finalize`: merge worker shards and run shared validation and summary steps

## Array-size rules

- `gg_input_generation_entrypoint.sh`:
  - `input_generation_mode=array_worker`: set the scheduler array range to the number of species tasks recorded in `workspace/output/input_generation/tmp/task_plan.json`
  - `input_generation_mode=array_prepare` and `array_finalize`: run once, not as an array
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

`workflow/support/gg_site_runtime.sh` automatically selects a site profile and
can apply behavior such as:

- changing into `PBS_O_WORKDIR`,
- adding site-specific container bind mounts,
- selecting a modified container shell command.

That logic is automatic; you do not normally need to edit scheduler headers just to get the runtime bindings right.

Current profiles cover SHIROKANE, the National Institute of Genetics, NHR@FAU,
and a default fallback. See [Site Runtime Profiles](site-runtime-profiles.md)
for detection and runtime details. `GG_SITE_PROFILE` can force a profile when
automatic hostname detection is not appropriate.

The SHIROKANE profile also initializes Environment Modules and loads
`apptainer` on AGE compute nodes. The login node is used for submission and
file management, not direct workflow execution.

## Debugging scheduler problems

The fastest places to look are:

- the wrapper stdout/stderr file declared in the scheduler header,
- the startup runtime summary printed by the wrapper,
- `workspace/output/versions/*.log` after successful jobs.
