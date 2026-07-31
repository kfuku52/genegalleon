# SHIROKANE AGE Guide

This guide covers running GeneGalleon on the Human Genome Center SHIROKANE
system with Altair Grid Engine (AGE) and Apptainer.

Official SHIROKANE references:

- [AGE basics](https://gc.hgc.jp/uge/uge-essential/)
- [AGE resources](https://gc.hgc.jp/uge/uge-resource/)
- [AGE queues](https://gc.hgc.jp/uge/queue/)
- [Multiple CPU cores and per-slot memory](https://gc.hgc.jp/uge/uge-def-slot/)
- [Installed software](https://gc.hgc.jp/util-info/installed/)
- [Home-disk file-count guidance](https://gc.hgc.jp/2026/03/file_limit/)

## Runtime model

GeneGalleon uses SHIROKANE as follows:

- AGE runs the workflow on a compute node.
- The SHIROKANE site profile loads the default `apptainer` module on that node.
- GeneGalleon normalizes `NSLOTS`, `JOB_ID`, `SGE_TASK_ID`, and `s_vmem` into
  the scheduler-neutral `GG_*` variables used inside the container.
- The repository, workspace, and SIF remain on a shared filesystem path that is
  visible from every compute node.

Do not run a GeneGalleon workflow directly on a login node. Apptainer is
provided through the compute-node module environment.

## 1. Place the repository on SHIROKANE

Clone or transfer the repository into the SHIROKANE home filesystem, then run
all commands below from the repository root:

```bash
cd /path/to/genegalleon
```

The default layout expects:

- `./genegalleon.sif`
- `./workflow/`
- `./workspace/`

## 2. Place the SIF

Build or obtain the SIF outside SHIROKANE, then transfer the completed file.
Use an amd64 image because SHIROKANE compute nodes are amd64. The workflow
expects the final file at `./genegalleon.sif`.

Example on the machine that holds the SIF:

```bash
sha256=$(shasum -a 256 genegalleon.sif | awk '{print $1}')
ssh shirokane-kf52 'mkdir -p "$HOME/genegalleon"'
rsync --partial --progress genegalleon.sif \
  "shirokane-kf52:genegalleon/genegalleon.sif.incoming"
ssh shirokane-kf52 \
  "cd \"\$HOME/genegalleon\" && \
printf '%s  %s\n' '${sha256}' genegalleon.sif.incoming | sha256sum -c - && \
mv -f genegalleon.sif.incoming genegalleon.sif"
```

## 3. Verify an AGE submission

The UGE header in each entrypoint already contains the SHIROKANE `ljob`
resource. Verification does not submit a job:

```bash
qsub -verify -t 1 workflow/gg_gene_evolution_entrypoint.sh
```

## 4. Submit a workflow

Array-aware example:

```bash
qsub -terse -t 1-N workflow/gg_gene_evolution_entrypoint.sh
```

Fixed single-task example:

```bash
qsub -terse workflow/gg_genome_evolution_entrypoint.sh
```

Each entrypoint now carries SHIROKANE defaults for `def_slot`, `s_vmem`, and
`ljob`. Override them with standard `qsub` options only when the data require
it:

```bash
qsub -terse -t 1-N -pe def_slot 4 -l s_vmem=32G,ljob \
  workflow/gg_transcriptome_generation_entrypoint.sh
```

On SHIROKANE, `s_vmem` is memory per slot. The last example requests
`4 * 32G = 128G` for each array task.

GeneGalleon forwards the allocated slot count to `OMP_NUM_THREADS`,
`OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, and `NUMEXPR_NUM_THREADS` inside the
container. This prevents numerical libraries from creating more worker threads
or reserving more virtual memory than AGE allocated.

## Resource selection

- `ljob`: default; Thin nodes, maximum runtime 62 days.
- `mjob`: Thin nodes, maximum runtime 2 days.
- `lmem`: Fat nodes, maximum runtime 14 days.
- `exclusive`: dedicated group queue when the account has an
  allocation.

Use the default `ljob` for GeneGalleon production work. Use `lmem` only when
the requested total memory or observed workload requires a Fat node. Specify an
alternative together with `s_vmem`, for example `-l s_vmem=64G,lmem`.

AGE chooses SHIROKANE queues through these resources. Do not add a direct
`-q ljobs.q` selection.

## Array sizing

Follow the rules in
[Scheduler and Array Semantics](scheduler-and-array-semantics.md). In
particular:

- query2family uses the number of files in `workspace/input/query_gene`;
- orthogroup mode uses the number of selected orthogroup rows;
- genome annotation uses the number of species CDS inputs;
- transcriptome generation uses the number of selected species/SRA input
  units.

Start with `qsub -verify -t 1`, then submit one real task before launching the
full range.

## Storage checks

SHIROKANE limits both storage capacity and inode count. Check the group inode
usage before and after large workflows:

```bash
lfsq
```

Avoid placing more than 5,000 files in one directory. Archive completed,
fine-grained intermediate outputs when they no longer need to be accessed
individually.

GeneGalleon outputs and `genegalleon.sif` remain on the shared filesystem.
Runtime workflows may use job-local `/tmp` for transient work.

## Expected startup evidence

A successful AGE launch should report values similar to:

```text
site profile = shirokane
container runtime = apptainer
scheduler=uge
GG_TASK_CPUS=4
GG_MEM_PER_CPU_GB=8
GG_MEM_TOTAL_GB=32
GG_ARRAY_TASK_ID=1
```

The exact CPU and memory values follow the submission request.
