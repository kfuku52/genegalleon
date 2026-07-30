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

## 2. Install a versioned SIF

Build or obtain the SIF outside SHIROKANE, then transfer the completed file.
SHIROKANE currently provides Apptainer 1.2.4, so the deployment path does not
convert the large OCI image on the cluster.

Example on the machine that holds the SIF:

```bash
sha256=$(shasum -a 256 genegalleon.sif | awk '{print $1}')
tag="sha256-${sha256:0:12}"
ssh shirokane-kf52 'mkdir -p "$HOME/genegalleon/incoming"'
rsync --partial --progress genegalleon.sif \
  "shirokane-kf52:genegalleon/incoming/genegalleon-${tag}.sif"
```

Submit the validation and installation job from the same local shell:

```bash
ssh shirokane-kf52 \
  "cd \"\$HOME/genegalleon\" && qsub -terse -v \
GG_SHIROKANE_PREBUILT_SIF=incoming/genegalleon-${tag}.sif,\
GG_SHIROKANE_SIF_TAG=${tag},\
GG_SHIROKANE_SIF_SHA256=${sha256} \
workflow/sites/shirokane/gg_shirokane_prepare_sif.sh"
```

The preparation job:

- requests the SHIROKANE `ljob` resource,
- validates the transferred image with `apptainer inspect` on a compute node,
- requires and verifies the expected SHA-256 before and after installation,
- executes `uname -m` inside the SIF and rejects non-amd64 images,
- serializes installations that use the same immutable tag,
- writes a versioned SIF under `containers/`,
- writes a relocatable SHA-256 file,
- updates `./genegalleon.sif` atomically as a relative symlink after validation.

It refuses to replace an existing regular `./genegalleon.sif`. Move or rename
such a file explicitly before switching to versioned symlinks.

Monitor the job with:

```bash
qstat
```

## 3. Verify an AGE submission

The SHIROKANE submit helper defaults to `-l ljob`. Verification does not submit
a job. Normal submissions also use `qsub -terse`, so successful output contains
the AGE job ID:

```bash
bash workflow/sites/shirokane/gg_shirokane_submit.sh \
  --entrypoint gg_gene_evolution_entrypoint.sh \
  --tasks 1 \
  --verify
```

Use `--dry-run` to print the `qsub` command without contacting AGE.
Task IDs, range endpoints, and range steps must be positive integers; decreasing
ranges and step zero are rejected before `qsub` is called.

## 4. Submit a workflow

Array-aware example:

```bash
bash workflow/sites/shirokane/gg_shirokane_submit.sh \
  --entrypoint gg_gene_evolution_entrypoint.sh \
  --tasks 1-N
```

Fixed single-task example:

```bash
bash workflow/sites/shirokane/gg_shirokane_submit.sh \
  --entrypoint gg_genome_evolution_entrypoint.sh
```

The helper preserves each entrypoint's default `def_slot` and `s_vmem`
directives while adding `ljob`. Override them only when the data require it:

```bash
bash workflow/sites/shirokane/gg_shirokane_submit.sh \
  --entrypoint gg_transcriptome_generation_entrypoint.sh \
  --tasks 1-N \
  --slots 4 \
  --mem-per-slot 32G
```

On SHIROKANE, `s_vmem` is memory per slot. The last example requests
`4 * 32G = 128G` for each array task.

GeneGalleon forwards the allocated slot count to `OMP_NUM_THREADS`,
`OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS`, and `NUMEXPR_NUM_THREADS` inside the
container. This prevents numerical libraries from creating more worker threads
or reserving more virtual memory than AGE allocated.

## Resource selection

The helper accepts:

- `--resource ljob`: default; Thin nodes, maximum runtime 62 days.
- `--resource mjob`: Thin nodes, maximum runtime 2 days.
- `--resource lmem`: Fat nodes, maximum runtime 14 days.
- `--resource exclusive`: dedicated group queue when the account has an
  allocation.

Use the default `ljob` for GeneGalleon production work. Use `lmem` only when
the requested total memory or observed workload requires a Fat node.

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

Start with `--tasks 1 --verify`, then submit one real task before launching the
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

GeneGalleon outputs remain in the shared workspace. SIF installation copies to
`containers/<name>.sif.partial.<JOB_ID>` and removes an incomplete partial file
on exit. Runtime workflows may use job-local `/tmp` for transient work.

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
