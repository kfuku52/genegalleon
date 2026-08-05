# Migrating legacy unzipped workspaces to ZIP storage

`gg_workspace_storage.sh` audits or converts every managed high-file-count
output store in one GeneGalleon workspace. It is the preferred command for an
existing project because it performs a complete preflight before changing any
target, verifies the result, and writes durable JSON and TSV reports.

This is the recommended migration procedure before running a current
ZIP-default GeneGalleon checkout against a workspace created by a release that
stored one ordinary file per artifact. The same command also audits an already
ZIP-backed workspace and converts it back to raw files when needed.

The short migration sequence is:

```bash
# Stop jobs that can write this workspace first.
bash workflow/gg_workspace_storage.sh audit \
  --workspace /path/to/legacy-project
bash workflow/gg_workspace_storage.sh convert \
  --workspace /path/to/legacy-project \
  --to zip \
  --dry-run \
  --available-bytes QUOTA_REMAINING_BYTES
bash workflow/gg_workspace_storage.sh convert \
  --workspace /path/to/legacy-project \
  --to zip \
  --workers 4 \
  --available-bytes QUOTA_REMAINING_BYTES
```

Review the dry-run JSON report before running the third command. The sections
below explain catalog overrides, unmatched files, quota estimates, target
selection, rollback, and verification in detail.

## Managed and excluded paths

The command detects these workspace layouts:

- a current workspace directory containing `output/` and `input/`;
- a project directory whose child `workspace/` contains `output/` and `input/`;
- a project directory containing `gfe_data/`;
- `workspace/output` or `gfe_data` supplied directly as `--workspace`.

Within the detected output root, it manages only `query2family/`,
`orthogroup/`, and the five documented high-file-count directories below
`species_tree/`: `single_copy_cds_fasta`, `single_copy_mafft`,
`single_copy_trimal`, `single_copy_iqtree_pep`, and
`single_copy_iqtree_dna`.

It never converts or edits `orthofinder/`, `genome_evolution/busco_*`,
unrecognized files or subdirectories below a gene-family output root, or small
shared species-tree summary directories outside the managed list. Excluded
top-level paths and unmatched gene-family file counts are recorded in the
report. At most 20 unmatched path examples are included so a report cannot
itself become a high-file-count artifact.

## Before running a conversion

Stop `gg_gene_evolution`, `gg_genome_evolution`, and progress-summary jobs that
can write the workspace. Current GeneGalleon processes participate in storage
locks, but older running code cannot be forced to honor locks introduced by a
newer checkout.

Run a read-only audit first:

```bash
bash workflow/gg_workspace_storage.sh audit \
  --workspace /path/to/project/workspace
```

The default is a quick verification: ZIP central directories, GeneGalleon
markers/manifests, indexes, member names, recorded sizes, and recorded CRC
values are checked without reading every payload byte. For a complete CRC and
SHA256 read of every gene-family member and a complete CRC read of every
species-tree member, use:

```bash
bash workflow/gg_workspace_storage.sh audit \
  --workspace /path/to/project/workspace \
  --verification deep
```

Deep verification can read hundreds of gigabytes. Use it periodically or when
storage corruption is suspected; it is not needed merely to repeat the full
payload verification that is already performed before every newly created ZIP
is committed.

## If a ZIP-default job is submitted before this migration

The workflow does not silently run this complete workspace migration before
starting the requested biological job. Instead, migration is incremental:

- `gg_gene_evolution` creates ZIP-store metadata, reuses existing raw files for
  normal stage skip/staleness decisions, and makes a successful family eligible
  for batched archiving; after a controlled failure it attempts to archive that
  family's new or changed partial outputs;
- its completion sweep may also archive other recognized legacy families that
  satisfy the legacy completion test; `gg_progress_summary` performs a similar
  sweep with a one-file minimum;
- remaining raw files continue to override older ZIP members, so a temporary
  mixed raw/ZIP store is supported by transparent readers;
- `gg_genome_evolution` packs each managed high-file-count species-tree stage
  when that stage reaches its normal archive point;
- unrelated or unrecognized files remain raw.

Each archive write is atomic and preserves the raw source when packing fails,
but this incremental path does not provide the workspace command's complete
catalog/conflict scan, workspace-wide quota plan, strict-unmatched option, or
durable before/after reports. It can therefore leave a partially migrated
workspace for a while. Run the audit and dry-run below first when practical,
especially for a large legacy project. If a job has already started this
incremental migration, stop all writers and run the same workspace command; it
accepts and completes a supported mixed raw/ZIP store.

## Dry-run conversion

Inventory every target, resolve its input catalog, estimate temporary space,
and write reports without changing managed output storage:

```bash
bash workflow/gg_workspace_storage.sh convert \
  --workspace /path/to/project/workspace \
  --to zip \
  --dry-run
```

For raw-to-ZIP conversion, `raw_zip_peak_new_bytes` is a conservative estimate
of the extra space needed while raw sources and in-progress ZIPs coexist. It
accounts for the configured number of concurrent writers and possible small
deflate expansion.

Filesystem free space on a shared filesystem may be much larger than the
applicable user or project quota. After checking the quota, pass the remaining
bytes explicitly:

```bash
bash workflow/gg_workspace_storage.sh convert \
  --workspace /path/to/project/workspace \
  --to zip \
  --dry-run \
  --available-bytes 500000000000
```

The effective available space is the smaller of filesystem free space and
`--available-bytes`. The workspace preflight evaluates the targets in their
actual sequential order: it carries forward a conservative maximum net growth
from earlier targets and rejects the complete run if a later target could
exceed the initial headroom. The JSON `quota_plan.operations` array shows this
calculation. A real conversion refuses to start when the workspace-wide
estimate is larger. The estimate is intentionally conservative but cannot
predict other jobs consuming the same quota after preflight.

## Convert the complete workspace to ZIP

```bash
bash workflow/gg_workspace_storage.sh convert \
  --workspace /path/to/project/workspace \
  --to zip \
  --workers 4 \
  --available-bytes 500000000000
```

The command has two phases:

1. It inventories every selected target and checks catalogs, path conflicts,
   symlinks, pending conversions, unmatched-file policy, and temporary space.
2. Only if all targets pass does it convert each target, verify the result, and
   audit the final physical/logical state.

If a preflight fails, no selected target is converted by default. If a failure
occurs after conversion has begun, the current atomic ZIP operation preserves
its raw source or committed ZIP, and remaining targets are marked `not-run`.
Rerunning the same command resumes a compatible interrupted gene-family
conversion. `--continue-on-error` is available for independent batch
maintenance, but is not recommended for a first conversion.

New ZIP creation always receives one full payload read before the raw source is
removed. The default workspace-level post-conversion verification is `quick`,
which avoids immediately rereading the same payload. Use `--verification deep`
when an additional independent full read is desired.

### Catalog discovery

For `query2family`, the command normally uses
`workspace/input/query_gene/`. Override it with `--query-dir`.

For `orthogroup`, it normally uses
`orthofinder/Orthogroups_filtered/Orthogroups.GeneCount.selected.tsv`.
Override it with `--genecount`.

If old raw outputs belong to a family that is no longer in the current input,
list one family ID per line and supplement the current catalog:

```bash
bash workflow/gg_workspace_storage.sh convert \
  --workspace /path/to/project/workspace \
  --to zip \
  --query-family-id-file old-query-ids.txt \
  --orthogroup-family-id-file old-orthogroup-ids.txt
```

Without `--strict-unmatched`, unrecognized files/directories are reported and
left unchanged. `--strict-unmatched` makes any such gene-family file a blocking
preflight issue.

### Select targets

All detected targets are selected by default. Repeat `--target` to restrict the
operation:

```bash
bash workflow/gg_workspace_storage.sh convert \
  --workspace /path/to/project/workspace \
  --target orthogroup \
  --target species_tree \
  --to zip
```

Valid values are `query2family`, `orthogroup`, and `species_tree`.

### Compression and shard settings

```bash
bash workflow/gg_workspace_storage.sh convert \
  --workspace /path/to/project/workspace \
  --to zip \
  --compression adaptive \
  --compression-level 6 \
  --max-files-per-shard 5000 \
  --workers 4
```

`adaptive` deflates ordinary text/FASTA/tree/table files and stores already
compressed members. `workers` is limited to 1–4 to avoid unbounded filesystem
load.

By default, a gene-family subdirectory is finalized as one visible
`<subdirectory>.zip`, regardless of size. A warning is included in reports
when its logical input exceeds 20 GiB. Change or disable the warning threshold:

```bash
--large-zip-warning-bytes 53687091200  # warn above 50 GiB
--large-zip-warning-bytes 0            # disable the warning
```

To avoid a very large single ZIP, set a positive final-ZIP limit:

```bash
--max-final-zip-bytes 21474836480
```

Above that logical size, GeneGalleon retains multiple visible, meaningful
archives such as:

```text
orthogroup/archives/stat_branch/stat_branch.part-000123.zip
```

Transparent readers treat the parts as one logical subdirectory. The default
remains an unlimited final ZIP because one `<subdirectory>.zip` is simpler to
move and extract on macOS. The same persistent workflow setting is available
as `GG_COMMON_GENE_FAMILY_FINAL_ZIP_MAX_BYTES` or the entrypoint variable
`gene_family_final_zip_max_bytes`.
The same limit is used when grouping members into parts, so each part stays at
or below the configured logical size unless one individual member is itself
larger than the limit.

## Convert the complete workspace back to raw files

```bash
bash workflow/gg_workspace_storage.sh convert \
  --workspace /path/to/project/workspace \
  --to raw \
  --available-bytes 500000000000
```

This restores current logical versions, including live overrides and logical
deletions. `--available-bytes` is optional but recommended on quota-managed
storage. ZIP-to-raw estimates use allocation-block-rounded file sizes and also
check inode availability; a preflight failure leaves every selected ZIP
unchanged. Archive control metadata is retained so the project can safely
return to ZIP mode later. To discard gene-family archive history/control
metadata after materialization:

```bash
bash workflow/gg_workspace_storage.sh convert \
  --workspace /path/to/project/workspace \
  --to raw \
  --pure-raw
```

For gene-family roots, `--pure-raw` also removes the generated visible README
and archive-status snapshot. It does not change the species-tree result because
its neighboring ZIP is consumed directly during materialization.

## Reports

Every successful dry-run, audit, or real conversion creates two visible files
under `storage_reports/` next to the detected output root:

```text
storage_reports/
  gg_storage_20260805T123456Z_12345.json
  gg_storage_20260805T123456Z_12345.tsv
```

Use `--report-dir` and `--report-prefix` to control their location/name.

The JSON report is the authoritative audit record. It contains workspace,
output, and input roots; the GeneGalleon Git commit and dirty-worktree flag; command options;
verification depth; timestamps; explicitly excluded paths; per-target
catalogs and preflight state; logical managed file/byte counts before
conversion; physical ZIP/metadata/raw file and byte counts afterward;
temporary-space estimates; unmatched examples; conversion/verification
results; warnings; and errors.

The TSV is a compact one-row-per-target summary plus a `TOTAL` workspace row,
suitable for spreadsheets. Its
before columns are logical managed contents; its after columns are physical
managed storage. Unchanged shared/unmatched files are deliberately excluded
from those comparison columns and remain described in JSON.

Reports are written atomically. A nonzero exit status means at least one target
failed or was not run. The terminal output always prints the report paths when
a report could be completed.

## Scheduler example

Run storage conversion as a single non-array job. Supply site-specific account,
partition, time, and memory settings rather than copying a cluster name from a
different host:

```bash
#!/usr/bin/env bash
#SBATCH --job-name=gg-storage
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=gg-storage-%j.out
#SBATCH --error=gg-storage-%j.err

set -euo pipefail

bash /path/to/genegalleon/workflow/gg_workspace_storage.sh convert \
  --workspace /path/to/project/workspace \
  --to zip \
  --workers "${SLURM_CPUS_PER_TASK:-1}" \
  --available-bytes 500000000000
```

Use an absolute GeneGalleon path or load the project runtime explicitly. The
storage command itself uses Python standard-library functionality and does not
require Apptainer for ZIP management.
