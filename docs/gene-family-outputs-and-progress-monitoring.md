# Gene-Family Outputs and Progress Monitoring

This page expands the short stage references for
`gg_gene_evolution_entrypoint.sh`, especially for users familiar with the old
`gfe_geneFamilyPhylogeny` workflow.

## Query2family input recap

In `mode_gene_evolution=query2family`:

- place one file per family under `workspace/input/query_gene/`,
- the file basename becomes the family/task ID,
- the scheduler array size is the number of files in that directory.

Each input file should contain query genes for one gene family only. If genes
from different families are combined in one file, GeneGalleon treats them as one
family-level task, which can artificially merge unrelated family phylogenies into
one unnatural tree and PDF output.

Accepted query-file forms are:

- plain gene ID list,
- in-frame CDS FASTA,
- protein FASTA.

See [Input Conventions](input-conventions.md) for concrete examples.

Family-scoped materialization checks a live filename's family before reading
its metadata or content. Other families' transient provenance files and shared
lockfiles cannot abort a selected family's restoration. Archived records retain
their explicit family identity; corruption or unreadability of a selected
artifact remains an error rather than being silently skipped.

## Main per-family outputs

If the input file is `workspace/input/query_gene/2_WOX`, the main outputs are
typically written under `workspace/output/query2family/` as:

- `query_gene/2_WOX`
- `query_aa_fasta/2_WOX_query.aa.fa.gz`
- `cds_fasta/2_WOX_cds.fa.gz`
- `mafft/`, `clipkit/`, and/or `trimal/` alignment products depending on the enabled cleaning steps
- `iqtree_tree/2_WOX_iqtree.nwk`
- `rooted_tree/2_WOX_root.nwk`
- `dated_tree/2_WOX_dated.nwk` and `dated_tree_log/2_WOX_dated.log.txt` when `run_tree_dating=1`
- `stat_branch/2_WOX_stat.branch.tsv`
- `stat_tree/2_WOX_stat.tree.tsv`
- `tree_plot/2_WOX_tree_plot.pdf`

When `run_orthogroup_extraction=1`, GeneGalleon extracts the seed-containing
orthogroup from the rooted homolog tree before GeneRax. The extracted tree and
FASTA are written below `orthogroup_extraction_nwk/`,
`orthogroup_extraction_rooted_nwk/`, and `orthogroup_extraction_fasta/` and are
then used as the GeneRax inputs in the same run when `run_generax=1`. A missing
required extraction input is an error rather than a silently skipped
refinement.

When `run_expression_trait_pgls=1`, the same family receives the methods named
by `pgls_methods` plus common comparison and audit files. RSC members are:

- `rsc_status/2_WOX_rsc.status.tsv`: `ok` or `not_estimable`, with the reason
  and counts used for screening,
- `rsc_regression/2_WOX_rsc.regression.tsv`: coefficients, uncertainty,
  inference status, evolutionary-model parameters, and event counts,
- `rsc_reconciliation/`, `rsc_gene_contrasts/`, and
  `rsc_species_contrasts/`: the event mapping and both sides of the contrast,
- `rsc_response_*` and `rsc_predictor_*`: replicate summaries and propagated
  sampling covariance,
- `rsc_random_effects/`, `rsc_sensitivity/`, and `rsc_trait_origins/`:
  hierarchical effects and optional lineage/origin diagnostics,
- `rsc_audit/` and `rsc_log/`: NWKIT provenance and execution diagnostics.

Species-tree comparison members are:

- `pgls_species_nwkit/`: NWKIT ordinary PGLS after within-sample paralog
  aggregation,
- `pgls_species_rphylopars/`: the matched Rphylopars comparator; unsupported
  categorical or covariance combinations remain header-only and are explained
  in the status table,
- `pgls_comparison/`: long-form coefficient rows from RSC and both species
  methods, with `analysis_method`, `aggregation`, `estimand`, and explicit
  comparability notes,
- `pgls_method_status/` and `pgls_method_audit/`: one place to distinguish
  successful, non-estimable, and unrequested methods and to record engine
  versions,
- `pgls_species_expression_summary/` and
  `pgls_species_expression_audit/`: the actual species/sample expression used,
  paralog coverage, aggregation rule, and number of declared cross-paralog
  sampling-covariance pairs used,
- `pgls_species_response_*` and `pgls_species_predictor_*`: species-tip
  replicate summaries and propagated sampling covariance.

Every combined bundle table carries `analysis_id`, so predictor-specific
shape-parameter estimates, contrasts, and sampling covariances remain
distinguishable when predictors are screened separately.

All members are written even when a family is not estimable; optional tables
then contain headers only. This makes family completion and ZIP-backed storage
unambiguous. The expression-trait bundle is staged completely and published as
a recoverable transaction; a failed move or caught interruption restores the
previous complete bundle before provenance is recorded. `stat_tree/*.tsv`
includes the family status and best overall RSC
summary plus namespaced per-response/per-term fields, all beginning with
`rsc_`. It also includes bounded `pgls_species_nwkit_*` and
`pgls_species_rphylopars_*` status/best-row fields. To keep `stat_tree` bounded
for large screens, it contains counts and
fields from the best usable (successfully converged) row only, under
`rsc_best_*`; it does not flatten every RSC result row into a new group of
columns. The full `rsc_regression` table remains the authoritative result when
several responses or predictors were fitted. Rows with failed inference or
failed optimizers cannot become the reported best row. The best row is still
chosen by its raw p-value, recorded as `rsc_best_p_value_raw` (and the analogous
species-method field). The unsuffixed `*_best_p_value` and `rsc_min_p_value`
fields are Holm-adjusted across every usable response, predictor term,
analysis, and paralog aggregation represented by that family/method. The same
summary also records the number and scope of tests plus Holm and
Benjamini-Hochberg values explicitly; use the full result tables for a
different prespecified multiplicity family.

In `mode_gene_evolution=orthogroup`, the same style of per-family outputs is
written under `workspace/output/orthogroup/`.

## ZIP-backed storage

`workflow/gg_common_params.sh` defaults
`GG_COMMON_GENE_FAMILY_OUTPUT_STORAGE` to `zip`. The alternative value
`files` preserves the historical one-file-per-artifact layout for newly
produced artifacts; `raw` is accepted as an alias for `files`. Existing ZIP
shards remain readable in either mode; a rerun in `files` mode restores the
selected family as needed and leaves it as normal files.

ZIP mode applies only below:

- `workspace/output/query2family/`
- `workspace/output/orthogroup/`

It does not archive `workspace/output/orthofinder/` or
`workspace/output/genome_evolution/busco_*`.

Routine batch archiving makes a family eligible when its gene-evolution task
reaches the end successfully and records a completion state. Current state is stored in at
most 256 bucketed JSON files instead of an ever-growing per-event log. This is
independent of optional outputs, so runs with `run_summary=0` or
`run_tree_plot=0` can still be archived. Outputs created before
completion-state tracking use the historical
`stat_branch` + `stat_tree` + `tree_plot` test as a one-time legacy fallback.
Starting a rerun immediately invalidates the previous completion state until
that rerun succeeds.

Storage conversion and controlled failed-run cleanup may archive an incomplete
family without marking it complete. Before a rerun starts, GeneGalleon
materializes only that family's archived artifacts at their historical paths,
preserves their mtimes, and then applies the normal stage skip/staleness tests.
Consequently, a legacy run converted after MAFFT, for example, can resume from
the later stages instead of starting over. A controlled failure records
`failed`, removes unchanged materialized copies, and archives new or changed
partial artifacts. After an abrupt kill, any remaining live files override the
older ZIP members and are consumed on the next rerun.

Workspaces from releases that used dot-separated subdirectories and filenames
(for example `amas.cleaned/HOG0000010.amas.cleaned.tsv`) are exposed through
the current underscore-based logical names. On rerun, only the selected
family is materialized under the corresponding current directory, and the
workflow accepts the historical uncompressed FASTA names as valid stage
inputs. This allows a partially completed pre-ZIP family to continue after a
raw-to-ZIP conversion without rewriting every archive member first.

The other output subdirectories do not have to contain every expected family.
For example, files for completed families can be moved from a partially
complete `mafft/` directory into ZIP storage while incomplete-family files
remain live. Each array task archives only its own family so it cannot hold the
global reader-maintenance lock while sweeping the entire catalog. The progress
summary separately flushes any remaining completed outputs.

`GG_COMMON_GENE_FAMILY_ZIP_MIN_BATCH_FILES` is retained as a deprecated
configuration compatibility setting and no longer controls array-task cleanup.

Gene-family tasks hold a shared lock for their family while they may read or
write live outputs. Lock files use 16 fixed hash stripes, so lock
metadata itself cannot grow with the number of families. Archive maintenance
takes only the relevant buckets exclusively and nonblockingly: an active
family is skipped while completed families in other buckets can still be
archived. A separate reader-maintenance lock prevents a ZIP shard from being
replaced while a downstream reader has it open.

The bounded internal metadata under `.gg_store/` costs at most 16 family-lock
files, 16 state-lock files, 256 family-index JSON files, and 256 state-index
JSON files, plus a small number of global/per-subdirectory files. Existing stores
created by an older build can remove unused 256-way lock files with
`optimize-metadata` after every job using that output root has stopped.

Immutable shards are compacted automatically before the number of referenced
shards in one logical subdirectory becomes large. Index updates and obsolete
shard reclamation are committed one subdirectory at a time, bounding peak
space during archive and compaction. A durable update marker makes readers
fail closed if a process stops between updates to the denormalized index
views; `repair` then rebuilds every view from the ZIP manifests. Family-bucket
and subdirectory indexes, plus an epoch used to invalidate long-lived reader
caches, avoid opening every ZIP manifest for normal listing, status, and
task-level materialization.
Repeated transparent reads reuse one ZIP reader per worker thread until that
worker moves to another archive. This avoids reparsing a large central
directory once per member in summaries and database generation; an index
epoch change closes the cached readers before refreshed artifacts are used.

Temporary task data remains under the historical
`workspace/output/{query2family,orthogroup}/tmp/` roots, avoiding an
unbounded burst in a node-wide system temporary filesystem. When a rerun
temporarily restores archived inputs into the historical paths, it writes a
receipt before computation begins. A failed or killed run subsequently
removes only receipt-listed files whose size and SHA-256 still match the
restored member. Files modified by the failed run or by a person remain live;
stale-task maintenance applies the same receipt cleanup.

Successful jobs
still remove their own task directory immediately. In ZIP mode, opportunistic
and progress-summary maintenance removes failed-job task directories older
than `GG_COMMON_GENE_FAMILY_TMP_RETENTION_DAYS` (default `7`) only after it
can obtain the affected family lock. It also retains at most
`GG_COMMON_GENE_FAMILY_TMP_MAX_DIRS` recent failed-task directories (default
`100`), `GG_COMMON_GENE_FAMILY_TMP_MAX_BYTES` bytes (default 100 GiB), and
`GG_COMMON_GENE_FAMILY_TMP_MAX_FILES` files (default `100000`). Set an
individual limit to `0` to disable it. Cleanup recognizes only task directory
names beginning with a numeric array index and an underscore.

Archives are ordinary, visible ZIP files. While a run is active, immutable
parts are stored under
`<gene-family-root>/archives/<subdirectory>/<subdirectory>.part-NNNNNN.zip`.
Once an orthogroup catalog is fully complete and no family task holds a lock,
progress-summary maintenance consolidates each logical output set once into
`<gene-family-root>/<subdirectory>.zip` and removes its parts. Raw-to-ZIP
conversion writes previously raw subdirectories directly to this finalized
layout, in parallel up to `--workers`, rather than creating parts and then
recompressing them. Each member keeps its logical
`<subdirectory>/<filename>` path. GeneGalleon summaries, database
generation, HGT plots, csubst site plots, and orthogroup-based GRAmPA read
these archives transparently or materialize only the family inputs required
by an external tool. HGT and csubst materialization directories carry
per-run locks; after SIGKILL, a later invocation reclaims only run directories
whose owner lock is no longer held.

`query2family` remains open for later inputs. A finalized
`<subdirectory>.zip` is an immutable base; newly completed queries go into
visible parts below `archives/`. The logical reader overlays the parts on the
base, and an explicit `finalize` safely rebuilds the single ZIP when desired.
GeneGalleon writes `README_GENE_FAMILY_OUTPUTS.txt` and `ARCHIVE_STATUS.tsv`
at the output root so a user browsing with Finder, FileZilla, or `ls` can see
where each logical output set is stored. Only locks, indexes, tombstones, and
transaction markers remain hidden under `.gg_store/`.
Rows whose ZIP storage also has visible overrides or shared files use a
`+live` storage suffix (for example `finalized+live`) and report the live count
separately.

The default `adaptive` compression deflates ordinary text and sequence files
while storing already compressed `.gz`, `.zip`, image, and PDF members without
redundant recompression. Configure routine workflow archiving with:

- `GG_COMMON_GENE_FAMILY_ZIP_COMPRESSION=adaptive|deflate|store`
- `GG_COMMON_GENE_FAMILY_ZIP_COMPRESSION_LEVEL=0..9`
- `GG_COMMON_GENE_FAMILY_ZIP_WORKERS=1..4`
- `GG_COMMON_GENE_FAMILY_FINAL_ZIP_MAX_BYTES=0..` (`0` keeps the one-final-ZIP behavior)

`store` prioritizes packing files into a small number of inodes over reducing
bytes. Worker concurrency is deliberately capped at four to avoid an
unbounded burst of metadata and read traffic on a shared filesystem.

### Adding, replacing, and deleting files manually

Do not edit `*.zip` in place. The archive manifest intentionally detects
direct ZIP edits as an error.

To add or replace an artifact, put it at its historical path. A live file
always overrides the archived version:

```bash
cp replacement.tsv \
  workspace/output/query2family/stat_branch/2_WOX_stat.branch.tsv
```

`ARCHIVE_STATUS.tsv` is a snapshot rather than a filesystem watcher; its file
modification time is the snapshot publication time. Per-task `archive-family`
does not refresh it or scan the whole store, even after successful cleanup.
Parent/batch archive operations still refresh it. After individual task
completions or manual `cp`, `mv`, or `rm` operations, refresh its counts once with:

```bash
bash workflow/gg_gene_family_archive.sh refresh-status \
  --root workspace/output/query2family
```

Perform direct `cp` changes only when the affected gene-family task is idle.
The managed `delete`, `undelete`, and `restore` commands take the affected
family lock before changing logical state. For a query2family live file that
has never been archived, supply `--family-id` explicitly because its family
cannot be recovered from archive metadata. The command fails without changing
the logical store when a family ID cannot be inferred safely or when an
explicit ID disagrees with existing archive metadata or the logical filename.

A later archive pass stores the replacement as a new immutable ZIP shard. To
inspect the logical view and see whether each artifact is live or archived:

```bash
bash workflow/gg_gene_family_archive.sh list \
  --root workspace/output/query2family \
  --subdir stat_branch
```

For deletion, use the management command instead of `rm`. It records a
tombstone so an older archived version cannot reappear:

```bash
bash workflow/gg_gene_family_archive.sh delete \
  --root workspace/output/query2family \
  --family-id 2_WOX \
  --path stat_branch/2_WOX_stat.branch.tsv
```

The operation can be reversed, or the archived artifact can be restored as a
normal file:

```bash
bash workflow/gg_gene_family_archive.sh undelete \
  --root workspace/output/query2family \
  --family-id 2_WOX \
  --path stat_branch/2_WOX_stat.branch.tsv

bash workflow/gg_gene_family_archive.sh restore \
  --root workspace/output/query2family \
  --family-id 2_WOX \
  --path stat_branch/2_WOX_stat.branch.tsv
```

Consolidate the current base and parts into one ordinary ZIP per logical
subdirectory only while the affected gene-family jobs are stopped:

```bash
bash workflow/gg_gene_family_archive.sh finalize \
  --root workspace/output/query2family \
  --mode query2family \
  --query-dir workspace/input/query_gene
```

Finalization refuses to replace an unrelated pre-existing
`<subdirectory>.zip`; move or rename that file explicitly after inspecting it.
Raw-to-ZIP conversion checks every such collision during preflight, before it
creates a conversion marker or writes any managed ZIP. Dry-run and
`storage-status` report both `conflicting_final_zip_files` and each
`conflicting-final-zip` path.

Verify every archive member with CRC and SHA-256 checks after copying a run:

```bash
bash workflow/gg_gene_family_archive.sh verify \
  --root workspace/output/query2family \
  --progress-interval 30
```

This deep verification reads every payload byte. Use `--quick` to validate ZIP
inventories, manifests, recorded CRC/size values, and indexes without rereading
every member. ZIP creation itself always performs one complete CRC read before
raw sources are removed.
Add `--json` to `verify`, `convert-storage`, `storage-status`, or
`conversion-status` for one machine-readable result object; progress remains
on standard error.

`verify` also cross-checks the family and subdirectory index views and rejects
a missing authoritative index, an index reference to a missing shard, an
unfinished index-update marker, mixed query2family/orthogroup shard modes, and
unreferenced/orphan ZIP shards. If an
interrupted archive, compaction, repair, or manual operation left valid ZIP
manifests but incomplete or missing indexes, rebuild them explicitly:

```bash
bash workflow/gg_gene_family_archive.sh repair \
  --root workspace/output/query2family \
  --progress-interval 10
```

Add `--remove-orphans` only after reviewing the reported orphan paths.
Both explicit `repair` and `finalize` commands emit the same periodic progress
fields as conversion. The orthogroup finalization automatically triggered by
`gg_progress_summary` does so as well.

### Converting an existing workspace

Inspect a query2family root before changing its storage:

```bash
bash workflow/gg_gene_family_archive.sh convert-storage \
  --root workspace/output/query2family \
  --mode query2family \
  --query-dir workspace/input/query_gene \
  --to zip \
  --dry-run
```

Remove `--dry-run` to archive every recognized family-owned live artifact,
including outputs of incomplete families. Conversion does not create a
completion state. Intentional shared files such as common files below
`parameters/` remain live and are reported as `shared_raw_files`; they do not
trigger `--strict-unmatched`. Other unrecognized files remain live and are
reported as `unmatched_live_files`; `--strict-unmatched` aborts before
conversion when those files exist. A parameter filename that identifies an
orthogroup outside the supplied catalog is likewise unmatched, rather than
being hidden in the shared count. If query inputs were
removed after their outputs were created, list the missing family IDs one per
line and pass `--family-id-file`.

The command does not edit `gg_common_params.sh`. Set
`GG_COMMON_GENE_FAMILY_OUTPUT_STORAGE=zip` before resubmitting converted
families. A rerun in `files`/`raw` mode remains readable, but its new outputs
stay as ordinary files and GeneGalleon prints a mode-mismatch warning.

For orthogroup output, use the selected gene-count table as a read-only family
catalog:

```bash
bash workflow/gg_gene_family_archive.sh convert-storage \
  --root workspace/output/orthogroup \
  --mode orthogroup \
  --genecount workspace/output/orthofinder/Orthogroups_filtered/Orthogroups.GeneCount.selected.tsv \
  --to zip
```

`--query-dir` unambiguously selects query2family mode and `--genecount`
selects orthogroup mode, so `--mode` can be omitted in these examples. The
first ZIP write records `store.json` with the mode, current catalog count and
family-ID digest, compression settings, and catalog source paths. Commands
against an existing ZIP store then infer the mode and reject an explicit
conflicting value. Adding query input files updates the catalog count and
digest without invalidating older shards.

Large conversions, repairs, and finalizations emit a flushed `progress` line
immediately, at each shard or subdirectory transition, and every 30 seconds
while a ZIP operation is still running. Direct final-ZIP writes report
`active_subdirs`, aggregate
`active_zip_bytes`, and the most recently updated `current_subdir` and
`current_zip_bytes`; completed ZIP writers waiting for index commit are not
counted as active. The mandatory post-write CRC read uses
`phase=verifying-final-zip` with completed/total member and logical-byte
counts, so a large final ZIP does not appear stalled after its size stops
growing. The final whole-store verification and the standalone `verify`
command report completed/total ZIP, member, and logical-byte counts under
`phase=verifying`. Finalization reports `finalized_subdirs` and
`total_subdirs`. Interrupted conversions use the explicit
`phase=repairing-index` phase with verified ZIP/member/byte counts instead of
remaining at `phase=starting`. Each line includes elapsed seconds and the
latest phase/counts. Use
`--progress-interval 10` for more frequent heartbeats or `0` for transition
lines only. Conversion writes its current phase, heartbeat, completed
subdirectories, and committed counts into the durable conversion marker.

Compression and bounded parallel writing can also be selected per conversion:

```bash
bash workflow/gg_gene_family_archive.sh convert-storage \
  --root workspace/output/orthogroup \
  --genecount workspace/output/orthofinder/Orthogroups_filtered/Orthogroups.GeneCount.selected.tsv \
  --to zip \
  --compression adaptive \
  --compression-level 1 \
  --workers 2
```

The command reads but never modifies `workspace/output/orthofinder/`.
Progress-summary wrappers likewise write their AMAS-augmented copy to the
summary output directory; `orthogroup_output_summary.py` accepts
`--updated-genecount-out` for callers that need an explicit destination and
otherwise writes beside `--out`, never beside the input gene-count table.

Stores created by the first experimental ZIP implementation can be moved out
of the hidden mixed `.gg_archives/` layout without changing their logical
contents. Stop every job using the output root, then run:

```bash
bash workflow/gg_gene_family_archive.sh migrate-layout \
  --root workspace/output/orthogroup

bash workflow/gg_gene_family_archive.sh finalize \
  --root workspace/output/orthogroup \
  --mode orthogroup \
  --genecount workspace/output/orthofinder/Orthogroups_filtered/Orthogroups.GeneCount.selected.tsv
```

The migration moves ZIP payload into visible `archives/`, moves internal
metadata into `.gg_store/`, rebuilds indexes from verified manifests, and
writes the user-facing README and status table. `finalize` then creates the
single `<subdirectory>.zip` files.

To restore the current logical view to ordinary files, first inspect the
uncompressed byte and file counts and then run the conversion:

```bash
bash workflow/gg_gene_family_archive.sh convert-storage \
  --root workspace/output/query2family \
  --to raw \
  --dry-run

bash workflow/gg_gene_family_archive.sh convert-storage \
  --root workspace/output/query2family \
  --to raw
```

ZIP-to-raw conversion materializes and verifies every non-deleted logical
artifact before removing ZIP payload. Tombstones become absence and live
overrides remain authoritative. Bounded family-state and lock metadata remain
by default; add `--pure-raw` only when an exact pre-ZIP physical layout is
required and discarding that control metadata, `README_GENE_FAMILY_OUTPUTS.txt`,
and `ARCHIVE_STATUS.tsv` is acceptable. The preflight
estimate rounds each restored file to the filesystem allocation-block size and
includes missing output directories. Filesystem-wide free-space and inode
values can still exceed a Shirokane user or project quota, so compare the
reported `raw_materialize_allocated_bytes` and `raw_peak_new_files` with the
applicable quota before a large conversion. `--available-bytes` applies to
ZIP-to-raw as well as raw-to-ZIP conversion; a quota shortfall is rejected
before a resumable conversion marker or raw output is created.

Raw-to-ZIP dry-runs report `raw_zip_peak_new_bytes`, including concurrent
writer and possible compression-expansion overhead. Pass the applicable quota
remainder with `--available-bytes`; conversion then uses the smaller of that
value and filesystem free space.

Subdirectories above `--large-zip-warning-bytes` are identified before
conversion. `--max-final-zip-bytes` can retain meaningful
`archives/<subdir>/<subdir>.part-N.zip` files instead of consolidating an
oversized `<subdir>.zip`. The default limit is 0 (unlimited).
Applying a positive limit to an already ZIP-backed store also splits an older
oversized final ZIP into the readable part names.

For complete existing projects, use the workspace-level command documented in
`docs/workspace-storage-management.md`; it audits all gene-family and
species-tree targets before changing any of them and writes JSON/TSV receipts.

Both directions write `storage-conversion.pending`. A stopped conversion is
resumed automatically by running the same command again. Add `--resume` when
you want the command to require an existing matching interrupted conversion;
it fails instead of silently starting a new conversion when no marker exists.
Before resuming raw-to-ZIP, GeneGalleon removes partial shards and rebuilds the
authoritative indexes from committed manifests, verifying CRC and SHA-256 in
one read pass and making shard/index crash windows recoverable.
`gg_gene_evolution` and progress
summary refuse to start while the marker remains. Stop old array jobs before
conversion; writers from code predating this lock protocol cannot be prevented
from creating new live files, although signature checks keep such files from
being deleted.

To inspect the current physical/logical counts independently:

```bash
bash workflow/gg_gene_family_archive.sh storage-status \
  --root workspace/output/query2family \
  --query-dir workspace/input/query_gene
```

The status distinguishes `zip_files`, bounded GeneGalleon
`metadata_files`, and their managed total in `physical_store_files`.
Intentional common files are counted separately as `shared_raw_files` and
`shared_raw_bytes`, while `unmatched_live_files` is reserved for files that
need review.

To inspect only the resumable marker and archive format without reading the
family catalog:

```bash
bash workflow/gg_gene_family_archive.sh conversion-status \
  --root workspace/output/query2family
```

After upgrading a project and stopping all old/new gene-family jobs, remove
obsolete lock files created by the former 256-way layout:

```bash
bash workflow/gg_gene_family_archive.sh optimize-metadata \
  --root workspace/output/query2family
```

To compact existing shards explicitly:

```bash
bash workflow/gg_gene_family_archive.sh compact \
  --root workspace/output/query2family \
  --mode query2family
```

Logical deletion uses a tombstone log; it is automatically reduced to the
latest record per path after it reaches its rotation threshold. To physically
rewrite shards so tombstoned members and archived versions hidden by live
overrides no longer occupy space, use the irreversible maintenance command:

```bash
bash workflow/gg_gene_family_archive.sh purge \
  --root workspace/output/query2family \
  --mode query2family \
  --query-dir workspace/input/query_gene
```

`purge --drop-unlisted` additionally removes archived members and family-state
records for query inputs no longer present (or orthogroups absent from the
specified gene-count table). Run `verify` and keep a backup before purge.

To produce a normal directory containing only the current logical versions
(live overrides applied and tombstoned files omitted), use `export-current`.
This is the safe form to transfer for ordinary Finder/Archive Utility use when
the store has replacements or deletions:

```bash
bash workflow/gg_gene_family_archive.sh export-current \
  --root workspace/output/query2family \
  --destination-root ../query2family-current
```

The export destination must not exist or must be empty. This prevents files
deleted from the logical store from surviving in a reused export directory.

Directly deleting a live override with `rm` is ambiguous: GeneGalleon will
then expose an older archived version if one exists. Use the managed `delete`
command when logical deletion is intended.

## How to read `tree_plot`

`tree_plot` is generated from `stat_branch/*.tsv` and summarizes many
gene-family attributes around the inferred tree. In the current default
GeneGalleon configuration, the panel order is:

1. gene tree
2. expression heatmap
3. expression point plot
4. cluster-membership panel
5. local synteny panel
6. tip labels
7. signal peptide summary
8. transmembrane-domain summary
9. intron-count panel
10. protein-domain panel from RPS-BLAST
11. alignment panel
12. promoter-motif panel from FIMO
13. MEME motif summary
14. ortholog-context panel

Practical interpretation:

- the leftmost tree is the anchor and carries the branch length, support, and
  branch-color settings chosen by `treevis_*` parameters,
- expression, cluster, and synteny panels help interpret lineage-specific
  expansions or expression shifts,
- structural panels summarize domains, introns, signal peptides, and membrane
  predictions on the same row order as the tree tips,
- if an upstream analysis was disabled, or the corresponding inputs were not
  available, the associated panel may be blank or minimally populated.

`stat.branch.tsv` is the master table that collects per-branch and per-tip
annotations for plotting. `stat.tree.tsv` is the paired tree-level summary.
Its `original_num_site` and `cleaned_num_site` fields are the alignment-column
counts before and after site trimming, respectively. The corresponding
`original_num_seq`, `original_len_*`, `cleaned_num_seq`, and `cleaned_len_*`
fields likewise summarize the untrimmed and trimmed analysis alignments.
See [Example Plots](example-plots.md) for a compact generated tree-plot
example from the bundled quick-start data.

## Gene-family progress summaries

GeneGalleon provides a built-in summary wrapper for orthogroup and
query2family runs, plus transcriptome-scale runs:

```bash
cd workflow
bash gg_progress_summary_entrypoint.sh
```

This writes:

- `workspace/orthogroup_summary.tsv`
- `workspace/query2family_summary.tsv`
- `workspace/transcriptome_assembly_summary.tsv`

For orthogroup runs, `orthogroup_summary.tsv` is useful because:

- it adds `GG_ARRAY_TASK_ID`, which is the row index to resubmit,
- it appends AMAS-derived alignment statistics such as
  `Parsimony_informative_sites_clean`,
- it adds one `1/0` completion column per visible output subdirectory under
  `workspace/output/orthogroup/`.

Practical reading tips:

- rows with `Parsimony_informative_sites_clean == 0` cannot produce normal
  IQ-TREE-based downstream outputs,
- rows with `0` in late-stage columns such as `stat_tree`, `stat_branch`, or
  `tree_plot` are the first candidates to inspect or rerun.

Transcriptome summary rows are species based rather than orthogroup based.
The wrapper also prints incomplete species IDs inferred from missing
`amalgkit_getfastq/<species>/*.safely_removed` markers.

## Query2family completion audit

For query2family runs, `query2family_summary.tsv` is generated when
`workspace/output/query2family` and `workspace/input/query_gene` exist.
Rows follow the same sorted input-file order used for query2family array
tasks, so `GG_ARRAY_TASK_ID` can be used directly for resubmission.

The summary adds one `1/0` completion column per visible output subdirectory
under `workspace/output/query2family/`. It also appends AMAS-derived alignment
statistics when `amas_original` or `amas_cleaned` outputs are present.

For large query2family runs, inspect late-stage completion markers such as
`tree_plot`, `stat_branch`, or `stat_tree`.

Archive-aware manual status:

```bash
bash workflow/gg_gene_family_archive.sh status \
  --root workspace/output/query2family \
  --mode query2family \
  --query-dir workspace/input/query_gene
```

For task-level resubmission decisions, use `query2family_summary.tsv`; it
contains the archive-aware `tree_plot`, `stat_branch`, and `stat_tree`
completion columns together with `GG_ARRAY_TASK_ID`.

## Rerunning incomplete tasks

Once you know the task IDs:

Slurm:

```bash
sbatch --array=17,42,105 workflow/gg_gene_evolution_entrypoint.sh
```

UGE:

```bash
qsub -t 17,42,105 workflow/gg_gene_evolution_entrypoint.sh
```

Local rerun of one task:

```bash
GG_ARRAY_TASK_ID=17 GG_TASK_CPUS=4 bash workflow/gg_gene_evolution_entrypoint.sh
```

If your site prefers one submission per task rather than an explicit list-style
array, submit a short loop around the same wrapper.
