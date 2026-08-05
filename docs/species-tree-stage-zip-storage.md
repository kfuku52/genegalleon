# Species-Tree Stage ZIP Storage

`gg_genome_evolution` can store its high-file-count species-tree intermediates
as ordinary, visible ZIP files. ZIP mode is the default and is controlled by:

```bash
GG_COMMON_SPECIES_TREE_OUTPUT_STORAGE=zip # zip|files|raw
GG_COMMON_SPECIES_TREE_ZIP_COMPRESSION=adaptive # adaptive|deflate|store
GG_COMMON_SPECIES_TREE_ZIP_COMPRESSION_LEVEL=6 # 0-9
```

`raw` is an alias for `files`. Project-specific values may instead be set as
`species_tree_output_storage`, `species_tree_zip_compression`, and
`species_tree_zip_compression_level` in
`workflow/gg_genome_evolution_entrypoint.sh`.

When `files`/`raw` mode is selected, all managed ZIPs are materialized near the
start of `gg_genome_evolution`, after stale-input invalidation but before later
workflow validation. Thus an early failure does not leave an only-partly
converted files-mode workspace. ZIP mode remains demand-driven and restores
only the directory needed by the current stage.

## Managed directories

Only these direct children of `workspace/output/species_tree` are managed:

- `single_copy_cds_fasta`
- `single_copy_mafft`
- `single_copy_trimal`
- `single_copy_iqtree_pep`
- `single_copy_iqtree_dna`

In ZIP mode, a stage that has reached its archive point is represented by a
neighboring file with the same name, for example:

```text
workspace/output/species_tree/
  single_copy_iqtree_dna.zip
```

The archive contains a top-level `single_copy_iqtree_dna/` directory, so normal
macOS extraction recreates the familiar directory name. GeneGalleon restores a
managed directory immediately before a stage needs it and archives it again
after its last consumer. This limits the number of simultaneously live stage
files. Completion checks count matching entries directly from the ZIP central
directory, so a rerun with no pending work does not expand and recompress
completed archives. An exit handler also attempts to archive managed raw
directories after a normal error exit. It deliberately leaves them raw if
background producers are still active. A non-catchable termination such as
`SIGKILL` can leave a raw or partial directory; the next managed operation
removes its own stale partial paths and resumes from the preserved raw files.
If a recoverable interrupted operation leaves both forms, workflow-controlled
packing materializes and merges them first; an existing live file wins.

Small, externally referenced directories such as `species_tree_summary`,
`single_copy_astral_*`, `concatenated_*`, `constrained_tree`, and `mcmctree_*`
remain raw. Archiving them would save few inodes while forcing unrelated
workflows to materialize common summary files.

## Existing workspace conversion

Convert all managed directories to ZIP:

```bash
bash workflow/gg_species_tree_archive.sh convert-storage \
  --root workspace/output/species_tree \
  --to zip
```

Restore them to the historical raw layout:

```bash
bash workflow/gg_species_tree_archive.sh convert-storage \
  --root workspace/output/species_tree \
  --to raw
```

Inspect or CRC-check archives:

```bash
bash workflow/gg_species_tree_archive.sh status \
  --root workspace/output/species_tree

bash workflow/gg_species_tree_archive.sh verify \
  --root workspace/output/species_tree
```

`verify` is deep by default. `--quick` validates the GeneGalleon marker and ZIP
member inventory without reading payload bytes. `convert-storage --dry-run`
reports raw/ZIP byte counts and a conservative temporary-space estimate.
`--available-bytes` supplies the applicable user/project quota remainder for
both packing and materialization. ZIP-to-raw preflight includes
allocation-block-rounded file bytes, directory overhead, and required inodes.
The option can also be supplied directly to `materialize`. A failure leaves
the source ZIP unchanged. `--progress-interval` controls JSON progress records
written to standard error during packing, materialization, and deep
verification.
Species-tree command results are JSON by default and also accept the common
`--json` flag for command-line symmetry with the gene-family manager.

To audit or convert gene-family and species-tree storage together, use
`workflow/gg_workspace_storage.sh` as documented in
`docs/workspace-storage-management.md`.

Add one or more `--directory NAME` options to limit any command to selected
managed directories.

## Manual edits

Before manually adding, deleting, or replacing files, restore the directory:

```bash
bash workflow/gg_species_tree_archive.sh materialize \
  --root workspace/output/species_tree \
  --directory single_copy_iqtree_dna
```

Edit the restored directory and either rerun `gg_genome_evolution` or run the
conversion command with `--to zip`. Packing refuses an ambiguous state in which
both the raw directory and its ZIP exist. Materialization can safely merge that
state after an interrupted extraction; an existing live file wins over the
older archived copy.

Do not run a manual conversion while `gg_genome_evolution` is actively writing
the same species-tree directory. Archive-manager commands serialize with each
other, and the workflow waits for all producers in a stage before packing, but
an unrelated process cannot participate in that lock. If a source changes
during packing, GeneGalleon refuses to replace it; avoiding concurrent manual
conversion also prevents a producer from recreating a raw directory just after
packing finishes.

ZIP creation is atomic. GeneGalleon verifies the new archive CRC and confirms
that the source inventory did not change before replacing an older ZIP and
removing the raw directory. ZIPs without GeneGalleon archive metadata are not
modified or extracted by this manager.

Standard ZIP stores regular files, directory entries, Unix permission bits,
and modification times (timestamps outside ZIP's 1980–2107 range are clamped).
File ownership, ACLs, hard-link identity, and extended attributes are not
portable ZIP metadata and are not retained. GeneGalleon-generated files do not
normally rely on those attributes.
