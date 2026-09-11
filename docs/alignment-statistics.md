# Alignment statistics

Gene-evolution workflows use `cdskit stats --mode alignment` with CDSKIT 0.31.0
or newer. Rebuild older containers before using this workflow. Container source
defaults continue to follow `container/source_branches.env`.

The switches are `run_alignment_stats_original` and
`run_alignment_stats_cleaned`. Outputs are stored in `alignment_stats_original`
and `alignment_stats_cleaned`, with filenames such as
`HOG0000010_alignment_stats.original.tsv`. The corresponding provenance steps
also use `alignment_stats_*`. There are no old-name aliases in normal workflows.

DNA (`--seq_type dna`) and protein (`--seq_type aa`) inputs are supported. seqkit
still decompresses alignments. The statistical column names are unchanged.
`GC_content` uses a 0–1 scale and per-taxon averaging, not a pooled percentage;
proteins report `NA`. See the
[CDSKIT statistics documentation](https://github.com/kfuku52/cdskit/wiki/cdskit-stats)
for character and site definitions. Summary readers reread current statistics
on each run, including when an earlier aggregate table already contains values.

## Migrating an existing workspace

This is a naming change: update custom configuration switches and scripts that
refer to old paths. Old `run_amas_*` settings are no longer recognized. Migration
of stored output is an explicit offline operation. Stop workflows using the
output root before running it, and use the updated GeneGalleon container.

From the repository root, for orthogroups:

```bash
python workflow/migrations/migrate_alignment_statistics.py \
  --root workspace/output/orthogroup \
  --workspace-root workspace \
  --mode orthogroup \
  --genecount workspace/output/orthofinder/Orthogroups_filtered/Orthogroups.GeneCount.selected.tsv \
  --summary-out workspace/output/orthogroup_summary.tsv
```

For query2family, use `--mode query2family`, its output root, and
`--dir-query-gene workspace/input/query_gene` instead of `--genecount`.
Choose `--summary-out` to match the aggregate table used by your workspace.

The migration:

1. Reads each old statistics manifest to identify the actual selected alignment
   and DNA/protein mode. It reads live files or ZIP members without rerunning
   alignment, trimming, or tree inference.
2. Generates and validates all replacement statistics before publishing any.
3. Writes the new files and provenance, then deletes the old statistics and
   old statistics provenance using the output store's deletion mechanism.
4. Purges deleted ZIP members and regenerates the requested summary and its
   statistics-augmented gene-count table. ZIP purge rewrites archives and can
   require additional disk space.

If a historical manifest is absent, supply `--alignments mappings.tsv` with
columns `family_id`, `stage`, `alignment`, and `seq_type`. Stage is `original`
or `cleaned`; sequence type is `dna` or `aa`. Relative alignment paths are
relative to `--root`. Paths inside the output root can refer to archived files.
The migration deliberately does not guess whether cleaned input came from
MaxAlign, trimAl, or ClipKIT. Missing or invalid input stops migration before
old statistics are removed. Rerunning the command completes an interrupted
migration or refreshes the aggregate summary.

Only this one-time migration and its fixtures need to know historical names.
Known old aggregate filenames beside the requested summary and source gene-count
table are removed after successful summary generation. Other copies and custom
configuration files are not automatically discovered or removed.

## Validation

The naming migration was validated in a Docker-backed GeneGalleon runtime with
CDSKIT 0.31.0: 417 focused tests passed. Coverage includes both workflow modes,
DNA/protein input, raw/ZIP storage, both historical filename layouts, explicit
alignment mappings, repeat migration, invalid-input preservation, and refreshed
aggregate values. Related output-store and shell-contract tests also passed.
Ruff and Bash syntax checks passed. SIF execution and a complete multi-platform
container rebuild were not performed.
