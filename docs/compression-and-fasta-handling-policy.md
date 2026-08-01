# Compression and FASTA Handling Policy

GeneGalleon is intentionally gzip-friendly. Current practice is:

- inputs may be plain or gzipped,
- retained pipeline outputs are standardized toward `.fa.gz`,
- plain FASTA is generated only when an intermediate third-party step requires it,
- query2family and orthogroup artifacts may be grouped into standard
  ZIP shards to reduce filesystem inode usage.

Gzip and ZIP serve different purposes here. Gzip remains the retained format
for an individual FASTA. ZIP-backed gene-family storage groups many otherwise
independent artifacts and presents them through a logical file view. Already
compressed members such as `.fa.gz`, `.pdf`, and nested `.zip` files are stored
without redundant recompression.

Set `GG_COMMON_GENE_FAMILY_OUTPUT_STORAGE=files` (or its `raw` alias) in
`workflow/gg_common_params.sh` to retain the historical physical-file layout,
or leave the experimental branch default `zip`. See
[Gene-Family Outputs and Progress Monitoring](gene-family-outputs-and-progress-monitoring.md#zip-backed-storage)
for archive scope and manual-management commands.

## Accepted input forms

Across major stages, input discovery supports both plain and gzipped FASTA extensions such as:

- `.fa`, `.fas`, `.fasta`, `.fna`
- `.fa.gz`, `.fas.gz`, `.fasta.gz`, `.fna.gz`

This applies to the main CDS/query discovery paths described in `docs/input-conventions.md`.

## Retained output convention

Pipeline-tracked FASTA outputs, especially `file_*` outputs in `workflow/core/gg_*_core.sh`, are standardized to:

- `.fa.gz`

This convention is used to:

- reduce storage footprint,
- avoid ambiguity between retained artifacts and temporary plain FASTA files,
- make downstream scripts prefer one predictable naming convention.

## Where plain FASTA still appears

Plain FASTA may still be created temporarily when:

- a tool does not read gzip-compressed input directly,
- a conversion step is used as a short-lived staging file,
- a plotting or reconciliation helper expects an uncompressed intermediate.

Those plain FASTA files are not the preferred retained format for long-lived outputs.

## Stage-specific implications

Current noteworthy behavior:

- species-tree intermediates retain primary FASTA/alignment artifacts as `.fa.gz`,
- gene-evolution intermediates follow the same retained-artifact convention,
- input-formatting helpers also normalize formatted CDS outputs to gzipped `.fa.gz`.

Legacy plain `.fasta` artifacts in some species-tree intermediate locations are auto-migrated by current code paths when those stages are rerun.

## Operational guidance

Recommended practice:

- keep your curated input datasets gzipped unless a provider/tool specifically requires plain text,
- do not rename retained `.fa.gz` outputs to plain `.fasta` just for convenience,
- inspect compressed FASTA with `gzip -cd` or tools that read gzip transparently.

Examples:

```bash
gzip -cd workspace/input/species_cds/Arabidopsis_thaliana.fa.gz | head
```

```bash
gzip -cd workspace/output/species_tree/some_alignment.fa.gz | seqkit stats
```
