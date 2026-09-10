# Reliable GIFT trait retrieval

`generate_species_trait.py` uses validated, resumable GIFT requests during the
normal input-generation workflow. Its source map remains
`workspace/input/input_generation/trait_database_sources.tsv`.

## Cache and interrupted runs

The default `gift_cache_mode=reuse` saves each successful JSON response beneath
`<downloads-dir>/gift/pages/` within a query generation. Each cache record contains its exact URL, retrieval
time and a SHA-256 checksum of the canonical JSON payload. The URL includes the
API endpoint, concrete release, trait, bias settings, offset and limit. Species
selection, agreement thresholds and aggregation are applied again after loading
these raw pages, so changing those settings does not reuse a stale final table.

If a request fails, rerun the same command: completed pages are reused and the
missing page is fetched. Temporary writes are installed atomically; truncated,
corrupt and API-error responses cannot become valid cached pages. Transient
network errors, HTTP 408/429 and selected 5xx errors receive up to two retries by
default. Permanent request errors fail immediately. A server retry delay above
30 seconds ends the attempt so the operator can retry later.

The following optional source-map columns control retrieval:

| Column | Default | Meaning |
| --- | --- | --- |
| `gift_cache_mode` | `reuse` | `reuse`, `refresh`, or `offline` |
| `gift_retries` | `2` | Number of retries after transient failures, from 0 to 5 |
| `gift_species_mapping_file` | empty | Additional reviewed mapping TSV; relative paths are relative to the source map |
| `gift_page_size` | `10000` | Page size, 1–10000; sent explicitly as the API `limit` |
| `gift_max_pages_per_trait` | unlimited | Safety cap; reaching it before completion fails instead of publishing a partial trait |
| `gift_agreement_min` | unset | Optional finite threshold from 0 to 1 |

`refresh` fetches a fresh copy of every required page. `offline` reads validated
cache pages only and fails on a missing/corrupt page; it requires an explicit
stable `gift_version`. Online `latest` discovers the current release on each
invocation, then reuses that release's pages. Mutable `beta` pages are always
fetched online. Stable-release pages remain cached until explicitly refreshed.

The workflow's existing final-output cache can skip the whole trait stage.
After an interrupted attempt no successful output receipt is written. To force
an already successful stage to execute again, use the normal `overwrite=1`
entrypoint setting (or invoke the Python helper directly). `overwrite` reruns the
stage; `gift_cache_mode=refresh` additionally refreshes the API pages.

GIFT's [documented traits endpoint](https://biogeomacro.github.io/GIFT/articles/web_only/GIFT_API.html)
supports trait and page parameters, but no species-ID filter. A cold acquisition
therefore still downloads global trait pages. Repeated runs share those pages
across species sets and trait plans. No parallel flood of API requests is needed.

## Reviewed scientific names

Exact species matches must resolve to one unique GIFT ID. The client never
chooses a non-exact synonym or an ambiguous hit merely because it has the largest
match score. Such rows remain unresolved and their candidates appear in the
retrieval report.

The bundled [reviewed mappings](../workflow/support/gift_species_mappings.tsv)
cover four source-verified synonyms in the official public GIFT 3.2 release:
`Mimulus_guttatus`, `Populus_trichocarpa`, `Dendrobium_catenatum` and
`Physcomitrium_patens`. Their evidence URLs are stored alongside the decisions.
Each run queries the reviewed accepted name and checks both its name and expected
GIFT ID before using it. Mappings never apply automatically to a different
release. `Citrus_sinensis` is explicitly excluded because the available
`Citrus × aurantium` aggregate is a broader taxonomic concept. The old
`Dendrobium moniliforme` interpretation is not used for `D. catenatum`.

A custom mapping TSV uses these required columns:

```text
species  gift_version  work_ID  work_species  decision  evidence_url
```

Use tabs, not spaces, between columns. `species` is the input base species label
with underscores; `work_species` is the exact GIFT binomial with a space.
`decision` is `map` or `exclude`. A `map` requires a numeric ID and an HTTPS source
supporting the species equivalence. An exclusion documents a known unsafe
mapping. Each species/release pair must be unique within a file. Custom rows
override bundled rows for the same species/release pair. The source map, custom
mapping files, bundled decisions and retrieval code are hashed in the workflow's
input provenance. Changing a mapping invalidates the corresponding output cache.
The default stale-artifact policy stops before replacing that output; choose
`artifact_stale_policy=rebuild` to regenerate it, or use `overwrite=1`.
Multiple input labels for one resolved taxon retain their own output rows.

## Reports and missing data

Every real acquisition attempt writes a uniquely named JSON report in
`<downloads-dir>/gift/runs/`, and prints its path. It records completion/failure,
resolved version, reviewed mapping hashes, every species decision and candidate,
every page's network/cache source and checksum, and matched rows per trait.
A complete download can still have no observations for a requested trait or
species; that is distinct from a failed/incomplete acquisition.

Empty values and the explicit missing tokens `NA`, `NaN`, `N/A`, `null`, `none`
and `unknown` remain missing, including when mapping a positive category to a
binary value. Missing binary observations never become zero. With `--strict`,
a requested output trait with no observations fails before replacing the output.
The stats JSON includes `num_observed_by_trait`. Missing public observations are
not imputed from related species or broader taxa.

## Additional integrity checks

Paged trait caches use a separate generation for each endpoint/release/trait/
bias/page-size query. A refresh switches to a new generation before its first
request. Resuming that refresh cannot reuse pages from the older generation;
a reader already in progress pins its own generation. Old generations remain
available for such readers. A killed process leaves a `running` report and
validated completed pages; the next invocation resumes the active generation.
No atomic snapshot of a changing remote database is implied, especially for beta.

Trait names must match an exact `Trait2`, or a `Trait1` identifying one unique
trait ID. Ambiguous group names such as `Plant_height` fail with candidate IDs;
use `Plant_height_max` or its explicit ID. Substring/popularity matching is not
used. Specifying one trait by both ID and name does not duplicate observations.
Explicitly unresolved taxonomic hits and qualified/infraspecific target labels
are not silently promoted to species-level matches.

Plan and source TSVs reject duplicate/empty column names, extra fields, empty
required fields, duplicate database definitions and reserved output names.
Optional trailing fields may be omitted. Allowed type/aggregation pairs are:

| `value_type` | `aggregation` | Default |
| --- | --- | --- |
| `numeric` | `median`, `mean`, `min`, `max` | `median` |
| `binary` | `any`, `all`, `min`, `max`, `sum`, `mean` | `any` |
| `categorical` | `mode`, `first` | `mode` |
| `text` | `unique`, `first` | `unique` |

Text `unique` writes a sorted JSON array of distinct strings into the TSV cell;
embedded separators are preserved. Only actual null/empty text values are missing,
so a literal textual `NA` is retained. Numeric parsing rejects malformed or
infinite observations and overflowing aggregate results; it does not silently
remove them. Strict mode fails the run; non-strict mode warns and leaves that
trait missing. Numeric TSV output preserves round-trip floating-point precision.
Non-strict outputs retain requested trait columns even when a source fails or
has no matching rows, and observation counts then report zero.

Input/output aliases (including existing hard links), output symlinks and output
paths inside the GIFT cache or species-input directory are rejected. The trait
TSV, its `<output>.schema.json` type metadata, and optional stats JSON are
prepared before installation; ordinary
publication errors roll back previously installed files. A hard process kill
between file replacements is not a transactional guarantee; the workflow must
complete successfully before its output receipt is recorded.

The schema records each output column's declared `value_type` and the TSV's
SHA-256 fingerprint. Multiple source rows may fill the same output column only
when their declared types agree. Both copy-number analysis stages use it to
exclude text/categorical columns from `all` and publish `trait_selection.tsv`
with exclusion reasons; numeric columns remain strictly validated. The workflow
tracks the schema as an acquisition output and an analysis input, so schema
changes or removal invalidate the corresponding cached analysis. Preserve the
sidecar when copying a generated table. See
[copy-number models](copy-number-trait-models.md) for legacy table behavior.

The [integrity audit](gift-trait-retrieval-audit.md) records reproduced defects,
regression coverage and the full saved-dataset replay.
