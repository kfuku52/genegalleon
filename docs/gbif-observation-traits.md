# GBIF observation traits

The `gbif_distribution` preset describes **retained GBIF presence records**.
Its extrema, means and occupied-cell area are not unbiased estimates of a
species' biological range. Acquisition completeness means completeness of a
query/download, not complete observation of the species in nature.

## Generate a bundle

Run entrypoints with a GeneGalleon container runtime as usual:

```bash
GG_INPUT_TRAIT_PROFILE=gbif_distribution \
bash workflow/gg_input_generation_entrypoint.sh
```

No GBIF account is needed for occurrence search. Each page contains at most 300
records and each species search retrieves at most 100,000 records. Records are
read in API page order; a capped acquisition is not a random or spatially
stratified sample. For larger datasets, use an existing GBIF download.
[GBIF occurrence API](https://techdocs.gbif.org/en/openapi/v1/occurrence).

For `species_trait.tsv`, the generator publishes one transaction containing:

| File | Meaning |
|---|---|
| `species_trait.tsv` | Ordinary traits and renamed GBIF observation metrics. Unresolved/incomplete GBIF values are missing |
| `species_trait.tsv.metadata.json` | Hash-bound table contract, trait roles/units/definitions, effective filters, acquisition and source identities, geometry, species-level quality and exclusion reasons |
| `species_trait.tsv.gbif-quality.tsv` | Counts, acquisition status, taxon match, date/missingness/dataset summaries; not a trait matrix |
| `species_trait.tsv.gbif-observations.tsv` | Descriptive summaries including available partial acquisitions; not an analysis-ready trait matrix |

The sidecars are also emitted for non-GBIF runs, with empty GBIF tables. Keep the
metadata next to the trait file when moving it. A modified table with stale
metadata is rejected rather than silently interpreted using the old contract.
For a documented subset, use `species_trait_contract.py` (below), which writes
a new contract bound to the selected table.

Normalized records and match/acquisition information are saved in immutable
`.records.jsonl.gz` snapshots under the trait download cache. The metadata
records the snapshot's path and SHA-256. These records retain the fields needed
for refiltering, dates, source attribution and observation-process audits; they
are not verbatim copies of every API response field. Keep snapshots for replay.
Move a snapshot only with an explicit update of its recorded path; its hash must
still match. A live paged search is not an atomic GBIF database snapshot.

## Definitions and migration

| Previous column | New location/name |
|---|---|
| `gbif_occurrence_count` | Quality `reported_count`: API query count before local filtering, or rows for the mapped species in a local file |
| `gbif_occurrence_used` | Quality `retained_count`, separately from `raw_fetched`, `unique_gbif_ids` and exclusions |
| `gbif_occurrence_truncated` | Quality `status` and `termination_reason`, covering caps, interrupted/invalid pages and inconsistent acquisition |
| `gbif_northern_limit_lat`, `gbif_southern_limit_lat` | `gbif_observed_northern_limit_lat`, `gbif_observed_southern_limit_lat` |
| `gbif_latitudinal_breadth_deg` | `gbif_observed_latitudinal_breadth_deg` |
| `gbif_western_limit_lon`, `gbif_eastern_limit_lon` | `gbif_observed_western_limit_lon`, `gbif_observed_eastern_limit_lon` |
| `gbif_longitudinal_breadth_deg` | `gbif_observed_longitudinal_breadth_deg` |
| `gbif_occupied_grid_area_km2` | `gbif_observed_occupied_grid_area_km2` |
| `gbif_centroid_lat` | `gbif_observed_record_mean_lat` |
| `gbif_centroid_lon` | `gbif_observed_record_circular_mean_lon` |
| `gbif_country_count` | `gbif_observed_country_count` |
| `gbif_convex_hull_area_km2` | Removed: the previous equirectangular plane approximation had no validated global/geodetic interpretation |

Regenerate old GBIF outputs and update explicit trait plans. Old quality columns
cannot be exported as traits. Custom aliases for new metrics keep their
observation role in metadata and cannot share a column with another source.
Metrics remain numeric; this preset does not define binary foregrounds.

The latitude mean uses one weight per retained GBIF record ID. Separate events
at one coordinate can therefore change the mean. It is not an area centroid or
an abundance-weighted population centre. The longitude mean is circular and is
missing for a near-zero resultant (threshold `1e-12`). The shortest longitude
arc handles the antimeridian; its starting longitude can exceed its ending
longitude in the usual numeric ordering. Tied shortest arcs use a deterministic
ordering of normalized longitudes, not a biological choice of eastern/western
boundaries.
Records exactly at the poles do not contribute to longitude summaries because
their longitude is undefined. If all retained records are polar, longitude
limits, breadth and mean are missing, with an explicit quality reason.

The default grid is 1 degree, anchored at latitude −90 and longitude −180.
Longitudes use `[-180, 180)`, so ±180 is one location. The grid includes its
lower edges and excludes upper edges, except that the poles are assigned to
the end latitude cells with canonical longitude −180. End cells are clipped
at +90/+180. Cell areas are summed once per occupied cell using
`R² × longitude_width_radians × (sin(latitude_max) − sin(latitude_min))`,
with `R = 6371.0088 km`. The whole cell contributes, including any sea or
unobserved habitat in that cell. There is no land/habitat mask. Valid grid
widths are finite numbers in `(0, 180]` degrees.

This is neither equal-area gridding nor IUCN AOO. IUCN's 2×2 km scale and its
biological/seasonal scope requirements are separate from this preset. Changing
the grid size alone cannot turn raw GBIF records into a Red List assessment.
[IUCN Guidelines §4.10](https://nc.iucnredlist.org/redlist/content/attachment_files/RedListGuidelines.pdf).

## Acquisition and filtering

Species search must resolve to a species-level match with a usable taxon key,
EXACT/FUZZY match type and confidence at least 90 by default. Missing confidence
or higher-rank matches are unresolved, not automatically accepted. Local
downloads use a reviewed species-level mapping supplied by the user.
Different target species that resolve to the same GBIF speciesKey are
ineligible for analysis (`shared_species_key`); a local map with such duplicate
keys is rejected. Synonyms must not become independent species replicates.

Search requests restrict records to PRESENT, coordinates and no GBIF
geospatial issue. Date/country/basis/establishment/uncertainty filters are then
applied locally, preserving the original acquisition for sensitivity analyses.
Local checks also reject failed/suspicious coordinate reprojections; the exact
excluded flag set is recorded in the metadata. Successfully reprojected
coordinates remain eligible. [GBIF issue definitions](https://techdocs.gbif.org/en/data-use/occurrence-issues-and-flags).
Thus the search cap applies **before** these local filters. Duplicate gbifIDs
are counted and removed; identical coordinates with different IDs are retained.
Duplicate or missing IDs during paging invalidate a complete-search claim.

| Status | Interpretation |
|---|---|
| `complete_search` | Count, page positions, unique IDs and end flags agree for the requested search |
| `complete_download` | Whole local file row count agrees with saved successful GBIF download metadata |
| `complete_local_file` | File read successfully, but download completeness has not been verified; descriptive output only |
| `capped_partial` | Requested species had more records than the acquisition cap |
| `interrupted` | Request failure or early empty/short page |
| `invalid_response` | Invalid count, page shape/position, or missing/mismatched speciesKey |
| `inconsistent_search` | Count, record IDs or end flags conflict during paging |
| `taxon_unresolved` | Species-level match was not accepted |

Incomplete/unresolved acquisitions are audited but not reused as complete
caches. A later workflow run retries an ineligible acquisition even if the
earlier run produced an all-missing trait table. `gbif_require_complete=yes` fails the run, preserving the acquisition
audit without replacing the previous trait bundle. The default `no` publishes
missing values for those species in the main trait table. Zero retained records
give missing observation metrics; they do not establish biological absence.
Check exclusions and unknown counts even for a complete acquisition.

Use these settings in the entrypoint, or uppercase them with the `GG_INPUT_`
prefix. Corresponding Python CLI flags replace underscores with hyphens.
They can also be columns in the GBIF row of `trait_database_sources.tsv`.
An explicitly supplied CLI value takes precedence over that row, including
values equal to the standard defaults.

| Setting | Default/meaning |
|---|---|
| `gbif_api` | GBIF v1 API base; no credentials/query string in the base URI |
| `gbif_page_size` | 300 maximum |
| `gbif_max_occurrences_per_species` | 100,000 maximum for search; not a sample-size standardization |
| `gbif_grid_degrees` | 1 |
| `gbif_min_match_confidence` | 90 |
| `gbif_year_min`, `gbif_year_max` | No period restriction; when set, the entire reported event-year interval must fit |
| `gbif_countries` | No country restriction; comma-separated country codes |
| `gbif_include_basis_of_record`, `gbif_exclude_basis_of_record` | No basis restriction; explicit comma-separated values |
| `gbif_include_establishment_means` | No origin restriction; when set, only specified values pass, not unknown |
| `gbif_max_coordinate_uncertainty_m` | Disabled |
| `gbif_min_distance_from_known_centroid_m` | Disabled; reject points nearer than the threshold to a known georeferencing centroid |
| `gbif_missing_date` | `exclude` with an active year window; optionally `keep` |
| `gbif_missing_uncertainty`, `gbif_missing_centroid_distance` | `keep` with an active threshold; optionally `exclude` |
| `gbif_use_cache` | `yes`; `no` requests a fresh acquisition and invalidates the workflow's trait-stage cache |
| `gbif_require_complete` | `no` |
| `gbif_occurrence_file`, `gbif_taxon_map`, `gbif_download_metadata` | Empty; local-download inputs described below |

`gbif_max_distance_from_centroid_m` and
`GG_INPUT_GBIF_MAX_DISTANCE_FROM_CENTROID_M` are removed and must not be carried
over unchanged. GBIF's field measures distance from a known georeferencing
centroid, not the species' observed centre. The new **minimum** distance
excludes proximity to those points. It is not an outlier filter for range
edges. Missing distances do not prove good coordinates; their handling is
explicit. [GBIF field definitions](https://techdocs.gbif.org/en/data-use/download-formats),
[GBIF filtering guide](https://data-blog.gbif.org/post/gbif-filtering-guide/).

For contemporary wild occurrences, decide whether introduced populations are
part of the target, and choose period and basis filters accordingly. No default
date, native-only filter or universal fossil exclusion is silently imposed.
Unknown establishment is never relabelled native. Excluding LIVING_SPECIMEN
alone does not remove all cultivation. Observation counts are not measured
effort; event/protocol/effort fields and their unknown counts are retained for
auditing, not collapsed into a supposedly comparable effort number.
[GBIF sampling-event documentation](https://techdocs.gbif.org/en/data-publishing/data-quality-recommendations).

## Import an existing download without network access

Provide the entire interpreted SIMPLE_CSV download table (GBIF's `.csv` is
tab-delimited), its gzip equivalent, or a ZIP containing exactly one `.csv`/`.tsv`
table. For other archive formats, extract the interpreted table explicitly.
Required headers include gbifID, speciesKey, decimalLatitude, decimalLongitude,
occurrenceStatus and issue(s). Do not use a species-list or a verbatim table.

The reviewed mapping must contain:

```tsv
species	taxon_key	scientific_name
Arabidopsis_thaliana	3052436	Arabidopsis thaliana
```

Verify keys/names against the taxonomy of your download; the example is only a
format illustration. `species` must match the requested GeneGalleon labels.
Species not in the map cause an error. Records are selected by speciesKey so
infraspecific records within the mapped species are included.

Save the official response from `/occurrence/download/{key}` as JSON. Its
SUCCEEDED status, download key, DOI, request and totalRecords establish the
download scope and row-count check. A table without this evidence remains
`complete_local_file`; a count mismatch is an error. Matching counts are an
integrity check, not proof that a publisher sampled the true range completely.
Local files, mappings and metadata are content-hashed, including when their
paths are supplied indirectly through the database-source table. Array preparation also fingerprints these source files, so changing them requires a new preparation.

```bash
GG_INPUT_TRAIT_PROFILE=gbif_distribution \
GG_INPUT_GBIF_OCCURRENCE_FILE=/data/gbif/download.zip \
GG_INPUT_GBIF_TAXON_MAP=/data/gbif/reviewed_taxa.tsv \
GG_INPUT_GBIF_DOWNLOAD_METADATA=/data/gbif/download.json \
GG_INPUT_GBIF_REQUIRE_COMPLETE=yes \
bash workflow/gg_input_generation_entrypoint.sh
```

The GBIF acquisition mode requires `database=gbif`; custom database names cannot strip the observation/quality contract. Custom output column aliases remain supported.

This import never submits a new download or makes a taxonomy API request. A
future authenticated download-submission feature would be a separate adapter.

## Use in related analyses

`trait=all` in copy-number PGLS / nested-CV selection and `rsc_predictors=all` exclude observation and
quality roles, including custom aliases. Explicitly select an observation
metric to analyze it. The metadata is checked and ineligible species are
masked with recorded reasons. Legacy GBIF columns without a contract cannot be
selected explicitly. Copy-number results include the selected table, its
metadata and `species_trait_input.json`; prepared RSC and species PGLS keep the
input interpretation in their metadata/audits. Nested-CV selection also saves the selected table and audit, and carries response meaning into each trait result metadata file. It remains exploratory, without post-selection inference. Quality columns are not valid
responses. Unrelated ordinary user trait tables do not require a sidecar.

The null hypothesis for such a regression concerns the selected **observed
record summary**, conditional on the model's specified predictors and
phylogenetic covariance. It is not a test of adaptation or an unbiased
biological-range estimate. Restricting to complete acquisitions can itself
select a nonrepresentative set of species. Study effort and genome/annotation
quality may confound observed associations; simply adjusting for record count
does not guarantee removal of that confounding.

GBIF observations are not automatically thresholded into codeml/HyPhy/CSUBST
foregrounds. Those stages use ordinary foreground traits, excluding GBIF
observations/quality; supplying only GBIF fields does not define a usable
foreground. When ordinary traits are binarized, missing/unrecognized values remain missing instead of becoming background 0. A biologically justified foreground definition is a separate
research decision.

For a portable, explicit subset, run inside the GeneGalleon runtime:

```bash
python workflow/support/species_trait_contract.py \
  --input workspace/input/species_trait/species_trait.tsv \
  --select gbif_observed_latitudinal_breadth_deg \
  --output /data/analysis/selected_traits.tsv \
  --report /data/analysis/trait_input.json
```

## Replay sensitivity analyses

This command uses only the saved, hash-verified records. It does not refetch the
API, correct for missing regions or provide confidence intervals for true
ranges. All comparisons remain within the original download/search universe.

```bash
python workflow/support/analyze_gbif_sensitivity.py \
  --metadata workspace/input/species_trait/species_trait.tsv.metadata.json \
  --grid-degrees 0.5,1,2 \
  --year-windows source,2000:2025 \
  --countries 'source;JP;US,CA' --missing-date source,keep,exclude \
  --establishment-means source,native \
  --selection records,unique_location,one_per_cell \
  --sample-sizes all,100,1000 \
  --replicates 20 --seed 123 \
  --output /data/analysis/gbif_sensitivity.tsv
```

Unique-location and one-per-cell selection choose an actual retained record
per group using a fixed seed. Fixed-size sampling is without replacement after
filtering/ID deduplication and any spatial selection. Insufficient sample sizes
are explicitly missing, not silently reduced. `--exclude-datasets` compares the
baseline with exclusion of each listed dataset in turn. The output identifies
every scenario, replicate, seed, acquisition status, available sampling units,
number of selected records and hash of the selected gbifID set, with a separate
provenance JSON. Country sets are separated with semicolons; codes within a
set use commas. Missing-date comparisons apply only with an active period.
Do not pick the most favourable
scenario's association P value as a confirmatory result.

The number of output rows grows with species × grid sizes × periods × origin
filters × regions × missing-date policies × dataset exclusions × selections × sample sizes × replicates (the
unsampled record baseline is emitted once). Measure time/RSS on a small subset
before expanding. Full unfiltered record storage is O(N), grid accumulation
O(N), and sorting for longitude intervals is O(N log N). Real processing costs
include record fields and sensitivity repetitions, not just coordinates.

Independent survey/effort data and study-specific validation are still needed
for a natural-distribution model. Target-group backgrounds and joint survey/
collection models are possible research approaches, not corrections implicitly
performed by this preset. [Phillips et al. 2009](https://pubmed.ncbi.nlm.nih.gov/19323182/),
[Fithian et al.](https://www.stat.berkeley.edu/~wfithian/biasCorrection.pdf).


The generator publishes the table, quality and descriptive sidecars, metadata,
and optional statistics JSON together, restoring earlier outputs if publication
fails. RSC preparation applies the same rollback contract to its prepared input
bundle. Outputs cannot alias source files or each other. Cache manifests carry
a digest of the full saved bundle as well as the record-snapshot digest; a
corrupt manifest or altered cached summaries trigger reacquisition.

### 形質型との接続

コピー数PGLSと同時選択は、観測値の利用制約と `.schema.json` の型検証を両方適用する。
`all` では通常のnumeric/binary列を選び、GBIF観測列・text/categorical列は理由付きで除外する。
観測列の明示選択でも、取得品質による種の欠測化と型検証を省略しない。
`trait_selection.tsv` に型と選択理由、`species_trait_input.json` に観測定義・品質・入力hashを残す。
解析用 `selected_species_traits.tsv` にもmetadataと型schemaを保存する。
