# Public plant trait databases

Input generation supports anonymous retrieval from BROT 2.0 (`brot`), China
Plant Trait Database v2 (`cpt`), AlgaeTraits/WoRMS (`algaetraits`), and the BIEN
CSV species download service (`bien`). GIFT and AusTraits retain their existing adapters.
TRY remains **local export input only**; it is not an anonymous download adapter.
EOL's registered generic `species_api` entry is not a configured TraitBank client.

Use the bundled source map and select an explicit trait plan, for example:

```bash
export GG_INPUT_RUN_GENERATE_SPECIES_TRAIT=1
export GG_INPUT_TRAIT_PROFILE=none
export GG_INPUT_TRAIT_PLAN=workspace/input/input_generation/trait_plan_public_plants.tsv
export GG_INPUT_TRAIT_DATABASES=brot,cpt,algaetraits,bien
bash workflow/gg_input_generation_entrypoint.sh
```

The example is deliberately small. Every retrieved observation is also retained
in `<trait_download_dir>/<database>/runs/<run>.tsv`, with `species`, `trait_key`,
`value`, `unit`, `source_id`, `source_reference`, and original `context`. Inspect
its distinct keys to extend the plan; use `source_column=value` and
`trait_key_column=trait_key`. Select `numeric/median`, `categorical/mode`, or
`text/unique` according to the source definition. Never infer the type from a
numeric-looking category. Database-prefixed output names prevent accidental
cross-database pooling. Missing observations remain missing.

## Acquisition and interpretation

| ID | Default anonymous source | Interpretation |
| --- | --- | --- |
| `brot` | [BROT 2.0 data CSV](https://figshare.com/articles/dataset/BROT_plant_functional_trait_database_Data_file/5280868) | Key is `Trait\|DataType\|Units`; preserves quantitative, categorical, boolean and conditional representations separately. SourceID links to the release bibliography. |
| `cpt` | [CPT v2 release](https://figshare.com/articles/dataset/The_China_Plant_Trait_Database_Version_2_0/19448219) | Sample traits join Species translations by SAMPLE ID; Photo Pathway joins Taxonomic standardisation by SPECIES ID. Keys are `table:column`. Units require the [published dictionary](https://www.nature.com/articles/s41597-022-01884-4); `see_CPT_v2_dictionary` is not a physical unit. |
| `algaetraits` | [WoRMS REST](https://www.marinespecies.org/rest/) | Exact accepted species only. Excludes inherited attributes from other taxon IDs. Keys are `measurementTypeID\|measurementType`. Attributes with nested qualifiers retain their value and qualifiers as JSON and must be selected as text. |
| `bien` | `https://mint-pheasant.nceas.ucsb.edu:5775/api/download/traits?species=...` | Uses the [official backend CSV contract](https://github.com/EnquistLab/Biendata-Backend-Express/blob/main/controllers/traitDownloadController.js), preserving units in `trait_name\|unit` keys and source URLs (full original fields remain in context). Non-public records are excluded. |

BROT `RespFire` can be either a yes/no observation or a percentage. These must
not share one numerical column. BROT conditional values and ranges are preserved
as published; request them as text or curate a separate conversion.
CPT rows with a nonempty outlier flag (except `0`/`NA`) are conservatively excluded
as entire records; raw files retain them. Climatic/site covariates are not exposed
as species traits. Qualified names, unidentified species and subspecies are not
silently reduced to binomials. CPT names are the release's accepted names;
no additional synonym inference occurs. AlgaeTraits also contains non-green
algae: the genome manifest remains responsible for the target taxonomic scope.

The source-map acquisition mode is `public_plant_traits`, with
`species_column=species` and `trait_key_column=trait_key`. An optional `uri`
overrides the BROT CSV URL, CPT file-ID base URL, WoRMS REST base URL, or BIEN
species download endpoint respectively. No login, tokens, or account creation are used.

## Availability, provenance and failures

On 2026-09-14 anonymous BROT, CPT, AlgaeTraits and BIEN CSV downloads were
verified from the GeneGalleon Docker runtime. BIEN's separate JSON route
`/api/traits/species` was unavailable (HTML 404); the supported default uses its
working official `/api/download/traits` CSV route. An explicitly selected JSON
service can use `response_format=json`. A missing route is never interpreted as
missing biological observations. BIEN requests are spaced at least 9.1 seconds
apart within an acquisition to respect the documented 100 requests/15 minutes;
run one BIEN acquisition at a time per public IP.

Each acquisition writes content-addressed raw responses and a uniquely named JSON
receipt with URLs, retrieval times, SHA-256 hashes, completion/failure status and
the normalized table hash. Files are installed atomically. An error does not
publish a normalized success table. A documented BIEN JSON no-data 404 is an
empty result; HTML 404, malformed payloads, missing columns and unsafe joins fail.
Requests use the workflow's guarded network transport. HTTP service limits and
failures remain visible. Raw downloads are evidence, not repository fixtures.
Acquisitions currently refetch sources; the normal successful workflow artifact
cache still skips an unchanged completed stage. Use a fresh output workspace for
new snapshots and retain the source receipts with the final dataset.

Use `strict=1` for production to prevent source failures becoming partial trait
outputs. Raw data and normalized observations retain source identifiers; merging
across databases still requires explicit unit/trait harmonization and original
study deduplication. The adapter does not claim that DB records are independent
measurements or that a genome accession and trait observation share a genotype.

## Validation

Validated on 2026-09-14 using `local/genegalleon:dev` with the current workspace
mounted into the Docker container:

- 67 adapter/generator/template/GIFT/schema tests passed.
- 29 input-generation core, end-to-end and species-array tests passed.
- Live strict generator run: BROT, CPT and AlgaeTraits produced a five-species
  typed trait table; all five species had at least one selected observation.
- Live strict BIEN generator run: Quercus robur returned 704 public observations
  across 34 unit-qualified keys; the selected wood-density column was generated.

These are Docker checks, not SIF or audrey1 deployment validation. Live downloads
and receipts are retained outside the repository; no fetched datasets are tracked.
