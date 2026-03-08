# Input Conventions

### Species naming

Species IDs are inferred from filename prefixes such as `Genus_species`.
Many joins/merges across pipeline stages rely on this convention.

### FASTA file extensions

Major scripts accept:

- `.fa`, `.fas`, `.fasta`, `.fna`
- `.fa.gz`, `.fas.gz`, `.fasta.gz`, `.fna.gz`

### `workspace/input/species_cds`

- one CDS FASTA per species,
- sequence IDs should be unique,
- avoid spaces and `|` in IDs when possible,
- species prefix on IDs (`Genus_species_...`) is strongly recommended.

### `workspace/input/query_gene`

Each file is one family-level task in `mode_query2family=1`.
Accepted forms:

- amino-acid FASTA,
- in-frame CDS FASTA,
- plain gene ID list (one ID per line).

`gg_gene_evolution_core.sh` auto-detects query type:

- FASTA with DNA sequences: translated to AA query,
- FASTA with protein sequences: used directly,
- non-FASTA text: treated as gene IDs to extract CDS then translate.

### Transcriptome assembly input modes

`gg_transcriptome_generation_core.sh` supports three mutually exclusive modes:

- `mode_sraid=1`
  - input: `workspace/input/query_sra_id/GENUS_SPECIES.txt`
  - one SRA/BioProject ID per line
- `mode_fastq=1`
  - input: `workspace/input/species_rnaseq/GENUS_SPECIES/*.fastq.gz`
- `mode_metadata=1`
  - input: `workspace/input/amalgkit_metadata/GENUS_SPECIES_metadata.tsv`

For `mode_sraid=1` and `mode_fastq=1`, auto-generated amalgkit metadata is written to:

- `workspace/output/transcriptome_assembly/amalgkit_metadata/GENUS_SPECIES_metadata.tsv`

This keeps `workspace/input/amalgkit_metadata` reserved for explicit `mode_metadata=1` inputs.

### Automated provider formatting helper

Manual formatting can be replaced with `workflow/support/format_species_inputs.py`, which ports key rules from legacy raw genome formatting scripts such as:

- `data_formatting_ensemblplants.sh`
- `data_formatting_phycocosm.sh`
- `data_formatting_phytozome.sh`

Example (single provider, explicit input directory):

```bash
python workflow/support/format_species_inputs.py \
  --provider ensemblplants \
  --input-dir "/path/to/raw_genomes/20230216_EnsemblPlants/original_files" \
  --species-cds-dir workspace/output/input_generation/species_cds \
  --species-gff-dir workspace/output/input_generation/species_gff
```

Example (all providers, provider-root input directory):

```bash
python workflow/support/format_species_inputs.py \
  --provider all \
  --input-dir "/path/to/raw_genomes" \
  --species-cds-dir workspace/output/input_generation/species_cds \
  --species-gff-dir workspace/output/input_generation/species_gff
```

Notes:

- output filenames are normalized to start with `Genus_species_...`,
- formatted CDS outputs are always gzipped with `.fa.gz` extension,
- formatted GFF outputs are always gzipped with `.gff.gz` extension,
- CDS IDs are prefixed with `Genus_species_...` and aggregated to one representative CDS per gene,
- common historical replacements are applied to CDS/GFF text,
- CDS are padded to codon-length multiples and transcript-level redundancies are collapsed at gene level.

Download-first workflow (manifest driven):

```bash
python workflow/support/format_species_inputs.py \
  --provider all \
  --download-manifest /path/to/download_plan.xlsx \
  --download-dir workspace/output/input_generation/tmp/input_download_cache \
  --species-cds-dir workspace/output/input_generation/species_cds \
  --species-gff-dir workspace/output/input_generation/species_gff
```

Manifest required columns:

- `provider` (required; first column in XLSX templates)
- `id` (required; second column in XLSX templates)
- `species_key` is optional
- `cds_url` and `gff_url` (or `id` to auto-resolve provider-specific URLs when supported)
  - for `provider=ncbi`, when `species_key` is omitted and `id` is given, `species_key` is inferred from NCBI species metadata (e.g. `Homo_sapiens`).
  - `provider=ncbi` accepts both `GCF_*` and `GCA_*` assembly accessions and auto-resolves NCBI assembly URLs.
  - for `provider=coge`, `id` must be CoGe `genome_id` (numeric `gid`), and CDS/GFF/Genome URLs are auto-built.
  - for `provider=cngb`, built-in inference resolves CNGB assembly IDs (`CNA...`, `cngb:...`) or linked `GCA/GCF` accessions and maps to downloadable assembly files.
  - for `provider=flybase`, `provider=wormbase`, `provider=vectorbase`, and `provider=fernbase`, `id` can be resolved via explicit URL columns or `GG_<PROVIDER>_*_URL_TEMPLATE`.
  - for `provider=ensembl` and `provider=ensemblplants`, `id`-only URL inference is supported via:
    - provider defaults (for example, Ensembl/EnsemblPlants index discovery),
    - or env templates: `GG_<PROVIDER>_CDS_URL_TEMPLATE`, `GG_<PROVIDER>_GFF_URL_TEMPLATE`, `GG_<PROVIDER>_GENOME_URL_TEMPLATE`,
    - and optional page template: `GG_<PROVIDER>_ID_URL_TEMPLATE`.
  - for `provider=fernbase`, `id`-only URL inference can also resolve from FernBase FTP directory discovery.
  - for `provider=local`, `id` can point to a local species directory (or set explicit local file paths) and formatting starts from local files.
  - `provider=phycocosm` and `provider=phytozome` are not supported in `--download-manifest`.
    Use `--input-dir` for local raw-file formatting only.

Optional columns:

- `species_key`
- `cds_filename`
- `gff_filename`
- `genome_filename`
- `cds_url_template`, `gff_url_template`, `genome_url_template`
  (`{id}`, `{species_key}`, `{provider}` placeholders)
- `local_cds_path`, `local_gff_path`, `local_genome_path`
  (for `provider=local`)

XLSX template notes:

- `download_plan.xlsx` includes drop-downs for `provider` and `id`.
- `id` drop-down values are provider-specific.
- provider drop-down order is fixed, and `local` is always listed last.
- for large provider (`ncbi`), five model-organism IDs are shown as examples (mixed `GCF_*`/`GCA_*` formats).
- for `coge` and `cngb`, IDs are also example-based by default (both show five model-organism examples).
- for `ensembl`, `ensemblplants`, `flybase`, `wormbase`, `vectorbase`, `fernbase`, and `local`,
  IDs can be supplied from a prebuilt `id_options_snapshot.json`.
- when no snapshot is supplied, non-large providers fall back to IDs discovered from `--input-dir`.
- drop-down IDs are shown as `ID (Species name)` for non-`local` providers.
- for `provider=local`, drop-down IDs are plain local directory/path-style IDs.
- at runtime, `id` parsing uses the token before the first space (for non-`local` providers), so labels like `GCF_000001405.40 (Homo sapiens)` are accepted.
- any other value can still be typed manually.

Use `--download-only` to fetch raw files and stop before formatting.
Use `--dry-run` to preview downloads and formatting outputs without writing files.
Use `--jobs` to set download parallelism (defaults to `GG_TASK_CPUS`, fallback `1`).
Resolved manifest rows (for example URL/template/local-path auto-resolution and filename filling)
are written to:
- `workspace/output/input_generation/download_plan.resolved.tsv`

Download concurrency safeguards:

- per-provider caps are applied even when `--jobs` is large.
- override per-provider cap with:
  `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_<PROVIDER>` (e.g. `..._NCBI=2`).
- NCBI E-utilities calls are throttled to documented limits by default:
  - `3 req/s` without API key,
  - `10 req/s` with `GG_NCBI_API_KEY`.
- override NCBI E-utilities rate with `GG_NCBI_EUTILS_MAX_RPS`.

Optional download authentication:

- `--auth-bearer-token-env ENV_NAME`
  - adds `Authorization: Bearer <token>` from an environment variable.
- `--http-header "Key: Value"` (repeatable)
  - add arbitrary request headers (for provider-specific requirements).

Build a manifest automatically from an existing local dataset tree:

```bash
python workflow/support/build_download_manifest.py \
  --provider all \
  --input-dir "/path/to/raw_genomes" \
  --id-options-snapshot /path/to/id_options_snapshot.json \
  --output workspace/input/input_generation/download_plan.xlsx
```

This writes `file://` URLs, so you can test the full download+format pipeline locally before replacing them with remote URLs.

Use `workspace/input/input_generation/download_plan.xlsx` as the runtime manifest input file.
Resolved/filled rows are written to `workspace/output/input_generation/download_plan.resolved.tsv`.

Latest template distribution:

- On each push to the default branch, GitHub Actions runs `build_download_manifest.py`
  and publishes the latest `download_plan.xlsx`.
- The workflow first builds `id_options_snapshot.json` (remote provider ID choices),
  then generates `download_plan.xlsx` from that snapshot.
- If remote ID fetch fails for some providers, the workflow reuses provider entries
  from the previous `id_options_snapshot.json` release asset.
- Release asset URL (rolling latest):
  `https://github.com/<owner>/<repo>/releases/download/download-plan-latest/download_plan.xlsx`
- Snapshot asset URL (rolling latest):
  `https://github.com/<owner>/<repo>/releases/download/download-plan-latest/id_options_snapshot.json`
- Workflow definition:
  `.github/workflows/download-plan.yml`

`gg_input_generation_entrypoint.sh` wrapper:

- runs `gg_input_generation_core.sh` inside the container as a single entrypoint,
- can download provider files from a manifest and format inputs in one run,
- can validate produced `species_cds` naming, `species_gff` consistency, and CDS-to-GFF mapping compatibility,
- CDS-to-GFF mapping validation is species-parallel and accepts `validate_cds_gff_mapping.py --nthreads N` (`--ncpu` remains as a compatibility alias),
- can optionally generate `workspace/input/species_trait/species_trait.tsv`
  from configured trait databases.

Configuration:

- edit the `### Start: Modify this block ... ###` section in
  `workflow/gg_input_generation_entrypoint.sh`,
- keep `workspace/input/input_generation/download_plan.xlsx` as the runtime manifest input file.
- check resolved rows at `workspace/output/input_generation/download_plan.resolved.tsv`.

Alternative runtime overrides (without editing files) via env vars:

- `GG_INPUT_PROVIDER`, `GG_INPUT_STRICT`, `GG_INPUT_OVERWRITE`,
- `GG_INPUT_DOWNLOAD_ONLY`, `GG_INPUT_DRY_RUN`,
- `GG_INPUT_DOWNLOAD_TIMEOUT`,
- `GG_INPUT_DOWNLOAD_MANIFEST`, `GG_INPUT_INPUT_DIR`, `GG_INPUT_DOWNLOAD_DIR`,
- `GG_INPUT_SPECIES_CDS_DIR`, `GG_INPUT_SPECIES_GFF_DIR`, `GG_INPUT_SPECIES_GENOME_DIR`,
- `GG_INPUT_SPECIES_SUMMARY_OUTPUT`,
- `GG_INPUT_RESOLVED_MANIFEST_OUTPUT`,
- `GG_INPUT_AUTH_BEARER_TOKEN_ENV`, `GG_INPUT_HTTP_HEADER`,
- `GG_INPUT_SUMMARY_OUTPUT`,
- `GG_INPUT_RUN_FORMAT_INPUTS`, `GG_INPUT_RUN_VALIDATE_INPUTS`,
- trait generation:
  `GG_INPUT_TRAIT_PROFILE` (`none` or `gift_starter`; default `none`),
  `GG_INPUT_RUN_GENERATE_SPECIES_TRAIT`,
  `GG_INPUT_TRAIT_SPECIES_SOURCE` (`download_manifest` or `species_cds`; default `download_manifest`),
  `GG_INPUT_SPECIES_TRAIT_OUTPUT`,
  `GG_INPUT_TRAIT_PLAN`,
  `GG_INPUT_TRAIT_DATABASE_SOURCES`,
  `GG_INPUT_TRAIT_DATABASES`,
  `GG_INPUT_TRAIT_DOWNLOAD_DIR`,
  `GG_INPUT_TRAIT_DOWNLOAD_TIMEOUT`.
- per-provider download caps:
  `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_ENSEMBL`,
  `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_ENSEMBLPLANTS`,
  `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_PHYCOCOSM`,
  `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_PHYTOZOME`,
  `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_NCBI`,
  `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_REFSEQ`,
  `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_GENBANK`,
  `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_COGE`,
  `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_CNGB`,
  `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_FLYBASE`,
  `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_WORMBASE`,
  `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_VECTORBASE`.

Quick preset example (enable trait stage with GIFT starter config):

```bash
GG_INPUT_TRAIT_PROFILE=gift_starter bash workflow/gg_input_generation_entrypoint.sh
```
