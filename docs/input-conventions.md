# Input Conventions

### Species naming

Species IDs are inferred from filename prefixes such as `GENUS_SPECIES`.
Many joins/merges across pipeline stages rely on this convention.
`GENUS_SPECIES` is a placeholder for the actual species label, for example
`Arabidopsis_thaliana` or `Oryza_sativa`.

For most species-wise inputs, use either:

- a filename that starts with `GENUS_SPECIES_...`
- or a directory named `GENUS_SPECIES`

Typical examples:

- `Arabidopsis_thaliana_Araport11.fa.gz`
- `Capsella_rubella_annotation.tsv`
- `Utricularia_gibba/`

### FASTA file extensions

Major scripts accept:

- `.fa`, `.fas`, `.fasta`, `.fna`
- `.fa.gz`, `.fas.gz`, `.fasta.gz`, `.fna.gz`

### `workspace/input/species_cds`

- one CDS FASTA per species,
- accepted extensions follow the FASTA list above,
- filenames must start with the species label, such as `GENUS_SPECIES_...`,
- sequence IDs must be unique,
- every sequence ID must use the required `GENUS_SPECIES_GENEID` format:
  `GENUS_SPECIES_` followed by the gene identifier.

Examples of valid `species_cds` FASTA headers:

```fasta
>Arabidopsis_thaliana_AT1G08465.1
ATGGAGAGAAATCTTCTCTCT...
>Oryza_sativa_Os01g0100100
ATGGCGTCGACGATCGAGAT...
>Selaginella_moellendorffii_GeneID112340394
ATGGCAGCTGCTGATCGTAA...
```

The `GENUS_SPECIES` prefix in each sequence ID must match the species prefix
parsed from that FASTA filename. Avoid whitespace and special punctuation in
IDs; use letters, numbers, `_`, `.`, or `-` for portable downstream joins.

Examples that fail validation:

- `>AT1G08465.1`: missing the required species prefix,
- `>Arabidopsis_thaliana|AT1G08465.1`: contains the prohibited `|` character,
- `>Arabidopsis_thaliana_AT1G08465.1` in `Oryza_sativa_IRGSP.fa.gz`:
  sequence prefix does not match the FASTA filename prefix.

### `workspace/input/species_protein`

This directory is usually not needed for CDS-first runs. When `species_cds`
inputs are available, GeneGalleon can translate them into temporary protein
FASTA files for species-tree, orthogroup, and OMArk steps.

- optional curated input used by `gg_genome_evolution_entrypoint.sh` only when
  `input_sequence_mode="protein"`,
- one protein FASTA per species,
- accepted extensions follow the same FASTA list as `species_cds`,
- sequence IDs should stay consistent with downstream orthogroup and annotation joins,
- one gene per protein is strongly recommended when using OMArk quality assessment,
- species prefix on IDs (`GENUS_SPECIES_...`) is strongly recommended.

Provide `species_protein` only when you want GeneGalleon to use a protein set
directly instead of translating `species_cds`. Common cases:

- you have curated or official protein FASTA files that should define the
  species-tree and orthogroup inputs,
- CDS files are unavailable, incomplete, or not trustworthy enough for
  translation,
- you want OMArk or other protein-based checks to evaluate the provided
  protein products rather than CDS-derived translations,
- you need to include lineages with different genetic codes and provide
  proteins that were already translated with the correct code upstream,
- you intentionally want a protein-only genome-evolution run and accept that
  CDS-only and codon-sequence-based downstream analyses are disabled.

Unlike `species_cds`, protein-mode inputs may come from external protein sets,
so the `GENUS_SPECIES_GENEID` pattern is recommended for join compatibility but
is not enforced as a CDS validation rule.

Important behavior:

- in `input_sequence_mode="protein"`, GeneGalleon uses `species_protein`
  directly when files are present,
- in `input_sequence_mode="cds"`, GeneGalleon ignores `species_protein` and
  always generates temporary proteins from `species_cds`.
- OMArk always runs against the effective protein set; in CDS mode that means
  temporary proteins translated from `species_cds`, not the raw CDS files.
- if you want `species_genetic_code.tsv` to affect translation, leave
  `species_protein` absent so GeneGalleon translates from `species_cds`.
- using `species_protein` can include species with different genetic codes
  because GeneGalleon does not translate CDS in that path; the trade-off is
  that codon-sequence-based analyses are unavailable.

### `workspace/input/species_genetic_code/species_genetic_code.tsv`

This file is optional and is consulted only when GeneGalleon translates CDS to proteins.
The current intended use is mixed-code `gg_genome_evolution` runs with `input_sequence_mode="protein"` and no prebuilt `species_protein`.

Expected columns:

- `species`
- `genetic_code`

Rules:

- `species` must match the species filename prefix such as `Tetrahymena_thermophila`,
- `genetic_code` must be an NCBI translation-table integer such as `1` or `6`,
- species absent from the table fall back to the run's default `genetic_code`,
- rows for species not present in `species_cds` are ignored with a warning.

Minimal example:

```tsv
species	genetic_code
Tetrahymena_thermophila	6
```

### `workspace/input/query_gene`

Each file is one family-level task in `mode_gene_evolution=query2family`.
Accepted forms:

- amino-acid FASTA,
- in-frame CDS FASTA,
- plain gene ID list (one ID per line).

`gg_gene_evolution_core.sh` auto-detects query type:

- FASTA with DNA sequences: translated to AA query,
- FASTA with protein sequences: used directly,
- non-FASTA text: treated as gene IDs to extract CDS then translate.

Practical notes:

- the file basename becomes the family/task ID under `workspace/output/query2family`,
- one file must correspond to one intended gene family or family-level analysis unit,
- do not put query genes from unrelated families into the same file; otherwise,
  homolog search and tree-building steps can merge different families into one
  artificial, biologically unnatural phylogeny,
- in query2family array runs, the number of files here is the scheduler array size.

If you want to analyze multiple families, create separate input files, for
example `WOX`, `MADS`, and `YABBY`, each containing only the query genes for
that family.

Gene-list example:

```text
Arabidopsis_thaliana_AT1G08465.1
AT4G00180.1
AT1G23420.2
```

CDS FASTA example:

```fasta
>Arabidopsis_thaliana_AT1G08465.1
ATGGAGAGAAATCTTCTCTCT...
>Arabidopsis_thaliana_AT4G00180.1
ATGGCGAAGAAGATCAAGAA...
```

Protein FASTA example:

```fasta
>Arabidopsis_thaliana_AT1G08465.1
MERNLLS...
>Arabidopsis_thaliana_AT4G00180.1
MAKKIK...
```

When using gene IDs rather than FASTA, keep the identifiers consistent with
the headers in `workspace/input/species_cds`.

### `workspace/input/species_expression`

These files are optional, but they enable expression-aware downstream analyses
such as expression matrix generation, expression panels in `tree_plot`, and
expression-evolution models in `gg_gene_evolution_entrypoint.sh`.

Expected layout:

- one TSV/TSV.GZ per species,
- filename starts with the species label such as `Arabidopsis_thaliana_...`,
- first column is gene ID,
- remaining columns are tissues, stages, treatments, or other conditions,
- pre-aggregated values are recommended when you want one value per condition.

Minimal example:

```tsv
gene-id	root	leaf	flower
Arabidopsis_thaliana_AT1G08465.1	0.0000	5.3219	1.2741
Arabidopsis_thaliana_AT4G00180.1	3.1920	0.0000	0.0000
Arabidopsis_thaliana_AT1G23420.2	0.2124	2.9133	4.1170
```

Related runtime setting:

- `exp_value_type` in `workflow/gg_gene_evolution_entrypoint.sh` tells
  downstream plots and models how to interpret the input scale
  (`log2p1` by default).

### Transcriptome assembly input modes

`gg_transcriptome_generation_core.sh` defaults to `mode_transcriptome_assembly="auto"` and supports three explicit modes:

- `mode_transcriptome_assembly="sraid"`
  - input: `workspace/input/query_sra_id/GENUS_SPECIES.txt`
  - one SRA/BioProject ID per line
  - by default, GeneGalleon appends `amalgkit_sra_strategy_query` to keep transcriptome-oriented SRA strategies (`RNA-seq`, `EST`, `CLONE`); set that variable empty or broaden it if you intentionally need other strategies
- `mode_transcriptome_assembly="fastq"`
  - input: `workspace/input/species_rnaseq/GENUS_SPECIES/*.fastq.gz`
- `mode_transcriptome_assembly="metadata"`
  - input: `workspace/input/amalgkit_metadata/GENUS_SPECIES_metadata.tsv`
- `mode_transcriptome_assembly="auto"` (default)
  - auto-selects the single available input layout
  - exits with an error when multiple layouts are simultaneously present

For `mode_transcriptome_assembly="sraid"` and `mode_transcriptome_assembly="fastq"`, auto-generated amalgkit metadata is written to:

- `workspace/output/transcriptome_assembly/amalgkit_metadata/GENUS_SPECIES_metadata.tsv`

This keeps `workspace/input/amalgkit_metadata` reserved for explicit `mode_transcriptome_assembly="metadata"` inputs.

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

- output filenames are normalized to start with `GENUS_SPECIES_...`,
- formatted CDS outputs are always gzipped with `.fa.gz` extension,
- formatted GFF outputs are always gzipped with `.gff.gz` extension,
- CDS IDs are normalized into the required `GENUS_SPECIES_GENEID` pattern and
  aggregated to one representative CDS per gene,
- common historical replacements are applied to CDS/GFF text,
- CDS are padded to codon-length multiples and transcript-level redundancies are collapsed at gene level.
- when taxonomy cache preparation succeeds, the generated
  `gg_input_generation_species.tsv` also includes:
  - `taxid`
  - `nuclear_genetic_code_id` / `nuclear_genetic_code_name`
  - `mitochondrial_genetic_code_id` / `mitochondrial_genetic_code_name`
  - `plastid_genetic_code_id` / `plastid_genetic_code_name`

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
- either `cds_url`, `gbff_url`, or the pair `gff_url` + `genome_url` (or `id` to auto-resolve provider-specific URLs when supported)
  - for `provider=ncbi`, when `species_key` is omitted and `id` is given, `species_key` is inferred from NCBI species metadata (e.g. `Homo_sapiens`).
  - `provider=ncbi` accepts both `GCF_*` and `GCA_*` assembly accessions and auto-resolves NCBI assembly URLs.
  - for `provider=ddbj`, `id` can be a DDBJ BioProject accession (for example `PRJDB15739`), a WGS master accession (for example `BAAHMP000000000`), or a DDBJ BioProject URL.
    - the resolver follows the public DDBJ Search API to find the `insdc-master` accession, then downloads the anonymous-FTP WGS flatfile (`GBFF`) and derives genome/GFF/CDS from it.
  - for `provider=coge`, `id` must be CoGe `genome_id` (numeric `gid`), and CDS/GFF/Genome URLs are auto-built.
  - for `provider=cngb`, built-in inference resolves CNGB assembly IDs (`CNA...`, `cngb:...`) or linked `GCA/GCF` accessions and maps to downloadable assembly files.
  - for `provider=gwh`, `id` can be a `GWH...` accession (for example `GWHIGRM00000000.1`), a GWH assembly show URL, or a GWH folder/index URL.
    - accession-only inputs first try the public GWH download tree and then fall back to `gwhSearch/api -> /Assembly/.../show` when directory listing is unavailable.
  - for `provider=citrusgenomedb`, `id` can be a Citrus Genome Database organism page URL (for example `/organism/5799`) or an analysis page URL (for example `/Analysis/2530647`); the resolver discovers public genome/GFF/CDS downloads from the selected page.
  - for `provider=figshare`, `id` can be a public figshare article URL or numeric article id; the resolver uses the figshare article API to discover file downloads.
    - for multi-file articles, set one or more of `cds_filename`, `gff_filename`, or `genome_filename` to select the exact files within the article.
  - for `provider=plantgarden`, `id` can be a PlantGARDEN species/genome/download URL (for example `/en/list/t64480/genome/t64480.G001`); the resolver follows the PlantGARDEN assembly page to the public download index and discovers CDS or transcript FASTA together with GFF/genome files when available.
  - for `provider=plantaedb`, `id` can be a PlantaeDB species page URL (or a relative PlantaeDB `/taxa/...` path); the resolver extracts the linked NCBI `GCA/GCF` assembly accession and delegates to the NCBI resolver.
  - for `provider=flybase`, `provider=wormbase`, `provider=vectorbase`, `provider=fernbase`, `provider=veupathdb`, and `provider=dictybase`, `id` can be resolved via explicit URL columns or `GG_<PROVIDER>_*_URL_TEMPLATE`.
  - for `provider=insectbase`, `id` can be an `IBG_*` genome identifier (for example `IBG_00001`), and the resolver uses the InsectBase genome detail API to derive CDS/GFF/genome downloads.
  - for `provider=oryza_minuta`, set `id=gramene_tetraploids`; the resolver downloads the public Gramene `oryza_minutabb` and `oryza_minutacc` CDS/GFF/genome bundles and merges them into one species-level bundle automatically.
  - for `provider=direct`, set explicit `cds_url`, `gbff_url`, or `gff_url` plus `genome_url`, or provide an index-style URL in `id` that exposes downloadable files.
    When a URL points to a shared archive, use `cds_archive_member`, `gff_archive_member`, and/or `genome_archive_member` to name the member to extract.
  - for `provider=ensembl` and `provider=ensemblplants`, `id`-only URL inference is supported via:
    - provider defaults (for example, Ensembl/EnsemblPlants index discovery),
    - or env templates: `GG_<PROVIDER>_CDS_URL_TEMPLATE`, `GG_<PROVIDER>_GFF_URL_TEMPLATE`, `GG_<PROVIDER>_GENOME_URL_TEMPLATE`,
    - and optional page template: `GG_<PROVIDER>_ID_URL_TEMPLATE`.
  - for `provider=fernbase`, `id`-only URL inference can also resolve from FernBase FTP directory discovery.
  - for `provider=veupathdb`, `id`-only URL inference resolves from the VEuPathDB `GenomeDataTypes` service and derives CDS from the corresponding `AnnotatedCDSs.fasta` download.
  - for `provider=dictybase`, `id`-only resolution can use `GG_DICTYBASE_*_URL_TEMPLATE`; explicit URL columns are also supported.
  - for `provider=local`, `id` can point to a local species directory (or set explicit local file paths) and formatting starts from local files.
  - when `cds_url` is blank but both `gff_url` and `genome_url` are available, `format_species_inputs.py` derives CDS sequences from the annotation and genome FASTA during formatting.
  - when only `cds_url` is available, formatting proceeds with CDS output only and leaves GFF/genome outputs empty.
  - `provider=phycocosm` and `provider=phytozome` are not supported in `--download-manifest`.
    Use `--input-dir` for local raw-file formatting only.

Optional columns:

- `species_key`
- `cds_archive_member`, `gff_archive_member`, `gbff_archive_member`, `genome_archive_member`
  (optional paths inside `.zip` or `.tar.*` archives referenced by the corresponding `*_url`)
- `cds_filename`
- `gff_filename`
- `gbff_filename`
- `genome_filename`
- `fernbase_confidence_mode`
  (`high-confidence only` or `high-low combined`; only used for legacy FernBase releases
  that publish separate `highconfidence` and `lowconfidence` CDS/GFF files)
- `cds_url_template`, `gff_url_template`, `genome_url_template`
  (`{id}`, `{species_key}`, `{provider}` placeholders)
- `local_cds_path`, `local_gff_path`, `local_gbff_path`, `local_genome_path`
  (for `provider=local`)

XLSX template notes:

- `download_plan.xlsx` includes drop-downs for `provider` and `id`.
- `id` drop-down values are provider-specific.
- provider drop-down order is fixed, and `local` is always listed last.
- repo-curated public `provider=direct` entries can be versioned in
  `workspace/input/input_generation/direct_catalog_curated.tsv`.
- use curated direct entries for public bundles that still need explicit filenames
  or hand-vetted URLs, instead of falling back to `provider=local`.
- for large provider (`ncbi`), five model-organism IDs are shown as examples (mixed `GCF_*`/`GCA_*` formats).
- for `ddbj`, `coge`, and `cngb`, IDs are example-based by default.
- for `ddbj`, `gwh`, `citrusgenomedb`, `figshare`, `plantgarden`, `ensembl`, `ensemblplants`, `flybase`, `wormbase`, `vectorbase`, `fernbase`, `veupathdb`, `dictybase`, `insectbase`, `oryza_minuta`, `direct`, and `local`,
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
- built-in defaults stay conservative by provider; FernBase and InsectBase currently default to `2`.
- override per-provider cap with:
  `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_<PROVIDER>` (e.g. `..._NCBI=2`).
- NCBI E-utilities calls are throttled to documented limits by default:
  - `3 req/s` without API key,
  - `10 req/s` with `GG_NCBI_API_KEY`.
- override NCBI E-utilities rate with `GG_NCBI_EUTILS_MAX_RPS`.

Optional download authentication:

- `--auth-bearer-token-env ENV_NAME`
  - adds `Authorization: Bearer <token>` from an environment variable.
- `--auth-cookie-env ENV_NAME`
  - adds `Cookie: <value>` from an environment variable.
  - use this when you already have a copied session cookie for a portal.
- `--jgi-login-env ENV_NAME --jgi-password-env ENV_NAME`
  - logs into JGI Genome Portal via the public sign-on form and derives a session cookie automatically.
  - this is suitable for non-interactive Genome Portal downloads when direct anonymous access is not available.
- `--http-header "Key: Value"` (repeatable)
  - add arbitrary request headers (for provider-specific requirements).
  - this can be used for session-cookie based portals such as JGI Genome Portal, for example
    `--http-header "Cookie: <copied_session_cookie>"`.
  - local redirect tests confirmed custom headers are preserved by the current downloader across HTTP redirects.

Practical notes for JGI Genome Portal:

- if you use `--auth-cookie-env`, set the environment variable to the cookie value only.
  - do not include the literal `Cookie:` prefix.
- copied cookies should be treated as short-lived session artifacts.
  - lifetime depends on JGI server-side session policy and may be affected by idle timeout or re-login.
- copied cookies may not be reusable across different machines or browser contexts.
  - server-side session binding can depend on more than the raw cookie value.
- before launching a large download, run a small probe first.
  - use `--download-only --dry-run` when URL resolution is enough, or a single small direct URL when you need to validate authenticated download itself.

Build a manifest automatically from an existing local dataset tree:

```bash
python workflow/support/build_download_manifest.py \
  --provider all \
  --input-dir "/path/to/raw_genomes" \
  --id-options-snapshot /path/to/id_options_snapshot.json \
  --direct-catalog-manifest workspace/input/input_generation/direct_catalog_curated.tsv \
  --output workspace/input/input_generation/download_plan.xlsx
```

This writes `file://` URLs, so you can test the full download+format pipeline locally before replacing them with remote URLs.

Use `workspace/input/input_generation/download_plan.xlsx` as the runtime manifest input file.
Resolved/filled rows are written to `workspace/output/input_generation/download_plan.resolved.tsv`.

Latest template distribution:

- On each push to the default branch, GitHub Actions runs `build_download_manifest.py`
  and publishes the latest `download_plan.xlsx`.
- The workflow first builds `id_options_snapshot.json` (remote provider ID choices),
  then generates `download_plan.xlsx` from that snapshot plus
  `workspace/input/input_generation/direct_catalog_curated.tsv` when present.
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
  `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_VECTORBASE`,
  `GG_INPUT_MAX_CONCURRENT_DOWNLOADS_FERNBASE`.

Quick preset example (enable trait stage with GIFT starter config):

```bash
GG_INPUT_TRAIT_PROFILE=gift_starter bash workflow/gg_input_generation_entrypoint.sh
```
