#!/usr/bin/env bash

# Scheduler header notes:
# - Keep sections ordered as SLURM -> UGE -> PBS across entrypoints.
# - Update job name, CPU count, memory, walltime, log paths, and array size together.
# - UGE resources default to SHIROKANE AGE; SLURM/PBS site-specific lines remain examples.

# SLURM
# Common parameters: job name, cores per task, memory per core, walltime, log files, and working directory.
#SBATCH -J gg_input_generation
#SBATCH -c 4
#SBATCH --mem-per-cpu=8G
#SBATCH -t 3-00:00:00
#SBATCH --output=gg_input_generation_entrypoint.sh_%A_%a.out
#SBATCH --error=gg_input_generation_entrypoint.sh_%A_%a.err
#SBATCH --chdir=.
#SBATCH --ignore-pbs
# Array example for array-aware entrypoints.
#SBATCH -a 1
# Site-specific partition example.
#SBATCH -p epyc,rome,medium
# Optional notifications and single-node examples.
##SBATCH --mail-type=ALL
##SBATCH --mail-user=<aaa@bbb.com>

## UGE
# SHIROKANE AGE defaults: shell, working directory, slot count, memory per slot, and ljob.
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 4
#$ -l s_vmem=8G
#$ -l ljob
# Array example for array-aware entrypoints.
#$ -t 1

## PBS
# Common parameters: shell, CPU count, total memory, and exported environment.
#PBS -S /bin/bash
#PBS -l ncpus=4
#PBS -l mem=16G
# Array example for array-aware entrypoints.
#PBS -J 1
# Site-specific queue example.
#PBS -q small
#PBS -V

# Number of parallel batch jobs ("-t 1-N" in SGE or "--array 1-N" in SLURM):
# N = Number of species tasks in workspace/output/input_generation/tmp/task_plan.json when input_generation_mode=array_worker

set -euo pipefail

echo "$(date): Starting"

# Resolve workflow paths for local and scheduler-spooled execution.
gg_bootstrap_script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
gg_bootstrap_checked_bases=""
if [[ -n "${KFAUTO_RUNTIME_RELEASE_ROOT:-}" ]]; then
  gg_bootstrap_bases=("${KFAUTO_RUNTIME_RELEASE_ROOT}")
else
  gg_bootstrap_bases=(
    "${SLURM_SUBMIT_DIR:-}"
    "${PBS_O_WORKDIR:-}"
    "${PWD:-}"
    "${gg_bootstrap_script_dir}"
  )
fi
for gg_bootstrap_base in "${gg_bootstrap_bases[@]}"
do
  [[ -n "${gg_bootstrap_base}" ]] || continue
  case ":${gg_bootstrap_checked_bases}:" in
    *":${gg_bootstrap_base}:"*) continue ;;
  esac
  gg_bootstrap_checked_bases="${gg_bootstrap_checked_bases:+${gg_bootstrap_checked_bases}:}${gg_bootstrap_base}"
  for bootstrap_path in \
    "${gg_bootstrap_base}/support/gg_entrypoint_bootstrap.sh" \
    "${gg_bootstrap_base}/workflow/support/gg_entrypoint_bootstrap.sh"
  do
    if [[ -s "${bootstrap_path}" ]]; then
      # shellcheck disable=SC1090
      source "${bootstrap_path}"
      break
    fi
  done
  if declare -F gg_entrypoint_initialize >/dev/null 2>&1; then
    break
  fi
done
unset gg_bootstrap_base gg_bootstrap_bases gg_bootstrap_checked_bases gg_bootstrap_script_dir bootstrap_path
if ! declare -F gg_entrypoint_initialize >/dev/null 2>&1; then
  echo "Failed to locate gg_entrypoint_bootstrap.sh from BASH_SOURCE[0]=${BASH_SOURCE[0]}" >&2
  exit 1
fi
if ! gg_entrypoint_initialize "${BASH_SOURCE[0]}" 1 "gg_input_generation"; then
  exit 1
fi
gg_entrypoint_name="gg_input_generation_entrypoint.sh"

### Start: Modify this block to tailor your analysis ###

# Workflow flags
run_gene_model_refinement=0 # Opt-in all-isoform selection, existing-model revision and isoform addition.
gene_model_refinement_dir="" # Blank uses output/input_generation/gene_model_refinement.
gene_model_refinement_policy="conserved" # longest|conserved; representative ranking policy.
gene_model_refinement_mode="conservative" # off|audit|conservative; audit publishes proposals without accepting predictions.
gene_model_refinement_inputs="" # Optional raw CDS/GFF/genome TSV; requires an explicit frozen correspondence table.
gene_model_refinement_edges="" # Optional trusted synteny locus correspondence TSV.
gene_model_refinement_rescue_dir="" # Existing frozen rescue plan; blank uses gene_model_rescue_dir.
gene_model_refinement_rna="" # Optional complete coding RNA paths TSV, with zero-based half-open CDS blocks.
gene_model_refinement_min_margin=0.10 # Experimental minimum candidate score margin; not a calibrated probability.
gene_model_refinement_min_support=2 # Minimum independent donor species for accepting homology predictions.
gene_model_refinement_candidate_limit=32 # Maximum candidate paths per locus for bounded selection/prediction.
gene_model_refinement_padding=2000 # Local prediction padding in genomic bases.
run_gene_model_rescue=0 # Opt-in sparse synteny/model rescue after initial formatting, BUSCO and taxonomy.
gene_model_rescue_tree="auto" # External initial species tree, or auto for output/species_taxonomy/taxonomy_tree.nwk.
gene_model_rescue_dir="" # Rescue outputs; blank uses output/input_generation/gene_model_rescue.
gene_model_rescue_comparison_cache="" # Shared hashed comparison cache; blank uses the rescue output's parent.
gene_model_rescue_interval_workers=0 # Independent interval workers; 0 uses task CPUs, total miniprot threads stay within task CPUs.
gene_model_rescue_common_references=5 # Common high-BUSCO species selected by NWKIT max-pd.
gene_model_rescue_nearest_references=3 # Additional nearest species per target; redundant pairs are computed once.
gene_model_rescue_minimum_busco=90 # Minimum complete BUSCO percentage for common references; S+D, same dataset/version/mode.
gene_model_rescue_minimum_coverage=0.95 # Minimum donor protein coverage for automatic model acceptance.
gene_model_rescue_minimum_identity=0.5 # Minimum protein identity for automatic model acceptance.
gene_model_rescue_max_interval=200000 # Maximum interval between target flanking anchors, in bp.
gene_model_rescue_max_intron=20000 # Maximum intron size for miniprot, in bp.
gene_model_rescue_genome_fallback=1 # Search unresolved syntenic candidates across the target genome; outside-block hits remain unresolved.
gene_model_rescue_gemoma_jar="" # Optional GeMoMa refinement jar; blank uses miniprot alone.
gene_model_rescue_gemoma_java="java" # Java executable compatible with the supplied GeMoMa jar; GeMoMa 1.9 requires a JavaScript engine.
run_species_taxonomy=1 # Resolve NCBI taxonomic ranks and plot them on the available species tree or an NCBI taxonomy tree.
run_format_inputs=1 # Format local inputs or download-manifest targets into workspace layout.
run_validate_inputs=1 # Validate formatted inputs before downstream workflows use them.
run_cds_fx2tab=1 # Run seqkit fx2tab for formatted species CDS files.
run_species_busco=1 # Run BUSCO for formatted species CDS files.
busco_timeout_seconds=0 # 0 disables the per-species BUSCO wall-time limit; set a positive number for large arrays.
run_multispecies_summary=1 # Generate multi-species BUSCO summary plots and tables from species BUSCO outputs.
run_generate_species_trait=0 # Generate species_trait.tsv from downloaded or local metadata sources.

# Taxonomic annotation parameters
taxonomy_species_tree="auto" # Species-tree Newick path, or auto to discover a selected tree before using NCBI taxonomy.
taxonomy_ranks="all" # all includes every available NCBI lineage rank; alternatively use a comma-separated list. Per-species missing ranks remain blank.
taxonomy_plot_clades=0 # Set to 1 to draw clade columns; tables and NHX retain clades either way.
taxonomy_taxid_map="" # Optional TSV with species and taxid columns for explicit taxonomy corrections.
taxonomy_taxid_override="" # One explicit species:TaxID correction for a scheduled run.

# Shared parameters
provider="all" # all|ensembl|ensemblplants|ensemblmetazoa|ensemblprotists|phycocosm|phytozome|ncbi|ddbj|refseq|genbank|coge|cngb|flybase|wormbase|vectorbase|fernbase|insectbase|local; selects which provider-specific local layout or download-manifest rows are formatted, with all scanning every supported provider directory.
input_generation_mode="single" # single | array_prepare | array_worker | array_finalize | rescue_prepare | rescue_synteny | rescue_models | rescue_finalize | refinement_prepare | refinement_catalog | refinement_correspondence | refinement_predict | refinement_finalize. Workers use GG_ARRAY_TASK_ID.
species_busco_parallel_jobs="auto" # In single mode, auto runs up to four species BUSCO jobs within GG_TASK_CPUS; array_worker still runs one species per task.
species_busco_memory_gb_per_job=4 # Minimum tool-memory budget per concurrent BUSCO species job; parallelism is capped by GG_MEM_TOOL_GB / this value.
trait_profile="none" # none|gift_starter|gbif_distribution; optional preset for generating species_trait.tsv from external trait databases.
busco_lineage="${GG_COMMON_BUSCO_LINEAGE:-auto}" # BUSCO lineage dataset name, or auto to infer a shared dataset from the discovered species set.
strict=0 # Treat input formatting and validation warnings as fatal errors.
require_cds=0 # 0|1; require a nonempty formatted CDS FASTA for every selected species (including CDS derived from GFF+genome or GBFF).
require_gff=0 # 0|1; require a formatted GFF with at least one feature for every selected species (including GFF derived from GBFF).
require_genome=0 # 0|1; require a nonempty formatted genome FASTA for every selected species. Default off for annotation-only projects.
overwrite=0 # Regenerate formatted/downloaded outputs even when existing non-empty outputs are present.
download_only=0 # Single mode only: stop after manifest downloads, before formatting or downstream processing.
dry_run=0 # Print planned downloads/formatting actions without writing formatted outputs.
download_timeout=120 # Per-request timeout in seconds for remote downloads.
gene_grouping_mode="rescue_overlap" # strict|rescue_overlap; strict keeps provider gene models as-is, while rescue_overlap merges likely fragmented/overlapping CDS records into a gene-level representative when possible.
gff_repair_mode="safe" # off|safe|strict; safe repairs only unique collision-free GFF gene IDs against final CDS IDs, while strict rejects ambiguous repair candidates.
trait_species_source="download_manifest" # download_manifest|species_cds; source used to decide which species are included during trait-table generation.
trait_databases="auto" # auto|all|comma-separated IDs; trait databases queried by the selected trait_profile, with auto choosing profile defaults.

# Request parameters
auth_bearer_token_env="" # Environment variable name containing a bearer token for authenticated download-manifest URLs, e.g., GFE_DOWNLOAD_BEARER_TOKEN.
http_header="" # Extra HTTP header forwarded to download requests, e.g., "User-Agent: genegalleon-input-generation".

# Path and output parameters
input_dir="" # Local raw input directory to ingest instead of downloading.
download_manifest="" # Path to the download manifest file.
download_dir="" # Directory for downloaded raw files.
download_limit_dir="" # Shared atomic-namespace directory for database request limits; blank uses workspace/.gg_cache/input_download_limits. Stop old clients before migrating protocols.
summary_output="" # Output path for the run summary table.
species_cds_dir="" # Output directory for formatted CDS FASTA files.
species_cds_fx2tab_dir="" # Output directory for CDS fx2tab TSV files.
species_busco_full_dir="" # Output directory for BUSCO full tables under output/input_generation/.
species_busco_short_dir="" # Output directory for BUSCO short summaries under output/input_generation/.
species_gff_dir="" # Output directory for formatted GFF files.
species_genome_dir="" # Output directory for formatted genome FASTA files.
species_summary_output="" # Output path for the species-level summary table.
resolved_manifest_output="" # Output path for the resolved download-manifest TSV.
species_trait_output="" # Output path for the generated species trait table.
task_plan_output="" # Output path for the discovered array task-plan JSON.
resume_from_task_plan="" # Optional frozen donor plan for importing verified completed formatting/validation into a new workspace.
resume_from_task_plan_sha256="" # Required SHA-256 of resume_from_task_plan.
resume_from_input_generation_root="" # Required donor output/input_generation directory; donor must have no active workers.
resume_fallback_task_plan="" # Optional second donor, tried only when the primary has no reusable species format.
resume_fallback_task_plan_sha256="" # Required SHA-256 of the fallback donor plan.
resume_fallback_input_generation_root="" # Required fallback output/input_generation directory; prepare verifies both inactive donors.
trait_plan="" # Optional trait plan file describing requested traits.
trait_database_sources="" # Optional mapping file that defines trait database sources.
trait_download_dir="" # Directory for cached or raw trait database downloads.
trait_download_timeout=120 # Per-request timeout in seconds for trait database downloads.
gbif_api="" # Optional GBIF API base URI override for trait_profile=gbif_distribution.
gbif_page_size="" # Number of occurrence records requested per GBIF API page.
gbif_max_occurrences_per_species="" # Maximum no-login GBIF occurrence records fetched per species before summarizing distribution traits.
gbif_grid_degrees="" # Grid size in degrees for observed occupied-cell area (not IUCN AOO).
gbif_min_match_confidence="" # Minimum GBIF species-match confidence required before using occurrence records.
gbif_max_coordinate_uncertainty_m="" # Optional maximum GBIF coordinate uncertainty in meters; blank keeps GBIF records regardless of uncertainty.
gbif_min_distance_from_known_centroid_m="" # Minimum distance from a known georeferencing centroid in meters; blank disables proximity filtering.
gbif_year_min="" # Minimum event year; date intervals must lie inside the selected window.
gbif_year_max="" # Maximum event year.
gbif_countries="" # Comma-separated country codes to retain.
gbif_include_basis_of_record="" # Comma-separated basisOfRecord values to retain.
gbif_exclude_basis_of_record="" # Comma-separated basisOfRecord values to exclude.
gbif_include_establishment_means="" # Retain specified establishmentMeans values; unknown is not native.
gbif_missing_date="" # keep|exclude for unknown dates with an active year filter; default exclude.
gbif_missing_uncertainty="" # keep|exclude for unknown uncertainty with an active threshold; default keep.
gbif_missing_centroid_distance="" # keep|exclude for missing known-centroid distance; default keep.
gbif_use_cache="" # yes|no; reuse verified records or acquire a new snapshot; default yes.
gbif_require_complete="" # yes|no; fail on incomplete acquisition instead of publishing NA; default no.
gbif_occurrence_file="" # Local GBIF SIMPLE_CSV table (.tsv, .csv, .gz or single-table .zip).
gbif_taxon_map="" # Reviewed TSV with species, taxon_key and scientific_name for local records.
gbif_download_metadata="" # Saved official GBIF download metadata JSON to verify completeness.

### End: Modify this block to tailor your analysis ###

source "${gg_support_dir}/gg_util.sh" # loading utility functions

if [[ -n "${GG_INPUT_GBIF_MAX_DISTANCE_FROM_CENTROID_M:-}" || -n "${gbif_max_distance_from_centroid_m:-}" ]]; then
  echo "Removed GBIF maximum-centroid option: use gbif_min_distance_from_known_centroid_m / GG_INPUT_GBIF_MIN_DISTANCE_FROM_KNOWN_CENTROID_M."
  exit 1
fi
# Apply documented GG_INPUT_* overrides, then forward canonical config variables.
gg_apply_registered_env_overrides "${gg_entrypoint_name}"
forward_config_vars_to_container_env "${gg_entrypoint_name}"

# Keep per-provider input downloads modest by default; callers can override these
# environment variables for sites with different network limits.
: "${GG_INPUT_MAX_CONCURRENT_DOWNLOADS_COGE:=2}"
: "${GG_INPUT_MAX_CONCURRENT_DOWNLOADS_GWH:=2}"
: "${GG_INPUT_MAX_CONCURRENT_DOWNLOADS_CNGB:=1}"
: "${GG_INPUT_MAX_CONCURRENT_DOWNLOADS_DIRECT:=2}"
export GG_INPUT_MAX_CONCURRENT_DOWNLOADS_COGE
export GG_INPUT_MAX_CONCURRENT_DOWNLOADS_GWH
export GG_INPUT_MAX_CONCURRENT_DOWNLOADS_CNGB
export GG_INPUT_MAX_CONCURRENT_DOWNLOADS_DIRECT

# Provider-specific download caps are consumed directly downstream.
gg_forward_env_vars_with_prefix_to_container_env "GG_INPUT_MAX_CONCURRENT_DOWNLOADS_"
gg_forward_env_vars_with_prefix_to_container_env "GG_INPUT_REQUEST_INTERVAL_"
gg_forward_env_vars_with_prefix_to_container_env "GG_INPUT_DOWNLOAD_LIMIT_"
gg_forward_env_vars_with_prefix_to_container_env "GG_DOWNLOAD_"

if ! gg_entrypoint_prepare_container_runtime 0; then
  exit 1
fi
gg_entrypoint_activate_container_runtime

gg_entrypoint_enter_workspace
gg_runtime_core_script="$(gg_prepare_entrypoint_runtime_snapshot "${gg_entrypoint_name}" "${gg_core_dir}/gg_input_generation_core.sh")"
gg_run_container_shell_script "${gg_container_image_path}" "${gg_runtime_core_script}"
gg_require_versions_dump "${gg_entrypoint_name}"

echo "$(date): Ending"
