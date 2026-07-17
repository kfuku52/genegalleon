#!/usr/bin/env bash

# Scheduler header notes:
# - Keep sections ordered as SLURM -> UGE -> PBS across entrypoints.
# - Update job name, CPU count, memory, walltime, log paths, and array size together.
# - Site-specific partition/queue/resource lines stay commented examples by default.

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
# Common parameters: shell, working directory, slot count, memory per slot, and runtime limits.
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 4
#$ -l s_vmem=8G
#$ -l mem_req=8G
# Site-specific resource example.
#$ -l epyc
#$ -l d_rt=124:00:00:00
#$ -l s_rt=124:00:00:00
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
gg_bootstrap_submit_dir="${SLURM_SUBMIT_DIR:-${PBS_O_WORKDIR:-${PWD:-}}}"
gg_bootstrap_script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
for bootstrap_path in \
  "${gg_bootstrap_submit_dir}/support/gg_entrypoint_bootstrap.sh" \
  "${gg_bootstrap_submit_dir}/workflow/support/gg_entrypoint_bootstrap.sh" \
  "${gg_bootstrap_script_dir}/support/gg_entrypoint_bootstrap.sh"
do
  if [[ -s "${bootstrap_path}" ]]; then
    # shellcheck disable=SC1090
    source "${bootstrap_path}"
    break
  fi
done
unset gg_bootstrap_submit_dir gg_bootstrap_script_dir bootstrap_path
if ! declare -F gg_entrypoint_initialize >/dev/null 2>&1; then
  echo "Failed to locate gg_entrypoint_bootstrap.sh from BASH_SOURCE[0]=${BASH_SOURCE[0]}" >&2
  exit 1
fi
if ! gg_entrypoint_initialize "${BASH_SOURCE[0]}" 1; then
  exit 1
fi
gg_entrypoint_name="gg_input_generation_entrypoint.sh"

### Start: Modify this block to tailor your analysis ###

# Workflow flags
run_format_inputs=1 # Format local inputs or download-manifest targets into workspace layout.
run_validate_inputs=1 # Validate formatted inputs before downstream workflows use them.
run_cds_fx2tab=1 # Run seqkit fx2tab for formatted species CDS files.
run_species_busco=1 # Run BUSCO for formatted species CDS files.
run_multispecies_summary=1 # Generate multi-species BUSCO summary plots and tables from species BUSCO outputs.
run_generate_species_trait=0 # Generate species_trait.tsv from downloaded or local metadata sources.

# Shared parameters
provider="all" # all|ensembl|ensemblplants|ensemblmetazoa|ensemblprotists|phycocosm|phytozome|ncbi|ddbj|refseq|genbank|coge|cngb|flybase|wormbase|vectorbase|fernbase|insectbase|local; selects which provider-specific local layout or download-manifest rows are formatted, with all scanning every supported provider directory.
input_generation_mode="single" # single=all stages in one run | array_prepare=build task plan | array_worker=run one species task per GG_ARRAY_TASK_ID | array_finalize=merge shards and run shared validation/summaries; use array_* for large downloads/formatting across many species.
trait_profile="none" # none|gift_starter|gbif_distribution; optional preset for generating species_trait.tsv from external trait databases.
busco_lineage="${GG_COMMON_BUSCO_LINEAGE:-auto}" # BUSCO lineage dataset name, or auto to infer a shared dataset from the discovered species set.
strict=0 # Treat input formatting and validation warnings as fatal errors.
overwrite=0 # Regenerate formatted/downloaded outputs even when existing non-empty outputs are present.
download_only=0 # Stop after provider/download-manifest download and formatting; skip validation, BUSCO, summaries, and trait generation.
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
trait_plan="" # Optional trait plan file describing requested traits.
trait_database_sources="" # Optional mapping file that defines trait database sources.
trait_download_dir="" # Directory for cached or raw trait database downloads.
trait_download_timeout=120 # Per-request timeout in seconds for trait database downloads.
gbif_api="" # Optional GBIF API base URI override for trait_profile=gbif_distribution.
gbif_page_size="" # Number of occurrence records requested per GBIF API page.
gbif_max_occurrences_per_species="" # Maximum no-login GBIF occurrence records fetched per species before summarizing distribution traits.
gbif_grid_degrees="" # Latitude/longitude grid size in degrees used to estimate occupied GBIF area.
gbif_min_match_confidence="" # Minimum GBIF species-match confidence required before using occurrence records.
gbif_max_coordinate_uncertainty_m="" # Optional maximum GBIF coordinate uncertainty in meters; blank keeps GBIF records regardless of uncertainty.
gbif_max_distance_from_centroid_m="" # Optional maximum distance from the species occurrence centroid in meters; blank disables centroid-distance filtering.

### End: Modify this block to tailor your analysis ###

source "${gg_support_dir}/gg_util.sh" # loading utility functions

# Apply documented GG_INPUT_* overrides, then forward canonical config variables.
gg_apply_registered_env_overrides "${gg_entrypoint_name}"
forward_config_vars_to_container_env "${gg_entrypoint_name}"

# Keep per-provider input downloads modest by default; callers can override these
# environment variables for sites with different network limits.
: "${GG_INPUT_MAX_CONCURRENT_DOWNLOADS_COGE:=2}"
: "${GG_INPUT_MAX_CONCURRENT_DOWNLOADS_GWH:=2}"
: "${GG_INPUT_MAX_CONCURRENT_DOWNLOADS_CNGB:=2}"
: "${GG_INPUT_MAX_CONCURRENT_DOWNLOADS_DIRECT:=2}"
export GG_INPUT_MAX_CONCURRENT_DOWNLOADS_COGE
export GG_INPUT_MAX_CONCURRENT_DOWNLOADS_GWH
export GG_INPUT_MAX_CONCURRENT_DOWNLOADS_CNGB
export GG_INPUT_MAX_CONCURRENT_DOWNLOADS_DIRECT

# Provider-specific download caps are consumed directly downstream.
gg_forward_env_vars_with_prefix_to_container_env "GG_INPUT_MAX_CONCURRENT_DOWNLOADS_"

if ! gg_entrypoint_prepare_container_runtime 0; then
  exit 1
fi
gg_entrypoint_activate_container_runtime

gg_entrypoint_enter_workspace
gg_runtime_core_script="$(gg_prepare_entrypoint_runtime_snapshot "${gg_entrypoint_name}" "${gg_core_dir}/gg_input_generation_core.sh")"
gg_run_container_shell_script "${gg_container_image_path}" "${gg_runtime_core_script}"
gg_require_versions_dump "${gg_entrypoint_name}"

echo "$(date): Ending"
