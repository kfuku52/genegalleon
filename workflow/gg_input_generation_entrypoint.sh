#!/usr/bin/env bash

# SLURM in NIG supercomputer
#SBATCH -J gg_input_generation
#SBATCH -c 2
#SBATCH --mem-per-cpu=8G
#SBATCH -t 2976:00:00
#SBATCH --output=gg_input_generation_%j.out
#SBATCH --error=gg_input_generation_%j.err
#SBATCH -p medium
#SBATCH --chdir=.
##SBATCH --mail-type=ALL
##SBATCH --mail-user=<aaa@bbb.com>
##SBATCH -N 1-1 # Number of nodes
##SBATCH -n 1 # Number of tasks

#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 2
#$ -l s_vmem=8G
#$ -l mem_req=8G
##$ -l epyc
#$ -l d_rt=124:00:00:00
#$ -l s_rt=124:00:00:00

# Number of parallel batch jobs ("-t 1-N" in SGE or "--array 1-N" in SLURM):
# Fixed to 1

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
if ! gg_entrypoint_initialize "${BASH_SOURCE[0]}" 0; then
  exit 1
fi

### Start: Modify this block to tailor your analysis ###

run_format_inputs=1
run_validate_inputs=1
run_generate_species_trait=0

provider="all" # all|ensembl|ensemblplants|phycocosm|phytozome|ncbi|refseq|genbank|coge|cngb|flybase|wormbase|vectorbase|local
trait_profile="none" # none|gift_starter
strict=0
overwrite=0
download_only=0
dry_run=0
download_timeout=120
trait_species_source="download_manifest" # download_manifest|species_cds
trait_databases="auto" # auto|all|comma-separated IDs

# Optional download request settings.
auth_bearer_token_env="" # e.g., GFE_DOWNLOAD_BEARER_TOKEN
http_header="" # e.g., "User-Agent: genegalleon-input-generation"

# Optional local raw input directory and manifest.
input_dir=""
download_manifest=""
download_dir=""
summary_output=""
species_cds_dir=""
species_gff_dir=""
species_genome_dir=""
species_summary_output=""
resolved_manifest_output=""
species_trait_output=""
trait_plan=""
trait_database_sources=""
trait_download_dir=""
trait_download_timeout=120

### End: Modify this block to tailor your analysis ###

source "${gg_support_dir}/gg_util.sh" # loading utility functions

# Apply documented GG_INPUT_* overrides before forwarding canonical config vars.
gg_apply_named_env_overrides \
  provider GG_INPUT_PROVIDER \
  strict GG_INPUT_STRICT \
  overwrite GG_INPUT_OVERWRITE \
  download_only GG_INPUT_DOWNLOAD_ONLY \
  dry_run GG_INPUT_DRY_RUN \
  download_timeout GG_INPUT_DOWNLOAD_TIMEOUT \
  input_dir GG_INPUT_INPUT_DIR \
  download_manifest GG_INPUT_DOWNLOAD_MANIFEST \
  download_dir GG_INPUT_DOWNLOAD_DIR \
  summary_output GG_INPUT_SUMMARY_OUTPUT \
  auth_bearer_token_env GG_INPUT_AUTH_BEARER_TOKEN_ENV \
  http_header GG_INPUT_HTTP_HEADER \
  run_format_inputs GG_INPUT_RUN_FORMAT_INPUTS \
  run_validate_inputs GG_INPUT_RUN_VALIDATE_INPUTS \
  run_generate_species_trait GG_INPUT_RUN_GENERATE_SPECIES_TRAIT \
  trait_profile GG_INPUT_TRAIT_PROFILE \
  species_cds_dir GG_INPUT_SPECIES_CDS_DIR \
  species_gff_dir GG_INPUT_SPECIES_GFF_DIR \
  species_genome_dir GG_INPUT_SPECIES_GENOME_DIR \
  species_summary_output GG_INPUT_SPECIES_SUMMARY_OUTPUT \
  resolved_manifest_output GG_INPUT_RESOLVED_MANIFEST_OUTPUT \
  species_trait_output GG_INPUT_SPECIES_TRAIT_OUTPUT \
  trait_plan GG_INPUT_TRAIT_PLAN \
  trait_database_sources GG_INPUT_TRAIT_DATABASE_SOURCES \
  trait_download_dir GG_INPUT_TRAIT_DOWNLOAD_DIR \
  trait_download_timeout GG_INPUT_TRAIT_DOWNLOAD_TIMEOUT \
  trait_species_source GG_INPUT_TRAIT_SPECIES_SOURCE \
  trait_databases GG_INPUT_TRAIT_DATABASES

# Forward canonical config variables into the container environment.
forward_config_vars_to_container_env "${BASH_SOURCE[0]}"

# Provider-specific download caps are consumed directly downstream.
gg_forward_env_vars_with_prefix_to_container_env "GG_INPUT_MAX_CONCURRENT_DOWNLOADS_"

if ! gg_entrypoint_prepare_container_runtime 0; then
  exit 1
fi
gg_entrypoint_activate_container_runtime

gg_entrypoint_enter_workspace
${singularity_command} "${gg_container_image_path}" < "${gg_core_dir}/gg_input_generation_core.sh"
if ! gg_trigger_versions_dump "$(basename "${BASH_SOURCE[0]}")"; then
  echo "Warning: gg_versions trigger failed."
fi

echo "$(date): Ending"
