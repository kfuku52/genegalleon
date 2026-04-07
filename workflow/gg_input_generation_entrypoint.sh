#!/usr/bin/env bash

# Scheduler header notes:
# - Keep sections ordered as SLURM -> UGE -> PBS across entrypoints.
# - Update job name, CPU count, memory, walltime, log paths, and array size together.
# - Site-specific partition/queue/resource lines stay commented examples by default.

# SLURM
# Common parameters: job name, cores per task, memory per core, walltime, log files, and working directory.
#SBATCH -J gg_input_generation
#SBATCH -c 2
#SBATCH --mem-per-cpu=8G
#SBATCH -t 2976:00:00
#SBATCH --output=gg_input_generation_entrypoint.sh_%A_%a.out
#SBATCH --error=gg_input_generation_entrypoint.sh_%A_%a.err
#SBATCH --chdir=.
#SBATCH --ignore-pbs
# Array example for array-aware entrypoints.
#SBATCH -a 1
# Site-specific partition example.
#SBATCH -p epyc
# Optional notifications and single-node examples.
##SBATCH --mail-type=ALL
##SBATCH --mail-user=<aaa@bbb.com>

## UGE
# Common parameters: shell, working directory, slot count, memory per slot, and runtime limits.
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 2
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
#PBS -l ncpus=2
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
provider="all" # all|ensembl|ensemblplants|phycocosm|phytozome|ncbi|ddbj|refseq|genbank|coge|cngb|flybase|wormbase|vectorbase|fernbase|insectbase|local
input_generation_mode="single" # single=all stages in one run | array_prepare=build task plan | array_worker=run one species task per GG_ARRAY_TASK_ID | array_finalize=merge shards and run shared validation/summaries
trait_profile="none" # none|gift_starter
busco_lineage="${GG_COMMON_BUSCO_LINEAGE:-auto}" # BUSCO lineage dataset name, or auto to infer a shared dataset from the discovered species set.
strict=0 # Treat input formatting and validation warnings as errors.
overwrite=0 # Overwrite existing formatted or downloaded outputs.
download_only=0 # Stop after download/format steps without validation-only downstream work.
dry_run=0 # Print planned actions without downloading or writing outputs.
download_timeout=120 # Per-request timeout in seconds for remote downloads.
gene_grouping_mode="rescue_overlap" # strict|rescue_overlap
trait_species_source="download_manifest" # download_manifest|species_cds
trait_databases="auto" # auto|all|comma-separated IDs

# Request parameters
auth_bearer_token_env="" # e.g., GFE_DOWNLOAD_BEARER_TOKEN
http_header="" # e.g., "User-Agent: genegalleon-input-generation"

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

### End: Modify this block to tailor your analysis ###

source "${gg_support_dir}/gg_util.sh" # loading utility functions

# Apply documented GG_INPUT_* overrides before forwarding canonical config vars.
gg_apply_named_env_overrides \
  provider GG_INPUT_PROVIDER \
  input_generation_mode GG_INPUT_INPUT_GENERATION_MODE \
  strict GG_INPUT_STRICT \
  overwrite GG_INPUT_OVERWRITE \
  download_only GG_INPUT_DOWNLOAD_ONLY \
  dry_run GG_INPUT_DRY_RUN \
  download_timeout GG_INPUT_DOWNLOAD_TIMEOUT \
  gene_grouping_mode GG_INPUT_GENE_GROUPING_MODE \
  busco_lineage GG_INPUT_BUSCO_LINEAGE \
  input_dir GG_INPUT_INPUT_DIR \
  download_manifest GG_INPUT_DOWNLOAD_MANIFEST \
  download_dir GG_INPUT_DOWNLOAD_DIR \
  summary_output GG_INPUT_SUMMARY_OUTPUT \
  auth_bearer_token_env GG_INPUT_AUTH_BEARER_TOKEN_ENV \
  http_header GG_INPUT_HTTP_HEADER \
  run_format_inputs GG_INPUT_RUN_FORMAT_INPUTS \
  run_validate_inputs GG_INPUT_RUN_VALIDATE_INPUTS \
  run_cds_fx2tab GG_INPUT_RUN_CDS_FX2TAB \
  run_species_busco GG_INPUT_RUN_SPECIES_BUSCO \
  run_multispecies_summary GG_INPUT_RUN_MULTISPECIES_SUMMARY \
  run_generate_species_trait GG_INPUT_RUN_GENERATE_SPECIES_TRAIT \
  trait_profile GG_INPUT_TRAIT_PROFILE \
  species_cds_dir GG_INPUT_SPECIES_CDS_DIR \
  species_cds_fx2tab_dir GG_INPUT_SPECIES_CDS_FX2TAB_DIR \
  species_busco_full_dir GG_INPUT_SPECIES_BUSCO_FULL_DIR \
  species_busco_short_dir GG_INPUT_SPECIES_BUSCO_SHORT_DIR \
  species_gff_dir GG_INPUT_SPECIES_GFF_DIR \
  species_genome_dir GG_INPUT_SPECIES_GENOME_DIR \
  species_summary_output GG_INPUT_SPECIES_SUMMARY_OUTPUT \
  resolved_manifest_output GG_INPUT_RESOLVED_MANIFEST_OUTPUT \
  species_trait_output GG_INPUT_SPECIES_TRAIT_OUTPUT \
  task_plan_output GG_INPUT_TASK_PLAN_OUTPUT \
  trait_plan GG_INPUT_TRAIT_PLAN \
  trait_database_sources GG_INPUT_TRAIT_DATABASE_SOURCES \
  trait_download_dir GG_INPUT_TRAIT_DOWNLOAD_DIR \
  trait_download_timeout GG_INPUT_TRAIT_DOWNLOAD_TIMEOUT \
  trait_species_source GG_INPUT_TRAIT_SPECIES_SOURCE \
  trait_databases GG_INPUT_TRAIT_DATABASES

# Forward canonical config variables into the container environment.
forward_config_vars_to_container_env "${gg_entrypoint_name}"

# Provider-specific download caps are consumed directly downstream.
gg_forward_env_vars_with_prefix_to_container_env "GG_INPUT_MAX_CONCURRENT_DOWNLOADS_"

if ! gg_entrypoint_prepare_container_runtime 0; then
  exit 1
fi
gg_entrypoint_activate_container_runtime

gg_entrypoint_enter_workspace
gg_run_container_shell_script "${gg_container_image_path}" "${gg_core_dir}/gg_input_generation_core.sh"
gg_require_versions_dump "${gg_entrypoint_name}"

echo "$(date): Ending"
