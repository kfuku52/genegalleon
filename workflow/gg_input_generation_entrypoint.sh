#!/usr/bin/env bash
set -euo pipefail

# SLURM in NIG supercomputer
#SBATCH -J gg_input_generation
#SBATCH -c 2 # Number of CPUs
#SBATCH --mem-per-cpu=8G # RAM per CPU in GB
#SBATCH -t 2976:00:00 # maximum time in d-hh:mm:ss format. NIG supercomputer epyc/medium MaxTime=2976:00:00
#SBATCH --output=gg_input_generation_%j.out
#SBATCH --error=gg_input_generation_%j.err
#SBATCH -p medium # partition name, cluster environment specific
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

echo "$(date): Starting"

# Change these directories for your custom-made analysis
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
dir_pg="${script_dir}/../workspace" # pg input and output directory
dir_script="${script_dir}" # directory where gg_util.sh and core/gg_*_core.sh locate
gg_image="${script_dir}/../genegalleon.sif" # path to the singularity image

### Start: Modify this block to tailor your analysis ###

run_format_inputs=1
run_validate_inputs=1

provider="all" # all|ensembl|ensemblplants|phycocosm|phytozome|ncbi|refseq|genbank|coge|cngb|flybase|wormbase|vectorbase|local
strict=0
overwrite=0
download_only=0
dry_run=0
download_timeout=120

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

### End: Modify this block to tailor your analysis ###

source "${dir_script}/support/gg_util.sh" # loading utility functions
# Forward config variables (including external overrides) into container environment.
forward_config_vars_to_container_env "${BASH_SOURCE[0]}"

# Forward GG_INPUT_* runtime overrides consumed by gg_input_generation_core.sh.
for gg_input_var_name in ${!GG_INPUT_@}; do
  gg_export_var_to_container_env_if_set "${gg_input_var_name}"
done

gg_scheduler_runtime_prelude
unset_singularity_envs
if ! set_singularity_command; then
  exit 1
fi
variable_SGEnizer
set_singularityenv
gg_print_scheduler_runtime_summary

mkdir -p "${dir_pg}"
cd "${dir_pg}"
${singularity_command} "${gg_image}" < "${dir_script}/core/gg_input_generation_core.sh"
if ! gg_trigger_versions_dump "$(basename "${BASH_SOURCE[0]}")"; then
  echo "Warning: gg_versions trigger failed."
fi

echo "$(date): Ending"
