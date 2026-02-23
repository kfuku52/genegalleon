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

run_build_manifest="${run_build_manifest:-1}"
run_format_inputs="${run_format_inputs:-1}"
run_validate_inputs="${run_validate_inputs:-1}"

provider="${provider:-all}" # all|ensemblplants|phycocosm|phytozome
strict="${strict:-0}"
overwrite="${overwrite:-0}"
download_only="${download_only:-0}"
dry_run="${dry_run:-0}"
download_timeout="${download_timeout:-120}"

# Optional download request settings.
auth_bearer_token_env="${auth_bearer_token_env:-}" # e.g., GFE_DOWNLOAD_BEARER_TOKEN
http_header="${http_header:-}" # e.g., "User-Agent: genegalleon-inputprep"

# Optional runtime config file.
config_file="${config_file:-}"

# Optional local raw dataset roots and manifests.
dataset_root="${dataset_root:-}"
input_dir="${input_dir:-}"
download_manifest="${download_manifest:-}"
download_dir="${download_dir:-}"
manifest_output="${manifest_output:-}"
summary_output="${summary_output:-}"

### End: Modify this block to tailor your analysis ###

# Forward config variables (including external overrides) into container environment.
forward_config_vars_to_container_env() {
  local job_script=$1
  local var_name
  local var_names
  var_names=$(
    sed -n '/^### Start: Modify this block to tailor your analysis ###$/,/^### End: Modify this block to tailor your analysis ###$/p' "${job_script}" \
      | sed -n -E 's/^[[:space:]]*([A-Za-z_][A-Za-z0-9_]*)=.*/\1/p' \
      | sort -u
  )
  while IFS= read -r var_name; do
    if [[ -z "${var_name}" ]]; then
      continue
    fi
    if [[ -n "${!var_name+x}" ]]; then
      export "${var_name}"
      export "SINGULARITYENV_${var_name}=${!var_name}"
      export "APPTAINERENV_${var_name}=${!var_name}"
    fi
  done <<< "${var_names}"

  # gg_debug_mode is read by gg_*_core.sh scripts but defined outside the config block.
  if [[ -n "${gg_debug_mode+x}" ]]; then
    export gg_debug_mode
    export "SINGULARITYENV_gg_debug_mode=${gg_debug_mode}"
    export "APPTAINERENV_gg_debug_mode=${gg_debug_mode}"
  fi

  # Forward GG_INPUT_* runtime overrides consumed by gg_input_generation_core.sh.
  local gg_input_var_name
  for gg_input_var_name in ${!GG_INPUT_@}; do
    export "${gg_input_var_name}"
    export "SINGULARITYENV_${gg_input_var_name}=${!gg_input_var_name}"
    export "APPTAINERENV_${gg_input_var_name}=${!gg_input_var_name}"
  done

  # Optional env-driven defaults consumed by gg_input_generation_core.sh.
  if [[ -n "${GG_INPUTPREP_CONFIG+x}" ]]; then
    export GG_INPUTPREP_CONFIG
    export "SINGULARITYENV_GG_INPUTPREP_CONFIG=${GG_INPUTPREP_CONFIG}"
    export "APPTAINERENV_GG_INPUTPREP_CONFIG=${GG_INPUTPREP_CONFIG}"
  fi
  if [[ -n "${GG_DATASET_ROOT+x}" ]]; then
    export GG_DATASET_ROOT
    export "SINGULARITYENV_GG_DATASET_ROOT=${GG_DATASET_ROOT}"
    export "APPTAINERENV_GG_DATASET_ROOT=${GG_DATASET_ROOT}"
  fi
}

forward_config_vars_to_container_env "${BASH_SOURCE[0]}"
unset -f forward_config_vars_to_container_env

source "${dir_script}/support/gg_util.sh" # loading utility functions
gg_scheduler_runtime_prelude
unset_singularity_envs
if ! set_singularity_command; then
  exit 1
fi
variable_SGEnizer
set_singularityenv
gg_print_scheduler_runtime_summary

if [[ -n "${GG_INPUT_DATASET_HOST_PATH:-}" ]]; then
  if [[ ! -d "${GG_INPUT_DATASET_HOST_PATH}" ]]; then
    echo "GG_INPUT_DATASET_HOST_PATH was set but directory was not found: ${GG_INPUT_DATASET_HOST_PATH}"
    exit 1
  fi
  dataset_host_path=$(cd "$(dirname "${GG_INPUT_DATASET_HOST_PATH}")" && pwd -P)/$(basename "${GG_INPUT_DATASET_HOST_PATH}")
  gg_add_container_bind_mount "${dataset_host_path}:/external/gfe_dataset"
  export SINGULARITYENV_GG_DATASET_ROOT="/external/gfe_dataset"
  export APPTAINERENV_GG_DATASET_ROOT="/external/gfe_dataset"
  echo "Binding external dataset root: ${dataset_host_path} -> /external/gfe_dataset"
fi

mkdir -p "${dir_pg}"
cd "${dir_pg}"
${singularity_command} "${gg_image}" < "${dir_script}/core/gg_input_generation_core.sh"
if ! gg_trigger_versions_dump "$(basename "${BASH_SOURCE[0]}")"; then
  echo "Warning: gg_versions trigger failed."
fi

echo "$(date): Ending"
