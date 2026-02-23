#!/usr/bin/env bash
set -euo pipefail

# SLURM in NIG supercomputer
#SBATCH -J gg_gene_convergence
#SBATCH -c 8 # Number of CPUs
#SBATCH --mem-per-cpu=32G # RAM per CPU in GB
#SBATCH -t 60-00:00:00 # maximum time in d-hh:mm:ss format
#SBATCH --output=gg_gene_convergence_%j.out
#SBATCH --error=gg_gene_convergence_%j.err
#SBATCH -p medium # partition name, cluster environment specific
#SBATCH --chdir=.
##SBATCH --mail-type=ALL
##SBATCH --mail-user=<aaa@bbb.com>
##SBATCH -N 1-1 # Number of nodes
##SBATCH -n 1 # Number of tasks

## UGE
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 8
#$ -l s_vmem=32G
#$ -l mem_req=32G
#$ -l medium-ubuntu
#$ -l d_rt=60:00:00:00
#$ -l s_rt=60:00:00:00

## PBS in BIAS at NIBB
##PBS -S /bin/bash
##PBS -l ncpus=8
##PBS -l mem=256G
##PBS -q small
##PBS -V

# Number of parallel batch jobs ("-t 1-N" in SGE or "--array 1-N" in SLURM):
# Fixed to 1

echo "$(date): Starting"

# Change these directories for your custom-made analysis
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
dir_pg="${script_dir}/../workspace" # pg input and output directory
dir_script="${script_dir}" # directory where core/gg_*_core.sh and support/ locate
gg_image="${script_dir}/../genegalleon.sif" # path to the singularity image

### Start: Modify this block to tailor your analysis ###

arity_range="${arity_range:-2-10}"
trait="${trait:-all}"
skip_lower_order="${skip_lower_order:-yes}"
min_fg_stem_ratio="${min_fg_stem_ratio:-0.5}"
min_OCNany2spe="${min_OCNany2spe:-1.8}"
min_omegaCany2spe="${min_omegaCany2spe:-3.0}"
min_OCNCoD="${min_OCNCoD:-0}"
max_per_K="${max_per_K:-100}"
file_trait="${file_trait:-/workspace/input/species_trait/species_trait.tsv}"
dir_orthogroup="${dir_orthogroup:-/workspace/output/orthogroup}"
dir_orthofinder="${dir_orthofinder:-/workspace/output/orthofinder}"
dir_out="${dir_out:-/workspace/output/csubst_site}"

### End: Modify this block to tailor your analysis ###

# Misc
exit_if_running=0

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

  # Optional plotting/runtime env overrides consumed by gg_gene_convergence_core.sh.
  if [[ -n "${PYMOL_HEADLESS+x}" ]]; then
    export PYMOL_HEADLESS
    export "SINGULARITYENV_PYMOL_HEADLESS=${PYMOL_HEADLESS}"
    export "APPTAINERENV_PYMOL_HEADLESS=${PYMOL_HEADLESS}"
  fi
  if [[ -n "${QT_QPA_PLATFORM+x}" ]]; then
    export QT_QPA_PLATFORM
    export "SINGULARITYENV_QT_QPA_PLATFORM=${QT_QPA_PLATFORM}"
    export "APPTAINERENV_QT_QPA_PLATFORM=${QT_QPA_PLATFORM}"
  fi
}

forward_config_vars_to_container_env "${BASH_SOURCE[0]}"
unset -f forward_config_vars_to_container_env

source "${dir_script}/support/gg_util.sh" # loading utility functions
gg_scheduler_runtime_prelude
unset_singularity_envs
exit_if_running_qstat
if ! set_singularity_command; then
  exit 1
fi
variable_SGEnizer
set_singularityenv
gg_print_scheduler_runtime_summary

mkdir -p "${dir_pg}"
cd "${dir_pg}"
${singularity_command} "${gg_image}" < "${dir_script}/core/gg_gene_convergence_core.sh"
if ! gg_trigger_versions_dump "$(basename "${BASH_SOURCE[0]}")"; then
  echo "Warning: gg_versions trigger failed."
fi

echo "$(date): Ending"
