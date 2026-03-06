#!/usr/bin/env bash

# SLURM in NIG supercomputer
#SBATCH -J gg_gene_convergence
#SBATCH -c 8
#SBATCH --mem-per-cpu=32G
#SBATCH -t 60-00:00:00
#SBATCH --output=gg_gene_convergence_%j.out
#SBATCH --error=gg_gene_convergence_%j.err
#SBATCH -p medium
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

arity_range="2-10"
trait="all"
skip_lower_order="yes"
min_fg_stem_ratio="0.5"
min_OCNany2spe="1.8"
min_omegaCany2spe="3.0"
min_OCNCoD=0
max_per_K=100
file_trait="auto"
dir_orthogroup="auto"
dir_orthofinder="auto"
dir_out="auto"

### End: Modify this block to tailor your analysis ###

# Misc
exit_if_running=0

source "${gg_support_dir}/gg_util.sh" # loading utility functions
# Forward config variables (including external overrides) into container environment.
forward_config_vars_to_container_env "${BASH_SOURCE[0]}"
gg_export_var_to_container_env_if_set "PYMOL_HEADLESS"
gg_export_var_to_container_env_if_set "QT_QPA_PLATFORM"
if ! gg_entrypoint_prepare_container_runtime 1; then
  exit 1
fi
gg_entrypoint_activate_container_runtime

gg_entrypoint_enter_workspace
${singularity_command} "${gg_container_image_path}" < "${gg_core_dir}/gg_gene_convergence_core.sh"
if ! gg_trigger_versions_dump "$(basename "${BASH_SOURCE[0]}")"; then
  echo "Warning: gg_versions trigger failed."
fi

echo "$(date): Ending"
