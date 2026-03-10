#!/usr/bin/env bash

# SLURM
#SBATCH -J gg_progress_summary
#SBATCH -c 4
#SBATCH --mem-per-cpu=8G
#SBATCH -t 2976:00:00
#SBATCH --output=gg_progress_summary_%j.out
#SBATCH --error=gg_progress_summary_%j.err
#SBATCH -p epyc
#SBATCH --chdir=.
##SBATCH --mail-type=ALL
##SBATCH --mail-user=<aaa@bbb.com>
##SBATCH -N 1-1 # Number of nodes
##SBATCH -n 1 # Number of tasks

## UGE in NIG supercomputer
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 4
#$ -l s_vmem=8G
#$ -l mem_req=8G
#$ -l epyc
#$ -l d_rt=124:00:00:00
#$ -l s_rt=124:00:00:00

## PBS in BIAS at NIBB
##PBS -S /bin/bash
##PBS -l ncpus=4
##PBS -l mem=32G
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
gg_entrypoint_name="gg_progress_summary_entrypoint.sh"

### Start: Modify this block to tailor your analysis ###
mode_transcriptome_assembly="auto" # {"auto", "sraid", "fastq", "metadata"}
ncpu_progress_summary="" # Number of CPU threads used by summary scripts; empty falls back to GG_TASK_CPUS.
### End: Modify this block to tailor your analysis ###

source "${gg_support_dir}/gg_util.sh"
if ! gg_entrypoint_prepare_container_runtime 0; then
  exit 1
fi
: "${ncpu_progress_summary:=${GG_TASK_CPUS:-1}}"
gg_entrypoint_activate_container_runtime

forward_config_vars_to_container_env "${gg_entrypoint_name}"

gg_entrypoint_enter_workspace
gg_run_container_shell_script "${gg_container_image_path}" "${gg_core_dir}/gg_progress_summary_core.sh"
gg_require_versions_dump "${gg_entrypoint_name}"

echo "$(date): Ending"
