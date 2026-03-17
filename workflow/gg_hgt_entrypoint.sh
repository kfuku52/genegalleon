#!/usr/bin/env bash

# Scheduler header notes:
# - Keep sections ordered as SLURM -> UGE -> PBS across entrypoints.
# - Update job name, CPU count, memory, walltime, log paths, and array size together.
# - Site-specific partition/queue/resource lines stay commented examples by default.

# SLURM
# Common parameters: job name, cores per task, memory per core, walltime, log files, and working directory.
#SBATCH -J gg_hgt
#SBATCH -c 1
#SBATCH --mem-per-cpu=32G
#SBATCH -t 2976:00:00
#SBATCH --output=gg_hgt_entrypoint.sh_%j.out
#SBATCH --error=gg_hgt_entrypoint.sh_%j.err
#SBATCH --chdir=.
#SBATCH --ignore-pbs
# Site-specific partition example.
#SBATCH -p epyc
# Optional notifications and single-node examples.
##SBATCH --mail-type=ALL
##SBATCH --mail-user=<aaa@bbb.com>

## UGE
# Common parameters: shell, working directory, slot count, memory per slot, and runtime limits.
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 1
#$ -l s_vmem=32G
#$ -l mem_req=32G
# Site-specific resource example.
#$ -l epyc
#$ -l d_rt=124:00:00:00
#$ -l s_rt=124:00:00:00

## PBS
# Common parameters: shell, CPU count, total memory, and exported environment.
#PBS -S /bin/bash
#PBS -l ncpus=1
#PBS -l mem=32G
# Site-specific queue example.
#PBS -q small
#PBS -V

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
gg_entrypoint_name="gg_hgt_entrypoint.sh"

### Start: Modify this block to tailor your analysis ###

run_hgt_eval=1 # Score GeneRax-first HGT candidates from gg_orthogroup.db.
hgt_min_branch_score="0.45" # Minimum composite HGT score retained in the branch candidate table.
hgt_use_taxonomy_db=1 # Resolve UniProt best-hit taxonomic distances with the local ETE taxonomy DB when available.
hgt_contamination_dir="" # Optional directory containing species_cds_contamination_removal_tsv files; empty auto-detects the workspace default.

### End: Modify this block to tailor your analysis ###

source "${gg_support_dir}/gg_util.sh"
forward_config_vars_to_container_env "${gg_entrypoint_name}"
if ! gg_entrypoint_prepare_container_runtime 0; then
  exit 1
fi
gg_entrypoint_activate_container_runtime

gg_entrypoint_enter_workspace
gg_run_container_shell_script "${gg_container_image_path}" "${gg_core_dir}/gg_hgt_core.sh"
gg_require_versions_dump "${gg_entrypoint_name}"

echo "$(date): Ending"
