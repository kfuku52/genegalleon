#!/usr/bin/env bash
set -euo pipefail

# SLURM in NIG supercomputer
#SBATCH -J gg_genome_annotation
#SBATCH -c 4 # Number of CPUs
#SBATCH --mem-per-cpu=8G # RAM per CPU in GB
#SBATCH -t 2976:00:00 # maximum time in d-hh:mm:ss format. NIG supercomputer epyc/medium MaxTime=2976:00:00
#SBATCH --output=gg_genome_annotation_%A_%a.out
#SBATCH --error=gg_genome_annotation_%A_%a.err
#SBATCH -p medium # partition name, cluster environment specific
#SBATCH --chdir=.
#SBATCH -a 1 # Array job, 1-N
##SBATCH --mail-type=ALL
##SBATCH --mail-user=<aaa@bbb.com>
##SBATCH -N 1-1 # Number of nodes
##SBATCH -n 1 # Number of tasks

#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 4
#$ -l s_vmem=8G
#$ -l mem_req=8G
#$ -l epyc
#$ -l d_rt=124:00:00:00
#$ -l s_rt=124:00:00:00
#$ -t 1

# Number of parallel batch jobs ("-t 1-N" in SGE or "--array 1-N" in SLURM):
# N = Number of fasta files in workspace/input/species_cds

echo "$(date): Starting"

# Change these directories for your custom-made analysis
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
dir_pg="${script_dir}/../workspace" # pg input and output directory
dir_script="${script_dir}" # directory where gg_util.sh and core/gg_*_core.sh locate
gg_image="${script_dir}/../genegalleon.sif" # path to the singularity image

# Load shared defaults when available.
if [[ -s "${script_dir}/gg_common_params.sh" ]]; then
  # shellcheck disable=SC1091
  source "${script_dir}/gg_common_params.sh"
fi

### Start: Modify this block to tailor your analysis ###

# CDS analysis
run_get_gff_info="${run_get_gff_info:-0}" # Collect gene information from workspace/input/species_gff
run_busco_cds="${run_busco_cds:-0}" # Gene set completeness analysis
run_uniprot_annotation="${run_uniprot_annotation:-0}" # DIAMOND-based CDS annotation against UniProt Swiss-Prot.
run_cds_fx2tab="${run_cds_fx2tab:-0}" # Sequence stats of CDS sequences
run_cds_mmseqs2taxonomy="${run_cds_mmseqs2taxonomy:-0}" # Taxonomic assignment of CDS sequences
run_cds_contamination_removal="${run_cds_contamination_removal:-0}" # Removal of contaminated sequences
run_annotation="${run_annotation:-0}" # Per-gene annotation summary
run_wgd_ksd="${run_wgd_ksd:-0}" # WGD inference by dS distribution

# Genome analysis
run_busco_genome="${run_busco_genome:-0}" # Gene set completeness analysis
run_subphaser="${run_subphaser:-0}" # Subgenome structure inference
run_genome_fx2tab="${run_genome_fx2tab:-0}" # Sequence stats of reference genome
run_scaffold_histogram="${run_scaffold_histogram:-0}" # Scaffold histogram
run_genome_mmseqs2taxonomy="${run_genome_mmseqs2taxonomy:-0}" # Taxonomic assignment of genome assembly
run_genome_contamination_removal="${run_genome_contamination_removal:-0}" # Removal of contaminated sequences
run_jcvi_dotplot="${run_jcvi_dotplot:-0}" # Self-self synteny dotplot

# DNA-seq analysis
run_genomescope="${run_genomescope:-0}" # GenomeScope

# Summary
run_multispecies_summary="${run_multispecies_summary:-1}" # Multi-species summary plots and tables

### End: Modify this block to tailor your analysis ###

# Misc
exit_if_running=0 # Exit without main analysis if the same SGE_TASK_ID is already running.
delete_tmp_dir="${delete_tmp_dir:-1}" # After this run, delete tmp directory created for each job. Set 0 when debugging.

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

  # Cleanup controls are defined outside the config block but consumed by core script.
  if [[ -n "${delete_tmp_dir+x}" ]]; then
    export delete_tmp_dir
    export "SINGULARITYENV_delete_tmp_dir=${delete_tmp_dir}"
    export "APPTAINERENV_delete_tmp_dir=${delete_tmp_dir}"
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
${singularity_command} "${gg_image}" < "${dir_script}/core/gg_genome_annotation_core.sh"
if ! gg_trigger_versions_dump "$(basename "${BASH_SOURCE[0]}")"; then
  echo "Warning: gg_versions trigger failed."
fi

echo "$(date): Ending"
