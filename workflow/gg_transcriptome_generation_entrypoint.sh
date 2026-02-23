#!/usr/bin/env bash
set -euo pipefail

# SLURM in NIG supercomputer
#SBATCH -J gg_transcriptome_generation
#SBATCH -c 4 # Number of CPUs
#SBATCH --mem-per-cpu=32G # RAM per CPU in GB
#SBATCH -t 2976:00:00 # maximum time in d-hh:mm:ss format. NIG supercomputer epyc/medium MaxTime=2976:00:00
#SBATCH --output=gg_transcriptome_generation_%A_%a.out
#SBATCH --error=gg_transcriptome_generation_%A_%a.err
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
#$ -l s_vmem=32G
#$ -l mem_req=32G
#$ -l epyc
#$ -l d_rt=124:00:00:00
#$ -l s_rt=124:00:00:00
#$ -t 1

# Number of parallel batch jobs ("-t 1-N" in SGE or "--array 1-N" in SLURM):
# N = Number of directories in workspace/input/species_rnaseq for mode_fastq=1
# N = Number of files in workspace/input/query_sra_id for mode_sraid=1

echo "$(date): Starting"

# In many cases, 4 cores and mem_req=16G (=64G) worked well in per-SRA Trinity assembly.
# Increase per-core RAM to 32G if Trinity fails.

# Change these directories for your custom-made analysis
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
dir_pg="${script_dir}/../workspace" # pg input and output directory
dir_script="${script_dir}" # directory core/gg_*_core.sh locate
gg_image="${script_dir}/../genegalleon.sif" # path to the singularity image

# Load shared defaults when available.
if [[ -s "${script_dir}/gg_common_params.sh" ]]; then
  # shellcheck disable=SC1091
  source "${script_dir}/gg_common_params.sh"
fi

### Start: Modify this block to tailor your analysis ###

mode_sraid="${mode_sraid:-1}" # Need input at workspace/input/query_sra_id.
mode_fastq="${mode_fastq:-0}" # Need input at workspace/input/species_rnaseq.
mode_metadata="${mode_metadata:-0}" # Need input at workspace/input/transcriptome_assembly/amalgkit_metadata.

run_amalgkit_metadata_or_integrate="${run_amalgkit_metadata_or_integrate:-1}" # Metadata retrieval.
run_amalgkit_getfastq="${run_amalgkit_getfastq:-1}" # fastq generation from NCBI SRA.
run_rRNA_contamination_report="${run_rRNA_contamination_report:-1}" # rRNA contamination report.
run_assembly="${run_assembly:-1}" # Transcriptome assembly with Trinity or rnaSPAdes.
run_longestcds="${run_longestcds:-1}" # Longest CDS extraction.
run_longestcds_fx2tab="${run_longestcds_fx2tab:-1}" # Sequence stats for longest CDS.
run_longestcds_mmseqs2taxonomy="${run_longestcds_mmseqs2taxonomy:-0}" # MMseqs2 taxonomy.
run_longestcds_contamination_removal="${run_longestcds_contamination_removal:-0}" # Contamination removal.
run_busco1="${run_busco1:-1}" # BUSCO for transcriptome assembly with isoforms.
run_busco2="${run_busco2:-1}" # BUSCO for longest CDS.
run_busco3="${run_busco3:-0}" # BUSCO for contamination-removed longest CDS.
run_assembly_stat="${run_assembly_stat:-1}" # seqkit stat for assembly and CDS files.
run_amalgkit_quant="${run_amalgkit_quant:-1}" # Expression quantification.
run_amalgkit_merge="${run_amalgkit_merge:-1}" # Expression merge.
run_multispecies_summary="${run_multispecies_summary:-1}" # Multi-species summary.

remove_amalgkit_fastq_after_completion="${remove_amalgkit_fastq_after_completion:-1}"
max_assembly_input_fastq_size="${max_assembly_input_fastq_size:-30,000,000,000}"
assembly_method="${assembly_method:-rnaSPAdes}" # {Trinity,rnaSPAdes}
protocol_rna_seq="${protocol_rna_seq:-mixed}" # {same,mixed}
kallisto_reference="${kallisto_reference:-longest_cds}" # {species_cds,longest_transcript,longest_cds,contamination_removed_longest_cds}
orf_aggregation_level="${orf_aggregation_level:-i}" # {c,g,i,p}
assembly_cpu_offset="${assembly_cpu_offset:-0}"
assembly_ram_offset="${assembly_ram_offset:-4}"

### End: Modify this block to tailor your analysis ###

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
if ! set_singularity_command; then
  exit 1
fi
variable_SGEnizer
set_singularityenv
gg_print_scheduler_runtime_summary

mkdir -p "${dir_pg}"
cd "${dir_pg}"
${singularity_command} "${gg_image}" < "${dir_script}/core/gg_transcriptome_generation_core.sh"
if ! gg_trigger_versions_dump "$(basename "${BASH_SOURCE[0]}")"; then
  echo "Warning: gg_versions trigger failed."
fi

echo "$(date): Ending"
