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

mode_sraid=1 # Need input at workspace/input/query_sra_id.
mode_fastq=0 # Need input at workspace/input/species_rnaseq.
mode_metadata=0 # Need input at workspace/input/amalgkit_metadata.

run_amalgkit_metadata_or_integrate=1 # Metadata retrieval.
run_amalgkit_getfastq=1 # fastq generation from NCBI SRA.
run_assembly=1 # Transcriptome assembly with Trinity or rnaSPAdes.
run_longestcds=1 # Longest CDS extraction.
run_longestcds_fx2tab=1 # Sequence stats for longest CDS.
run_longestcds_mmseqs2taxonomy=0 # MMseqs2 taxonomy.
run_longestcds_contamination_removal=0 # Contamination removal.
run_busco1=1 # BUSCO for transcriptome assembly with isoforms.
run_busco2=1 # BUSCO for longest CDS.
run_busco3=0 # BUSCO for contamination-removed longest CDS.
run_assembly_stat=1 # seqkit stat for assembly and CDS files.
run_amalgkit_quant=1 # Expression quantification.
run_amalgkit_merge=1 # Expression merge.
run_multispecies_summary=1 # Multi-species summary.

amalgkit_rrna_filter="yes" # read-level rRNA removal in amalgkit getfastq.
amalgkit_contam_filter="yes" # read-level contamination removal in amalgkit getfastq.
amalgkit_contam_filter_rank="phylum" # taxonomy rank for read-level contamination removal.
amalgkit_filter_order="fastp_first" # {fastp_first,rrna_first}
remove_amalgkit_fastq_after_completion=1
max_assembly_input_fastq_size="30,000,000,000"
assembly_method="rnaSPAdes" # {Trinity,rnaSPAdes}
protocol_rna_seq="mixed" # {same,mixed}
kallisto_reference="longest_cds" # {species_cds,longest_transcript,longest_cds,contamination_removed_longest_cds}
orf_aggregation_level="i" # {c,g,i,p}
assembly_cpu_offset=0
assembly_ram_offset=4

### End: Modify this block to tailor your analysis ###

delete_tmp_dir=1 # After this run, delete tmp directory created for each job. Set 0 when debugging.

source "${dir_script}/support/gg_util.sh" # loading utility functions
# Forward config variables (including external overrides) into container environment.
forward_config_vars_to_container_env "${BASH_SOURCE[0]}" "delete_tmp_dir"
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
