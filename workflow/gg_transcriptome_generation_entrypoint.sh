#!/usr/bin/env bash

# Scheduler header notes:
# - Keep sections ordered as SLURM -> UGE -> PBS across entrypoints.
# - Update job name, CPU count, memory, walltime, log paths, and array size together.
# - Site-specific partition/queue/resource lines stay commented examples by default.

# SLURM
# Common parameters: job name, cores per task, memory per core, walltime, log files, and working directory.
#SBATCH -J gg_transcriptome_generation
#SBATCH -c 4
#SBATCH --mem-per-cpu=32G
#SBATCH -t 2976:00:00
#SBATCH --output=gg_transcriptome_generation_entrypoint.sh_%A_%a.out
#SBATCH --error=gg_transcriptome_generation_entrypoint.sh_%A_%a.err
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
#$ -pe def_slot 4
#$ -l s_vmem=32G
#$ -l mem_req=32G
# Site-specific resource example.
#$ -l epyc
#$ -l d_rt=124:00:00:00
#$ -l s_rt=124:00:00:00
# Array example for array-aware entrypoints.
#$ -t 1

## PBS
# Common parameters: shell, CPU count, total memory, and exported environment.
#PBS -S /bin/bash
#PBS -l ncpus=4
#PBS -l mem=128G
# Array example for array-aware entrypoints.
#PBS -J 1
# Site-specific queue example.
#PBS -q small
#PBS -V

# Number of parallel batch jobs ("-t 1-N" in SGE or "--array 1-N" in SLURM):
# N = Number of directories in workspace/input/species_rnaseq for mode_transcriptome_assembly=fastq
# N = Number of files in workspace/input/query_sra_id for mode_transcriptome_assembly=sraid
# N = Number of files in workspace/input/amalgkit_metadata for mode_transcriptome_assembly=metadata
# In many cases, 4 cores and mem_req=16G (=64G) worked well in per-SRA Trinity assembly.
# Increase per-core RAM to 32G if Trinity fails.

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
gg_entrypoint_name="gg_transcriptome_generation_entrypoint.sh"

### Start: Modify this block to tailor your analysis ###

# Mode
mode_transcriptome_assembly="${mode_transcriptome_assembly:-sraid}" # {"auto", "sraid", "fastq", "metadata"}

# Workflow flags
run_amalgkit_metadata_or_integrate=1 # Metadata retrieval.
run_amalgkit_getfastq=1 # fastq generation from NCBI SRA.
run_assembly=1 # Transcriptome assembly with Trinity or rnaSPAdes.
run_longestcds=1 # Longest CDS extraction.
run_longestcds_fx2tab=1 # Sequence stats for longest CDS.
run_longestcds_mmseqs2taxonomy=0 # MMseqs2 taxonomy.
run_longestcds_contamination_removal=0 # Contamination removal.
run_busco_isoforms=1 # BUSCO for transcriptome assembly with isoforms.
run_busco_longest_cds=1 # BUSCO for longest CDS.
run_busco_contamination_removed_longest_cds=0 # BUSCO for contamination-removed longest CDS.
run_assembly_stat=1 # seqkit stat for assembly and CDS files.
run_amalgkit_quant=1 # Expression quantification.
run_amalgkit_merge=1 # Expression merge.
run_multispecies_summary=1 # Multi-species summary.

# Input-download parameters
amalgkit_rrna_filter="yes" # read-level rRNA removal in amalgkit getfastq. Observed to finish for ~4.2 Gbp total / ~14 million reads within a 350G job before contamination filtering; exact peak RAM and elapsed time for this step alone were not logged.
amalgkit_contam_filter="no" # read-level contamination removal in amalgkit getfastq. Rank follows contamination_removal_rank below. Setting yes may require >350G RAM for large private FASTQ inputs (for example, ~4.2 Gbp total, ~14 million reads).
amalgkit_metadata_max_concurrent_jobs="${amalgkit_metadata_max_concurrent_jobs:-10}" # Maximum number of concurrent array tasks allowed to call NCBI-backed amalgkit metadata; set 0 to disable gg-side throttling.
amalgkit_getfastq_max_concurrent_jobs="${amalgkit_getfastq_max_concurrent_jobs:-10}" # Maximum number of concurrent array tasks allowed to call NCBI-backed amalgkit getfastq or fallback FASTQ recovery; set 0 to disable gg-side throttling.
amalgkit_sra_strategy_query="${amalgkit_sra_strategy_query:-\"RNA-seq\"[Strategy] OR \"EST\"[Strategy] OR \"CLONE\"[Strategy]}" # Entrez strategy clause appended in mode_transcriptome_assembly=sraid; include CLONE so capillary/Sanger cDNA libraries are eligible. Set empty to disable strategy filtering.
amalgkit_long_read_instrument_pattern="${amalgkit_long_read_instrument_pattern:-pacbio|smrt|nanopore|minion|gridion|promethion|flongle|sequel|revio}" # Case-insensitive regex used to exclude long-read instruments after amalgkit metadata retrieval; leave empty to disable post-filtering.
remove_amalgkit_fastq_after_completion=1 # Delete per-species amalgkit FASTQ files after downstream completion.

# Assembly and quantification parameters
max_assembly_input_fastq_size="30,000,000,000" # Maximum total FASTQ length in bp used for transcriptome assembly.
assembly_method="rnaSPAdes" # {Trinity,rnaSPAdes}
protocol_rna_seq="mixed" # {same,mixed}
kallisto_reference="longest_cds" # {species_cds,longest_transcript,longest_cds,contamination_removed_longest_cds}
orf_aggregation_level="i" # {c,g,i,p}
assembly_cpu_offset=0 # Number of CPU cores reserved from GG_TASK_CPUS before launching the assembler.
assembly_ram_offset=4 # Amount of RAM in GB reserved from GG_MEM_TOTAL_GB before launching the assembler.

# Contamination-removal parameters
contamination_removal_rank="domain" # Taxonomic rank for contamination removal. Canonical value is domain; GeneGalleon normalizes tool-specific synonyms automatically.
contamination_removal_target_taxon="${contamination_removal_target_taxon:-}" # Optional NCBI taxon name used as the lineage anchor for contamination removal (for example, Eukaryota when the sample species name is unknown).

### End: Modify this block to tailor your analysis ###

delete_tmp_dir=1 # After this run, delete tmp directory created for each job. Set 0 when debugging.

source "${gg_support_dir}/gg_util.sh" # loading utility functions
# Forward config variables (including external overrides) into container environment.
forward_config_vars_to_container_env "${gg_entrypoint_name}" "delete_tmp_dir"
if ! gg_entrypoint_prepare_container_runtime 0; then
  exit 1
fi
gg_entrypoint_activate_container_runtime

gg_entrypoint_enter_workspace
gg_run_container_shell_script "${gg_container_image_path}" "${gg_core_dir}/gg_transcriptome_generation_core.sh"
gg_require_versions_dump "${gg_entrypoint_name}"

echo "$(date): Ending"
