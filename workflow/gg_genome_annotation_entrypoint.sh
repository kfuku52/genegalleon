#!/usr/bin/env bash

# Scheduler header notes:
# - Keep sections ordered as SLURM -> UGE -> PBS across entrypoints.
# - Update job name, CPU count, memory, walltime, log paths, and array size together.
# - UGE resources default to SHIROKANE AGE; SLURM/PBS site-specific lines remain examples.

# SLURM
# Common parameters: job name, cores per task, memory per core, walltime, log files, and working directory.
#SBATCH -J gg_genome_annotation
#SBATCH -c 4
#SBATCH --mem-per-cpu=8G
#SBATCH -t 3-00:00:00
#SBATCH --output=gg_genome_annotation_entrypoint.sh_%A_%a.out
#SBATCH --error=gg_genome_annotation_entrypoint.sh_%A_%a.err
#SBATCH --chdir=.
#SBATCH --ignore-pbs
# Array example for array-aware entrypoints.
#SBATCH -a 1
# Site-specific partition example.
#SBATCH -p epyc,rome,medium
# Optional notifications and single-node examples.
##SBATCH --mail-type=ALL
##SBATCH --mail-user=<aaa@bbb.com>

## UGE
# SHIROKANE AGE defaults: shell, working directory, slot count, memory per slot, and ljob.
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 4
#$ -l s_vmem=8G
#$ -l ljob
# Array example for array-aware entrypoints.
#$ -t 1

## PBS
# Common parameters: shell, CPU count, total memory, and exported environment.
#PBS -S /bin/bash
#PBS -l ncpus=4
#PBS -l mem=32G
# Array example for array-aware entrypoints.
#PBS -J 1
# Site-specific queue example.
#PBS -q small
#PBS -V

# Number of parallel batch jobs ("-t 1-N" in SGE or "--array 1-N" in SLURM):
# N = Number of fasta files in workspace/input/species_cds

set -euo pipefail

echo "$(date): Starting"

# Resolve workflow paths for local and scheduler-spooled execution.
gg_bootstrap_script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
gg_bootstrap_checked_bases=""
for gg_bootstrap_base in \
  "${SLURM_SUBMIT_DIR:-}" \
  "${PBS_O_WORKDIR:-}" \
  "${PWD:-}" \
  "${gg_bootstrap_script_dir}"
do
  [[ -n "${gg_bootstrap_base}" ]] || continue
  case ":${gg_bootstrap_checked_bases}:" in
    *":${gg_bootstrap_base}:"*) continue ;;
  esac
  gg_bootstrap_checked_bases="${gg_bootstrap_checked_bases:+${gg_bootstrap_checked_bases}:}${gg_bootstrap_base}"
  for bootstrap_path in \
    "${gg_bootstrap_base}/support/gg_entrypoint_bootstrap.sh" \
    "${gg_bootstrap_base}/workflow/support/gg_entrypoint_bootstrap.sh"
  do
    if [[ -s "${bootstrap_path}" ]]; then
      # shellcheck disable=SC1090
      source "${bootstrap_path}"
      break
    fi
  done
  if declare -F gg_entrypoint_initialize >/dev/null 2>&1; then
    break
  fi
done
unset gg_bootstrap_base gg_bootstrap_checked_bases gg_bootstrap_script_dir bootstrap_path
if ! declare -F gg_entrypoint_initialize >/dev/null 2>&1; then
  echo "Failed to locate gg_entrypoint_bootstrap.sh from BASH_SOURCE[0]=${BASH_SOURCE[0]}" >&2
  exit 1
fi
if ! gg_entrypoint_initialize "${BASH_SOURCE[0]}" 1; then
  exit 1
fi
gg_entrypoint_name="gg_genome_annotation_entrypoint.sh"

### Start: Modify this block to tailor your analysis ###

# CDS workflow flags
run_collect_gff_info=0 # Collect per-gene coordinates, exon/intron structure, and gene stats from workspace/input/species_gff.
run_busco_cds=0 # BUSCO completeness analysis on species CDS inputs.
run_uniprot_annotation=0 # CDS annotation against UniProt Swiss-Prot.
run_cdskit_localize=0 # Predict targeting-peptide and peroxisome localization signals with cdskit localize.
run_cds_fx2tab=0 # Sequence-length and composition stats for species CDS FASTA files.
run_cds_mmseqs2taxonomy=0 # MMseqs2 taxonomy assignment for CDS sequences.
run_cds_contamination_removal=0 # Remove CDS sequences assigned outside the expected lineage.
run_annotation=0 # Combine GFF, UniProt, domain, taxonomy, and sequence-derived signals into per-gene annotation summaries.
run_wgd_ksd=0 # WGD inference by dS distribution

# Genome workflow flags
run_busco_genome=0 # BUSCO completeness analysis on genome assembly inputs.
run_subphaser=0 # Subgenome structure inference
run_genome_fx2tab=0 # Sequence-length and composition stats for genome assembly FASTA files.
run_scaffold_histogram=0 # Plot scaffold/contig length distributions for genome assemblies.
run_genome_mmseqs2taxonomy=0 # MMseqs2 taxonomy assignment for genome assembly sequences.
run_genome_contamination_removal=0 # Remove genome scaffolds/contigs assigned outside the expected lineage.
run_jcvi_dotplot=0 # Self-self synteny dotplot

# DNA-seq workflow flags
run_genomescope=0 # GenomeScope

# Summary workflow flags
run_multispecies_summary=1 # Multi-species summary plots and tables

# Annotation parameters
uniprot_annotation_method="mmseqs2" # blastp|mmseqs2 for UniProt Swiss-Prot annotation search engine.
cdskit_localize_model="${cdskit_localize_model:-targeting5-perox-deeploc21-et-v1}" # cdskit localize model path or alias; default includes the peroxisome head.
cdskit_localize_organism_group="${cdskit_localize_organism_group:-auto}" # auto|unknown|plant|non_plant; auto infers from GG_COMMON_BUSCO_LINEAGE/busco_lineage.
cdskit_localize_include_features=0 # Include internal cdskit localize sequence features in the output TSV.
cdskit_localize_no_model_download=0 # Set 1 to require the localize model to already exist in the cdskit model cache.

# Contamination-removal parameters
contamination_removal_rank="domain" # Taxonomic rank for contamination removal. Canonical value is domain; GeneGalleon normalizes tool-specific synonyms automatically.
contamination_removal_target_taxon="${contamination_removal_target_taxon:-}" # Optional NCBI taxon name used as the lineage anchor for contamination removal (for example, Eukaryota when the sample species name is unknown).

### End: Modify this block to tailor your analysis ###

# Misc
exit_if_running=0 # Exit without main analysis if the same GG_ARRAY_TASK_ID is already running.
delete_tmp_dir=1 # After this run, delete tmp directory created for each job. Set 0 when debugging.

source "${gg_support_dir}/gg_util.sh" # loading utility functions
# Forward config variables (including external overrides) into container environment.
gg_apply_registered_env_overrides "${gg_entrypoint_name}" "delete_tmp_dir"
forward_config_vars_to_container_env "${gg_entrypoint_name}" "delete_tmp_dir"
if ! gg_entrypoint_prepare_container_runtime 1; then
  exit 1
fi
gg_entrypoint_activate_container_runtime

gg_entrypoint_enter_workspace
gg_runtime_core_script="$(gg_prepare_entrypoint_runtime_snapshot "${gg_entrypoint_name}" "${gg_core_dir}/gg_genome_annotation_core.sh")"
gg_run_container_shell_script "${gg_container_image_path}" "${gg_runtime_core_script}"
gg_require_versions_dump "${gg_entrypoint_name}"

echo "$(date): Ending"
