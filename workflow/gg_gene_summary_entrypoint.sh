#!/usr/bin/env bash

# Scheduler header notes:
# - Keep sections ordered as SLURM -> UGE -> PBS across entrypoints.
# - Update job name, CPU count, memory, walltime, log paths, and working directory together.
# - Site-specific partition/queue/resource lines stay commented examples by default.

# SLURM
# Common parameters: job name, cores per task, memory per core, walltime, log files, and working directory.
#SBATCH -J gg_gene_summary
#SBATCH -c 4
#SBATCH --mem-per-cpu=16G
#SBATCH -t 2976:00:00
#SBATCH --output=gg_gene_summary_entrypoint.sh_%j.out
#SBATCH --error=gg_gene_summary_entrypoint.sh_%j.err
#SBATCH --chdir=.
#SBATCH --ignore-pbs
# Site-specific partition example.
#SBATCH -p epyc,rome,medium
# Optional notifications and single-node examples.
##SBATCH --mail-type=ALL
##SBATCH --mail-user=<aaa@bbb.com>

## UGE
# Common parameters: shell, working directory, slot count, memory per slot, and runtime limits.
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 4
#$ -l s_vmem=16G
#$ -l mem_req=16G
# Site-specific resource example.
#$ -l epyc
#$ -l d_rt=124:00:00:00
#$ -l s_rt=124:00:00:00

## PBS
# Common parameters: shell, CPU count, total memory, and exported environment.
#PBS -S /bin/bash
#PBS -l ncpus=4
#PBS -l mem=64G
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
gg_entrypoint_name="gg_gene_summary_entrypoint.sh"

### Start: Modify this block to tailor your analysis ###

# Mode
mode_gene_summary="${mode_gene_summary:-query2family}" # query2family|orthogroup

# Output and shared summary parameters
dir_gene_summary="${dir_gene_summary:-auto}" # Output directory for mode-specific gene summaries.
run_gene_completion_summary="${run_gene_completion_summary:-1}" # Generate per-family completion/alignment summary tables.

# Gene-family species-tree presence/absence summary
run_gene_presence_absence="${run_gene_presence_absence:-1}" # Generate species x gene-family presence/absence and copy-number matrices.
gene_summary_species_tree="${gene_summary_species_tree:-auto}" # Species tree path, or auto to use query2family/species_tree outputs.
gene_summary_species_tree_ci="${gene_summary_species_tree_ci:-auto}" # Dated species-tree CI path, or auto to use mcmctree_95CI.nwk when available.
gene_summary_species_tree_support="${gene_summary_species_tree_support:-auto}" # Species tree with numeric branch-support labels, or auto to use species_tree outputs.
gene_summary_busco_table="${gene_summary_busco_table:-auto}" # BUSCO full-table directory or species summary table for the right-side stacked bars, or auto to use workspace outputs.
gene_summary_include_incomplete="${gene_summary_include_incomplete:-0}" # Include query files without stat_branch as NA columns.
gene_summary_heatmap_value="${gene_summary_heatmap_value:-presence}" # presence|copy_number for the PDF/SVG heatmap fill.
gene_summary_plot_width="${gene_summary_plot_width:-7.2}" # Width in inches for the presence/absence summary plot.
gene_summary_max_families="${gene_summary_max_families:-auto}" # Maximum plotted families; auto means all query2family queries and first 100 orthogroups. Use 0/all for no cap.
gene_summary_family_ids="${gene_summary_family_ids:-}" # Optional comma/space-separated family IDs to plot, overriding the default subset.
gene_summary_family_file="${gene_summary_family_file:-}" # Optional file listing family IDs to plot, one ID per line or first column.

# Database summary
run_database_prep="${run_database_prep:-0}" # Prepare or refresh gg_orthogroup.db for the selected mode.

# HGT summary
run_hgt_eval="${run_hgt_eval:-0}" # Summarize GeneRax-first HGT candidate evidence from gg_orthogroup.db.
run_hgt_plot="${run_hgt_plot:-0}" # Generate HGT overview, taxonomy-flow, and per-family HGT tree plots.
hgt_use_taxonomy_db="${hgt_use_taxonomy_db:-1}" # Resolve UniProt best-hit taxonomic distances with the local ETE taxonomy DB when available.
hgt_contamination_dir="${hgt_contamination_dir:-}" # Optional directory containing species_cds_contamination_removal_tsv files; empty auto-detects the workspace default.
hgt_taxonomy_flow_rank="${hgt_taxonomy_flow_rank:-phylum}" # Taxonomic rank used to collapse recipient/best-hit lineages in the flow plot.
hgt_taxonomy_flow_max_categories="${hgt_taxonomy_flow_max_categories:-12}" # Maximum recipient and best-hit categories retained before collapsing to Other.
hgt_tree_plot_width="${hgt_tree_plot_width:-24}" # Width in inches for HGT-specific tree plot PDFs.
hgt_promoter_bp="${hgt_promoter_bp:-2000}" # Promoter length used when re-rendering FIMO panels in HGT tree plots.
hgt_fimo_qvalue="${hgt_fimo_qvalue:-0.05}" # FIMO q-value threshold used when re-rendering HGT tree plots.
dir_hgt="${dir_hgt:-auto}" # HGT output directory; auto is mode-specific.

# Convergent-site summary
run_convergent_sites="${run_convergent_sites:-0}" # Run site-level convergence screening for the selected mode.
arity_range="${arity_range:-2-10}" # Branch-combination arity range to scan.
trait="${trait:-all}" # Trait column name(s) to analyze, or "all".
skip_lower_order="${skip_lower_order:-yes}" # Skip lower-order combinations once higher-order hits are retained.
min_fg_stem_ratio="${min_fg_stem_ratio:-0.5}" # Minimum fraction of foreground branches that must be stem branches.
min_OCNany2spe="${min_OCNany2spe:-1.8}" # Minimum OCNany2spe cutoff for candidate combinations.
min_omegaCany2spe="${min_omegaCany2spe:-3.0}" # Minimum omegaCany2spe cutoff for candidate combinations.
min_OCNCoD="${min_OCNCoD:-0}" # Minimum OCNCoD cutoff for candidate combinations.
max_per_K="${max_per_K:-100}" # Maximum number of combinations retained per arity K.
file_trait="${file_trait:-auto}" # Trait table path, or auto to use the workspace default.
dir_orthofinder="${dir_orthofinder:-auto}" # OrthoFinder result directory, or auto to detect it from the workspace.
dir_convergent_sites="${dir_convergent_sites:-auto}" # Output directory for convergence results; auto is mode-specific.

### End: Modify this block to tailor your analysis ###

source "${gg_support_dir}/gg_util.sh"
gg_apply_registered_env_overrides "${gg_entrypoint_name}"
forward_config_vars_to_container_env "${gg_entrypoint_name}"
gg_export_var_to_container_env_if_set "PYMOL_HEADLESS"
gg_export_var_to_container_env_if_set "QT_QPA_PLATFORM"
if ! gg_entrypoint_prepare_container_runtime 0; then
  exit 1
fi
gg_entrypoint_activate_container_runtime

gg_entrypoint_enter_workspace
gg_runtime_core_script="$(gg_prepare_entrypoint_runtime_snapshot "${gg_entrypoint_name}" "${gg_core_dir}/gg_gene_summary_core.sh")"
gg_run_container_shell_script "${gg_container_image_path}" "${gg_runtime_core_script}"
gg_require_versions_dump "${gg_entrypoint_name}"

echo "$(date): Ending"
