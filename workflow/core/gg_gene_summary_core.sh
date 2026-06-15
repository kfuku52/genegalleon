#!/usr/bin/env bash
set -euo pipefail

gg_core_bootstrap="/script/support/gg_core_bootstrap.sh"
if [[ ! -s "${gg_core_bootstrap}" ]]; then
  gg_core_bootstrap="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../support/gg_core_bootstrap.sh"
fi
# shellcheck disable=SC1090
source "${gg_core_bootstrap}"
unset gg_core_bootstrap

### Start: Job-supplied configuration ###
# Configuration variables are provided by gg_gene_summary_entrypoint.sh.
### End: Job-supplied configuration ###

gg_bootstrap_core_runtime "${BASH_SOURCE[0]:-$0}" "base" 1 1

mode_gene_summary=$(echo "${mode_gene_summary:-query2family}" | tr '[:upper:]' '[:lower:]')
run_gene_completion_summary="${run_gene_completion_summary:-1}"
run_gene_presence_absence="${run_gene_presence_absence:-1}"
run_database_prep="${run_database_prep:-0}"
run_hgt_eval="${run_hgt_eval:-0}"
run_hgt_plot="${run_hgt_plot:-0}"
run_convergent_sites="${run_convergent_sites:-0}"
gene_summary_include_incomplete="${gene_summary_include_incomplete:-0}"
gene_summary_heatmap_value=$(echo "${gene_summary_heatmap_value:-presence}" | tr '[:upper:]' '[:lower:]')
gene_summary_species_tree="${gene_summary_species_tree:-auto}"
gene_summary_species_tree_ci="${gene_summary_species_tree_ci:-auto}"
gene_summary_species_tree_support="${gene_summary_species_tree_support:-auto}"
gene_summary_busco_table="${gene_summary_busco_table:-auto}"
gene_summary_plot_width="${gene_summary_plot_width:-7.2}"
gene_summary_max_families="${gene_summary_max_families:-auto}"
gene_summary_family_ids="${gene_summary_family_ids:-}"
gene_summary_family_file="${gene_summary_family_file:-}"
dir_gene_summary="${dir_gene_summary:-auto}"
dir_hgt="${dir_hgt:-auto}"
dir_convergent_sites="${dir_convergent_sites:-auto}"
dir_orthofinder="${dir_orthofinder:-auto}"
file_trait="${file_trait:-auto}"

enable_all_run_flags_for_debug_mode

validate_binary_flag() {
  local name=$1
  local value=$2
  if [[ "${value}" != "0" && "${value}" != "1" ]]; then
    echo "Invalid binary flag value: ${name}=${value} (expected 0 or 1)"
    exit 1
  fi
}

validate_binary_flag "run_gene_completion_summary" "${run_gene_completion_summary}"
validate_binary_flag "run_gene_presence_absence" "${run_gene_presence_absence}"
validate_binary_flag "run_database_prep" "${run_database_prep}"
validate_binary_flag "run_hgt_eval" "${run_hgt_eval}"
validate_binary_flag "run_hgt_plot" "${run_hgt_plot}"
validate_binary_flag "run_convergent_sites" "${run_convergent_sites}"
validate_binary_flag "gene_summary_include_incomplete" "${gene_summary_include_incomplete}"

case "${mode_gene_summary}" in
  query2family)
    dir_gene_family="${gg_workspace_output_dir}/query2family"
    dir_query_gene="${gg_workspace_input_dir}/query_gene"
    default_dir_summary="${gg_workspace_output_dir}/gene_summary/query2family"
    default_dir_hgt="${gg_workspace_output_dir}/query2family_hgt"
    default_dir_convergent_sites="${gg_workspace_output_dir}/query2family_csubst_site"
    ;;
  orthogroup)
    dir_gene_family="${gg_workspace_output_dir}/orthogroup"
    default_dir_summary="${gg_workspace_output_dir}/gene_summary/orthogroup"
    default_dir_hgt="${gg_workspace_output_dir}/hgt"
    default_dir_convergent_sites="${gg_workspace_output_dir}/csubst_site"
    ;;
  *)
    echo "Invalid mode_gene_summary: ${mode_gene_summary}"
    echo 'mode_gene_summary must be either "query2family" or "orthogroup". Exiting.'
    exit 1
    ;;
esac

case "${gene_summary_heatmap_value}" in
  presence|copy_number)
    ;;
  *)
    echo "Invalid gene_summary_heatmap_value: ${gene_summary_heatmap_value}"
    echo 'gene_summary_heatmap_value must be either "presence" or "copy_number". Exiting.'
    exit 1
    ;;
esac

if [[ "${dir_gene_summary}" == "auto" ]]; then
  dir_gene_summary="${default_dir_summary}"
fi
if [[ "${dir_hgt}" == "auto" ]]; then
  dir_hgt="${default_dir_hgt}"
fi
if [[ "${dir_convergent_sites}" == "auto" ]]; then
  dir_convergent_sites="${default_dir_convergent_sites}"
fi

mkdir -p "${dir_gene_summary}"

file_gene_family_db="${dir_gene_family}/gg_orthogroup.db"
resolve_orthogroup_genecount_selected() {
  local candidates=(
    "${gg_workspace_output_dir}/orthofinder/Orthogroups_filtered/Orthogroups.GeneCount.selected.tsv"
    "${gg_workspace_output_dir}/orthofinder/Orthogroups/Orthogroups.GeneCount.selected.tsv"
    "${gg_workspace_output_dir}/orthofinder/hog2og/Orthogroups.GeneCount.selected.tsv"
    "${gg_workspace_output_dir}/orthofinder/Orthogroups.GeneCount.selected.tsv"
  )
  local candidate
  for candidate in "${candidates[@]}"; do
    if [[ -s "${candidate}" ]]; then
      printf '%s\n' "${candidate}"
      return 0
    fi
  done
  return 1
}

file_orthogroup_genecount_selected=""
if file_orthogroup_genecount_selected=$(resolve_orthogroup_genecount_selected); then
  :
else
  file_orthogroup_genecount_selected=""
fi

resolve_gene_summary_species_tree() {
  if [[ "${gene_summary_species_tree}" != "auto" ]]; then
    if [[ -s "${gene_summary_species_tree}" ]]; then
      printf '%s\n' "${gene_summary_species_tree}"
      return 0
    fi
    echo "Configured gene_summary_species_tree was not found: ${gene_summary_species_tree}" >&2
    return 1
  fi

  local candidates=(
    "${dir_gene_family}/parameters/dated_species_tree.pruned.nwk"
    "${dir_gene_family}/parameters/undated_species_tree.pruned.nwk"
    "${dir_gene_family}/parameters/dated_species_tree.nwk"
    "${dir_gene_family}/parameters/undated_species_tree.nwk"
    "${gg_workspace_output_dir}/species_tree/species_tree_summary/dated_species_tree.nwk"
    "${gg_workspace_output_dir}/species_tree/species_tree_summary/undated_species_tree.nwk"
    "${gg_workspace_output_dir}/species_tree/dated_species_tree.pruned.nwk"
    "${gg_workspace_output_dir}/species_tree/undated_species_tree.pruned.nwk"
    "${gg_workspace_output_dir}/species_tree/dated_species_tree.nwk"
    "${gg_workspace_output_dir}/species_tree/undated_species_tree.nwk"
  )
  local candidate
  for candidate in "${candidates[@]}"; do
    if [[ -s "${candidate}" ]]; then
      printf '%s\n' "${candidate}"
      return 0
    fi
  done
  return 1
}

tree_has_numeric_internal_labels() {
  local tree_file=$1
  grep -Eq '\)[[:space:]]*[-+]?[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?[[:space:]]*:' "${tree_file}"
}

resolve_gene_summary_species_tree_ci() {
  local file_species_tree=$1
  if [[ "${gene_summary_species_tree_ci}" != "auto" ]]; then
    if [[ -s "${gene_summary_species_tree_ci}" ]]; then
      printf '%s\n' "${gene_summary_species_tree_ci}"
      return 0
    fi
    echo "Configured gene_summary_species_tree_ci was not found: ${gene_summary_species_tree_ci}" >&2
    return 1
  fi

  case "${file_species_tree}" in
    *dated*|*mcmctree*)
      ;;
    *)
      return 1
      ;;
  esac

  local candidates=(
    "${dir_gene_family}/parameters/mcmctree_95CI.nwk"
    "${gg_workspace_output_dir}/species_tree/mcmctree_main/mcmctree_95CI.nwk"
  )
  local candidate
  for candidate in "${candidates[@]}"; do
    if [[ -s "${candidate}" ]]; then
      printf '%s\n' "${candidate}"
      return 0
    fi
  done
  return 1
}

resolve_gene_summary_species_tree_support() {
  if [[ "${gene_summary_species_tree_support}" != "auto" ]]; then
    if [[ -s "${gene_summary_species_tree_support}" ]]; then
      printf '%s\n' "${gene_summary_species_tree_support}"
      return 0
    fi
    echo "Configured gene_summary_species_tree_support was not found: ${gene_summary_species_tree_support}" >&2
    return 1
  fi

  local candidates=(
    "${gg_workspace_output_dir}/species_tree/concatenated_iqtree_dna/concat.dna.rooted.nwk"
    "${gg_workspace_output_dir}/species_tree/concatenated_iqtree_dna/concat.dna.iqtree.nwk"
    "${gg_workspace_output_dir}/species_tree/concatenated_iqtree_pep/concat.pep.rooted.nwk"
    "${gg_workspace_output_dir}/species_tree/concatenated_iqtree_pep/concat.pep.iqtree.nwk"
    "${gg_workspace_output_dir}/species_tree/single_copy_astral_dna/single_copy.astral.dna.pp1.nwk"
    "${gg_workspace_output_dir}/species_tree/single_copy_astral_pep/single_copy.astral.pep.pp1.nwk"
    "${dir_gene_family}/parameters/undated_species_tree.pruned.nwk"
    "${dir_gene_family}/parameters/undated_species_tree.nwk"
    "${gg_workspace_output_dir}/species_tree/species_tree_summary/undated_species_tree.nwk"
    "${gg_workspace_output_dir}/species_tree/undated_species_tree.pruned.nwk"
    "${gg_workspace_output_dir}/species_tree/undated_species_tree.nwk"
  )
  local candidate
  for candidate in "${candidates[@]}"; do
    if [[ -s "${candidate}" ]] && tree_has_numeric_internal_labels "${candidate}"; then
      printf '%s\n' "${candidate}"
      return 0
    fi
  done
  return 1
}

resolve_gene_summary_busco_table() {
  busco_summary_input_exists() {
    local path=$1
    if [[ -s "${path}" ]]; then
      return 0
    fi
    if [[ -d "${path}" ]] && [[ -n "$(find "${path}" -maxdepth 1 -type f -name "*busco.full.tsv" -print -quit 2> /dev/null)" ]]; then
      return 0
    fi
    return 1
  }

  if [[ "${gene_summary_busco_table}" != "auto" ]]; then
    if busco_summary_input_exists "${gene_summary_busco_table}"; then
      printf '%s\n' "${gene_summary_busco_table}"
      return 0
    fi
    echo "Configured gene_summary_busco_table was not found or contains no BUSCO full tables: ${gene_summary_busco_table}" >&2
    return 1
  fi

  local candidates=(
    "${gg_workspace_output_dir}/input_generation/species_cds_busco_full"
    "${gg_workspace_output_dir}/species_cds_busco_full"
    "${gg_workspace_output_dir}/species_protein_busco_full"
    "${gg_workspace_output_dir}/species_genome_busco_full"
    "${gg_workspace_output_dir}/species_tree/busco_summary_table/busco_summary.tsv"
    "${gg_workspace_output_dir}/genome_evolution/busco_summary_table/busco_summary.tsv"
    "${gg_workspace_output_dir}/annotation_summary/busco_cds.tsv"
  )
  local candidate
  for candidate in "${candidates[@]}"; do
    if busco_summary_input_exists "${candidate}"; then
      printf '%s\n' "${candidate}"
      return 0
    fi
  done
  return 1
}

run_completion_summary_for_mode() {
  if [[ ${run_gene_completion_summary} -ne 1 ]]; then
    echo "Skipping gene completion summary because run_gene_completion_summary=0."
    return 0
  fi
  if [[ "${mode_gene_summary}" == "query2family" ]]; then
    echo "Generating query2family completion summary."
    if [[ ! -d "${dir_gene_family}" ]]; then
      echo "Skipping query2family summary because directory was not found: ${dir_gene_family}"
      return 0
    fi
    if [[ ! -d "${dir_query_gene}" ]]; then
      echo "Skipping query2family summary because query_gene input directory was not found: ${dir_query_gene}"
      return 0
    fi
    python "${gg_support_dir}/query2family_output_summary.py" \
      --dir_query2family "${dir_gene_family}" \
      --dir_query_gene "${dir_query_gene}" \
      --ncpu "${GG_TASK_CPUS:-1}" \
      --out "${dir_gene_summary}/query2family_summary.tsv"
  else
    echo "Generating orthogroup completion summary."
    if [[ ! -d "${dir_gene_family}" ]]; then
      echo "Skipping orthogroup summary because directory was not found: ${dir_gene_family}"
      return 0
    fi
    if [[ -z "${file_orthogroup_genecount_selected}" || ! -s "${file_orthogroup_genecount_selected}" ]]; then
      echo "Skipping orthogroup summary because no selected gene-count table was found under workspace/output/orthofinder."
      return 0
    fi
    if [[ ! -d "${dir_gene_family}/amas_original" || ! -d "${dir_gene_family}/amas_cleaned" ]]; then
      echo "Skipping orthogroup summary because required AMAS directories were not found under: ${dir_gene_family}"
      return 0
    fi
    python "${gg_support_dir}/orthogroup_output_summary.py" \
      --dir_og "${dir_gene_family}" \
      --genecount "${file_orthogroup_genecount_selected}" \
      --ncpu "${GG_TASK_CPUS:-1}" \
      --out "${dir_gene_summary}/orthogroup_summary.tsv"
  fi
}

run_gene_presence_absence_summary() {
  if [[ ${run_gene_presence_absence} -ne 1 ]]; then
    echo "Skipping gene-family presence/absence summary because run_gene_presence_absence=0."
    return 0
  fi
  if [[ ! -d "${dir_gene_family}/stat_branch" ]]; then
    echo "Skipping gene-family presence/absence summary because stat_branch directory was not found: ${dir_gene_family}/stat_branch"
    return 0
  fi
  if [[ -z "$(find "${dir_gene_family}/stat_branch" -mindepth 1 -maxdepth 1 -type f -name '*_stat.branch.tsv' -print -quit 2>/dev/null)" ]]; then
    echo "Skipping gene-family presence/absence summary because no stat_branch files were found under: ${dir_gene_family}/stat_branch"
    return 0
  fi
  if [[ "${mode_gene_summary}" == "query2family" && ! -d "${dir_query_gene}" ]]; then
    echo "Skipping query2family presence/absence summary because query_gene input directory was not found: ${dir_query_gene}"
    return 0
  fi

  local file_species_tree
  if ! file_species_tree=$(resolve_gene_summary_species_tree); then
    echo "Skipping gene-family presence/absence summary because no species tree was found."
    return 0
  fi
  echo "Using species tree for ${mode_gene_summary} presence/absence summary: ${file_species_tree}"

  local file_species_tree_ci=""
  if file_species_tree_ci=$(resolve_gene_summary_species_tree_ci "${file_species_tree}"); then
    echo "Using dated species-tree CI annotations: ${file_species_tree_ci}"
  else
    file_species_tree_ci=""
  fi

  local file_species_tree_support=""
  if file_species_tree_support=$(resolve_gene_summary_species_tree_support); then
    echo "Using species-tree branch support labels: ${file_species_tree_support}"
  else
    file_species_tree_support=""
  fi

  local file_busco_table=""
  if file_busco_table=$(resolve_gene_summary_busco_table); then
    echo "Using BUSCO input for species stacked bars: ${file_busco_table}"
  else
    file_busco_table=""
  fi

  local outbase="${mode_gene_summary}_presence_absence"
  local file_presence="${dir_gene_summary}/${outbase}.tsv"
  local file_copy_number="${dir_gene_summary}/${mode_gene_summary}_copy_number.tsv"
  local file_long="${dir_gene_summary}/${outbase}.long.tsv"
  local file_plot_presence="${dir_gene_summary}/${outbase}.plot.tsv"
  local file_plot_copy_number="${dir_gene_summary}/${mode_gene_summary}_copy_number.plot.tsv"
  local file_plot_long="${dir_gene_summary}/${outbase}.plot.long.tsv"
  local file_selection="${dir_gene_summary}/${outbase}.plot_selection.tsv"
  local file_pdf="${dir_gene_summary}/${outbase}.pdf"
  local file_svg="${dir_gene_summary}/${outbase}.svg"

  local collect_args=(
    --mode "${mode_gene_summary}"
    --dir_gene_family "${dir_gene_family}"
    --species_tree "${file_species_tree}"
    --out_presence "${file_presence}"
    --out_copy_number "${file_copy_number}"
    --out_long "${file_long}"
    --out_plot_presence "${file_plot_presence}"
    --out_plot_copy_number "${file_plot_copy_number}"
    --out_plot_long "${file_plot_long}"
    --out_selection "${file_selection}"
    --include_incomplete "${gene_summary_include_incomplete}"
    --max_families "${gene_summary_max_families}"
  )
  if [[ "${mode_gene_summary}" == "query2family" ]]; then
    collect_args+=(--dir_query_gene "${dir_query_gene}")
  elif [[ -n "${file_orthogroup_genecount_selected}" ]]; then
    collect_args+=(--orthogroup_genecount "${file_orthogroup_genecount_selected}")
  fi
  if [[ -n "${gene_summary_family_ids}" ]]; then
    collect_args+=(--family_ids "${gene_summary_family_ids}")
  fi
  if [[ -n "${gene_summary_family_file}" ]]; then
    collect_args+=(--family_file "${gene_summary_family_file}")
  fi

  python "${gg_support_dir}/gene_family_presence_absence.py" "${collect_args[@]}"

  local plot_args=(
    --species_tree="${file_species_tree}" \
    --long_table="${file_plot_long}" \
    --value="${gene_summary_heatmap_value}" \
    --width="${gene_summary_plot_width}" \
    --out_pdf="${file_pdf}" \
    --out_svg="${file_svg}"
  )
  if [[ -n "${file_species_tree_ci}" ]]; then
    plot_args+=(--species_tree_ci="${file_species_tree_ci}")
  fi
  if [[ -n "${file_species_tree_support}" ]]; then
    plot_args+=(--support_tree="${file_species_tree_support}")
  fi
  if [[ -n "${file_busco_table}" ]]; then
    plot_args+=(--busco_table="${file_busco_table}")
  fi
  Rscript "${gg_support_dir}/plot_query2family_presence_absence.R" "${plot_args[@]}"
}

run_database_summary_for_mode() {
  if [[ ${run_database_prep} -ne 1 ]]; then
    echo "Skipping database prep because run_database_prep=0."
    return 0
  fi
  local missing_input=0
  for required_dir in \
    "${dir_gene_family}/stat_tree" \
    "${dir_gene_family}/stat_branch"
  do
    if [[ ! -d "${required_dir}" ]]; then
      echo "Skipping database prep because required directory is missing: ${required_dir}"
      missing_input=1
    fi
  done
  if [[ ${missing_input} -ne 0 ]]; then
    return 0
  fi
  echo "Generating gene-family database for mode_gene_summary=${mode_gene_summary}: ${file_gene_family_db}"
  python "${gg_support_dir}/generate_orthogroup_database.py" \
    --overwrite 1 \
    --dbpath "${file_gene_family_db}" \
    --dir_stat_tree "${dir_gene_family}/stat_tree" \
    --dir_stat_branch "${dir_gene_family}/stat_branch" \
    --dir_csubst_cb_prefix "${dir_gene_family}/csubst_cb_" \
    --row_threshold 8000 \
    --cutoff_stat "OCNany2spe,0.8" \
    --ncpu "${GG_TASK_CPUS:-1}"
}

run_hgt_summary_for_mode() {
  if [[ ${run_hgt_eval} -ne 1 && ${run_hgt_plot} -ne 1 ]]; then
    echo "Skipping HGT summary because run_hgt_eval=0 and run_hgt_plot=0."
    return 0
  fi
  echo "Running HGT summary for mode_gene_summary=${mode_gene_summary}."
  hgt_gene_family_dir="${dir_gene_family}" \
  hgt_db_path="${file_gene_family_db}" \
  hgt_output_dir="${dir_hgt}" \
  run_hgt_eval="${run_hgt_eval}" \
  run_hgt_plot="${run_hgt_plot}" \
  hgt_use_taxonomy_db="${hgt_use_taxonomy_db:-1}" \
  hgt_contamination_dir="${hgt_contamination_dir:-}" \
  hgt_taxonomy_flow_rank="${hgt_taxonomy_flow_rank:-phylum}" \
  hgt_taxonomy_flow_max_categories="${hgt_taxonomy_flow_max_categories:-12}" \
  hgt_tree_plot_width="${hgt_tree_plot_width:-24}" \
  hgt_promoter_bp="${hgt_promoter_bp:-2000}" \
  hgt_fimo_qvalue="${hgt_fimo_qvalue:-0.05}" \
    bash "${gg_core_dir}/gg_hgt_core.sh"
}

run_convergent_site_summary_for_mode() {
  if [[ ${run_convergent_sites} -ne 1 ]]; then
    echo "Skipping convergent-site summary because run_convergent_sites=0."
    return 0
  fi
  echo "Running convergent-site summary for mode_gene_summary=${mode_gene_summary}."
  dir_orthogroup="${dir_gene_family}" \
  dir_orthofinder="${dir_orthofinder}" \
  dir_out="${dir_convergent_sites}" \
  file_trait="${file_trait}" \
  arity_range="${arity_range:-2-10}" \
  trait="${trait:-all}" \
  skip_lower_order="${skip_lower_order:-yes}" \
  min_fg_stem_ratio="${min_fg_stem_ratio:-0.5}" \
  min_OCNany2spe="${min_OCNany2spe:-1.8}" \
  min_omegaCany2spe="${min_omegaCany2spe:-3.0}" \
  min_OCNCoD="${min_OCNCoD:-0}" \
  max_per_K="${max_per_K:-100}" \
    bash "${gg_core_dir}/gg_convergent_sites_core.sh"
}

echo "mode_gene_summary=${mode_gene_summary}"
echo "dir_gene_family=${dir_gene_family}"
echo "dir_gene_summary=${dir_gene_summary}"

run_completion_summary_for_mode
run_gene_presence_absence_summary
run_database_summary_for_mode
run_hgt_summary_for_mode
run_convergent_site_summary_for_mode

echo "$(date): Exiting Singularity environment"
