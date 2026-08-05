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

gene_family_source=$(echo "${gene_family_source:-query2family}" | tr '[:upper:]' '[:lower:]')
run_family_completion_summary="${run_family_completion_summary:-1}"
run_presence_absence_summary="${run_presence_absence_summary:-1}"
run_gene_family_database_build="${run_gene_family_database_build:-0}"
run_csubst_scan_aa_change_summary="${run_csubst_scan_aa_change_summary:-0}"
run_hgt_candidate_summary="${run_hgt_candidate_summary:-0}"
run_hgt_summary_plots="${run_hgt_summary_plots:-0}"
run_csubst_site_convergence_summary="${run_csubst_site_convergence_summary:-0}"
csubst_scan_aa_change_top_n="${csubst_scan_aa_change_top_n:-30}"
csubst_site_nonsyn_recode=$(echo "${csubst_site_nonsyn_recode:-${GG_COMMON_CSUBST_NONSYN_RECODE:-no}}" | tr '[:upper:]' '[:lower:]')
presence_absence_include_incomplete="${presence_absence_include_incomplete:-0}"
presence_absence_heatmap_value=$(echo "${presence_absence_heatmap_value:-presence}" | tr '[:upper:]' '[:lower:]')
presence_absence_species_tree="${presence_absence_species_tree:-auto}"
presence_absence_species_tree_ci="${presence_absence_species_tree_ci:-auto}"
presence_absence_species_tree_support="${presence_absence_species_tree_support:-auto}"
presence_absence_busco_table="${presence_absence_busco_table:-auto}"
presence_absence_plot_width="${presence_absence_plot_width:-7.2}"
presence_absence_max_families="${presence_absence_max_families:-auto}"
presence_absence_family_ids="${presence_absence_family_ids:-}"
presence_absence_family_file="${presence_absence_family_file:-}"
summary_output_dir="${summary_output_dir:-auto}"
hgt_summary_output_dir="${hgt_summary_output_dir:-auto}"
csubst_site_output_dir="${csubst_site_output_dir:-auto}"
csubst_site_orthofinder_dir="${csubst_site_orthofinder_dir:-auto}"
csubst_site_trait_file="${csubst_site_trait_file:-auto}"
csubst_site_arity_range="${csubst_site_arity_range:-2-10}"
csubst_site_trait="${csubst_site_trait:-all}"
csubst_site_skip_lower_order="${csubst_site_skip_lower_order:-yes}"
csubst_site_min_fg_stem_ratio="${csubst_site_min_fg_stem_ratio:-0.5}"
csubst_site_min_ocn_any2spe="${csubst_site_min_ocn_any2spe:-1.8}"
csubst_site_min_omega_c_any2spe="${csubst_site_min_omega_c_any2spe:-3.0}"
csubst_site_min_ocn_cod="${csubst_site_min_ocn_cod:-0}"
csubst_site_max_candidates_per_arity="${csubst_site_max_candidates_per_arity:-100}"

enable_all_run_flags_for_debug_mode

validate_binary_flag() {
  local name=$1
  local value=$2
  if [[ "${value}" != "0" && "${value}" != "1" ]]; then
    echo "Invalid binary flag value: ${name}=${value} (expected 0 or 1)"
    exit 1
  fi
}

validate_binary_flag "run_family_completion_summary" "${run_family_completion_summary}"
validate_binary_flag "run_presence_absence_summary" "${run_presence_absence_summary}"
validate_binary_flag "run_gene_family_database_build" "${run_gene_family_database_build}"
validate_binary_flag "run_csubst_scan_aa_change_summary" "${run_csubst_scan_aa_change_summary}"
validate_binary_flag "run_hgt_candidate_summary" "${run_hgt_candidate_summary}"
validate_binary_flag "run_hgt_summary_plots" "${run_hgt_summary_plots}"
validate_binary_flag "run_csubst_site_convergence_summary" "${run_csubst_site_convergence_summary}"
validate_binary_flag "presence_absence_include_incomplete" "${presence_absence_include_incomplete}"

case "${gene_family_source}" in
  query2family)
    dir_gene_family="${gg_workspace_output_dir}/query2family"
    dir_query_gene="${gg_workspace_input_dir}/query_gene"
    default_summary_output_dir="${gg_workspace_output_dir}/gene_summary/query2family"
    default_dir_hgt="${gg_workspace_output_dir}/query2family_hgt"
    default_dir_convergent_sites="${gg_workspace_output_dir}/query2family_csubst_site"
    ;;
  orthogroup)
    dir_gene_family="${gg_workspace_output_dir}/orthogroup"
    default_summary_output_dir="${gg_workspace_output_dir}/gene_summary/orthogroup"
    default_dir_hgt="${gg_workspace_output_dir}/hgt"
    default_dir_convergent_sites="${gg_workspace_output_dir}/csubst_site"
    ;;
  *)
    echo "Invalid gene_family_source: ${gene_family_source}"
    echo 'gene_family_source must be either "query2family" or "orthogroup". Exiting.'
    exit 1
    ;;
esac

case "${presence_absence_heatmap_value}" in
  presence|copy_number)
    ;;
  *)
    echo "Invalid presence_absence_heatmap_value: ${presence_absence_heatmap_value}"
    echo 'presence_absence_heatmap_value must be either "presence" or "copy_number". Exiting.'
    exit 1
    ;;
esac

if [[ "${summary_output_dir}" == "auto" ]]; then
  summary_output_dir="${default_summary_output_dir}"
fi
if [[ "${hgt_summary_output_dir}" == "auto" ]]; then
  hgt_summary_output_dir="${default_dir_hgt}"
fi
if [[ "${csubst_site_output_dir}" == "auto" ]]; then
  csubst_site_output_dir="${default_dir_convergent_sites}"
fi

mkdir -p "${summary_output_dir}"

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

resolve_presence_absence_species_tree() {
  if [[ "${presence_absence_species_tree}" != "auto" ]]; then
    if [[ -s "${presence_absence_species_tree}" ]]; then
      printf '%s\n' "${presence_absence_species_tree}"
      return 0
    fi
    echo "Configured presence_absence_species_tree was not found: ${presence_absence_species_tree}" >&2
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

resolve_presence_absence_species_tree_ci() {
  local file_species_tree=$1
  if [[ "${presence_absence_species_tree_ci}" != "auto" ]]; then
    if [[ -s "${presence_absence_species_tree_ci}" ]]; then
      printf '%s\n' "${presence_absence_species_tree_ci}"
      return 0
    fi
    echo "Configured presence_absence_species_tree_ci was not found: ${presence_absence_species_tree_ci}" >&2
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

resolve_presence_absence_species_tree_support() {
  if [[ "${presence_absence_species_tree_support}" != "auto" ]]; then
    if [[ -s "${presence_absence_species_tree_support}" ]]; then
      printf '%s\n' "${presence_absence_species_tree_support}"
      return 0
    fi
    echo "Configured presence_absence_species_tree_support was not found: ${presence_absence_species_tree_support}" >&2
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

resolve_presence_absence_busco_table() {
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

  if [[ "${presence_absence_busco_table}" != "auto" ]]; then
    if busco_summary_input_exists "${presence_absence_busco_table}"; then
      printf '%s\n' "${presence_absence_busco_table}"
      return 0
    fi
    echo "Configured presence_absence_busco_table was not found or contains no BUSCO full tables: ${presence_absence_busco_table}" >&2
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

run_family_completion_summary_for_source() {
  if [[ ${run_family_completion_summary} -ne 1 ]]; then
    echo "Skipping gene completion summary because run_family_completion_summary=0."
    return 0
  fi
  if [[ "${gene_family_source}" == "query2family" ]]; then
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
      --out "${summary_output_dir}/query2family_summary.tsv"
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
    python "${gg_support_dir}/orthogroup_output_summary.py" \
      --dir_og "${dir_gene_family}" \
      --genecount "${file_orthogroup_genecount_selected}" \
      --ncpu "${GG_TASK_CPUS:-1}" \
      --updated-genecount-out "${summary_output_dir}/orthogroup_genecount.amas.tsv" \
      --out "${summary_output_dir}/orthogroup_summary.tsv"
  fi
}

run_presence_absence_summary_for_source() {
  local store_status
  if [[ ${run_presence_absence_summary} -ne 1 ]]; then
    echo "Skipping gene-family presence/absence summary because run_presence_absence_summary=0."
    return 0
  fi
  if python "${gg_support_dir}/gene_family_output_store.py" has-files \
    --root "${dir_gene_family}" \
    --subdir stat_branch \
    --suffix "_stat.branch.tsv"
  then
    :
  else
    store_status=$?
    if [[ ${store_status} -eq 1 ]]; then
      echo "Skipping gene-family presence/absence summary because no live or ZIP-backed stat_branch files were found under: ${dir_gene_family}"
      return 0
    fi
    return "${store_status}"
  fi
  if [[ "${gene_family_source}" == "query2family" && ! -d "${dir_query_gene}" ]]; then
    echo "Skipping query2family presence/absence summary because query_gene input directory was not found: ${dir_query_gene}"
    return 0
  fi

  local file_species_tree
  if ! file_species_tree=$(resolve_presence_absence_species_tree); then
    echo "Skipping gene-family presence/absence summary because no species tree was found."
    return 0
  fi
  echo "Using species tree for ${gene_family_source} presence/absence summary: ${file_species_tree}"

  local file_species_tree_ci=""
  if file_species_tree_ci=$(resolve_presence_absence_species_tree_ci "${file_species_tree}"); then
    echo "Using dated species-tree CI annotations: ${file_species_tree_ci}"
  else
    file_species_tree_ci=""
  fi

  local file_species_tree_support=""
  if file_species_tree_support=$(resolve_presence_absence_species_tree_support); then
    echo "Using species-tree branch support labels: ${file_species_tree_support}"
  else
    file_species_tree_support=""
  fi

  local file_busco_table=""
  if file_busco_table=$(resolve_presence_absence_busco_table); then
    echo "Using BUSCO input for species stacked bars: ${file_busco_table}"
  else
    file_busco_table=""
  fi

  local outbase="${gene_family_source}_presence_absence"
  local file_presence="${summary_output_dir}/${outbase}.tsv"
  local file_copy_number="${summary_output_dir}/${gene_family_source}_copy_number.tsv"
  local file_long="${summary_output_dir}/${outbase}.long.tsv"
  local file_plot_presence="${summary_output_dir}/${outbase}.plot.tsv"
  local file_plot_copy_number="${summary_output_dir}/${gene_family_source}_copy_number.plot.tsv"
  local file_plot_long="${summary_output_dir}/${outbase}.plot.long.tsv"
  local file_selection="${summary_output_dir}/${outbase}.plot_selection.tsv"
  local file_pdf="${summary_output_dir}/${outbase}.pdf"
  local file_svg="${summary_output_dir}/${outbase}.svg"

  local collect_args=(
    --mode "${gene_family_source}"
    --dir_gene_family "${dir_gene_family}"
    --species_tree "${file_species_tree}"
    --out_presence "${file_presence}"
    --out_copy_number "${file_copy_number}"
    --out_long "${file_long}"
    --out_plot_presence "${file_plot_presence}"
    --out_plot_copy_number "${file_plot_copy_number}"
    --out_plot_long "${file_plot_long}"
    --out_selection "${file_selection}"
    --include_incomplete "${presence_absence_include_incomplete}"
    --max_families "${presence_absence_max_families}"
  )
  if [[ "${gene_family_source}" == "query2family" ]]; then
    collect_args+=(--dir_query_gene "${dir_query_gene}")
  elif [[ -n "${file_orthogroup_genecount_selected}" ]]; then
    collect_args+=(--orthogroup_genecount "${file_orthogroup_genecount_selected}")
  fi
  if [[ -n "${presence_absence_family_ids}" ]]; then
    collect_args+=(--family_ids "${presence_absence_family_ids}")
  fi
  if [[ -n "${presence_absence_family_file}" ]]; then
    collect_args+=(--family_file "${presence_absence_family_file}")
  fi

  python "${gg_support_dir}/gene_family_presence_absence.py" "${collect_args[@]}"

  local plot_args=(
    --species_tree="${file_species_tree}" \
    --long_table="${file_plot_long}" \
    --value="${presence_absence_heatmap_value}" \
    --width="${presence_absence_plot_width}" \
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

run_gene_family_database_for_source() {
  if [[ ${run_gene_family_database_build} -ne 1 ]]; then
    echo "Skipping database prep because run_gene_family_database_build=0."
    return 0
  fi
  local required_subdir
  local store_status
  for required_subdir in stat_tree stat_branch; do
    if python "${gg_support_dir}/gene_family_output_store.py" has-files \
      --root "${dir_gene_family}" \
      --subdir "${required_subdir}"
    then
      :
    else
      store_status=$?
      if [[ ${store_status} -eq 1 ]]; then
        echo "Skipping database prep because no live or ZIP-backed files were found in logical subdirectory: ${required_subdir}"
        return 0
      fi
      return "${store_status}"
    fi
  done
  echo "Generating gene-family database for gene_family_source=${gene_family_source}: ${file_gene_family_db}"
  python "${gg_support_dir}/generate_orthogroup_database.py" \
    --overwrite 1 \
    --dbpath "${file_gene_family_db}" \
    --dir_gene_family "${dir_gene_family}" \
    --dir_stat_tree "${dir_gene_family}/stat_tree" \
    --dir_stat_branch "${dir_gene_family}/stat_branch" \
    --dir_csubst_cb_prefix "${dir_gene_family}/csubst_cb_" \
    --dir_csubst_aa_change "${dir_gene_family}/csubst_scan" \
    --dir_csubst_aa_change_unit "${dir_gene_family}/csubst_scan_units" \
    --row_threshold 8000 \
    --cutoff_stat "OCNany2spe,0.8" \
    --ncpu "${GG_TASK_CPUS:-1}"
}

run_csubst_scan_aa_change_summary_for_source() {
  if [[ ${run_csubst_scan_aa_change_summary} -ne 1 ]]; then
    echo "Skipping CSUBST scan summary because run_csubst_scan_aa_change_summary=0."
    return 0
  fi
  if [[ ! -s "${file_gene_family_db}" ]]; then
    echo "Skipping CSUBST scan summary because gene-family database was not found: ${file_gene_family_db}"
    echo "Set run_gene_family_database_build=1 with run_csubst_scan_aa_change_summary=1 to refresh the database and generate CSUBST scan summaries in one run."
    return 0
  fi
  echo "Generating CSUBST scan AA-change summary for gene_family_source=${gene_family_source}: ${file_gene_family_db}"
  python "${gg_support_dir}/plot_csubst_aa_change_summary.py" \
    --dbpath "${file_gene_family_db}" \
    --out_prefix "${summary_output_dir}/${gene_family_source}_csubst_aa_change" \
    --out_tsv "${summary_output_dir}/${gene_family_source}_csubst_aa_change_summary.tsv" \
    --top_n "${csubst_scan_aa_change_top_n}"
}

run_hgt_summary_for_source() {
  if [[ ${run_hgt_candidate_summary} -ne 1 && ${run_hgt_summary_plots} -ne 1 ]]; then
    echo "Skipping HGT summary because run_hgt_candidate_summary=0 and run_hgt_summary_plots=0."
    return 0
  fi
  echo "Running HGT summary for gene_family_source=${gene_family_source}."
  hgt_gene_family_dir="${dir_gene_family}" \
  hgt_gene_family_mode="${gene_family_source}" \
  hgt_db_path="${file_gene_family_db}" \
  hgt_output_dir="${hgt_summary_output_dir}" \
  run_hgt_eval="${run_hgt_candidate_summary}" \
  run_hgt_plot="${run_hgt_summary_plots}" \
  hgt_use_taxonomy_db="${hgt_summary_use_taxonomy_db:-1}" \
  hgt_contamination_dir="${hgt_summary_contamination_dir:-}" \
  hgt_taxonomy_flow_rank="${hgt_summary_taxonomy_flow_rank:-phylum}" \
  hgt_taxonomy_flow_max_categories="${hgt_summary_taxonomy_flow_max_categories:-12}" \
  hgt_tree_plot_width="${hgt_summary_tree_plot_width:-24}" \
  hgt_promoter_bp="${hgt_summary_promoter_bp:-2000}" \
  hgt_fimo_qvalue="${hgt_summary_fimo_qvalue:-0.05}" \
    bash "${gg_core_dir}/gg_hgt_core.sh"
}

run_csubst_site_convergence_summary_for_source() {
  if [[ ${run_csubst_site_convergence_summary} -ne 1 ]]; then
    echo "Skipping convergent-site summary because run_csubst_site_convergence_summary=0."
    return 0
  fi
  echo "Running convergent-site summary for gene_family_source=${gene_family_source}."
  dir_orthogroup="${dir_gene_family}" \
  dir_orthofinder="${csubst_site_orthofinder_dir}" \
  dir_out="${csubst_site_output_dir}" \
  file_trait="${csubst_site_trait_file}" \
  arity_range="${csubst_site_arity_range}" \
  trait="${csubst_site_trait}" \
  skip_lower_order="${csubst_site_skip_lower_order}" \
  min_fg_stem_ratio="${csubst_site_min_fg_stem_ratio}" \
  min_OCNany2spe="${csubst_site_min_ocn_any2spe}" \
  min_omegaCany2spe="${csubst_site_min_omega_c_any2spe}" \
  min_OCNCoD="${csubst_site_min_ocn_cod}" \
  max_per_K="${csubst_site_max_candidates_per_arity}" \
  csubst_nonsyn_recode="${csubst_site_nonsyn_recode}" \
    bash "${gg_core_dir}/gg_convergent_sites_core.sh"
}

echo "gene_family_source=${gene_family_source}"
echo "dir_gene_family=${dir_gene_family}"
echo "summary_output_dir=${summary_output_dir}"

run_family_completion_summary_for_source
run_presence_absence_summary_for_source
run_gene_family_database_for_source
run_csubst_scan_aa_change_summary_for_source
run_hgt_summary_for_source
run_csubst_site_convergence_summary_for_source

echo "$(date): Exiting Singularity environment"
