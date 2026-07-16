# shellcheck shell=bash
# Sourced by gg_genome_evolution_core.sh.

task="OrthoFinder"
if [[ ! -s "${file_orthofinder_done_marker}" && ${run_orthofinder} -eq 1 ]]; then
  dir_sp_protein_orthofinder="${dir_sp_protein}_orthofinder"
  gg_step_start "${task}"
  prepare_species_protein_tmp
  prepare_species_protein_orthofinder_dir "${dir_sp_protein}" "${dir_sp_protein_orthofinder}"
  ensure_dir "${dir_orthofinder}"
  ensure_dir "${dir_orthofinder_hog2og}"

  orthofinder_algorithm_threads=$((${GG_TASK_CPUS} / 4))
  if [[ ${orthofinder_algorithm_threads} -eq 0 ]]; then
    orthofinder_algorithm_threads=1
  fi
  param_species_tree=()
  species_tree=""
  if [[ -s "${dir_species_tree_summary}/dated_species_tree.nwk" ]]; then
    species_tree="${dir_species_tree_summary}/dated_species_tree.nwk"
  elif [[ -s "${dir_species_tree_summary}/undated_species_tree.nwk" ]]; then
    species_tree="${dir_species_tree_summary}/undated_species_tree.nwk"
  fi
  if [[ -n "${species_tree}" && -s "${species_tree}" ]]; then
    echo "OrthoFinder will use the species tree: ${species_tree}"
    param_species_tree=(-s "${species_tree}")
  elif [[ ${species_tree_requested_for_orthofinder} -eq 1 ]]; then
    echo "Refusing to run OrthoFinder without a species tree."
    echo "Species-tree generation was requested, but no summary tree is available."
    echo "Expected one of:"
    echo "  ${dir_species_tree_summary}/dated_species_tree.nwk"
    echo "  ${dir_species_tree_summary}/undated_species_tree.nwk"
    echo "Please check the species-tree stage logs or disable species-tree generation flags intentionally before running OrthoFinder without a species tree."
    exit 1
  else
    echo "No species tree summary was found. Species-tree generation flags are disabled, so OrthoFinder will run without species tree constraints."
  fi
  echo "OrthoFinder will use ${GG_TASK_CPUS} threads for diamond search."
  echo "OrthoFinder will use ${orthofinder_algorithm_threads} threads for the OrthoFinder algorithm."

  protein_files=()
  mapfile -t protein_files < <(find "${dir_sp_protein_orthofinder}" -maxdepth 1 -type f ! -name '.*' \( -name "*.fa" -o -name "*.fa.gz" -o -name "*.fas" -o -name "*.fas.gz" -o -name "*.fasta" -o -name "*.fasta.gz" -o -name "*.fna" -o -name "*.fna.gz" \) | sort)
  num_sp=${#protein_files[@]}
  if [[ ${num_sp} -eq 0 ]]; then
    echo "No protein FASTA files were found in: ${dir_sp_protein_orthofinder}. Exiting."
    exit 1
  fi
  if [[ ${#param_species_tree[@]} -gt 0 ]]; then
    species_ids=()
    for protein_file in "${protein_files[@]}"; do
      species_base=$(basename "${protein_file}")
      species_base=${species_base%.gz}
      species_base=${species_base%.*}
      species_ids+=("${species_base}")
    done
    mapfile -t missing_species < <(
      python - "${species_tree}" "${species_ids[@]}" << 'PY'
import re
import sys

tree_file = sys.argv[1]
species = sys.argv[2:]
with open(tree_file) as f:
    tree = f.read().strip()

leaves = set(re.findall(r'(?<=[(,])([^:(),]+?)(?=[:),])', tree))
missing = [sp for sp in species if sp not in leaves]
for sp in missing:
    print(sp)
PY
    )
    if [[ ${#missing_species[@]} -gt 0 ]]; then
      echo "Species tree is missing ${#missing_species[@]} species: ${missing_species[*]}"
      echo "Refusing to run OrthoFinder without species tree constraints because a species tree was found but does not match the current OrthoFinder species set."
      echo "Please regenerate workspace/output/species_tree for the current inputs or update the species inputs to match the existing tree."
      exit 1
    fi
  fi
  if [[ ${num_sp} -gt ${max_orthofinder_core_species} ]]; then
    echo "The number of species (${num_sp}) is greater than the maximum number of core species (${max_orthofinder_core_species}) for OrthoFinder."
    echo "OrthoFinder will be run for 2 rounds (using --assign): For details, see https://github.com/davidemms/OrthoFinder"
    echo "OrthoFinder core-species filters: ${orthofinder_core_filters}"
    echo "OrthoFinder core-species rank keys: ${orthofinder_core_rank}"
    echo "OrthoFinder core-species sampling method: ${orthofinder_core_method}"

    species_cds_core=()
    select_orthofinder_core_args=(
      --protein-dir "${dir_sp_protein_orthofinder}"
      --busco-short-dir "${dir_species_busco_short}"
      --max-core-species "${max_orthofinder_core_species}"
      --filters "${orthofinder_core_filters}"
      --rank "${orthofinder_core_rank}"
      --method "${orthofinder_core_method}"
      --candidates-table "${file_orthofinder_core_candidates}"
      --selected-table "${file_orthofinder_core_selected}"
      --selected-list "${file_orthofinder_core_selected_list}"
      --core-tree "${file_orthofinder_core_species_tree}"
    )
    if [[ ${#param_species_tree[@]} -gt 0 ]]; then
      select_orthofinder_core_args+=(--species-tree "${species_tree}")
    fi
    if ! python "${gg_support_dir}/select_orthofinder_core_species.py" "${select_orthofinder_core_args[@]}"; then
      echo "Failed to select OrthoFinder core species. Exiting."
      exit 1
    fi
    mapfile -t species_cds_core < "${file_orthofinder_core_selected_list}"
    printf 'Core CDS files: %s\n' "${species_cds_core[*]}"
    if [[ ${#species_cds_core[@]} -eq 0 ]]; then
      echo "No OrthoFinder core species were selected. Exiting."
      exit 1
    fi
    if [[ -e "${dir_sp_protein}_core" ]]; then
      rm -rf -- "${dir_sp_protein}_core"
    fi
    mkdir -p "${dir_sp_protein}_core"
    for sp_cds_core in "${species_cds_core[@]}"; do
      cp_out "${dir_sp_protein_orthofinder}"/"${sp_cds_core}" "${dir_sp_protein}"_core
    done

    species_cds_additional=()
    mapfile -t species_cds_additional < <(comm -23 <(printf "%s\n" "${protein_files[@]##*/}" | sort) <(printf "%s\n" "${species_cds_core[@]}" | sort))
    printf 'Additional CDS files: %s\n' "${species_cds_additional[*]}"
    if [[ -e "${dir_sp_protein}_additional" ]]; then
      rm -rf -- "${dir_sp_protein}_additional"
    fi
    mkdir -p "${dir_sp_protein}_additional"
    for sp_cds_additional in "${species_cds_additional[@]}"; do
      cp_out "${dir_sp_protein_orthofinder}"/"${sp_cds_additional}" "${dir_sp_protein}"_additional
    done

    if [[ -e "${dir_orthofinder}"/core ]]; then
      rm -rf -- "${dir_orthofinder}/core"
    fi

    orthofinder_core_species_tree_args=()
    if [[ ${#param_species_tree[@]} -gt 0 ]]; then
      if [[ ! -s "${file_orthofinder_core_species_tree}" ]]; then
        echo "Core-species tree was not created: ${file_orthofinder_core_species_tree}"
        exit 1
      fi
      orthofinder_core_species_tree_args=(-s "${file_orthofinder_core_species_tree}")
    fi

    if orthofinder \
      -t "${GG_TASK_CPUS}" \
      -a "${orthofinder_algorithm_threads}" \
      -M "msa" \
      -S "diamond" \
      -f "${dir_sp_protein}_core" \
      -n "core" \
      -o "${dir_orthofinder}/core" \
      "${orthofinder_core_species_tree_args[@]}"; then
      orthofinder_core_exit_code=0
    else
      orthofinder_core_exit_code=$?
    fi

    if [[ ${orthofinder_core_exit_code} -ne 0 ]]; then
      echo "OrthoFinder failed in the core-species run. Exiting."
      exit 1
    fi

    if orthofinder \
      -t "${GG_TASK_CPUS}" \
      -a "${orthofinder_algorithm_threads}" \
      -M "msa" \
      -S "diamond" \
      -n "all" \
      --core "${dir_orthofinder}/core/Results_core" \
      --assign "${dir_sp_protein}_additional" \
      "${param_species_tree[@]}"; then
      orthofinder_main_exit_code=0
    else
      orthofinder_main_exit_code=$?
    fi

    if [[ ${orthofinder_main_exit_code} -ne 0 ]]; then
      echo "OrthoFinder failed in the all-species run. Exiting."
      exit 1
    fi

    shopt -s nullglob
    orthofinder_all_outputs=("${dir_orthofinder}"/core/Results_all/*)
    orthofinder_core_outputs=("${dir_orthofinder}"/core/Results_core/*)
    shopt -u nullglob
    if [[ ${#orthofinder_all_outputs[@]} -eq 0 || ${#orthofinder_core_outputs[@]} -eq 0 ]]; then
      echo "OrthoFinder core/main output files were expected but not found after completion."
      exit 1
    fi
    mv_out "${orthofinder_all_outputs[@]}" "${dir_orthofinder}"
    mv_out "${orthofinder_core_outputs[@]}" "${dir_orthofinder}/core"
    shopt -s nullglob
    orthofinder_result_dirs=("${dir_orthofinder}/core/Results_"*)
    shopt -u nullglob
    if [[ ${#orthofinder_result_dirs[@]} -gt 0 ]]; then
      rm -rf -- "${orthofinder_result_dirs[@]}"
    fi
    orthofinder_output_directory_cleanup "${dir_orthofinder}/core" "${GG_TASK_CPUS}"
  else
    echo "The number of species (${num_sp}) is less than or equal to the maximum number of core species (${max_orthofinder_core_species}) for OrthoFinder."
    echo "OrthoFinder will be run for 1 round."

    if orthofinder \
      -t "${GG_TASK_CPUS}" \
      -a "${orthofinder_algorithm_threads}" \
      -M "msa" \
      -S "diamond" \
      -f "${dir_sp_protein_orthofinder}" \
      -n "main" \
      -o "${dir_orthofinder}/main" \
      "${param_species_tree[@]}"; then
      orthofinder_main_exit_code=0
    else
      orthofinder_main_exit_code=$?
    fi

    if [[ ${orthofinder_main_exit_code} -ne 0 ]]; then
      echo "OrthoFinder failed in the all-species run. Exiting."
      exit 1
    fi

    shopt -s nullglob
    orthofinder_main_outputs=("${dir_orthofinder}"/main/Results_main/*)
    shopt -u nullglob
    if [[ ${#orthofinder_main_outputs[@]} -eq 0 ]]; then
      echo "OrthoFinder main output files were expected but not found after completion."
      exit 1
    fi
    mv_out "${orthofinder_main_outputs[@]}" "${dir_orthofinder}"
    rm -rf -- "${dir_orthofinder}/main"
  fi

  echo "OrthoFinder finished successfully."
  orthofinder_output_directory_cleanup "${dir_orthofinder}" "${GG_TASK_CPUS}"

  orthofinder_version=$(detect_orthofinder_version)
  if [[ -n "${orthofinder_version}" ]]; then
    echo "Detected OrthoFinder version: ${orthofinder_version}"
  else
    echo "Could not detect OrthoFinder version from orthofinder -v."
  fi

  hog_table="${dir_orthofinder}/Phylogenetic_Hierarchical_Orthogroups/N0.tsv"
  if orthofinder_supports_root_hog_equivalent "${orthofinder_version}"; then
    if ! copy_root_hog_equivalent_from_orthogroups "${dir_orthofinder}" "${dir_orthofinder_hog2og}" "${orthofinder_version}"; then
      if [[ -s "${hog_table}" ]]; then
        echo "Using detected N0.tsv because OrthoFinder v3.1+ Orthogroups.tsv files were incomplete: ${hog_table}"
        python "${gg_support_dir}/orthogroup_table_formatter.py" \
          --file_orthogroup_table "${hog_table}" \
          --dir_out "${dir_orthofinder_hog2og}" \
          --mode "hog2og"
      else
        echo "Neither OrthoFinder v3.1+ Orthogroups.tsv nor root-level HOG table N0.tsv was available. Exiting."
        exit 1
      fi
    fi
  elif [[ -s "${hog_table}" ]]; then
    python "${gg_support_dir}/orthogroup_table_formatter.py" \
      --file_orthogroup_table "${hog_table}" \
      --dir_out "${dir_orthofinder_hog2og}" \
      --mode "hog2og"
  else
    hog_candidates=()
    mapfile -t hog_candidates < <(find "${dir_orthofinder}/Phylogenetic_Hierarchical_Orthogroups" -maxdepth 1 -type f -name 'N*.tsv' | sort -V)
    echo "Required root-level HOG table was not found: ${hog_table}"
    echo "For OrthoFinder versions earlier than v3.1, genegalleon requires Phylogenetic_Hierarchical_Orthogroups/N0.tsv for whole-dataset HOG selection."
    if [[ ${#hog_candidates[@]} -gt 0 ]]; then
      echo "Other HOG tables were found, but they correspond to individual species-tree clades and will not be selected automatically:"
      printf '  %s\n' "${hog_candidates[@]}"
      printf '%s\n' 'Use OrthoFinder species-tree node labels to choose a clade-specific N*.tsv manually, upgrade to OrthoFinder v3.1-compatible root-HOG handling, or set orthogroup_table="OG" to use standard OrthoFinder orthogroups.'
    else
      echo "No N*.tsv HOG tables were found in ${dir_orthofinder}/Phylogenetic_Hierarchical_Orthogroups."
    fi
    exit 1
  fi
else
  gg_step_skip "${task}"
fi

task="OMArk analysis of species-wise protein input files"
run_shared_species_omark_stage

task="Summarizing OMArk species quality results"
run_shared_omark_summary_stage

task="Selecting orthogroups based on gene and species numbers"
if [[ ! -s "${file_orthogroup_selection}" && ${run_og_selection} -eq 1 ]]; then
  gg_step_start "${task}"
  prepare_species_protein_tmp
  ensure_dir "${dir_orthofinder_filtered}"
  if [[ "${orthogroup_annotation_method}" == "blastp" ]]; then
    if ! uniprot_db_prefix=$(ensure_uniprot_sprot_blast_db "${gg_workspace_dir}"); then
      echo "Failed to prepare UniProt Swiss-Prot BLASTP DB. Exiting."
      exit 1
    fi
    if ! validate_uniprot_sprot_db_prefix "${uniprot_db_prefix}" "blastp"; then
      echo "Invalid UniProt Swiss-Prot BLASTP DB prefix. Exiting."
      exit 1
    fi
  else
    if ! uniprot_db_prefix=$(ensure_uniprot_sprot_mmseqs_db "${gg_workspace_dir}"); then
      echo "Failed to prepare UniProt Swiss-Prot MMseqs2 DB. Exiting."
      exit 1
    fi
    if ! validate_uniprot_sprot_db_prefix "${uniprot_db_prefix}" "mmseqs2"; then
      echo "Invalid UniProt Swiss-Prot MMseqs2 DB prefix. Exiting."
      exit 1
    fi
  fi

  if [[ ${orthogroup_table} == "OG" ]]; then
    dir_orthogroup_selection_input="${dir_orthofinder_og}"
  elif [[ ${orthogroup_table} == "HOG" ]]; then
    dir_orthogroup_selection_input="${dir_orthofinder_hog2og}"
  else
    echo "Unsupported orthogroup_table: ${orthogroup_table}. Allowed values are OG or HOG."
    exit 1
  fi

  if [[ ! -s "${dir_orthogroup_selection_input}/Orthogroups.tsv" || ! -s "${dir_orthogroup_selection_input}/Orthogroups.GeneCount.tsv" ]]; then
    echo "Orthogroup source files were not found in ${dir_orthogroup_selection_input}. Exiting."
    exit 1
  fi
  cp_out "${dir_orthogroup_selection_input}/Orthogroups.tsv" "${dir_orthofinder_filtered}/Orthogroups.tsv"
  cp_out "${dir_orthogroup_selection_input}/Orthogroups.GeneCount.tsv" "${dir_orthofinder_filtered}/Orthogroups.GeneCount.tsv"

  if python "${gg_support_dir}/orthogroup_selection.py" \
    --dir_orthofinder_og "${dir_orthofinder_filtered}" \
    --dir_species_protein "${dir_sp_protein}" \
    --min_gene_num "${min_num_gene}" \
    --max_gene_num "${max_num_gene}" \
    --min_species_num "${min_num_species}" \
    --min_percent_species_coverage "${min_percent_species_coverage}" \
    --remove_unannotated 'yes' \
    --gene_size_quantiles '0.05,0.25,0.5,0.75,0.95' \
    --annotation_search_method "${orthogroup_annotation_method}" \
    --path_search_db "${uniprot_db_prefix}" \
    --evalue '1e-2' \
    --mmseqs_split_memory_limit "$(gg_memory_fraction_gb "${GG_MEM_TOOL_GB}" 3 4)G" \
    --ncpu "${GG_TASK_CPUS}"; then
    exit_code=0
  else
    exit_code=$?
  fi

  if [[ ${exit_code} -eq 0 ]]; then
    echo "Orthogroup selection finished successfully."
  else
    echo "Orthogroup selection failed. Exiting."
    exit 1
  fi
else
  gg_step_skip "${task}"
fi

task="Orthogroup method comparison"
disable_if_no_input_file "run_orthogroup_method_comparison" "${file_orthofinder_done_marker}"
if [[ ! -s "${file_orthogroup_method_comparison}" && ${run_orthogroup_method_comparison} -eq 1 ]]; then
  gg_step_start "${task}"

  if python "${gg_support_dir}/orthogroup_method_comparison.py" \
    --orthofinder_og_genecount "${dir_orthofinder_og}/Orthogroups.GeneCount.tsv" \
    --orthofinder_hog_genecount "${dir_orthofinder_hog2og}/Orthogroups.GeneCount.tsv"; then
    exit_code=0
  else
    exit_code=$?
  fi
  if [[ ${exit_code} -eq 0 ]]; then
    echo "Orthogroup method comparison finished successfully."
    mv_out orthogroup_histogram.pdf "${file_orthogroup_method_comparison}"
  else
    echo "Orthogroup method comparison failed. Exiting."
    exit 1
  fi
else
  gg_step_skip "${task}"
fi

task="Single-copy ortholog decay plot"
orthogroup_decay_genecount=""
if [[ "${orthogroup_table}" == "OG" ]]; then
  orthogroup_decay_genecount="${dir_orthofinder_og}/Orthogroups.GeneCount.tsv"
elif [[ "${orthogroup_table}" == "HOG" ]]; then
  orthogroup_decay_genecount="${dir_orthofinder_hog2og}/Orthogroups.GeneCount.tsv"
elif [[ ${run_single_copy_ortholog_decay_plot} -eq 1 ]]; then
  echo "Unsupported orthogroup_table for single-copy ortholog decay plot: ${orthogroup_table}. Allowed values are OG or HOG."
  exit 1
fi
disable_if_no_input_file "run_single_copy_ortholog_decay_plot" "${orthogroup_decay_genecount}"
if [[ (! -s "${file_single_copy_ortholog_decay_plot}" || ! -s "${file_single_copy_ortholog_decay_summary}") && ${run_single_copy_ortholog_decay_plot} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_orthogroup_decay}"

  if python "${gg_support_dir}/single_copy_ortholog_decay_plot.py" \
    --orthogroup-genecount "${orthogroup_decay_genecount}" \
    --outdir "${dir_orthogroup_decay}" \
    --replicates "${orthogroup_decay_replicates}" \
    --species-counts "${orthogroup_decay_species_counts}" \
    --seed "${orthogroup_decay_seed}" \
    --formats "pdf,svg"; then
    exit_code=0
  else
    exit_code=$?
  fi
  if [[ ${exit_code} -eq 0 ]]; then
    echo "Single-copy ortholog decay plot finished successfully."
  else
    echo "Single-copy ortholog decay plot failed. Exiting."
    exit 1
  fi
else
  gg_step_skip "${task}"
fi

cleanup_species_protein_tmp
trap - EXIT

# Genome evolution
if [[ ${run_busco_dupaware_notung_root_dna} -eq 1 || ${run_busco_dupaware_notung_root_pep} -eq 1 || ${run_busco_dupaware_root_dna} -eq 1 || ${run_busco_dupaware_root_pep} -eq 1 || ${run_busco_dupaware_grampa_dna} -eq 1 || ${run_busco_dupaware_grampa_pep} -eq 1 ]]; then
  if [[ ! -s "${notung_jar}" ]]; then
    echo "Notung jar was not found: ${notung_jar}"
    echo "Disabling NOTUNG-dependent duplicate-aware BUSCO tasks: run_busco_dupaware_notung_root_dna, run_busco_dupaware_notung_root_pep, run_busco_dupaware_root_dna, run_busco_dupaware_root_pep, run_busco_dupaware_grampa_dna, run_busco_dupaware_grampa_pep"
    run_busco_dupaware_notung_root_dna=0
    run_busco_dupaware_notung_root_pep=0
    run_busco_dupaware_root_dna=0
    run_busco_dupaware_root_pep=0
    run_busco_dupaware_grampa_dna=0
    run_busco_dupaware_grampa_pep=0
  fi
fi

if [[ ${run_orthogroup_grampa} -eq 1 ]]; then
  if [[ ! -s "${file_orthogroup_genecount_selected}" ]]; then
    echo "Disabling run_orthogroup_grampa because required file is missing: ${file_orthogroup_genecount_selected}"
    run_orthogroup_grampa=0
  elif [[ ! -d "${dir_og_rooted_tree}" ]]; then
    echo "Disabling run_orthogroup_grampa because required directory is missing: ${dir_og_rooted_tree}"
    run_orthogroup_grampa=0
  fi
fi

if [[ -z "${grampa_h1}" ]]; then
  if [[ ${run_busco_dupaware_grampa_dna} -eq 1 || ${run_busco_dupaware_grampa_pep} -eq 1 || ${run_orthogroup_grampa} -eq 1 ]]; then
    echo "Disabling GRAMPA tasks because grampa_h1 is empty. Set grampa_h1 in gg_genome_evolution_entrypoint.sh to enable them."
  fi
  run_busco_dupaware_grampa_dna=0
  run_busco_dupaware_grampa_pep=0
  run_orthogroup_grampa=0
fi

if [[ -z "${target_branch_go}" && ${run_go_enrichment} -eq 1 ]]; then
  echo "Disabling run_go_enrichment because target_branch_go is empty. Set target_branch_go in gg_genome_evolution_entrypoint.sh to enable it."
  run_go_enrichment=0
fi

ensure_dir "${dir_tmp}"
cd "${dir_tmp}" || exit 1
