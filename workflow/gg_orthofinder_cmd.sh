#!/usr/bin/env bash

### Start: Modify this block to tailor your analysis ###

run_cds_translation=1 # Translating species_cds to species_protein. Necessary for OrthoFinder
run_orthofinder=1 # OrthoFinder run
run_sonicparanoid=1 # SonicParanoid run
run_og_selection=1 # Selecting orthogroups for downstream analyses including gg_geneFamilyPhylogeny
run_orthogroup_method_comparison=1 # Plotting comparison of orthogroup selection methods

# Orthogroup table for downstream analyses
orthogroup_table="HOG" # "OG" for OrthoFinder's classical orthogroups, "HOG" for OrthoFinder's hierarchical orthogroups, "SP" for SonicParanoid's orthogroups

# Codon table
genetic_code=1

# Parameters for orthogroup selection
min_num_gene=4
min_num_species=2
max_orthofinder_core_species=50
min_percent_species_coverage=50
max_num_gene=1000

### End: Modify this block to tailor your analysis ###

### Modify below if you need to add a new analysis or need to fix some bugs ###

dir_pg="/workspace"
dir_myscript="/script/script"
source "${dir_myscript}/gg_util.sh" # Load utility functions
gg_source_home_bashrc
gg_prepare_cmd_runtime "${dir_pg}" "base" 1 1

enable_all_run_flags_for_debug_mode

dir_sp_cds="${dir_pg_input}/species_cds"
dir_sp_protein="${dir_pg_output}/species_protein"
dir_orthofinder="${dir_pg_output}/orthofinder"
dir_sonicparanoid="${dir_pg_output}/sonicparanoid"
dir_orthogroup_method_comparison="${dir_pg_output}/orthogroup_method_comparison"
dir_tmp="${dir_pg_output}/tmp"

dir_orthofinder_og="${dir_orthofinder}/Orthogroups"
dir_orthofinder_og_classical="${dir_orthofinder}/Orthogroups_classical"
dir_orthofinder_hog2og="${dir_orthofinder}/hog2og"
dir_sonicparanoid_sp2og="${dir_sonicparanoid}/sp2og"

file_orthofinder="${dir_orthofinder_hog2og}/README.txt"
file_sonicparanoid="${dir_sonicparanoid_sp2og}/README.txt"
file_orthogroup_selection="${dir_orthofinder_og}/Orthogroups.selected.tsv"
file_orthogroup_method_comparison="${dir_orthogroup_method_comparison}/orthogroup_method_comparison.pdf"

ensure_dir "${dir_tmp}"
cd "${dir_tmp}"
check_species_cds ${dir_pg}
check_if_species_files_unique ${dir_sp_cds}

task="Translate CDS sequences"
if [[ ${run_cds_translation} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_sp_protein}"
  cds_files=()
  mapfile -t cds_files < <(gg_find_fasta_files "${dir_sp_cds}" 1)
  echo "${#cds_files[@]} fasta files were detected."
  for cds_path in "${cds_files[@]}"; do
    cds=$(basename "${cds_path}")
    flag_seqkit=1
    if zgrep -e "^>" "${cds_path}" | grep -q -e "[[:blank:]]"; then
      echo "Space (\" \") is detected. Please remove all annotation info after spaces in sequence names. Exiting: ${cds}"
      exit 1
    fi
    sp_ub=$(gg_species_name_from_path "${cds}")
    translated_file="${sp_ub}.fasta"
    if [[ -s ${dir_sp_protein}/${translated_file} ]]; then
      echo "Translated file detected: ${cds}"
      flag_seqkit=0
    fi
    if [[ -s ${dir_sp_protein}/${translated_file} && ${dir_sp_protein}/${translated_file} -ot "${cds_path}" ]]; then
      echo "Translated file older than CDS file: ${cds}"
      flag_seqkit=1
    fi
    if [[ ${flag_seqkit} -eq 1 ]]; then
        echo "Translation started: ${cds}"
        touch "${cds_path}"

        seqkit seq --remove-gaps --threads ${NSLOTS} "${cds_path}" \
        | gg_prepare_cds_fasta_stream "${NSLOTS}" "${genetic_code}" \
        | seqkit translate --transl-table ${genetic_code} --threads ${NSLOTS} \
        | sed -e "s/^>${sp_ub}[-_\.]/>/" -e "s/^>/>${sp_ub}_/" \
        | sed -e '/^1 1$/d' -e 's/_frame=1[[:space:]]*//' \
        > ${dir_sp_protein}/${translated_file}
    else
        echo "Translation skipped: ${cds}"
    fi
  done
else
  gg_step_skip "${task}"
fi

task="OrthoFinder"
if [[ ! -s ${file_orthofinder} && ${run_orthofinder} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_orthofinder}"
  ensure_dir "${dir_orthofinder_hog2og}"

  nslots_of=$((${NSLOTS}/4))
  if [[ ${nslots_of} -eq 0 ]]; then
      nslots_of=1
  fi
  param_species_tree=''
  if [[ -s "${dir_pg_output}/species_tree/dated_species_tree.nwk" ]]; then
    file_species_tree_base="dated_species_tree"
  elif [[ -s "${dir_pg_output}/species_tree/undated_species_tree.nwk" ]]; then
    file_species_tree_base="undated_species_tree"
  fi
  species_tree="${dir_pg_output}/species_tree/${file_species_tree_base}.nwk"
  if [[ -s ${species_tree} ]]; then
    echo "OrthoFinder will use the species tree: ${species_tree}"
    param_species_tree="-s ${species_tree}"
  fi
  echo "OrthoFinder will use ${NSLOTS} threads for diamond search."
  echo "OrthoFinder will use ${nslots_of} threads for the OrthoFinder algorithm."

    protein_files=()
    mapfile -t protein_files < <(find "${dir_sp_protein}" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fa.gz" -o -name "*.fas" -o -name "*.fas.gz" -o -name "*.fasta" -o -name "*.fasta.gz" -o -name "*.fna" -o -name "*.fna.gz" \) | sort)
    num_sp=${#protein_files[@]}
    if [[ ${num_sp} -eq 0 ]]; then
      echo "No protein FASTA files were found in: ${dir_sp_protein}. Exiting."
      exit 1
    fi
  if [[ -n "${param_species_tree}" ]]; then
    species_ids=()
    for protein_file in "${protein_files[@]}"; do
      species_base=$(basename "${protein_file}")
      species_base=${species_base%.gz}
      species_base=${species_base%.*}
      species_ids+=( "${species_base}" )
    done
    mapfile -t missing_species < <(python - "${species_tree}" "${species_ids[@]}" << 'PY'
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
      echo "Running OrthoFinder without species tree constraints."
      param_species_tree=''
    fi
  fi
    if [[ ${num_sp} -gt ${max_orthofinder_core_species} ]]; then
    echo "The number of species (${num_sp}) is greater than the maximum number of core species (${max_orthofinder_core_species}) for OrthoFinder."
    echo "OrthoFinder will be run for 2 rounds (using --assign): For details, see https://github.com/davidemms/OrthoFinder"

      species_cds_core=()
      mapfile -t species_cds_core < <(printf "%s\n" "${protein_files[@]##*/}" | shuf -n ${max_orthofinder_core_species})
      echo "Core CDS files: ${species_cds_core[@]}"
    if [[ -e "${dir_sp_protein}_core" ]]; then
      rm -r "${dir_sp_protein}_core"
    fi
    mkdir -p "${dir_sp_protein}_core"
    for sp_cds_core in "${species_cds_core[@]}"; do
      cp_out ${dir_sp_protein}/${sp_cds_core} ${dir_sp_protein}_core
    done

      species_cds_additional=()
      mapfile -t species_cds_additional < <(printf "%s\n" "${protein_files[@]##*/}" | grep -v -F -x -f <(printf "%s\n" "${species_cds_core[@]}"))
      echo "Additional CDS files: ${species_cds_additional[@]}"
    if [[ -e "${dir_sp_protein}_additional" ]]; then
      rm -r "${dir_sp_protein}_additional"
    fi
    mkdir -p "${dir_sp_protein}_additional"
    for sp_cds_additional in "${species_cds_additional[@]}"; do
      cp_out ${dir_sp_protein}/${sp_cds_additional} ${dir_sp_protein}_additional
    done

    if [[ -e ${dir_orthofinder}/core ]]; then
      rm -r ${dir_orthofinder}/core
    fi

    if [[ -e "${dir_orthofinder}/species_tree_core.nwk" ]]; then
      rm "${dir_orthofinder}/species_tree_core.nwk"
    fi
      core_species_names=()
      mapfile -t core_species_names < <(
        find "${dir_sp_protein}_core" -maxdepth 1 -type f ! -name '.*' | sort \
        | awk '{name=$0; sub(/^.*\//, "", name); sub(/\.[^.]*$/, "", name); print name}'
      )
      core_species_regex=$(printf '%s|' "${core_species_names[@]}")
      core_species_regex=${core_species_regex%|}
      if [[ -z "${core_species_regex}" ]]; then
        echo "Failed to build core species regex. No core protein files were detected. Exiting."
        exit 1
      fi
      nwkit prune --invert_match yes --pattern "${core_species_regex}" --infile ${species_tree} --outfile "${dir_orthofinder}/species_tree_core.nwk"

    orthofinder \
    -t ${NSLOTS} \
    -a ${nslots_of} \
    -M "msa" \
    -S "diamond" \
    -f "${dir_sp_protein}_core" \
    -n "core" \
    -o "${dir_orthofinder}/core" \
    -s "${dir_orthofinder}/species_tree_core.nwk"; orthofinder_core_exit_code=$?

    if [[ ${orthofinder_core_exit_code} -ne 0 ]]; then
      echo "OrthoFinder failed in the core-species run. Exiting."
      exit 1
    fi

    orthofinder \
    -t ${NSLOTS} \
    -a ${nslots_of} \
    -M "msa" \
    -S "diamond" \
    -n "all" \
    --core "${dir_orthofinder}/core/Results_core" \
    --assign "${dir_sp_protein}_additional" \
    ${param_species_tree}; orthofinder_main_exit_code=$?

    if [[ ${orthofinder_main_exit_code} -ne 0 ]]; then
      echo "OrthoFinder failed in the all-species run. Exiting."
      exit 1
    fi

    mv_out ${dir_orthofinder}/core/Results_all/* ${dir_orthofinder}
    mv_out ${dir_orthofinder}/core/Results_core/* ${dir_orthofinder}/core
    rm -r ${dir_orthofinder}/core/Results_*
    orthofinder_output_directory_cleanup ${dir_orthofinder}/core ${NSLOTS}
  else
    echo "The number of species (${num_sp}) is less than or equal to the maximum number of core species (${max_orthofinder_core_species}) for OrthoFinder."
    echo "OrthoFinder will be run for 1 round."

    orthofinder \
    -t ${NSLOTS} \
    -a ${nslots_of} \
    -M "msa" \
    -S "diamond" \
    -f ${dir_sp_protein} \
    -n "main" \
    -o ${dir_orthofinder}/main \
    ${param_species_tree}; orthofinder_main_exit_code=$?

    if [[ ${orthofinder_main_exit_code} -ne 0 ]]; then
      echo "OrthoFinder failed in the all-species run. Exiting."
      exit 1
    fi

    mv_out ${dir_orthofinder}/main/Results_main/* ${dir_orthofinder}
    rm -r ${dir_orthofinder}/main
  fi

  echo "OrthoFinder finished successfully."
  orthofinder_output_directory_cleanup ${dir_orthofinder} ${NSLOTS}
  mv_out ${dir_orthofinder_og} ${dir_orthofinder_og_classical}

  hog_table="${dir_orthofinder}/Phylogenetic_Hierarchical_Orthogroups/N0.tsv"
  if [[ ! -s "${hog_table}" ]]; then
    hog_candidates=()
    mapfile -t hog_candidates < <(find "${dir_orthofinder}/Phylogenetic_Hierarchical_Orthogroups" -maxdepth 1 -type f -name 'N*.tsv' | sort -V)
    if [[ ${#hog_candidates[@]} -eq 0 ]]; then
      echo "No HOG tables were found in ${dir_orthofinder}/Phylogenetic_Hierarchical_Orthogroups. Exiting."
      exit 1
    fi
    hog_table="${hog_candidates[0]}"
    echo "N0.tsv was not found. Falling back to HOG table: ${hog_table}"
  fi

  python ${dir_myscript}/orthogroup_table_formatter.py \
  --file_orthogroup_table "${hog_table}" \
  --dir_out "${dir_orthofinder_hog2og}" \
  --mode "hog2og"
else
  gg_step_skip "${task}"
fi

task="SonicParanoid"
if [[ ! -s ${file_sonicparanoid} && ${run_sonicparanoid} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_sonicparanoid}"
  ensure_dir "${dir_sonicparanoid_sp2og}"

  sonicparanoid_runner=(sonicparanoid)
  if [[ "$(uname -m)" == "aarch64" && -x "/opt/conda/envs/base/bin/python" ]]; then
    # SonicParanoid's install-type detection can misclassify /opt/conda/bin/python
    # as "Python" on arm64, which triggers replacement with bundled x86 binaries.
    sonicparanoid_runner=(/opt/conda/envs/base/bin/python -m sonicparanoid.sonic_paranoid)
    echo "Using ARM-safe SonicParanoid runner: ${sonicparanoid_runner[*]}"
  fi

  "${sonicparanoid_runner[@]}" \
  --input-directory ${dir_sp_protein} \
  --output-directory "${dir_sonicparanoid}" \
  --project-id "gg_run" \
  --graph-only \
  --overwrite \
  --threads ${NSLOTS}
  sonicparanoid_exit_code=$?
  if [[ ${sonicparanoid_exit_code} -ne 0 ]]; then
    echo "SonicParanoid failed with exit code ${sonicparanoid_exit_code}. Exiting."
    exit 1
  fi

  mv_out "${dir_sonicparanoid}/runs/gg_run/ortholog_groups" "${dir_sonicparanoid}"
  rm -r "${dir_sonicparanoid}/alignments"
  rm -r "${dir_sonicparanoid}/orthologs_db"
  rm -r "${dir_sonicparanoid}/runs"

  python ${dir_myscript}/orthogroup_table_formatter.py \
  --mode "sp2og" \
  --file_orthogroup_table "${dir_sonicparanoid}/ortholog_groups/flat.ortholog_groups.tsv" \
  --dir_out "${dir_sonicparanoid}/sp2og"

  if [[ -s "${dir_sonicparanoid}/sp2og/README.txt" ]]; then
    echo "Output file detected for the task: ${task}"
  else
    echo "No output file detected for the task: ${task}"
    exit 1
  fi
else
  gg_step_skip "${task}"
fi

task="Selecting orthogroups based on gene and species numbers"
if [[ ! -s ${file_orthogroup_selection} && ${run_og_selection} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_orthofinder_og}"
  uniprot_db_prefix=$(ensure_uniprot_sprot_db "${dir_pg}")
  if [[ $? -ne 0 ]]; then
    echo "Failed to prepare UniProt Swiss-Prot DB. Exiting."
    exit 1
  fi

  if [[ ${orthogroup_table} == "OG" ]]; then
    dir_orthogroup_selection_input=${dir_orthofinder_og_classical}
  elif [[ ${orthogroup_table} == "HOG" ]]; then
    dir_orthogroup_selection_input=${dir_orthofinder_hog2og}
  elif [[ ${orthogroup_table} == "SP" ]]; then
    dir_orthogroup_selection_input="${dir_sonicparanoid_sp2og}"
  fi

  python ${dir_myscript}/orthogroup_selection.py \
  --dir_orthofinder_og ${dir_orthogroup_selection_input} \
  --dir_species_protein ${dir_sp_protein} \
  --min_gene_num ${min_num_gene} \
  --max_gene_num ${max_num_gene} \
  --min_species_num ${min_num_species} \
  --min_percent_species_coverage ${min_percent_species_coverage} \
  --remove_unannotated 'yes' \
  --gene_size_quantiles '0.05,0.25,0.5,0.75,0.95' \
  --path_diamond_db "${uniprot_db_prefix}" \
  --evalue '1e-2' \
  --ncpu ${NSLOTS}; exit_code=$?

  if [[ ${exit_code} -eq 0 ]]; then
    echo "Orthogroup selection finished successfully."
    echo "Copying files from ${dir_orthogroup_selection_input} to ${dir_orthofinder_og}"
    cp_out -r ${dir_orthogroup_selection_input} ${dir_orthofinder_og}
  else
    echo "Orthogroup selection failed. Exiting."
    exit 1
  fi
else
  gg_step_skip "${task}"
fi

task="Orthogroup method comparison"
disable_if_no_input_file "run_orthogroup_method_comparison" ${file_orthofinder}
disable_if_no_input_file "run_orthogroup_method_comparison" ${file_sonicparanoid}
if [[ ! -s ${file_orthogroup_method_comparison} && ${run_orthogroup_method_comparison} -eq 1 ]]; then
  gg_step_start "${task}"

  python ${dir_myscript}/orthogroup_method_comparison.py \
  --orthofinder_og_genecount ${dir_orthofinder_og_classical}/Orthogroups.GeneCount.tsv \
  --orthofinder_hog_genecount ${dir_orthofinder_hog2og}/Orthogroups.GeneCount.tsv \
  --sonicparanoid_genecount ${dir_sonicparanoid_sp2og}/Orthogroups.GeneCount.tsv
  exit_code=$?
  if [[ ${exit_code} -eq 0 ]]; then
    echo "Orthogroup method comparison finished successfully."
    mv_out orthogroup_histogram.pdf ${file_orthogroup_method_comparison}
  else
    echo "Orthogroup method comparison failed. Exiting."
    exit 1
  fi
else
  gg_step_skip "${task}"
fi

echo "$(date): Exiting Singularity environment"
