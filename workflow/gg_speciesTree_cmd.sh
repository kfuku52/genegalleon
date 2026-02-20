#!/usr/bin/env bash

#run a busco analysis for all species in species_cds and extract shared complete BUSCO genes among all species.
#then make cds fasta files including all species' cds sequenes for each BUSCO genes and make alignment
#then concatenate all fasta files for further downstream analysis -> task = trimal in gg_species_tree.

### Start: Modify this block to tailor your analysis ###

# Mode
# Values should be boolean (0 or 1).
# Only 1 mode should be activated.
mode_orthogroup=0 # Analyze OrthoFinder-based single-copy genes, not supported yet.
mode_busco=1 # Analyze BUSCO-based single-copy genes

# Workflow
# Values should be boolean (0 or 1).
run_busco=1
run_get_busco_summary=1
run_individual_get_fasta=1
run_individual_mafft=1
run_individual_trimal=1
# Species tree with concatenated sequences
run_concat_alignment=1
run_concat_iqtree_protein=1
run_concat_iqtree_dna=1
# Species tree with ASTRAL
run_individual_iqtree_pep=1
run_astral_pep=1
run_individual_iqtree_dna=1
run_astral_dna=1
run_plot_species_trees=1 # Plotting 4 species trees
# Divergence time estimation (IQ2MC workflow)
run_constrained_tree=1 # Introduce divergence time constraints for IQ2MC input
run_plot_constrained_tree=1 # Plot node-wise constrained ranges in constrained.nwk
run_mcmctree1=1 # IQ2MC step 2: Hessian estimation + mcmctree control generation
run_mcmctree2=1 # IQ2MC step 3: MCMC run by modified mcmctree
run_convert_tree_format=1
run_plot_mcmctreer=1 # Plotting dated species tree using the local R plotting script

# Other parameters
strictly_single_copy_only=0 # To remove BUSCO gene groups where some species have missing or duplicated BUSCO genes.
genetic_code=1 # Integer. See here https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes
busco_lineage="embryophyta_odb12" # See here for available datasets: https://busco-data.ezlab.org/v5/data/lineages/
bootstrap_params="-bb 1000 -bnni" # IQ-TREE's bootstrap parameters
nucleotide_model='GTR+R4' # IQ-TREE's model
protein_model='LG+R4' # IQ-TREE's model
undated_species_tree='astral_pep' # {iqtree_dna,iqtree_pep,astral_dna,astral_pep}
astral_min_tips=4 # Minimum number of tips per gene tree to include in ASTRAL.

# Constraints for divergence time estimation
timetree_constraint=1 # Whether to automatically retrieve divergence time constraints from timetree.org using `nwkit mcmctree --timetree`. If set to 1, mcmctree_divergence_time_constraints will not be used.
#outgroup_labels='Amborella_trichopoda'
outgroup_labels='Oryza_sativa' # Comma-separated list of outgroup species labels
# left_species,right_species,lower_bound,upper_bound
#mcmctree_divergence_time_constraints=(
#  "Arabidopsis_thaliana,Amborella_trichopoda,-,200"
#  "Arabidopsis_thaliana,Populus_trichocarpa,100,120"
#  "Arabidopsis_thaliana,Oryza_sativa,130,-"
#)
mcmctree_divergence_time_constraints=(
  "Arabidopsis_thaliana,Oryza_sativa,130,-"
)
mcmc_burnin=20000 # https://link.springer.com/protocol/10.1007/978-1-4939-9074-0_10
mcmc_sampfreq=100 # https://link.springer.com/protocol/10.1007/978-1-4939-9074-0_10
mcmc_nsample=20000 # https://link.springer.com/protocol/10.1007/978-1-4939-9074-0_10
mcmc_birth_death_sampling="1,1,0.5" # birth,death,sampling_fraction for --mcmc-bds
mcmc_clock_model="IND" # Clock model for IQ2MC: {EQUAL, IND, CORR}

### End: Modify this block to tailor your analysis ###

### ----------------------------------------------------------------------- ###

### Modify below if you need to add a new analysis or need to fix some bugs ###

dir_pg="/workspace"
dir_myscript="/script/script"
source "${dir_myscript}/gg_util.sh" # Load utility functions
gg_source_home_bashrc
gg_prepare_cmd_runtime "${dir_pg}" "base" 1 1

copy_busco_tables() {
  local busco_root_dir="$1"
  local lineage="$2"
  local file_full="$3"
  local file_short="$4"
  local run_dir="${busco_root_dir}/run_${lineage}"
  local full_src=""
  local short_src=""
  local -a short_candidates=()

  if [[ -s "${run_dir}/full_table.tsv" ]]; then
    full_src="${run_dir}/full_table.tsv"
  elif [[ -s "${run_dir}/full_table.tsv.gz" ]]; then
    gzip -cd "${run_dir}/full_table.tsv.gz" > "./tmp.busco.full_table.tsv"
    full_src="./tmp.busco.full_table.tsv"
  fi

  if [[ -s "${run_dir}/short_summary.txt" ]]; then
    short_src="${run_dir}/short_summary.txt"
  else
    mapfile -t short_candidates < <(find "${run_dir}" -maxdepth 1 -type f -name "short_summary*.txt" | sort)
    if [[ ${#short_candidates[@]} -gt 0 ]]; then
      short_src="${short_candidates[0]}"
    fi
  fi

  if [[ -z "${full_src}" || -z "${short_src}" ]]; then
    echo "BUSCO outputs were not found under ${run_dir}."
    echo "Expected full_table.tsv(.gz) and short_summary*.txt."
    rm -f "./tmp.busco.full_table.tsv"
    return 1
  fi

  cp_out "${full_src}" "${file_full}"
  cp_out "${short_src}" "${file_short}"
  rm -f "./tmp.busco.full_table.tsv"
}

dir_sp_cds="${dir_pg_input}/species_cds"
dir_og="${dir_pg_output}/orthogroup"
dir_og_rooted_tree="${dir_og}/rooted_tree"
dir_species_tree="${dir_pg_output}/species_tree"
dir_species_busco_full="${dir_pg_output}/species_cds_busco_full"
dir_species_busco_short="${dir_pg_output}/species_cds_busco_short"
dir_busco_summary_table="${dir_species_tree}/busco_summary_table"
dir_single_copy_fasta="${dir_species_tree}/single_copy_cds_fasta"
dir_single_copy_mafft="${dir_species_tree}/single_copy_mafft"
dir_single_copy_trimal="${dir_species_tree}/single_copy_trimal"
dir_single_copy_iqtree_dna="${dir_species_tree}/single_copy_iqtree_dna"
dir_single_copy_iqtree_pep="${dir_species_tree}/single_copy_iqtree_pep"
dir_astral_dna="${dir_species_tree}/single_copy_astral_dna"
dir_astral_pep="${dir_species_tree}/single_copy_astral_pep"
dir_concat_fasta="${dir_species_tree}/concatenated_alignment"
dir_concat_iqtree_dna="${dir_species_tree}/concatenated_iqtree_dna"
dir_concat_iqtree_pep="${dir_species_tree}/concatenated_iqtree_pep"
dir_constrained_tree="${dir_species_tree}/constrained_tree"
dir_mcmctree1="${dir_species_tree}/mcmctree_parameter_estimation"
dir_mcmctree2="${dir_species_tree}/mcmctree_main"
dir_tmp="${dir_species_tree}/tmp"

check_species_cds ${dir_pg}
check_if_species_files_unique ${dir_sp_cds}
memory_notung=${MEM_PER_SLOT}

file_busco_summary_table="${dir_busco_summary_table}/busco_summary.tsv"
file_astral_tree_dna_q1="${dir_astral_dna}/single_copy.astral.dna.q1.nwk" # Quartet supports for the main topology; The lengths of terminal branches are not computed.
file_astral_tree_dna="${dir_astral_dna}/single_copy.astral.dna.optimized.nwk" # ASTRAL topology with branch lengths optimized by IQ-TREE.
file_astral_log_dna="${dir_astral_dna}/single_copy.astral.dna.log"
file_astral_tree_pep_q1="${dir_astral_pep}/single_copy.astral.pep.q1.nwk" # Quartet supports for the main topology; The lengths of terminal branches are not computed.
file_astral_tree_pep="${dir_astral_pep}/single_copy.astral.pep.optimized.nwk" # ASTRAL topology with branch lengths optimized by IQ-TREE.
file_astral_log_pep="${dir_astral_pep}/single_copy.astral.pep.log"
file_concat_cds="${dir_concat_fasta}/concat.cds.trimal.fa.gz"
file_concat_pep="${dir_concat_fasta}/concat.pep.trimal.fa.gz"
file_concat_iqtree_dna="${dir_concat_iqtree_dna}/concat.dna.iqtree.nwk"
file_concat_iqtree_pep="${dir_concat_iqtree_pep}/concat.pep.iqtree.nwk"
file_concat_iqtree_dna_root="${dir_concat_iqtree_dna}/concat.dna.rooted.nwk"
file_concat_iqtree_pep_root="${dir_concat_iqtree_pep}/concat.pep.rooted.nwk"
file_constrained_tree="${dir_constrained_tree}/constrained.nwk"
file_plot_constrained_tree="${dir_constrained_tree}/constrained_tree_constraints.pdf"
file_iq2mc_prefix="${dir_mcmctree1}/iq2mc"
file_iq2mc_ctl="${file_iq2mc_prefix}.mcmctree.ctl"
file_iq2mc_hessian="${file_iq2mc_prefix}.mcmctree.hessian"
file_iq2mc_rooted_tree="${file_iq2mc_prefix}.rooted.nwk"
file_iq2mc_dummy_phy="${file_iq2mc_prefix}.dummy.phy"
file_mcmctree1="${file_iq2mc_ctl}" # Backward-compatible marker name
file_mcmctree2_raw="${dir_mcmctree2}/iq2mc.mcmctree.out"
file_mcmctree2="${dir_mcmctree2}/FigTree.tre"
file_mcmctree2_nwk="${dir_mcmctree2}/dated_species_tree.nwk"
file_plot_species_trees="${dir_species_tree}/species_trees.pdf"
file_plot_mcmctreer="${file_mcmctree2_nwk%.*}.pdf"
file_dated_species_tree="${dir_species_tree}/dated_species_tree.nwk"
file_undated_species_tree="${dir_species_tree}/undated_species_tree.nwk"

ensure_dir "${dir_tmp}"
cd ${dir_tmp}

enable_all_run_flags_for_debug_mode

root_tree_with_outgroup () {
  local infile=$1
  local outfile=$2
  local tree_description=$3
  local root_log="${dir_tmp}/tmp.nwkit.root.$$.log"
  local outgroup_label_list=()
  ensure_parent_dir "${outfile}"
  rm -f "${outfile}" "${root_log}"

  local missing_outgroup=0
  mapfile -t outgroup_label_list < <(printf '%s' "${outgroup_labels}" | tr ',' '\n' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' -e '/^$/d')
  for outgroup_label in "${outgroup_label_list[@]}"; do
    if ! grep -q -F "${outgroup_label}" "${infile}"; then
      missing_outgroup=1
      break
    fi
  done
  if [[ ${missing_outgroup} -eq 1 ]]; then
    echo "Error: Outgroup labels (${outgroup_labels}) are not present in ${tree_description}."
    return 1
  fi

  nwkit root \
  --method outgroup \
  --outgroup "${outgroup_labels}" \
  --infile "${infile}" \
  --outfile "${outfile}" \
  2> "${root_log}"
  local root_exit_code=$?

  if [[ ${root_exit_code} -eq 0 && -s "${outfile}" ]]; then
    rm -f "${root_log}"
    return 0
  fi

  echo "Error: Failed to root ${tree_description}."
  [[ -s "${root_log}" ]] && cat "${root_log}"
  rm -f "${root_log}"
  return 1
}

normalize_iq2mc_constraint_tree () {
  local infile=$1
  local tmpfile

  if [[ ! -s "${infile}" ]]; then
    echo "Error: Constraint tree file is missing or empty: ${infile}"
    return 1
  fi

  tmpfile=$(mktemp "${dir_tmp}/tmp.iq2mc.constraint.XXXXXX")
  # IQ-TREE's dated-tree parser rejects spaces after commas inside B(...) labels.
  tr -d '\r' < "${infile}" | sed -E 's/,[[:space:]]+/,/g; s/[[:space:]]+$//' > "${tmpfile}"

  if [[ ! -s "${tmpfile}" ]]; then
    echo "Error: Failed to normalize IQ2MC constraint tree: ${infile}"
    rm -f "${tmpfile}"
    return 1
  fi

  mv_out "${tmpfile}" "${infile}"
  return 0
}

iq2mc_option_supported () {
  local candidate=$1
  local resolved_candidate
  resolved_candidate=$(command -v "${candidate}" 2>/dev/null || true)
  if "${candidate}" -h 2>&1 | grep -q -- "--mcmc-bds"; then
    return 0
  fi
  if command -v strings >/dev/null 2>&1 && [[ -n "${resolved_candidate}" ]]; then
    if strings "${resolved_candidate}" 2>/dev/null | grep -q -- "--mcmc-bds"; then
      return 0
    fi
  fi
  return 1
}

resolve_iq2mc_binary () {
  local candidate
  for candidate in iqtree3 iqtree; do
    if ! command -v "${candidate}" >/dev/null 2>&1; then
      continue
    fi
    if iq2mc_option_supported "${candidate}"; then
      echo "${candidate}"
      return 0
    fi
  done
  return 1
}

count_newick_tips () {
  local tree_file=$1
  python - "${tree_file}" <<'PY'
import sys
from Bio import Phylo
try:
    tree = Phylo.read(sys.argv[1], "newick")
    print(len(tree.get_terminals()))
except Exception:
    print(0)
PY
}

build_astral_input () {
  local tree_dir=$1
  local merged_file=$2
  local min_tips=$3
  local kept skipped ntips tree_file
  local -a tree_files

  kept=0
  skipped=0
  : > "${merged_file}"

  mapfile -t tree_files < <(find "${tree_dir}" -maxdepth 1 -type f -name "*.nwk" | sort)
  if [[ ${#tree_files[@]} -eq 0 ]]; then
    echo "No gene-tree files were found for ASTRAL in: ${tree_dir}"
    rm -f "${merged_file}"
    return 1
  fi

  for tree_file in "${tree_files[@]}"; do
    ntips=$(count_newick_tips "${tree_file}")
    if [[ "${ntips}" -ge "${min_tips}" ]]; then
      cat "${tree_file}" >> "${merged_file}"
      echo "" >> "${merged_file}"
      kept=$((kept + 1))
    else
      skipped=$((skipped + 1))
      echo "Skipped for ASTRAL (tips=${ntips} < ${min_tips}): $(basename "${tree_file}")"
    fi
  done

  echo "ASTRAL input tree filtering summary: kept=${kept}, skipped=${skipped}, min_tips=${min_tips}"
  if [[ ${kept} -eq 0 ]]; then
    rm -f "${merged_file}"
    return 1
  fi
  return 0
}

sanitize_newick_topology_for_iqtree () {
  local infile=$1
  local outfile=$2
  python - "${infile}" "${outfile}" <<'PY'
import sys
from Bio import Phylo

try:
    tree = Phylo.read(sys.argv[1], "newick")
except Exception as exc:
    sys.stderr.write(f"Failed to read tree file for IQ-TREE topology sanitization: {exc}\n")
    sys.exit(1)

for clade in tree.find_clades():
    if not clade.is_terminal():
        clade.name = None
        clade.confidence = None

Phylo.write(tree, sys.argv[2], "newick")
PY
}

optimize_astral_tree_branch_lengths () {
  local astral_support_tree=$1
  local concat_alignment=$2
  local model=$3
  local optimized_outfile=$4
  local tag=$5
  local tmp_topology="tmp.astral.${tag}.topology.nwk"
  local iqtree_prefix="tmp.astral.${tag}.optimized"
  local tmp_concat_alignment="tmp.astral.${tag}.concat.fasta"
  local tmp_rooted_optimized_tree="tmp.astral.${tag}.optimized.rooted.nwk"
  local tmp_optimized_tree_with_support="tmp.astral.${tag}.optimized.with_support.nwk"

  if [[ ! -s "${astral_support_tree}" ]]; then
    echo "Warning: Cannot optimize ASTRAL tree branch lengths. Missing tree: ${astral_support_tree}"
    return 1
  fi
  if [[ ! -s "${concat_alignment}" ]]; then
    echo "Warning: Cannot optimize ASTRAL tree branch lengths. Missing alignment: ${concat_alignment}"
    return 1
  fi

  rm -f "${tmp_topology}" "${iqtree_prefix}".* "${tmp_concat_alignment}" "${tmp_rooted_optimized_tree}" "${tmp_optimized_tree_with_support}"
  if ! sanitize_newick_topology_for_iqtree "${astral_support_tree}" "${tmp_topology}"; then
    echo "Warning: Failed to sanitize ASTRAL topology for IQ-TREE branch-length optimization."
    rm -f "${tmp_topology}" "${iqtree_prefix}".* "${tmp_concat_alignment}" "${tmp_rooted_optimized_tree}" "${tmp_optimized_tree_with_support}"
    return 1
  fi
  seqkit seq --threads 1 "${concat_alignment}" --out-file "${tmp_concat_alignment}"
  if [[ ! -s "${tmp_concat_alignment}" ]]; then
    echo "Warning: Failed to prepare alignment for IQ-TREE branch-length optimization: ${concat_alignment}"
    rm -f "${tmp_topology}" "${iqtree_prefix}".* "${tmp_concat_alignment}" "${tmp_rooted_optimized_tree}" "${tmp_optimized_tree_with_support}"
    return 1
  fi

  iqtree \
  -s "${tmp_concat_alignment}" \
  -te "${tmp_topology}" \
  -m "${model}" \
  -n 0 \
  -T ${NSLOTS} \
  --prefix "${iqtree_prefix}" \
  --seed 12345 \
  --redo

  if [[ ! -s "${iqtree_prefix}.treefile" ]]; then
    echo "Warning: IQ-TREE optimization did not produce a tree for ${tag}."
    rm -f "${tmp_topology}" "${iqtree_prefix}".* "${tmp_concat_alignment}" "${tmp_rooted_optimized_tree}" "${tmp_optimized_tree_with_support}"
    return 1
  fi

  root_tree_with_outgroup \
  "${iqtree_prefix}.treefile" \
  "${tmp_rooted_optimized_tree}" \
  "optimized ASTRAL ${tag} tree"
  local root_exit_code=$?
  if [[ ${root_exit_code} -ne 0 || ! -s "${tmp_rooted_optimized_tree}" ]]; then
    echo "Warning: Failed to save optimized ASTRAL ${tag} tree."
    rm -f "${tmp_topology}" "${iqtree_prefix}".* "${tmp_concat_alignment}" "${tmp_rooted_optimized_tree}" "${tmp_optimized_tree_with_support}"
    return 1
  fi

  nwkit transfer \
  --infile "${tmp_rooted_optimized_tree}" \
  --infile2 "${astral_support_tree}" \
  --target intnode \
  --support yes \
  --name no \
  --length no \
  --outfile "${tmp_optimized_tree_with_support}"
  local transfer_exit_code=$?
  if [[ ${transfer_exit_code} -eq 0 && -s "${tmp_optimized_tree_with_support}" ]]; then
    cp_out "${tmp_optimized_tree_with_support}" "${optimized_outfile}"
  else
    echo "Warning: Failed to transfer ASTRAL support values to optimized ${tag} tree. Keeping branch-length-optimized tree without support transfer."
    cp_out "${tmp_rooted_optimized_tree}" "${optimized_outfile}"
  fi
  local copy_exit_code=$?
  rm -f "${tmp_topology}" "${iqtree_prefix}".* "${tmp_concat_alignment}" "${tmp_rooted_optimized_tree}" "${tmp_optimized_tree_with_support}"

  if [[ ${copy_exit_code} -ne 0 || ! -s "${optimized_outfile}" ]]; then
    echo "Warning: Failed to write optimized ASTRAL ${tag} tree."
    return 1
  fi
  return 0
}

iq2mc_binary=""
if [[ ${run_mcmctree1} -eq 1 || ${run_mcmctree2} -eq 1 ]]; then
  echo "IQ2MC dating will run. Checking tools and calibration species."
  if ! iq2mc_binary=$(resolve_iq2mc_binary); then
    echo "Error: IQ2MC-capable IQ-TREE was not found."
    echo "Please install an IQ-TREE build that supports --dating mcmctree (and --mcmc-bds)."
    exit 1
  fi
  if ! iq2mc_option_supported "${iq2mc_binary}"; then
    echo "Warning: Could not verify --mcmc-bds support in ${iq2mc_binary} by help/strings checks."
    echo "Warning: Proceeding with ${iq2mc_binary}; IQ2MC step 2 will fail if this binary is incompatible."
  fi
  if ! command -v mcmctree >/dev/null 2>&1; then
    echo "Error: mcmctree command was not found."
    echo "Please install the IQ2MC-compatible mcmctree binary."
    exit 1
  fi
  echo "Using IQ2MC binary: ${iq2mc_binary}"
  outgroup_label_list=()
  mapfile -t outgroup_label_list < <(printf '%s' "${outgroup_labels}" | tr ',' '\n' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' -e '/^$/d')
  for outgroup_label in "${outgroup_label_list[@]}"; do
    stop_if_species_not_found_in ${dir_sp_cds} ${outgroup_label}
  done
  if [[ ${timetree_constraint} -eq 0 ]]; then
    mcmctree_params=(${mcmctree_divergence_time_constraints//,/ })
    stop_if_species_not_found_in ${dir_sp_cds} ${mcmctree_params[0]}
    stop_if_species_not_found_in ${dir_sp_cds} ${mcmctree_params[1]}
  fi
  echo ""
fi

task="BUSCO analysis of species-wise CDS files"
if [[ ${run_busco} -eq 1 ]]; then
  ensure_dir "${dir_species_busco_full}"
  ensure_dir "${dir_species_busco_short}"
  species_cds_fasta=()
  mapfile -t species_cds_fasta < <(gg_find_fasta_files "${dir_sp_cds}" 1)
  echo "Number of CDS files for BUSCO: ${#species_cds_fasta[@]}"
  if [[ ${#species_cds_fasta[@]} -eq 0 ]]; then
    echo "No CDS file found. Exiting."
    exit 1
  fi
  legacy_busco_files=()
  mapfile -t legacy_busco_files < <(find "${dir_species_busco_full}" "${dir_species_busco_short}" -maxdepth 1 -type f \( -name "*_busco.full.tsv" -o -name "*_busco.short.txt" \) 2>/dev/null | sort)
  for legacy_busco_file in "${legacy_busco_files[@]}"; do
    legacy_base=$(basename "${legacy_busco_file}")
    normalized_base=$(printf '%s\n' "${legacy_base}" | sed -E 's/_busco\.full\.tsv$/.busco.full.tsv/; s/_busco\.short\.txt$/.busco.short.txt/')
    normalized_path="$(dirname "${legacy_busco_file}")/${normalized_base}"
    if [[ "${legacy_busco_file}" == "${normalized_path}" ]]; then
      continue
    fi
    if [[ -e "${normalized_path}" ]]; then
      echo "Removing legacy BUSCO duplicate file: ${legacy_busco_file}"
      rm -f "${legacy_busco_file}"
    else
      echo "Renaming legacy BUSCO file: ${legacy_busco_file} -> ${normalized_path}"
      mv_out "${legacy_busco_file}" "${normalized_path}"
    fi
  done
  # Remove stale BUSCO outputs from species not present in current species_cds inputs.
  input_species_set=$(
    for cds_full in "${species_cds_fasta[@]}"; do
      gg_species_name_from_path_or_dot "$(basename "${cds_full}")"
    done | sort -u
  )
  busco_output_files=()
  mapfile -t busco_output_files < <(find "${dir_species_busco_full}" "${dir_species_busco_short}" -maxdepth 1 -type f \( -name "*.busco.full.tsv" -o -name "*.busco.short.txt" \) 2>/dev/null | sort)
  for busco_file in "${busco_output_files[@]}"; do
    busco_base=$(basename "${busco_file}")
    busco_species=${busco_base%.busco.full.tsv}
    if [[ "${busco_species}" == "${busco_base}" ]]; then
      busco_species=${busco_base%.busco.short.txt}
    fi
    if ! printf '%s\n' "${input_species_set}" | grep -qx "${busco_species}"; then
      echo "Removing stale BUSCO output for species not in current input: ${busco_file}"
      rm -f "${busco_file}"
    fi
  done
  for cds_full in "${species_cds_fasta[@]}"; do
    cds=$(basename "${cds_full}")
    sp_ub=$(gg_species_name_from_path "${cds}")
    file_sp_busco_full="${dir_species_busco_full}/${sp_ub}.busco.full.tsv"
    file_sp_busco_short="${dir_species_busco_short}/${sp_ub}.busco.short.txt"
    if [[ ! -s ${file_sp_busco_full} || ! -s ${file_sp_busco_short} ]]; then
      gg_step_start "${task}: ${cds}"

      if [[ "${cds}" == *gz ]]; then # Do not quote like "*gz"
        echo "Decompressing gzipped CDS fasta for BUSCO: ${cds}"
        seqkit seq --threads ${NSLOTS} "${dir_sp_cds}/${cds}" --out-file "tmp.busco_input.cds.fasta"
      else
        cp_out ${dir_sp_cds}/${cds} ./tmp.busco_input.cds.fasta
      fi

      dir_busco_db=$(ensure_busco_download_path "${dir_pg}" "${busco_lineage}")
      if [[ $? -ne 0 ]]; then
        echo "Failed to prepare BUSCO dataset: ${busco_lineage}"
        exit 1
      fi
      dir_busco_lineage="${dir_busco_db}/lineages/${busco_lineage}"

      busco \
      --in "tmp.busco_input.cds.fasta" \
      --mode "transcriptome" \
      --out "busco_tmp" \
      --cpu ${NSLOTS} \
      --force \
      --evalue 1e-03 \
      --limit 20 \
      --lineage_dataset ${dir_busco_lineage} \
      --download_path ${dir_busco_db} \
      --offline

      if [[ $? -eq 0 ]]; then
        if copy_busco_tables "./busco_tmp" "${busco_lineage}" "${file_sp_busco_full}" "${file_sp_busco_short}"; then
          rm -r './busco_tmp'
        else
          echo "Failed to locate normalized BUSCO outputs for ${sp_ub}. Exiting."
          exit 1
        fi
      fi
    else
      echo "Skipped BUSCO: ${cds}"
    fi
  done
  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

task="Collecting IDs of common BUSCO genes"
if ! is_species_set_identical "${dir_sp_cds}" "${dir_species_busco_full}"; then
  echo "Exiting due to species-set mismatch between ${dir_sp_cds} and ${dir_species_busco_full}"
  exit 1
fi
if [[ ! -s ${file_busco_summary_table} && ${run_get_busco_summary} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_parent_dir "${file_busco_summary_table}"

  python ${dir_myscript}/collect_common_BUSCO_genes.py \
  --busco_outdir ${dir_species_busco_full} \
  --ncpu ${NSLOTS} \
  --outfile "tmp.busco_summary_table.tsv"
  mv_out "tmp.busco_summary_table.tsv" ${file_busco_summary_table}

  num_busco_ids=$(( $(wc -l < "${file_busco_summary_table}") - 1 ))
  echo "Number of BUSCO genes: ${num_busco_ids}"
  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

task="Generating fasta files for individual single-copy genes"
ensure_dir "${dir_single_copy_fasta}"
legacy_singlecopy_files=()
mapfile -t legacy_singlecopy_files < <(find "${dir_single_copy_fasta}" -maxdepth 1 -type f -name "*.cds.fasta" | sort)
for legacy_singlecopy_file in "${legacy_singlecopy_files[@]}"; do
  migrated_singlecopy_file="${legacy_singlecopy_file%.fasta}.fa.gz"
  if [[ ! -s "${migrated_singlecopy_file}" ]]; then
    seqkit seq --threads 1 "${legacy_singlecopy_file}" --out-file "${migrated_singlecopy_file}"
  fi
  rm -f "${legacy_singlecopy_file}"
done
num_busco_ids=$(( $(wc -l < "${file_busco_summary_table}") - 1 ))
singlecopy_fasta_files=()
mapfile -t singlecopy_fasta_files < <(gg_find_file_basenames "${dir_single_copy_fasta}" "*.cds.fa.gz")
num_singlecopy_fasta=${#singlecopy_fasta_files[@]}
if [[ ! ${num_busco_ids} -eq ${num_singlecopy_fasta} && ${run_individual_get_fasta} -eq 1 ]]; then
  gg_step_start "${task}"

  generate_single_copy_fasta() {
    local busco_id
    busco_id=$(awk -v row="$1" 'NR==row {print $1; exit}' "${file_busco_summary_table}")
    local remove_nonsingle=$2
    local outfile1="${dir_single_copy_fasta}/${busco_id}.cds.fa.gz"
    if [[ -s ${outfile1} ]]; then
      return 0
    fi
    local genes=()
    IFS=$'\t' read -r -a genes <<< "$(sed -n "${1}P" "${file_busco_summary_table}" | cut -f 4-)"
    echo "busco_id: ${busco_id}"
    if [[ ${strictly_single_copy_only} -eq 1 ]]; then
      if [[ "${genes[@]}" == *" -"* || "${genes[@]}" == *"- "* ]]; then
        echo "Skipping. ${busco_id} has missing gene(s)."
        return 0
      fi
      if [[ "${genes[@]}" == *","* ]]; then
        echo "Skipping. ${busco_id} has duplicated gene(s)."
        return 0
      fi
    fi
    local genes1=()
    read -r -a genes1 <<< "$(printf "%s " "${genes[@]}" | sed -e "s/[[:space:]]-[[:space:]]/ /g" -e "s/[^[[:space:]]]*,[^[[:space:]]]*//g" -e "s/[[:space:]]+/ /g")"
    local num_gene=${#genes1[@]}
    if [[ ${num_gene} -eq 0 ]]; then
      echo "Skipping. ${busco_id} has no genes in the selected species."
      return 0
    fi
    if [[ ! -s ${outfile1} ]]; then
      local pattern_args=()
      for gene in "${genes1[@]}"; do
        pattern_args+=(--pattern "${gene}")
      done
      seqkit grep --threads 1 "${pattern_args[@]}" --infile-list "species_cds_fasta_list.txt" \
      | seqkit replace --pattern X --replacement N --by-seq --ignore-case --threads 1 \
      | seqkit replace --pattern " .*" --replacement "" --ignore-case --threads 1 \
      | cdskit pad \
      | sed -e "s/_/|/" -e "s/_.*//" -e "s/|/_/" \
      | seqkit seq --threads 1 --out-file "${outfile1}"
      if [[ ! -s ${outfile1} ]]; then
        echo "File is empty. Removing: ${outfile1}"
        rm ${outfile1}
      fi
    fi
    local fasta_genes=()
    if [[ -s ${outfile1} ]]; then
      mapfile -t fasta_genes < <(seqkit seq --name --threads 1 "${outfile1}")
    fi
    local num_seq=${#fasta_genes[@]}
    # this block needs to be disabeld for ${strictly_single_copy_only} -eq 1, because orthogroups won't be complete
    if [[ ${num_gene} -ne ${num_seq} && ${strictly_single_copy_only} -eq 1 ]]; then
      echo "${busco_id}: Error. Number of genes and sequences did not match."
      echo "Genes in the orthogroup or BLAST hit:"
      echo -e "${genes1[@]}"
      echo ""
      echo "Genes in the generated FASTA:"
      echo -e "${fasta_genes[@]}"
      echo ""
      echo "Check duplicated sequence names in the species_cds. Exiting."
      rm ${outfile1}
      exit 1
    fi
  }

  gg_find_fasta_files "${dir_sp_cds}" 1 > species_cds_fasta_list.txt
  num_busco_ids=$(( $(wc -l < "${file_busco_summary_table}") - 1 ))
  for (( i=2; i<=num_busco_ids+1; i++ )); do # starting from 2 because the line 1 is header.
    wait_until_jobn_le ${NSLOTS}
    generate_single_copy_fasta ${i} ${strictly_single_copy_only} &
  done
  wait
  rm -f tmp.*
else
  gg_step_skip "${task}"
fi

task="In-frame mafft alignment"
ensure_dir "${dir_single_copy_mafft}"
legacy_mafft_files=()
mapfile -t legacy_mafft_files < <(find "${dir_single_copy_mafft}" -maxdepth 1 -type f -name "*.cds.aln.fasta" | sort)
for legacy_mafft_file in "${legacy_mafft_files[@]}"; do
  migrated_mafft_file="${legacy_mafft_file%.fasta}.fa.gz"
  if [[ ! -s "${migrated_mafft_file}" ]]; then
    seqkit seq --threads 1 "${legacy_mafft_file}" --out-file "${migrated_mafft_file}"
  fi
  rm -f "${legacy_mafft_file}"
done
num_busco_ids=$(( $(wc -l < "${file_busco_summary_table}") - 1 ))
mafft_fasta_files=()
mapfile -t mafft_fasta_files < <(gg_find_file_basenames "${dir_single_copy_mafft}" "*.cds.aln.fa.gz")
num_mafft_fasta=${#mafft_fasta_files[@]}
if [[ ! ${num_busco_ids} -eq ${num_mafft_fasta} && ${run_individual_mafft} -eq 1 ]]; then
  gg_step_start "${task}"

  run_mafft() {
    local infile=$1
    local infile_base=${infile%%.*}
    local outfile=${dir_single_copy_mafft}/${infile_base}.cds.aln.fa.gz
    if [[ -s ${outfile} ]]; then
      return 0
    fi
    local infile_path="${dir_single_copy_fasta}/${infile}"
    local num_seq=$(gg_count_fasta_records "${infile_path}")
    if [[ ${num_seq} -lt 2 ]]; then
      echo "Skipped. At least 2 sequences are necessary for MAFFT: ${infile}"
      return 0
    fi
    echo "$(date): start mafft: ${infile_base}"
    seqkit seq --threads 1 "${infile_path}" --out-file "tmp.${infile_base}.input.cds.fasta"
    cdskit mask \
    --seqfile "tmp.${infile_base}.input.cds.fasta" \
    --outfile tmp.${infile_base}.cds.fasta
    seqkit translate \
    --allow-unknown-codon \
    --transl-table ${genetic_code} \
    --threads 1 \
    tmp.${infile_base}.cds.fasta \
    > tmp.${infile_base}.pep.fasta
    mafft \
    --auto \
    --thread 1 \
    tmp.${infile_base}.pep.fasta \
    > tmp.${infile_base}.pep.aln.fasta
    cdskit backalign \
    --seqfile tmp.${infile_base}.cds.fasta \
    --aa_aln tmp.${infile_base}.pep.aln.fasta \
    --codontable ${genetic_code} \
    --outfile tmp.${infile_base}.cds.aln.fasta
    if [[ -s tmp.${infile_base}.cds.aln.fasta ]]; then
      seqkit seq --threads 1 tmp.${infile_base}.cds.aln.fasta --out-file "tmp.${infile_base}.cds.aln.out.fa.gz"
      mv_out "tmp.${infile_base}.cds.aln.out.fa.gz" ${outfile}
    fi
    rm tmp.${infile_base}*
  }

  input_alignment_files=()
  mapfile -t input_alignment_files < <(gg_find_file_basenames "${dir_single_copy_fasta}" "*.cds.fa.gz")
  echo "Number of input alignments: ${#input_alignment_files[@]}"
  for input_alignment_file in "${input_alignment_files[@]}"; do
    wait_until_jobn_le ${NSLOTS}
    run_mafft ${input_alignment_file} &
  done
  wait
else
  gg_step_skip "${task}"
fi

task="TrimAl"
ensure_dir "${dir_single_copy_trimal}"
legacy_trimal_files=()
mapfile -t legacy_trimal_files < <(find "${dir_single_copy_trimal}" -maxdepth 1 -type f -name "*.trimal.fasta" | sort)
for legacy_trimal_file in "${legacy_trimal_files[@]}"; do
  migrated_trimal_file="${legacy_trimal_file%.fasta}.fa.gz"
  if [[ ! -s "${migrated_trimal_file}" ]]; then
    seqkit seq --threads 1 "${legacy_trimal_file}" --out-file "${migrated_trimal_file}"
  fi
  rm -f "${legacy_trimal_file}"
done
trimal_fasta_files=()
mapfile -t trimal_fasta_files < <(gg_find_file_basenames "${dir_single_copy_trimal}" "*.trimal.fa.gz")
num_trimal_fasta=${#trimal_fasta_files[@]}
if [[ ! ${num_busco_ids} -eq ${num_trimal_fasta} && ${run_individual_trimal} -eq 1 ]]; then
  gg_step_start "${task}"

  run_trimal() {
    local infile=$1
    local infile_base=${infile%%.*}
    local outfile="${dir_single_copy_trimal}/${infile_base}.trimal.fa.gz"
    if [[ ! -s ${outfile} ]]; then
      seqkit seq --remove-gaps --threads 1 ${dir_single_copy_mafft}/${infile} > tmp.${infile_base}.degap.fasta
      seqkit translate --transl-table ${genetic_code} --threads 1 ${dir_single_copy_mafft}/${infile} > tmp.${infile_base}.pep.fasta
      trimal \
      -in tmp.${infile_base}.pep.fasta \
      -backtrans tmp.${infile_base}.degap.fasta \
      -out tmp.${infile_base}.trimal.fasta \
      -ignorestopcodon \
      -automated1
      if [[ -s tmp.${infile_base}.trimal.fasta ]]; then
        seqkit seq --threads 1 tmp.${infile_base}.trimal.fasta --out-file "tmp.${infile_base}.trimal.out.fa.gz"
        mv_out "tmp.${infile_base}.trimal.out.fa.gz" ${outfile}
      fi
      rm tmp.${infile_base}.*
    fi
  }

  input_alignment_files=()
  mapfile -t input_alignment_files < <(gg_find_file_basenames "${dir_single_copy_mafft}" "*.cds.aln.fa.gz")
  echo "Number of input alignments: ${#input_alignment_files[@]}"
  for input_alignment_file in "${input_alignment_files[@]}"; do
    wait_until_jobn_le ${NSLOTS}
    run_trimal ${input_alignment_file} &
  done
  wait
else
  gg_step_skip "${task}"
fi

task="Concatenating single-copy CDS fasta files"
if [[ ( ! -s ${file_concat_cds} || ! -s ${file_concat_pep} ) && ${run_concat_alignment} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_parent_dir "${file_concat_cds}"
  ensure_parent_dir "${file_concat_pep}"
  ensure_dir "${dir_concat_fasta}"
  mapfile -t trimal_files < <(find "${dir_single_copy_trimal}" -maxdepth 1 -type f -name "*.trimal.fa.gz" | sort)
  if [[ ${#trimal_files[@]} -eq 0 ]]; then
    echo "No trimmed single-copy CDS fasta files were found in: ${dir_single_copy_trimal}"
    exit 1
  fi

  # Build partition ranges from per-locus alignment lengths.
  partition_start=1
  > "${dir_concat_fasta}/partitions.txt"
  for trimal_file in "${trimal_files[@]}"; do
    aln_len=$(seqkit fx2tab -l "${trimal_file}" | awk 'NR==1 {print $NF; exit}')
    if [[ -z "${aln_len}" ]]; then
      echo "Failed to read alignment length from: ${trimal_file}"
      exit 1
    fi
    partition_end=$((partition_start + aln_len - 1))
    echo "DNA, $(basename "${trimal_file}") = ${partition_start}-${partition_end}" >> "${dir_concat_fasta}/partitions.txt"
    partition_start=$((partition_end + 1))
  done

  if [[ ${strictly_single_copy_only} -eq 0 ]]; then
    concat_cds_tmp="tmp.concat.cds.fa.gz"
    seqkit concat --full --fill "-" --threads ${NSLOTS} "${trimal_files[@]}" \
    | seqkit sort --threads ${NSLOTS} \
    | seqkit seq --threads ${NSLOTS} --out-file "${concat_cds_tmp}"
  else
    concat_cds_tmp="tmp.concat.cds.fa.gz"
    seqkit concat --threads ${NSLOTS} "${trimal_files[@]}" \
    | seqkit sort --threads ${NSLOTS} \
    | seqkit seq --threads ${NSLOTS} --out-file "${concat_cds_tmp}"
  fi
  mv_out "${concat_cds_tmp}" "${file_concat_cds}"

  seqkit translate --transl-table ${genetic_code} --threads ${NSLOTS} "${file_concat_cds}" \
  | seqkit seq --threads ${NSLOTS} --out-file "tmp.concat.pep.fa.gz"
  mv_out "tmp.concat.pep.fa.gz" "${file_concat_pep}"
else
  gg_step_skip "${task}"
fi

task="IQ-TREE of the concatenated alignment with a Protein evolution model"
if [[ ! -s ${file_concat_iqtree_pep} && ${run_concat_iqtree_protein} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_concat_iqtree_pep}"

  ntaxa=$(gg_count_fasta_records "${file_concat_pep}")
  if [[ ${ntaxa} -lt 4 ]]; then
    bootstrap_params=''
  fi
  iqtree_mem_arg=''
  if [[ -n "${MEM_PER_HOST:-}" ]]; then
    iqtree_mem_arg="-mem ${MEM_PER_HOST}G"
  fi

  cd ${dir_concat_iqtree_pep}
  concat_pep_local="tmp.concat.pep.input.fasta"
  seqkit seq --threads 1 "${file_concat_pep}" --out-file "./${concat_pep_local}"
  iqtree \
  -s "${concat_pep_local}" \
  -st AA \
  -m ${protein_model} \
  -pre "${concat_pep_local}" \
  -nt ${NSLOTS} \
  ${iqtree_mem_arg} \
  -seed 12345 \
  ${bootstrap_params}; iqtree_exit_status=$?

  if [[ ${iqtree_exit_status} -ne 0 ]]; then
    echo "Error. IQ-TREE exited with non-zero status: ${iqtree_exit_status}"
  else
    echo "IQ-TREE successfully exited with zero status: ${iqtree_exit_status}"
    if [[ ${ntaxa} -lt 4 ]]; then
      out_nwk="${concat_pep_local}.treefile"
    else
      out_nwk="${concat_pep_local}.contree"
    fi
    cp_out ${out_nwk} ${file_concat_iqtree_pep}
    rm "./${concat_pep_local}"
  fi
  cd ${dir_tmp}
else
  gg_step_skip "${task}"
fi

task="Rooting of IQ-TREE's protein tree"
if [[ ! -s ${file_concat_iqtree_pep_root} && ${run_concat_iqtree_protein} -eq 1 ]]; then
  gg_step_start "${task}"

  root_tree_with_outgroup \
  "${file_concat_iqtree_pep}" \
  "${file_concat_iqtree_pep_root}" \
  "concatenated protein tree"

  if [[ -s ${file_concat_iqtree_pep_root} ]]; then
    Rscript ${dir_myscript}/nwk2pdf.r --underbar2space=yes --italic=yes --infile=${file_concat_iqtree_pep_root}
  fi

  if [[ -s ${file_concat_iqtree_pep_root} && ${undated_species_tree} == 'iqtree_pep' ]]; then
    echo "Undated species tree is copied,"
    echo "from: ${file_concat_iqtree_pep_root}"
    echo "to: ${file_undated_species_tree}"
    cp_out "${file_concat_iqtree_pep_root}" "${file_undated_species_tree}"
    Rscript ${dir_myscript}/nwk2pdf.r --underbar2space=yes --italic=yes --infile=${file_undated_species_tree}
  fi
else
  gg_step_skip "${task}"
fi

task="IQ-TREE of the concatenated alignment with a DNA evolution model"
if [[ ! -s ${file_concat_iqtree_dna} && ${run_concat_iqtree_dna} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_concat_iqtree_dna}"

  ntaxa=$(gg_count_fasta_records "${file_concat_cds}")
  if [[ ${ntaxa} -lt 4 ]]; then
    bootstrap_params=''
  fi
  iqtree_mem_arg=''
  if [[ -n "${MEM_PER_HOST:-}" ]]; then
    iqtree_mem_arg="-mem ${MEM_PER_HOST}G"
  fi

  cd ${dir_concat_iqtree_dna}
  concat_cds_local="tmp.concat.cds.input.fasta"
  seqkit seq --threads 1 "${file_concat_cds}" --out-file "./${concat_cds_local}"
  iqtree \
  -s "${concat_cds_local}" \
  -st DNA \
  -m ${nucleotide_model} \
  -pre "${concat_cds_local}" \
  -nt ${NSLOTS} \
  ${iqtree_mem_arg} \
  -seed 12345 \
  ${bootstrap_params}; iqtree_exit_status=$?

  if [[ ${iqtree_exit_status} -ne 0 ]]; then
    echo "Error. IQ-TREE exited with non-zero status: ${iqtree_exit_status}"
  else
    echo "IQ-TREE successfully exited with zero status: ${iqtree_exit_status}"
    if [[ ${ntaxa} -lt 4 ]]; then
      out_nwk="${concat_cds_local}.treefile"
    else
      out_nwk="${concat_cds_local}.contree"
    fi
    cp_out ${out_nwk} ${file_concat_iqtree_dna}
    rm "./${concat_cds_local}"
  fi
  cd ${dir_tmp}
else
  gg_step_skip "${task}"
fi

task="Rooting of IQ-TREE's DNA tree"
if [[ ! -s ${file_concat_iqtree_dna_root} && ${run_concat_iqtree_dna} -eq 1 ]]; then
  gg_step_start "${task}"

  root_tree_with_outgroup \
  "${file_concat_iqtree_dna}" \
  "${file_concat_iqtree_dna_root}" \
  "concatenated DNA tree"

  if [[ -s ${file_concat_iqtree_dna_root} ]]; then
    Rscript ${dir_myscript}/nwk2pdf.r --underbar2space=yes --italic=yes --infile=${file_concat_iqtree_dna_root}
  fi

  if [[ -s ${file_concat_iqtree_dna_root} && ${undated_species_tree} == 'iqtree_dna' ]]; then
    echo "Undated species tree is copied,"
    echo "from: ${file_concat_iqtree_dna_root}"
    echo "to: ${file_undated_species_tree}"
    cp_out "${file_concat_iqtree_dna_root}" "${file_undated_species_tree}"
    Rscript ${dir_myscript}/nwk2pdf.r --underbar2space=yes --italic=yes --infile=${file_undated_species_tree}
  fi
else
  gg_step_skip "${task}"
fi

task="IQ-TREE for individual single-copy protein trees"
ensure_dir "${dir_single_copy_iqtree_pep}"
iqtree_pep_tree_files=()
mapfile -t iqtree_pep_tree_files < <(gg_find_file_basenames "${dir_single_copy_iqtree_pep}" "*.nwk")
num_iqtree_pep=${#iqtree_pep_tree_files[@]}
if [[ ! ${num_busco_ids} -eq ${num_iqtree_pep} && ${run_individual_iqtree_pep} -eq 1 ]]; then
  gg_step_start "${task}"

  run_iqtree_pep() {
    local infile=$1
    local infile_base=${infile%%.*}
    local outfile="${dir_single_copy_iqtree_pep}/${infile_base}.pep.nwk"
    if [[ -s ${outfile} ]]; then
      return 0
    fi
    local num_seq=$(gg_count_fasta_records "${dir_single_copy_trimal}/${infile}")
    if [[ ${num_seq} -lt 3 ]]; then
      echo "Skipped. At least 3 sequences are necessary for IQ-TREE: ${infile}"
      return 0
    fi
    seqkit translate --transl-table ${genetic_code} --threads 1 ${dir_single_copy_trimal}/${infile} \
    > "tmp.${infile_base}.pep.fasta"
    iqtree \
    -s "tmp.${infile_base}.pep.fasta" \
    -m ${protein_model} \
    -T 1 \
    --prefix tmp.${infile_base} \
    --seed 12345 \
    --redo
    mv_out tmp.${infile_base}.treefile ${outfile}
    rm tmp.${infile_base}.*
  }

  input_alignment_files=()
  mapfile -t input_alignment_files < <(gg_find_file_basenames "${dir_single_copy_trimal}" "*.trimal.fa.gz")
  echo "Number of input alignments: ${#input_alignment_files[@]}"
  for input_alignment_file in "${input_alignment_files[@]}"; do
    wait_until_jobn_le ${NSLOTS}
    run_iqtree_pep ${input_alignment_file} &
  done
  wait
else
  gg_step_skip "${task}"
fi

task="ASTRAL of individual single-copy protein trees"
if [[ ( ! -s ${file_astral_tree_pep} || ! -s ${file_astral_log_pep} ) && ${run_astral_pep} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_astral_pep}"

  if compgen -G "tmp.astral.*" > /dev/null; then
    rm tmp.astral.*
  fi

  if ! build_astral_input "${dir_single_copy_iqtree_pep}" "tmp.astral.merged.iqtree.nwk" "${astral_min_tips}"; then
    echo "Skipped. No eligible protein gene trees for ASTRAL after filtering."
  else
    astral-hybrid \
    --input "tmp.astral.merged.iqtree.nwk" \
    --output "tmp.astral.out.tree" \
    --mode 3 \
    --support 2 \
    --thread ${NSLOTS} \
    2> "tmp.astral.log.txt"

    labels=("pp1" "pp2" "pp3" "f1" "f2" "f3" "q1" "q2" "q3") # https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md
    for i in "${!labels[@]}"; do
      tmp_label_tree="tmp.astral.pep.${labels[i]}.tmp.nwk"
      python ${dir_myscript}/extract_astral_support_label.py \
      --infile "tmp.astral.out.tree" \
      --outfile "${tmp_label_tree}" \
      --label_key "${labels[i]}"
      root_tree_with_outgroup \
      "${tmp_label_tree}" \
      "single_copy.astral.pep.${labels[i]}.nwk" \
      "ASTRAL protein tree (${labels[i]})"
      if [[ -s "single_copy.astral.pep.${labels[i]}.nwk" ]]; then
        Rscript ${dir_myscript}/nwk2pdf.r --underbar2space=yes --italic=yes --infile="single_copy.astral.pep.${labels[i]}.nwk"
      fi
      rm -f "${tmp_label_tree}"
    done

    if [[ -s "single_copy.astral.pep.q1.nwk" ]]; then
      mv_out single_copy.astral.pep.* "${dir_astral_pep}"
      mv_out "tmp.astral.log.txt" "${file_astral_log_pep}"
      echo "For more information on support values (e.g., f1, f2, pp1, q1, ...), please refer to: https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md" > "${dir_astral_pep}/README.txt"
    fi

    if [[ -s "${file_astral_tree_pep_q1}" ]]; then
      if optimize_astral_tree_branch_lengths "${file_astral_tree_pep_q1}" "${file_concat_pep}" "${protein_model}" "${file_astral_tree_pep}" "pep"; then
        Rscript ${dir_myscript}/nwk2pdf.r --underbar2space=yes --italic=yes --infile="${file_astral_tree_pep}"
      else
        echo "Warning: Falling back to unoptimized ASTRAL protein tree."
        cp_out "${file_astral_tree_pep_q1}" "${file_astral_tree_pep}"
      fi
    fi
  fi
  rm tmp.astral.*

  if [[ -s ${file_astral_tree_pep} && ${undated_species_tree} == 'astral_pep' ]]; then
    echo "Undated species tree is copied,"
    echo "from: ${file_astral_tree_pep}"
    echo "to: ${file_undated_species_tree}"
    nwkit drop --name yes --target intnode --infile "${file_astral_tree_pep}" | nwkit label --target intnode --force yes --outfile "tmp.undated_species_tree.nwk"
    if [[ -s "tmp.undated_species_tree.nwk" ]]; then
      mv_out "tmp.undated_species_tree.nwk" "${file_undated_species_tree}"
    fi
    if [[ ! -s "${file_undated_species_tree}" ]]; then
      echo "Warning: Failed to convert optimized ASTRAL protein tree with nwkit drop|label. Copying optimized tree instead."
      cp_out "${file_astral_tree_pep}" "${file_undated_species_tree}"
    fi
    if [[ -s "${file_undated_species_tree}" ]]; then
      Rscript ${dir_myscript}/nwk2pdf.r --underbar2space=yes --italic=yes --infile=${file_undated_species_tree}
    fi
  fi
else
  gg_step_skip "${task}"
fi

task="IQ-TREE for individual single-copy DNA trees"
ensure_dir "${dir_single_copy_iqtree_dna}"
iqtree_dna_tree_files=()
mapfile -t iqtree_dna_tree_files < <(gg_find_file_basenames "${dir_single_copy_iqtree_dna}" "*.nwk")
num_iqtree_dna=${#iqtree_dna_tree_files[@]}
if [[ ! ${num_busco_ids} -eq ${num_iqtree_dna} && ${run_individual_iqtree_dna} -eq 1 ]]; then
	gg_step_start "${task}"

  run_iqtree_dna() {
    local infile=$1
    local infile_base=${infile%%.*}
    local outfile="${dir_single_copy_iqtree_dna}/${infile_base}.dna.nwk"
    if [[ -s ${outfile} ]]; then
      return 0
    fi
    local num_seq=$(gg_count_fasta_records "${dir_single_copy_trimal}/${infile}")
    if [[ ${num_seq} -lt 3 ]]; then
      echo "Skipped. At least 3 sequences are necessary for IQ-TREE: ${infile}"
      return 0
    fi
    cp_out "${dir_single_copy_trimal}/${infile}" "./tmp.${infile_base}.input.fasta"
    iqtree \
    -s "./tmp.${infile_base}.input.fasta" \
    -m ${nucleotide_model} \
    -T 1 \
    --prefix tmp.${infile_base} \
    --seed 12345 \
    --redo
    mv_out tmp.${infile_base}.treefile ${outfile}
    rm tmp.${infile_base}.*
  }

  input_alignment_files=()
  mapfile -t input_alignment_files < <(gg_find_file_basenames "${dir_single_copy_trimal}" "*.trimal.fa.gz")
  echo "Number of input alignments: ${#input_alignment_files[@]}"
  for input_alignment_file in "${input_alignment_files[@]}"; do
    wait_until_jobn_le ${NSLOTS}
    run_iqtree_dna ${input_alignment_file} &
  done
  wait
else
	gg_step_skip "${task}"
fi

task="ASTRAL of individual single-copy DNA trees"
if [[ ( ! -s ${file_astral_tree_dna} || ! -s ${file_astral_log_dna} ) && ${run_astral_dna} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_astral_dna}"

  if compgen -G "tmp.astral.*" > /dev/null; then
    rm tmp.astral.*
  fi

  if ! build_astral_input "${dir_single_copy_iqtree_dna}" "tmp.astral.merged.iqtree.nwk" "${astral_min_tips}"; then
    echo "Skipped. No eligible DNA gene trees for ASTRAL after filtering."
  else
    astral-hybrid \
    --input "tmp.astral.merged.iqtree.nwk" \
    --output "tmp.astral.out.tree" \
    --mode 3 \
    --support 2 \
    --thread ${NSLOTS} \
    2> "tmp.astral.log.txt"

    labels=("pp1" "pp2" "pp3" "f1" "f2" "f3" "q1" "q2" "q3") # https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md
    for i in "${!labels[@]}"; do
      tmp_label_tree="tmp.astral.dna.${labels[i]}.tmp.nwk"
      python ${dir_myscript}/extract_astral_support_label.py \
      --infile "tmp.astral.out.tree" \
      --outfile "${tmp_label_tree}" \
      --label_key "${labels[i]}"
      root_tree_with_outgroup \
      "${tmp_label_tree}" \
      "single_copy.astral.dna.${labels[i]}.nwk" \
      "ASTRAL DNA tree (${labels[i]})"
      if [[ -s "single_copy.astral.dna.${labels[i]}.nwk" ]]; then
        Rscript ${dir_myscript}/nwk2pdf.r --underbar2space=yes --italic=yes --infile="single_copy.astral.dna.${labels[i]}.nwk"
      fi
      rm -f "${tmp_label_tree}"
    done

    if [[ -s "single_copy.astral.dna.q1.nwk" ]]; then
      mv_out single_copy.astral.dna.* "${dir_astral_dna}"
      mv_out "tmp.astral.log.txt" "${file_astral_log_dna}"
      echo "For more information on support values (e.g., f1, f2, pp1, q1, ...), please refer to: https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md" > "${dir_astral_dna}/README.txt"
    fi

    if [[ -s "${file_astral_tree_dna_q1}" ]]; then
      if optimize_astral_tree_branch_lengths "${file_astral_tree_dna_q1}" "${file_concat_cds}" "${nucleotide_model}" "${file_astral_tree_dna}" "dna"; then
        Rscript ${dir_myscript}/nwk2pdf.r --underbar2space=yes --italic=yes --infile="${file_astral_tree_dna}"
      else
        echo "Warning: Falling back to unoptimized ASTRAL DNA tree."
        cp_out "${file_astral_tree_dna_q1}" "${file_astral_tree_dna}"
      fi
    fi
  fi
  rm tmp.astral.*

  if [[ -s ${file_astral_tree_dna} && ${undated_species_tree} == 'astral_dna' ]]; then
    echo "Undated species tree is copied,"
    echo "from: ${file_astral_tree_dna}"
    echo "to: ${file_undated_species_tree}"
    nwkit drop --name yes --target intnode --infile "${file_astral_tree_dna}" | nwkit label --target intnode --force yes --outfile "tmp.undated_species_tree.nwk"
    if [[ -s "tmp.undated_species_tree.nwk" ]]; then
      mv_out "tmp.undated_species_tree.nwk" "${file_undated_species_tree}"
    fi
    if [[ ! -s "${file_undated_species_tree}" ]]; then
      echo "Warning: Failed to convert optimized ASTRAL DNA tree with nwkit drop|label. Copying optimized tree instead."
      cp_out "${file_astral_tree_dna}" "${file_undated_species_tree}"
    fi
    if [[ -s "${file_undated_species_tree}" ]]; then
      Rscript ${dir_myscript}/nwk2pdf.r --underbar2space=yes --italic=yes --infile=${file_undated_species_tree}
    fi
  fi
else
  gg_step_skip "${task}"
fi

task="Species tree plotting"
disable_if_no_input_file "run_plot_species_trees" ${file_concat_iqtree_dna_root} ${file_concat_iqtree_pep_root} ${file_astral_tree_dna} ${file_astral_tree_pep}
if [[ ! -s ${file_plot_species_trees} && ${run_plot_species_trees} -eq 1 ]]; then
  gg_step_start "${task}"

  Rscript ${dir_myscript}/plot_species_trees.r \
  --iqtree_dna_nwk=${file_concat_iqtree_dna_root} \
  --iqtree_pep_nwk=${file_concat_iqtree_pep_root} \
  --iqtree_dna_log=${dir_concat_iqtree_dna}/"tmp.concat.cds.input.fasta.log" \
  --iqtree_pep_log=${dir_concat_iqtree_pep}/"tmp.concat.pep.input.fasta.log" \
  --astral_dna_nwk=${file_astral_tree_dna} \
  --astral_pep_nwk=${file_astral_tree_pep} \
  --astral_dna_log=${file_astral_log_dna} \
  --astral_pep_log=${file_astral_log_pep}

  if [[ -s "species_trees.pdf" ]]; then
    echo "Output file found for the task: ${task}"
    mv_out "species_trees.pdf" ${file_plot_species_trees}
  fi
else
  gg_step_skip "${task}"
fi

task="Time-constrained tree preparation"
disable_if_no_input_file "run_constrained_tree" ${file_undated_species_tree}
if [[ ! -s ${file_constrained_tree} && ${run_constrained_tree} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_parent_dir "${file_constrained_tree}"
  if [[ ${timetree_constraint} -eq 1 ]]; then
    if ! ensure_ete_taxonomy_db "${dir_pg}"; then
      echo "Error: Failed to prepare ETE taxonomy database under workspace/db."
      echo "Error: Network access is required to download/update ETE taxonomy data."
      exit 1
    fi
    nwkit mcmctree \
    --infile ${file_undated_species_tree} \
    --timetree "ci" \
    --min_clade_prop 0.2 \
    --outfile "tmp.constrained.tree.nwk"
    if [[ -s "tmp.constrained.tree.nwk" ]]; then
      mv_out "tmp.constrained.tree.nwk" "${file_constrained_tree}"
    fi
  else
    tree_string=$(< "${file_undated_species_tree}")
    for i in $(seq 0 $((${#mcmctree_divergence_time_constraints[@]} - 1))); do
      mdtc=${mcmctree_divergence_time_constraints[${i}]}
      echo "applying ${mdtc}"
      mcmctree_params=(${mdtc//,/ })
      bound_params=""
      if [[ "${mcmctree_params[2]}" != "-" ]]; then
        bound_params="${bound_params} --lower_bound ${mcmctree_params[2]}"
      fi
      if [[ "${mcmctree_params[3]}" != "-" ]]; then
        bound_params="${bound_params} --upper_bound ${mcmctree_params[3]}"
      fi
      left_right="--left_species ${mcmctree_params[0]} --right_species ${mcmctree_params[1]}"
      echo "nwkit mcmctree params: ${left_right} ${bound_params}"
      tree_string=$(printf '%s\n' "${tree_string}" | nwkit mcmctree ${left_right} ${bound_params})
    done
    echo -e "${tree_string}" > "tmp.constrained.tree.nwk"
    if [[ -s "tmp.constrained.tree.nwk" ]]; then
      mv_out "tmp.constrained.tree.nwk" "${file_constrained_tree}"
    fi
  fi
  if [[ -s ${file_constrained_tree} ]]; then
    if ! normalize_iq2mc_constraint_tree ${file_constrained_tree}; then
      echo "Error: Failed to normalize the constrained tree for IQ2MC."
      exit 1
    fi
  fi
else
  gg_step_skip "${task}"
fi

task="Constrained range plotting"
disable_if_no_input_file "run_plot_constrained_tree" ${file_constrained_tree}
if [[ ! -s ${file_plot_constrained_tree} && ${run_plot_constrained_tree} -eq 1 ]]; then
  gg_step_start "${task}"
  Rscript ${dir_myscript}/plot_constrained_tree.r \
  --infile="${file_constrained_tree}" \
  --outfile="tmp.constrained_tree_plot.pdf"
  if [[ -s "tmp.constrained_tree_plot.pdf" ]]; then
    mv_out "tmp.constrained_tree_plot.pdf" "${file_plot_constrained_tree}"
  fi
  if [[ -s ${file_plot_constrained_tree} ]]; then
    echo "Output file found for the task: ${task}"
  fi
else
  gg_step_skip "${task}"
fi

task="IQ2MC step 2 (IQ-TREE Hessian/control generation)"
disable_if_no_input_file "run_mcmctree1" ${file_constrained_tree} ${file_concat_cds}
if [[ ( ! -s ${file_iq2mc_ctl} || ! -s ${file_iq2mc_hessian} || ! -s ${file_iq2mc_rooted_tree} || ! -s ${file_iq2mc_dummy_phy} ) && ${run_mcmctree1} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_mcmctree1}"

  if ! normalize_iq2mc_constraint_tree ${file_constrained_tree}; then
    echo "Error: Failed to normalize the constrained tree for IQ2MC."
    exit 1
  fi
  find "${dir_mcmctree1}" -mindepth 1 -maxdepth 1 -exec rm -rf {} +
  cd ${dir_mcmctree1}
  seqkit seq --threads 1 "${file_concat_cds}" --out-file "./tmp.iq2mc.concat.cds.fasta"

  ${iq2mc_binary} \
  -s "./tmp.iq2mc.concat.cds.fasta" \
  -m ${nucleotide_model} \
  -te ${file_constrained_tree} \
  --dating mcmctree \
  --mcmc-bds ${mcmc_birth_death_sampling} \
  --mcmc-clock ${mcmc_clock_model} \
  --mcmc-iter ${mcmc_burnin},${mcmc_sampfreq},${mcmc_nsample} \
  -T ${NSLOTS} \
  --prefix iq2mc
  if [[ $? -ne 0 ]]; then
    echo "Error: IQ2MC step 2 failed. Deleting generated files."
    rm -f "${file_iq2mc_prefix}".*
  elif [[ ! -s ${file_iq2mc_ctl} || ! -s ${file_iq2mc_hessian} || ! -s ${file_iq2mc_rooted_tree} || ! -s ${file_iq2mc_dummy_phy} ]]; then
    echo "Error: IQ2MC step 2 did not generate all expected files."
    rm -f "${file_iq2mc_prefix}".*
  fi
  rm -f "./tmp.iq2mc.concat.cds.fasta"
  cd ${dir_tmp}
else
  gg_step_skip "${task}"
fi

task="IQ2MC step 3 (MCMCtree dating run)"
disable_if_no_input_file "run_mcmctree2" ${file_iq2mc_ctl} ${file_iq2mc_hessian} ${file_iq2mc_rooted_tree} ${file_iq2mc_dummy_phy}
if [[ ! -s ${file_mcmctree2_raw} && ${run_mcmctree2} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_mcmctree2}"

  find "${dir_mcmctree2}" -mindepth 1 -maxdepth 1 -exec rm -rf {} +
  cd ${dir_mcmctree2}
  cp_out ${file_iq2mc_ctl} ./
  cp_out ${file_iq2mc_hessian} ./
  cp_out ${file_iq2mc_rooted_tree} ./
  cp_out ${file_iq2mc_dummy_phy} ./

  # Ensure MCMCtree emits CI-rich summaries (including 95% HPD annotations in FigTree output).
  ctl_basename="$(basename "${file_iq2mc_ctl}")"
  sed -i -e "s/FossilErrprint/FossilErr\\nprint/g" "${ctl_basename}"
  if ! grep -q "^print[[:space:]]*=" "${ctl_basename}"; then
    echo "print = 1" >> "${ctl_basename}"
  fi

  mcmctree "${ctl_basename}"
  if [[ $? -ne 0 ]]; then
    echo "Error: IQ2MC step 3 failed."
    rm -f ${file_mcmctree2_raw}
  fi
  cd ${dir_tmp}
else
  gg_step_skip "${task}"
fi

if [[ ! -s ${file_mcmctree2} && -s ${file_mcmctree2_raw} ]]; then
  echo "Generating ${file_mcmctree2} from ${file_mcmctree2_raw}"
  awk '
  /Species tree for FigTree/ {print; in_figtree=1; next}
  in_figtree && /^\(\(/ {print; count++; if (count >= 3) exit}
  ' "${file_mcmctree2_raw}" > "tmp.mcmctree2.txt"
  if [[ -s "tmp.mcmctree2.txt" ]]; then
    mv_out "tmp.mcmctree2.txt" "${file_mcmctree2}"
  fi
  if [[ ! -s "${file_mcmctree2}" ]]; then
    echo "Warning: Failed to extract FigTree content from ${file_mcmctree2_raw}. Copying raw file instead."
    cp_out "${file_mcmctree2_raw}" "${file_mcmctree2}"
  fi
fi

task="Convert tree format"
disable_if_no_input_file "run_convert_tree_format" ${file_mcmctree2}
if [[ ! -s ${file_mcmctree2_nwk} && ${run_convert_tree_format} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_parent_dir "${file_mcmctree2_nwk}"

  if grep -q -e "UTREE" "${file_mcmctree2}"; then
    grep -e "UTREE" "${file_mcmctree2}" |
      sed -e "s/.*UTREE 1 = //" -e "s/;.*/;/" \
        > "${dir_mcmctree2}/mcmctree_95CI.nwk"

    grep -e "UTREE" "${file_mcmctree2}" |
      sed -e "s/.*UTREE 1 = //" -e "s/;.*/;/" -e "s/[[:space:]]*\[&95%={[0-9.]*,[[:space:]][0-9.]*}\][[:space:]]*//g" -e "s/:[[:space:]]/:/g" \
        > "${dir_mcmctree2}/mcmctree_no95CI.nwk"
  else
    tree_line="$(awk '/^\(\(/ {line=$0} END {print line}' "${file_mcmctree2}")"
    if [[ -n "${tree_line}" ]]; then
      echo "${tree_line}" > "${dir_mcmctree2}/mcmctree_95CI.nwk"
      echo "${tree_line}" |
        sed -e "s/[[:space:]]*\[&95%={[0-9.]*,[[:space:]][0-9.]*}\][[:space:]]*//g" -e "s/:[[:space:]]/:/g" \
        > "${dir_mcmctree2}/mcmctree_no95CI.nwk"
    else
      echo "Error: Failed to detect a tree string in ${file_mcmctree2}"
      rm -f "${dir_mcmctree2}/mcmctree_95CI.nwk" "${dir_mcmctree2}/mcmctree_no95CI.nwk"
    fi
  fi

  if [[ -s "${dir_mcmctree2}/mcmctree_no95CI.nwk" ]]; then
    Rscript -e "library(ape); t=read.tree(\"${dir_mcmctree2}/mcmctree_no95CI.nwk\"); \
    t[['node.label']]=paste0('s',1:(length(t[['tip.label']])-1)); \
    write.tree(t, \"${file_mcmctree2_nwk}\")"
  else
    echo "Error: Missing mcmctree_no95CI.nwk. Skipping tree conversion."
  fi

  if [[ ! -s ${file_dated_species_tree} ]]; then
    echo "Dated species tree for gg_pipeline is not placed yet: ${file_dated_species_tree}"
    if [[ -s ${file_mcmctree2_nwk} ]]; then
      echo "Coping from: ${file_mcmctree2_nwk}"
      echo "Coping to: ${file_dated_species_tree}"
      cp_out ${file_mcmctree2_nwk} ${file_dated_species_tree}
      echo "Please manually check whether the species tree is valid."
    fi
  else
    echo "Dated species tree for gg_pipeline is already placed: ${file_dated_species_tree}"
    echo "If necessary, please replace the file with: ${file_mcmctree2_nwk}"
  fi
else
  gg_step_skip "${task}"
fi

task="Dated species tree plotting"
disable_if_no_input_file "run_plot_mcmctreer" ${file_mcmctree2_nwk}
if [[ ! -s ${file_plot_mcmctreer} && ${run_plot_mcmctreer} -eq 1 ]]; then
  gg_step_start "${task}"

  Rscript ${dir_myscript}/plot_mcmctreer.r \
  --infile="${file_mcmctree2_nwk}" \
  --outfile="tmp.plot_mcmctreer.pdf"
  if [[ -s "tmp.plot_mcmctreer.pdf" ]]; then
    mv_out "tmp.plot_mcmctreer.pdf" "${file_plot_mcmctreer}"
  fi

  if [[ -s "${file_plot_mcmctreer}" ]]; then
    echo "Output file found for the task: ${task}"
    echo "Output file: ${file_plot_mcmctreer}"
  fi
else
  gg_step_skip "${task}"
fi

remove_empty_subdirs ${dir_species_tree}
if [[ ${delete_tmp_dir} -eq 1 ]]; then
  if [[ -d "${dir_tmp}" ]]; then
    echo "Removing tmp directory: ${dir_tmp}"
    rm -rf -- "${dir_tmp}"
  fi
fi

echo "$(date): Exiting Singularity environment"
