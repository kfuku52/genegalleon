#!/usr/bin/env bash
set -eo pipefail

# Unified pipeline entrypoint (inlined stages):
# 1) gg_speciesTree_cmd.sh
# 2) gg_orthofinder_cmd.sh
# 3) gg_genomeEvolution_cmd.sh

run_stage_inline() {
  local stage_name=$1
  local stage_rc=0
  local tmp_script

  tmp_script=$(mktemp)
  cat > "${tmp_script}"
  chmod +x "${tmp_script}"

  echo "$(date): Starting stage: ${stage_name}"
  set +e
  bash "${tmp_script}"
  stage_rc=$?
  set -e

  rm -f "${tmp_script}"

  if [[ ${stage_rc} -eq 8 ]]; then
    echo "Skipped stage: ${stage_name} (exit code 8: outputs already exist)"
  elif [[ ${stage_rc} -ne 0 ]]; then
    echo "Failed stage: ${stage_name} (exit code ${stage_rc})"
    return "${stage_rc}"
  fi

  echo "$(date): Ending stage: ${stage_name}"
  echo ""
}

run_stage_inline "gg_speciesTree_cmd.sh" <<'__GG_INLINE_STAGE_SPECIESTREE_20260221__'
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
__GG_INLINE_STAGE_SPECIESTREE_20260221__

run_stage_inline "gg_orthofinder_cmd.sh" <<'__GG_INLINE_STAGE_ORTHOFINDER_20260221__'
#!/usr/bin/env bash

### Start: Modify this block to tailor your analysis ###

run_cds_translation=1 # Translating species_cds to species_protein. Necessary for OrthoFinder
run_orthofinder=1 # OrthoFinder run
run_sonicparanoid=1 # SonicParanoid run
run_og_selection=1 # Selecting orthogroups for downstream analyses including gg_gene_evolution
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

orthofinder_output_directory_cleanup() {
  local target_dir=$1
  local _threads=${2:-1}
  if [[ -d "${target_dir}" ]]; then
    remove_empty_subdirs "${target_dir}"
  fi
}

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
    cp_out -r "${dir_orthogroup_selection_input}/." "${dir_orthofinder_og}"
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
__GG_INLINE_STAGE_ORTHOFINDER_20260221__

run_stage_inline "gg_genomeEvolution_cmd.sh" <<'__GG_INLINE_STAGE_GENOMEEVOLUTION_20260221__'
#!/usr/bin/env bash

#run a busco analysis for all species in species_cds and extract shared complete BUSCO genes among all species.
#then make cds fasta files including all species' cds sequenes for each BUSCO genes and make alignment
#then concatenate all fasta files for further downstream analysis -> task = trimal in gg_species_tree.

### Start: Modify this block to tailor your analysis ###

# Workflow
# Values should be boolean (0 or 1).
run_busco=1
run_get_busco_summary=1
# Polyploidization history estimation
run_busco_getfasta=1
run_busco_mafft=1
run_busco_trimal=1
run_busco_iqtree_dna=1
run_busco_iqtree_pep=1
run_busco_notung_root_dna=1
run_busco_notung_root_pep=1
run_busco_root_dna=1
run_busco_root_pep=1
run_busco_grampa_dna=1
run_busco_grampa_pep=1
run_orthogroup_grampa=1 # Requires outputs of run_tree_root=1 of gg_gene_evolution with mode_orthogroup=1
# Gene family evolution
run_cafe=0
run_go_enrichment=0

# Other parameters
genetic_code=1 # Integer. See here https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes
busco_lineage="embryophyta_odb12" # See here for available datasets: https://busco-data.ezlab.org/v5/data/lineages/
nucleotide_model='GTR+R4' # IQ-TREE's model
protein_model='LG+R4' # IQ-TREE's model
grampa_h1='' # '' for exploratory analysis. If an internal node should be specified, run grampa once with grampa_h1='' to obtain internal node IDs.
min_gene_orthogroup_grampa=5 # Lower threshold for the orthogroup size to be included in the orthogroup-based grampa analysis.
max_gene_orthogroup_grampa=50 # Upper threshold for the orthogroup size to be included in the orthogroup-based grampa analysis.
notung_jar="/usr/local/bin/Notung.jar"
max_size_differential_cafe=9999999 # Maximum allowed size differential among species for each orthogroup to be included in the CAFE5 analysis.
n_gamma_cats_cafe=4 # Number of gamma rate categories to use in CAFE5.
target_branch_go="<25>" # Target branch for GO enrichment analysis (species name or branch ID; see workspace/genome_evolution/cafe/summary_plot/branch_id.pdf)
change_direction_go="increase" # Direction of gene change to focus on ("increase" or "decrease") in GO enrichment analysis.
species_go_annotation="Arabidopsis_thaliana" # Species used for GO annotation
go_category="BP,MF,CC" # GO categories to analyze (BP: Biological Process, MF: Molecular Function, CC: Cellular Component)

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

migrate_legacy_fasta_outputs() {
  local target_dir=$1
  local legacy_pattern=$2
  local legacy_files=()

  if [[ -z "${target_dir}" || -z "${legacy_pattern}" || ! -d "${target_dir}" ]]; then
    return 0
  fi

  mapfile -t legacy_files < <(find "${target_dir}" -maxdepth 1 -type f -name "${legacy_pattern}" | sort)
  for legacy_file in "${legacy_files[@]}"; do
    local migrated_file="${legacy_file%.fasta}.fa.gz"
    if [[ ! -s "${migrated_file}" ]]; then
      seqkit seq --threads 1 "${legacy_file}" --out-file "${migrated_file}"
    fi
    rm -f "${legacy_file}"
  done
}
dir_sp_cds="${dir_pg_input}/species_cds"
dir_og="${dir_pg_output}/orthogroup"
dir_og_rooted_tree="${dir_og}/rooted_tree"
dir_species_tree="${dir_pg_output}/species_tree"
dir_genome_evolution="${dir_pg_output}/genome_evolution"
dir_species_busco_full="${dir_pg_output}/species_cds_busco_full"
dir_species_busco_short="${dir_pg_output}/species_cds_busco_short"
dir_busco_summary_table="${dir_genome_evolution}/busco_summary_table"
dir_busco_fasta="${dir_genome_evolution}/busco_cds_fasta"
dir_busco_mafft="${dir_genome_evolution}/busco_mafft"
dir_busco_trimal="${dir_genome_evolution}/busco_trimal"
dir_busco_iqtree_dna="${dir_genome_evolution}/busco_iqtree_dna"
dir_busco_notung_dna="${dir_genome_evolution}/busco_notung_dna"
dir_busco_rooted_nwk_dna="${dir_genome_evolution}/busco_rooted_nwk_dna"
dir_busco_rooted_txt_dna="${dir_genome_evolution}/busco_rooted_txt_dna"
dir_busco_grampa_dna="${dir_genome_evolution}/grampa_busco_dna"
dir_busco_iqtree_pep="${dir_genome_evolution}/busco_iqtree_pep"
dir_busco_notung_pep="${dir_genome_evolution}/busco_notung_pep"
dir_busco_rooted_nwk_pep="${dir_genome_evolution}/busco_rooted_nwk_pep"
dir_busco_rooted_txt_pep="${dir_genome_evolution}/busco_rooted_txt_pep"
dir_busco_grampa_pep="${dir_genome_evolution}/grampa_busco_pep"
dir_orthogroup_grampa="${dir_genome_evolution}/grampa_orthogroup"
dir_cafe="${dir_genome_evolution}/cafe"
dir_cafe_orthogroup_selection="${dir_cafe}/orthogroup_selection"
dir_cafe_output="${dir_cafe}/cafe_output"
dir_cafe_summary_plot="${dir_cafe}/summary_plot"
dir_cafe_each_family_plot="${dir_cafe}/each_family_plot"
dir_species_cds_annotation="${dir_pg_output}/species_cds_annotation"
dir_go_enrichment="${dir_cafe}/go_enrichment"
dir_tmp="${dir_genome_evolution}/tmp"

check_species_cds ${dir_pg}
check_if_species_files_unique ${dir_sp_cds}
memory_notung=${MEM_PER_SLOT}

file_genecount="${dir_pg_output}/orthofinder/Orthogroups/Orthogroups.GeneCount.selected.tsv"
file_dated_species_tree="${dir_species_tree}/dated_species_tree.nwk"
file_busco_summary_table="${dir_busco_summary_table}/busco_summary.tsv"
file_busco_grampa_dna="${dir_busco_grampa_dna}/grampa_summary.tsv"
file_busco_grampa_pep="${dir_busco_grampa_pep}/grampa_summary.tsv"
file_orthogroup_grampa="${dir_orthogroup_grampa}/grampa_summary.tsv"
file_gene_id="${dir_pg_output}/orthofinder/Orthogroups/Orthogroups.selected.tsv"
file_go_annotation="${dir_species_cds_annotation}/${species_go_annotation}_annotation.tsv"

if [[ ! -s "${file_genecount}" && -s "${dir_pg_output}/orthofinder/Orthogroups/hog2og/Orthogroups.GeneCount.selected.tsv" ]]; then
  file_genecount="${dir_pg_output}/orthofinder/Orthogroups/hog2og/Orthogroups.GeneCount.selected.tsv"
fi
if [[ ! -s "${file_gene_id}" && -s "${dir_pg_output}/orthofinder/Orthogroups/hog2og/Orthogroups.selected.tsv" ]]; then
  file_gene_id="${dir_pg_output}/orthofinder/Orthogroups/hog2og/Orthogroups.selected.tsv"
fi

enable_all_run_flags_for_debug_mode

if [[ ${run_busco_notung_root_dna} -eq 1 || ${run_busco_notung_root_pep} -eq 1 || ${run_busco_root_dna} -eq 1 || ${run_busco_root_pep} -eq 1 || ${run_busco_grampa_dna} -eq 1 || ${run_busco_grampa_pep} -eq 1 ]]; then
  if [[ ! -s "${notung_jar}" ]]; then
    echo "Notung jar was not found: ${notung_jar}"
    echo "Disabling NOTUNG-dependent BUSCO tasks: run_busco_notung_root_dna, run_busco_notung_root_pep, run_busco_root_dna, run_busco_root_pep, run_busco_grampa_dna, run_busco_grampa_pep"
    run_busco_notung_root_dna=0
    run_busco_notung_root_pep=0
    run_busco_root_dna=0
    run_busco_root_pep=0
    run_busco_grampa_dna=0
    run_busco_grampa_pep=0
  fi
fi

if [[ ${run_orthogroup_grampa} -eq 1 ]]; then
  if [[ ! -s "${file_genecount}" ]]; then
    echo "Disabling run_orthogroup_grampa because required file is missing: ${file_genecount}"
    run_orthogroup_grampa=0
  elif [[ ! -d "${dir_og_rooted_tree}" ]]; then
    echo "Disabling run_orthogroup_grampa because required directory is missing: ${dir_og_rooted_tree}"
    run_orthogroup_grampa=0
  fi
fi

ensure_dir "${dir_tmp}"
cd ${dir_tmp}

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
if [[ ${run_busco_getfasta} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_busco_fasta}"
  migrate_legacy_fasta_outputs "${dir_busco_fasta}" "*.busco.cds.fasta"
  busco_rows=()
  mapfile -t busco_rows < <(tail -n +2 "${file_busco_summary_table}")
  num_busco_ids=${#busco_rows[@]}

  generate_single_copy_fasta() {
    local busco_idx=$1
    local busco_row="${busco_rows[${busco_idx}]}"
    local busco_id
    local outfile2
    local cols=()
    local genes=()
    local genes2=()
    IFS=$'\t' read -r -a cols <<< "${busco_row}"
    busco_id="${cols[0]:-}"
    if [[ -z "${busco_id}" ]]; then
      return 0
    fi
    if [[ ${#cols[@]} -gt 3 ]]; then
      genes=( "${cols[@]:3}" )
    fi
    outfile2="${dir_busco_fasta}/${busco_id}.busco.cds.fa.gz"
    if [[ -s ${outfile2} ]]; then
      return 0
    fi
    echo "busco_id: ${busco_id}"
    if [[ ! -s ${outfile2} ]]; then
      mapfile -t genes2 < <(gg_busco_gene_tokens "split_duplicated" "${genes[@]}")
      if [[ ${#genes2[@]} -eq 0 ]]; then
        echo "Skipping. ${busco_id} has no genes in the selected species."
        return 0
      fi
      gg_seqkit_grep_by_patterns_from_infile_list 1 "species_cds_fasta_list.txt" "${genes2[@]}" \
      | gg_prepare_cds_fasta_stream 1 \
      | seqkit seq --threads 1 --out-file "${outfile2}"
      if [[ ! -s ${outfile2} ]]; then
        echo "File is empty. Removing: ${outfile2}"
        rm ${outfile2}
      fi
    fi
  }

  gg_find_fasta_files "${dir_sp_cds}" 1 > species_cds_fasta_list.txt
  for (( busco_idx=0; busco_idx<num_busco_ids; busco_idx++ )); do
    wait_until_jobn_le ${NSLOTS}
    generate_single_copy_fasta "${busco_idx}" &
  done
  wait
  find . -maxdepth 1 -type f -name 'tmp.*' -delete
else
  gg_step_skip "${task}"
fi

task="In-frame mafft alignment of duplicate-containing BUSCO genes"
if [[ ${run_busco_mafft} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_busco_mafft}"
  migrate_legacy_fasta_outputs "${dir_busco_fasta}" "*.busco.cds.fasta"
  migrate_legacy_fasta_outputs "${dir_busco_mafft}" "*.busco.cds.aln.fasta"

  run_mafft() {
    infile=$1
    infile_base=${infile%%.*}
    outfile=${dir_busco_mafft}/${infile_base}.busco.cds.aln.fa.gz
    if [[ -s ${outfile} ]]; then
      return 0
    fi
    echo "$(date): start mafft: ${infile_base}"

    seqkit seq --threads 1 "${dir_busco_fasta}/${infile}" --out-file "tmp.${infile_base}.input.cds.fasta"
    cdskit mask \
    --seqfile "tmp.${infile_base}.input.cds.fasta" \
    --outfile tmp.${infile_base}.cds.fasta

    num_seq=$(gg_count_fasta_records "tmp.${infile_base}.cds.fasta")
    if [[ ${num_seq} -lt 2 ]]; then
      echo "Skipped MAFFT/backalign because fewer than 2 sequences were found: ${infile}"
      seqkit seq --threads 1 "tmp.${infile_base}.cds.fasta" --out-file "tmp.${infile_base}.cds.out.fa.gz"
      mv_out "tmp.${infile_base}.cds.out.fa.gz" "${outfile}"
      rm tmp.${infile_base}*
      return 0
    fi

    seqkit translate \
    --allow-unknown-codon \
    --transl-table ${genetic_code} \
    --threads 1 \
    tmp.${infile_base}.cds.fasta \
    > tmp.${infile_base}.pep.fasta

    mafft \
    --auto \
    --amino \
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
  mapfile -t input_alignment_files < <(gg_find_file_basenames "${dir_busco_fasta}" "*.busco.cds.fa.gz")
  echo "Number of input alignments: ${#input_alignment_files[@]}"
  for input_alignment_file in "${input_alignment_files[@]}"; do
    wait_until_jobn_le ${NSLOTS}
    run_mafft ${input_alignment_file} &
  done
  wait
else
  gg_step_skip "${task}"
fi

task="TrimAl of duplicate-containing BUSCO genes"
if [[ ${run_busco_trimal} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_busco_trimal}"
  migrate_legacy_fasta_outputs "${dir_busco_mafft}" "*.busco.cds.aln.fasta"
  migrate_legacy_fasta_outputs "${dir_busco_trimal}" "*.busco.trimal.fasta"

  run_trimal() {
    infile=$1
    infile_base=${infile%%.*}
    outfile="${dir_busco_trimal}/${infile_base}.busco.trimal.fa.gz"
    if [[ -s ${outfile} ]]; then
      return 0
    fi
    seqkit seq --remove-gaps --threads 1 ${dir_busco_mafft}/${infile} > tmp.${infile_base}.degap.fasta
    seqkit translate --transl-table ${genetic_code} --threads 1 ${dir_busco_mafft}/${infile} > tmp.${infile_base}.pep.fasta

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
  }

  input_alignment_files=()
  mapfile -t input_alignment_files < <(gg_find_file_basenames "${dir_busco_mafft}" "*.busco.cds.aln.fa.gz")
  echo "Number of input alignments: ${#input_alignment_files[@]}"
  for input_alignment_file in "${input_alignment_files[@]}"; do
    wait_until_jobn_le ${NSLOTS}
    run_trimal ${input_alignment_file} &
  done
  wait
else
  gg_step_skip "${task}"
fi

task="IQ-TREE for duplicate-containing BUSCO DNA trees"
if [[ ${run_busco_iqtree_dna} -eq 1 ]]; then
    gg_step_start "${task}"
  ensure_dir "${dir_busco_iqtree_dna}"
  migrate_legacy_fasta_outputs "${dir_busco_trimal}" "*.busco.trimal.fasta"

  busco_iqtree_dna() {
    infile=$1
    indir=$2
    outdir=$3
    infile_base=${infile%%.*}
    outfile="${outdir}/${infile_base}.busco.nwk"
    if [[ -s ${outfile} ]]; then
      return 0
    fi
    num_seq=$(gg_count_fasta_records "${dir_busco_trimal}/${infile}")
    if [[ ${num_seq} -lt 3 ]]; then
      echo "Skipped. At least 3 sequences are necessary for IQ-TREE: ${infile}"
      return 0
    fi

    seqkit seq --threads 1 "${indir}/${infile}" --out-file "./tmp.${infile_base}.input.fasta"

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
  mapfile -t input_alignment_files < <(gg_find_file_basenames "${dir_busco_trimal}" "*.busco.trimal.fa.gz")
  echo "Number of input alignments: ${#input_alignment_files[@]}"
  for input_alignment_file in "${input_alignment_files[@]}"; do
    wait_until_jobn_le ${NSLOTS}
    busco_iqtree_dna ${input_alignment_file} ${dir_busco_trimal} ${dir_busco_iqtree_dna} &
  done
  wait
else
  gg_step_skip "${task}"
fi

task="IQ-TREE for duplicate-containing BUSCO protein trees"
if [[ ${run_busco_iqtree_pep} -eq 1 ]]; then
    gg_step_start "${task}"
  ensure_dir "${dir_busco_iqtree_pep}"
  migrate_legacy_fasta_outputs "${dir_busco_trimal}" "*.busco.trimal.fasta"

  busco_iqtree_pep() {
    infile=$1
    indir=$2
    outdir=$3
    infile_base=${infile%%.*}
    outfile="${outdir}/${infile_base}.busco.nwk"
    if [[ -s ${outfile} ]]; then
      return 0
    fi
    num_seq=$(gg_count_fasta_records "${dir_busco_trimal}/${infile}")
    if [[ ${num_seq} -lt 3 ]]; then
      echo "Skipped. At least 3 sequences are necessary for IQ-TREE: ${infile}"
      return 0
    fi

    seqkit translate --transl-table ${genetic_code} --threads 1 "${indir}/${infile}" > "./tmp.${infile_base}.input.fasta"

    iqtree \
    -s "./tmp.${infile_base}.input.fasta" \
    -m ${protein_model} \
    -st AA \
    -T 1 \
    --prefix tmp.${infile_base} \
    --seed 12345 \
    --redo

    mv_out tmp.${infile_base}.treefile ${outfile}
    rm tmp.${infile_base}.*
  }

  input_alignment_files=()
  mapfile -t input_alignment_files < <(gg_find_file_basenames "${dir_busco_trimal}" "*.busco.trimal.fa.gz")
  echo "Number of input alignments: ${#input_alignment_files[@]}"
  for input_alignment_file in "${input_alignment_files[@]}"; do
    wait_until_jobn_le ${NSLOTS}
    busco_iqtree_pep ${input_alignment_file} ${dir_busco_trimal} ${dir_busco_iqtree_pep} &
  done
  wait
else
  gg_step_skip "${task}"
fi

busco_notung () {
  infile=$1
  indir=$2
  outdir=$3
  busco_id=${infile%%.*}
  outfile="${outdir}/${busco_id}.busco.notung.root.zip"
  if [[ -s "${outfile}" ]]; then
    return 0
  fi
  if [[ -e "./${busco_id}.notung.root" ]]; then
    rm -r "./${busco_id}.notung.root"
  fi
  java -jar -Xmx${memory_notung}g ${notung_jar} \
  -s "${file_dated_species_tree}" \
  -g "${indir}/${infile}" \
  --root \
  --infertransfers "false" \
  --treeoutput newick \
  --log \
  --treestats \
  --events \
  --parsable \
  --speciestag prefix \
  --allopt \
  --maxtrees 1000 \
  --nolosses \
  --outputdir "./${busco_id}.notung.root"
  if [[ -e "${busco_id}.notung.root/${busco_id}.busco.nwk.rooting.0" ]]; then
    zip -rq "${busco_id}.busco.notung.root.zip" "${busco_id}.notung.root"
    mv_out "${busco_id}.busco.notung.root.zip" "${outfile}"
    rm -r "./${busco_id}.notung.root"
  fi
}

task="NOTUNG rooting of duplicate-containing BUSCO DNA trees"
if [[ ${run_busco_notung_root_dna} -eq 1 ]]; then
  gg_step_start "${task}"

  infiles=()
  mapfile -t infiles < <(gg_find_file_basenames "${dir_busco_iqtree_dna}")
  for infile in "${infiles[@]}"; do
    wait_until_jobn_le $((${NSLOTS}/2))
    busco_notung ${infile} ${dir_busco_iqtree_dna} ${dir_busco_notung_dna} &
  done
  wait
else
  gg_step_skip "${task}"
fi

task="NOTUNG rooting of duplicate-containing BUSCO protein trees"
if [[ ${run_busco_notung_root_pep} -eq 1 ]]; then
  gg_step_start "${task}"

  infiles=()
  mapfile -t infiles < <(gg_find_file_basenames "${dir_busco_iqtree_pep}")
  for infile in "${infiles[@]}"; do
    wait_until_jobn_le $((${NSLOTS}/2))
    busco_notung ${infile} ${dir_busco_iqtree_pep} ${dir_busco_notung_pep} &
  done
  wait
else
  gg_step_skip "${task}"
fi

busco_species_tree_assisted_gene_tree_rooting () {
  infile=$1
  indir=$2
  intreedir=$3
  outdir_txt=$4
  outdir_nwk=$5
  busco_id=${infile%%.*}
  intree="${intreedir}/${busco_id}.busco.nwk"
  outfile_txt="${outdir_txt}/${busco_id}.busco.root.txt"
  outfile_nwk="${outdir_nwk}/${busco_id}.busco.root.nwk"
  if [[ -s ${outfile_txt} && -s ${outfile_nwk} ]]; then
    return 0
  fi
  echo "Start NOTUNG root: ${busco_id}"
  if [[ -e ./${busco_id}.notung.root ]]; then
    rm -r ./${busco_id}.notung.root
  fi
  cp_out ${indir}/${infile} .
  unzip -q ${infile}

  Rscript ${dir_myscript}/species_tree_guided_gene_tree_rooting.r \
  "--notung_root_zip=./${infile}" \
  "--in_tree=${intree}" \
  "--out_tree=${busco_id}.root.nwk" \
  "--ncpu=${NSLOTS}" \
  2>&1 | tee ${busco_id}.root.txt

  if [[ -s ${busco_id}.root.nwk ]]; then
    mv_out ${busco_id}.root.txt ${outfile_txt}
    mv_out ${busco_id}.root.nwk ${outfile_nwk}
    rm ${infile}
    rm -r ${busco_id}.notung.root
  fi
}

task="Species-tree-guided gene tree rooting of duplicate-containing BUSCO DNA trees"
if [[ ${run_busco_root_dna} -eq 1 ]]; then
    gg_step_start "${task}"

  infiles=()
  mapfile -t infiles < <(gg_find_file_basenames "${dir_busco_notung_dna}")
  for infile in "${infiles[@]}"; do
    wait_until_jobn_le ${NSLOTS}
    busco_species_tree_assisted_gene_tree_rooting ${infile} ${dir_busco_notung_dna} ${dir_busco_iqtree_dna} ${dir_busco_rooted_txt_dna} ${dir_busco_rooted_nwk_dna} &
  done
  wait
else
  gg_step_skip "${task}"
fi

task="Species-tree-guided gene tree rooting of duplicate-containing BUSCO protein trees"
if [[ ${run_busco_root_pep} -eq 1 ]]; then
    gg_step_start "${task}"

  infiles=()
  mapfile -t infiles < <(gg_find_file_basenames "${dir_busco_notung_pep}")
  for infile in "${infiles[@]}"; do
    wait_until_jobn_le ${NSLOTS}
    busco_species_tree_assisted_gene_tree_rooting ${infile} ${dir_busco_notung_pep} ${dir_busco_iqtree_pep} ${dir_busco_rooted_txt_pep} ${dir_busco_rooted_nwk_pep} &
  done
  wait
else
  gg_step_skip "${task}"
fi

busco_grampa () {
  indir=$1
  outdir=$2
  outfile=$3
  ensure_dir "${outdir}"
  local nwk_files=()
  shopt -s nullglob
  nwk_files=( "${indir}"/*.nwk )
  shopt -u nullglob
  if [[ ${#nwk_files[@]} -eq 0 ]]; then
    echo "Skipping Grampa because no rooted gene trees were found in: ${indir}"
    return 0
  fi
  if [[ -e "./grampa_out" ]]; then
    rm -r "./grampa_out"
  fi
  nwkit drop --infile ${file_dated_species_tree} --target 'intnode' --name 'yes' \
  | sed -e "s/_/-/g" \
  > "grampa_input_species_tree.nwk"

  cat "${nwk_files[@]}" \
  | sed -E "s/([(,])([^_)(,:]+)_([^_)(,:]+)_([^)(,:]+)([)(,:])/\1\4|||\2\-\3\5/g" \
  | sed -e "s/_/-/g" -e "s/|||/_/g" \
  > "grampa_input_gene_trees.nwk"

  printf "%s\n" "${nwk_files[@]##*/}" > "busco_genetree_filenames.txt"

  if [[ ${grampa_h1} == "" ]]; then
    h1_param=""
  else
    local grampa_h1_normalized=${grampa_h1//_/-}
    grampa_h1_normalized=${grampa_h1_normalized//[[:space:]]/-}
    h1_param="-h1 ${grampa_h1_normalized}"
  fi

  grampa.py \
  -s "grampa_input_species_tree.nwk" \
  -g "grampa_input_gene_trees.nwk" \
  -o "./grampa_out" \
  -p ${NSLOTS} \
  -v -1 \
  --maps \
  ${h1_param}

  grampa_filtered_file=""
  if [[ -s "./grampa_out/grampa_trees_filtered.txt" ]]; then
    grampa_filtered_file="./grampa_out/grampa_trees_filtered.txt"
  elif [[ -s "./grampa_out/grampa-trees-filtered.txt" ]]; then
    grampa_filtered_file="./grampa_out/grampa-trees-filtered.txt"
  fi

  if [[ -n "${grampa_filtered_file}" ]]; then
    # Filtered gene trees are not analyzed by Grampa, but it does not disturb gene tree IDs.
    # For example, if the 136th gene tree is filtered out, it is replaced with a placeholder text in grampa_trees_filtered.txt.
    # And GT-136 does not appear in the Grampa outputs. Still, the 137th gene tree is labeled correctly as GT-137.
    sed -e "s/$/;/" "${grampa_filtered_file}" > "grampa_input_gene_trees_filtered.nwk"
  fi

  grampa_det_file=""
  grampa_out_file=""
  grampa_checknums_file=""
  for candidate in "./grampa_out/grampa_det.txt" "./grampa_out/grampa-detailed.txt"; do
    if [[ -s "${candidate}" ]]; then
      grampa_det_file="${candidate}"
      break
    fi
  done
  for candidate in "./grampa_out/grampa_out.txt" "./grampa_out/grampa-scores.txt"; do
    if [[ -s "${candidate}" ]]; then
      grampa_out_file="${candidate}"
      break
    fi
  done
  for candidate in "./grampa_out/grampa_checknums.txt" "./grampa_out/grampa-checknums.txt"; do
    if [[ -s "${candidate}" ]]; then
      grampa_checknums_file="${candidate}"
      break
    fi
  done

  if [[ -z "${grampa_det_file}" || -z "${grampa_out_file}" || -z "${grampa_checknums_file}" ]]; then
    echo "Grampa output files are missing. Ending Grampa."
    return 0
  fi

  python ${dir_myscript}/parse_grampa.py \
  --grampa_det "${grampa_det_file}" \
  --grampa_out "${grampa_out_file}" \
  --gene_trees "./grampa_input_gene_trees.nwk" \
  --species_tree "./grampa_input_species_tree.nwk" \
  --ncpu ${NSLOTS} \
  --sorted_gene_tree_file_names "./busco_genetree_filenames.txt"

  if [[ -s "${grampa_checknums_file}" && -s "${grampa_det_file}" && -s "${grampa_out_file}" && -s "grampa_summary.tsv" ]]; then
    if ! awk -F'\t' '/^The MUL-tree with the minimum parsimony score/ {print $NF; found=1} END{exit(found?0:1)}' "${grampa_out_file}" > "${outdir}/best_mul_tree.nwk"; then
      awk -F'\t' '
        /^MT-/ {
          score = $NF + 0
          if (!seen || score < best_score) {
            best_score = score
            best_tree = $(NF-1)
            seen = 1
          }
        }
        END {
          if (seen) {
            print best_tree
          }
        }
      ' "${grampa_out_file}" > "${outdir}/best_mul_tree.nwk"
    fi
    cp_out "${grampa_checknums_file}" "${outdir}/grampa_checknums.txt"
    cp_out "${grampa_det_file}" "${outdir}/grampa_det.txt"
    cp_out "${grampa_out_file}" "${outdir}/grampa_out.txt"
    if [[ -s "./grampa_input_gene_trees_filtered.nwk" ]]; then
      mv_out "./grampa_input_gene_trees_filtered.nwk" "${outdir}"
    fi
    mv_out "./grampa_input_species_tree.nwk" "${outdir}"
    mv_out "./grampa_input_gene_trees.nwk" "${outdir}"
    mv_out "./grampa_summary.tsv" "${outfile}"
    rm -r "./grampa_out"
  else
    echo "Grampa output files are missing. Ending Grampa."
  fi
}

task="BUSCO-based Grampa analysis for the polyploidization history with BUSCO DNA trees"
if [[ ! -s ${file_busco_grampa_dna} && ${run_busco_grampa_dna} -eq 1 ]]; then
  gg_step_start "${task}"

  busco_grampa ${dir_busco_rooted_nwk_dna} ${dir_busco_grampa_dna} ${file_busco_grampa_dna}
else
  gg_step_skip "${task}"
fi

task="BUSCO-based Grampa analysis for the polyploidization history with BUSCO protein trees"
if [[ ! -s ${file_busco_grampa_pep} && ${run_busco_grampa_pep} -eq 1 ]]; then
  gg_step_start "${task}"

  busco_grampa ${dir_busco_rooted_nwk_pep} ${dir_busco_grampa_pep} ${file_busco_grampa_pep}
else
  gg_step_skip "${task}"
fi

task="Orthogroup-based Grampa analysis for polyploidization history"
if [[ ! -s ${file_orthogroup_grampa} && ${run_orthogroup_grampa} -eq 1 ]]; then
  gg_step_start "${task}"

  og_ids=()
  mapfile -t og_ids < <(python -c 'import pandas,sys; d=pandas.read_csv(sys.argv[1],sep="\t",header=0); ids=d.loc[(d["Total"]>=int(sys.argv[2]))&(d["Total"]<=int(sys.argv[3])),"Orthogroup"].astype(str).tolist(); print("\n".join(ids))' ${file_genecount} ${min_gene_orthogroup_grampa} ${max_gene_orthogroup_grampa})
  file_names=()
  mapfile -t file_names < <(gg_find_file_basenames "${dir_og_rooted_tree}")
  echo "Number of files in ${dir_og_rooted_tree}: ${#file_names[@]}"
  echo "Number of selected orthogroups with ${min_gene_orthogroup_grampa}<=gene number<=${max_gene_orthogroup_grampa}: ${#og_ids[@]}"
  if [[ -e ./tmp.orthogroup_grampa_indir ]]; then
    rm -r ./tmp.orthogroup_grampa_indir
  fi
  mkdir ./tmp.orthogroup_grampa_indir
  for file_name in "${file_names[@]}"; do
    for og_id in "${og_ids[@]}"; do
      if [[ ${file_name} == ${og_id}* ]]; then
        cp_out ${dir_og_rooted_tree}/${file_name} ./tmp.orthogroup_grampa_indir
        mapfile -t og_ids < <(printf "%s\n" "${og_ids[@]}" | grep -v -Fx "${og_id}" || true)
        break
      fi
    done
  done

  busco_grampa ./tmp.orthogroup_grampa_indir ${dir_orthogroup_grampa} ${file_orthogroup_grampa}
else
  gg_step_skip "${task}"
fi

task='CAFE analysis'
disable_if_no_input_file "run_cafe" ${file_genecount} ${file_dated_species_tree}
if [[ (! -s "${dir_cafe_summary_plot}/summary_all.pdf" || ! -s "${dir_cafe_summary_plot}/summary_significant.pdf") && ${run_cafe} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_cafe}"

  if [[ ! -s "${dir_cafe_output}/Gamma_asr.tre" || ! -s "${dir_cafe_output}/Gamma_count.tab" || ! -s "${dir_cafe_output}/Gamma_change.tab" ]]; then
    if [[ -d ${dir_cafe_orthogroup_selection} ]]; then
      rm -r ${dir_cafe_orthogroup_selection}
    fi
    python ${dir_myscript}/cafe_orthogroup_selection.py \
    --genecount ${file_genecount} \
    --dated_species_tree ${file_dated_species_tree} \
    --output_dir ${dir_cafe_orthogroup_selection} \
    --max_size_differential ${max_size_differential_cafe}

    if [[ -d ${dir_cafe_output} ]]; then
      rm -r ${dir_cafe_output}
    fi
    cafe5 \
    --infile ${dir_cafe_orthogroup_selection}/cafe_input.tsv \
    --tree ${file_dated_species_tree} \
    --n_gamma_cats ${n_gamma_cats_cafe} \
    --pvalue 0.05 \
    --cores ${NSLOTS} \
    --output_prefix ${dir_cafe_output}
  else
    echo "CAFE output files already exist. Skipping CAFE run."
  fi

  if [[ -s "${dir_cafe_output}/Gamma_asr.tre" && -s "${dir_cafe_output}/Gamma_count.tab" && -s "${dir_cafe_output}/Gamma_change.tab" ]]; then
    if [[ -d ${dir_cafe_each_family_plot} ]]; then
      rm -r ${dir_cafe_each_family_plot}
    fi
    Rscript ${dir_myscript}/cafe_plot_each_family.r \
    ${dir_cafe_output}/Gamma_asr.tre \
    ${dir_cafe_output}/Gamma_count.tab \
    ${dir_cafe_output}/Gamma_change.tab \
    ${dir_cafe_each_family_plot} \
    ${NSLOTS}
    status=$?
    if [[ ${status} -ne 0 ]]; then
      echo "Error in Rscript cafe_plot_each_family.r. Exiting."
      exit 1
    fi

    if [[ -d ${dir_cafe_summary_plot} ]]; then
      rm -r ${dir_cafe_summary_plot}
    fi
    Rscript ${dir_myscript}/cafe_plot_summary.r \
    ${dir_cafe_output}/Gamma_asr.tre \
    ${dir_cafe_output}/Gamma_change.tab \
    ${dir_cafe_summary_plot}

    Rscript ${dir_myscript}/cafe_plot_branch_id.r \
    ${dir_cafe_output}/Gamma_asr.tre \
    ${dir_cafe_summary_plot}
  else
    echo "CAFE did not finish successfully. Exiting."
    exit 1
  fi

  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

task="GO enrichment analysis"
disable_if_no_input_file "run_go_enrichment" "${dir_cafe_output}/Gamma_change.tab" "${dir_cafe_output}/Gamma_branch_probabilities.tab" ${file_gene_id} ${file_go_annotation}
if [[(! -s "${dir_go_enrichment}/enrichment_significant_${change_direction_go}_${target_branch_go}_significant_go.tsv") && ${run_go_enrichment} -eq 1 ]]; then
  gg_step_start "${task}"
  ensure_dir "${dir_go_enrichment}"
  Rscript ${dir_myscript}/cafe_go_enrichment.r \
  "${dir_cafe_output}/Gamma_change.tab" \
  "${dir_cafe_output}/Gamma_branch_probabilities.tab" \
  "${file_gene_id}" \
  "${file_go_annotation}" \
  "${dir_go_enrichment}" \
  "${target_branch_go}" \
  "${change_direction_go}" \
  "${go_category}"
  status=$?
  if [[ ${status} -ne 0 ]]; then
    echo "Error in Rscript cafe_go_enrichment.r. Exiting."
    exit 1
  fi
  echo "$(date): End: ${task}"
else
  gg_step_skip "${task}"
fi

if [[ ${delete_tmp_dir} -eq 1 ]]; then
  if [[ -d "${dir_tmp}" ]]; then
    echo "Removing tmp directory: ${dir_tmp}"
    rm -rf -- "${dir_tmp}"
  fi
fi

echo $(date): end
__GG_INLINE_STAGE_GENOMEEVOLUTION_20260221__

echo "$(date): Completed unified genome evolution pipeline"
