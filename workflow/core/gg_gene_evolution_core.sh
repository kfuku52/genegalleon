#!/usr/bin/env bash
set -euo pipefail

# Load shared defaults when available.
gg_core_self="${BASH_SOURCE[0]:-/script/core/gg_gene_evolution_core.sh}"
gg_common_params_file=""
for gg_common_params_candidate in \
  "$(cd "$(dirname "${gg_core_self}")" && pwd)/gg_common_params.sh" \
  "$(cd "$(dirname "${gg_core_self}")" && pwd)/../gg_common_params.sh" \
  "/script/gg_common_params.sh"
do
  if [[ -s "${gg_common_params_candidate}" ]]; then
    gg_common_params_file="${gg_common_params_candidate}"
    break
  fi
done
if [[ -n "${gg_common_params_file}" ]]; then
  # shellcheck disable=SC1090
  source "${gg_common_params_file}"
fi
unset gg_common_params_file gg_common_params_candidate gg_core_self

### Start: Job-supplied configuration ###
# Configuration variables are provided by gg_gene_evolution_entrypoint.sh.
genetic_code="${genetic_code:-${GG_COMMON_GENETIC_CODE:-1}}"
treevis_clade_ortholog_prefix="${treevis_clade_ortholog_prefix:-${GG_COMMON_TREEVIS_CLADE_ORTHOLOG_PREFIX:-Arabidopsis_thaliana_}}"

# Substitution model in CSUBST and mapdNdS
if [[ ${genetic_code} -eq 1 ]]; then
  codon_model="ECMK07+F+R4" # Standard genetic code
else
  codon_model="GY+F+R4" # Non-standard genetic code
fi
### End: Job-supplied configuration ###

### ----------------------------------------------------------------------- ###

### Modify below if you need to add a new analysis or need to fix some bugs ###

dir_pg="/workspace"
dir_script="/script/support"
source "${dir_script}/gg_util.sh" # Load utility functions
gg_source_home_bashrc
gg_prepare_cmd_runtime "${dir_pg}" "base" 1 1
delete_preexisting_tmp_dir=${delete_preexisting_tmp_dir:-1}
delete_tmp_dir=${delete_tmp_dir:-1}

build_iqtree_mem_args() {
  IQTREE_MEM_ARGS=()
  if [[ -n "${MEM_PER_HOST:-}" ]]; then
    IQTREE_MEM_ARGS=(--mem "${MEM_PER_HOST}G")
  fi
}

binarize_species_trait() {
  local file_in="$1"
  local file_out="$2"
  python - "${file_in}" "${file_out}" <<'PY'
import sys
import numpy
import pandas

file_in, file_out = sys.argv[1], sys.argv[2]
df = pandas.read_csv(file_in, sep="\t", header=0, dtype=str)
if df.shape[1] < 2:
    raise ValueError(f"Trait file must have at least 2 columns: {file_in}")

species_col = df.columns[0]
out = pandas.DataFrame()
out[species_col] = df[species_col].astype(str)

na_tokens = {"", "NA", "NaN", "nan", "."}
truthy = {"1", "true", "yes", "y", "t"}
falsey = {"0", "false", "no", "n", "f"}

for col in df.columns[1:]:
    raw = df[col].astype(str).str.strip()
    raw = raw.where(~raw.isin(na_tokens), numpy.nan)
    numeric = pandas.to_numeric(raw, errors="coerce")
    valid = numeric.dropna()

    if valid.size > 0:
        uniq = set(valid.astype(float).tolist())
        if uniq.issubset({0.0, 1.0}):
            binary = numeric.fillna(0).astype(float).clip(0, 1).round().astype(int)
        else:
            threshold = float(valid.median())
            binary = (numeric > threshold).astype(int)
            # Keep the split informative when many values are tied at median.
            if int(binary.sum()) in {0, int(binary.shape[0])}:
                threshold = float(valid.mean())
                binary = (numeric > threshold).astype(int)
            binary.loc[numeric.isna()] = 0
    else:
        lowered = raw.fillna("").str.lower()
        mapped = lowered.map(lambda x: 1 if x in truthy else (0 if x in falsey else 0))
        binary = mapped.astype(int)

    out[col] = binary

out.to_csv(file_out, sep="\t", index=False)
PY
}

normalize_mapnh_params_for_mapnh_v1() {
  local file_param="$1"
  local default_genetic_code="$2"

  if [[ ! -s "${file_param}" ]]; then
    echo "MapNH parameter file is missing: ${file_param}"
    return 1
  fi

  if grep -q '^input.sequence.file=' "${file_param}"; then
    return 0
  fi

  if ! grep -q '^input.data1=' "${file_param}"; then
    echo "Unknown MapNH parameter format (keeping original): ${file_param}"
    return 0
  fi

  local kappa omega alpha gamma_n genetic_code_local
  kappa=$(awk -F'kappa=|,omega=' '/kappa=/{ print $2; exit }' "${file_param}")
  omega=$(awk -F'omega=|\\)' '/omega=/{ print $2; exit }' "${file_param}")
  alpha=$(awk -F'alpha=|,Gamma.beta=' '/alpha=/{ print $2; exit }' "${file_param}")
  gamma_n=$(awk -F'Gamma\\(n=|,alpha=' '/Gamma\\(n=/{ print $2; exit }' "${file_param}")
  genetic_code_local=$(awk -F'=' '/^genetic_code=/{ print $2; exit }' "${file_param}")

  [[ -z "${kappa}" ]] && kappa="2.0"
  [[ -z "${omega}" ]] && omega="0.2"
  [[ -z "${alpha}" ]] && alpha="0.5"
  [[ -z "${gamma_n}" ]] && gamma_n="4"
  [[ -z "${genetic_code_local}" ]] && genetic_code_local="${default_genetic_code}"

  cat > "${file_param}" <<EOF
alphabet=Codon(letter=DNA)
genetic_code=${genetic_code_local}
input.sequence.file=\$(SEQ)
input.sequence.format=Fasta
input.sequence.remove_stop_codons=yes
input.tree.file=\$(TREE)
input.tree.format=Newick
model=YN98(kappa=${kappa},omega=${omega},initFreqs=observed)
rate_distribution=Gamma(n=${gamma_n},alpha=${alpha})
map.type=DnDs
output.counts=PerBranch(prefix=\$(OUT).)
output.tree_with_id.file=\$(OUT).with_id.nwk
EOF
}

assert_strictly_bifurcating_tree() {
  local file_tree="$1"
  local context="${2:-tree}"
  python - "${file_tree}" "${context}" <<'PY'
import sys
from Bio import Phylo

tree_file = sys.argv[1]
context = sys.argv[2]

try:
    tree = Phylo.read(tree_file, "newick")
except Exception as exc:  # pragma: no cover
    print(f"{context}: failed to parse Newick: {tree_file}")
    print(f"{exc}")
    sys.exit(1)

violations = []
for clade in tree.get_nonterminals(order="level"):
    n_children = len(clade.clades)
    if n_children != 2:
        label = clade.name if clade.name else "<unnamed>"
        violations.append((label, n_children))

if violations:
    print(f"{context}: non-bifurcating tree detected: {tree_file}")
    print("Resolve polytomies upstream before running downstream analyses.")
    for label, n_children in violations[:10]:
        print(f"  node={label}, children={n_children}")
    sys.exit(1)
PY
}

# Setting modes
if [[ ${gg_debug_mode:-0} -eq 1 ]]; then
  enable_all_run_flags_for_debug_mode "gg debug mode: All run_* variables are forced to set 1, except for too-time-consuming tasks."
  #run_orthogroup_extraction=0; echo "gg debug mode: run_orthogroup_extraction=${run_orthogroup_extraction}"
  run_codeml_two_ratio=0; echo "gg debug mode: run_codeml_two_ratio=${run_codeml_two_ratio}"
  run_hyphy_relax=0; echo "gg debug mode: run_hyphy_relax=${run_hyphy_relax}"
  run_hyphy_relax_reversed=0; echo "gg debug mode: run_hyphy_relax_reversed=${run_hyphy_relax_reversed}"
  csubst_cutoff_stat="OCNany2spe,0|omegaCany2spe,1"; echo "gg debug mode: csubst_cutoff_stat=${csubst_cutoff_stat}"
fi
query_blast_method=$(echo "${query_blast_method}" | tr '[:upper:]' '[:lower:]')
tree_rooting_method=$(echo "${tree_rooting_method}" | tr '[:upper:]' '[:lower:]')
if [[ "${tree_rooting_method}" != "notung" && "${tree_rooting_method}" != "midpoint" && "${tree_rooting_method}" != "mad" && "${tree_rooting_method}" != "md" ]]; then
  echo "Invalid tree_rooting_method: ${tree_rooting_method}"
  echo "tree_rooting_method must be one of notung, midpoint, mad, md. Exiting."
  exit 1
fi
if [[ ${mode_query2family} -eq 1 && ${run_query_blast} -eq 1 ]]; then
  if [[ "${query_blast_method}" != "tblastn" && "${query_blast_method}" != "diamond" ]]; then
		  echo "Invalid query_blast_method: ${query_blast_method}"
		  echo 'query_blast_method must be either "tblastn" or "diamond". Exiting.'
		  exit 1
	  fi
fi
case "${mode_orthogroup}:${mode_query2family}" in
		"1:0")
			dir_output_active="${dir_pg_output}/orthogroup"
			file_orthogroup_genecount_selected="${dir_pg_output}/orthofinder/Orthogroups_filtered/Orthogroups.GeneCount.selected.tsv"
			if [[ ! -s "${file_orthogroup_genecount_selected}" ]]; then
			  echo "Orthogroup gene-count table not found: ${file_orthogroup_genecount_selected}"
			  exit 1
			fi
			if [[ ! "${SGE_TASK_ID}" =~ ^[0-9]+$ ]] || [[ ${SGE_TASK_ID} -lt 1 ]]; then
			  echo "Invalid SGE_TASK_ID value (must be a positive integer): ${SGE_TASK_ID}"
			  exit 1
			fi
			num_orthogroups=$(awk 'END { print (NR > 0 ? NR - 1 : 0) }' "${file_orthogroup_genecount_selected}")
			if [[ ${num_orthogroups} -le 0 ]]; then
			  echo "No orthogroup rows were found in: ${file_orthogroup_genecount_selected}"
			  exit 1
			fi
			if [[ ${SGE_TASK_ID} -gt ${num_orthogroups} ]]; then
			  echo "SGE_TASK_ID=${SGE_TASK_ID} is out of range for ${num_orthogroups} orthogroups in: ${file_orthogroup_genecount_selected}"
			  exit 1
			fi
			og_id=$(awk -F'\t' -v row="${SGE_TASK_ID}" 'NR == (row + 1) { print $1; exit }' "${file_orthogroup_genecount_selected}")
			if [[ -z "${og_id}" ]]; then
			  echo "Failed to read OrthoGroup ID at row ${SGE_TASK_ID} from: ${file_orthogroup_genecount_selected}"
			  exit 1
			fi
			echo "OrthoGroup ID: ${og_id}"
			run_get_query_fasta=0
			run_query_blast=0
			run_orthogroup_extraction=0
			;;
		"0:1")
			dir_output_active="${dir_pg_output}/query2family"
			dir_genelist="${dir_pg_input}/query_gene"
			if [[ ! -d "${dir_genelist}" ]]; then
			  echo "Input directory does not exist: ${dir_genelist}"
			  exit 1
			fi
			files=()
			mapfile -t files < <(find "${dir_genelist}" -mindepth 1 -maxdepth 1 -type f ! -name '.*' | sort)
			if [[ ${#files[@]} -eq 0 ]]; then
			  echo "Input directory is empty: ${dir_genelist}"
			  exit 1
			fi
			if [[ ! "${SGE_TASK_ID}" =~ ^[0-9]+$ ]] || [[ ${SGE_TASK_ID} -lt 1 ]]; then
			  echo "Invalid SGE_TASK_ID value (must be a positive integer): ${SGE_TASK_ID}"
			  exit 1
			fi
			idx=$((SGE_TASK_ID-1))
			if [[ ${idx} -ge ${#files[@]} ]]; then
			  echo "Input genelist file not found, probably due to out-of-range SGE_TASK_ID=${SGE_TASK_ID}. Number of files: ${#files[@]}"
			  exit 1
			fi
			file_query_gene="${files[${idx}]}"
			og_id="$(basename "${file_query_gene}")"
			echo "mode_query2family: ${#files[@]} input genelist file(s) were detected in ${dir_genelist}/"
		echo "mode_query2family: Input genelist file = ${file_query_gene}"
		echo "output file prefix: ${og_id}"
		if [[ ! -f "${file_query_gene}" ]]; then
		  echo "Input genelist file not found, probably due to out-of-range SGE_TASK_ID=${SGE_TASK_ID}: ${file_query_gene}"
		  exit 1
		fi
		;;
	*)
		echo 'Exactly one of ${mode_orthogroup} or ${mode_query2family} should be set to 1. Exiting.'
		exit 1
		;;
esac
if [[ -z "$og_id" ]]; then
    echo "og_id is empty. Exiting."
    exit 1
fi

dir_sp_genome="${dir_pg_input}/species_genome"
dir_sp_gff="${dir_pg_input}/species_gff"
dir_sp_expression="${dir_pg_input}/species_expression"
dir_sp_cds="${dir_pg_input}/species_cds"
dir_sp_blastdb="${dir_pg_output}/species_cds_blastdb"
file_sp_trait="${dir_pg_input}/species_trait/species_trait.tsv"
file_og="${dir_pg_output}/orthofinder/Orthogroups_filtered/Orthogroups.selected.tsv"
file_og_parameters_dir="${dir_output_active}/parameters"
dir_species_tree="${dir_pg_output}/species_tree"
dir_species_tree_summary="${dir_species_tree}/species_tree_summary"
if [[ -s "${dir_species_tree_summary}/dated_species_tree.nwk" ]]; then
	species_tree_basename="dated_species_tree"
	species_tree="${dir_species_tree_summary}/${species_tree_basename}.nwk"
elif [[ -s "${dir_species_tree_summary}/undated_species_tree.nwk" ]]; then
	species_tree_basename="undated_species_tree"
	species_tree="${dir_species_tree_summary}/${species_tree_basename}.nwk"
else
	species_tree_basename="undated_species_tree"
	species_tree="${dir_species_tree_summary}/${species_tree_basename}.nwk"
fi
species_tree_generax="${file_og_parameters_dir}/${species_tree_basename}.generax.nwk" # generated later
species_tree_pruned="${file_og_parameters_dir}/${species_tree_basename}.pruned.nwk"
ensure_dir "${file_og_parameters_dir}"
notung_jar="/usr/local/bin/Notung.jar"
dir_rpsblastdb="/usr/local/db/Pfam_LE"

# Directory PATHs
# Directories for temporary files
dir_tmp="${dir_output_active}/tmp/${SGE_TASK_ID}_${og_id}" #_${RANDOM}

# File PATHs
# Alignment and gene tree preparation and others
file_og_query_aa_fasta="${dir_output_active}/query_aa_fasta/${og_id}_query.aa.fa.gz"
file_og_query_blast="${dir_output_active}/query_blast/${og_id}_query_blast.tsv"
file_og_cds_fasta="${dir_output_active}/cds_fasta/${og_id}_cds.fa.gz"
file_og_rpsblast="${dir_output_active}/rpsblast/${og_id}_rpsblast.tsv"
file_og_uniprot_annotation="${dir_output_active}/uniprot_annotation/${og_id}_uniprot.tsv"
file_og_mafft="${dir_output_active}/mafft/${og_id}_cds.aln.fa.gz"
file_og_maxalign="${dir_output_active}/maxalign/${og_id}_cds.maxalign.fa.gz"
file_og_trimal="${dir_output_active}/trimal/${og_id}_cds.trimal.fa.gz"
file_og_clipkit="${dir_output_active}/clipkit/${og_id}_cds.clipkit.fa.gz"
file_og_clipkit_log="${dir_output_active}/clipkit_log/${og_id}_cds.clipkit.log"
file_og_iqtree_tree="${dir_output_active}/iqtree_tree/${og_id}_iqtree.nwk"
file_og_iqtree_generax_ufboot="${dir_output_active}/generax_ufboot_tree/${og_id}_generax.ufboot.nwk"
file_og_orthogroup_extraction_nwk="${dir_output_active}/orthogroup_extraction_nwk/${og_id}_orthogroup_extraction.nwk"
file_og_orthogroup_extraction_fasta="${dir_output_active}/orthogroup_extraction_fasta/${og_id}_orthogroup_extraction.fa.gz"
file_og_generax_nhx="${dir_output_active}/generax_tree/${og_id}_generax.nhx"
file_og_generax_nwk="${dir_output_active}/generax_nwk/${og_id}_generax.nwk"
file_og_generax_xml="${dir_output_active}/generax_xml/${og_id}_generax.xml"
file_og_rooted_tree="${dir_output_active}/rooted_tree/${og_id}_root.nwk"
file_og_rooted_log="${dir_output_active}/rooted_tree_log/${og_id}_root.txt"
file_og_notung_reconcil="${dir_output_active}/notung_reconcile/${og_id}_notung_reconcile.zip"
file_og_dated_tree="${dir_output_active}/dated_tree/${og_id}_dated.nwk"
file_og_dated_tree_log="${dir_output_active}/dated_tree_log/${og_id}_dated.log.txt"
file_og_mapdnds_parameter="${dir_output_active}/mapdnds_parameter/${og_id}_parameter.zip"
file_og_mapdnds_dn="${dir_output_active}/mapdnds_dn_tree/${og_id}_mapdNdS.dN.nwk"
file_og_mapdnds_ds="${dir_output_active}/mapdnds_ds_tree/${og_id}_mapdNdS.dS.nwk"
file_og_codeml_two_ratio="${dir_output_active}/codeml_two_ratio/${og_id}_codeml.two_ratio.tsv"
file_og_hyphy_dnds="${dir_output_active}/hyphy_dnds/${og_id}_hyphy.dnds.json"
file_og_hyphy_relax="${dir_output_active}/hyphy_relax/${og_id}_hyphy.relax.json"
file_og_hyphy_relax_reversed="${dir_output_active}/hyphy_relax_reversed/${og_id}_hyphy.relax.reversed.json"
file_og_expression="${dir_output_active}/character_expression/${og_id}_expression.tsv"
file_og_gff_info="${dir_output_active}/character_gff_info/${og_id}_gff.tsv"
file_og_scm_intron_summary="${dir_output_active}/scm_intron_summary/${og_id}_scm.intron.tsv"
file_og_scm_intron_plot="${dir_output_active}/scm_intron_plot/${og_id}_scm.intron.pdf"
# Cis-regulatory motif
file_og_promoter_fasta="${dir_output_active}/promoter_fasta/${og_id}_promoter.fa.gz"
file_og_meme="${dir_output_active}/meme/${og_id}_meme.xml"
file_og_fimo="${dir_output_active}/fimo/${og_id}_fimo.tsv"
file_og_fimo_collapsed="${dir_output_active}/fimo_collapsed/${og_id}_fimo.collapsed.tsv"
# OU expression modeling
file_og_pem_rdata="${dir_output_active}/phylogeneticem_rdata/${og_id}_PhylogeneticEM.RData"
file_og_pem_tree="${dir_output_active}/phylogeneticem_tree/${og_id}_PhylogeneticEM.tree.tsv"
file_og_pem_regime="${dir_output_active}/phylogeneticem_regime/${og_id}_PhylogeneticEM.regime.tsv"
file_og_pem_leaf="${dir_output_active}/phylogeneticem_leaf/${og_id}_PhylogeneticEM.leaf.tsv"
file_og_pem_plot="${dir_output_active}/phylogeneticem_plot/${og_id}_PhylogeneticEM.pdf"
file_og_l1ou_fit_rdata="${dir_output_active}/l1ou_fit_rdata/${og_id}_l1ou.RData"
file_og_l1ou_fit_conv_rdata="${dir_output_active}/l1ou_fit_conv_rdata/${og_id}_l1ou.conv.RData"
file_og_l1ou_fit_tree="${dir_output_active}/l1ou_fit_tree/${og_id}_l1ou.tree.tsv"
file_og_l1ou_fit_regime="${dir_output_active}/l1ou_fit_regime/${og_id}_l1ou.regime.tsv"
file_og_l1ou_fit_leaf="${dir_output_active}/l1ou_fit_leaf/${og_id}_l1ou.leaf.tsv"
file_og_l1ou_fit_plot="${dir_output_active}/l1ou_fit_plot/${og_id}_l1ou.pdf"
# Protein convergence analysis
file_og_iqtree_anc="${dir_output_active}/iqtree_anc/${og_id}_iqtree.anc.zip"
file_og_csubst_b="${dir_output_active}/csubst_b/${og_id}_csubst_b.tsv"
file_og_csubst_cb_2="${dir_output_active}/csubst_cb_2/${og_id}_csubst_cb_2.tsv"
file_og_csubst_cb_stats="${dir_output_active}/csubst_cb_stats/${og_id}_csubst_cb_stats.tsv"
if [[ ${csubst_max_arity} -gt 2 ]]; then
  for (( i=3; i<=csubst_max_arity; i++ )); do
    declare file_og_csubst_cb_${i}="${dir_output_active}/csubst_cb_${i}/${og_id}.csubst_cb_${i}.tsv"
  done
fi
# PGLS output
file_og_gene_pgls="${dir_output_active}/pgls_gene_tree/${og_id}_gene_tree_PGLS.tsv"
file_og_gene_pgls_plot="${dir_output_active}/pgls_gene_tree_plot/${og_id}_gene_PGLS.barplot.pdf"
file_og_species_pgls="${dir_output_active}/pgls_species_tree/${og_id}_species_PGLS.tsv"
file_og_species_pgls_plot="${dir_output_active}/pgls_species_tree_plot/${og_id}_species_PGLS.barplot.pdf"
# Summary
file_og_stat_branch="${dir_output_active}/stat_branch/${og_id}_stat.branch.tsv"
file_og_stat_tree="${dir_output_active}/stat_tree/${og_id}_stat.tree.tsv"
file_og_amas_original="${dir_output_active}/amas_original/${og_id}_amas.original.tsv"
file_og_amas_cleaned="${dir_output_active}/amas_cleaned/${og_id}_amas.cleaned.tsv"
file_og_tree_plot="${dir_output_active}/tree_plot/${og_id}_tree_plot.pdf"
file_og_synteny="${dir_output_active}/synteny/${og_id}_synteny.tsv"
# Pruned datasets
file_og_untrimmed_aln_pruned="${dir_output_active}/pruned_untrimmed_alignment/${og_id}_cds.untrimmed.pruned.fa.gz"
file_og_trimmed_aln_pruned="${dir_output_active}/pruned_trimmed_alignment/${og_id}_cds.trimmed.pruned.fa.gz"
file_og_unrooted_tree_pruned="${dir_output_active}/pruned_unrooted_tree/${og_id}_unrooted.pruned.nwk"
file_og_rooted_tree_pruned="${dir_output_active}/pruned_rooted_tree/${og_id}_rooted.pruned.nwk"
file_og_dated_tree_pruned="${dir_output_active}/pruned_dated_tree/${og_id}_dated.pruned.nwk"

# Define intermediate files for downstream analysis.
# These variables are updated via helper functions to keep routing logic explicit.
set_analysis_file() {
  local slot=$1
  local path=$2
  case "${slot}" in
    untrimmed_aln) file_og_untrimmed_aln_analysis=${path} ;;
    trimmed_aln) file_og_trimmed_aln_analysis=${path} ;;
    unrooted_tree) file_og_unrooted_tree_analysis=${path} ;;
    rooted_tree) file_og_rooted_tree_analysis=${path} ;;
    dated_tree) file_og_dated_tree_analysis=${path} ;;
    *)
      echo "Error: Unknown analysis file slot: ${slot}"
      exit 1
      ;;
  esac
}

set_default_analysis_files() {
  set_analysis_file untrimmed_aln "${file_og_mafft}"
  set_analysis_file trimmed_aln "${file_og_mafft}"
  set_analysis_file unrooted_tree "${file_og_iqtree_tree}"
  set_analysis_file rooted_tree "${file_og_rooted_tree}"
  set_analysis_file dated_tree "${file_og_dated_tree}"
  if [[ ${run_generax} -eq 1 ]]; then
    if [[ -s "${species_tree_pruned}" ]]; then
      set_analysis_file rooted_tree "${file_og_generax_nhx}"
    else
      echo "run_generax is deactivated: missing species tree for GeneRax (${species_tree_pruned})"
      run_generax=0
    fi
  fi
}

switch_alignment_analysis_source() {
  local infile=$1
  set_analysis_file untrimmed_aln "${infile}"
  set_analysis_file trimmed_aln "${infile}"
}

prepare_species_tree_pruned() {
  local task_local="Species tree pruning"
  if [[ ! -s "${species_tree}" ]]; then
    echo "$(date): Warning: ${task_local}: source species tree was not found."
    echo "Missing: ${species_tree}"
    return 1
  fi

  if [[ -s "${species_tree_pruned}" ]]; then
    return 0
  fi

  ensure_parent_dir "${species_tree_pruned}"

  local cds_files=()
  mapfile -t cds_files < <(gg_find_fasta_files "${dir_sp_cds}" 1)
  if [[ ${#cds_files[@]} -eq 0 ]]; then
    echo "$(date): ${task_local}: no species CDS files detected. Copying source species tree as-is."
    cp_out "${species_tree}" "${species_tree_pruned}"
    return 0
  fi

  local cds_spp=()
  local cds_file
  for cds_file in "${cds_files[@]}"; do
    cds_spp+=( "$(gg_species_name_from_path "${cds_file}")" )
  done
  mapfile -t cds_spp < <(printf '%s\n' "${cds_spp[@]}" | sed -e '/^[[:space:]]*$/d' | sort -u)
  if [[ ${#cds_spp[@]} -eq 0 ]]; then
    echo "$(date): ${task_local}: species names could not be parsed from species CDS files. Copying source tree."
    cp_out "${species_tree}" "${species_tree_pruned}"
    return 0
  fi

  local keep_pattern
  keep_pattern=$(
    printf '%s\n' "${cds_spp[@]}" \
    | sed -e 's/[][(){}.^$+*?|\\-]/\\&/g' \
    | paste -sd'|' -
  )
  if [[ -z "${keep_pattern}" ]]; then
    echo "$(date): ${task_local}: keep-pattern is empty. Copying source tree."
    cp_out "${species_tree}" "${species_tree_pruned}"
    return 0
  fi

  local tmp_pruned="${species_tree_pruned}.tmp.$$"
  if nwkit prune \
    --infile "${species_tree}" \
    --pattern "^(${keep_pattern})$" \
    --invert_match yes \
    --outfile "${tmp_pruned}"; then
    if [[ -s "${tmp_pruned}" ]]; then
      mv_out "${tmp_pruned}" "${species_tree_pruned}"
    else
      echo "$(date): ${task_local}: pruned tree is empty. Copying source tree."
      rm -f -- "${tmp_pruned}"
      cp_out "${species_tree}" "${species_tree_pruned}"
    fi
  else
    echo "$(date): ${task_local}: nwkit prune failed. Copying source tree."
    rm -f -- "${tmp_pruned}"
    cp_out "${species_tree}" "${species_tree_pruned}"
  fi
}

prepare_species_tree_pruned || true
set_default_analysis_files

memory_notung=$((MEM_PER_HOST/4))

echo "Checking parameter conflicts..."
if [[ ${run_trimal} -eq 1 && ${run_clipkit} -eq 1 ]]; then
  echo 'Both of ${run_trimal} and ${run_clipkit} are set to 1.'
	echo '${run_trimal} is deactivated. ${run_clipkit} is still active.'
	run_trimal=0
fi
if [[ ! -d "${dir_sp_gff}" ]] || [[ -z "$(find "${dir_sp_gff}" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]]; then
  if [[ ${run_get_gff_info} -eq 1 ]]; then
    echo "\${run_get_gff_info} is deactivated. Empty input: ${dir_sp_gff}"
    run_get_gff_info=0
  fi
  if [[ ${run_scm_intron} -eq 1 ]]; then
    echo "\${run_scm_intron} is deactivated. Empty input: ${dir_sp_gff}"
    run_scm_intron=0
  fi
fi
if [[ -d "${dir_sp_expression}" ]] && [[ -n "$(find "${dir_sp_expression}" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]]; then
  echo "\${dir_sp_expression} is not empty. Continued: ${dir_sp_expression}"
else
  echo "\${dir_sp_expression} is empty: ${dir_sp_expression}"
  echo '${run_get_expression_matrix}, ${run_tree_pruning}, and other options are deactivated.'
  run_tree_pruning=0
  run_get_expression_matrix=0
  run_phylogeneticem=0
  run_l1ou=0
fi
if [[ -d "${dir_sp_genome}" ]] && [[ -n "$(find "${dir_sp_genome}" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]]; then
  echo "\${dir_sp_genome} is not empty. Continued: ${dir_sp_genome}"
else
  echo "\${dir_sp_genome} is empty: ${dir_sp_genome}"
  echo '${run_get_promoter_fasta} and ${run_fimo} are deactivated.'
  run_get_promoter_fasta=0
  run_fimo=0
fi
echo "Checking preexisting tmp directory."
if [[ -e "${dir_tmp}" && ${delete_preexisting_tmp_dir} -eq 1 ]]; then
	echo "$(date): Deleting preexisting ${dir_tmp}"
	shopt -s nullglob
	stale_tmp_paths=( "${dir_output_active}/tmp/${SGE_TASK_ID}_"* )
	shopt -u nullglob
	if [[ ${#stale_tmp_paths[@]} -gt 0 ]]; then
		rm -rf -- "${stale_tmp_paths[@]}"
	fi
fi
if [[ ! -e "${dir_tmp}" ]]; then
	echo "Making ${dir_tmp}"
	mkdir -p "${dir_tmp}"
fi
cd "${dir_tmp}"
echo "Working at: $(pwd)"

cleanup_tmp_dir_on_normal_exit() {
  local exit_code=$?
  if [[ ${delete_tmp_dir} -eq 1 && ( ${exit_code} -eq 0 || ${exit_code} -eq 8 ) ]]; then
    if [[ -n "${dir_tmp:-}" && -d "${dir_tmp}" && "${dir_tmp}" != "/" ]]; then
      echo "Deleting ${dir_tmp}"
      rm -rf -- "${dir_tmp}"
    elif [[ -n "${dir_tmp:-}" ]]; then
      echo "Refusing to delete unsafe tmp directory: ${dir_tmp}"
    fi
  fi
  return ${exit_code}
}
trap cleanup_tmp_dir_on_normal_exit EXIT

task="Species tree availability check"
if [[ ! -s "${species_tree_pruned}" ]]; then
  echo "$(date): Warning: ${task}: species tree file was not found."
  echo "Missing: ${species_tree_pruned}"
else
  gg_step_skip "${task}"
fi

task="Query fasta generation"
if [[ ! -s "${file_og_query_aa_fasta}" && ${run_get_query_fasta} -eq 1 ]]; then
    gg_step_start "${task}"
    if [[ "$(head -c 1 "${file_query_gene}")" == ">" ]]; then
        seqtype=$(seqkit stats --tabular "${file_query_gene}" | awk 'NR>1 {print $3}')
        if [[ ${seqtype} == "DNA" ]]; then
            echo "DNA sequences were detected. The file will be treated as in-frame CDS sequences, translated into amino acids, and used as a ${query_blast_method} query: ${file_query_gene}"
            seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" --threads "${NSLOTS}" "${file_query_gene}" > "${og_id}.query.aa.tmp.fasta"
            seqkit seq --threads "${NSLOTS}" "${og_id}.query.aa.tmp.fasta" --out-file "${og_id}.query.aa.out.fa.gz"
            mv_out "${og_id}.query.aa.out.fa.gz" "${file_og_query_aa_fasta}"
            rm -f -- "${og_id}.query.aa.tmp.fasta"
        elif [[ ${seqtype} == "Protein" ]]; then
            echo "Amino acid sequences were detected. The file will be used as a ${query_blast_method} query: ${file_query_gene}"
            seqkit seq --threads "${NSLOTS}" "${file_query_gene}" --out-file "${og_id}.query.aa.out.fa.gz"
            mv_out "${og_id}.query.aa.out.fa.gz" "${file_og_query_aa_fasta}"
        else
            echo "Unsupported sequence type '${seqtype}' in '${file_query_gene}'. Only \"DNA\" or \"Protein\" are allowed. Exiting."
            exit 1
        fi
    else
        echo "Gene IDs were detected. Extracting in-frame CDS sequences from species_cds: ${file_query_gene}"
        cp_out "${file_query_gene}" "${dir_output_active}/query_gene/$(basename "${file_query_gene}")"
        mapfile -t genes < <(sed -e '/^[[:space:]]*$/d' "${file_query_gene}")
	        mapfile -t cds_files < <(gg_find_fasta_files "${dir_sp_cds}" 1)
	        if [[ -e pattern.txt ]]; then
	            rm -f -- pattern.txt
	        fi
	        touch pattern.txt
        for gene in "${genes[@]}"; do
            echo "${gene}" >> pattern.txt
            if [[ "${gene}" == *"−"* ]]; then
                echo "Query sequence name contains minus sign. Searching the sequence name with hyphen as well: ${gene}"
                echo "${gene//−/-}" >> pattern.txt # Replace minus signs ("−") with hyphens ("-") and add to pattern.txt
            fi
        done
	        if [[ -e "${og_id}.query.cds.fasta" ]]; then
	            rm -f -- "${og_id}.query.cds.fasta"
	        fi
	        if [[ -e "${og_id}.query.cds.2.fasta" ]]; then
	            rm -f -- "${og_id}.query.cds.2.fasta"
	        fi
	        touch "${og_id}.query.cds.fasta"
	        query_hits_tmp_dir="./tmp.query_hits"
	        if [[ -d "${query_hits_tmp_dir}" ]]; then
	            rm -rf -- "${query_hits_tmp_dir}"
	        fi
	        mkdir -p "${query_hits_tmp_dir}"
	        for file_cds in "${cds_files[@]}"; do
            wait_until_jobn_le ${NSLOTS}
            (
                sp_ub=$(gg_species_name_from_path "${file_cds}")
                query_hits_tmp_file="${query_hits_tmp_dir}/$(basename "${file_cds}").hits.fasta"
                seqkit grep --threads "${NSLOTS}" --ignore-case --pattern-file <(awk -v sp="${sp_ub}" '{print $0; print sp "_" $0}' pattern.txt) "${file_cds}" \
                | sed -e "s/^>${sp_ub}_/>/" -e "s/^>${sp_ub}-/>/" -e "s/^>${sp_ub}[[:space:]]/>/" -e "s/^>${sp_ub}\./>/" -e "s/^>/>${sp_ub}_/" \
	                > "${query_hits_tmp_file}"
	            ) &
	        done
	        wait_for_background_jobs
	        shopt -s nullglob
	        query_hits_tmp_files=( "${query_hits_tmp_dir}"/*.hits.fasta )
	        shopt -u nullglob
	        for query_hits_tmp_file in "${query_hits_tmp_files[@]}"; do
	            cat "${query_hits_tmp_file}" >> "${og_id}.query.cds.fasta"
	        done
	        rm -rf -- "${query_hits_tmp_dir}"
	        gg_prepare_cds_fasta_stream "${NSLOTS}" "${genetic_code}" < "${og_id}.query.cds.fasta" \
	        | sed -e '/^1 1$/d' -e 's/_frame=1[[:space:]]*//' \
	        > "${og_id}.query.cds.2.fasta"
        num_query=${#genes[@]}
        num_result=$(grep -c -e "^>" "${og_id}.query.cds.2.fasta")
        echo "Number of gene names in query: ${num_query}"
        echo "Number of gene names in extracted fasta: ${num_result}"
	        if [[ ${num_query} -ne ${num_result} ]]; then
	            echo "Some gene names were not found in species_cds."
	            for gene_name in "${genes[@]}"; do
	                if ! awk -v gene="${gene_name}" '
	                  /^>/ {
	                    header=$0
	                    sub(/^>/, "", header)
	                    sub(/[[:space:]].*$/, "", header)
	                    if (header == gene) {
	                      found=1
	                      exit
	                    }
	                  }
	                  END { exit(found ? 0 : 1) }
	                ' "${og_id}.query.cds.2.fasta"; then
	                    echo "Query gene not found in species_cds: ${gene_name}"
	                fi
	            done
            echo "Exiting."
            exit 1
        fi
		        if [[ -s "${og_id}.query.cds.2.fasta" ]]; then
		            echo "Translating in-frame CDS sequences to amino acid sequences: ${og_id}.query.cds.2.fasta"
		            seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" --threads "${NSLOTS}" "${og_id}.query.cds.2.fasta" > "${og_id}.query.aa.tmp.fasta"
		            seqkit seq --threads "${NSLOTS}" "${og_id}.query.aa.tmp.fasta" --out-file "${og_id}.query.aa.out.fa.gz"
	            mv_out "${og_id}.query.aa.out.fa.gz" "${file_og_query_aa_fasta}"
	            rm -f -- "${og_id}.query.aa.tmp.fasta"
	            rm -f -- "${og_id}.query.cds.2.fasta"
	        fi
	    fi
else
    gg_step_skip "${task}"
fi

task="In-frame query BLAST (${query_blast_method})"
if [[ ! -s "${file_og_query_blast}" && ${run_query_blast} -eq 1 && ${mode_query2family} -eq 1 ]]; then
  gg_step_start "${task}"

  if [[ ${query_blast_method} == "tblastn" ]]; then
    if ! type makeblastdb >/dev/null 2>&1; then
      echo "makeblastdb was not found but query_blast_method=tblastn. Exiting."
      exit 1
    fi
    if ! type tblastn >/dev/null 2>&1; then
      echo "tblastn was not found but query_blast_method=tblastn. Exiting."
      exit 1
    fi
  elif [[ ${query_blast_method} == "diamond" ]]; then
    if ! type diamond >/dev/null 2>&1; then
      echo "diamond was not found but query_blast_method=diamond. Exiting."
      exit 1
    fi
    echo "DIAMOND mode selected: species CDS will be translated to proteins because diamond makedb/blastp use protein reference databases."
  fi

  export BLASTDB_LMDB_MAP_SIZE=100000000
  check_species_cds "${dir_pg}"
  check_if_species_files_unique "${dir_sp_cds}"

  if [[ -e "${og_id}".blastQuery.fasta ]]; then
    rm -f -- "${og_id}.blastQuery.fasta"
  fi
  touch "${og_id}.blastQuery.fasta"

  db_files=()
  ensure_dir "${dir_sp_blastdb}"
  cds_files=()
  mapfile -t cds_files < <(gg_find_fasta_files "${dir_sp_cds}" 1)
  cds_spp=()
  for cds_file in "${cds_files[@]}"; do
    cds_spp+=( "$(gg_species_name_from_path "${cds_file}")" )
  done
  mapfile -t cds_spp < <(printf "%s\n" "${cds_spp[@]}" | sort -u)
  for sp in "${cds_spp[@]}"; do
    wait_until_jobn_le ${NSLOTS}
    echo "sp: ${sp}"
    sp_cds_candidates=()
    for cds_file in "${cds_files[@]}"; do
      if [[ "$(gg_species_name_from_path "${cds_file}")" == "${sp}" ]]; then
        sp_cds_candidates+=( "${cds_file}" )
      fi
    done
    if [[ ${#sp_cds_candidates[@]} -eq 0 ]]; then
      echo "No CDS file was found for species: ${sp}. Skipping."
      continue
    fi
    mapfile -t sp_cds_candidates < <(printf "%s\n" "${sp_cds_candidates[@]}" | sort)
    sp_cds=${sp_cds_candidates[0]}
    sp_cds_blastdb="${dir_sp_blastdb}/$(basename "${sp_cds}")"
    db_files+=("${sp_cds_blastdb}")
    if [[ ${query_blast_method} == "tblastn" ]]; then
      echo "makeblastdb input CDS file: ${sp_cds}"
      echo "makeblastdb output database file: ${sp_cds_blastdb}"
      if [[ ! -e "${sp_cds_blastdb}".nhr || ! -e "${sp_cds_blastdb}".nin || ! -e "${sp_cds_blastdb}".nsq || ! -e "${sp_cds_blastdb}".ndb ]]; then
        db_lock_file="${sp_cds_blastdb}.tblastn.lock"
        (
          flock 9
          if [[ ! -e "${sp_cds_blastdb}".nhr || ! -e "${sp_cds_blastdb}".nin || ! -e "${sp_cds_blastdb}".nsq || ! -e "${sp_cds_blastdb}".ndb ]]; then
            if zgrep -q -e "^>.*[[:blank:]]" "${sp_cds}"; then
              echo "Space is detected. Please remove all annotation info after spaces in sequence names. Exiting: ${sp_cds}"
              exit 1
            fi
            if zgrep -q -e "^>.*[|]" "${sp_cds}"; then
              echo "Bar (|) is detected. Bars in sequence names will be replaced with underlines (_): ${sp_cds}"
            fi
            echo "Generating BLAST database: ${sp_cds}"
            echo "Generating BLAST database: ${sp_cds}" >&2
            if [[ ${sp_cds} =~ ".gz" ]]; then
	              seqkit seq --threads "${NSLOTS}" "${sp_cds}" | makeblastdb -dbtype nucl -title "${sp_cds}" -out "${sp_cds_blastdb}"
	            else
	              makeblastdb -dbtype nucl -in "${sp_cds}" -out "${sp_cds_blastdb}"
	            fi
          fi
        ) 9>"${db_lock_file}"
      fi
    elif [[ ${query_blast_method} == "diamond" ]]; then
      sp_cds_diamond_fasta="${sp_cds_blastdb}.diamond.fasta"
      echo "diamond input CDS file: ${sp_cds}"
      echo "diamond translated protein file: ${sp_cds_diamond_fasta}"
      echo "diamond database file: ${sp_cds_blastdb}.dmnd"
      if [[ ! -e "${sp_cds_blastdb}".dmnd ]]; then
        db_lock_file="${sp_cds_blastdb}.diamond.lock"
        (
          flock 9
          if [[ ! -e "${sp_cds_blastdb}".dmnd ]]; then
            if zgrep -q -e "^>.*[[:blank:]]" "${sp_cds}"; then
              echo "Space is detected. Please remove all annotation info after spaces in sequence names. Exiting: ${sp_cds}"
              exit 1
            fi
            if zgrep -q -e "^>.*[|]" "${sp_cds}"; then
              echo "Bar (|) is detected. Bars in sequence names will be replaced with underlines (_): ${sp_cds}"
            fi
            echo "Generating DIAMOND database: ${sp_cds}"
            echo "Generating DIAMOND database: ${sp_cds}" >&2
            if [[ ${sp_cds} =~ ".gz" ]]; then
              seqkit seq --remove-gaps --threads "${NSLOTS}" "${sp_cds}" \
              | seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" --threads "${NSLOTS}" \
              | sed '/^>/! s/\*//g' \
              > "${sp_cds_diamond_fasta}"
            else
              seqkit seq --remove-gaps --threads "${NSLOTS}" "${sp_cds}" \
              | seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" --threads "${NSLOTS}" \
              | sed '/^>/! s/\*//g' \
              > "${sp_cds_diamond_fasta}"
            fi
            if [[ "$(head -c 1 "${sp_cds_diamond_fasta}")" != '>' ]]; then
              sed -e "1d" "${sp_cds_diamond_fasta}" > "${sp_cds_diamond_fasta}.tmp"
              mv_out "${sp_cds_diamond_fasta}.tmp" "${sp_cds_diamond_fasta}"
            fi
            if [[ ! -s "${sp_cds_diamond_fasta}" ]]; then
              echo "Translated FASTA for DIAMOND is empty: ${sp_cds_diamond_fasta}. Exiting."
              exit 1
            fi
            if ! diamond makedb --in "${sp_cds_diamond_fasta}" --db "${sp_cds_blastdb}"; then
              echo "diamond makedb failed for ${sp_cds}. Exiting."
              exit 1
            fi
            rm -f -- "${sp_cds_diamond_fasta}"
          fi
        ) 9>"${db_lock_file}"
      fi
    fi
  done
  wait_for_background_jobs
  echo "db_files: ${db_files[*]}"
  query_aa_local="${og_id}.query.aa.tmp.for_blast.fasta"
  seqkit seq --threads "${NSLOTS}" "${file_og_query_aa_fasta}" --out-file "${query_aa_local}"

  outfmt="qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore frames qlen slen"
  if [[ ${query_blast_method} == "tblastn" ]]; then
    db_files_str=$(printf " %s" "${db_files[@]}")
    db_files_str="${db_files_str# }"
    echo "Running tblastn."
    if ! tblastn \
    -query "${query_aa_local}" \
    -db "${db_files_str}" \
    -out blast_out.tsv \
    -db_gencode "${genetic_code}" \
    -evalue "${query_blast_evalue}" \
    -max_target_seqs 50000 \
    -outfmt "6 ${outfmt}" \
    -num_threads "${NSLOTS}"; then
      echo "tblastn failed. Exiting."
      exit 1
    fi
  elif [[ ${query_blast_method} == "diamond" ]]; then
    echo "Running diamond blastp."
    rm -f -- blast_out.tsv
    touch blast_out.tsv
    for db_file in "${db_files[@]}"; do
      if [[ ! -e "${db_file}".dmnd ]]; then
        echo "DIAMOND database file is missing: ${db_file}.dmnd. Exiting."
        exit 1
      fi
      tmp_diamond_out="$(basename "${db_file}").diamond.out.tsv"
      if ! diamond blastp \
      --query "${query_aa_local}" \
      --db "${db_file}" \
      --out "${tmp_diamond_out}" \
      --evalue "${query_blast_evalue}" \
      --max-target-seqs 50000 \
      --threads "${NSLOTS}" \
      --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen; then
        echo "diamond blastp failed for database: ${db_file}. Exiting."
        exit 1
      fi
      if [[ -s "${tmp_diamond_out}" ]]; then
        awk -F '\t' 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,"0/1",$13,$14}' "${tmp_diamond_out}" >> blast_out.tsv
      fi
      rm -f -- "${tmp_diamond_out}"
    done
  fi
  rm -f -- "${query_aa_local}"

  python "${dir_script}/annotate_blast_coverage.py" \
  --in blast_out.tsv \
  --ncpu "${NSLOTS}" \
  --outfmt-columns "${outfmt}" \
  --frame-filter "0/1" \
  --out blast_out_inframe.tmp3.tsv

  if [[ -s blast_out_inframe.tmp3.tsv ]]; then
    mv_out blast_out_inframe.tmp3.tsv "${file_og_query_blast}"
  else
    echo "No query BLAST hits were detected after in-frame filtering. Exiting."
    exit 1
  fi
else
	gg_step_skip "${task}"
fi

task="Fasta generation"
if [[ ! -s "${file_og_cds_fasta}" && ${run_get_fasta} -eq 1 ]]; then
	gg_step_start "${task}"

	if [[ ${mode_orthogroup} -eq 1 ]]; then
    genes=()
	    read -r -a genes <<< "$(awk -v og="${og_id}" '$1==og {$1=""; sub(/^[[:space:]]+/, "", $0); gsub(",", "", $0); gsub(/\t/, " ", $0); sub(/[[:space:]]*$/, "", $0); gsub(/\047|"/, "", $0); print; exit}' "${file_og}")"
	elif [[ ${mode_query2family} -eq 1 ]]; then
	  python "${dir_script}/extract_gene_id_from_blast_table.py" \
	  --infile "${file_og_query_blast}" \
	  --outfile gene_id_list.txt \
	  --min_query_blast_coverage "${query_blast_coverage}" \
	  --max_num_gene_blast_hit_retrieval "${max_num_gene_blast_hit_retrieval}"
		mapfile -t genes < gene_id_list.txt
	fi
  cds_files=()
  mapfile -t cds_files < <(gg_find_fasta_files "${dir_sp_cds}" 1)
  echo "Number of CDS files in ${dir_sp_cds}: ${#cds_files[@]}"
  spp=()
	for file_cds in "${cds_files[@]}"; do
	    sp_ub=$(gg_species_name_from_path "${file_cds}")
	    spp+=("${sp_ub}")
	  done
		if [[ -e pattern.txt ]]; then
			rm -f -- pattern.txt
		fi
		touch pattern.txt
	for gene in "${genes[@]}"; do
    echo "${gene}" >> pattern.txt
	done
		if [[ -e "${og_id}.cds.fasta" ]]; then
			rm -f -- "${og_id}.cds.fasta"
		fi
		touch "${og_id}.cds.fasta"
	  for file_cds in "${cds_files[@]}"; do
	    sp_ub=$(gg_species_name_from_path "${file_cds}")
	    seqkit grep --threads "${NSLOTS}" --pattern-file pattern.txt "${file_cds}" \
	    >> "${og_id}.cds.fasta"
	  done

  seqkit replace --pattern "X" --replacement "N" --by-seq --ignore-case --threads "${NSLOTS}" "${og_id}.cds.fasta" \
  | seqkit replace --pattern " .*" --replacement "" --ignore-case --threads "${NSLOTS}" \
	  | seqkit replace --pattern "\+" --replacement "_" --ignore-case --threads "${NSLOTS}" \
	  | cdskit pad --codontable "${genetic_code}" \
	  | sed -e '/^1 1$/d' -e 's/_frame=1[[:space:]]*//' \
	  > "${og_id}.cds.2.fasta"

  num_gene=${#genes[@]}
  fasta_genes=()
  mapfile -t fasta_genes < <(awk '/^>/ {sub(/^>/, "", $0); print}' "${og_id}.cds.2.fasta")
  num_seq=${#fasta_genes[@]}
  echo "Number of genes in the orthogroup or BLAST hit: ${num_gene}"
  echo "Number of sequences in the fasta: ${num_seq}"
  if [[ ${num_gene} -eq ${num_seq} ]]; then
      echo "Number of genes and sequences matched. Fasta generation completed!"
      seqkit seq --threads "${NSLOTS}" "${og_id}.cds.2.fasta" --out-file "${og_id}.cds.out.fa.gz"
      mv_out "${og_id}.cds.out.fa.gz" "${file_og_cds_fasta}"
  else
      echo "Number of genes and sequences did not match."
      echo "Genes in the orthogroup or BLAST hit:"
      printf '%s\n' "${genes[@]}"
      echo ""
      echo "Genes in the generated FASTA:"
			printf '%s\n' "${fasta_genes[@]}" | sort | tr '\n' ' '
      echo ""
      echo "There may be duplicated or missing sequences."
      echo "If you have recently replaced species_cds files, please make sure to remove species_cds_blastdb before rerunning."
      echo "Exiting."
      exit 1
  fi
else
	gg_step_skip "${task}"
fi

task="Protein RPS-BLAST"
disable_if_no_input_file "run_rps_blast" "${file_og_cds_fasta}"
if [[ ! -s "${file_og_rpsblast}" && ${run_rps_blast} -eq 1 ]]; then
    gg_step_start "${task}"
    if ! dir_rpsblastdb=$(ensure_pfam_le_db "${dir_pg}"); then
        echo "Failed to prepare Pfam_LE DB. Exiting."
        exit 1
    fi
    if [[ -s "${dir_rpsblastdb}/Pfam.pal" ]]; then
        db_rpsblast="${dir_rpsblastdb}/Pfam"
    else
        rps_db_candidates=()
        mapfile -t rps_db_candidates < <(find "${dir_rpsblastdb}" -maxdepth 1 -type f -name "*.loo" | sort)
        if [[ ${#rps_db_candidates[@]} -eq 0 ]]; then
            echo "No RPS-BLAST DB index (*.loo) was found in: ${dir_rpsblastdb}. Exiting."
            exit 1
        fi
        if [[ ${#rps_db_candidates[@]} -gt 1 ]]; then
            echo "Multiple RPS-BLAST DB indices were found. Using the first sorted entry: ${rps_db_candidates[0]}"
        fi
        db_rpsblast="${rps_db_candidates[0]%.loo}"
    fi
    echo "db_rpsblast: ${db_rpsblast}"

	    if [[ -e "${og_id}.rpsblast.tmp.tsv" ]]; then
	        rm -f -- "${og_id}.rpsblast.tmp.tsv"
	    fi

    seqkit seq --remove-gaps --threads "${NSLOTS}" "${file_og_cds_fasta}" \
    | seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" --threads "${NSLOTS}" \
    > ungapped_translated_cds.fas

    if [[ "$(head -c 1 ungapped_translated_cds.fas)" != '>' ]]; then
        sed -e "1d" ungapped_translated_cds.fas > ungapped_translated_cds2.fas
        mv_out ungapped_translated_cds2.fas ungapped_translated_cds.fas
    fi

    outfmt="qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle"

	    if ! rpsblast \
	    -query ungapped_translated_cds.fas \
	    -db "${db_rpsblast}" \
	    -out "${og_id}.rpsblast.tmp.tsv" \
	    -evalue 0.01 \
    -outfmt "6 ${outfmt}" \
    -num_threads "${NSLOTS}"; then
        echo "RPS-BLAST failed. Exiting."
        exit 1
    fi

    genes=()
    mapfile -t genes < <(awk '/^>/ {sub(/^>/, "", $0); sub(/^[[:space:]]*/, "", $0); sub(/[[:space:]].*$/, "", $0); print}' ungapped_translated_cds.fas)
    for gene in "${genes[@]}"; do
        if ! awk -F '\t' -v gene="${gene}" '$1 == gene {found=1; exit} END {exit(found ? 0 : 1)}' "${og_id}.rpsblast.tmp.tsv"; then
            echo "${gene}: no hit in RPS-BLAST. Appending qlen to output tsv."
            qlen=$(seqkit fx2tab --length ungapped_translated_cds.fas | awk -F '\t' -v gene="${gene}" '$1 == gene {print $NF}')
	            printf '%s\t\t\t\t\t\t\t\t\t\t\t\t%s\t\t\n' "${gene}" "${qlen}" >> "${og_id}.rpsblast.tmp.tsv"
        else
            echo "${gene}: RPS-BLAST hit found."
        fi
    done
	    {
	      printf '%s\n' "${outfmt}" | tr ' ' '\t'
	      cat "${og_id}.rpsblast.tmp.tsv"
	    } > "${og_id}.rpsblast.tsv"

	    cp_out "${og_id}.rpsblast.tsv" "${file_og_rpsblast}"
else
	gg_step_skip "${task}"
fi

task="Gene trait extraction from gff files"
disable_if_no_input_file "run_get_gff_info" "${file_og_cds_fasta}"
	if [[ ! -s "${file_og_gff_info}" && ${run_get_gff_info} -eq 1 ]]; then
	  gg_step_start "${task}"
	  if [[ -e gff2genestat.tsv ]]; then
	    rm -f -- gff2genestat.tsv
	  fi
  seqkit seq --threads "${NSLOTS}" "${file_og_cds_fasta}" --out-file "${og_id}.gff2genestat_input.fasta"

  python "${dir_script}/gff2genestat.py" \
  --dir_gff "${dir_sp_gff}" \
  --feature "CDS" \
  --multiple_hits "longest" \
  --seqfile "${og_id}.gff2genestat_input.fasta" \
  --ncpu "${NSLOTS}" \
  --outfile gff2genestat.tsv
  rm -f -- "${og_id}.gff2genestat_input.fasta"

	  if [[ -s gff2genestat.tsv ]]; then
	    mv_out gff2genestat.tsv "${file_og_gff_info}"
	  fi
else
	gg_step_skip "${task}"
fi

task="UniProt annotation by DIAMOND"
disable_if_no_input_file "run_uniprot_annotation" "${file_og_cds_fasta}"
if [[ ! -s "${file_og_uniprot_annotation}" && ${run_uniprot_annotation} -eq 1 ]]; then
    gg_step_start "${task}"
    if ! uniprot_db_prefix=$(ensure_uniprot_sprot_db "${dir_pg}"); then
      echo "Failed to prepare UniProt Swiss-Prot DB. Exiting."
      exit 1
    fi

    seqkit seq --remove-gaps --only-id --threads "${NSLOTS}" "${file_og_cds_fasta}" \
    | seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" --threads "${NSLOTS}" \
    > uniprot.query.pep.fas

	    diamond blastp \
	    --query uniprot.query.pep.fas \
	    --threads "${NSLOTS}" \
	    --db "${uniprot_db_prefix}" \
    --very-sensitive \
    --out uniprot.diamond.tsv \
    --outfmt 6 qseqid sseqid pident length evalue bitscore qlen \
    --max-target-seqs 1 \
    --evalue 1e-2

	    python "${dir_script}/reformat_uniprot_diamond.py" \
	    --diamond_tsv uniprot.diamond.tsv \
	    --query_fasta uniprot.query.pep.fas \
	    --uniprot_fasta "${uniprot_db_prefix}.pep" \
	    --outfile uniprot.annotation.tsv

    cp_out uniprot.annotation.tsv "${file_og_uniprot_annotation}"
else
		gg_step_skip "${task}"
fi

task="In-frame mafft alignment"
disable_if_no_input_file "run_mafft" "${file_og_cds_fasta}"
if [[ ! -s "${file_og_mafft}" && ${run_mafft} -eq 1 ]]; then
  gg_step_start "${task}"

  seqkit seq --threads "${NSLOTS}" "${file_og_cds_fasta}" --out-file tmp.cds.input.fasta
  cdskit mask --seqfile tmp.cds.input.fasta --codontable "${genetic_code}" --outfile tmp.cds.fasta

	seqkit translate \
	--allow-unknown-codon \
	--transl-table "${genetic_code}" \
	--threads "${NSLOTS}" \
	tmp.cds.fasta \
	> tmp.pep.fasta

	mafft \
	--auto \
	--amino \
	--thread "${NSLOTS}" \
	--quiet \
	tmp.pep.fasta \
	> tmp.pep.aln.fasta

	cdskit backalign \
	--seqfile tmp.cds.fasta \
	--aa_aln tmp.pep.aln.fasta \
	--codontable "${genetic_code}" \
	--outfile "${og_id}.cds.aln.fasta"

  seqkit seq --threads "${NSLOTS}" "${og_id}.cds.aln.fasta" --out-file "${og_id}.cds.aln.out.fa.gz"
  mv_out "${og_id}.cds.aln.out.fa.gz" "${file_og_mafft}"
  rm -f -- tmp.cds.input.fasta
else
	gg_step_skip "${task}"
fi

task="AMAS for original alignment"
disable_if_no_input_file "run_amas_original" "${file_og_untrimmed_aln_analysis}"
if [[ ! -s "${file_og_amas_original}" && ${run_amas_original} -eq 1 ]]; then
	gg_step_start "${task}"
  seqkit seq --threads "${NSLOTS}" "${file_og_untrimmed_aln_analysis}" --out-file "${og_id}.amas.original.input.fasta"

	AMAS.py summary \
	--in-format fasta \
	--data-type dna \
	--in-files "${og_id}.amas.original.input.fasta"

	mv_out summary.txt "${file_og_amas_original}"
  rm -f -- "${og_id}.amas.original.input.fasta"
else
	gg_step_skip "${task}"
fi

task="MaxAlign"
disable_if_no_input_file "run_maxalign" "${file_og_untrimmed_aln_analysis}"
if [[ ! -s "${file_og_maxalign}" && ${run_maxalign} -eq 1 ]]; then
	gg_step_start "${task}"

	seqkit seq --threads "${NSLOTS}" "${file_og_untrimmed_aln_analysis}" --out-file "${og_id}.cds.aln.fasta"
	maxalign_keep_regex=""
	if [[ ${mode_query2family} -eq 1 && ${retain_query_in_maxalign} -eq 0 ]]; then
		echo "Query sequence(s) is NOT necessarily retained in MaxAlign."
	elif [[ ${mode_query2family} -eq 1 && ${retain_query_in_maxalign} -eq 1 ]]; then
		echo "Query sequence(s) is retained in MaxAlign."
		maxalign_keep_regex=$(python - "${file_query_gene}" <<'PY'
import re
import sys

infile = sys.argv[1]
gene_ids = []
with open(infile, 'r', encoding='utf-8', errors='replace') as handle:
    first_char = handle.read(1)
    handle.seek(0)
    if first_char == '>':
        for line in handle:
            if not line.startswith('>'):
                continue
            gene = line[1:].strip().split()[0]
            if gene:
                gene_ids.append(gene)
    else:
        for line in handle:
            gene = line.strip()
            if gene:
                gene_ids.append(gene)

normalized_ids = []
seen = set()
for gene in gene_ids:
    for candidate in (gene, gene.replace('−', '-')):
        if candidate and candidate not in seen:
            seen.add(candidate)
            normalized_ids.append(candidate)

patterns = [f"(?i:{re.escape(gene)}.*)" for gene in normalized_ids]
print(','.join(patterns))
PY
)
		if [[ -z "${maxalign_keep_regex}" ]]; then
			echo "Warning: No query IDs were parsed for MaxAlign --keep. Running without keep constraints."
		fi
	else
		maxalign_keep_regex=""
	fi

	maxalign_cmd=( \
		cdskit maxalign \
		--seqfile "${og_id}.cds.aln.fasta" \
		--outfile "${og_id}.maxalign.output.fasta" \
	)
	if [[ -n "${maxalign_keep_regex}" ]]; then
		maxalign_cmd+=(--keep "${maxalign_keep_regex}")
	fi
	"${maxalign_cmd[@]}"

	echo Number of sequences before MaxAlign: $(gg_count_fasta_records "${og_id}.cds.aln.fasta")
	echo Number of sequences after MaxAlign: $(gg_count_fasta_records "${og_id}.maxalign.output.fasta")

	seqkit seq --threads "${NSLOTS}" "${og_id}.maxalign.output.fasta" --out-file "${og_id}.maxalign.out.fa.gz"
	mv_out "${og_id}.maxalign.out.fa.gz" "${file_og_maxalign}"
	rm -f -- "${og_id}.maxalign.output.fasta"
else
	gg_step_skip "${task}"
fi
if [[ ${run_maxalign} -eq 1 ]]; then
    switch_alignment_analysis_source "${file_og_maxalign}"
fi

task="TrimAl"
disable_if_no_input_file "run_trimal" "${file_og_untrimmed_aln_analysis}"
if [[ ! -s "${file_og_trimal}" && ${run_trimal} -eq 1 ]]; then
	gg_step_start "${task}"

		seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" --threads "${NSLOTS}" "${file_og_untrimmed_aln_analysis}" \
		| sed -e '/^1 1$/d' -e 's/_frame=1[[:space:]]*//' \
		> untrimmed.pep.fasta

		seqkit seq --remove-gaps --threads "${NSLOTS}" \
		"${file_og_untrimmed_aln_analysis}" \
		> untrimmed.cds.degap.fasta

	trimal \
	-in untrimmed.pep.fasta \
	-backtrans untrimmed.cds.degap.fasta \
	-out "${og_id}.cds.trimal.tmp1.fasta" \
	-ignorestopcodon \
	-automated1

  cdskit rmseq \
  --seqfile "${og_id}.cds.trimal.tmp1.fasta" \
  --problematic_percent 100 \
  | cdskit hammer \
  --seqfile "-" \
  --codontable "${genetic_code}" \
  --nail 4 \
  --outfile "${og_id}.cds.trimal.tmp2.fasta"

  if [[ -s "${og_id}.cds.trimal.tmp2.fasta" ]]; then
    echo "Copying. Output file detected for the task: ${task}"
    seqkit seq --threads "${NSLOTS}" "${og_id}.cds.trimal.tmp2.fasta" --out-file "${og_id}.cds.trimal.out.fa.gz"
    mv_out "${og_id}.cds.trimal.out.fa.gz" "${file_og_trimal}"
  fi
else
	gg_step_skip "${task}"
fi
if [[ ${run_trimal} -eq 1 ]]; then
    set_analysis_file trimmed_aln "${file_og_trimal}"
fi

task="ClipKIT"
disable_if_no_input_file "run_clipkit" "${file_og_untrimmed_aln_analysis}"
if [[ ! -s "${file_og_clipkit}" && ${run_clipkit} -eq 1 ]]; then
	gg_step_start "${task}"

  seqkit seq --threads "${NSLOTS}" "${file_og_untrimmed_aln_analysis}" --out-file "${og_id}.cds.clipkit.input.fasta"

	clipkit \
	"${og_id}.cds.clipkit.input.fasta" \
	--mode smart-gap \
	--sequence_type nt \
	--codon \
	--input_file_format "fasta" \
	--output_file_format "fasta" \
	--output "${og_id}.cds.clipkit.tmp.fasta" \
	--log

	cdskit hammer \
  --codontable "${genetic_code}" \
  --nail 4 \
  --seqfile "${og_id}.cds.clipkit.tmp.fasta" \
  | cdskit rmseq \
  --problematic_percent 100 \
  --outfile "${og_id}.cds.clipkit.hammer.fasta"

  if [[ -s "${og_id}.cds.clipkit.hammer.fasta" ]]; then
    echo "Copying. Output file detected for the task: ${task}"
    seqkit seq --threads "${NSLOTS}" "${og_id}.cds.clipkit.hammer.fasta" --out-file "${og_id}.cds.clipkit.out.fa.gz"
    mv_out "${og_id}.cds.clipkit.out.fa.gz" "${file_og_clipkit}"
    cp_out "${og_id}.cds.clipkit.tmp.fasta.log" "${file_og_clipkit_log}"
  fi
  rm -f -- "${og_id}.cds.clipkit.input.fasta"
else
	gg_step_skip "${task}"
fi
if [[ ${run_clipkit} -eq 1 ]]; then
    set_analysis_file trimmed_aln "${file_og_clipkit}"
fi

task="AMAS for cleaned alignment"
disable_if_no_input_file "run_amas_cleaned" "${file_og_trimmed_aln_analysis}"
if [[ ! -s "${file_og_amas_cleaned}" && ${run_amas_cleaned} -eq 1 ]]; then
	gg_step_start "${task}"
  seqkit seq --threads "${NSLOTS}" "${file_og_trimmed_aln_analysis}" --out-file "${og_id}.amas.cleaned.input.fasta"

	AMAS.py summary \
	--in-format fasta \
	--data-type dna \
	--in-files "${og_id}.amas.cleaned.input.fasta"

	mv_out summary.txt "${file_og_amas_cleaned}"
  rm -f -- "${og_id}.amas.cleaned.input.fasta"
else
	gg_step_skip "${task}"
fi
if [[ ${run_maxalign} -eq 1 ]]; then
	# This code block should be placed immediately after "AMAS for cleaned alignment".
	# orthogroup_summary.tsv will not include necessary info otherwise.
    num_gene_before_maxalign=$(gg_count_fasta_records "${file_og_mafft}")
    num_gene_after_maxalign=$(gg_count_fasta_records "${file_og_maxalign}")
    echo "Number of genes before MaxAlign: ${num_gene_before_maxalign}"
    echo "Number of genes after MaxAlign: ${num_gene_after_maxalign}"
    if [[ ${num_gene_after_maxalign} -lt 3 ]]; then
        echo "This is not sufficient for tree-based analysis (<3). Exiting."
        exit 1
    fi
fi

task="IQ-TREE"
disable_if_no_input_file "run_iqtree" "${file_og_trimmed_aln_analysis}"
if [[ ! -s "${file_og_iqtree_tree}" && ${run_iqtree} -eq 1 ]]; then
	gg_step_start "${task}"
	num_seq=$(gg_count_fasta_records "${file_og_trimmed_aln_analysis}")
	if [[ ${num_seq} -ge 4 ]]; then
		if [[ ${run_generax} -eq 1 ]]; then
			other_iqtree_params=()
			file_tree="${og_id}.treefile"
			echo "run_generax=1: disabling UFBOOT in the initial IQ-TREE run. Support will be computed on the GeneRax topology."
		else
			other_iqtree_params=( --ufboot 1000 --bnni )
			file_tree="${og_id}.contree"
		fi
	else
		other_iqtree_params=()
		file_tree="${og_id}.treefile"
	fi
	if [[ ${num_seq} -gt ${iqtree_fast_mode_gt} ]]; then
	  other_iqtree_params+=( --fast )
	fi

		if [[ ${run_generax} -eq 1 ]]; then
			base_model=${generax_model%%+*}
		aa_models=( Blosum62 cpREV Dayhoff DCMut DEN FLU HIVb HIVw JTT JTT-DCMut LG mtART mtMAM mtREV mtZOA PMB rtREV stmtREV VT WAG LG4M LG4X PROTGTR )
		is_aa_model=0
		for aa_model in "${aa_models[@]}"; do
		  if [[ "${aa_model}" == "${base_model}" ]]; then
		    is_aa_model=1
		    break
		  fi
		done
		if [[ ${is_aa_model} -eq 1 ]]; then
			echo "Specified substitution model was interpreted as an amino acid model (base model = ${base_model})."
				seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" "${file_og_trimmed_aln_analysis}" \
				| sed -e '/^1 1$/d' -e 's/_frame=1[[:space:]]*//' \
				> iqtree_input.fa
		else
			echo "Specified substitution model was interpreted as a nucleotide model (base model = ${base_model})."
			seqkit seq --threads "${NSLOTS}" "${file_og_trimmed_aln_analysis}" --out-file iqtree_input.fa
		fi
	else
		seqkit seq --threads "${NSLOTS}" "${file_og_trimmed_aln_analysis}" --out-file iqtree_input.fa
	fi
	echo "IQ-TREE starting..."

	iqtree_model_string=${generax_model}
	build_iqtree_mem_args
	
	iqtree \
	-s iqtree_input.fa \
	-m "${iqtree_model_string}" \
	-T AUTO \
	--threads-max "${NSLOTS}" \
	--prefix "${og_id}" \
	"${IQTREE_MEM_ARGS[@]}" \
	--seed 12345 \
	--redo \
	"${other_iqtree_params[@]}"

	cp_out "${file_tree}" "${file_og_iqtree_tree}"
else
	gg_step_skip "${task}"
fi

task="Gene tree rooting"
disable_if_no_input_file "run_tree_root" "${file_og_unrooted_tree_analysis}"
if [[ ( ! -s "${file_og_rooted_tree}" || ! -s "${file_og_rooted_log}" ) && ${run_tree_root} -eq 1 ]]; then
	gg_step_start "${task}"

  if [[ "${tree_rooting_method}" == "notung" ]]; then
    if [[ ! -s "${species_tree_pruned}" ]]; then
      echo "tree_rooting_method=notung requires species tree: ${species_tree_pruned}"
      exit 1
    fi
    if [[ -e "./${og_id}.notung.root" ]]; then
      rm -rf -- "./${og_id}.notung.root"
    fi

    echo "memory_notung: ${memory_notung}"
    java -jar -Xmx${memory_notung}g "${notung_jar}" \
    -s "${species_tree_pruned}" \
    -g "${file_og_unrooted_tree_analysis}" \
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
    --outputdir "./${og_id}.notung.root"

    rooted_candidates=()
    mapfile -t rooted_candidates < <(find "./${og_id}.notung.root" -maxdepth 1 -type f -name "${og_id}.iqtree.nwk.rooting.*" | sort -V)
    selected_rooted_tree=""
    for candidate in "${rooted_candidates[@]}"; do
      if [[ "${candidate}" =~ \.rooting\.[0-9]+$ ]]; then
        selected_rooted_tree="${candidate}"
        break
      fi
    done
    if [[ -z "${selected_rooted_tree}" ]]; then
      echo "NOTUNG did not generate rooted-tree candidates in ./${og_id}.notung.root"
      exit 1
    fi

    nwkit label --target intnode --force yes --infile "${selected_rooted_tree}" --outfile "${og_id}.root.tmp.nwk"
    mv_out "${og_id}.root.tmp.nwk" "${file_og_rooted_tree}"
    {
      echo "tree_rooting_method=notung"
      echo "selected_rooting=${selected_rooted_tree}"
    } > "${og_id}.root.txt"
    mv_out "${og_id}.root.txt" "${file_og_rooted_log}"
  else
    nwkit_root_method="${tree_rooting_method}"
    if [[ "${nwkit_root_method}" == "md" ]]; then
      nwkit_root_method="mv"
    fi
    nwkit root --method "${nwkit_root_method}" --infile "${file_og_unrooted_tree_analysis}" \
    | nwkit label --target intnode --force yes --outfile "${og_id}.root.tmp.nwk"
    mv_out "${og_id}.root.tmp.nwk" "${file_og_rooted_tree}"
    {
      echo "tree_rooting_method=${tree_rooting_method}"
      echo "nwkit_method=${nwkit_root_method}"
    } > "${og_id}.root.txt"
    mv_out "${og_id}.root.txt" "${file_og_rooted_log}"
  fi
else
	gg_step_skip "${task}"
fi

task="Orthogroup extraction with NWKIT"
run_orthogroup_extraction_original="${run_orthogroup_extraction}" # This variable may be disabled by disable_if_no_input_file but the original value is necessary to properly update file_og_*_analysis
disable_if_no_input_file "run_orthogroup_extraction" "${file_query_gene:-}" "${file_og_trimmed_aln_analysis}" "${file_og_rooted_tree_analysis}"
if [[ ( ! -s "${file_og_orthogroup_extraction_nwk}" || ! -s "${file_og_orthogroup_extraction_fasta}" ) && ${run_orthogroup_extraction} -eq 1 ]]; then
	gg_step_start "${task}"

	  if [[ "$(head -c 1 "${file_query_gene}")" == ">" ]]; then
    echo "Fasta format was detected. Running run_orthogroup_extraction but gene names in the input fasta may not be compatible with this task."
    comma_separated_genes=$(awk '/^>/ {sub(/^>/, "", $0); sub(/[[:space:]].*$/, "", $0); gsub(/−/, "-", $0); print}' "${file_query_gene}" | paste -sd, -)
	  else
	    echo "Gene IDs were detected."
	    cp_out "${file_query_gene}" "${dir_output_active}/query_gene/$(basename "${file_query_gene}")"
		    comma_separated_genes=$(tr '\n' ',' < "${file_query_gene}" | sed -e 's/,$//' | tr '−' '-')
	  fi
  echo "Seed genes for orthogroup extraction: ${comma_separated_genes}"

  run_nwkit_subtree () {
    local infile=$1
    echo "Running nwkit subtree for ${infile}"
  	local info_txt=$(nwkit subtree --infile "${infile}" --leaves "${comma_separated_genes}" --orthogroup "yes" --dup_conf_score_threshold 0 2> /dev/null | nwkit info 2> /dev/null)
	  	local num_leaf
	  	num_leaf=$(awk -F': *' '/Number of leaves/ {print $2; exit}' <<< "${info_txt}")
    printf '%s\t%s\n' "${num_leaf}" "${infile}" >> tmp_num_leaf.tsv
  }

  subtree_infiles=()
  if [[ "${tree_rooting_method}" == "notung" && -d "./${og_id}.notung.root" ]]; then
    mapfile -t subtree_infiles < <(
      find "./${og_id}.notung.root" -maxdepth 1 -type f \
      | awk -v og="${og_id}" '$0 ~ (og "\\.iqtree\\.nwk\\.rooting\\.[0-9]+$") {print}' \
      | sort -V
    )
  fi
  if [[ ${#subtree_infiles[@]} -eq 0 ]]; then
    if [[ -s "${file_og_rooted_tree_analysis}" ]]; then
      subtree_infiles=( "${file_og_rooted_tree_analysis}" )
    else
      echo "No rooted tree is available for orthogroup extraction."
      exit 1
    fi
  fi

  printf 'num_leaf\tfile\n' > tmp_num_leaf.tsv
  for subtree_infile in "${subtree_infiles[@]}"; do
    wait_until_jobn_le ${NSLOTS}
    run_nwkit_subtree "${subtree_infile}"
  done

  if ! IFS=$'\t' read -r min_leaf_num min_leaf_file max_leaf_num max_leaf_file < <(
    awk -F'\t' '
      NR==1 {next}
      NR==2 {min=$1; max=$1; min_file=$2; max_file=$2; next}
      {
        if ($1 < min) {min=$1; min_file=$2}
        if ($1 > max) {max=$1; max_file=$2}
      }
      END {
        if (NR < 2) exit 1
        printf "%s\t%s\t%s\t%s\n", min, min_file, max, max_file
      }
    ' tmp_num_leaf.tsv
  ); then
    echo "Failed to parse tmp_num_leaf.tsv."
    exit 1
  fi
  echo "Minimum number of orthogroup subtree leaves after checking all rooting positions: ${min_leaf_num} in ${min_leaf_file} (will be used for orthogroup extraction)"
  echo "Maximum number of orthogroup subtree leaves after checking all rooting positions: ${max_leaf_num} in ${max_leaf_file} (shown just as a reference)"

	  nwkit subtree --infile "${min_leaf_file}" --leaves "${comma_separated_genes}" --orthogroup "yes" --dup_conf_score_threshold 0 \
	  --outfile "${og_id}.orthogroup_seed.tmp.nwk"

  seqkit seq --threads "${NSLOTS}" "${file_og_trimmed_aln_analysis}" --out-file tmp.trimmed.input.fasta
  nwkit intersection \
  --infile "${og_id}.orthogroup_seed.tmp.nwk" \
  --outfile /dev/null \
  --seqin tmp.trimmed.input.fasta \
	--seqout tmp.fasta \
	--match "complete"
		rm -f -- tmp.trimmed.input.fasta

	  # Preserve IQ-TREE support values in the extracted unrooted tree.
	  nwkit intersection \
	  --infile "${file_og_iqtree_tree}" \
	  --outfile "${og_id}.orthogroup_extraction.tmp.nwk" \
	  --seqin tmp.fasta \
	  --seqout /dev/null \
	  --match "complete"
	  mv_out "${og_id}.orthogroup_extraction.tmp.nwk" "${file_og_orthogroup_extraction_nwk}"
	  rm -f -- "${og_id}.orthogroup_seed.tmp.nwk"

	  cdskit hammer --nail 4 -s tmp.fasta -o "${og_id}.orthogroup_extraction.tmp.fasta"
	  seqkit seq --threads "${NSLOTS}" "${og_id}.orthogroup_extraction.tmp.fasta" --out-file "${og_id}.orthogroup_extraction.out.fa.gz"
	  mv_out "${og_id}.orthogroup_extraction.out.fa.gz" "${file_og_orthogroup_extraction_fasta}"
	  rm -f -- "${og_id}.orthogroup_extraction.tmp.fasta"
else
	gg_step_skip "${task}"
fi
if [[ ${run_orthogroup_extraction_original} -eq 1 && -s "${file_og_orthogroup_extraction_nwk}" && -s "${file_og_orthogroup_extraction_fasta}" ]]; then
  set_analysis_file unrooted_tree "${file_og_orthogroup_extraction_nwk}"
  set_analysis_file trimmed_aln "${file_og_orthogroup_extraction_fasta}"
fi

task="GeneRax"
disable_if_no_input_file "run_generax" "${file_og_trimmed_aln_analysis}" "${file_og_unrooted_tree_analysis}" "${species_tree_pruned}"
if [[ ! -s "${file_og_generax_nhx}" && ${run_generax} -eq 1 ]]; then
	gg_step_start "${task}"

	base_model=${generax_model%%+*}
	aa_models=( Blosum62 cpREV Dayhoff DCMut DEN FLU HIVb HIVw JTT JTT-DCMut LG mtART mtMAM mtREV mtZOA PMB rtREV stmtREV VT WAG LG4M LG4X PROTGTR )
	is_aa_model=0
	for aa_model in "${aa_models[@]}"; do
	  if [[ "${aa_model}" == "${base_model}" ]]; then
	    is_aa_model=1
	    break
	  fi
	done
	if [[ ${is_aa_model} -eq 1 ]]; then
		echo "Specified substitution model was interpreted as an amino acid model (base model = ${base_model})."
			seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" "${file_og_trimmed_aln_analysis}" \
			| sed -e '/^1 1$/d' -e 's/_frame=1[[:space:]]*//' \
			> generax_input_alignment.fas
	else
		echo "Specified substitution model was interpreted as a nucleotide model (base model = ${base_model})."
		seqkit seq --threads "${NSLOTS}" "${file_og_trimmed_aln_analysis}" --out-file generax_input_alignment.fas
	fi

  nwkit drop --target intnode --support yes --name yes \
	--infile "${file_og_unrooted_tree_analysis}" \
	--outfile generax_input_gene_tree.nwk

	#avoid multifurcating tree
	  R -q -e "library(ape); t=read.tree(\"generax_input_gene_tree.nwk\"); t=multi2di(t,random=FALSE); write.tree(t, \"generax_input_gene_tree_bi.nwk\")"

	generate_generax_mapfile () {
		# https://github.com/BenoitMorel/GeneRax/wiki/Gene-to-species-mapping
			my_aln_file=$1
				awk '/^>/ {sub(/^>/, "", $0); print}' "${my_aln_file}" > tmp.gene_names.txt
			while IFS= read -r gene_name; do gg_species_name_from_path "${gene_name}"; done < tmp.gene_names.txt > tmp.species_names.txt
			paste tmp.gene_names.txt tmp.species_names.txt > generax_map.txt
			rm -f -- tmp.gene_names.txt tmp.species_names.txt
	}
	generate_generax_mapfile generax_input_alignment.fas

		printf '%s\n' \
		  '[FAMILIES]' \
		  '- family_1' \
		  'starting_gene_tree = generax_input_gene_tree_bi.nwk' \
		  'alignment = generax_input_alignment.fas' \
		  'mapping = generax_map.txt' \
		  "subst_model = ${generax_model}" \
		  > generax_families.txt

  #${host_mpiexec_path} -np ${NSLOTS} generax \
  mpiexec_args=(mpiexec -oversubscribe -np ${NSLOTS})
  mpi_env_args=()
  running_under_scheduler=0
  if [[ -n "${SLURM_JOB_ID:-}" || -n "${PBS_JOBID:-}" || -n "${PE_HOSTFILE:-}" || -n "${LSB_JOBID:-}" ]]; then
    running_under_scheduler=1
  fi
  if [[ ${running_under_scheduler} -eq 0 ]]; then
    # Local/container runs may not have ssh/rsh and often probe unavailable OpenIB transports.
    mpi_env_args=(env OMPI_MCA_plm=isolated OMPI_MCA_plm_rsh_agent=/bin/false OMPI_MCA_btl=^openib)
  fi
  if [[ "$(id -u)" -eq 0 ]]; then
    mpiexec_args+=(--allow-run-as-root)
  fi
	"${mpi_env_args[@]}" "${mpiexec_args[@]}" generax \
	--species-tree "${species_tree_pruned}" \
	--families generax_families.txt \
	--strategy "SPR" \
	--rec-model "${generax_rec_model}" \
	--prefix "generax_${og_id}" \
	--per-family-rates \
  --skip-family-filtering \
	--mad-rooting \
	--seed 12345 </dev/null

  echo "GeneRax exit code = $?"

	generax_out_sptree="./generax_${og_id}/species_trees/starting_species_tree.newick" # generax v2.0
  if [[ -s "${generax_out_sptree}" ]]; then
    lock_file="${species_tree_generax}.lock"
	    if command -v flock >/dev/null 2>&1; then
	      exec 9> "${lock_file}"
	      flock 9
	      if [[ ! -s "${species_tree_generax}" ]]; then
	        echo "copying GeneRax output species tree (first writer only)."
	        cp_out "${generax_out_sptree}" "${species_tree_generax}"
	      fi
	      flock -u 9
	      exec 9>&-
	    else
	      echo "Error: flock command is required but not available."
	      exit 1
	    fi
	  elif [[ ! -s "${species_tree_generax}" ]]; then
	    echo "GeneRax species tree file was not found yet: ${generax_out_sptree}"
	  fi
	echo "copying GeneRax output gene tree."
	reconciled_base="./generax_${og_id}/reconciliations/family_1_reconciliated"
	reconciled_xml="${reconciled_base}.xml"
	reconciled_nhx="${reconciled_base}.nhx"
		if [[ -e "${reconciled_nhx}" ]]; then
			echo "GeneRax outfile was found. Copying."
			nwkit nhx2nwk --infile "${reconciled_nhx}" --outfile "${og_id}.generax.tmp.nwk"
			mv_out "${og_id}.generax.tmp.nwk" "${file_og_generax_nwk}"
			cp_out "${reconciled_xml}" "${file_og_generax_xml}"
			cp_out "${reconciled_nhx}" "${file_og_generax_nhx}"
	else
		echo "GeneRax outfile was not found. Exiting."
		exit 1
	fi
else
	gg_step_skip "${task}"
fi

task="IQ-TREE UFBOOT on GeneRax topology"
if [[ ${run_generax} -eq 1 ]]; then
	if [[ ! -s "${file_og_generax_nwk}" ]]; then
		echo "Skipped: ${task}. Missing GeneRax output tree: ${file_og_generax_nwk}"
	elif [[ ! -s "${file_og_trimmed_aln_analysis}" ]]; then
		echo "Skipped: ${task}. Missing alignment: ${file_og_trimmed_aln_analysis}"
	elif [[ ! -s "${file_og_iqtree_generax_ufboot}" || "${file_og_generax_nwk}" -nt "${file_og_iqtree_generax_ufboot}" || "${file_og_trimmed_aln_analysis}" -nt "${file_og_iqtree_generax_ufboot}" ]]; then
		gg_step_start "${task}"
		num_seq=$(gg_count_fasta_records "${file_og_trimmed_aln_analysis}")
		if [[ ${num_seq} -lt 4 ]]; then
			echo "UFBOOT requires >=4 sequences. Using the GeneRax topology without bootstrap support."
			nwkit drop --target root --length yes --infile "${file_og_generax_nwk}" --outfile "${og_id}.generax_ufboot.tmp.nwk"
			mv_out "${og_id}.generax_ufboot.tmp.nwk" "${file_og_iqtree_generax_ufboot}"
		else
			base_model=${generax_model%%+*}
			aa_models=( Blosum62 cpREV Dayhoff DCMut DEN FLU HIVb HIVw JTT JTT-DCMut LG mtART mtMAM mtREV mtZOA PMB rtREV stmtREV VT WAG LG4M LG4X PROTGTR )
			is_aa_model=0
			for aa_model in "${aa_models[@]}"; do
			  if [[ "${aa_model}" == "${base_model}" ]]; then
			    is_aa_model=1
			    break
			  fi
			done
			if [[ ${is_aa_model} -eq 1 ]]; then
				echo "Specified substitution model was interpreted as an amino acid model (base model = ${base_model})."
					seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" "${file_og_trimmed_aln_analysis}" \
					| sed -e '/^1 1$/d' -e 's/_frame=1[[:space:]]*//' \
					> "${og_id}.generax_ufboot.input.fa"
			else
				echo "Specified substitution model was interpreted as a nucleotide model (base model = ${base_model})."
				seqkit seq --threads "${NSLOTS}" "${file_og_trimmed_aln_analysis}" --out-file "${og_id}.generax_ufboot.input.fa"
			fi

			nwkit drop --target intnode --support yes --name yes --infile "${file_og_generax_nwk}" \
			| nwkit drop --target root --length yes --outfile "${og_id}.generax_ufboot.constraint.nwk"

			other_iqtree_params=( --ufboot 1000 --bnni )
			if [[ ${num_seq} -gt ${iqtree_fast_mode_gt} ]]; then
			  other_iqtree_params+=( --fast )
			fi

			build_iqtree_mem_args
			iqtree \
			-s "${og_id}.generax_ufboot.input.fa" \
			-g "${og_id}.generax_ufboot.constraint.nwk" \
			-m "${generax_model}" \
			-T AUTO \
			--threads-max "${NSLOTS}" \
			--prefix "${og_id}.generax_ufboot" \
			"${IQTREE_MEM_ARGS[@]}" \
			--seed 12345 \
			--redo \
			"${other_iqtree_params[@]}"

			if [[ -s "${og_id}.generax_ufboot.contree" ]]; then
				cp_out "${og_id}.generax_ufboot.contree" "${file_og_iqtree_generax_ufboot}"
			elif [[ -s "${og_id}.generax_ufboot.treefile" ]]; then
				cp_out "${og_id}.generax_ufboot.treefile" "${file_og_iqtree_generax_ufboot}"
			else
				echo "IQ-TREE UFBOOT on GeneRax topology failed to generate a tree file."
				exit 1
			fi
			rm -f -- "${og_id}.generax_ufboot.input.fa" "${og_id}.generax_ufboot.constraint.nwk"
		fi
	else
		gg_step_skip "${task}"
	fi
else
	gg_step_skip "${task}"
fi
if [[ ${run_generax} -eq 1 && -s "${file_og_iqtree_generax_ufboot}" ]]; then
	set_analysis_file unrooted_tree "${file_og_iqtree_generax_ufboot}"
fi

task="NOTUNG reconciliation"
disable_if_no_input_file "run_notung_reconcil" "${file_og_rooted_tree}" "${species_tree_pruned}"
if [[ ! -s "${file_og_notung_reconcil}" && ${run_notung_reconcil} -eq 1 ]]; then
	gg_step_start "${task}"

	echo "memory_notung: ${memory_notung}"

		if [[ -s "./${og_id}.root.nwk" ]]; then
			rm -f -- "${og_id}.root.nwk"
		fi
		if [[ -e "./${og_id}.notung_reconcile" ]]; then
			rm -rf -- "${og_id}.notung_reconcile"
		fi

  nwkit drop --target intnode --support yes --name yes \
	--infile "${file_og_rooted_tree}" \
	--outfile "${og_id}.root.nwk"

	java -jar -Xmx${memory_notung}g "${notung_jar}" \
	-s "${species_tree_pruned}" \
	-g "${og_id}.root.nwk" \
	--reconcile \
	--infertransfers "false" \
	--treeoutput newick \
	--log \
	--treestats \
	--events \
	--parsable \
	--speciestag prefix \
	--maxtrees 1 \
	--nolosses \
	--outputdir ./${og_id}.notung_reconcile

		if [[ -s "${og_id}.notung_reconcile/${og_id}.root.nwk.reconciled.parsable.txt" || -s "${og_id}.notung_reconcile/${og_id}.root.nwk.reconciled.0.parsable.txt" ]]; then
			zip -rq "${og_id}.notung_reconcile.zip" "${og_id}.notung_reconcile"
			cp_out "${og_id}.notung_reconcile.zip" "${file_og_notung_reconcil}"
		fi
else
	gg_step_skip "${task}"
fi

task="Species-tree-guided divergence time estimation"
disable_if_no_input_file "run_tree_dating" "${species_tree_pruned}" "${file_og_unrooted_tree_analysis}"
	if [[ ( ! -s "${file_og_dated_tree}" || ! -s "${file_og_dated_tree_log}" ) && ${run_tree_dating} -eq 1 ]]; then
		gg_step_start "${task}"
    radte_args=()

		if [[ ${run_generax} -eq 1 ]]; then
			radte_args+=("--species_tree=${species_tree_generax}")
			radte_args+=("--generax_nhx=${file_og_generax_nhx}")
		else
				if [[ -e "./${og_id}.notung_reconcile" ]]; then
					rm -rf -- "./${og_id}.notung_reconcile"
				fi
		cp_out "${file_og_notung_reconcil}" .
		unzip -q "$(basename "${file_og_notung_reconcil}")"
		if [[ -s ./${og_id}.notung_reconcile/${og_id}.root.nwk.reconciled.0 ]]; then
				cp_out ./"${og_id}".notung_reconcile/"${og_id}".root.nwk.reconciled.0 ./"${og_id}".notung_reconcile/"${og_id}".root.nwk.reconciled
				cp_out ./"${og_id}".notung_reconcile/"${og_id}".root.nwk.reconciled.0.parsable.txt ./"${og_id}".notung_reconcile/"${og_id}".root.nwk.reconciled.parsable.txt
			fi
				radte_args+=("--species_tree=${species_tree_pruned}")
				radte_args+=("--gene_tree=./${og_id}.notung_reconcile/${og_id}.root.nwk.reconciled")
				radte_args+=("--notung_parsable=./${og_id}.notung_reconcile/${og_id}.root.nwk.reconciled.parsable.txt")
		fi

		radte.r \
		"${radte_args[@]}" \
		--max_age="${radte_max_age}" \
		--chronos_lambda=1 \
		--chronos_model=discrete \
		--pad_short_edge=0.001 \
	2>&1 | tee radte.log

	constrained_node=$(awk -F': *' '/^Calibrated nodes:/ {gsub(/[[:space:]]/, "", $2); print $2; exit}' radte.log)
		echo "${constrained_node}" > "${og_id}.dated.log.txt"

	if grep -q ":-" radte_gene_tree_output.nwk; then
      contain_negative_bl=1
    else
      contain_negative_bl=0
    fi
	if [[ ${contain_negative_bl} -eq 1 ]]; then
			echo "Dated tree has negative branch length. Deleting output files depending on the tree file."
			for key in l1ou pem scm dated stat tree_plot; do
				files=()
				mapfile -t files < <(compgen -A variable "file_og_${key}")
				for f in "${files[@]}"; do
					target_file="${!f}"
					if [[ -e "${target_file}" ]]; then
						echo "deleting: ${target_file}"
							rm -f -- "${target_file}"
					fi
				done
			done
		else
			echo "Dated tree has no negative branch length. Continue."
			cp_out radte_calibrated_nodes.txt "${file_og_dated_tree_log}"
			cp_out radte_gene_tree_output.nwk "${file_og_dated_tree}"
		fi
else
	gg_step_skip "${task}"
fi

task="Expression matrix preparation"
disable_if_no_input_file "run_get_expression_matrix" "${file_og_trimmed_aln_analysis}"
if [[ ! -s "${file_og_expression}" && ${run_get_expression_matrix} -eq 1 ]]; then
	gg_step_start "${task}"
  seqkit seq --threads "${NSLOTS}" "${file_og_trimmed_aln_analysis}" --out-file "${og_id}.trait_matrix_input.fasta"

		python "${dir_script}/get_trait_matrix.py" \
		--dir_trait "${dir_sp_expression}" \
		--seqfile "${og_id}.trait_matrix_input.fasta" \
		--ncpu "${NSLOTS}" \
		--outfile expression_matrix.tsv
  rm -f -- "${og_id}.trait_matrix_input.fasta"
  if [[ -s expression_matrix.tsv ]]; then
  	mv_out expression_matrix.tsv "${file_og_expression}"
  fi
else
	gg_step_skip "${task}"
fi

task="Promoter fasta generation"
disable_if_no_input_file "run_get_promoter_fasta" "${file_og_gff_info}"
if [[ ! -s "${file_og_promoter_fasta}" && ${run_get_promoter_fasta} -eq 1 ]]; then
	gg_step_start "${task}"

	  python "${dir_script}/get_promoter_fasta.py" \
	  --dir_genome "${dir_sp_genome}" \
	  --geneinfo_tsv "${file_og_gff_info}" \
	  --seqkit_exe "seqkit" \
	  --outfile "${og_id}.promoter.tmp.fasta" \
	  --promoter_bp "${promoter_bp}" \
	  --ncpu "${NSLOTS}"
  if [[ -s "${og_id}.promoter.tmp.fasta" ]]; then
    seqkit seq --threads "${NSLOTS}" "${og_id}.promoter.tmp.fasta" --out-file "${og_id}.promoter.out.fa.gz"
    mv_out "${og_id}.promoter.out.fa.gz" "${file_og_promoter_fasta}"
    rm -f -- "${og_id}.promoter.tmp.fasta"
  fi
else
	gg_step_skip "${task}"
fi

task="fimo"
jaspar_path=""
if [[ ${run_fimo} -eq 1 ]]; then
  if ! jaspar_path=$(ensure_jaspar_file "${dir_pg}" "${jaspar_file}") || [[ -z "${jaspar_path}" ]]; then
    echo "Failed to prepare JASPAR motif file (${jaspar_file}). Exiting."
    exit 1
  fi
fi
disable_if_no_input_file "run_fimo" "${file_og_promoter_fasta}" "${jaspar_path}"
if [[ ! -s "${file_og_fimo}" && ${run_fimo} -eq 1 ]]; then
  gg_step_start "${task}"
  seqkit seq --threads "${NSLOTS}" "${file_og_promoter_fasta}" --out-file "${og_id}.fimo.input.fasta"

  fimo \
  --oc "fimo_out" \
  "${jaspar_path}" \
  "${og_id}.fimo.input.fasta"
  rm -f -- "${og_id}.fimo.input.fasta"

  fimo_result_table=""
  if [[ -s "./fimo_out/fimo.tsv" ]]; then
    fimo_result_table="./fimo_out/fimo.tsv"
  elif [[ -s "./fimo_out/fimo.txt" ]]; then
    fimo_result_table="./fimo_out/fimo.txt"
  fi

  if [[ -n "${fimo_result_table}" ]]; then
    mv_out "${fimo_result_table}" "${file_og_fimo}"
    rm -rf -- "./fimo_out"
  else
    echo "FIMO result table was not detected (expected fimo.tsv or fimo.txt). Keeping fimo_out for inspection."
  fi
else
	gg_step_skip "${task}"
fi

task="Tree pruning"
disable_if_no_input_file "run_tree_pruning" "${file_og_expression}" "${file_og_untrimmed_aln_analysis}" "${file_og_trimmed_aln_analysis}" "${file_og_unrooted_tree_analysis}" "${file_og_rooted_tree_analysis}"
if [[ ( ! -s "${file_og_untrimmed_aln_pruned}" || ! -s "${file_og_trimmed_aln_pruned}" || ! -s "${file_og_unrooted_tree_pruned}" || ! -s "${file_og_rooted_tree_pruned}" ) ]]; then
  is_all_outputs_exist=0
else
  if [[ ${run_tree_dating} -eq 1 && ! -s "${file_og_dated_tree_pruned}" ]]; then
    is_all_outputs_exist=0
  else
    is_all_outputs_exist=1
  fi
fi
if [[ ${is_all_outputs_exist} -eq 0 && ${run_tree_pruning} -eq 1 ]]; then
	gg_step_start "${task}"

	cut -f 1 "${file_og_expression}" | tail -n +2 > target_genes.txt

  if [[ -s "${file_og_untrimmed_aln_analysis}" ]]; then
    seqkit seq --threads "${NSLOTS}" "${file_og_untrimmed_aln_analysis}" --out-file "${og_id}.untrimmed.input.fasta"
    awk '
    NR==FNR {
      sub(/\r$/, "", $1)
      if ($1 != "") {
        keep[$1] = 1
      }
      next
    }
    /^>/ {
      id = substr($0, 2)
      sub(/[[:space:]].*$/, "", id)
      keep_seq = (id in keep)
    }
    keep_seq {
      print
    }
    ' target_genes.txt "${og_id}.untrimmed.input.fasta" > "${og_id}.untrimmed.pruned.tmp.fasta"
    rm -f -- "${og_id}.untrimmed.input.fasta"
    seqkit seq --threads "${NSLOTS}" "${og_id}.untrimmed.pruned.tmp.fasta" --out-file "${og_id}.untrimmed.pruned.out.fa.gz"
    mv_out "${og_id}.untrimmed.pruned.out.fa.gz" "${file_og_untrimmed_aln_pruned}"
    rm -f -- "${og_id}.untrimmed.pruned.tmp.fasta"
  fi

  if [[ -s "${file_og_trimmed_aln_analysis}" ]]; then
    seqkit seq --threads "${NSLOTS}" "${file_og_trimmed_aln_analysis}" --out-file "${og_id}.trimmed.input.fasta"
    awk '
    NR==FNR {
      sub(/\r$/, "", $1)
      if ($1 != "") {
        keep[$1] = 1
      }
      next
    }
    /^>/ {
      id = substr($0, 2)
      sub(/[[:space:]].*$/, "", id)
      keep_seq = (id in keep)
    }
    keep_seq {
      print
    }
    ' target_genes.txt "${og_id}.trimmed.input.fasta" > "${og_id}.trimmed.pruned.tmp.fasta"
    rm -f -- "${og_id}.trimmed.input.fasta"
    seqkit seq --threads "${NSLOTS}" "${og_id}.trimmed.pruned.tmp.fasta" --out-file "${og_id}.trimmed.pruned.out.fa.gz"
    mv_out "${og_id}.trimmed.pruned.out.fa.gz" "${file_og_trimmed_aln_pruned}"
    rm -f -- "${og_id}.trimmed.pruned.tmp.fasta"
  fi

  mapfile -t prune_genes < <(sed -e '/^[[:space:]]*$/d' target_genes.txt)
  prune_pattern=""
  if [[ ${#prune_genes[@]} -gt 0 ]]; then
    prune_pattern=$(
      printf '%s\n' "${prune_genes[@]}" \
      | sed -e 's/[][(){}.^$+*?|\\-]/\\&/g' \
      | paste -sd'|' -
    )
  fi

  if [[ -s "${file_og_unrooted_tree_analysis}" ]]; then
    if [[ ${#prune_genes[@]} -eq 0 ]]; then
      cat "${file_og_unrooted_tree_analysis}"
    else
      nwkit prune \
      --infile "${file_og_unrooted_tree_analysis}" \
      --pattern "^(${prune_pattern})$" \
      --invert_match yes
    fi \
    | nwkit drop --target root --length yes \
    | nwkit sanitize --remove_singleton yes --resolve_polytomy no \
    > "${og_id}.unrooted.pruned.tmp.nwk"
    mv_out "${og_id}.unrooted.pruned.tmp.nwk" "${file_og_unrooted_tree_pruned}"
  fi

  if [[ -s "${file_og_rooted_tree_analysis}" ]]; then
    if [[ ${#prune_genes[@]} -eq 0 ]]; then
      cat "${file_og_rooted_tree_analysis}"
    else
      nwkit prune \
      --infile "${file_og_rooted_tree_analysis}" \
      --pattern "^(${prune_pattern})$" \
      --invert_match yes
    fi \
    | nwkit drop --target root --length yes \
    | nwkit sanitize --remove_singleton yes --resolve_polytomy no \
    > "${og_id}.rooted.pruned.tmp.nwk"
    mv_out "${og_id}.rooted.pruned.tmp.nwk" "${file_og_rooted_tree_pruned}"
  fi

  if [[ -s "${file_og_dated_tree_analysis}" ]]; then
    if [[ ${#prune_genes[@]} -eq 0 ]]; then
      cat "${file_og_dated_tree_analysis}"
    else
      nwkit prune \
      --infile "${file_og_dated_tree_analysis}" \
      --pattern "^(${prune_pattern})$" \
      --invert_match yes
    fi \
    | nwkit drop --target root --length yes \
    | nwkit sanitize --remove_singleton yes --resolve_polytomy no \
    > "${og_id}.dated.pruned.tmp.nwk"
    mv_out "${og_id}.dated.pruned.tmp.nwk" "${file_og_dated_tree_pruned}"
  fi
else
	gg_step_skip "${task}"
fi
if [[ ${run_tree_pruning} -eq 1 ]]; then
	num_gene_before_pruning=$(gg_count_fasta_records "${file_og_trimmed_aln_analysis}")
	num_gene_after_pruning=$(gg_count_fasta_records "${file_og_trimmed_aln_pruned}")
    echo "Number of genes before pruning: ${num_gene_before_pruning}"
    echo "Number of genes after pruning: ${num_gene_after_pruning}"
	if [[ ${num_gene_after_pruning} -lt 3 ]]; then
			echo 'This is not sufficient for tree-based analysis (<3). Exiting.'
		exit 0
	fi
		set_analysis_file untrimmed_aln "${file_og_untrimmed_aln_pruned}"
		set_analysis_file trimmed_aln "${file_og_trimmed_aln_pruned}"
		set_analysis_file unrooted_tree "${file_og_unrooted_tree_pruned}"
		set_analysis_file rooted_tree "${file_og_rooted_tree_pruned}"
		set_analysis_file dated_tree "${file_og_dated_tree_pruned}"
    if [[ -s "${file_og_rooted_tree_analysis}" ]]; then
      assert_strictly_bifurcating_tree "${file_og_rooted_tree_analysis}" "Rooted analysis tree"
    fi
    if [[ -s "${file_og_dated_tree_analysis}" ]]; then
      assert_strictly_bifurcating_tree "${file_og_dated_tree_analysis}" "Dated analysis tree"
    fi
fi
if [[ -s "${file_og_expression}" && ( ${run_l1ou} -eq 1 || ${run_phylogeneticem} -eq 1 ) ]]; then
  # This block should be run after tree pruning.
  num_gene_trait=$(( $(wc -l < "${file_og_expression}") - 1 )) # -1 for header
  num_gene_tree=$(gg_count_fasta_records "${file_og_trimmed_aln_analysis}")
  if [[ ${num_gene_trait} -eq ${num_gene_tree} ]]; then
      echo "num_gene_trait (${num_gene_trait}) and num_gene_tree (${num_gene_tree}) matched."
  else
      echo "num_gene_trait (${num_gene_trait}) and num_gene_tree (${num_gene_tree}) did not match."
      if [[ ${run_tree_pruning} -ne 1 && ( ${run_phylogeneticem} -eq 1 || ${run_l1ou} -eq 1 ) ]]; then
        echo "Set run_tree_pruning=1 to run phylogenetic comparative analysis. Exiting."
        exit 1
    fi
  fi
fi

task="Checking whether downstream output files are correctly unpruned/pruned."
if [[ ${check_pruned} -eq 1 ]]; then
    gg_step_start "${task}"
    file_to_remove=(
        "${file_og_mapdnds_dn}"
        "${file_og_mapdnds_ds}"
        "${file_og_hyphy_dnds}"
        "${file_og_gff_info}"
        "${file_og_scm_intron_summary}"
        "${file_og_scm_intron_plot}"
        "${file_og_pem_rdata}"
        "${file_og_pem_tree}"
        "${file_og_pem_regime}"
        "${file_og_pem_leaf}"
        "${file_og_pem_plot}"
        "${file_og_l1ou_fit_rdata}"
        "${file_og_l1ou_fit_conv_rdata}"
        "${file_og_l1ou_fit_tree}"
        "${file_og_l1ou_fit_regime}"
        "${file_og_l1ou_fit_leaf}"
        "${file_og_l1ou_fit_plot}"
        "${file_og_iqtree_anc}"
        "${file_og_csubst_b}"
        "${file_og_csubst_cb_stats}"
        "${file_og_gene_pgls}"
        "${file_og_gene_pgls_plot}"
        "${file_og_species_pgls}"
        "${file_og_species_pgls_plot}"
    )
    for ((i=2; i<=csubst_max_arity; i++)); do
        varname="file_og_csubst_cb_${i}"
        if [[ -n "${!varname:-}" ]]; then
            file_to_remove+=( "${!varname}" )
        fi
    done
    remove_flag=0
    if [[ ${run_tree_pruning} -eq 0 ]]; then
        if [[ -s "${file_og_untrimmed_aln_pruned}" ]]; then
            remove_flag=1
        fi
    fi
    if [[ ${run_tree_pruning} -eq 1 ]]; then
        if [[ -s "${file_og_stat_tree}" ]]; then
            num_gene_pruned=$(gg_count_fasta_records "${file_og_trimmed_aln_pruned}")
            num_stat_branch_row=$(wc -l < "${file_og_stat_branch}")
            if [[ $((${num_gene_pruned}*2)) -lt ${num_stat_branch_row} ]]; then
                remove_flag=1
            fi
        fi
    fi
    if [[ ${remove_flag} -eq 1 ]]; then
        echo "Downstream output files are inconsistent with the run_tree_pruning setting."
        for file in "${file_to_remove[@]}"; do
            echo "Not found: ${file}"
	            if [[ -e "${file}" ]]; then
	                echo "Deleting: ${file}"
	                rm -f -- "${file}"
	            fi
	        done
        echo "Completed the deletion of inconsistent files."
    else
        echo "Downstream output files are consistent with the run_tree_pruning setting."
    fi
else
	gg_step_skip "${task}"
fi

task="Parameter estimation for mapdNdS"
disable_if_no_input_file "run_mapdnds_parameter_estimation" "${file_og_rooted_tree_analysis}" "${file_og_trimmed_aln_analysis}"
if [[ ! -s "${file_og_mapdnds_parameter}" && ${run_mapdnds_parameter_estimation} -eq 1 ]]; then
	gg_step_start "${task}"

  nwkit drop --target intnode --support yes --name yes \
	--infile "${file_og_rooted_tree_analysis}" \
  | nwkit sanitize --remove_singleton yes --resolve_polytomy no \
  > mapdnds_input.nwk

	seqkit seq --threads "${NSLOTS}" "${file_og_trimmed_aln_analysis}" --out-file ./mapdnds_input.fasta

  # F3X4+G4 shouldn not be changed otherwise iqtree2mapnh.py has to be updated.
	build_iqtree_mem_args
	iqtree \
	-s mapdnds_input.fasta \
	-m "GY+F3X4+G4" \
	-te mapdnds_input.nwk \
	-T AUTO \
	--threads-max "${NSLOTS}" \
	--seqtype "CODON${genetic_code}" \
	--prefix "${og_id}.iqtree2mapdNdS" \
	"${IQTREE_MEM_ARGS[@]}" \
	--ancestral \
	--seed 12345 \
	--redo

	python "${dir_script}/iqtree2mapnh.py" \
	--iqtree "${og_id}.iqtree2mapdNdS.iqtree" \
	--log "${og_id}.iqtree2mapdNdS.log" \
	--state "${og_id}.iqtree2mapdNdS.state" \
	--alignment mapdnds_input.fasta \
	--treefile "${og_id}.iqtree2mapdNdS.treefile" \
	--rooted_tree mapdnds_input.nwk \
	--genetic_code "${genetic_code}"

	if [[ -s "iqtree2mapnh.params" && -s "iqtree2mapnh.nwk" ]]; then
	  echo "iqtree2mapnh was successfully completed."
			  if [[ -e "${og_id}.mapdnds_parameter" ]]; then
			    rm -rf -- "${og_id}.mapdnds_parameter"
			  fi
			  mkdir -p "${og_id}.mapdnds_parameter"
		    mv_out "iqtree2mapnh.params" ./"${og_id}".mapdnds_parameter
		    mv_out "iqtree2mapnh.nwk" ./"${og_id}".mapdnds_parameter
		    zip -r "${og_id}.mapdnds_parameter.zip" "${og_id}.mapdnds_parameter"
			    mv_out "${og_id}.mapdnds_parameter.zip" "${file_og_mapdnds_parameter}"
  else
    echo "iqtree2mapnh.params was not generated."
  fi
else
	gg_step_skip "${task}"
fi

task="mapdNdS main run"
disable_if_no_input_file "run_mapdnds" "${file_og_mapdnds_parameter}" "${file_og_trimmed_aln_analysis}"
if [[ ( ! -s "${file_og_mapdnds_dn}" || ! -s "${file_og_mapdnds_ds}" ) && ${run_mapdnds} -eq 1 ]]; then
	gg_step_start "${task}"

		unzip -o "${file_og_mapdnds_parameter}"
		cd "${dir_tmp}/${og_id}.mapdnds_parameter"
	seqkit seq --threads "${NSLOTS}" "${file_og_trimmed_aln_analysis}" --out-file ./mapdnds_input.fasta
  normalize_mapnh_params_for_mapnh_v1 "iqtree2mapnh.params" "${genetic_code}"

  mapnh_exit_code=0
  if ! mapnh \
    SEQ=mapdnds_input.fasta \
    TREE=iqtree2mapnh.nwk \
    OUT=${og_id} \
    param=iqtree2mapnh.params \
    2>&1 | tee mapnh.log.txt; then
    mapnh_exit_code=${PIPESTATUS[0]:-1}
  fi

	  if [[ -s "${og_id}.dN.dnd" && -s "${og_id}.dS.dnd" ]]; then
	    echo "mapnh successfully generated dN and dS trees."
	    mv_out "${og_id}.dN.dnd" "${file_og_mapdnds_dn}"
	  	mv_out "${og_id}.dS.dnd" "${file_og_mapdnds_ds}"
  else
    echo "mapnh did not generate dN/dS trees (exit code: ${mapnh_exit_code})."
    echo "mapnh output and HyPhy output are managed separately; no cross-substitution is applied."
  fi
	  cd "${dir_tmp}"
else
	gg_step_skip "${task}"
fi

task="CodeML two-ratio model"
disable_if_no_input_file "run_codeml_two_ratio" "${file_og_rooted_tree_analysis}" "${file_og_trimmed_aln_analysis}" "${file_sp_trait}"
if [[ ! -s "${file_og_codeml_two_ratio}" && "${run_codeml_two_ratio}" -eq 1 ]]; then
	gg_step_start "${task}"

  binarize_species_trait "${file_sp_trait}" species_trait_binary.tsv
  sed '2,$ s/\t/_.*\t/' species_trait_binary.tsv > foreground.tsv
  IFS=$'\t' read -r -a colname_array < foreground.tsv

  if [[ ${#colname_array[@]} -le 1 ]]; then
    echo "No trait columns were detected in ${file_sp_trait}. Skipping ${task}."
  else
  for ((i=1; i<${#colname_array[@]}; i++)); do
    trait="${colname_array[$i]}"
    echo "Processing trait: ${trait}"
    awk -F'\t' -v trait_col="$((i+1))" 'NR>1 && $trait_col == 1 { print $1 }' foreground.tsv > "foreground_${trait}.txt"

	    target_spnode=$(paste -sd'|' "foreground_${trait}.txt") # Bar-separated list of target species nodes
    echo "Regular expression for CodeML foreground node search: ${target_spnode}"

    nwkit drop \
    --infile "${file_og_rooted_tree_analysis}" \
    --target "intnode" \
    --support "yes" \
    --name "yes" \
    | nwkit mark \
    --pattern "${target_spnode}" \
    --insert_txt "#1" \
    --insert_sep "" \
    --target "mrca" \
    --target_only_clade "yes" \
    --outfile "codeml_input_${trait}.nwk"
    echo "CodeML input tree: $(< "codeml_input_${trait}.nwk")"

    bash "${dir_script}/shorten_fasta_newick_names.sh" \
    "${file_og_trimmed_aln_analysis}" "codeml_input2_${trait}.fasta" "codeml_input_${trait}.nwk" "codeml_input2_${trait}.nwk" 90

    if grep -q "#1:0;$" "codeml_input_${trait}.nwk"; then
      exit_code1=0
    else
      exit_code1=1
    fi
    if grep -q "#1" "codeml_input_${trait}.nwk"; then
      exit_code2=0
    else
      exit_code2=1
    fi
    flag_unanalyzable=0
    if [[ ${exit_code1} -eq 0 || ${exit_code2} -eq 1 ]]; then
      flag_unanalyzable=1
    fi
    if [[ ${flag_unanalyzable} -eq 1 ]]; then
      echo "Target species tree node (${target_spnode}) is the gene tree root node. Generating an empty output file."
      codeml_out_treelength="NA"
      codeml_out_treelength_dn="NA"
      codeml_out_treelength_ds="NA"
      codeml_out_kappa="NA"
      codeml_out_background_omega="NA"
      codeml_out_foreground_omega="NA"
      codeml_out_time="NA"
    else
      python -c 'from pathlib import Path; import sys; template = Path(sys.argv[1]).read_text(encoding="utf-8"); rendered = template.replace("__SEQFILE__", sys.argv[2]).replace("__TREEFILE__", sys.argv[3]).replace("__ICODE__", sys.argv[4]); Path(sys.argv[5]).write_text(rendered, encoding="utf-8")' \
        "${dir_script}/codeml/codeml_two_ratio.ctl.template" \
        "codeml_input2_${trait}.fasta" \
        "codeml_input2_${trait}.nwk" \
        "$((genetic_code - 1))" \
        "my_codeml_${trait}.ctl"

      codeml "my_codeml_${trait}.ctl"
	      codeml_out_treelength=$(awk '/^tree length =/ {sub(/^tree length =[[:space:]]*/, "", $0); print; exit}' mlc)
	      codeml_out_treelength_dn=$(awk '/^tree length for dN:/ {sub(/^tree length for dN:[[:space:]]*/, "", $0); print; exit}' mlc)
	      codeml_out_treelength_ds=$(awk '/^tree length for dS:/ {sub(/^tree length for dS:[[:space:]]*/, "", $0); print; exit}' mlc)
	      codeml_out_kappa=$(awk '/^kappa \(ts\/tv\) =/ {sub(/^kappa \(ts\/tv\) =[[:space:]]*/, "", $0); print; exit}' mlc)
	      read -r -a codeml_out_omegas <<< "$(awk '/^w \(dN\/dS\) for branches:/ {sub(/^w \(dN\/dS\) for branches:[[:space:]]*/, "", $0); print; exit}' mlc)"
	      codeml_out_background_omega=${codeml_out_omegas[0]:-}
	      codeml_out_foreground_omega=${codeml_out_omegas[1]:-}
	      codeml_out_time=$(awk '/^Time used:/ {sub(/^Time used:[[:space:]]*/, "", $0); print; exit}' mlc)
    fi
    if [[ -n "${codeml_out_background_omega}" && -n "${codeml_out_foreground_omega}" ]]; then
      echo "The task '${task}' has completed successfully for trait '${trait}'."
      printf 'tree_length_%s\ttree_length_dn_%s\ttree_length_ds_%s\tkappa_%s\tbackground_omega_%s\tforeground_omega_%s\tcodeml_time_%s\n' \
        "${trait}" "${trait}" "${trait}" "${trait}" "${trait}" "${trait}" "${trait}" > "file_og_codeml_two_ratio_${trait}.tsv"
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${codeml_out_treelength}" "${codeml_out_treelength_dn}" "${codeml_out_treelength_ds}" "${codeml_out_kappa}" "${codeml_out_background_omega}" "${codeml_out_foreground_omega}" "${codeml_out_time}" >> "file_og_codeml_two_ratio_${trait}.tsv"
    else
      echo "The task '${task}' failed for trait '${trait}'."
    fi
  done

  # Combine all codeml output files
  missing_files=()
  codeml_output_files=()
  for ((i=1; i<${#colname_array[@]}; i++)); do
    trait="${colname_array[$i]}"
    codeml_output_file="file_og_codeml_two_ratio_${trait}.tsv"
    if [[ -s "${codeml_output_file}" ]]; then
      codeml_output_files+=("${codeml_output_file}")
    else
      missing_files+=("${codeml_output_file}")
    fi
  done
  if [[ ${#missing_files[@]} -gt 0 ]]; then
    echo "The following codeml output files are missing:"
    for f in "${missing_files[@]}"; do
      echo "${f}"
    done
    echo "The task has failed: ${task}"
  else
    echo "All codeml two-ratio model output files are generated. Combining them into a single tsv file: ${file_og_codeml_two_ratio}"
    header_files=()
    data_files=()
    for f in "${codeml_output_files[@]}"; do
      base=$(basename "$f")
      head -n1 "$f" > "header_${base}"
      tail -n1 "$f" > "data_${base}"
      header_files+=("header_${base}")
      data_files+=("data_${base}")
    done
	    paste -d$'\t' "${header_files[@]}" > combined_header.tsv
	    paste -d$'\t' "${data_files[@]}" > combined_data.tsv
	    cat combined_header.tsv combined_data.tsv > "${og_id}.codeml.two_ratio.tmp.tsv"
	    mv_out "${og_id}.codeml.two_ratio.tmp.tsv" "${file_og_codeml_two_ratio}"
	    echo "The task has completed successfully: ${task}"
	  fi
  fi
else
	gg_step_skip "${task}"
fi

task="HyPhy dN-dS estimation"
disable_if_no_input_file "run_hyphy_dnds" "${file_og_rooted_tree_analysis}" "${file_og_trimmed_aln_analysis}"
if [[ ! -s "${file_og_hyphy_dnds}" && ${run_hyphy_dnds} -eq 1 ]]; then
	gg_step_start "${task}"

	nwkit drop --target intnode --support yes --name yes \
	--infile "${file_og_rooted_tree_analysis}" \
  | nwkit sanitize --remove_singleton yes --resolve_polytomy no \
  > "hyphy_input.nwk"

	seqkit seq --threads "${NSLOTS}" "${file_og_trimmed_aln_analysis}" --out-file "hyphy_input.fasta"

	      hyphy_genetic_code=$(get_hyphy_genetic_code "${genetic_code}")


  hyphy "${dir_script}/hyphy/FitMG94.bf" \
  --alignment "hyphy_input.fasta" \
  --tree "hyphy_input.nwk" \
  --code "${hyphy_genetic_code}" \
  --frequencies "CF3x4" \
  --type "local" \
  --lrt "No" \
  --rooted "Yes" \
  --CPU "${NSLOTS}"

  # --lrt "Yes" took too long time for some genes. 20 sec vs 10 min in a small tree.

  hyphy_dnds_json=""
  if [[ -s "hyphy_input.fasta.FITTER.json" ]]; then
    hyphy_dnds_json="hyphy_input.fasta.FITTER.json"
  else
    hyphy_dnds_candidates=()
    mapfile -t hyphy_dnds_candidates < <(find "." -maxdepth 1 -type f -name "hyphy_input.fasta*.json" | sort)
    if [[ ${#hyphy_dnds_candidates[@]} -gt 0 ]]; then
      hyphy_dnds_json="${hyphy_dnds_candidates[0]}"
    fi
  fi
  if [[ -z "${hyphy_dnds_json}" ]]; then
    echo "HyPhy FitMG94 output JSON was not detected. Exiting."
    exit 1
  fi
	mv_out "${hyphy_dnds_json}" "${file_og_hyphy_dnds}"
else
	gg_step_skip "${task}"
fi


# HyPhy RELAX
run_hyphy_relax_for_all_traits() {
  local foreground="$1"
  local out_json="$2"
  local relax_multiple_hits_value=""

  relax_multiple_hits_value="$(detect_hyphy_relax_multiple_hits_off_value)"
  if [[ -n "${relax_multiple_hits_value}" ]]; then
    echo "HyPhy RELAX multiple-hits mode (default-off): ${relax_multiple_hits_value}"
  else
    echo "HyPhy RELAX multiple-hits mode could not be detected. Running without --multiple-hits."
  fi

  binarize_species_trait "${file_sp_trait}" species_trait_binary.tsv
  sed '2,$ s/\t/_.*\t/' species_trait_binary.tsv > foreground.tsv
  IFS=$'\t' read -r -a colname_array < foreground.tsv

  local reversed_mark=""
  if [[ "${foreground}" == "1" ]]; then
    echo "Running HyPhy RELAX for foreground=1"
  elif [[ "${foreground}" == "0" ]]; then
    echo "Running HyPhy RELAX for foreground=0 (reversed)"
    reversed_mark="_reversed"
  else
    echo "Error: foreground must be 1 or 0"
    echo "Exiting."
    exit 1
  fi

  for ((i=1; i<${#colname_array[@]}; i++)); do
    trait="${colname_array[$i]}"
    hyphy_tree_file="hyphy_input_${trait}${reversed_mark}.nwk"
    echo "Processing trait: ${trait}"
    awk -F'\t' -v trait_col="$((i+1))" -v foreground="${foreground}" 'NR>1 && $trait_col == foreground { print $1 }' foreground.tsv > "foreground_${trait}${reversed_mark}.txt"

	    fg_regex=$(paste -sd'|' "foreground_${trait}${reversed_mark}.txt")
    echo "Foreground node search pattern: ${fg_regex}"
    nwkit drop --target intnode --support yes --name yes --infile "${file_og_rooted_tree_analysis}" \
    | nwkit mark --pattern "${fg_regex}" --target "clade" --target_only_clade "yes" --insert_txt "{Foreground}" --outformat 1 \
    | nwkit sanitize --remove_singleton yes --resolve_polytomy no \
    > "${hyphy_tree_file}"

    if grep -q "{Foreground}" "${hyphy_tree_file}"; then
      seqkit seq --threads "${NSLOTS}" "${file_og_trimmed_aln_analysis}" --out-file "hyphy_input_${trait}${reversed_mark}.fasta"
	      hyphy_genetic_code=$(get_hyphy_genetic_code "${genetic_code}")
      relax_multiple_hits_args=()
      if [[ -n "${relax_multiple_hits_value}" ]]; then
        relax_multiple_hits_args=(--multiple-hits "${relax_multiple_hits_value}")
      fi
      hyphy_relax_common_args=(
        --alignment "hyphy_input_${trait}${reversed_mark}.fasta"
        --tree "${hyphy_tree_file}"
        --code "${hyphy_genetic_code}"
        --mode "Classic mode"
        --test "Foreground"
        --models "Minimal"
        --srv "No"
        --rooted "Yes"
        --CPU "${NSLOTS}"
      )
      if ! hyphy relax "${hyphy_relax_common_args[@]}" "${relax_multiple_hits_args[@]}"; then
        if [[ ${#relax_multiple_hits_args[@]} -gt 0 ]]; then
          echo "HyPhy RELAX failed with --multiple-hits=${relax_multiple_hits_value}. Retrying without --multiple-hits."
          hyphy relax "${hyphy_relax_common_args[@]}"
        else
          echo "HyPhy RELAX failed. Exiting."
          exit 1
        fi
      fi
    else
      echo "No foreground lineage is included in this gene tree. Generating an empty JSON file."
      echo "{}" > "hyphy_input_${trait}${reversed_mark}.fasta.RELAX.json"
    fi
  done

  # Combine all hyphy relax output files
  missing_files=()
  relax_output_files=()
  for ((i=1; i<${#colname_array[@]}; i++)); do
    trait="${colname_array[$i]}"
    relax_output_file="hyphy_input_${trait}${reversed_mark}.fasta.RELAX.json"
    if [[ ! -s "${relax_output_file}" ]]; then
      relax_candidates=()
      mapfile -t relax_candidates < <(find "." -maxdepth 1 -type f -name "hyphy_input_${trait}${reversed_mark}.fasta*.json" | sort)
      if [[ ${#relax_candidates[@]} -gt 0 ]]; then
        relax_output_file="${relax_candidates[0]}"
      fi
    fi
    if [[ -s "${relax_output_file}" ]]; then
      relax_output_files+=("${relax_output_file}")
    else
      missing_files+=("${relax_output_file}")
    fi
  done
  if [[ ${#missing_files[@]} -gt 0 ]]; then
    echo "The following HyPhy RELAX output files are missing:"
    for f in "${missing_files[@]}"; do
      echo "${f}"
    done
  else
    echo "All HyPhy RELAX output files are generated. Combining them into a single JSON file: ${out_json}"
    echo "{}" > "combined_relax_output${reversed_mark}.json"
    for ((i=0; i<${#relax_output_files[@]}; i++)); do
      file=${relax_output_files[$i]}
      trait=${colname_array[$((i+1))]}
      jq --arg key "${trait}${reversed_mark}" --slurpfile value "${file}" '. + {($key): $value[0]}' "combined_relax_output${reversed_mark}.json" > "tmp${reversed_mark}.json"
      mv_out "tmp${reversed_mark}.json" "combined_relax_output${reversed_mark}.json"
    done
	    jq . "combined_relax_output${reversed_mark}.json" > "combined_relax_output${reversed_mark}.tmp.json"
	    mv_out "combined_relax_output${reversed_mark}.tmp.json" "${out_json}"
	  fi
}

detect_hyphy_relax_multiple_hits_off_value() {
  local help_text
  help_text="$(hyphy relax --help 2>&1 || true)"
  if [[ "${help_text}" == *"default value: None"* ]]; then
    echo "None"
    return 0
  fi
  if [[ "${help_text}" == *"default value: No"* ]]; then
    echo "No"
    return 0
  fi
  if [[ "${help_text}" == *"Double+Triple"* ]]; then
    echo "None"
    return 0
  fi
  if [[ "${help_text}" == *"multiple-hits"* ]]; then
    echo "None"
    return 0
  fi
  echo ""
}

task="HyPhy RELAX"
disable_if_no_input_file "run_hyphy_relax" "${file_og_rooted_tree_analysis}" "${file_og_trimmed_aln_analysis}" "${file_sp_trait}"
if [[ ! -s "${file_og_hyphy_relax}" && ${run_hyphy_relax} -eq 1 ]]; then
	gg_step_start "${task}"
  run_hyphy_relax_for_all_traits 1 "${file_og_hyphy_relax}"
  if [[ -s "${file_og_hyphy_relax}" ]]; then
    echo "The task has completed successfully: ${task}"
  else
    echo "The task has failed: ${task}"
  fi
else
	gg_step_skip "${task}"
fi

task="HyPhy RELAX with reversed foreground/background"
disable_if_no_input_file "run_hyphy_relax_reversed" "${file_og_rooted_tree_analysis}" "${file_og_trimmed_aln_analysis}" "${file_sp_trait}"
if [[ ! -s "${file_og_hyphy_relax_reversed}" && ${run_hyphy_relax_reversed} -eq 1 ]]; then
  gg_step_start "${task}"
  run_hyphy_relax_for_all_traits 0 "${file_og_hyphy_relax_reversed}"
  if [[ -s "${file_og_hyphy_relax_reversed}" ]]; then
    echo "The task has completed successfully: ${task}"
  else
    echo "The task has failed: ${task}"
  fi
else
  gg_step_skip "${task}"
fi

task="Stochastic character mapping of intron evolution"
disable_if_no_input_file "run_scm_intron" "${file_og_gff_info}" "${file_og_dated_tree_analysis}"
if [[ ! -s "${file_og_scm_intron_summary}" && ${run_scm_intron} -eq 1 ]]; then
	gg_step_start "${task}"


	  Rscript "${dir_script}/scm_intron_evolution.r" \
	  --tree_file="${file_og_dated_tree_analysis}" \
	  --trait_file="${file_og_gff_info}" \
	  --intron_gain_rate="${intron_gain_rate}" \
	  --retrotransposition_rate="${retrotransposition_rate}" \
	  --nrep=1000 \
	  --nslots="${NSLOTS}"

  cp_out intron_evolution_summary.tsv "${file_og_scm_intron_summary}"
  if [[ -e intron_evolution_plot.pdf ]]; then
    cp_out intron_evolution_plot.pdf "${file_og_scm_intron_plot}"
  fi
else
	gg_step_skip "${task}"
fi

task="PhylogeneticEM"
disable_if_no_input_file "run_phylogeneticem" "${file_og_expression}" "${file_og_dated_tree_analysis}"
if [[ ${run_phylogeneticem} -eq 1 && ( ! -s "${file_og_pem_rdata}" || ! -s "${file_og_pem_tree}" || ! -s "${file_og_pem_regime}" || ! -s "${file_og_pem_leaf}" ) && ${num_gene_after_maxalign:-0} -gt 3 ]]; then
	gg_step_start "${task}"

	pem_fit_file=''
	if [[ ${phylogeneticem_use_fit_file} -eq 1 && -s "${file_og_pem_rdata}" ]]; then
		pem_fit_file=${file_og_pem_rdata}
	fi

		Rscript "${dir_script}/detect_OU_shift_PhylogeneticEM.r" \
		--tree_file="${file_og_dated_tree_analysis}" \
		--trait_file="${file_og_expression}" \
		--nslots="${NSLOTS}" \
		--fit_file="${pem_fit_file}" \
		--require_internal_node_labels="${require_internal_node_labels:-1}" \
		--clade_collapse_similarity_method="${clade_collapse_similarity_method}" \
		--clade_collapse_similarity_threshold="${clade_collapse_similarity_threshold}" \
		--ceil_negative=0 \
	  --replicate_sep="_"

	mv_out PhylogeneticEM.tree.tsv "${file_og_pem_tree}"
	mv_out PhylogeneticEM.regime.tsv "${file_og_pem_regime}"
	mv_out PhylogeneticEM.leaf.tsv "${file_og_pem_leaf}"
	mv_out PhylogeneticEM.plot.pdf "${file_og_pem_plot}"
	mv_out PhylogeneticEM.out.RData "${file_og_pem_rdata}"

else
	gg_step_skip "${task}"
fi

task="l1ou"
disable_if_no_input_file "run_l1ou" "${file_og_trimmed_aln_analysis}" "${file_og_expression}" "${file_og_dated_tree_analysis}"
if [[ ( ! -s "${file_og_l1ou_fit_rdata}" || ! -s "${file_og_l1ou_fit_tree}" || ! -s "${file_og_l1ou_fit_regime}" || ! -s "${file_og_l1ou_fit_leaf}" ) && ${run_l1ou} -eq 1 ]]; then
	gg_step_start "${task}"

  num_gene=$(gg_count_fasta_records "${file_og_trimmed_aln_analysis}")
  if [[ ${num_gene} -ge ${large_tree_num_gene} ]]; then
    max_nshift=${large_tree_max_nshift}
  else
    max_nshift=0
  fi

  # taskset is needed because l1ou is not good at handling multiple CPUs.
  if command -v nproc >/dev/null 2>&1; then
    CPU_PER_HOST=$(nproc)
  elif command -v getconf >/dev/null 2>&1; then
    CPU_PER_HOST=$(getconf _NPROCESSORS_ONLN 2>/dev/null || true)
  elif [[ -r /proc/cpuinfo ]]; then
    CPU_PER_HOST=$(awk '/^processor[[:space:]]*:/ {n++} END {print n+0}' /proc/cpuinfo)
  else
    CPU_PER_HOST=1
  fi
  if [[ -z "${CPU_PER_HOST}" || "${CPU_PER_HOST}" -lt 1 ]]; then
    CPU_PER_HOST=1
  fi
  cpu_pick="${NSLOTS}"
  if [[ "${cpu_pick}" -gt "${CPU_PER_HOST}" ]]; then
    cpu_pick="${CPU_PER_HOST}"
  fi
  if [[ "${cpu_pick}" -lt 1 ]]; then
    cpu_pick=1
  fi
		echo "CPU_PER_HOST: ${CPU_PER_HOST}"
			cpu_id=$(python -c 'import sys; from numpy import random; a = random.choice(range(int(sys.argv[2])), int(sys.argv[1]), replace=False); print(",".join([str(b) for b in a]))' "${cpu_pick}" "${CPU_PER_HOST}")
			echo "CPU IDs for l1ou: ${cpu_id}"


	fit_ind_file=''
	if [[ ${l1ou_use_fit_file} -eq 1 && -s "${file_og_l1ou_fit_rdata}" ]]; then
		fit_ind_file=${file_og_l1ou_fit_rdata}
	fi

			l1ou_cmd=(
			  Rscript "${dir_script}/detect_OU_shift_l1ou.r"
			  --max_nshift="${max_nshift}"
			  --tree_file="${file_og_dated_tree_analysis}"
			  --trait_file="${file_og_expression}"
			  --nslots="${NSLOTS}"
			  --require_internal_node_labels="${require_internal_node_labels:-1}"
			  --clade_collapse_similarity_method="${clade_collapse_similarity_method}"
			  --clade_collapse_similarity_threshold="${clade_collapse_similarity_threshold}"
			  --ceil_negative=0
			  --criterion="${l1ou_criterion}"
			  --nbootstrap="${l1ou_nbootstrap}"
			  --fit_ind_file="${fit_ind_file}"
			  --fit_conv_file=''
			  --alpha_upper="${l1ou_alpha_upper}"
			  --detect_convergence="${l1ou_convergence}"
			  --replicate_sep="_"
			)
			if command -v taskset >/dev/null 2>&1 && [[ -n "${cpu_id}" ]]; then
			  taskset -c "${cpu_id}" "${l1ou_cmd[@]}"
			else
			  echo "taskset not available. Running l1ou without CPU pinning."
			  "${l1ou_cmd[@]}"
			fi

	mv_out fit_ind.RData "${file_og_l1ou_fit_rdata}"
	mv_out l1ou_tree.tsv "${file_og_l1ou_fit_tree}"
	mv_out l1ou_regime.tsv "${file_og_l1ou_fit_regime}"
	mv_out l1ou_leaf.tsv "${file_og_l1ou_fit_leaf}"
	mv_out l1ou_plot.pdf "${file_og_l1ou_fit_plot}"
	if [[ ${l1ou_nbootstrap} -gt 0 ]]; then
		cp_out l1ou_bootstrap.tsv "${dir_l1ou_bootstrap}"/"${l1ou_bootstrap}"
	fi
	if [[ ${l1ou_convergence} -eq 1 ]]; then
		cp_out fit_conv.RData "${file_og_l1ou_fit_conv_rdata}"
	fi

else
	gg_step_skip "${task}"
fi

task="Gene tree E-PGLS analysis"
gg_step_skip "${task}"

task="Species tree PGLS analysis"
disable_if_no_input_file "run_pgls_species_tree" "${file_sp_trait}" "${species_tree_pruned}" "${file_og_expression}"
if [[ ! -s "${file_og_species_pgls}" && ${run_pgls_species_tree} -eq 1 ]]; then
	gg_step_start "${task}"
	pgls_merge_replicates="yes"
	if [[ ${pgls_use_phenocov} -eq 1 ]]; then
		pgls_merge_replicates="no"
	fi

		Rscript "${dir_script}/species_tree_pgls.r" \
		--file_sptree="${species_tree_pruned}" \
		--file_exp="${file_og_expression}" \
		--file_trait="${file_sp_trait}" \
		--replicate_sep="_" \
		--exp_value_type="${exp_value_type}" \
		--merge_replicates="${pgls_merge_replicates}" \
		2>&1 | tee pgls.log


	mv_out species_tree_PGLS.tsv "${file_og_species_pgls}"
	mv_out species_tree_PGLS.barplot.pdf "${file_og_species_pgls_plot}"
else
	gg_step_skip "${task}"
fi

task="IQ-TREE ancestral codon sequence reconstruction for CSUBST"
disable_if_no_input_file "run_iqtree_anc" "${file_og_trimmed_aln_analysis}" "${file_og_rooted_tree_analysis}"
if [[ ! -s "${file_og_iqtree_anc}" && ${run_iqtree_anc} -eq 1 ]]; then
		gg_step_start "${task}"

		shopt -s nullglob
		csubst_cleanup_paths=( csubst.* )
		shopt -u nullglob
		if [[ ${#csubst_cleanup_paths[@]} -gt 0 ]]; then
			rm -rf -- "${csubst_cleanup_paths[@]}"
		fi
	  seqkit seq --threads "${NSLOTS}" "${file_og_trimmed_aln_analysis}" --out-file "tmp.csubst.fasta"

  nwkit drop --target intnode --support yes --name no --infile "${file_og_rooted_tree_analysis}" \
	| nwkit intersection --seqin tmp.csubst.fasta --seqout csubst.fasta \
  | nwkit sanitize --remove_singleton yes --resolve_polytomy no \
  > csubst.nwk

	  rm -f -- tmp.csubst.fasta

	build_iqtree_mem_args
	iqtree \
	-s csubst.fasta \
	-te csubst.nwk \
	-m "${codon_model}" \
	-T AUTO \
	--threads-max "${NSLOTS}" \
	--seqtype "CODON${genetic_code}" \
	--prefix csubst \
	--ancestral \
	--rate \
	"${IQTREE_MEM_ARGS[@]}" \
	--seed 12345 \
	--redo

		if [[ -s csubst.rate && -s csubst.state && -s csubst.treefile ]]; then
				if [[ -e "${og_id}.iqtree.anc" ]]; then
					rm -rf -- "${og_id}.iqtree.anc"
				fi
				mkdir -p "${og_id}.iqtree.anc"
				shopt -s nullglob
				csubst_outputs=( csubst.* )
				shopt -u nullglob
				if [[ ${#csubst_outputs[@]} -eq 0 ]]; then
					echo "Expected csubst output files were not found."
					exit 1
				fi
				mv_out "${csubst_outputs[@]}" "${og_id}.iqtree.anc"
				zip -rq "${og_id}.iqtree.anc.zip" "${og_id}.iqtree.anc"
				mv_out "${og_id}.iqtree.anc.zip" "${file_og_iqtree_anc}"
			rm -rf -- "${og_id}.iqtree.anc"
		fi
else
	gg_step_skip "${task}"
fi

task="CSUBST"
disable_if_no_input_file "run_csubst" "${file_og_iqtree_anc}"
if [[ ( ! -s "${file_og_csubst_b}" || ! -s "${file_og_csubst_cb_stats}" ) && ${run_csubst} -eq 1 ]]; then
	gg_step_start "${task}"

	if [[ -s "${file_sp_trait}" ]]; then
		echo "CSUBST foreground specification file: ${file_sp_trait}"
    first_trait_header=""
    IFS= read -r first_trait_header < "${file_sp_trait}" || true
    if [[ "${first_trait_header}" == *" "* ]]; then
      echo "Column names should not contain spaces: ${file_sp_trait}"
      echo "Exiting."
      exit 1
    fi
		sed '2,$ s/\t/_.*\t/' "${file_sp_trait}" > "foreground.tsv"
		foreground_params=( --foreground foreground.tsv --fg_format 2 )
	else
		echo 'Foreground specification file was not found. CSUBST will run without it.'
		foreground_params=()
	fi
  shopt -s nullglob
  csubst_cleanup_paths=( csubst.* )
  shopt -u nullglob
  if [[ ${#csubst_cleanup_paths[@]} -gt 0 ]]; then
    rm -rf -- "${csubst_cleanup_paths[@]}"
  fi
  if [[ -z "${og_id:-}" ]]; then
    echo "og_id is empty. Refusing to run cleanup pattern."
    exit 1
  fi
  shopt -s nullglob
  og_cleanup_paths=( "${og_id}".* )
  shopt -u nullglob
  if [[ ${#og_cleanup_paths[@]} -gt 0 ]]; then
    rm -rf -- "${og_cleanup_paths[@]}"
  fi
		  unzip -q "${file_og_iqtree_anc}"
  csubst_input_base="./${og_id}.iqtree.anc/csubst"

  csubst analyze \
  --genetic_code "${genetic_code}" \
  --infile_type "iqtree" \
  --alignment_file "${csubst_input_base}.fasta" \
  --rooted_tree_file "${csubst_input_base}.nwk" \
  --iqtree_treefile "${csubst_input_base}.treefile" \
  --iqtree_state "${csubst_input_base}.state" \
  --iqtree_rate "${csubst_input_base}.rate" \
  --iqtree_iqtree "${csubst_input_base}.iqtree" \
  --iqtree_log "${csubst_input_base}.log" \
  --iqtree_model "${codon_model}" \
  --iqtree_redo "no" \
  --max_arity "${csubst_max_arity}" \
  --exhaustive_until "${csubst_exhaustive_until}" \
  --cutoff_stat "${csubst_cutoff_stat}" \
  --max_combination "${csubst_max_combination}" \
  --fg_exclude_wg "${csubst_fg_exclude_wg}" \
  --fg_stem_only "${csubst_fg_stem_only}" \
  --mg_sister "no" \
  --exclude_sister_pair "yes" \
  --ml_anc "no" \
  --b "yes" \
  --s "no" \
  --cs "no" \
  --cb "yes" \
  --bs "no" \
  --cbs "no" \
  --calc_quantile "no" \
  --omegaC_method "submodel" \
  --asrv "each" \
  --threads "${NSLOTS}" \
  --calibrate_longtail "yes" \
  --float_type 32 \
  "${foreground_params[@]}"

	  if [[ -s csubst_cb_stats.tsv ]]; then
	    echo "CSUBST was successful."
	    mv_out csubst_b.tsv "${file_og_csubst_b}"
	    mv_out csubst_cb_stats.tsv "${file_og_csubst_cb_stats}"
	    if [[ ${csubst_max_arity} -gt 2 ]]; then
	      for (( i=2; i<=csubst_max_arity; i++ )); do
	        if [[ -e "csubst_cb_${i}.tsv" ]]; then
	          my_csubst_file=file_og_csubst_cb_${i}
	          mv_out csubst_cb_"${i}".tsv "${!my_csubst_file}"
	        fi
	      done
	    fi
  else
    echo "CSUBST failed."
  fi
else
	gg_step_skip "${task}"
fi

task="summary statistics"
if is_output_older_than_inputs "^file_og_" "${file_og_tree_plot}"; then
  summary_flag=0
else
  summary_flag=$?
fi
disable_if_no_input_file "run_summary" "${file_og_rooted_tree_analysis}"
if [[ ( ${summary_flag} -eq 1 || ! -s "${file_og_stat_branch}" || ! -s "${file_og_stat_tree}" ) && ${run_summary} -eq 1 ]]; then
		gg_step_start "${task}"

		if [[ -s "${file_og_notung_reconcil}" ]]; then
			unzip -qf "${file_og_notung_reconcil}"
		fi
  notung_root_log_for_summary="PLACEHOLDER"
  if [[ -d "./${og_id}.notung.root" ]]; then
    notung_log_candidates=()
    mapfile -t notung_log_candidates < <(find "./${og_id}.notung.root" -maxdepth 1 -type f -name "*.ntglog" | sort)
    if [[ ${#notung_log_candidates[@]} -gt 0 ]]; then
      notung_root_log_for_summary="${notung_log_candidates[0]}"
    fi
  fi
  notung_reconcil_stats_for_summary="PLACEHOLDER"
  if [[ -d "./${og_id}.notung_reconcile" ]]; then
    reconcil_stats_candidates=()
    mapfile -t reconcil_stats_candidates < <(find "./${og_id}.notung_reconcile" -maxdepth 1 -type f -name "*.reconciled*.parsable.txt" | sort)
    if [[ ${#reconcil_stats_candidates[@]} -gt 0 ]]; then
      notung_reconcil_stats_for_summary="${reconcil_stats_candidates[0]}"
    fi
  fi
  if [[ ${run_tree_pruning} -eq 1 ]]; then
    generax2orthogroup_statistics="PLACEHOLDER" # generax nhx should be pruned to get used here.
  else
    generax2orthogroup_statistics=${file_og_generax_nhx}
  fi
  summary_unaligned_fasta="PLACEHOLDER"
  summary_trimmed_fasta="PLACEHOLDER"
  summary_promoter_fasta="PLACEHOLDER"
  summary_tmp_files=()
  if [[ -s "${file_og_cds_fasta}" ]]; then
    summary_unaligned_fasta="${og_id}.summary.unaligned.fasta"
    seqkit seq --threads "${NSLOTS}" "${file_og_cds_fasta}" --out-file "${summary_unaligned_fasta}"
    summary_tmp_files+=( "${summary_unaligned_fasta}" )
  fi
  if [[ -s "${file_og_trimmed_aln_analysis}" ]]; then
    summary_trimmed_fasta="${og_id}.summary.trimmed.fasta"
    seqkit seq --threads "${NSLOTS}" "${file_og_trimmed_aln_analysis}" --out-file "${summary_trimmed_fasta}"
    summary_tmp_files+=( "${summary_trimmed_fasta}" )
  fi
  if [[ -s "${file_og_promoter_fasta}" ]]; then
    summary_promoter_fasta="${og_id}.summary.promoter.fasta"
    seqkit seq --threads "${NSLOTS}" "${file_og_promoter_fasta}" --out-file "${summary_promoter_fasta}"
    summary_tmp_files+=( "${summary_promoter_fasta}" )
  fi

		python "${dir_script}/orthogroup_statistics.py" \
		--species_tree "${species_tree_pruned}" \
		--unaligned_aln "${summary_unaligned_fasta}" \
		--trimal_aln "${summary_trimmed_fasta}" \
		--unrooted_tree "${file_og_unrooted_tree_analysis}" \
		--rooted_tree "${file_og_rooted_tree_analysis}" \
		--rooting_log "${file_og_rooted_log}" \
		--notung_root_log "${notung_root_log_for_summary}" \
		--notung_reconcil_stats "${notung_reconcil_stats_for_summary}" \
		--dated_tree "${file_og_dated_tree_analysis}" \
		--dated_log "${file_og_dated_tree_log}" \
		--generax_nhx "${generax2orthogroup_statistics}" \
		--hyphy_dnds_json "${file_og_hyphy_dnds}" \
		--hyphy_relax_json "${file_og_hyphy_relax}" \
		--hyphy_relax_reversed_json "${file_og_hyphy_relax_reversed}" \
		--l1ou_tree "${file_og_l1ou_fit_tree}" \
		--l1ou_regime "${file_og_l1ou_fit_regime}" \
		--l1ou_leaf "${file_og_l1ou_fit_leaf}" \
		--phylogeneticem_tree "${file_og_pem_tree}" \
		--phylogeneticem_regime "${file_og_pem_regime}" \
		--phylogeneticem_leaf "${file_og_pem_leaf}" \
		--expression "${file_og_expression}" \
		--mapdnds_tree_dn "${file_og_mapdnds_dn}" \
		--mapdnds_tree_ds "${file_og_mapdnds_ds}" \
		--codeml_tsv "${file_og_codeml_two_ratio}" \
		--character_gff "${file_og_gff_info}" \
		--fimo "${file_og_fimo}" \
		--promoter_fasta "${summary_promoter_fasta}" \
		--scm_intron "${file_og_scm_intron_summary}" \
		--csubst_b "${file_og_csubst_b}" \
		--gene_pgls_stats "${file_og_gene_pgls}" \
		--species_pgls_stats "${file_og_species_pgls}" \
		--rpsblast "${file_og_rpsblast}" \
		--uniprot "${file_og_uniprot_annotation}" \
		--ncpu "${NSLOTS}" \
		--clade_ortholog_prefix "${treevis_clade_ortholog_prefix}"
  if [[ ${#summary_tmp_files[@]} -gt 0 ]]; then
    rm -f -- "${summary_tmp_files[@]}"
  fi

	#--csubst_cb_stats ${file_og_csubst_cb_stats} \ # Does not support --arity 3 or larger

	cp_out orthogroup.branch.tsv "${file_og_stat_branch}"
	cp_out orthogroup.tree.tsv "${file_og_stat_tree}"

else
	gg_step_skip "${task}"
fi

task="Synteny neighborhood grouping"
if [[ ${treevis_synteny} -eq 1 && ${run_tree_plot} -eq 1 ]]; then
  synteny_needs_update=0
  if [[ ! -s "${file_og_synteny}" || "${file_og_cds_fasta}" -nt "${file_og_synteny}" ]]; then
    synteny_needs_update=1
  fi
  if [[ ${synteny_needs_update} -eq 1 ]]; then
    gg_step_start "${task}"
    if [[ ! -d "${dir_sp_cds}" ]]; then
      echo "species_cds directory not found. Skipping synteny panel input generation: ${dir_sp_cds}"
    elif [[ ! -d "${dir_sp_gff}" ]]; then
      echo "species_gff directory not found. Skipping synteny panel input generation: ${dir_sp_gff}"
    elif [[ ! -s "${file_og_cds_fasta}" ]]; then
      echo "Focal CDS fasta file not found. Skipping synteny panel input generation: ${file_og_cds_fasta}"
    else
      python "${dir_script}/synteny_neighbors.py" \
        --focal_cds_fasta "${file_og_cds_fasta}" \
        --dir_sp_cds "${dir_sp_cds}" \
        --dir_sp_gff "${dir_sp_gff}" \
        --cache_dir "${dir_pg_output}/species_gff_info" \
        --lock_dir "${file_og_parameters_dir}/synteny_locks" \
        --gff2genestat_script "${dir_script}/gff2genestat.py" \
        --window "${treevis_synteny_window}" \
        --evalue "${query_blast_evalue}" \
        --genetic_code "${genetic_code}" \
        --threads "${NSLOTS}" \
        --outfile "${file_og_synteny}"
      if [[ ! -s "${file_og_synteny}" ]]; then
        echo "No synteny links were generated: ${file_og_synteny}"
      fi
    fi
  else
    gg_step_skip "${task}"
  fi
else
  gg_step_skip "${task}"
fi

if [[ ${treevis_synteny} -eq 1 && -s "${file_og_synteny}" ]]; then
  if [[ ! -s "${file_og_tree_plot}" || "${file_og_synteny}" -nt "${file_og_tree_plot}" ]]; then
    summary_flag=1
  fi
fi

task="stat_branch2tree_plot"
disable_if_no_input_file "run_tree_plot" "${file_og_stat_branch}" "${file_og_stat_tree}"
if [[ ${run_tree_plot} -eq 1 ]]; then
  if ! Rscript -e "if (!requireNamespace('ggimage', quietly=TRUE)) quit(status=1)" >/dev/null 2>&1; then
    echo "ggimage package is unavailable. Disabling run_tree_plot."
    run_tree_plot=0
  fi
fi
if ( [[ ${summary_flag} -eq 1 || ! -s "${file_og_tree_plot}" ]] ) && [[ ${run_tree_plot} -eq 1 ]]; then
    gg_step_start "${task}"

    num_tip_treeplot=$(
      awk -F $'\t' '
        NR==1 {
          col=0
          for (i=1; i<=NF; i++) {
            if ($i=="so_event") {
              col=i
              break
            }
          }
          next
        }
        (col>0 && $col=="L") {n++}
        END {print n+0}
      ' "${file_og_stat_branch}"
    )
    panel11_trimmed_aln="${file_og_trimmed_aln_analysis}"
    panel11_trimmed_n=$(gg_count_fasta_records "${panel11_trimmed_aln}")
    if [[ ${panel11_trimmed_n} -lt ${num_tip_treeplot} ]]; then
      for candidate in \
        "${file_og_clipkit}" \
        "${file_og_orthogroup_extraction_fasta}" \
        "${file_og_maxalign}" \
        "${file_og_mafft}" \
        "${file_og_cds_fasta}"; do
        candidate_n=$(gg_count_fasta_records "${candidate}")
        if [[ ${candidate_n} -ge ${num_tip_treeplot} ]]; then
          panel11_trimmed_aln="${candidate}"
          panel11_trimmed_n=${candidate_n}
          break
        fi
      done
    fi
    panel11_untrimmed_aln="${file_og_untrimmed_aln_analysis}"
    panel11_untrimmed_n=$(gg_count_fasta_records "${panel11_untrimmed_aln}")
    if [[ ${panel11_untrimmed_n} -lt ${num_tip_treeplot} ]]; then
      for candidate in \
        "${file_og_orthogroup_extraction_fasta}" \
        "${file_og_mafft}" \
        "${file_og_cds_fasta}"; do
        candidate_n=$(gg_count_fasta_records "${candidate}")
        if [[ ${candidate_n} -ge ${num_tip_treeplot} ]]; then
          panel11_untrimmed_aln="${candidate}"
          panel11_untrimmed_n=${candidate_n}
          break
        fi
      done
    fi
    echo "Tree plot alignment inputs: tips=${num_tip_treeplot}, trimmed=${panel11_trimmed_n} (${panel11_trimmed_aln}), untrimmed=${panel11_untrimmed_n} (${panel11_untrimmed_aln})"

	    if [[ ${treevis_clade_ortholog} -eq 1 ]]; then
	        ortholog_prefix=${treevis_clade_ortholog_prefix}
	    else
	        ortholog_prefix=""
	    fi
    cb_path=${file_og_csubst_cb_2/cb_2/cb_ARITY}

    Rscript "${dir_script}/stat_branch2tree_plot.r" \
    --stat_branch="${file_og_stat_branch}" \
    --tree_annotation_dir="${dir_script}/tree_annotation" \
	    --max_delta_intron_present="${treevis_retrotransposition_delta_intron}" \
    --width="7.2" \
    --rel_widths="" \
    --panel1="tree,${treevis_branch_length},${treevis_support_value},${treevis_branch_color},L" \
    --panel2="heatmap,${treevis_heatmap_transform},abs,_,expression_" \
    --panel3="pointplot,no,rel,_,expression_" \
    --panel4="cluster_membership,${treevis_max_intergenic_dist}" \
	    --panel5="synteny,${file_og_synteny},${treevis_synteny_window}" \
	    --panel6="tiplabel" \
	    --panel7="signal_peptide" \
	    --panel8="transmembrane_domain" \
    --panel9="intron_number" \
    --panel10="domain,${file_og_rpsblast}" \
	    --panel11="alignment,${panel11_trimmed_aln},${panel11_untrimmed_aln}" \
	    --panel12="fimo,${promoter_bp},${fimo_qvalue}" \
	    --panel13="meme,${file_og_meme}" \
	    --panel14="ortholog,${ortholog_prefix},${file_og_dated_tree}" \
	    --show_branch_id="yes" \
		    --event_method="${treevis_event_method}" \
		    --species_color_table="PLACEHOLDER" \
	    --pie_chart_value_transformation="${treevis_pie_chart_value_transformation}" \
	    --long_branch_display="${treevis_long_branch_display}" \
	    --long_branch_ref_quantile="${treevis_long_branch_ref_quantile}" \
	    --long_branch_detect_ratio="${treevis_long_branch_detect_ratio}" \
	    --long_branch_cap_ratio="${treevis_long_branch_cap_ratio}" \
	    --long_branch_tail_shrink="${treevis_long_branch_tail_shrink}" \
	    --long_branch_max_fraction="${treevis_long_branch_max_fraction}" \
    --protein_convergence="100,100,yes,3-${csubst_max_arity},${cb_path},${csubst_cutoff_stat}"

	    if [[ -e "df_fimo.tsv" ]]; then
	    	mv_out "df_fimo.tsv" "${file_og_fimo_collapsed}"
	    fi
	    mv_out stat_branch2tree_plot.pdf "${file_og_tree_plot}"
else
	gg_step_skip "${task}"
fi

# Copy parameter files and codes to ${file_og_parameters_dir} for record
mkdir -p "${file_og_parameters_dir}"
file_params=(
	"${file_sp_trait}"
	"${species_tree}"
	"${species_tree_pruned}"
)
	for file_from in "${file_params[@]}"; do
		file_to="${file_og_parameters_dir}/$(basename "${file_from}")"
		if [[ ! -e "${file_from}" ]]; then
			continue
		fi
		lock_file="${file_to}.lock"
			if command -v flock >/dev/null 2>&1; then
				exec 7> "${lock_file}"
				flock 7
				filesize_from=$(stat -c%s "${file_from}")
				filesize_to=0
				if [[ -s "${file_to}" ]]; then
				filesize_to=$(stat -c%s "${file_to}")
		fi
		if [[ ${filesize_from} -ne ${filesize_to} ]]; then
			echo "Storing important files for record: ${file_to}"
				cp_out "${file_from}" "${file_to}"
			fi
				flock -u 7
				exec 7>&-
			else
				echo "Error: flock command is required but not available."
				exit 1
			fi
		done

cd "${dir_pg}"
remove_empty_subdirs "${dir_output_active}"

if [[ -s "${file_og_stat_branch}" && -s "${file_og_stat_tree}" && -s "${file_og_tree_plot}" && ${gg_debug_mode:-0} -eq 0 ]]; then
    echo "Output files detected."
    echo "${file_og_stat_branch}"
    echo "${file_og_stat_tree}"
    echo "${file_og_tree_plot}"
    echo "$(date): Exiting Singularity environment"
    exit 8
elif [[ -s "${file_og_stat_branch}" && -s "${file_og_stat_tree}" && -s "${file_og_tree_plot}" && ${gg_debug_mode:-0} -eq 1 ]]; then
    echo "Output files detected & debug mode."
else
    echo "Output files not found."
fi

###################
echo "$(date): Exiting Singularity environment"
