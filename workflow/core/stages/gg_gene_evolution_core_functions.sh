# Named stage functions for gg_gene_evolution_core.sh.
# This file is sourced by workflow/core/gg_gene_evolution_core.sh.

build_iqtree_mem_args() {
  IQTREE_MEM_ARGS=()
  if [[ -n "${GG_MEM_TOOL_GB:-}" ]]; then
    IQTREE_MEM_ARGS=(--mem "${GG_MEM_TOOL_GB}G")
  fi
}

resolve_query_blast_evalue() {
  local requested_evalue="$1"
  local query_fasta="$2"
  local cutoffs="$3"
  local requested_lower=""
  local query_blast_stats=""

  requested_lower=$(printf '%s' "${requested_evalue}" | tr '[:upper:]' '[:lower:]')
  effective_query_blast_evalue="${requested_evalue}"
  query_blast_query_num_seqs="NA"
  query_blast_query_min_aa_len="NA"
  query_blast_query_avg_aa_len="NA"
  query_blast_query_max_aa_len="NA"

  if [[ "${requested_lower}" != "auto" ]]; then
    return 0
  fi
  if [[ ! -s "${query_fasta}" ]]; then
    echo "query_blast_evalue=auto requires a non-empty query FASTA: ${query_fasta}"
    return 1
  fi
  if ! query_blast_stats=$(seqkit stats --tabular "${query_fasta}" | awk 'NR==2 {gsub(/,/, "", $4); gsub(/,/, "", $6); gsub(/,/, "", $7); gsub(/,/, "", $8); print $4, $6, $7, $8}'); then
    echo "Failed to calculate query FASTA length statistics for query_blast_evalue=auto: ${query_fasta}"
    return 1
  fi
  read -r query_blast_query_num_seqs query_blast_query_min_aa_len query_blast_query_avg_aa_len query_blast_query_max_aa_len <<< "${query_blast_stats}"
  if [[ ! "${query_blast_query_max_aa_len}" =~ ^[0-9]+$ ]]; then
    echo "Failed to parse max query amino-acid length for query_blast_evalue=auto: ${query_fasta}"
    return 1
  fi
  if ! effective_query_blast_evalue=$(awk -v max_len="${query_blast_query_max_aa_len}" -v cutoffs="${cutoffs}" '
    function trim(x) {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", x)
      return x
    }
    BEGIN {
      n = split(cutoffs, rules, ",")
      for (i = 1; i <= n; i++) {
        rule = trim(rules[i])
        if (rule == "") {
          continue
        }
        n_parts = split(rule, parts, ":")
        if (n_parts != 2) {
          exit 2
        }
        cutoff = trim(parts[1])
        evalue = trim(parts[2])
        if (evalue == "") {
          exit 2
        }
        if (tolower(cutoff) == "inf") {
          print evalue
          exit 0
        }
        if (cutoff !~ /^[0-9]+$/) {
          exit 2
        }
        if (max_len <= cutoff + 0) {
          print evalue
          exit 0
        }
      }
      exit 3
    }
  '); then
    echo "Invalid query_blast_auto_evalue_maxlen_cutoffs: ${cutoffs}"
    return 1
  fi
  if [[ -z "${effective_query_blast_evalue}" ]]; then
    echo "Failed to resolve query_blast_evalue=auto using cutoffs: ${cutoffs}"
    return 1
  fi
}

prepare_synteny_evalue_fasta() {
  local outfile="$1"

  if [[ -s "${file_og_query_aa_fasta}" ]]; then
    seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_query_aa_fasta}" --out-file "${outfile}"
    return
  fi
  if [[ "${input_sequence_mode}" == "protein" ]]; then
    seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_primary_fasta}" --out-file "${outfile}"
    return
  fi
  seqkit translate \
    --allow-unknown-codon \
    --transl-table "${genetic_code}" \
    --threads "${GG_TASK_CPUS}" \
    "${file_og_primary_fasta}" \
    --out-file "${outfile}"
}

binarize_species_trait() {
  local file_in="$1"
  local file_out="$2"
  python - "${file_in}" "${file_out}" << 'PY'
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

  cat > "${file_param}" << EOF
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
  python - "${file_tree}" "${context}" << 'PY'
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

set_profile_default_override() {
  local var_name=$1
  local default_value=$2
  local new_value=$3
  local current_value="${!var_name:-}"
  if [[ "${current_value}" == "${default_value}" ]]; then
    printf -v "${var_name}" '%s' "${new_value}"
    echo "gene_evolution_profile=${gene_evolution_profile}: ${var_name}=${new_value}"
  fi
}

apply_gene_evolution_profile() {
  case "${gene_evolution_profile}" in
    ""|default)
      gene_evolution_profile="default"
      ;;
    hgt)
      mode_gene_evolution="orthogroup"
      echo "gene_evolution_profile=hgt: mode_gene_evolution=${mode_gene_evolution}"
      set_profile_default_override run_uniprot_annotation "0" "1"
      set_profile_default_override run_generax "0" "1"
      set_profile_default_override generax_rec_model "UndatedDL" "UndatedDTL"
      set_profile_default_override run_collect_gff_info "0" "1"
      set_profile_default_override run_generate_expression_matrix "0" "1"
      set_profile_default_override run_scm_intron "0" "1"
      set_profile_default_override treevis_event_method "species_overlap" "auto"
      ;;
    *)
      echo "Invalid gene_evolution_profile: ${gene_evolution_profile}"
      echo 'gene_evolution_profile must be either "default" or "hgt". Exiting.'
      exit 1
      ;;
  esac
}

gene_evolution_model_is_aa() {
  local model_name=${1:-}
  local base_model=${model_name%%+*}
  case "${base_model}" in
    Blosum62|cpREV|Dayhoff|DCMut|DEN|FLU|HIVb|HIVw|JTT|JTT-DCMut|LG|mtART|mtMAM|mtREV|mtZOA|PMB|rtREV|stmtREV|VT|WAG|LG4M|LG4X|PROTGTR)
      return 0
      ;;
  esac
  return 1
}

assert_gene_evolution_aa_model_for_protein_mode() {
  local task_name=${1:-task}
  if [[ "${input_sequence_mode}" != "protein" ]]; then
    return 0
  fi
  if ! gene_evolution_model_is_aa "${generax_model}"; then
    echo "input_sequence_mode=protein requires an amino-acid substitution model for ${task_name}: ${generax_model}"
    echo 'Set generax_model to an amino-acid model such as LG+G4.'
    exit 1
  fi
}

apply_gene_evolution_input_sequence_mode() {
  if [[ "${input_sequence_mode}" != "protein" ]]; then
    return 0
  fi

  echo "input_sequence_mode=protein: gg_gene_evolution will keep only protein-compatible or output-only stages."

  if [[ "${mode_gene_evolution}" == "query2family" ]]; then
    echo 'input_sequence_mode=protein is not supported with mode_gene_evolution=query2family.'
    echo 'query2family currently requires species_cds-backed search and CDS extraction.'
    exit 1
  fi

  if [[ "${generax_model}" == "GTR+G4" ]]; then
    generax_model="LG+G4"
    echo "input_sequence_mode=protein: generax_model=${generax_model}"
  fi

  disable_flag_with_reason "run_maxalign" "input_sequence_mode=protein: the current MaxAlign path expects codon alignments."
  disable_flag_with_reason "run_mapdnds_parameter_estimation" "input_sequence_mode=protein: mapdNdS parameter estimation requires codon alignments."
  disable_flag_with_reason "run_mapdnds" "input_sequence_mode=protein: mapdNdS requires codon alignments."
  disable_flag_with_reason "run_codeml_two_ratio" "input_sequence_mode=protein: codeml two-ratio requires codon alignments."
  disable_flag_with_reason "run_hyphy_dnds" "input_sequence_mode=protein: HyPhy dN/dS requires codon alignments."
  disable_flag_with_reason "run_hyphy_relax" "input_sequence_mode=protein: HyPhy RELAX requires codon alignments."
  disable_flag_with_reason "run_hyphy_relax_reversed" "input_sequence_mode=protein: reversed HyPhy RELAX requires codon alignments."
  disable_flag_with_reason "run_iqtree_anc" "input_sequence_mode=protein: ancestral reconstruction for CSUBST runs in codon mode."
  disable_flag_with_reason "run_csubst" "input_sequence_mode=protein: CSUBST currently depends on codon-mode ancestral reconstruction."
  disable_flag_with_reason "run_csubst_scan" "input_sequence_mode=protein: CSUBST scan currently depends on codon-mode ancestral reconstruction."
}

write_species_trait_foreground_regex_table() {
  local input_tsv=$1
  local output_tsv=$2
  PYTHONPATH="${gg_support_dir}${PYTHONPATH:+:${PYTHONPATH}}" python - "${input_tsv}" "${output_tsv}" <<'PY'
import csv
import re
import sys

from species_labeling import extract_species_label

RANK_OR_QUALIFIER_TOKENS = (
    "cf", "aff", "nr", "sp", "spp",
    "subsp", "ssp", "subspecies", "var", "variety", "forma", "form", "f",
    "strain", "substrain", "serovar", "serotype", "serogroup",
    "pathovar", "pv", "biovar", "biotype", "chemovar", "morphovar",
    "cultivar", "cv", "isolate", "group", "subgroup", "complex",
    "clade", "lineage", "section", "series", "ecotype", "breed",
)


def species_foreground_regex(value):
    species_label = extract_species_label(value) or str(value or "").strip()
    if species_label == "":
        return ""
    token_pattern = "|".join(re.escape(token) for token in RANK_OR_QUALIFIER_TOKENS)
    return r"^{}_(?!(?:{})(?:\.|_)).*".format(re.escape(species_label), token_pattern)


infile, outfile = sys.argv[1:3]
with open(infile, "r", encoding="utf-8", errors="replace", newline="") as src, open(
    outfile, "w", encoding="utf-8", newline=""
) as dst:
    reader = csv.reader(src, delimiter="\t")
    writer = csv.writer(dst, delimiter="\t", lineterminator="\n")
    for row_index, row in enumerate(reader):
        if row_index > 0 and row:
            row[0] = species_foreground_regex(row[0])
        writer.writerow(row)
PY
}

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

species_protein_input_has_files() {
  local protein_files=()
  mapfile -t protein_files < <(gg_find_fasta_files "${dir_sp_protein_input}" 1)
  [[ ${#protein_files[@]} -gt 0 ]]
}

translate_orthogroup_cds_to_protein_fasta() {
  local cds_fasta=$1
  local protein_out=$2
  local table_path=$3
  local translated_tmp="${og_id}.translated.pep.tmp.fasta"
  local species_code=""
  local sp_ub=""
  local num_cds=0
  local num_protein=0
  local -a species_names=()

  rm -f -- "${translated_tmp}"
  touch "${translated_tmp}"
  mapfile -t species_names < <(
    seqkit seq --threads "${GG_TASK_CPUS}" --name "${cds_fasta}" |
      while IFS= read -r header; do
        gg_species_name_from_path_or_dot "${header}"
      done | sed -e '/^$/d' | sort -u
  )
  if [[ ${#species_names[@]} -eq 0 ]]; then
    echo "No species prefixes were detected in the focal CDS FASTA: ${cds_fasta}"
    exit 1
  fi

  for sp_ub in "${species_names[@]}"; do
    species_code=$(gg_lookup_species_genetic_code "${sp_ub}" "${table_path}" "${genetic_code}")
    echo "Translation started: ${sp_ub} (genetic_code=${species_code}) -> $(basename "${protein_out}")"
    seqkit grep --threads "${GG_TASK_CPUS}" --use-regexp --pattern "^${sp_ub}([_.-])" "${cds_fasta}" |
      seqkit seq --remove-gaps --threads "${GG_TASK_CPUS}" |
      gg_prepare_cds_fasta_stream "${GG_TASK_CPUS}" "${species_code}" |
      seqkit translate --allow-unknown-codon --transl-table "${species_code}" --threads "${GG_TASK_CPUS}" |
      sed -e '/^1 1$/d' -e 's/_frame=1[[:space:]]*//' \
        >> "${translated_tmp}"
  done
  num_cds=$(gg_count_fasta_records "${cds_fasta}")
  num_protein=$(gg_count_fasta_records "${translated_tmp}")
  if [[ ${num_cds} -ne ${num_protein} ]]; then
    echo "Protein translation produced a different number of sequences (${num_protein}) than the source CDS FASTA (${num_cds})."
    echo "Exiting."
    exit 1
  fi
  seqkit seq --threads "${GG_TASK_CPUS}" "${translated_tmp}" --out-file "${protein_out}"
  rm -f -- "${translated_tmp}"
}

prepare_species_tree_pruned() {
  local task_local="Species tree pruning"
  local species_sequence_dir="${dir_sp_cds}"
  if [[ ! -s "${species_tree}" ]]; then
    echo "$(date): Warning: ${task_local}: source species tree was not found."
    echo "Missing: ${species_tree}"
    return 1
  fi

  if [[ -s "${species_tree_pruned}" ]]; then
    return 0
  fi

  ensure_parent_dir "${species_tree_pruned}"

  if [[ "${input_sequence_mode}" == "protein" ]] && species_protein_input_has_files; then
    species_sequence_dir="${dir_sp_protein_input}"
  fi

  local sequence_files=()
  mapfile -t sequence_files < <(gg_find_fasta_files "${species_sequence_dir}" 1)
  if [[ ${#sequence_files[@]} -eq 0 ]]; then
    echo "$(date): ${task_local}: no species sequence files detected in ${species_sequence_dir}. Copying source species tree as-is."
    cp_out "${species_tree}" "${species_tree_pruned}"
    return 0
  fi

  local cds_spp=()
  local sequence_file
  for sequence_file in "${sequence_files[@]}"; do
    cds_spp+=("$(gg_species_name_from_path "${sequence_file}")")
  done
  mapfile -t cds_spp < <(printf '%s\n' "${cds_spp[@]}" | sed -e '/^[[:space:]]*$/d' | sort -u)
  if [[ ${#cds_spp[@]} -eq 0 ]]; then
    echo "$(date): ${task_local}: species names could not be parsed from species CDS files. Copying source tree."
    cp_out "${species_tree}" "${species_tree_pruned}"
    return 0
  fi

  local keep_pattern
  keep_pattern=$(
    printf '%s\n' "${cds_spp[@]}" |
      sed -e 's/[][(){}.^$+*?|\\-]/\\&/g' |
      paste -sd'|' -
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

cleanup_tmp_dir_on_normal_exit() {
  local exit_code=$?
  if [[ ${delete_tmp_dir} -eq 1 && (${exit_code} -eq 0 || ${exit_code} -eq 8) ]]; then
    if [[ -n "${dir_tmp:-}" && -d "${dir_tmp}" && "${dir_tmp}" != "/" ]]; then
      echo "Deleting ${dir_tmp}"
      rm -rf -- "${dir_tmp}"
    elif [[ -n "${dir_tmp:-}" ]]; then
      echo "Refusing to delete unsafe tmp directory: ${dir_tmp}"
    fi
  fi
  return ${exit_code}
}

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
  write_species_trait_foreground_regex_table species_trait_binary.tsv foreground.tsv
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

  for ((i = 1; i < ${#colname_array[@]}; i++)); do
    trait="${colname_array[$i]}"
    hyphy_tree_file="hyphy_input_${trait}${reversed_mark}.nwk"
    echo "Processing trait: ${trait}"
    awk -F'\t' -v trait_col="$((i + 1))" -v foreground="${foreground}" 'NR>1 && $trait_col == foreground { print $1 }' foreground.tsv > "foreground_${trait}${reversed_mark}.txt"

    fg_regex=$(paste -sd'|' "foreground_${trait}${reversed_mark}.txt")
    echo "Foreground node search pattern: ${fg_regex}"
    nwkit drop --target intnode --support yes --name yes --infile "${file_og_rooted_tree_analysis}" |
      nwkit mark --pattern "${fg_regex}" --target "clade" --target_only_clade "yes" --insert_txt "{Foreground}" --outformat 1 |
      nwkit sanitize --remove_singleton yes --resolve_polytomy no \
        > "${hyphy_tree_file}"

    if grep -q "{Foreground}" "${hyphy_tree_file}"; then
      seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file "hyphy_input_${trait}${reversed_mark}.fasta"
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
        --CPU "${GG_TASK_CPUS}"
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
  for ((i = 1; i < ${#colname_array[@]}; i++)); do
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
    for ((i = 0; i < ${#relax_output_files[@]}; i++)); do
      file=${relax_output_files[$i]}
      trait=${colname_array[$((i + 1))]}
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
