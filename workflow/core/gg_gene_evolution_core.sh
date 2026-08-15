#!/usr/bin/env bash
set -euo pipefail

gg_core_bootstrap="/script/support/gg_core_bootstrap.sh"
if [[ ! -s "${gg_core_bootstrap}" ]]; then
  gg_core_bootstrap="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)/../support/gg_core_bootstrap.sh"
fi
# shellcheck disable=SC1090
source "${gg_core_bootstrap}"
unset gg_core_bootstrap
gg_source_common_params_from_core "${BASH_SOURCE[0]:-$0}"

### Start: Job-supplied configuration ###
# Configuration variables are provided by gg_gene_evolution_entrypoint.sh.
busco_lineage="${busco_lineage:-${GG_COMMON_BUSCO_LINEAGE:-auto}}"
genetic_code="${genetic_code:-${GG_COMMON_GENETIC_CODE:-1}}"
annotation_species="${annotation_species:-${GG_COMMON_REFERENCE_SPECIES:-auto}}"
input_sequence_mode="${input_sequence_mode:-${GG_COMMON_INPUT_SEQUENCE_MODE:-cds}}"
species_label_parser="${species_label_parser:-${GG_COMMON_SPECIES_LABEL_PARSER:-taxonomic}}"
species_label_regex="${species_label_regex:-${GG_COMMON_SPECIES_LABEL_REGEX:-}}"
species_label_map_tsv="${species_label_map_tsv:-${GG_COMMON_SPECIES_LABEL_MAP_TSV:-}}"
cdskit_localize_model="${cdskit_localize_model:-targeting5-perox-deeploc21-et-v1}"
cdskit_localize_organism_group="${cdskit_localize_organism_group:-auto}"
cdskit_localize_include_features="${cdskit_localize_include_features:-0}"
cdskit_localize_no_model_download="${cdskit_localize_no_model_download:-0}"
run_extract_query_fasta="${run_extract_query_fasta:-1}"
run_extract_primary_fasta="${run_extract_primary_fasta:-1}"
run_generate_expression_matrix="${run_generate_expression_matrix:-0}"
run_collect_gff_info="${run_collect_gff_info:-0}"
run_extract_promoter_fasta="${run_extract_promoter_fasta:-0}"
run_expression_trait_pgls="${run_expression_trait_pgls:-0}"
treevis_query_marker="${treevis_query_marker:-1}"
query_blast_evalue="${query_blast_evalue:-auto}"
query_blast_auto_evalue_maxlen_cutoffs="${query_blast_auto_evalue_maxlen_cutoffs:-40:1000,80:100,150:10,300:1,inf:0.01}"

# Unified RSC/species-tree expression-trait PGLS. Defaults are repeated here so
# direct core-script invocations remain safe under set -u.
pgls_methods="${pgls_methods:-rsc}"
species_expression_aggregation="${species_expression_aggregation:-sum}"
species_paralog_missing="${species_paralog_missing:-error}"
species_paralog_sampling_covariance="${species_paralog_sampling_covariance:-}"
rphylopars_sampling_covariance="${rphylopars_sampling_covariance:-require-diagonal}"
rsc_responses="${rsc_responses:-all}"
rsc_predictors="${rsc_predictors:-all}"
rsc_predictor_mode="${rsc_predictor_mode:-separate}"
rsc_event_source="${rsc_event_source:-auto}"
rsc_speciation_coverage="${rsc_speciation_coverage:-complete}"
rsc_event_weighting="${rsc_event_weighting:-equal}"
rsc_model="${rsc_model:-hierarchical}"
rsc_gene_branch_length="${rsc_gene_branch_length:-original}"
rsc_gene_evolution_model="${rsc_gene_evolution_model:-brownian}"
rsc_gene_evolution_parameter="${rsc_gene_evolution_parameter:-auto}"
rsc_species_branch_length="${rsc_species_branch_length:-original}"
rsc_species_evolution_model="${rsc_species_evolution_model:-brownian}"
rsc_species_evolution_parameter="${rsc_species_evolution_parameter:-auto}"
rsc_inference="${rsc_inference:-wald}"
rsc_bootstrap_replicates="${rsc_bootstrap_replicates:-1000}"
rsc_seed="${rsc_seed:-1}"
rsc_confidence_level="${rsc_confidence_level:-0.95}"
rsc_reml="${rsc_reml:-yes}"
rsc_min_species_events="${rsc_min_species_events:-2}"
rsc_unmatched="${rsc_unmatched:-error}"
rsc_replicate_separator="${rsc_replicate_separator:-_}"
rsc_expression_sample_metadata="${rsc_expression_sample_metadata:-}"
rsc_within_variance="${rsc_within_variance:-pooled}"
rsc_technical_aggregation="${rsc_technical_aggregation:-error}"
rsc_predictor_biological_id="${rsc_predictor_biological_id:-}"
rsc_predictor_technical_id="${rsc_predictor_technical_id:-}"
rsc_predictor_batch="${rsc_predictor_batch:-}"
rsc_predictor_within_variance="${rsc_predictor_within_variance:-pooled}"
rsc_predictor_technical_aggregation="${rsc_predictor_technical_aggregation:-error}"
rsc_predictor_standard_error_columns="${rsc_predictor_standard_error_columns:-}"
rsc_predictor_sample_size_columns="${rsc_predictor_sample_size_columns:-}"
rsc_categorical_predictors="${rsc_categorical_predictors:-}"
rsc_ordered_predictors="${rsc_ordered_predictors:-}"
rsc_factor_reference="${rsc_factor_reference:-}"
rsc_factor_coding="${rsc_factor_coding:-treatment}"
rsc_categorical_replicate_policy="${rsc_categorical_replicate_policy:-latent}"
rsc_event_random_effect="${rsc_event_random_effect:-auto}"
rsc_lineage_random_slope="${rsc_lineage_random_slope:-auto}"
rsc_lineage_inference="${rsc_lineage_inference:-none}"
rsc_lineage_leave_one_out="${rsc_lineage_leave_one_out:-no}"
rsc_categorical_origin_diagnostics="${rsc_categorical_origin_diagnostics:-none}"
rsc_origin_map_replicates="${rsc_origin_map_replicates:-200}"
rsc_origin_map_threads="${rsc_origin_map_threads:-1}"
rsc_origin_min_posterior="${rsc_origin_min_posterior:-0.5}"
rsc_origin_leave_one_out="${rsc_origin_leave_one_out:-no}"
rsc_allow_large_dense="${rsc_allow_large_dense:-no}"

# Substitution model in CSUBST and mapdNdS
if [[ ${genetic_code} -eq 1 ]]; then
  codon_model="ECMK07+F+R4" # Standard genetic code
else
  codon_model="GY+F+R4" # Non-standard genetic code
fi
### End: Job-supplied configuration ###

### ----------------------------------------------------------------------- ###

### Modify below if you need to add a new analysis or need to fix some bugs ###

gg_bootstrap_core_runtime "${BASH_SOURCE[0]:-$0}" "base" 1 1
gg_extract_expected_zip_prefix() {
  local archive_path=$1
  local expected_prefix=$2
  python "${gg_support_dir}/safe_zip_extract.py" \
    --archive "${archive_path}" \
    --destination-root . \
    --expected-prefix "${expected_prefix}" >/dev/null
}

delete_preexisting_tmp_dir=${delete_preexisting_tmp_dir:-1}
delete_tmp_dir=${delete_tmp_dir:-1}
gene_family_run_token=""
gene_family_state_finalized=1
gene_family_run_succeeded=0
gene_family_materialization_receipt=""
gene_family_run_lock_path=""
gene_family_run_lock_heartbeat_pid=""
gene_family_output_storage=$(printf '%s' "${gene_family_output_storage:-${GG_COMMON_GENE_FAMILY_OUTPUT_STORAGE:-zip}}" | tr '[:upper:]' '[:lower:]')
gene_family_zip_min_batch_files="${gene_family_zip_min_batch_files:-${GG_COMMON_GENE_FAMILY_ZIP_MIN_BATCH_FILES:-100}}"
gene_family_zip_compression=$(printf '%s' "${gene_family_zip_compression:-${GG_COMMON_GENE_FAMILY_ZIP_COMPRESSION:-adaptive}}" | tr '[:upper:]' '[:lower:]')
gene_family_zip_compression_level="${gene_family_zip_compression_level:-${GG_COMMON_GENE_FAMILY_ZIP_COMPRESSION_LEVEL:-6}}"
gene_family_zip_workers="${gene_family_zip_workers:-${GG_COMMON_GENE_FAMILY_ZIP_WORKERS:-1}}"
gene_family_final_zip_max_bytes="${gene_family_final_zip_max_bytes:-${GG_COMMON_GENE_FAMILY_FINAL_ZIP_MAX_BYTES:-0}}"
gene_family_tmp_retention_days="${gene_family_tmp_retention_days:-${GG_COMMON_GENE_FAMILY_TMP_RETENTION_DAYS:-7}}"
gene_family_tmp_max_dirs="${gene_family_tmp_max_dirs:-${GG_COMMON_GENE_FAMILY_TMP_MAX_DIRS:-100}}"
gene_family_tmp_max_bytes="${gene_family_tmp_max_bytes:-${GG_COMMON_GENE_FAMILY_TMP_MAX_BYTES:-107374182400}}"
gene_family_tmp_max_files="${gene_family_tmp_max_files:-${GG_COMMON_GENE_FAMILY_TMP_MAX_FILES:-100000}}"
case "${gene_family_output_storage}" in
  zip|files)
    ;;
  raw)
    gene_family_output_storage="files"
    ;;
  *)
    echo "Invalid gene_family_output_storage: ${gene_family_output_storage}; expected zip, files, or raw." >&2
    exit 1
    ;;
esac
if [[ ! "${gene_family_zip_min_batch_files}" =~ ^[0-9]+$ || ${gene_family_zip_min_batch_files} -lt 1 ]]; then
  echo "Invalid gene_family_zip_min_batch_files: ${gene_family_zip_min_batch_files}; expected a positive integer." >&2
  exit 1
fi
case "${gene_family_zip_compression}" in
  adaptive|deflate|store)
    ;;
  *)
    echo "Invalid gene_family_zip_compression: ${gene_family_zip_compression}; expected adaptive, deflate, or store." >&2
    exit 1
    ;;
esac
if [[ ! "${gene_family_zip_compression_level}" =~ ^[0-9]+$ || ${gene_family_zip_compression_level} -gt 9 ]]; then
  echo "Invalid gene_family_zip_compression_level: ${gene_family_zip_compression_level}; expected an integer from 0 through 9." >&2
  exit 1
fi
if [[ ! "${gene_family_zip_workers}" =~ ^[0-9]+$ || ${gene_family_zip_workers} -lt 1 || ${gene_family_zip_workers} -gt 4 ]]; then
  echo "Invalid gene_family_zip_workers: ${gene_family_zip_workers}; expected an integer from 1 through 4." >&2
  exit 1
fi
if [[ ! "${gene_family_final_zip_max_bytes}" =~ ^[0-9]+$ ]]; then
  echo "Invalid gene_family_final_zip_max_bytes: ${gene_family_final_zip_max_bytes}; expected a non-negative integer." >&2
  exit 1
fi
gene_family_archive_write_args=(
  --compression "${gene_family_zip_compression}"
  --compression-level "${gene_family_zip_compression_level}"
  --workers "${gene_family_zip_workers}"
  --max-final-zip-bytes "${gene_family_final_zip_max_bytes}"
)
if [[ ! "${gene_family_tmp_retention_days}" =~ ^[0-9]+$ ]]; then
  echo "Invalid gene_family_tmp_retention_days: ${gene_family_tmp_retention_days}; expected a non-negative integer." >&2
  exit 1
fi
if [[ ! "${gene_family_tmp_max_dirs}" =~ ^[0-9]+$ ]]; then
  echo "Invalid gene_family_tmp_max_dirs: ${gene_family_tmp_max_dirs}; expected a non-negative integer." >&2
  exit 1
fi
if [[ ! "${gene_family_tmp_max_bytes}" =~ ^[0-9]+$ ]]; then
  echo "Invalid gene_family_tmp_max_bytes: ${gene_family_tmp_max_bytes}; expected a non-negative integer." >&2
  exit 1
fi
if [[ ! "${gene_family_tmp_max_files}" =~ ^[0-9]+$ ]]; then
  echo "Invalid gene_family_tmp_max_files: ${gene_family_tmp_max_files}; expected a non-negative integer." >&2
  exit 1
fi

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

write_csubst_foreground_table() {
  local input_tsv=$1
  local output_tsv=$2
  local resolved_tsv="${output_tsv}.binary-resolved.tsv"

  if [[ "${csubst_resolve_binary_foreground}" != "yes" ]]; then
    write_species_trait_foreground_regex_table "${input_tsv}" "${output_tsv}"
    return 0
  fi
  if [[ ! -s "${species_tree_pruned}" ]]; then
    echo "csubst_resolve_binary_foreground=yes requires the pruned species tree: ${species_tree_pruned}"
    exit 1
  fi

  python "${gg_support_dir}/resolve_binary_trait_foregrounds.py" \
    --input "${input_tsv}" \
    --species-tree "${species_tree_pruned}" \
    --output "${resolved_tsv}"
  write_species_trait_foreground_regex_table "${resolved_tsv}" "${output_tsv}"
  rm -f -- "${resolved_tsv}"
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

adopt_existing_output_path() {
  local variable_name=$1
  shift
  local configured_path=${!variable_name}
  local candidate
  if [[ -s "${configured_path}" ]]; then
    return 0
  fi
  for candidate in "$@"; do
    if [[ -s "${candidate}" ]]; then
      printf -v "${variable_name}" '%s' "${candidate}"
      echo "Using historical output path for ${variable_name}: ${candidate}"
      return 0
    fi
  done
}

adopt_historical_gene_family_outputs() {
  adopt_existing_output_path file_og_cds_fasta \
    "${dir_output_active}/cds_fasta/${og_id}_cds.fasta" \
    "${dir_output_active}/cds.fasta/${og_id}.cds.fasta"
  adopt_existing_output_path file_og_rpsblast \
    "${dir_output_active}/rpsblast/${og_id}_rpsblast.tsv" \
    "${dir_output_active}/rpsblast/${og_id}.rpsblast.tsv"
  adopt_existing_output_path file_og_mafft \
    "${dir_output_active}/mafft/${og_id}_cds.aln.fasta" \
    "${dir_output_active}/mafft/${og_id}.cds.aln.fasta"
  adopt_existing_output_path file_og_clipkit \
    "${dir_output_active}/clipkit/${og_id}_cds.clipkit.fasta" \
    "${dir_output_active}/clipkit/${og_id}.cds.clipkit.fasta"
  adopt_existing_output_path file_og_clipkit_log \
    "${dir_output_active}/clipkit_log/${og_id}_cds.clipkit.log" \
    "${dir_output_active}/clipkit.log/${og_id}.cds.clipkit.log"
  adopt_existing_output_path file_og_iqtree_tree \
    "${dir_output_active}/iqtree_tree/${og_id}_iqtree.nwk" \
    "${dir_output_active}/iqtree.tree/${og_id}.iqtree.nwk"
  adopt_existing_output_path file_og_generax_nhx \
    "${dir_output_active}/generax_tree/${og_id}_generax.nhx" \
    "${dir_output_active}/generax.tree/${og_id}.generax.nhx"
  adopt_existing_output_path file_og_generax_nwk \
    "${dir_output_active}/generax_nwk/${og_id}_generax.nwk" \
    "${dir_output_active}/generax.nwk/${og_id}.generax.nwk"
  adopt_existing_output_path file_og_generax_xml \
    "${dir_output_active}/generax_xml/${og_id}_generax.xml" \
    "${dir_output_active}/generax.xml/${og_id}.generax.xml"
  adopt_existing_output_path file_og_mapdnds_dn \
    "${dir_output_active}/mapdnds_dn_tree/${og_id}_mapdNdS.dN.nwk" \
    "${dir_output_active}/mapdNdS.dN.tree/${og_id}.mapdNdS.dN.nwk"
  adopt_existing_output_path file_og_mapdnds_ds \
    "${dir_output_active}/mapdnds_ds_tree/${og_id}_mapdNdS.dS.nwk" \
    "${dir_output_active}/mapdNdS.dS.tree/${og_id}.mapdNdS.dS.nwk"
  adopt_existing_output_path file_og_gff_info \
    "${dir_output_active}/character_gff_info/${og_id}_gff.tsv" \
    "${dir_output_active}/character.gff/${og_id}.gff.tsv"
  adopt_existing_output_path file_og_stat_branch \
    "${dir_output_active}/stat_branch/${og_id}_stat.branch.tsv" \
    "${dir_output_active}/stat.branch/${og_id}.stat.branch.tsv"
  adopt_existing_output_path file_og_stat_tree \
    "${dir_output_active}/stat_tree/${og_id}_stat.tree.tsv" \
    "${dir_output_active}/stat.tree/${og_id}.stat.tree.tsv"
  adopt_existing_output_path file_og_amas_original \
    "${dir_output_active}/amas_original/${og_id}_amas.original.tsv" \
    "${dir_output_active}/amas.original/${og_id}.amas.original.tsv"
  adopt_existing_output_path file_og_amas_cleaned \
    "${dir_output_active}/amas_cleaned/${og_id}_amas.cleaned.tsv" \
    "${dir_output_active}/amas.cleaned/${og_id}.amas.cleaned.tsv"
  adopt_existing_output_path file_og_tree_plot \
    "${dir_output_active}/tree_plot/${og_id}_tree_plot.pdf" \
    "${dir_output_active}/tree_plot/${og_id}.tree_plot.pdf"

  file_og_primary_fasta=${file_og_cds_fasta}
  if [[ "${input_sequence_mode}" == "protein" ]]; then
    file_og_primary_fasta=${file_og_pep_fasta}
  fi
  set_default_analysis_files
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
  local protein_part="${protein_out}.partial.$$"
  if [[ "${protein_out}" == *.gz ]]; then
    # seqkit selects output compression from the final filename suffix.  Keep
    # the staging file on the destination filesystem while retaining .gz as
    # its last suffix so a compressed public artifact cannot contain plain
    # FASTA bytes after the atomic rename.
    protein_part="${protein_out}.partial.$$.gz"
  fi
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
    rm -f -- "${translated_tmp}"
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
    rm -f -- "${translated_tmp}"
    exit 1
  fi
  ensure_parent_dir "${protein_out}"
  rm -f -- "${protein_part}"
  if ! seqkit seq --threads "${GG_TASK_CPUS}" "${translated_tmp}" --out-file "${protein_part}"; then
    rm -f -- "${protein_part}" "${translated_tmp}"
    echo "Failed to write translated protein FASTA: ${protein_out}" >&2
    return 1
  fi
  if [[ ! -s "${protein_part}" ]]; then
    rm -f -- "${protein_part}" "${translated_tmp}"
    echo "Translated protein FASTA is empty: ${protein_out}" >&2
    return 1
  fi
  if ! mv_out "${protein_part}" "${protein_out}"; then
    rm -f -- "${protein_part}" "${translated_tmp}"
    echo "Failed to publish translated protein FASTA: ${protein_out}" >&2
    return 1
  fi
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
    --invert-match yes \
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
  if [[ -n "${gene_family_run_token:-}" && ${gene_family_state_finalized:-1} -eq 0 ]]; then
    if ! python "${gene_family_store_script}" mark-failed \
      --root "${dir_output_active}" \
      --family-id "${og_id}" \
      --run-token "${gene_family_run_token}"
    then
      echo "Warning: Failed to record failed gene-family state for ${og_id}." >&2
    fi
    gene_family_state_finalized=1
  fi
  gg_advisory_shared_lock_release 2>/dev/null || true
  if [[
    ${gene_family_run_succeeded:-0} -eq 0
    && -n "${gene_family_materialization_receipt:-}"
    && -s "${gene_family_materialization_receipt}"
  ]]; then
    if ! python "${gene_family_store_script}" cleanup-materialized \
      --receipt "${gene_family_materialization_receipt}" \
      --nonblocking
    then
      echo "Warning: Failed to remove unchanged materialized artifacts for ${og_id}." >&2
    fi
  fi
  if [[
    ${gene_family_run_succeeded:-0} -eq 0
    && "${gene_family_output_storage}" == "zip"
    && ${gg_debug_mode:-0} -eq 0
    && -n "${gene_family_store_script:-}"
    && -n "${dir_output_active:-}"
    && -n "${og_id:-}"
  ]]; then
    if ! python "${gene_family_store_script}" archive-family \
      "${gene_family_store_context_args[@]}" \
      "${gene_family_archive_write_args[@]}" \
      --family-id "${og_id}" \
      --nonblocking
    then
      echo "Warning: Failed to archive partial outputs for ${og_id}; live files were preserved." >&2
    fi
  fi
  if [[ ${delete_tmp_dir} -eq 1 && (${exit_code} -eq 0 || ${exit_code} -eq 8) ]]; then
    if [[
      -n "${gene_family_materialization_receipt:-}"
      && -s "${gene_family_materialization_receipt}"
    ]]; then
      echo "Retaining ${dir_tmp}: materialization receipt still requires cleanup." >&2
    elif [[ -n "${dir_tmp:-}" && -d "${dir_tmp}" && "${dir_tmp}" != "/" ]]; then
      echo "Deleting ${dir_tmp}"
      rm -rf -- "${dir_tmp}"
    elif [[ -n "${dir_tmp:-}" ]]; then
      echo "Refusing to delete unsafe tmp directory: ${dir_tmp}"
    fi
  fi
  if [[
    "${gene_family_output_storage}" == "zip"
    && (
      ${gene_family_tmp_retention_days} -gt 0
      || ${gene_family_tmp_max_dirs} -gt 0
      || ${gene_family_tmp_max_bytes} -gt 0
      || ${gene_family_tmp_max_files} -gt 0
    )
  ]]; then
    if ! python "${gene_family_store_script}" cleanup-tmp \
      --root "${dir_output_active}" \
      --older-than-days "${gene_family_tmp_retention_days}" \
      --max-directories "${gene_family_tmp_max_dirs}" \
      --max-bytes "${gene_family_tmp_max_bytes}" \
      --max-files "${gene_family_tmp_max_files}" \
      --nonblocking
    then
      echo "Warning: Failed to clean stale or excess gene-family temporary directories." >&2
    fi
  fi
  gg_shared_lock_stop_heartbeat "${gene_family_run_lock_heartbeat_pid:-}"
  if [[ -n "${gene_family_run_lock_path:-}" ]]; then
    gg_shared_lock_release "${gene_family_run_lock_path}" 2>/dev/null || true
    gene_family_run_lock_path=""
  fi
  return ${exit_code}
}

finalize_gene_family_run_success() {
  if [[ -n "${gene_family_run_token:-}" && ${gene_family_state_finalized:-1} -eq 0 ]]; then
    if ! python "${gene_family_store_script}" mark-complete \
      --root "${dir_output_active}" \
      --family-id "${og_id}" \
      --run-token "${gene_family_run_token}"
    then
      return 1
    fi
    gene_family_state_finalized=1
  fi
  if ! gg_advisory_shared_lock_release; then
    return 1
  fi
  gene_family_run_succeeded=1
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
      nwkit mark --pattern "${fg_regex}" --target "clade" --target-only-clade "yes" --insert-text "{Foreground}" --outformat 1 |
      nwkit sanitize --remove-singleton yes --resolve-polytomy no \
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

rsc_validate_choice() {
  local variable_name=$1
  local value=$2
  shift 2
  local choice
  for choice in "$@"; do
    if [[ "${value}" == "${choice}" ]]; then
      return 0
    fi
  done
  echo "Invalid ${variable_name}: ${value}; expected one of: $*." >&2
  exit 1
}

rsc_validate_nonnegative_integer() {
  local variable_name=$1
  local value=$2
  if [[ ! "${value}" =~ ^[0-9]+$ ]]; then
    echo "Invalid ${variable_name}: ${value}; expected a non-negative integer." >&2
    exit 1
  fi
}

rsc_validate_positive_integer() {
  local variable_name=$1
  local value=$2
  rsc_validate_nonnegative_integer "${variable_name}" "${value}"
  if [[ ${value} -lt 1 ]]; then
    echo "Invalid ${variable_name}: ${value}; expected a positive integer." >&2
    exit 1
  fi
}

rsc_validate_probability() {
  local variable_name=$1
  local value=$2
  if [[ ! "${value}" =~ ^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$ ]] ||
    ! awk -v value="${value}" 'BEGIN { exit !(value > 0 && value < 1) }'
  then
    echo "Invalid ${variable_name}: ${value}; expected a number strictly between 0 and 1." >&2
    exit 1
  fi
}

rsc_validate_evolution_parameter() {
  local variable_name=$1
  local value=$2
  if [[ "${value}" == "auto" ]]; then
    return 0
  fi
  if [[ ! "${value}" =~ ^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$ ]]; then
    echo "Invalid ${variable_name}: ${value}; expected auto or a finite numeric value." >&2
    exit 1
  fi
}

rsc_read_metadata_value() {
  local metadata_path=$1
  local requested_key=$2
  awk -F '\t' -v key="${requested_key}" '$1 == key { print $2; exit }' "${metadata_path}"
}












# Setting modes
if [[ ${gg_debug_mode:-0} -eq 1 ]]; then
  enable_all_run_flags_for_debug_mode "gg debug mode: All run_* variables are forced to set 1, except for too-time-consuming tasks."
  #run_orthogroup_extraction=0; echo "gg debug mode: run_orthogroup_extraction=${run_orthogroup_extraction}"
  run_codeml_two_ratio=0
  echo "gg debug mode: run_codeml_two_ratio=${run_codeml_two_ratio}"
  run_hyphy_relax=0
  echo "gg debug mode: run_hyphy_relax=${run_hyphy_relax}"
  run_hyphy_relax_reversed=0
  echo "gg debug mode: run_hyphy_relax_reversed=${run_hyphy_relax_reversed}"
  run_expression_trait_pgls=0
  echo "gg debug mode: run_expression_trait_pgls=${run_expression_trait_pgls}"
  csubst_cutoff_stat="OCNany2spe,0|omegaCany2spe,1"
  echo "gg debug mode: csubst_cutoff_stat=${csubst_cutoff_stat}"
fi
query_blast_method=$(echo "${query_blast_method}" | tr '[:upper:]' '[:lower:]')
mode_gene_evolution=$(echo "${mode_gene_evolution:-query2family}" | tr '[:upper:]' '[:lower:]')
gene_evolution_profile=$(echo "${gene_evolution_profile:-default}" | tr '[:upper:]' '[:lower:]')
input_sequence_mode=$(gg_normalize_input_sequence_mode "${input_sequence_mode}")
csubst_nonsyn_recode=$(echo "${csubst_nonsyn_recode:-${GG_COMMON_CSUBST_NONSYN_RECODE:-no}}" | tr '[:upper:]' '[:lower:]')
csubst_resolve_binary_foreground=$(echo "${csubst_resolve_binary_foreground:-no}" | tr '[:upper:]' '[:lower:]')
csubst_scan_unit_mode=$(echo "${csubst_scan_unit_mode:-clade}" | tr '[:upper:]' '[:lower:]')
csubst_scan_match=$(echo "${csubst_scan_match:-any2spe}" | tr '[:upper:]' '[:lower:]')
csubst_scan_min_event_pp="${csubst_scan_min_event_pp:-0.5}"
csubst_scan_min_support="${csubst_scan_min_support:-2}"
csubst_scan_rate_event_mode=$(echo "${csubst_scan_rate_event_mode:-posterior_sum}" | tr '[:upper:]' '[:lower:]')
csubst_scan_rate_length=$(echo "${csubst_scan_rate_length:-n_rescaled}" | tr '[:upper:]' '[:lower:]')
csubst_scan_rate_exposure=$(echo "${csubst_scan_rate_exposure:-q_weighted}" | tr '[:upper:]' '[:lower:]')
csubst_scan_other_scope=$(echo "${csubst_scan_other_scope:-all}" | tr '[:upper:]' '[:lower:]')
csubst_scan_pvalue_calibration=$(echo "${csubst_scan_pvalue_calibration:-full_scan}" | tr '[:upper:]' '[:lower:]')
csubst_scan_n_permutations="${csubst_scan_n_permutations:-1000}"
csubst_scan_site_plot=$(echo "${csubst_scan_site_plot:-yes}" | tr '[:upper:]' '[:lower:]')
csubst_scan_tree_site_plot_format=$(echo "${csubst_scan_tree_site_plot_format:-pdf}" | tr '[:upper:]' '[:lower:]')
csubst_scan_tree_site_plot_max_sites="${csubst_scan_tree_site_plot_max_sites:-30}"
uniprot_annotation_method=$(echo "${uniprot_annotation_method:-mmseqs2}" | tr '[:upper:]' '[:lower:]')
tree_rooting_method=$(echo "${tree_rooting_method}" | tr '[:upper:]' '[:lower:]')
pgls_methods=$(printf '%s' "${pgls_methods}" | tr '[:upper:]' '[:lower:]' | tr -d '[:space:]')
species_expression_aggregation=$(printf '%s' "${species_expression_aggregation}" | tr '[:upper:]' '[:lower:]')
species_paralog_missing=$(printf '%s' "${species_paralog_missing}" | tr '[:upper:]' '[:lower:]')
rphylopars_sampling_covariance=$(printf '%s' "${rphylopars_sampling_covariance}" | tr '[:upper:]' '[:lower:]')
rsc_predictor_mode=$(printf '%s' "${rsc_predictor_mode}" | tr '[:upper:]' '[:lower:]')
rsc_event_source=$(printf '%s' "${rsc_event_source}" | tr '[:upper:]' '[:lower:]')
rsc_speciation_coverage=$(printf '%s' "${rsc_speciation_coverage}" | tr '[:upper:]' '[:lower:]')
rsc_event_weighting=$(printf '%s' "${rsc_event_weighting}" | tr '[:upper:]' '[:lower:]')
rsc_model=$(printf '%s' "${rsc_model}" | tr '[:upper:]' '[:lower:]')
rsc_gene_branch_length=$(printf '%s' "${rsc_gene_branch_length}" | tr '[:upper:]' '[:lower:]')
rsc_gene_evolution_model=$(printf '%s' "${rsc_gene_evolution_model}" | tr '[:upper:]' '[:lower:]')
rsc_gene_evolution_parameter=$(printf '%s' "${rsc_gene_evolution_parameter}" | tr '[:upper:]' '[:lower:]')
rsc_species_branch_length=$(printf '%s' "${rsc_species_branch_length}" | tr '[:upper:]' '[:lower:]')
rsc_species_evolution_model=$(printf '%s' "${rsc_species_evolution_model}" | tr '[:upper:]' '[:lower:]')
rsc_species_evolution_parameter=$(printf '%s' "${rsc_species_evolution_parameter}" | tr '[:upper:]' '[:lower:]')
rsc_inference=$(printf '%s' "${rsc_inference}" | tr '[:upper:]' '[:lower:]')
rsc_reml=$(printf '%s' "${rsc_reml}" | tr '[:upper:]' '[:lower:]')
rsc_unmatched=$(printf '%s' "${rsc_unmatched}" | tr '[:upper:]' '[:lower:]')
rsc_within_variance=$(printf '%s' "${rsc_within_variance}" | tr '[:upper:]' '[:lower:]')
rsc_technical_aggregation=$(printf '%s' "${rsc_technical_aggregation}" | tr '[:upper:]' '[:lower:]')
rsc_predictor_within_variance=$(printf '%s' "${rsc_predictor_within_variance}" | tr '[:upper:]' '[:lower:]')
rsc_predictor_technical_aggregation=$(printf '%s' "${rsc_predictor_technical_aggregation}" | tr '[:upper:]' '[:lower:]')
rsc_factor_coding=$(printf '%s' "${rsc_factor_coding}" | tr '[:upper:]' '[:lower:]')
rsc_categorical_replicate_policy=$(printf '%s' "${rsc_categorical_replicate_policy}" | tr '[:upper:]' '[:lower:]')
rsc_event_random_effect=$(printf '%s' "${rsc_event_random_effect}" | tr '[:upper:]' '[:lower:]')
rsc_lineage_random_slope=$(printf '%s' "${rsc_lineage_random_slope}" | tr '[:upper:]' '[:lower:]')
rsc_lineage_inference=$(printf '%s' "${rsc_lineage_inference}" | tr '[:upper:]' '[:lower:]')
rsc_lineage_leave_one_out=$(printf '%s' "${rsc_lineage_leave_one_out}" | tr '[:upper:]' '[:lower:]')
rsc_categorical_origin_diagnostics=$(printf '%s' "${rsc_categorical_origin_diagnostics}" | tr '[:upper:]' '[:lower:]')
rsc_origin_leave_one_out=$(printf '%s' "${rsc_origin_leave_one_out}" | tr '[:upper:]' '[:lower:]')
rsc_allow_large_dense=$(printf '%s' "${rsc_allow_large_dense}" | tr '[:upper:]' '[:lower:]')
apply_gene_evolution_profile
pgls_run_rsc=0
pgls_run_species_rphylopars=0
if [[ "${pgls_methods}" == "all" ]]; then
  pgls_methods="rsc,species-nwkit,species-rphylopars"
fi
if [[ -z "${pgls_methods}" ]]; then
  echo "pgls_methods must select at least one method." >&2
  exit 1
fi
pgls_seen_methods=","
IFS=',' read -r -a pgls_method_items <<< "${pgls_methods}"
for pgls_method in "${pgls_method_items[@]}"; do
  if [[ "${pgls_seen_methods}" == *",${pgls_method},"* ]]; then
    echo "pgls_methods contains a duplicate method: ${pgls_method}" >&2
    exit 1
  fi
  pgls_seen_methods+="${pgls_method},"
  case "${pgls_method}" in
    rsc) pgls_run_rsc=1 ;;
    species-nwkit) : ;;
    species-rphylopars) pgls_run_species_rphylopars=1 ;;
    *)
      echo "Invalid pgls_methods entry: ${pgls_method}. Use rsc, species-nwkit, species-rphylopars, or all." >&2
      exit 1
      ;;
  esac
done
unset pgls_method pgls_method_items pgls_seen_methods
if [[ ${run_expression_trait_pgls} -eq 1 ]]; then
  rsc_validate_choice species_expression_aggregation "${species_expression_aggregation}" sum mean max all
  rsc_validate_choice species_paralog_missing "${species_paralog_missing}" error ignore
  rsc_validate_choice rphylopars_sampling_covariance "${rphylopars_sampling_covariance}" require-diagonal diagonalize
  rsc_validate_choice rsc_predictor_mode "${rsc_predictor_mode}" separate joint
  rsc_validate_choice rsc_event_source "${rsc_event_source}" auto nhx lca species-overlap
  rsc_validate_choice rsc_speciation_coverage "${rsc_speciation_coverage}" complete any
  rsc_validate_choice rsc_event_weighting "${rsc_event_weighting}" equal observation
  rsc_validate_choice rsc_model "${rsc_model}" hierarchical replicate-reml legacy
  rsc_validate_choice rsc_gene_branch_length "${rsc_gene_branch_length}" original unit
  rsc_validate_choice rsc_species_branch_length "${rsc_species_branch_length}" original unit
  rsc_validate_choice rsc_gene_evolution_model "${rsc_gene_evolution_model}" brownian lambda ou kappa delta eb acdc independent
  rsc_validate_choice rsc_species_evolution_model "${rsc_species_evolution_model}" brownian lambda ou kappa delta eb acdc independent
  rsc_validate_choice rsc_inference "${rsc_inference}" wald parametric-bootstrap
  rsc_validate_choice rsc_reml "${rsc_reml}" yes no
  rsc_validate_choice rsc_unmatched "${rsc_unmatched}" error warn ignore
  rsc_validate_choice rsc_within_variance "${rsc_within_variance}" pooled leaf known-se
  rsc_validate_choice rsc_technical_aggregation "${rsc_technical_aggregation}" error mean
  rsc_validate_choice rsc_predictor_within_variance "${rsc_predictor_within_variance}" pooled leaf known-se
  rsc_validate_choice rsc_predictor_technical_aggregation "${rsc_predictor_technical_aggregation}" error mean
  rsc_validate_choice rsc_factor_coding "${rsc_factor_coding}" treatment sum
  rsc_validate_choice rsc_categorical_replicate_policy "${rsc_categorical_replicate_policy}" error latent
  rsc_validate_choice rsc_event_random_effect "${rsc_event_random_effect}" auto yes no
  rsc_validate_choice rsc_lineage_random_slope "${rsc_lineage_random_slope}" auto yes no
  rsc_validate_choice rsc_lineage_inference "${rsc_lineage_inference}" none likelihood-ratio parametric-bootstrap
  rsc_validate_choice rsc_lineage_leave_one_out "${rsc_lineage_leave_one_out}" yes no
  rsc_validate_choice rsc_categorical_origin_diagnostics "${rsc_categorical_origin_diagnostics}" none stochastic-map
  rsc_validate_choice rsc_origin_leave_one_out "${rsc_origin_leave_one_out}" yes no
  rsc_validate_choice rsc_allow_large_dense "${rsc_allow_large_dense}" yes no
  rsc_validate_evolution_parameter rsc_gene_evolution_parameter "${rsc_gene_evolution_parameter}"
  rsc_validate_evolution_parameter rsc_species_evolution_parameter "${rsc_species_evolution_parameter}"
  if [[ "${rsc_gene_evolution_model}" =~ ^(brownian|independent)$ && "${rsc_gene_evolution_parameter}" != "auto" ]]; then
    echo "${rsc_gene_evolution_model} has no gene-evolution shape parameter; set rsc_gene_evolution_parameter=auto." >&2
    exit 1
  fi
  if [[ "${rsc_species_evolution_model}" =~ ^(brownian|independent)$ && "${rsc_species_evolution_parameter}" != "auto" ]]; then
    echo "${rsc_species_evolution_model} has no species-evolution shape parameter; set rsc_species_evolution_parameter=auto." >&2
    exit 1
  fi
  rsc_validate_positive_integer rsc_bootstrap_replicates "${rsc_bootstrap_replicates}"
  if (( rsc_bootstrap_replicates < 2 )); then
    echo "rsc_bootstrap_replicates must be at least 2." >&2
    exit 1
  fi
  rsc_validate_nonnegative_integer rsc_seed "${rsc_seed}"
  rsc_validate_positive_integer rsc_min_species_events "${rsc_min_species_events}"
  rsc_validate_positive_integer rsc_origin_map_replicates "${rsc_origin_map_replicates}"
  rsc_validate_positive_integer rsc_origin_map_threads "${rsc_origin_map_threads}"
  rsc_validate_probability rsc_confidence_level "${rsc_confidence_level}"
  rsc_validate_probability rsc_origin_min_posterior "${rsc_origin_min_posterior}"
  if [[ -z "${rsc_responses//[[:space:]]/}" || -z "${rsc_predictors//[[:space:]]/}" ]]; then
    echo "rsc_responses and rsc_predictors must not be empty when run_expression_trait_pgls=1." >&2
    exit 1
  fi
  if [[ "${rsc_within_variance}" == "known-se" && -z "${rsc_expression_sample_metadata}" ]]; then
    echo "rsc_within_variance=known-se requires rsc_expression_sample_metadata." >&2
    exit 1
  fi
  if [[ -n "${species_paralog_sampling_covariance}" && "${rsc_within_variance}" != "known-se" ]]; then
    echo "species_paralog_sampling_covariance requires rsc_within_variance=known-se; raw paired biological replicates already preserve cross-paralog covariance after aggregation." >&2
    exit 1
  fi
  if [[ "${rsc_within_variance}" == "known-se" && "${rsc_technical_aggregation}" != "error" ]]; then
    echo "rsc_within_variance=known-se cannot use rsc_technical_aggregation; provide one summarized mean and SE per gene/response." >&2
    exit 1
  fi
  if [[ -n "${rsc_predictor_technical_id}${rsc_predictor_batch}" && -z "${rsc_predictor_biological_id}" ]]; then
    echo "RSC predictor technical-replicate and batch columns require rsc_predictor_biological_id." >&2
    exit 1
  fi
  if [[ "${rsc_predictor_within_variance}" == "known-se" && -z "${rsc_predictor_standard_error_columns}" ]]; then
    echo "rsc_predictor_within_variance=known-se requires rsc_predictor_standard_error_columns." >&2
    exit 1
  fi
  if [[ "${rsc_predictor_within_variance}" == "known-se" && -n "${rsc_predictor_biological_id}${rsc_predictor_technical_id}${rsc_predictor_batch}" ]]; then
    echo "rsc_predictor_within_variance=known-se requires one summarized row per species and cannot use predictor replicate IDs or batch." >&2
    exit 1
  fi
  if [[ "${rsc_predictor_within_variance}" == "known-se" && "${rsc_predictor_technical_aggregation}" != "error" ]]; then
    echo "rsc_predictor_within_variance=known-se cannot use rsc_predictor_technical_aggregation." >&2
    exit 1
  fi
  if [[ -n "${rsc_predictor_standard_error_columns}${rsc_predictor_sample_size_columns}" && "${rsc_predictor_within_variance}" != "known-se" ]]; then
    echo "RSC predictor standard-error/sample-size columns require rsc_predictor_within_variance=known-se." >&2
    exit 1
  fi
  if [[ "${rsc_origin_leave_one_out}" == "yes" && "${rsc_categorical_origin_diagnostics}" != "stochastic-map" ]]; then
    echo "rsc_origin_leave_one_out=yes requires rsc_categorical_origin_diagnostics=stochastic-map." >&2
    exit 1
  fi
  if [[ "${rsc_model}" == "legacy" && "${rsc_inference}" != "wald" ]]; then
    echo "rsc_model=legacy supports only rsc_inference=wald." >&2
    exit 1
  fi
  if [[ "${rsc_model}" == "legacy" && ! "${rsc_gene_evolution_model}" =~ ^(brownian|independent)$ && "${rsc_gene_evolution_parameter}" == "auto" ]]; then
    echo "Automatic gene evolution-parameter estimation requires a likelihood-based RSC model." >&2
    exit 1
  fi
  if [[ "${rsc_model}" != "hierarchical" && ( "${rsc_event_random_effect}" == "yes" || "${rsc_lineage_random_slope}" == "yes" ) ]]; then
    echo "Explicit RSC random effects require rsc_model=hierarchical." >&2
    exit 1
  fi
  if [[ "${rsc_model}" != "hierarchical" && "${rsc_lineage_inference}" != "none" ]]; then
    echo "rsc_lineage_inference requires rsc_model=hierarchical." >&2
    exit 1
  fi
  if [[ "${rsc_lineage_random_slope}" == "no" && "${rsc_lineage_inference}" != "none" ]]; then
    echo "rsc_lineage_inference requires lineage random slopes." >&2
    exit 1
  fi
  if [[ "${rsc_model}" == "legacy" && "${rsc_lineage_leave_one_out}" == "yes" ]]; then
    echo "rsc_lineage_leave_one_out requires a likelihood-based RSC model." >&2
    exit 1
  fi
  if [[ "${rsc_model}" == "legacy" && ( "${rsc_within_variance}" == "known-se" || "${rsc_predictor_within_variance}" == "known-se" || -n "${rsc_predictor_biological_id}" ) ]]; then
    echo "rsc_model=legacy cannot use response or predictor sampling uncertainty." >&2
    exit 1
  fi
fi
if [[ "${mode_gene_evolution}" != "orthogroup" && "${mode_gene_evolution}" != "query2family" ]]; then
  echo "Invalid mode_gene_evolution: ${mode_gene_evolution}"
  echo 'mode_gene_evolution must be either "orthogroup" or "query2family". Exiting.'
  exit 1
fi
if [[ "${tree_rooting_method}" != "notung" && "${tree_rooting_method}" != "midpoint" && "${tree_rooting_method}" != "mad" && "${tree_rooting_method}" != "md" ]]; then
  echo "Invalid tree_rooting_method: ${tree_rooting_method}"
  echo "tree_rooting_method must be one of notung, midpoint, mad, md. Exiting."
  exit 1
fi
if [[ "${uniprot_annotation_method}" != "blastp" && "${uniprot_annotation_method}" != "mmseqs2" ]]; then
  echo "Invalid uniprot_annotation_method: ${uniprot_annotation_method}"
  echo 'uniprot_annotation_method must be either "blastp" or "mmseqs2". Exiting.'
  exit 1
fi
case "${csubst_nonsyn_recode}" in
  no|3di20|dayhoff6|sr6|kgb6|sr4|dayhoff9|dayhoff12|dayhoff15|dayhoff18|srchisq6|kgbauto6)
    ;;
  *)
    echo "Invalid csubst_nonsyn_recode: ${csubst_nonsyn_recode}"
    echo 'csubst_nonsyn_recode must be one of no, 3di20, dayhoff6, sr6, kgb6, sr4, dayhoff9, dayhoff12, dayhoff15, dayhoff18, srchisq6, kgbauto6. Exiting.'
    exit 1
    ;;
esac
case "${csubst_resolve_binary_foreground}" in
  yes|no)
    ;;
  *)
    echo "Invalid csubst_resolve_binary_foreground: ${csubst_resolve_binary_foreground}"
    echo 'csubst_resolve_binary_foreground must be either "yes" or "no". Exiting.'
    exit 1
    ;;
esac
case "${csubst_scan_unit_mode}" in
  lineage|stem|clade)
    ;;
  *)
    echo "Invalid csubst_scan_unit_mode: ${csubst_scan_unit_mode}"
    echo 'csubst_scan_unit_mode must be one of lineage, stem, clade. Exiting.'
    exit 1
    ;;
esac
case "${csubst_scan_rate_event_mode}" in
  posterior_sum|called)
    ;;
  *)
    echo "Invalid csubst_scan_rate_event_mode: ${csubst_scan_rate_event_mode}"
    echo 'csubst_scan_rate_event_mode must be either "posterior_sum" or "called". Exiting.'
    exit 1
    ;;
esac
case "${csubst_scan_rate_length}" in
  raw|sn_rescaled|n_rescaled)
    ;;
  *)
    echo "Invalid csubst_scan_rate_length: ${csubst_scan_rate_length}"
    echo 'csubst_scan_rate_length must be one of raw, sn_rescaled, n_rescaled. Exiting.'
    exit 1
    ;;
esac
case "${csubst_scan_rate_exposure}" in
  q_weighted|state_aware|raw_branch_length)
    ;;
  *)
    echo "Invalid csubst_scan_rate_exposure: ${csubst_scan_rate_exposure}"
    echo 'csubst_scan_rate_exposure must be one of q_weighted, state_aware, raw_branch_length. Exiting.'
    exit 1
    ;;
esac
case "${csubst_scan_other_scope}" in
  all|sister)
    ;;
  *)
    echo "Invalid csubst_scan_other_scope: ${csubst_scan_other_scope}"
    echo 'csubst_scan_other_scope must be either "all" or "sister". Exiting.'
    exit 1
    ;;
esac
case "${csubst_scan_pvalue_calibration}" in
  none|candidate_fixed|full_scan)
    ;;
  *)
    echo "Invalid csubst_scan_pvalue_calibration: ${csubst_scan_pvalue_calibration}"
    echo 'csubst_scan_pvalue_calibration must be one of none, candidate_fixed, full_scan. Exiting.'
    exit 1
    ;;
esac
case "${csubst_scan_site_plot}" in
  yes|no)
    ;;
  *)
    echo "Invalid csubst_scan_site_plot: ${csubst_scan_site_plot}"
    echo 'csubst_scan_site_plot must be either "yes" or "no". Exiting.'
    exit 1
    ;;
esac
case "${csubst_scan_tree_site_plot_format}" in
  pdf|png|svg)
    ;;
  *)
    echo "Invalid csubst_scan_tree_site_plot_format: ${csubst_scan_tree_site_plot_format}"
    echo 'csubst_scan_tree_site_plot_format must be one of pdf, png, svg. Exiting.'
    exit 1
    ;;
esac
apply_gene_evolution_input_sequence_mode
if [[ "${mode_gene_evolution}" == "query2family" && ${run_query_blast} -eq 1 ]]; then
  if [[ "${query_blast_method}" != "tblastn" && "${query_blast_method}" != "diamond" ]]; then
    echo "Invalid query_blast_method: ${query_blast_method}"
    echo 'query_blast_method must be either "tblastn" or "diamond". Exiting.'
    exit 1
  fi
fi
case "${mode_gene_evolution}" in
  "orthogroup")
    dir_output_active="${gg_workspace_output_dir}/orthogroup"
    file_orthogroup_genecount_selected="${gg_workspace_output_dir}/orthofinder/Orthogroups_filtered/Orthogroups.GeneCount.selected.tsv"
    if [[ ! -s "${file_orthogroup_genecount_selected}" ]]; then
      echo "Orthogroup gene-count table not found: ${file_orthogroup_genecount_selected}"
      exit 1
    fi
    if [[ ! "${GG_ARRAY_TASK_ID}" =~ ^[0-9]+$ ]] || [[ ${GG_ARRAY_TASK_ID} -lt 1 ]]; then
      echo "Invalid GG_ARRAY_TASK_ID value (must be a positive integer): ${GG_ARRAY_TASK_ID}"
      exit 1
    fi
    num_orthogroups=$(awk 'END { print (NR > 0 ? NR - 1 : 0) }' "${file_orthogroup_genecount_selected}")
    if [[ ${num_orthogroups} -le 0 ]]; then
      echo "No orthogroup rows were found in: ${file_orthogroup_genecount_selected}"
      exit 1
    fi
    if [[ ${GG_ARRAY_TASK_ID} -gt ${num_orthogroups} ]]; then
      echo "GG_ARRAY_TASK_ID=${GG_ARRAY_TASK_ID} is out of range for ${num_orthogroups} orthogroups in: ${file_orthogroup_genecount_selected}"
      exit 1
    fi
    og_id=$(awk -F'\t' -v row="${GG_ARRAY_TASK_ID}" 'NR == (row + 1) { print $1; exit }' "${file_orthogroup_genecount_selected}")
    if [[ -z "${og_id}" ]]; then
      echo "Failed to read OrthoGroup ID at row ${GG_ARRAY_TASK_ID} from: ${file_orthogroup_genecount_selected}"
      exit 1
    fi
    echo "OrthoGroup ID: ${og_id}"
    run_extract_query_fasta=0
    run_query_blast=0
    run_orthogroup_extraction=0
    ;;
  "query2family")
    dir_output_active="${gg_workspace_output_dir}/query2family"
    dir_genelist="${gg_workspace_input_dir}/query_gene"
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
    if [[ ! "${GG_ARRAY_TASK_ID}" =~ ^[0-9]+$ ]] || [[ ${GG_ARRAY_TASK_ID} -lt 1 ]]; then
      echo "Invalid GG_ARRAY_TASK_ID value (must be a positive integer): ${GG_ARRAY_TASK_ID}"
      exit 1
    fi
    idx=$((GG_ARRAY_TASK_ID - 1))
    if [[ ${idx} -ge ${#files[@]} ]]; then
      echo "Input genelist file not found, probably due to out-of-range GG_ARRAY_TASK_ID=${GG_ARRAY_TASK_ID}. Number of files: ${#files[@]}"
      exit 1
    fi
    file_query_gene="${files[${idx}]}"
    og_id="$(basename "${file_query_gene}")"
    echo "mode_gene_evolution=query2family: ${#files[@]} input genelist file(s) were detected in ${dir_genelist}/"
    echo "mode_gene_evolution=query2family: Input genelist file = ${file_query_gene}"
    echo "output file prefix: ${og_id}"
    if [[ ! -f "${file_query_gene}" ]]; then
      echo "Input genelist file not found, probably due to out-of-range GG_ARRAY_TASK_ID=${GG_ARRAY_TASK_ID}: ${file_query_gene}"
      exit 1
    fi
    ;;
  *)
    echo "Invalid mode_gene_evolution: ${mode_gene_evolution}. Exiting."
    exit 1
    ;;
esac
if [[ -z "$og_id" ]]; then
  echo "og_id is empty. Exiting."
  exit 1
fi

gene_family_store_script="${gg_support_dir}/gene_family_output_store.py"
gene_family_store_context_args=(
  --root "${dir_output_active}"
  --mode "${mode_gene_evolution}"
)
if [[ "${mode_gene_evolution}" == "query2family" ]]; then
  gene_family_store_context_args+=(--query-dir "${dir_genelist}")
else
  gene_family_store_context_args+=(--genecount "${file_orthogroup_genecount_selected}")
fi
for gene_family_storage_conversion_marker in \
  "${dir_output_active}/.gg_store/storage-conversion.pending" \
  "${dir_output_active}/.gg_archives/storage-conversion.pending"
do
  if [[ -e "${gene_family_storage_conversion_marker}" || -L "${gene_family_storage_conversion_marker}" ]]; then
    echo "Gene-family storage conversion is in progress or needs to be resumed: ${gene_family_storage_conversion_marker}" >&2
    echo "Refusing to start gg_gene_evolution until convert-storage finishes." >&2
    exit 1
  fi
done
if [[
  "${gene_family_output_storage}" == "files"
  && (
    -d "${dir_output_active}/.gg_store"
    || -d "${dir_output_active}/.gg_archives"
  )
]]; then
  echo "Warning: ZIP-backed artifacts exist below ${dir_output_active}, but gene_family_output_storage=${gene_family_output_storage}." >&2
  echo "The family can resume, but new outputs will remain raw; use storage mode zip to return them to ZIP." >&2
fi
gene_family_run_lock_path="${dir_output_active}/.gg_run_locks/task.${GG_ARRAY_TASK_ID}.lock"
if ! gg_shared_lock_acquire \
  "${gene_family_run_lock_path}" \
  "gene-family producer (${og_id})"
then
  echo "Failed to acquire the exclusive gene-family producer lock." >&2
  exit 1
fi
trap cleanup_tmp_dir_on_normal_exit EXIT
gg_shared_lock_start_heartbeat "${gene_family_run_lock_path}"
gene_family_run_lock_heartbeat_pid=${GG_SHARED_LOCK_HEARTBEAT_PID:-}
gene_family_task_tmp_dir="${dir_output_active}/tmp/${GG_ARRAY_TASK_ID}_${og_id}"
if [[ -e "${gene_family_task_tmp_dir}" && ${delete_preexisting_tmp_dir} -eq 1 ]]; then
  echo "$(date): Deleting preexisting ${gene_family_task_tmp_dir}"
  stale_receipt="${gene_family_task_tmp_dir}/.gg_materialized.jsonl"
  if [[ -s "${stale_receipt}" ]]; then
    if ! python "${gene_family_store_script}" cleanup-materialized \
      --receipt "${stale_receipt}"
    then
      echo "Failed to clean the preexisting materialization receipt: ${stale_receipt}" >&2
      exit 1
    fi
    if [[ -s "${stale_receipt}" ]]; then
      echo "Refusing to remove tmp directory while its materialization receipt remains: ${stale_receipt}" >&2
      exit 1
    fi
  fi
  rm -rf -- "${gene_family_task_tmp_dir}"
fi
if [[
  "${gene_family_output_storage}" == "zip"
  || -d "${dir_output_active}/.gg_store"
  || -d "${dir_output_active}/.gg_archives"
]]; then
  gene_family_lock_path=$(python "${gene_family_store_script}" lock-path \
    --root "${dir_output_active}" \
    --family-id "${og_id}")
  if ! gg_advisory_shared_lock_acquire \
    "${gene_family_lock_path}"
  then
    echo "Failed to acquire the gene-family lock." >&2
    exit 1
  fi
fi
if [[ "${gene_family_output_storage}" == "zip" ]]; then
  gene_family_run_token="${GG_JOB_ID:-local}_${GG_ARRAY_TASK_ID:-1}_$$_$(date +%s)"
  gene_family_state_finalized=0
  if ! python "${gene_family_store_script}" mark-running \
    --root "${dir_output_active}" \
    --family-id "${og_id}" \
    --run-token "${gene_family_run_token}"
  then
    echo "Failed to record running gene-family state for ${og_id}." >&2
    exit 1
  fi
fi

dir_sp_genome="${gg_workspace_input_dir}/species_genome"
dir_sp_gff="${gg_workspace_input_dir}/species_gff"
dir_sp_expression="${gg_workspace_input_dir}/species_expression"
dir_sp_cds="${gg_workspace_input_dir}/species_cds"
dir_sp_protein_input="$(gg_species_protein_input_dir_path "${gg_workspace_input_dir}")"
dir_sp_blastdb="${gg_workspace_output_dir}/species_cds_blastdb"
dir_fasta_sequence_store="${gg_workspace_output_dir}/.gg_cache/fasta_sequence_store"
file_species_cds_store_db="${dir_fasta_sequence_store}/species_cds.sqlite3"
file_species_cds_store_manifest="${dir_fasta_sequence_store}/species_cds.json"
file_species_protein_store_db="${dir_fasta_sequence_store}/species_protein.sqlite3"
file_species_protein_store_manifest="${dir_fasta_sequence_store}/species_protein.json"
file_species_genetic_code="$(gg_species_genetic_code_table_path "${gg_workspace_input_dir}")"
file_species_genetic_code_resolved="${dir_output_active}/parameters/${og_id}_species_genetic_code.resolved.tsv"
annotation_species_resolved=""
treevis_clade_ortholog_prefix=""
annotation_species_candidates=()
if [[ "${input_sequence_mode}" == "protein" ]]; then
  mapfile -t annotation_species_candidates < <(gg_species_names_from_fasta_dir "${dir_sp_protein_input}")
fi
if [[ ${#annotation_species_candidates[@]} -eq 0 ]]; then
  mapfile -t annotation_species_candidates < <(gg_species_names_from_fasta_dir "${dir_sp_cds}")
fi
if annotation_species_resolved=$(gg_resolve_annotation_species "${annotation_species}" "${annotation_species_candidates[@]}"); then
  if [[ -n "${annotation_species_resolved}" ]]; then
    treevis_clade_ortholog_prefix="${annotation_species_resolved}_"
  fi
fi
file_sp_trait="${gg_workspace_input_dir}/species_trait/species_trait.tsv"
if [[ -n "${rsc_expression_sample_metadata}" && "${rsc_expression_sample_metadata}" != /* ]]; then
  rsc_expression_sample_metadata="${gg_workspace_dir}/${rsc_expression_sample_metadata}"
fi
if [[ -n "${species_paralog_sampling_covariance}" && "${species_paralog_sampling_covariance}" != /* ]]; then
  species_paralog_sampling_covariance="${gg_workspace_dir}/${species_paralog_sampling_covariance}"
fi
file_og="${gg_workspace_output_dir}/orthofinder/Orthogroups_filtered/Orthogroups.selected.tsv"
file_og_parameters_dir="${dir_output_active}/parameters"
dir_species_tree="${gg_workspace_output_dir}/species_tree"
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
dir_tmp="${gene_family_task_tmp_dir}" #_${RANDOM}
gene_family_materialization_receipt="${dir_tmp}/.gg_materialized.jsonl"

# File PATHs
# Alignment and gene tree preparation and others
file_og_query_aa_fasta="${dir_output_active}/query_aa_fasta/${og_id}_query.aa.fa.gz"
file_og_query_blast="${dir_output_active}/query_blast/${og_id}_query_blast.tsv"
file_og_cds_fasta="${dir_output_active}/cds_fasta/${og_id}_cds.fa.gz"
file_og_pep_fasta="${dir_output_active}/protein_fasta/${og_id}_pep.fa.gz"
file_og_rpsblast="${dir_output_active}/rpsblast/${og_id}_rpsblast.tsv"
file_og_uniprot_annotation="${dir_output_active}/uniprot_annotation/${og_id}_uniprot.tsv"
file_og_cdskit_localize="${dir_output_active}/cdskit_localize/${og_id}_cdskit_localize.tsv"
file_og_mafft="${dir_output_active}/mafft/${og_id}_cds.aln.fa.gz"
file_og_maxalign="${dir_output_active}/maxalign/${og_id}_cds.maxalign.fa.gz"
file_og_trimal="${dir_output_active}/trimal/${og_id}_cds.trimal.fa.gz"
file_og_clipkit="${dir_output_active}/clipkit/${og_id}_cds.clipkit.fa.gz"
file_og_clipkit_log="${dir_output_active}/clipkit_log/${og_id}_cds.clipkit.log"
file_og_iqtree_tree="${dir_output_active}/iqtree_tree/${og_id}_iqtree.nwk"
file_og_iqtree_generax_ufboot="${dir_output_active}/generax_ufboot_tree/${og_id}_generax.ufboot.nwk"
file_og_orthogroup_extraction_nwk="${dir_output_active}/orthogroup_extraction_nwk/${og_id}_orthogroup_extraction.nwk"
file_og_orthogroup_extraction_rooted_nwk="${dir_output_active}/orthogroup_extraction_rooted_nwk/${og_id}_orthogroup_extraction.rooted.nwk"
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
file_og_l1ou_fit_rdata="${dir_output_active}/l1ou_fit_rdata/${og_id}_l1ou.RData"
file_og_l1ou_fit_conv_rdata="${dir_output_active}/l1ou_fit_conv_rdata/${og_id}_l1ou.conv.RData"
file_og_l1ou_fit_tree="${dir_output_active}/l1ou_fit_tree/${og_id}_l1ou.tree.tsv"
file_og_l1ou_fit_regime="${dir_output_active}/l1ou_fit_regime/${og_id}_l1ou.regime.tsv"
file_og_l1ou_fit_leaf="${dir_output_active}/l1ou_fit_leaf/${og_id}_l1ou.leaf.tsv"
file_og_l1ou_fit_plot="${dir_output_active}/l1ou_fit_plot/${og_id}_l1ou.pdf"
file_og_l1ou_bootstrap="${dir_output_active}/l1ou_bootstrap/${og_id}_l1ou.bootstrap.tsv"
# Protein convergence analysis
file_og_iqtree_anc="${dir_output_active}/iqtree_anc/${og_id}_iqtree.anc.zip"
file_og_csubst_b="${dir_output_active}/csubst_b/${og_id}_csubst_b.tsv"
file_og_csubst_cb_2="${dir_output_active}/csubst_cb_2/${og_id}_csubst_cb_2.tsv"
file_og_csubst_cb_stats="${dir_output_active}/csubst_cb_stats/${og_id}_csubst_cb_stats.tsv"
file_og_csubst_scan="${dir_output_active}/csubst_scan/${og_id}_csubst_scan.tsv"
file_og_csubst_scan_units="${dir_output_active}/csubst_scan_units/${og_id}_csubst_scan_units.tsv"
file_og_csubst_scan_foreground_branch="${dir_output_active}/csubst_scan_foreground_branch/${og_id}_csubst_foreground_branch.txt"
file_og_csubst_scan_plot="${dir_output_active}/csubst_scan_plot/${og_id}_csubst_scan.tree_site.${csubst_scan_tree_site_plot_format}"
file_og_csubst_scan_log="${dir_output_active}/csubst_scan_log/${og_id}_csubst_scan.log"
file_og_iqtree_anc_provenance="${dir_output_active}/artifact_provenance/${og_id}.iqtree_anc.json"
file_og_csubst_provenance="${dir_output_active}/artifact_provenance/${og_id}.csubst.json"
file_og_csubst_scan_provenance="${dir_output_active}/artifact_provenance/${og_id}.csubst_scan.json"
file_og_summary_provenance="${dir_output_active}/artifact_provenance/${og_id}.summary_statistics.json"
file_og_tree_plot_provenance="${dir_output_active}/artifact_provenance/${og_id}.tree_plot.json"
if [[ ${csubst_max_arity} -gt 2 ]]; then
  for ((i = 3; i <= csubst_max_arity; i++)); do
    declare file_og_csubst_cb_${i}="${dir_output_active}/csubst_cb_${i}/${og_id}.csubst_cb_${i}.tsv"
  done
fi
# PGLS output
file_og_gene_pgls="${dir_output_active}/pgls_gene_tree/${og_id}_gene_tree_PGLS.tsv"
file_og_gene_pgls_plot="${dir_output_active}/pgls_gene_tree_plot/${og_id}_gene_PGLS.barplot.pdf"
file_og_species_nwkit_pgls="${dir_output_active}/pgls_species_nwkit/${og_id}_species_nwkit.pgls.tsv"
file_og_species_rphylopars_pgls="${dir_output_active}/pgls_species_rphylopars/${og_id}_species_rphylopars.pgls.tsv"
file_og_pgls_comparison="${dir_output_active}/pgls_comparison/${og_id}_pgls.comparison.tsv"
file_og_pgls_method_status="${dir_output_active}/pgls_method_status/${og_id}_pgls.method-status.tsv"
file_og_pgls_method_audit="${dir_output_active}/pgls_method_audit/${og_id}_pgls.method-audit.jsonl"
file_og_species_expression_summary="${dir_output_active}/pgls_species_expression_summary/${og_id}_species_expression.tsv"
file_og_species_expression_audit="${dir_output_active}/pgls_species_expression_audit/${og_id}_species_expression.audit.tsv"
file_og_species_response_tip_summary="${dir_output_active}/pgls_species_response_tip_summary/${og_id}_response-tip-summary.tsv"
file_og_species_response_sampling_covariance="${dir_output_active}/pgls_species_response_sampling_covariance/${og_id}_response-sampling-covariance.tsv"
file_og_species_predictor_tip_summary="${dir_output_active}/pgls_species_predictor_tip_summary/${og_id}_predictor-tip-summary.tsv"
file_og_species_predictor_sampling_covariance="${dir_output_active}/pgls_species_predictor_sampling_covariance/${og_id}_predictor-sampling-covariance.tsv"
file_og_rsc_status="${dir_output_active}/rsc_status/${og_id}_rsc.status.tsv"
file_og_rsc_pgls="${dir_output_active}/rsc_pgls/${og_id}_rsc.pgls.tsv"
file_og_rsc_reconciliation="${dir_output_active}/rsc_reconciliation/${og_id}_rsc.reconciliation.tsv"
file_og_rsc_gene_contrasts="${dir_output_active}/rsc_gene_contrasts/${og_id}_rsc.gene-contrasts.tsv"
file_og_rsc_species_contrasts="${dir_output_active}/rsc_species_contrasts/${og_id}_rsc.species-contrasts.tsv"
file_og_rsc_response_sampling_covariance="${dir_output_active}/rsc_response_sampling_covariance/${og_id}_rsc.response-sampling-covariance.tsv"
file_og_rsc_response_tip_summary="${dir_output_active}/rsc_response_tip_summary/${og_id}_rsc.response-tip-summary.tsv"
file_og_rsc_predictor_sampling_covariance="${dir_output_active}/rsc_predictor_sampling_covariance/${og_id}_rsc.predictor-sampling-covariance.tsv"
file_og_rsc_predictor_tip_summary="${dir_output_active}/rsc_predictor_tip_summary/${og_id}_rsc.predictor-tip-summary.tsv"
file_og_rsc_random_effects="${dir_output_active}/rsc_random_effects/${og_id}_rsc.random-effects.tsv"
file_og_rsc_sensitivity="${dir_output_active}/rsc_sensitivity/${og_id}_rsc.sensitivity.tsv"
file_og_rsc_trait_origins="${dir_output_active}/rsc_trait_origins/${og_id}_rsc.trait-origins.tsv"
file_og_rsc_audit="${dir_output_active}/rsc_audit/${og_id}_rsc.audit.jsonl"
file_og_rsc_log="${dir_output_active}/rsc_log/${og_id}_rsc.log"
file_og_expression_trait_pgls_log="${dir_output_active}/pgls_log/${og_id}_expression_trait_pgls.log"
file_og_expression_trait_pgls_provenance="${dir_output_active}/artifact_provenance/${og_id}.expression_trait_pgls.json"
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
file_og_primary_fasta="${file_og_cds_fasta}"
amas_data_type="dna"
if [[ "${input_sequence_mode}" == "protein" ]]; then
  file_og_primary_fasta="${file_og_pep_fasta}"
  amas_data_type="aa"
fi

# Define intermediate files for downstream analysis.
# These variables are updated via helper functions to keep routing logic explicit.






prepare_species_tree_pruned || true
set_default_analysis_files

memory_notung=$(gg_memory_fraction_gb "${GG_MEM_TOOL_GB}" 1 4)

echo "Checking parameter conflicts..."
if [[ ${run_trimal} -eq 1 && ${run_clipkit} -eq 1 ]]; then
  echo 'Both of ${run_trimal} and ${run_clipkit} are set to 1.'
  echo '${run_trimal} is deactivated. ${run_clipkit} is still active.'
  run_trimal=0
fi
if [[ ! -d "${dir_sp_gff}" ]] || [[ -z "$(find "${dir_sp_gff}" -mindepth 1 -maxdepth 1 -print -quit 2> /dev/null)" ]]; then
  if [[ ${run_collect_gff_info} -eq 1 ]]; then
    echo "\${run_collect_gff_info} is deactivated. Empty input: ${dir_sp_gff}"
    run_collect_gff_info=0
  fi
  if [[ ${run_scm_intron} -eq 1 ]]; then
    echo "\${run_scm_intron} is deactivated. Empty input: ${dir_sp_gff}"
    run_scm_intron=0
  fi
fi
if [[ -d "${dir_sp_expression}" ]] && [[ -n "$(find "${dir_sp_expression}" -mindepth 1 -maxdepth 1 -print -quit 2> /dev/null)" ]]; then
  echo "\${dir_sp_expression} is not empty. Continued: ${dir_sp_expression}"
else
  echo "\${dir_sp_expression} is empty: ${dir_sp_expression}"
  echo '${run_generate_expression_matrix}, ${run_tree_pruning}, and other options are deactivated.'
  run_tree_pruning=0
  run_generate_expression_matrix=0
  run_l1ou=0
fi
if [[ -d "${dir_sp_genome}" ]] && [[ -n "$(find "${dir_sp_genome}" -mindepth 1 -maxdepth 1 -print -quit 2> /dev/null)" ]]; then
  echo "\${dir_sp_genome} is not empty. Continued: ${dir_sp_genome}"
else
  echo "\${dir_sp_genome} is empty: ${dir_sp_genome}"
  echo '${run_extract_promoter_fasta} and ${run_fimo} are deactivated.'
  run_extract_promoter_fasta=0
  run_fimo=0
fi
echo "Checking preexisting tmp directory."
if [[ ! -e "${dir_tmp}" ]]; then
  echo "Making ${dir_tmp}"
  mkdir -p "${dir_tmp}"
fi
if [[ -d "${dir_output_active}/.gg_store" || -d "${dir_output_active}/.gg_archives" ]]; then
  materialize_args=(
    materialize-family
    "${gene_family_store_context_args[@]}"
    --family-id "${og_id}"
  )
  if [[ "${gene_family_output_storage}" == "zip" ]]; then
    materialize_args+=(
      --receipt "${gene_family_materialization_receipt}"
      --run-token "${gene_family_run_token}"
    )
  fi
  if ! python "${gene_family_store_script}" "${materialize_args[@]}"; then
    echo "Failed to materialize archived artifacts for ${og_id}." >&2
    exit 1
  fi
fi
adopt_historical_gene_family_outputs
cd "${dir_tmp}"
echo "Working at: $(pwd)"

trap cleanup_tmp_dir_on_normal_exit EXIT


ensure_species_fasta_sequence_store() {
  local sequence_kind=$1
  local source_dir=""
  local database=""
  local manifest=""
  local source_list="${dir_tmp}/.${sequence_kind}.sequence_store.sources.$$.tsv"
  local source_path=""
  local species=""
  case "${sequence_kind}" in
    cds)
      source_dir="${dir_sp_cds}"
      database="${file_species_cds_store_db}"
      manifest="${file_species_cds_store_manifest}"
      ;;
    protein)
      source_dir="${dir_sp_protein_input}"
      database="${file_species_protein_store_db}"
      manifest="${file_species_protein_store_manifest}"
      ;;
    *)
      echo "Invalid FASTA sequence-store kind: ${sequence_kind}" >&2
      return 1
      ;;
  esac
  : > "${source_list}"
  while IFS= read -r source_path; do
    species=$(gg_species_name_from_path "${source_path}")
    printf '%s\t%s\n' "${source_path}" "${species}" >> "${source_list}"
  done < <(gg_find_fasta_files "${source_dir}" 1)
  if [[ ! -s "${source_list}" ]]; then
    echo "No FASTA inputs were found for the ${sequence_kind} sequence store: ${source_dir}" >&2
    return 1
  fi
  ensure_dir "${dir_fasta_sequence_store}"
  if ! python "${gg_support_dir}/fasta_sequence_store.py" ensure \
    --database "${database}" \
    --manifest "${manifest}" \
    --source-list "${source_list}" \
    --digest-cache "${gg_workspace_dir}/.gg_cache/content_digests.sqlite3"
  then
    return 1
  fi
  rm -f -- "${source_list}"
}


# shellcheck shell=bash
# Sourced by gg_gene_evolution_core.sh.

task="Species tree availability check"
if [[ ! -s "${species_tree_pruned}" ]]; then
  echo "$(date): Warning: ${task}: species tree file was not found."
  echo "Missing: ${species_tree_pruned}"
else
  gg_step_skip "${task}"
fi

task="Query fasta generation"
query_fasta_needs_update=0
query_fasta_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.query_fasta.json"
  --step "query_fasta"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "query=${file_query_gene}"
  --output "query_fasta=${file_og_query_aa_fasta}"
  --parameter "genetic_code=${genetic_code}"
)
if [[ "$(head -c 1 "${file_query_gene}")" != ">" ]]; then
  ensure_species_fasta_sequence_store cds || exit 1
  query_fasta_provenance_args+=(--input "species_cds_index=${file_species_cds_store_manifest}")
  query_fasta_provenance_args+=(--parameter "fasta_sequence_store_schema=1")
fi
gg_artifact_prepare_stage query_fasta_needs_update run_extract_query_fasta "${query_fasta_provenance_args[@]}" || exit $?
if [[ ${query_fasta_needs_update} -eq 1 && ${run_extract_query_fasta} -eq 1 ]]; then
  gg_step_start "${task}"
  if [[ "$(head -c 1 "${file_query_gene}")" == ">" ]]; then
    seqtype=$(seqkit stats --tabular "${file_query_gene}" | awk 'NR>1 {print $3}')
    if [[ ${seqtype} == "DNA" ]]; then
      echo "DNA sequences were detected. The file will be treated as in-frame CDS sequences, translated into amino acids, and used as a ${query_blast_method} query: ${file_query_gene}"
      seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" --threads "${GG_TASK_CPUS}" "${file_query_gene}" > "${og_id}.query.aa.tmp.fasta"
      seqkit seq --threads "${GG_TASK_CPUS}" "${og_id}.query.aa.tmp.fasta" --out-file "${og_id}.query.aa.out.fa.gz"
      mv_out "${og_id}.query.aa.out.fa.gz" "${file_og_query_aa_fasta}"
      rm -f -- "${og_id}.query.aa.tmp.fasta"
    elif [[ ${seqtype} == "Protein" ]]; then
      echo "Amino acid sequences were detected. The file will be used as a ${query_blast_method} query: ${file_query_gene}"
      seqkit seq --threads "${GG_TASK_CPUS}" "${file_query_gene}" --out-file "${og_id}.query.aa.out.fa.gz"
      mv_out "${og_id}.query.aa.out.fa.gz" "${file_og_query_aa_fasta}"
    else
      echo "Unsupported sequence type '${seqtype}' in '${file_query_gene}'. Only \"DNA\" or \"Protein\" are allowed. Exiting."
      exit 1
    fi
  else
    echo "Gene IDs were detected. Extracting in-frame CDS sequences from species_cds: ${file_query_gene}"
    cp_out "${file_query_gene}" "${dir_output_active}/query_gene/$(basename "${file_query_gene}")"
    mapfile -t genes < <(sed -e '/^[[:space:]]*$/d' "${file_query_gene}")
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
    python "${gg_support_dir}/fasta_sequence_store.py" extract \
      --database "${file_species_cds_store_db}" \
      --pattern-file pattern.txt \
      --output "${og_id}.query.cds.fasta" \
      --ignore-case \
      --query-variants \
      --prefix-species
    gg_prepare_cds_fasta_stream "${GG_TASK_CPUS}" "${genetic_code}" < "${og_id}.query.cds.fasta" |
      sed -e '/^1 1$/d' -e 's/_frame=1[[:space:]]*//' \
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
      seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" --threads "${GG_TASK_CPUS}" "${og_id}.query.cds.2.fasta" > "${og_id}.query.aa.tmp.fasta"
      seqkit seq --threads "${GG_TASK_CPUS}" "${og_id}.query.aa.tmp.fasta" --out-file "${og_id}.query.aa.out.fa.gz"
      mv_out "${og_id}.query.aa.out.fa.gz" "${file_og_query_aa_fasta}"
      rm -f -- "${og_id}.query.aa.tmp.fasta"
      rm -f -- "${og_id}.query.cds.2.fasta"
    fi
  fi
  gg_artifact_record "${query_fasta_provenance_args[@]}"
else
  gg_step_skip "${task}"
fi

task="In-frame query BLAST (${query_blast_method})"
query_blast_needs_update=0
query_blast_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.query_blast.json"
  --step "query_blast"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "query_fasta=${file_og_query_aa_fasta}"
  --output "query_blast=${file_og_query_blast}"
  --parameter "method=${query_blast_method}"
  --parameter "genetic_code=${genetic_code}"
  --parameter "evalue=${query_blast_evalue}"
  --parameter "auto_evalue_cutoffs=${query_blast_auto_evalue_maxlen_cutoffs}"
)
if [[ "${mode_gene_evolution}" == "query2family" ]]; then
  ensure_species_fasta_sequence_store cds || exit 1
  query_blast_provenance_args+=(--input "species_cds_index=${file_species_cds_store_manifest}")
  query_blast_provenance_args+=(--parameter "fasta_sequence_store_schema=1")
  gg_artifact_prepare_stage query_blast_needs_update run_query_blast "${query_blast_provenance_args[@]}" || exit $?
fi
if [[ ${query_blast_needs_update} -eq 1 && ${run_query_blast} -eq 1 && "${mode_gene_evolution}" == "query2family" ]]; then
  gg_step_start "${task}"

  if [[ ${query_blast_method} == "tblastn" ]]; then
    if ! type makeblastdb > /dev/null 2>&1; then
      echo "makeblastdb was not found but query_blast_method=tblastn. Exiting."
      exit 1
    fi
    if ! type tblastn > /dev/null 2>&1; then
      echo "tblastn was not found but query_blast_method=tblastn. Exiting."
      exit 1
    fi
  elif [[ ${query_blast_method} == "diamond" ]]; then
    if ! type diamond > /dev/null 2>&1; then
      echo "diamond was not found but query_blast_method=diamond. Exiting."
      exit 1
    fi
    echo "DIAMOND mode selected: species CDS will be translated to proteins because diamond makedb/blastp use protein reference databases."
  fi

  export BLASTDB_LMDB_MAP_SIZE=100000000
  check_species_cds "${gg_workspace_dir}"
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
    cds_spp+=("$(gg_species_name_from_path "${cds_file}")")
  done
  mapfile -t cds_spp < <(printf "%s\n" "${cds_spp[@]}" | sort -u)
  db_build_jobs=${GG_TASK_CPUS}
  [[ ${db_build_jobs} -gt ${#cds_spp[@]} ]] && db_build_jobs=${#cds_spp[@]}
  [[ ${db_build_jobs} -lt 1 ]] && db_build_jobs=1
  db_threads_per_job=$((GG_TASK_CPUS / db_build_jobs))
  [[ ${db_threads_per_job} -lt 1 ]] && db_threads_per_job=1
  species_cds_digest_table="${dir_tmp}/species_cds.sequence_digests.$$.tsv"
  python - "${file_species_cds_store_manifest}" > "${species_cds_digest_table}" <<'PY'
import json
import pathlib
import sys

payload = json.loads(pathlib.Path(sys.argv[1]).read_text(encoding="utf-8"))
for entry in payload.get("sources", []):
    print(f"{entry['species']}\t{entry['sha256']}")
PY
  filter_translated_fasta_for_diamond() {
    awk '
      BEGIN {
        seen = 0
        dropped = 0
        header = ""
        seq = ""
      }
      /^>/ {
        if (seen) {
          if (seq != "") {
            print header
            print seq
          } else {
            dropped++
          }
        }
        header = $0
        seq = ""
        seen = 1
        next
      }
      {
        gsub(/\*/, "", $0)
        if ($0 != "") {
          seq = seq $0
        }
      }
      END {
        if (seen) {
          if (seq != "") {
            print header
            print seq
          } else {
            dropped++
          }
        }
        if (dropped > 0) {
          printf("Dropped %d translated protein records with empty sequence after stop-codon removal.\n", dropped) > "/dev/stderr"
        }
      }
    '
  }
  for sp in "${cds_spp[@]}"; do
    wait_until_jobn_le "${db_build_jobs}"
    echo "sp: ${sp}"
    sp_cds_candidates=()
    for cds_file in "${cds_files[@]}"; do
      if [[ "$(gg_species_name_from_path "${cds_file}")" == "${sp}" ]]; then
        sp_cds_candidates+=("${cds_file}")
      fi
    done
    if [[ ${#sp_cds_candidates[@]} -eq 0 ]]; then
      echo "No CDS file was found for species: ${sp}. Skipping."
      continue
    fi
    mapfile -t sp_cds_candidates < <(printf "%s\n" "${sp_cds_candidates[@]}" | sort)
    sp_cds=${sp_cds_candidates[0]}
    sp_cds_blastdb="${dir_sp_blastdb}/$(basename "${sp_cds}")"
    sp_cds_source_digest=$(awk -F '\t' -v species="${sp}" '$1 == species {print $2; exit}' "${species_cds_digest_table}")
    if [[ ! "${sp_cds_source_digest}" =~ ^[0-9a-f]{64}$ ]]; then
      echo "No indexed content digest was found for species CDS: ${sp_cds}. Exiting." >&2
      exit 1
    fi
    db_source_signature="source_sha256=${sp_cds_source_digest};method=${query_blast_method}"
    if [[ ${query_blast_method} == "diamond" ]]; then
      db_source_signature="${db_source_signature};genetic_code=${genetic_code};translation_filter=1"
    fi
    db_signature_file="${sp_cds_blastdb}.${query_blast_method}.build.signature"
    db_files+=("${sp_cds_blastdb}")
    if [[ ${query_blast_method} == "tblastn" ]]; then
      echo "makeblastdb input CDS file: ${sp_cds}"
      echo "makeblastdb output database file: ${sp_cds_blastdb}"
      if [[ ! -e "${sp_cds_blastdb}".nhr || ! -e "${sp_cds_blastdb}".nin || ! -e "${sp_cds_blastdb}".nsq || ! -e "${sp_cds_blastdb}".ndb || ! -s "${db_signature_file}" || "$(< "${db_signature_file}")" != "${db_source_signature}" ]]; then
        db_lock_file="${sp_cds_blastdb}.tblastn.build.lock"
        (
          if ! gg_shared_lock_acquire "${db_lock_file}" "TBLASTN database build (${sp})"; then
            exit 1
          fi
          gg_shared_lock_start_heartbeat "${db_lock_file}"
          heartbeat_pid=${GG_SHARED_LOCK_HEARTBEAT_PID:-}
          cleanup_tblastn_db_lock() {
            gg_shared_lock_stop_heartbeat "${heartbeat_pid}"
            gg_shared_lock_release "${db_lock_file}"
          }
          trap cleanup_tblastn_db_lock EXIT
          if [[ ! -e "${sp_cds_blastdb}".nhr || ! -e "${sp_cds_blastdb}".nin || ! -e "${sp_cds_blastdb}".nsq || ! -e "${sp_cds_blastdb}".ndb || ! -s "${db_signature_file}" || "$(< "${db_signature_file}")" != "${db_source_signature}" ]]; then
            if zgrep -q -e "^>.*[[:blank:]]" "${sp_cds}"; then
              echo "Space is detected. Please remove all annotation info after spaces in sequence names. Exiting: ${sp_cds}"
              exit 1
            fi
            if zgrep -q -e "^>.*[|]" "${sp_cds}"; then
              echo "Bar (|) is detected. Bars in sequence names will be replaced with underlines (_): ${sp_cds}"
            fi
            echo "Generating BLAST database: ${sp_cds}"
            echo "Generating BLAST database: ${sp_cds}" >&2
            if [[ ${sp_cds} == *.gz ]]; then
              seqkit seq --threads "${db_threads_per_job}" "${sp_cds}" | makeblastdb -dbtype nucl -title "${sp_cds}" -out "${sp_cds_blastdb}"
            else
              makeblastdb -dbtype nucl -in "${sp_cds}" -out "${sp_cds_blastdb}"
            fi
            printf '%s\n' "${db_source_signature}" > "${db_signature_file}"
          fi
        ) &
        gg_background_register "$!"
      fi
    elif [[ ${query_blast_method} == "diamond" ]]; then
      sp_cds_diamond_fasta="${sp_cds_blastdb}.diamond.fasta"
      echo "diamond input CDS file: ${sp_cds}"
      echo "diamond translated protein file: ${sp_cds_diamond_fasta}"
      echo "diamond database file: ${sp_cds_blastdb}.dmnd"
      if [[ ! -e "${sp_cds_blastdb}".dmnd || ! -s "${db_signature_file}" || "$(< "${db_signature_file}")" != "${db_source_signature}" ]]; then
        db_lock_file="${sp_cds_blastdb}.diamond.build.lock"
        (
          if ! gg_shared_lock_acquire "${db_lock_file}" "DIAMOND database build (${sp})"; then
            exit 1
          fi
          gg_shared_lock_start_heartbeat "${db_lock_file}"
          heartbeat_pid=${GG_SHARED_LOCK_HEARTBEAT_PID:-}
          cleanup_diamond_db_lock() {
            gg_shared_lock_stop_heartbeat "${heartbeat_pid}"
            gg_shared_lock_release "${db_lock_file}"
          }
          trap cleanup_diamond_db_lock EXIT
          if [[ ! -e "${sp_cds_blastdb}".dmnd || ! -s "${db_signature_file}" || "$(< "${db_signature_file}")" != "${db_source_signature}" ]]; then
            if zgrep -q -e "^>.*[[:blank:]]" "${sp_cds}"; then
              echo "Space is detected. Please remove all annotation info after spaces in sequence names. Exiting: ${sp_cds}"
              exit 1
            fi
            if zgrep -q -e "^>.*[|]" "${sp_cds}"; then
              echo "Bar (|) is detected. Bars in sequence names will be replaced with underlines (_): ${sp_cds}"
            fi
            echo "Generating DIAMOND database: ${sp_cds}"
            echo "Generating DIAMOND database: ${sp_cds}" >&2
            if [[ ${sp_cds} == *.gz ]]; then
              seqkit seq --remove-gaps --threads 1 "${sp_cds}" |
                seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" --threads "${db_threads_per_job}" |
                filter_translated_fasta_for_diamond \
                  > "${sp_cds_diamond_fasta}"
            else
              seqkit seq --remove-gaps --threads 1 "${sp_cds}" |
                seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" --threads "${db_threads_per_job}" |
                filter_translated_fasta_for_diamond \
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
            if ! diamond makedb --in "${sp_cds_diamond_fasta}" --db "${sp_cds_blastdb}" --threads "${db_threads_per_job}"; then
              echo "diamond makedb failed for ${sp_cds}. Exiting."
              exit 1
            fi
            rm -f -- "${sp_cds_diamond_fasta}"
            printf '%s\n' "${db_source_signature}" > "${db_signature_file}"
          fi
        ) &
        gg_background_register "$!"
      fi
    fi
  done
  wait_for_background_jobs
  rm -f -- "${species_cds_digest_table}"
  echo "db_files: ${db_files[*]}"
  query_aa_local="${og_id}.query.aa.tmp.for_blast.fasta"
  seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_query_aa_fasta}" --out-file "${query_aa_local}"
  resolve_query_blast_evalue "${query_blast_evalue}" "${query_aa_local}" "${query_blast_auto_evalue_maxlen_cutoffs}"
  if [[ "$(printf '%s' "${query_blast_evalue}" | tr '[:upper:]' '[:lower:]')" == "auto" ]]; then
    echo "query BLAST auto E-value: query_count=${query_blast_query_num_seqs} min_aa_len=${query_blast_query_min_aa_len} avg_aa_len=${query_blast_query_avg_aa_len} max_aa_len=${query_blast_query_max_aa_len}"
    echo "query BLAST auto E-value: cutoffs=${query_blast_auto_evalue_maxlen_cutoffs} effective_query_blast_evalue=${effective_query_blast_evalue}"
  else
    echo "query BLAST E-value: effective_query_blast_evalue=${effective_query_blast_evalue}"
  fi

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
      -evalue "${effective_query_blast_evalue}" \
      -max_target_seqs 50000 \
      -outfmt "6 ${outfmt}" \
      -num_threads "${GG_TASK_CPUS}"; then
      echo "tblastn failed. Exiting."
      exit 1
    fi
  elif [[ ${query_blast_method} == "diamond" ]]; then
    echo "Running diamond blastp."
    rm -f -- blast_out.tsv
    touch blast_out.tsv
    diamond_search_jobs=${GG_TASK_CPUS}
    [[ ${diamond_search_jobs} -gt ${#db_files[@]} ]] && diamond_search_jobs=${#db_files[@]}
    [[ ${diamond_search_jobs} -lt 1 ]] && diamond_search_jobs=1
    diamond_threads_per_job=$((GG_TASK_CPUS / diamond_search_jobs))
    [[ ${diamond_threads_per_job} -lt 1 ]] && diamond_threads_per_job=1
    for db_file in "${db_files[@]}"; do
      if [[ ! -e "${db_file}".dmnd ]]; then
        echo "DIAMOND database file is missing: ${db_file}.dmnd. Exiting."
        exit 1
      fi
      tmp_diamond_out="$(basename "${db_file}").diamond.out.tsv"
      wait_until_jobn_le "${diamond_search_jobs}"
      (
        if ! diamond blastp \
          --query "${query_aa_local}" \
          --db "${db_file}" \
          --out "${tmp_diamond_out}" \
          --evalue "${effective_query_blast_evalue}" \
          --max-target-seqs 50000 \
          --threads "${diamond_threads_per_job}" \
          --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen
        then
          echo "diamond blastp failed for database: ${db_file}." >&2
          exit 1
        fi
      ) &
      gg_background_register "$!"
    done
    wait_for_background_jobs
    for db_file in "${db_files[@]}"; do
      tmp_diamond_out="$(basename "${db_file}").diamond.out.tsv"
      if [[ -s "${tmp_diamond_out}" ]]; then
        awk -F '\t' 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,"0/1",$13,$14}' "${tmp_diamond_out}" >> blast_out.tsv
      fi
      rm -f -- "${tmp_diamond_out}"
    done
  fi
  rm -f -- "${query_aa_local}"

  python "${gg_support_dir}/annotate_blast_coverage.py" \
    --in blast_out.tsv \
    --ncpu "${GG_TASK_CPUS}" \
    --outfmt-columns "${outfmt}" \
    --frame-filter "0/1" \
    --out blast_out_inframe.tmp3.tsv

  if [[ -s blast_out_inframe.tmp3.tsv ]]; then
    mv_out blast_out_inframe.tmp3.tsv "${file_og_query_blast}"
  else
    echo "No query BLAST hits were detected after in-frame filtering. Exiting."
    exit 1
  fi
  gg_artifact_record "${query_blast_provenance_args[@]}"
else
  gg_step_skip "${task}"
fi

task="Fasta generation"
primary_fasta_needs_update=0
primary_fasta_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.primary_fasta.json"
  --step "primary_fasta"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --output "primary_fasta=${file_og_primary_fasta}"
  --parameter "mode=${mode_gene_evolution}"
  --parameter "input_sequence_mode=${input_sequence_mode}"
  --parameter "genetic_code=${genetic_code}"
  --parameter "query_blast_coverage=${query_blast_coverage}"
  --parameter "max_query_hits=${max_num_gene_blast_hit_retrieval}"
)
if [[ "${mode_gene_evolution}" == "orthogroup" ]]; then
  primary_fasta_provenance_args+=(--input "orthogroups=${file_og}")
else
  primary_fasta_provenance_args+=(--input "query_blast=${file_og_query_blast}")
fi
if [[ "${input_sequence_mode}" == "protein" ]] && species_protein_input_has_files; then
  ensure_species_fasta_sequence_store protein || exit 1
  primary_fasta_provenance_args+=(--input "species_protein_index=${file_species_protein_store_manifest}")
else
  ensure_species_fasta_sequence_store cds || exit 1
  primary_fasta_provenance_args+=(--input "species_cds_index=${file_species_cds_store_manifest}")
  if [[ -s "${file_species_genetic_code}" ]]; then
    primary_fasta_provenance_args+=(--input "species_genetic_code=${file_species_genetic_code}")
  fi
fi
primary_fasta_provenance_args+=(--parameter "fasta_sequence_store_schema=1")
gg_artifact_prepare_stage primary_fasta_needs_update run_extract_primary_fasta "${primary_fasta_provenance_args[@]}" || exit $?
if [[ ${primary_fasta_needs_update} -eq 1 && ${run_extract_primary_fasta} -eq 1 ]]; then
  gg_step_start "${task}"

  if [[ "${mode_gene_evolution}" == "orthogroup" ]]; then
    genes=()
    read -r -a genes <<< "$(awk -v og="${og_id}" '$1==og {$1=""; sub(/^[[:space:]]+/, "", $0); gsub(",", "", $0); gsub(/\t/, " ", $0); sub(/[[:space:]]*$/, "", $0); gsub(/\047|"/, "", $0); print; exit}' "${file_og}")"
  elif [[ "${mode_gene_evolution}" == "query2family" ]]; then
    python "${gg_support_dir}/extract_gene_id_from_blast_table.py" \
      --infile "${file_og_query_blast}" \
      --outfile gene_id_list.txt \
      --min_query_blast_coverage "${query_blast_coverage}" \
      --max_num_gene_blast_hit_retrieval "${max_num_gene_blast_hit_retrieval}"
    mapfile -t genes < gene_id_list.txt
  fi
  if [[ -e pattern.txt ]]; then
    rm -f -- pattern.txt
  fi
  touch pattern.txt
  for gene in "${genes[@]}"; do
    echo "${gene}" >> pattern.txt
  done

  num_gene=${#genes[@]}
  if [[ "${input_sequence_mode}" == "protein" ]] && species_protein_input_has_files; then
    check_species_protein_dir "${dir_sp_protein_input}"
    check_if_species_files_unique "${dir_sp_protein_input}"
    if [[ -s "${file_species_genetic_code}" ]]; then
      echo "species_genetic_code.tsv is ignored because species_protein inputs are provided: ${file_species_genetic_code}"
    fi
    rm -f -- "${og_id}.pep.fasta"
    python "${gg_support_dir}/fasta_sequence_store.py" extract \
      --database "${file_species_protein_store_db}" \
      --pattern-file pattern.txt \
      --output "${og_id}.pep.fasta" \
      --require-all
    seqkit replace --pattern " .*" --replacement "" --ignore-case --threads "${GG_TASK_CPUS}" "${og_id}.pep.fasta" |
      seqkit replace --pattern "\+" --replacement "_" --ignore-case --threads "${GG_TASK_CPUS}" |
      sed -e '/^1 1$/d' -e 's/_frame=1[[:space:]]*//' \
        > "${og_id}.pep.2.fasta"
    fasta_genes=()
    mapfile -t fasta_genes < <(awk '/^>/ {sub(/^>/, "", $0); print}' "${og_id}.pep.2.fasta")
    num_seq=${#fasta_genes[@]}
    echo "Number of genes in the orthogroup or BLAST hit: ${num_gene}"
    echo "Number of sequences in the protein fasta: ${num_seq}"
    if [[ ${num_gene} -eq ${num_seq} ]]; then
      echo "Number of genes and sequences matched. Protein fasta generation completed!"
      seqkit seq --threads "${GG_TASK_CPUS}" "${og_id}.pep.2.fasta" --out-file "${og_id}.pep.out.fa.gz"
      mv_out "${og_id}.pep.out.fa.gz" "${file_og_pep_fasta}"
    else
      echo "Number of genes and sequences did not match."
      echo "Genes in the orthogroup or BLAST hit:"
      printf '%s\n' "${genes[@]}"
      echo ""
      echo "Genes in the generated protein FASTA:"
      printf '%s\n' "${fasta_genes[@]}" | sort | tr '\n' ' '
      echo ""
      echo "There may be duplicated or missing sequences."
      echo "Exiting."
      exit 1
    fi
  else
    if [[ ! -s "${file_og_cds_fasta}" ]]; then
      if [[ -e "${og_id}.cds.fasta" ]]; then
        rm -f -- "${og_id}.cds.fasta"
      fi
      python "${gg_support_dir}/fasta_sequence_store.py" extract \
        --database "${file_species_cds_store_db}" \
        --pattern-file pattern.txt \
        --output "${og_id}.cds.fasta" \
        --require-all

      seqkit replace --pattern "X" --replacement "N" --by-seq --ignore-case --threads "${GG_TASK_CPUS}" "${og_id}.cds.fasta" |
        seqkit replace --pattern " .*" --replacement "" --ignore-case --threads "${GG_TASK_CPUS}" |
        seqkit replace --pattern "\+" --replacement "_" --ignore-case --threads "${GG_TASK_CPUS}" |
        cdskit pad --codon_table "${genetic_code}" |
        sed -e '/^1 1$/d' -e 's/_frame=1[[:space:]]*//' \
          > "${og_id}.cds.2.fasta"

      fasta_genes=()
      mapfile -t fasta_genes < <(awk '/^>/ {sub(/^>/, "", $0); print}' "${og_id}.cds.2.fasta")
      num_seq=${#fasta_genes[@]}
      echo "Number of genes in the orthogroup or BLAST hit: ${num_gene}"
      echo "Number of sequences in the fasta: ${num_seq}"
      if [[ ${num_gene} -eq ${num_seq} ]]; then
        echo "Number of genes and sequences matched. Fasta generation completed!"
        seqkit seq --threads "${GG_TASK_CPUS}" "${og_id}.cds.2.fasta" --out-file "${og_id}.cds.out.fa.gz"
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
    fi
    if [[ "${input_sequence_mode}" == "protein" ]]; then
      gg_prepare_species_genetic_code_table "${dir_sp_cds}" "${genetic_code}" "${file_species_genetic_code_resolved}" "${file_species_genetic_code}"
      translate_orthogroup_cds_to_protein_fasta "${file_og_cds_fasta}" "${file_og_pep_fasta}" "${file_species_genetic_code_resolved}"
    fi
  fi
  gg_artifact_record "${primary_fasta_provenance_args[@]}"
else
  gg_step_skip "${task}"
fi
# shellcheck shell=bash
# Sourced by gg_gene_evolution_core.sh.

task="Protein RPS-BLAST"
disable_if_no_input_file "run_rps_blast" "${file_og_primary_fasta}"
rpsblast_needs_update=0
rpsblast_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.rpsblast.json"
  --step "rpsblast"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "primary_fasta=${file_og_primary_fasta}"
  --output "rpsblast=${file_og_rpsblast}"
  --parameter "input_sequence_mode=${input_sequence_mode}"
  --parameter "genetic_code=${genetic_code}"
  --parameter "evalue=0.01"
)
gg_artifact_prepare_stage rpsblast_needs_update run_rps_blast "${rpsblast_provenance_args[@]}" || exit $?
if [[ ${rpsblast_needs_update} -eq 1 && ${run_rps_blast} -eq 1 ]]; then
  gg_step_start "${task}"
  if ! dir_rpsblastdb=$(ensure_pfam_le_db "${gg_workspace_dir}"); then
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

  if [[ "${input_sequence_mode}" == "protein" ]]; then
    seqkit seq --remove-gaps --threads "${GG_TASK_CPUS}" "${file_og_pep_fasta}" > ungapped_translated_cds.fas
  else
    seqkit seq --remove-gaps --threads "${GG_TASK_CPUS}" "${file_og_cds_fasta}" |
      seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" --threads "${GG_TASK_CPUS}" \
        > ungapped_translated_cds.fas
  fi

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
    -num_threads "${GG_TASK_CPUS}"; then
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
  gg_artifact_record "${rpsblast_provenance_args[@]}"
else
  gg_step_skip "${task}"
fi

task="Gene trait extraction from gff files"
disable_if_no_input_file "run_collect_gff_info" "${file_og_primary_fasta}"
gff_info_needs_update=0
gff_info_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.gff_info.json"
  --step "gff_info"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "primary_fasta=${file_og_primary_fasta}"
  --input "species_gff=${dir_sp_gff}"
  --output "gff_info=${file_og_gff_info}"
  --parameter "feature=CDS"
  --parameter "multiple_hits=longest"
)
gg_artifact_prepare_stage gff_info_needs_update run_collect_gff_info "${gff_info_provenance_args[@]}" || exit $?
if [[ ${gff_info_needs_update} -eq 1 && ${run_collect_gff_info} -eq 1 ]]; then
  gg_step_start "${task}"
  if [[ -e gff2genestat.tsv ]]; then
    rm -f -- gff2genestat.tsv
  fi
  seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_primary_fasta}" --out-file "${og_id}.gff2genestat_input.fasta"

  python "${gg_support_dir}/gff2genestat.py" \
    --dir_gff "${dir_sp_gff}" \
    --feature "CDS" \
    --multiple_hits "longest" \
    --seqfile "${og_id}.gff2genestat_input.fasta" \
    --ncpu "${GG_TASK_CPUS}" \
    --outfile gff2genestat.tsv
  rm -f -- "${og_id}.gff2genestat_input.fasta"

  if [[ -s gff2genestat.tsv ]]; then
    mv_out gff2genestat.tsv "${file_og_gff_info}"
  fi
  gg_artifact_record "${gff_info_provenance_args[@]}"
else
  gg_step_skip "${task}"
fi

task="UniProt annotation (${uniprot_annotation_method})"
disable_if_no_input_file "run_uniprot_annotation" "${file_og_primary_fasta}"
uniprot_annotation_needs_update=0
uniprot_annotation_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.uniprot_annotation.json"
  --step "uniprot_annotation"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "primary_fasta=${file_og_primary_fasta}"
  --output "uniprot_annotation=${file_og_uniprot_annotation}"
  --parameter "method=${uniprot_annotation_method}"
  --parameter "input_sequence_mode=${input_sequence_mode}"
  --parameter "genetic_code=${genetic_code}"
  --parameter "evalue=1e-2"
)
gg_artifact_prepare_stage uniprot_annotation_needs_update run_uniprot_annotation "${uniprot_annotation_provenance_args[@]}" || exit $?
if [[ ${uniprot_annotation_needs_update} -eq 1 && ${run_uniprot_annotation} -eq 1 ]]; then
  gg_step_start "${task}"

  if [[ "${input_sequence_mode}" == "protein" ]]; then
    seqkit seq --remove-gaps --only-id --threads "${GG_TASK_CPUS}" "${file_og_pep_fasta}" > uniprot.query.pep.fas
  else
    seqkit seq --remove-gaps --only-id --threads "${GG_TASK_CPUS}" "${file_og_cds_fasta}" |
      seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" --threads "${GG_TASK_CPUS}" \
        > uniprot.query.pep.fas
  fi

  if [[ "${uniprot_annotation_method}" == "blastp" ]]; then
    if ! uniprot_db_prefix=$(ensure_uniprot_sprot_blast_db "${gg_workspace_dir}"); then
      echo "Failed to prepare UniProt Swiss-Prot BLASTP DB. Exiting."
      exit 1
    fi
    if ! validate_uniprot_sprot_db_prefix "${uniprot_db_prefix}" "blastp"; then
      echo "Invalid UniProt Swiss-Prot BLASTP DB prefix. Exiting."
      exit 1
    fi

    blastp \
      -query uniprot.query.pep.fas \
      -num_threads "${GG_TASK_CPUS}" \
      -db "${uniprot_db_prefix}" \
      -out uniprot.search.tsv \
      -outfmt "6 qseqid sseqid pident length evalue bitscore qlen" \
      -max_target_seqs 1 \
      -evalue 1e-2
  else
    if ! uniprot_db_prefix=$(ensure_uniprot_sprot_mmseqs_db "${gg_workspace_dir}"); then
      echo "Failed to prepare UniProt Swiss-Prot MMseqs2 DB. Exiting."
      exit 1
    fi
    if ! validate_uniprot_sprot_db_prefix "${uniprot_db_prefix}" "mmseqs2"; then
      echo "Invalid UniProt Swiss-Prot MMseqs2 DB prefix. Exiting."
      exit 1
    fi

    mmseqs createdb "uniprot.query.pep.fas" "uniprot.queryDB"
    mmseqs search "uniprot.queryDB" "${uniprot_db_prefix}.mmseqs" "uniprot.resultDB" "tmp_mmseqs2_uniprot" \
      --threads "${GG_TASK_CPUS}" \
      --split-memory-limit "$(gg_memory_fraction_gb "${GG_MEM_TOOL_GB}" 3 4)G" \
      --max-seqs 1 \
      -e 1e-2 \
      -s 7.5
    mmseqs convertalis "uniprot.queryDB" "${uniprot_db_prefix}.mmseqs" "uniprot.resultDB" "uniprot.search.tsv" \
      --threads "${GG_TASK_CPUS}" \
      --format-output "query,target,pident,alnlen,evalue,bits,qlen"
    rm -f -- uniprot.queryDB* uniprot.resultDB*
    rm -rf -- "tmp_mmseqs2_uniprot"
  fi

  uniprot_meta_tsv=""
  if ! uniprot_meta_tsv=$(ensure_uniprot_sprot_metadata_tsv "${gg_workspace_dir}" "${uniprot_db_prefix}" 2>/dev/null); then
    echo "Warning: UniProt Swiss-Prot metadata TSV is unavailable for prefix: ${uniprot_db_prefix}" >&2
    uniprot_meta_tsv=""
  fi

  python "${gg_support_dir}/reformat_uniprot_diamond.py" \
    --diamond_tsv uniprot.search.tsv \
    --query_fasta uniprot.query.pep.fas \
    --uniprot_fasta "${uniprot_db_prefix}.pep" \
    --uniprot_meta_tsv "${uniprot_meta_tsv}" \
    --outfile uniprot.annotation.tsv

  cp_out uniprot.annotation.tsv "${file_og_uniprot_annotation}"
  gg_artifact_record "${uniprot_annotation_provenance_args[@]}"
else
  gg_step_skip "${task}"
fi

task="cdskit localize"
disable_if_no_input_file "run_cdskit_localize" "${file_og_primary_fasta}"
cdskit_localize_needs_update=0
cdskit_localize_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.cdskit_localize.json"
  --step "cdskit_localize"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "primary_fasta=${file_og_primary_fasta}"
  --output "cdskit_localize=${file_og_cdskit_localize}"
  --parameter "model=${cdskit_localize_model}"
  --parameter "organism_group=${cdskit_localize_organism_group}"
  --parameter "include_features=${cdskit_localize_include_features}"
  --parameter "input_sequence_mode=${input_sequence_mode}"
  --parameter "genetic_code=${genetic_code}"
  --parameter "busco_lineage=${busco_lineage}"
)
gg_artifact_prepare_stage cdskit_localize_needs_update run_cdskit_localize "${cdskit_localize_provenance_args[@]}" || exit $?
if [[ ${cdskit_localize_needs_update} -eq 1 && ${run_cdskit_localize} -eq 1 ]]; then
  gg_step_start "${task}"

  cdskit_localize_species_dir="${dir_sp_cds}"
  if [[ "${input_sequence_mode}" == "protein" ]] && species_protein_input_has_files; then
    cdskit_localize_species_dir="${dir_sp_protein_input}"
  fi
  cdskit_localize_species_names=()
  mapfile -t cdskit_localize_species_names < <(gg_species_names_from_fasta_dir "${cdskit_localize_species_dir}")
  cdskit_localize_group_resolved=$(
    gg_resolve_cdskit_localize_organism_group \
      "${cdskit_localize_organism_group}" \
      "${gg_workspace_dir}" \
      "${busco_lineage}" \
      "${cdskit_localize_species_names[@]}"
  )

  if [[ "${input_sequence_mode}" == "protein" ]]; then
    seqkit seq --remove-gaps --only-id --threads "${GG_TASK_CPUS}" "${file_og_pep_fasta}" > "cdskit_localize.input.pep.fasta"
    cdskit_localize_input="cdskit_localize.input.pep.fasta"
    cdskit_localize_seqtype="protein"
  else
    gg_prepare_cdskit_localize_cds_input \
      "${file_og_cds_fasta}" \
      "cdskit_localize.input.cds.fasta" \
      "${GG_TASK_CPUS}" \
      "${genetic_code}"
    cdskit_localize_input="cdskit_localize.input.cds.fasta"
    cdskit_localize_seqtype="dna"
  fi

  gg_run_cdskit_localize \
    "${cdskit_localize_input}" \
    "${cdskit_localize_seqtype}" \
    "cdskit_localize.tsv" \
    "${cdskit_localize_model}" \
    "${cdskit_localize_group_resolved}" \
    "${cdskit_localize_include_features}" \
    "${cdskit_localize_no_model_download}" \
    "${GG_TASK_CPUS}" \
    "${genetic_code}"
  if [[ -s "cdskit_localize.tsv" ]]; then
    mv_out "cdskit_localize.tsv" "${file_og_cdskit_localize}"
  fi
  rm -f -- "cdskit_localize.input.cds.fasta" "cdskit_localize.input.pep.fasta"
  gg_artifact_record "${cdskit_localize_provenance_args[@]}"
else
  gg_step_skip "${task}"
fi

task="In-frame mafft alignment"
disable_if_no_input_file "run_mafft" "${file_og_primary_fasta}"
mafft_needs_update=0
mafft_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.mafft.json"
  --step "mafft"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "primary_fasta=${file_og_primary_fasta}"
  --output "mafft=${file_og_mafft}"
  --parameter "input_sequence_mode=${input_sequence_mode}"
  --parameter "genetic_code=${genetic_code}"
)
gg_artifact_prepare_stage mafft_needs_update run_mafft "${mafft_provenance_args[@]}" || exit $?
if [[ ${mafft_needs_update} -eq 1 && ${run_mafft} -eq 1 ]]; then
  gg_step_start "${task}"

  if [[ "${input_sequence_mode}" == "protein" ]]; then
    seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_pep_fasta}" --out-file tmp.pep.input.fasta
    mafft \
      --auto \
      --amino \
      --thread "${GG_TASK_CPUS}" \
      --quiet \
      tmp.pep.input.fasta \
      > "${og_id}.cds.aln.fasta"
  else
    seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_cds_fasta}" --out-file tmp.cds.input.fasta
    cdskit mask --seq_file tmp.cds.input.fasta --codon_table "${genetic_code}" --out_file tmp.cds.fasta

    seqkit translate \
      --allow-unknown-codon \
      --transl-table "${genetic_code}" \
      --threads "${GG_TASK_CPUS}" \
      tmp.cds.fasta \
      > tmp.pep.fasta

    mafft \
      --auto \
      --amino \
      --thread "${GG_TASK_CPUS}" \
      --quiet \
      tmp.pep.fasta \
      > tmp.pep.aln.fasta

    cdskit backalign \
      --seq_file tmp.cds.fasta \
      --aa_aln tmp.pep.aln.fasta \
      --codon_table "${genetic_code}" \
      --out_file "${og_id}.cds.aln.fasta"
  fi

  seqkit seq --threads "${GG_TASK_CPUS}" "${og_id}.cds.aln.fasta" --out-file "${og_id}.cds.aln.out.fa.gz"
  mv_out "${og_id}.cds.aln.out.fa.gz" "${file_og_mafft}"
  rm -f -- tmp.cds.input.fasta tmp.cds.fasta tmp.pep.fasta tmp.pep.aln.fasta tmp.pep.input.fasta
  gg_artifact_record "${mafft_provenance_args[@]}"
else
  gg_step_skip "${task}"
fi

task="AMAS for original alignment"
disable_if_no_input_file "run_amas_original" "${file_og_untrimmed_aln_analysis}"
amas_original_needs_update=0
amas_original_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.amas_original.json"
  --step "amas_original"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "alignment=${file_og_untrimmed_aln_analysis}"
  --output "amas=${file_og_amas_original}"
  --parameter "data_type=${amas_data_type}"
)
gg_artifact_prepare_stage amas_original_needs_update run_amas_original "${amas_original_provenance_args[@]}" || exit $?
if [[ ${amas_original_needs_update} -eq 1 && ${run_amas_original} -eq 1 ]]; then
  gg_step_start "${task}"
  seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_untrimmed_aln_analysis}" --out-file "${og_id}.amas.original.input.fasta"

  AMAS.py summary \
    --in-format fasta \
    --data-type "${amas_data_type}" \
    --in-files "${og_id}.amas.original.input.fasta"

  mv_out summary.txt "${file_og_amas_original}"
  rm -f -- "${og_id}.amas.original.input.fasta"
  gg_artifact_record "${amas_original_provenance_args[@]}"
else
  gg_step_skip "${task}"
fi

task="MaxAlign"
disable_if_no_input_file "run_maxalign" "${file_og_untrimmed_aln_analysis}"
maxalign_needs_update=0
maxalign_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.maxalign.json"
  --step "maxalign"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "untrimmed_alignment=${file_og_untrimmed_aln_analysis}"
  --output "maxalign=${file_og_maxalign}"
  --parameter "mode=${mode_gene_evolution}"
  --parameter "retain_query=${retain_query_in_maxalign}"
)
if [[ "${mode_gene_evolution}" == "query2family" ]]; then
  maxalign_provenance_args+=(--input "query=${file_query_gene}")
fi
gg_artifact_prepare_stage maxalign_needs_update run_maxalign "${maxalign_provenance_args[@]}" || exit $?
if [[ ${maxalign_needs_update} -eq 1 && ${run_maxalign} -eq 1 ]]; then
  gg_step_start "${task}"

  seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_untrimmed_aln_analysis}" --out-file "${og_id}.cds.aln.fasta"
  maxalign_keep_regex=""
  if [[ "${mode_gene_evolution}" == "query2family" && ${retain_query_in_maxalign} -eq 0 ]]; then
    echo "Query sequence(s) is NOT necessarily retained in MaxAlign."
  elif [[ "${mode_gene_evolution}" == "query2family" && ${retain_query_in_maxalign} -eq 1 ]]; then
    echo "Query sequence(s) is retained in MaxAlign."
    maxalign_keep_regex=$(
      python - "${file_query_gene}" << 'PY'
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

patterns = [f"(?i:^{re.escape(gene)}$)" for gene in normalized_ids]
print(','.join(patterns))
PY
    )
    if [[ -z "${maxalign_keep_regex}" ]]; then
      echo "Warning: No query IDs were parsed for MaxAlign --keep_seq_name_regex. Running without keep constraints."
    fi
  else
    maxalign_keep_regex=""
  fi

  maxalign_cmd=(
    cdskit maxalign
    --seq_file "${og_id}.cds.aln.fasta"
    --out_file "${og_id}.maxalign.output.fasta"
  )
  if [[ -n "${maxalign_keep_regex}" ]]; then
    maxalign_cmd+=(--keep_seq_name_regex "${maxalign_keep_regex}")
  fi
  "${maxalign_cmd[@]}"

  echo "Number of sequences before MaxAlign: $(gg_count_fasta_records "${og_id}.cds.aln.fasta")"
  echo "Number of sequences after MaxAlign: $(gg_count_fasta_records "${og_id}.maxalign.output.fasta")"

  seqkit seq --threads "${GG_TASK_CPUS}" "${og_id}.maxalign.output.fasta" --out-file "${og_id}.maxalign.out.fa.gz"
  mv_out "${og_id}.maxalign.out.fa.gz" "${file_og_maxalign}"
  rm -f -- "${og_id}.maxalign.output.fasta"
  gg_artifact_record "${maxalign_provenance_args[@]}"
else
  gg_step_skip "${task}"
fi
if [[ ${run_maxalign} -eq 1 ]]; then
  switch_alignment_analysis_source "${file_og_maxalign}"
fi

task="TrimAl"
disable_if_no_input_file "run_trimal" "${file_og_untrimmed_aln_analysis}"
trimal_needs_update=0
trimal_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.trimal.json"
  --step "trimal"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "untrimmed_alignment=${file_og_untrimmed_aln_analysis}"
  --output "trimal=${file_og_trimal}"
  --parameter "input_sequence_mode=${input_sequence_mode}"
  --parameter "genetic_code=${genetic_code}"
)
gg_artifact_prepare_stage trimal_needs_update run_trimal "${trimal_provenance_args[@]}" || exit $?
if [[ ${trimal_needs_update} -eq 1 && ${run_trimal} -eq 1 ]]; then
  gg_step_start "${task}"

  if [[ "${input_sequence_mode}" == "protein" ]]; then
    seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_untrimmed_aln_analysis}" --out-file untrimmed.pep.fasta
    trimal \
      -in untrimmed.pep.fasta \
      -out "${og_id}.cds.trimal.tmp2.fasta" \
      -automated1
  else
    seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" --threads "${GG_TASK_CPUS}" "${file_og_untrimmed_aln_analysis}" |
      sed -e '/^1 1$/d' -e 's/_frame=1[[:space:]]*//' \
        > untrimmed.pep.fasta

    seqkit seq --remove-gaps --threads "${GG_TASK_CPUS}" \
      "${file_og_untrimmed_aln_analysis}" \
      > untrimmed.cds.degap.fasta

    trimal \
      -in untrimmed.pep.fasta \
      -backtrans untrimmed.cds.degap.fasta \
      -out "${og_id}.cds.trimal.tmp1.fasta" \
      -ignorestopcodon \
      -automated1

    cdskit rmseq \
      --seq_file "${og_id}.cds.trimal.tmp1.fasta" \
      --problematic_percent 100 |
      cdskit hammer \
        --seq_file "-" \
        --codon_table "${genetic_code}" \
        --nail 4 \
        --out_file "${og_id}.cds.trimal.tmp2.fasta"
  fi

  if [[ -s "${og_id}.cds.trimal.tmp2.fasta" ]]; then
    echo "Copying. Output file detected for the task: ${task}"
    seqkit seq --threads "${GG_TASK_CPUS}" "${og_id}.cds.trimal.tmp2.fasta" --out-file "${og_id}.cds.trimal.out.fa.gz"
    mv_out "${og_id}.cds.trimal.out.fa.gz" "${file_og_trimal}"
  fi
  gg_artifact_record "${trimal_provenance_args[@]}"
else
  gg_step_skip "${task}"
fi
if [[ ${run_trimal} -eq 1 ]]; then
  set_analysis_file trimmed_aln "${file_og_trimal}"
fi

task="ClipKIT"
disable_if_no_input_file "run_clipkit" "${file_og_untrimmed_aln_analysis}"
clipkit_needs_update=0
clipkit_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.clipkit.json"
  --step "clipkit"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "untrimmed_alignment=${file_og_untrimmed_aln_analysis}"
  --output "clipkit=${file_og_clipkit}"
  --optional-output "clipkit_log=${file_og_clipkit_log}"
  --parameter "input_sequence_mode=${input_sequence_mode}"
  --parameter "genetic_code=${genetic_code}"
  --parameter "mode=smart-gap"
)
gg_artifact_prepare_stage clipkit_needs_update run_clipkit "${clipkit_provenance_args[@]}" || exit $?
if [[ ${clipkit_needs_update} -eq 1 && ${run_clipkit} -eq 1 ]]; then
  gg_step_start "${task}"
  rm -f -- "${file_og_clipkit}" "${file_og_clipkit_log}"
  rm -f -- "${og_id}.cds.clipkit.tmp.fasta.log" "${og_id}.cds.clipkit.hammer.fasta.log"

  seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_untrimmed_aln_analysis}" --out-file "${og_id}.cds.clipkit.input.fasta"

  if [[ "${input_sequence_mode}" == "protein" ]]; then
    clipkit \
      "${og_id}.cds.clipkit.input.fasta" \
      --mode smart-gap \
      --sequence_type aa \
      --input_file_format "fasta" \
      --output_file_format "fasta" \
      --output "${og_id}.cds.clipkit.hammer.fasta" \
      --log
  else
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
      --codon_table "${genetic_code}" \
      --nail 4 \
      --seq_file "${og_id}.cds.clipkit.tmp.fasta" |
      cdskit rmseq \
        --problematic_percent 100 \
        --out_file "${og_id}.cds.clipkit.hammer.fasta"
  fi

  if [[ -s "${og_id}.cds.clipkit.hammer.fasta" ]]; then
    echo "Copying. Output file detected for the task: ${task}"
    seqkit seq --threads "${GG_TASK_CPUS}" "${og_id}.cds.clipkit.hammer.fasta" --out-file "${og_id}.cds.clipkit.out.fa.gz"
    mv_out "${og_id}.cds.clipkit.out.fa.gz" "${file_og_clipkit}"
    if [[ -s "${og_id}.cds.clipkit.tmp.fasta.log" ]]; then
      cp_out "${og_id}.cds.clipkit.tmp.fasta.log" "${file_og_clipkit_log}"
    elif [[ -s "${og_id}.cds.clipkit.hammer.fasta.log" ]]; then
      cp_out "${og_id}.cds.clipkit.hammer.fasta.log" "${file_og_clipkit_log}"
    fi
  fi
  rm -f -- "${og_id}.cds.clipkit.input.fasta"
  gg_artifact_record "${clipkit_provenance_args[@]}"
else
  gg_step_skip "${task}"
fi
if [[ ${run_clipkit} -eq 1 ]]; then
  set_analysis_file trimmed_aln "${file_og_clipkit}"
fi

task="AMAS for cleaned alignment"
disable_if_no_input_file "run_amas_cleaned" "${file_og_trimmed_aln_analysis}"
amas_cleaned_needs_update=0
amas_cleaned_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.amas_cleaned.json"
  --step "amas_cleaned"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "alignment=${file_og_trimmed_aln_analysis}"
  --output "amas=${file_og_amas_cleaned}"
  --parameter "data_type=${amas_data_type}"
)
gg_artifact_prepare_stage amas_cleaned_needs_update run_amas_cleaned "${amas_cleaned_provenance_args[@]}" || exit $?
if [[ ${amas_cleaned_needs_update} -eq 1 && ${run_amas_cleaned} -eq 1 ]]; then
  gg_step_start "${task}"
  seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file "${og_id}.amas.cleaned.input.fasta"

  AMAS.py summary \
    --in-format fasta \
    --data-type "${amas_data_type}" \
    --in-files "${og_id}.amas.cleaned.input.fasta"

  mv_out summary.txt "${file_og_amas_cleaned}"
  rm -f -- "${og_id}.amas.cleaned.input.fasta"
  gg_artifact_record "${amas_cleaned_provenance_args[@]}"
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
# shellcheck shell=bash
# Sourced by gg_gene_evolution_core.sh.

task="IQ-TREE"
disable_if_no_input_file "run_iqtree" "${file_og_trimmed_aln_analysis}"
iqtree_needs_update=0
iqtree_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.iqtree.json"
  --step "iqtree"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "trimmed_alignment=${file_og_trimmed_aln_analysis}"
  --output "iqtree=${file_og_iqtree_tree}"
  --parameter "model=${generax_model}"
  --parameter "input_sequence_mode=${input_sequence_mode}"
  --parameter "genetic_code=${genetic_code}"
  --parameter "generax_enabled=${run_generax}"
  --parameter "fast_mode_threshold=${iqtree_fast_mode_gt}"
)
gg_artifact_prepare_stage iqtree_needs_update run_iqtree "${iqtree_provenance_args[@]}" || exit $?
if [[ ${iqtree_needs_update} -eq 1 && ${run_iqtree} -eq 1 ]]; then
  gg_step_start "${task}"
  num_seq=$(gg_count_fasta_records "${file_og_trimmed_aln_analysis}")
  if [[ ${num_seq} -ge 4 ]]; then
    if [[ ${run_generax} -eq 1 ]]; then
      other_iqtree_params=()
      file_tree="${og_id}.treefile"
      use_ufboot=0
      echo "run_generax=1: disabling UFBOOT in the initial IQ-TREE run. Support will be computed on the GeneRax topology."
    else
      other_iqtree_params=(--ufboot 1000 --bnni)
      file_tree="${og_id}.contree"
      use_ufboot=1
    fi
  else
    other_iqtree_params=()
    file_tree="${og_id}.treefile"
    use_ufboot=0
  fi
  if [[ ${num_seq} -gt ${iqtree_fast_mode_gt} ]]; then
    if [[ ${use_ufboot} -eq 1 ]]; then
      echo "Disabling IQ-TREE UFBOOT because fast mode is enabled for large alignments (${num_seq} > ${iqtree_fast_mode_gt})."
      other_iqtree_params=()
      file_tree="${og_id}.treefile"
      use_ufboot=0
    fi
    other_iqtree_params+=(--fast)
  fi

  if [[ "${input_sequence_mode}" == "protein" ]]; then
    assert_gene_evolution_aa_model_for_protein_mode "${task}"
    seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file iqtree_input.fa
  elif [[ ${run_generax} -eq 1 ]]; then
    if gene_evolution_model_is_aa "${generax_model}"; then
      echo "Specified substitution model was interpreted as an amino acid model (base model = ${generax_model%%+*})."
      seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" "${file_og_trimmed_aln_analysis}" |
        sed -e '/^1 1$/d' -e 's/_frame=1[[:space:]]*//' \
          > iqtree_input.fa
    else
      echo "Specified substitution model was interpreted as a nucleotide model (base model = ${generax_model%%+*})."
      seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file iqtree_input.fa
    fi
  else
    seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file iqtree_input.fa
  fi
  echo "IQ-TREE starting..."

  iqtree_model_string=${generax_model}
  build_iqtree_mem_args

  iqtree \
    -s iqtree_input.fa \
    -m "${iqtree_model_string}" \
    -T AUTO \
    --threads-max "${GG_TASK_CPUS}" \
    --prefix "${og_id}" \
    "${IQTREE_MEM_ARGS[@]}" \
    --seed 12345 \
    --redo \
    "${other_iqtree_params[@]}"

  cp_out "${file_tree}" "${file_og_iqtree_tree}"
  gg_artifact_record "${iqtree_provenance_args[@]}"
else
  gg_step_skip "${task}"
fi

task="Gene tree rooting"
disable_if_no_input_file "run_tree_root" "${file_og_unrooted_tree_analysis}"
tree_root_needs_update=0
tree_root_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.tree_root.json"
  --step "tree_root"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "unrooted_tree=${file_og_unrooted_tree_analysis}"
  --output "rooted_tree=${file_og_rooted_tree}"
  --output "rooting_log=${file_og_rooted_log}"
  --parameter "method=${tree_rooting_method}"
  --parameter "species_parser=${species_label_parser}"
  --parameter "species_regex=${species_label_regex}"
  --parameter "species_map_present=$([[ -n "${species_label_map_tsv}" ]] && echo 1 || echo 0)"
)
if [[ "${tree_rooting_method}" == "notung" ]]; then
  tree_root_provenance_args+=(--input "species_tree=${species_tree_pruned}")
fi
if [[ -n "${species_label_map_tsv}" ]]; then
  tree_root_provenance_args+=(--input "species_map=${species_label_map_tsv}")
fi
gg_artifact_prepare_stage tree_root_needs_update run_tree_root "${tree_root_provenance_args[@]}" || exit $?
if [[ ${tree_root_needs_update} -eq 1 && ${run_tree_root} -eq 1 ]]; then
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
    nwkit_root_args=(root --method "${nwkit_root_method}" --infile "${file_og_unrooted_tree_analysis}")
    if [[ "${nwkit_root_method}" == "taxonomy" ]]; then
      nwkit_root_args+=(--species-parser "${species_label_parser}")
      if [[ -n "${species_label_regex}" ]]; then
        nwkit_root_args+=(--species-regex "${species_label_regex}")
      fi
      if [[ -n "${species_label_map_tsv}" ]]; then
        nwkit_root_args+=(--species-map-tsv "${species_label_map_tsv}")
      fi
    fi
    nwkit "${nwkit_root_args[@]}" |
      nwkit label --target intnode --force yes --outfile "${og_id}.root.tmp.nwk"
    mv_out "${og_id}.root.tmp.nwk" "${file_og_rooted_tree}"
    {
      echo "tree_rooting_method=${tree_rooting_method}"
      echo "nwkit_method=${nwkit_root_method}"
    } > "${og_id}.root.txt"
    mv_out "${og_id}.root.txt" "${file_og_rooted_log}"
  fi
  gg_artifact_record "${tree_root_provenance_args[@]}"
else
  gg_step_skip "${task}"
fi

task="Orthogroup extraction with NWKIT"
run_orthogroup_extraction_original="${run_orthogroup_extraction}" # This variable may be disabled by disable_if_no_input_file but the original value is necessary to properly update file_og_*_analysis
disable_if_no_input_file "run_orthogroup_extraction" "${file_query_gene:-}" "${file_og_trimmed_aln_analysis}" "${file_og_rooted_tree_analysis}"
orthogroup_extraction_needs_update=0
orthogroup_extraction_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.orthogroup_extraction.json"
  --step "orthogroup_extraction"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "query=${file_query_gene:-}"
  --input "trimmed_alignment=${file_og_trimmed_aln_analysis}"
  --input "rooted_tree=${file_og_rooted_tree_analysis}"
  --output "extracted_unrooted_tree=${file_og_orthogroup_extraction_nwk}"
  --output "extracted_rooted_tree=${file_og_orthogroup_extraction_rooted_nwk}"
  --output "extracted_fasta=${file_og_orthogroup_extraction_fasta}"
  --parameter "input_sequence_mode=${input_sequence_mode}"
  --parameter "genetic_code=${genetic_code}"
  --parameter "rooting_method=${tree_rooting_method}"
)
if [[ -s "${file_og_query_blast}" ]]; then
  orthogroup_extraction_provenance_args+=(--input "query_blast=${file_og_query_blast}")
fi
gg_artifact_prepare_stage orthogroup_extraction_needs_update run_orthogroup_extraction "${orthogroup_extraction_provenance_args[@]}" || exit $?
if [[ ${run_orthogroup_extraction} -eq 1 ]]; then
  run_orthogroup_extraction_original=1
fi
if [[ ${orthogroup_extraction_needs_update} -eq 1 && ${run_orthogroup_extraction} -eq 1 ]]; then
  gg_step_start "${task}"

  if [[ "$(head -c 1 "${file_query_gene}")" == ">" ]]; then
    echo "Fasta format was detected. Query IDs absent from the tree will be replaced by their best tree-backed query BLAST hit."
  else
    echo "Gene IDs were detected."
    cp_out "${file_query_gene}" "${dir_output_active}/query_gene/$(basename "${file_query_gene}")"
  fi

  run_nwkit_subtree() {
    local infile=$1
    local result_file=$2
    echo "Running nwkit subtree for ${infile}"
    local info_txt
    info_txt=$(nwkit subtree --infile "${infile}" --leaves "${comma_separated_genes}" --orthogroup "yes" --dup-conf-score-threshold 0 2> /dev/null | nwkit info 2> /dev/null)
    local num_leaf
    num_leaf=$(awk -F': *' '/Number of leaves/ {print $2; exit}' <<< "${info_txt}")
    printf '%s\t%s\n' "${num_leaf}" "${infile}" > "${result_file}"
  }

  subtree_infiles=()
  if [[ "${tree_rooting_method}" == "notung" && -d "./${og_id}.notung.root" ]]; then
    mapfile -t subtree_infiles < <(
      find "./${og_id}.notung.root" -maxdepth 1 -type f |
        awk -v og="${og_id}" '$0 ~ (og "\\.iqtree\\.nwk\\.rooting\\.[0-9]+$") {print}' |
        sort -V
    )
  fi
  if [[ ${#subtree_infiles[@]} -eq 0 ]]; then
    if [[ -s "${file_og_rooted_tree_analysis}" ]]; then
      subtree_infiles=("${file_og_rooted_tree_analysis}")
    else
      echo "No rooted tree is available for orthogroup extraction."
      exit 1
    fi
  fi
  comma_separated_genes=$(
    python - "${file_query_gene}" "${file_og_query_blast:-}" "${subtree_infiles[0]}" << 'PY'
import csv
import math
import os
import sys

query_gene_path = sys.argv[1]
query_blast_path = sys.argv[2]
tree_path = sys.argv[3]


def normalize_id(value):
    return value.strip().replace("−", "-")


def parse_query_ids(path):
    query_ids = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        first_char = handle.read(1)
        handle.seek(0)
        if first_char == ">":
            for line in handle:
                if not line.startswith(">"):
                    continue
                query_id = normalize_id(line[1:].split()[0])
                if query_id:
                    query_ids.append(query_id)
        else:
            for line in handle:
                query_id = normalize_id(line)
                if query_id:
                    query_ids.append(query_id)
    return query_ids


def parse_newick_leaves(path):
    text = open(path, "r", encoding="utf-8", errors="replace").read()
    leaves = set()
    index = 0
    while index < len(text):
        if text[index] not in "(,":
            index += 1
            continue
        index += 1
        while index < len(text) and text[index].isspace():
            index += 1
        if index >= len(text) or text[index] == ")":
            continue
        if text[index] == "'":
            index += 1
            token = []
            while index < len(text):
                char = text[index]
                if char == "'":
                    if index + 1 < len(text) and text[index + 1] == "'":
                        token.append("'")
                        index += 2
                        continue
                    index += 1
                    break
                token.append(char)
                index += 1
            label = "".join(token)
        else:
            start = index
            while index < len(text) and text[index] not in ":,();":
                index += 1
            label = text[start:index]
        label = normalize_id(label)
        if label:
            leaves.add(label)
    return leaves


def to_float(value, default):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def best_tree_backed_hits(path, leaves):
    best = {}
    if not path or not os.path.exists(path):
        return best
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row_index, row in enumerate(reader):
            qacc = normalize_id(row.get("qacc", ""))
            sacc = normalize_id(row.get("sacc", ""))
            if not qacc or not sacc or sacc not in leaves:
                continue
            evalue = to_float(row.get("min_evalue") or row.get("evalue"), math.inf)
            bitscore = to_float(row.get("max_bitscore") or row.get("bitscore"), -math.inf)
            rank = (evalue, -bitscore, row_index)
            current = best.get(qacc)
            if current is None or rank < current["rank"]:
                best[qacc] = {
                    "sacc": sacc,
                    "evalue": evalue,
                    "bitscore": bitscore,
                    "rank": rank,
                }
    return best


query_ids = parse_query_ids(query_gene_path)
tree_leaves = parse_newick_leaves(tree_path)
best_hits = best_tree_backed_hits(query_blast_path, tree_leaves)
resolved = []
seen = set()

for query_id in query_ids:
    if query_id in tree_leaves:
        seed = query_id
    elif query_id in best_hits:
        hit = best_hits[query_id]
        seed = hit["sacc"]
        print(
            "Orthogroup extraction seed fallback: "
            f"query {query_id} was not found in the tree; using best tree-backed "
            f"query BLAST hit {seed} (evalue={hit['evalue']}, bitscore={hit['bitscore']}).",
            file=sys.stderr,
        )
    else:
        print(
            "Warning: Orthogroup extraction seed skipped because neither the query "
            f"nor a tree-backed query BLAST hit was found: {query_id}",
            file=sys.stderr,
        )
        continue
    if seed in seen:
        print(
            f"Orthogroup extraction seed duplicate skipped: {seed} "
            f"(from query {query_id}).",
            file=sys.stderr,
        )
        continue
    seen.add(seed)
    resolved.append(seed)

if not resolved:
    print(
        "No orthogroup extraction seed genes were found in the rooted tree. Exiting.",
        file=sys.stderr,
    )
    sys.exit(1)

print(",".join(resolved))
PY
  )
  echo "Seed genes for orthogroup extraction: ${comma_separated_genes}"

  printf 'num_leaf\tfile\n' > tmp_num_leaf.tsv
  rm -rf -- tmp_num_leaf.parts
  mkdir -p tmp_num_leaf.parts
  subtree_index=0
  for subtree_infile in "${subtree_infiles[@]}"; do
    wait_until_jobn_le "${GG_TASK_CPUS}"
    run_nwkit_subtree "${subtree_infile}" "tmp_num_leaf.parts/${subtree_index}.tsv" &
    gg_background_register "$!"
    subtree_index=$((subtree_index + 1))
  done
  wait_for_background_jobs
  for ((subtree_index = 0; subtree_index < ${#subtree_infiles[@]}; subtree_index++)); do
    cat "tmp_num_leaf.parts/${subtree_index}.tsv" >> tmp_num_leaf.tsv
  done
  rm -rf -- tmp_num_leaf.parts

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

  nwkit subtree --infile "${min_leaf_file}" --leaves "${comma_separated_genes}" --orthogroup "yes" --dup-conf-score-threshold 0 \
    --outfile "${og_id}.orthogroup_seed.tmp.nwk"

  seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file tmp.trimmed.input.fasta
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
  mv_out "${og_id}.orthogroup_seed.tmp.nwk" "${file_og_orthogroup_extraction_rooted_nwk}"

  if [[ "${input_sequence_mode}" == "protein" ]]; then
    seqkit seq --threads "${GG_TASK_CPUS}" tmp.fasta --out-file "${og_id}.orthogroup_extraction.tmp.fasta"
  else
    cdskit hammer --nail 4 -s tmp.fasta -o "${og_id}.orthogroup_extraction.tmp.fasta"
  fi
  seqkit seq --threads "${GG_TASK_CPUS}" "${og_id}.orthogroup_extraction.tmp.fasta" --out-file "${og_id}.orthogroup_extraction.out.fa.gz"
  mv_out "${og_id}.orthogroup_extraction.out.fa.gz" "${file_og_orthogroup_extraction_fasta}"
  rm -f -- "${og_id}.orthogroup_extraction.tmp.fasta"
  gg_artifact_record "${orthogroup_extraction_provenance_args[@]}"
else
  gg_step_skip "${task}"
fi
if [[ ${run_orthogroup_extraction_original} -eq 1 && -s "${file_og_orthogroup_extraction_nwk}" && -s "${file_og_orthogroup_extraction_rooted_nwk}" && -s "${file_og_orthogroup_extraction_fasta}" ]]; then
  set_analysis_file unrooted_tree "${file_og_orthogroup_extraction_nwk}"
  set_analysis_file rooted_tree "${file_og_orthogroup_extraction_rooted_nwk}"
  set_analysis_file trimmed_aln "${file_og_orthogroup_extraction_fasta}"
fi

task="GeneRax"
disable_if_no_input_file "run_generax" "${file_og_trimmed_aln_analysis}" "${file_og_unrooted_tree_analysis}" "${species_tree_pruned}"
generax_needs_update=0
generax_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.generax.json"
  --step "generax"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "trimmed_alignment=${file_og_trimmed_aln_analysis}"
  --input "unrooted_tree=${file_og_unrooted_tree_analysis}"
  --input "species_tree=${species_tree_pruned}"
  --output "generax_nwk=${file_og_generax_nwk}"
  --output "generax_xml=${file_og_generax_xml}"
  --output "generax_nhx=${file_og_generax_nhx}"
  --parameter "model=${generax_model}"
  --parameter "reconciliation_model=${generax_rec_model}"
  --parameter "input_sequence_mode=${input_sequence_mode}"
  --parameter "genetic_code=${genetic_code}"
)
gg_artifact_prepare_stage generax_needs_update run_generax "${generax_provenance_args[@]}" || exit $?
if [[ ${generax_needs_update} -eq 1 && ${run_generax} -eq 1 ]]; then
  gg_step_start "${task}"

  if [[ "${input_sequence_mode}" == "protein" ]]; then
    assert_gene_evolution_aa_model_for_protein_mode "${task}"
    seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file generax_input_alignment.fas
  elif gene_evolution_model_is_aa "${generax_model}"; then
    echo "Specified substitution model was interpreted as an amino acid model (base model = ${generax_model%%+*})."
    seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" "${file_og_trimmed_aln_analysis}" |
      sed -e '/^1 1$/d' -e 's/_frame=1[[:space:]]*//' \
        > generax_input_alignment.fas
  else
    echo "Specified substitution model was interpreted as a nucleotide model (base model = ${generax_model%%+*})."
    seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file generax_input_alignment.fas
  fi

  nwkit drop --target intnode --support yes --name yes \
    --infile "${file_og_unrooted_tree_analysis}" \
    --outfile generax_input_gene_tree.nwk

  #avoid multifurcating tree
  R -q -e "library(ape); t=read.tree(\"generax_input_gene_tree.nwk\"); t=multi2di(t,random=FALSE); write.tree(t, \"generax_input_gene_tree_bi.nwk\")"

  generate_generax_mapfile() {
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

  # GeneRax runs within one scheduler task; ignore host Grid Engine discovery and use local launch inside containers.
  mpiexec_args=(mpiexec -oversubscribe -np "${GG_TASK_CPUS}")
  mpi_env_args=(env OMPI_MCA_ras=^gridengine OMPI_MCA_plm=isolated OMPI_MCA_plm_rsh_agent=/bin/false OMPI_MCA_btl=^openib)
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
    --seed 12345 < /dev/null

  echo "GeneRax exit code = $?"

  generax_out_sptree="./generax_${og_id}/species_trees/starting_species_tree.newick" # generax v2.0
  if [[ -s "${generax_out_sptree}" ]]; then
    lock_file="${species_tree_generax}.lock"
    (
      if ! gg_shared_lock_acquire "${lock_file}" "GeneRax species tree copy"; then
        exit 1
      fi
      gg_shared_lock_start_heartbeat "${lock_file}"
      heartbeat_pid=${GG_SHARED_LOCK_HEARTBEAT_PID:-}
      cleanup_generax_tree_lock() {
        gg_shared_lock_stop_heartbeat "${heartbeat_pid}"
        gg_shared_lock_release "${lock_file}"
      }
      trap cleanup_generax_tree_lock EXIT
      if [[ ! -s "${species_tree_generax}" ]]; then
        echo "copying GeneRax output species tree (first writer only)."
        cp_out "${generax_out_sptree}" "${species_tree_generax}"
      fi
    ) || exit 1
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
  gg_artifact_record "${generax_provenance_args[@]}"
else
  gg_step_skip "${task}"
fi

task="IQ-TREE UFBOOT on GeneRax topology"
generax_ufboot_needs_update=0
generax_ufboot_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.generax_ufboot.json"
  --step "generax_ufboot"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "generax_tree=${file_og_generax_nwk}"
  --input "trimmed_alignment=${file_og_trimmed_aln_analysis}"
  --output "ufboot_tree=${file_og_iqtree_generax_ufboot}"
  --parameter "model=${generax_model}"
  --parameter "input_sequence_mode=${input_sequence_mode}"
  --parameter "genetic_code=${genetic_code}"
  --parameter "ufboot=1000"
  --parameter "bnni=yes"
  --parameter "seed=12345"
)
if [[ ${run_generax} -eq 1 ]]; then
  if [[ ! -s "${file_og_generax_nwk}" ]]; then
    echo "Skipped: ${task}. Missing GeneRax output tree: ${file_og_generax_nwk}"
  elif [[ ! -s "${file_og_trimmed_aln_analysis}" ]]; then
    echo "Skipped: ${task}. Missing alignment: ${file_og_trimmed_aln_analysis}"
  else
    gg_artifact_prepare_stage generax_ufboot_needs_update run_generax "${generax_ufboot_provenance_args[@]}" || exit $?
  fi
  if [[ ${generax_ufboot_needs_update} -eq 1 && ${run_generax} -eq 1 ]]; then
    gg_step_start "${task}"
    num_seq=$(gg_count_fasta_records "${file_og_trimmed_aln_analysis}")
    if [[ ${num_seq} -lt 4 ]]; then
      echo "UFBOOT requires >=4 sequences. Using the GeneRax topology without bootstrap support."
      nwkit drop --target root --length yes --infile "${file_og_generax_nwk}" --outfile "${og_id}.generax_ufboot.tmp.nwk"
      mv_out "${og_id}.generax_ufboot.tmp.nwk" "${file_og_iqtree_generax_ufboot}"
    else
      if [[ "${input_sequence_mode}" == "protein" ]]; then
        assert_gene_evolution_aa_model_for_protein_mode "${task}"
        seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file "${og_id}.generax_ufboot.input.fa"
      elif gene_evolution_model_is_aa "${generax_model}"; then
        echo "Specified substitution model was interpreted as an amino acid model (base model = ${generax_model%%+*})."
        seqkit translate --allow-unknown-codon --transl-table "${genetic_code}" "${file_og_trimmed_aln_analysis}" |
          sed -e '/^1 1$/d' -e 's/_frame=1[[:space:]]*//' \
            > "${og_id}.generax_ufboot.input.fa"
      else
        echo "Specified substitution model was interpreted as a nucleotide model (base model = ${generax_model%%+*})."
        seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file "${og_id}.generax_ufboot.input.fa"
      fi

      nwkit drop --target intnode --support yes --name yes --infile "${file_og_generax_nwk}" |
        nwkit drop --target all --length yes --outformat 9 --outfile "${og_id}.generax_ufboot.constraint.nwk"

      other_iqtree_params=(--ufboot 1000 --bnni)
      if [[ ${num_seq} -gt ${iqtree_fast_mode_gt} ]]; then
        echo "Skipping IQ-TREE --fast in UFBOOT-on-GeneRax mode because the options are incompatible."
      fi

      build_iqtree_mem_args
      iqtree \
        -s "${og_id}.generax_ufboot.input.fa" \
        -g "${og_id}.generax_ufboot.constraint.nwk" \
        -m "${generax_model}" \
        -T AUTO \
        --threads-max "${GG_TASK_CPUS}" \
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
    gg_artifact_record "${generax_ufboot_provenance_args[@]}"
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
notung_reconcil_needs_update=0
notung_reconcil_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.notung_reconciliation.json"
  --step "notung_reconciliation"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "rooted_tree=${file_og_rooted_tree}"
  --input "species_tree=${species_tree_pruned}"
  --output "notung_reconciliation=${file_og_notung_reconcil}"
  --parameter "species_parser=${species_label_parser}"
)
gg_artifact_prepare_stage notung_reconcil_needs_update run_notung_reconcil "${notung_reconcil_provenance_args[@]}" || exit $?
if [[ ${notung_reconcil_needs_update} -eq 1 && ${run_notung_reconcil} -eq 1 ]]; then
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
    rm -f -- "${og_id}.notung_reconcile.zip"
    zip -rq "${og_id}.notung_reconcile.zip" "${og_id}.notung_reconcile"
    python "${gg_support_dir}/atomic_zip_publish.py" \
      --source "${og_id}.notung_reconcile.zip" \
      --destination "${file_og_notung_reconcil}" \
      --expected-prefix "${og_id}.notung_reconcile"
  fi
  gg_artifact_record "${notung_reconcil_provenance_args[@]}"
else
  gg_step_skip "${task}"
fi

task="Species-tree-guided divergence time estimation"
disable_if_no_input_file "run_tree_dating" "${species_tree_pruned}" "${file_og_unrooted_tree_analysis}"
tree_dating_needs_update=0
tree_dating_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.tree_dating.json"
  --step "tree_dating"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "species_tree=${species_tree_pruned}"
  --input "unrooted_tree=${file_og_unrooted_tree_analysis}"
  --output "dated_tree=${file_og_dated_tree}"
  --output "dating_log=${file_og_dated_tree_log}"
  --parameter "generax_enabled=${run_generax}"
  --parameter "max_age=${radte_max_age}"
  --parameter "species_parser=${species_label_parser}"
  --parameter "species_regex=${species_label_regex}"
)
if [[ ${run_generax} -eq 1 ]]; then
  tree_dating_provenance_args+=(
    --input "generax_species_tree=${species_tree_generax}"
    --input "generax_nhx=${file_og_generax_nhx}"
  )
else
  tree_dating_provenance_args+=(--input "notung_reconciliation=${file_og_notung_reconcil}")
fi
if [[ -n "${species_label_map_tsv}" ]]; then
  tree_dating_provenance_args+=(--input "species_map=${species_label_map_tsv}")
fi
gg_artifact_prepare_stage tree_dating_needs_update run_tree_dating "${tree_dating_provenance_args[@]}" || exit $?
if [[ ${tree_dating_needs_update} -eq 1 && ${run_tree_dating} -eq 1 ]]; then
  gg_step_start "${task}"
  radte_args=()

  if [[ ${run_generax} -eq 1 ]]; then
    radte_args+=("--species_tree=${species_tree_generax}")
    radte_args+=("--generax_nhx=${file_og_generax_nhx}")
  else
    gg_extract_expected_zip_prefix \
      "${file_og_notung_reconcil}" \
      "${og_id}.notung_reconcile"
    if [[ -s ./${og_id}.notung_reconcile/${og_id}.root.nwk.reconciled.0 ]]; then
      cp_out ./"${og_id}".notung_reconcile/"${og_id}".root.nwk.reconciled.0 ./"${og_id}".notung_reconcile/"${og_id}".root.nwk.reconciled
      cp_out ./"${og_id}".notung_reconcile/"${og_id}".root.nwk.reconciled.0.parsable.txt ./"${og_id}".notung_reconcile/"${og_id}".root.nwk.reconciled.parsable.txt
    fi
    radte_args+=("--species_tree=${species_tree_pruned}")
    radte_args+=("--gene_tree=./${og_id}.notung_reconcile/${og_id}.root.nwk.reconciled")
    radte_args+=("--notung_parsable=./${og_id}.notung_reconcile/${og_id}.root.nwk.reconciled.parsable.txt")
  fi
  radte_args+=("--species-parser=${species_label_parser}")
  if [[ -n "${species_label_regex}" ]]; then
    radte_args+=("--species-regex=${species_label_regex}")
  fi
  if [[ -n "${species_label_map_tsv}" ]]; then
    radte_args+=("--species-map-tsv=${species_label_map_tsv}")
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
  gg_artifact_record "${tree_dating_provenance_args[@]}"
else
  gg_step_skip "${task}"
fi
# shellcheck shell=bash
# Sourced by gg_gene_evolution_core.sh.

task="Expression matrix preparation"
disable_if_no_input_file "run_generate_expression_matrix" "${file_og_trimmed_aln_analysis}"
expression_needs_update=0
expression_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.expression_matrix.json"
  --step "expression_matrix"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "trimmed_alignment=${file_og_trimmed_aln_analysis}"
  --input "species_expression=${dir_sp_expression}"
  --output "expression_matrix=${file_og_expression}"
  --parameter "expression_value_type=${exp_value_type}"
)
gg_artifact_prepare_stage expression_needs_update run_generate_expression_matrix "${expression_provenance_args[@]}" || exit $?
if [[ ${expression_needs_update} -eq 1 && ${run_generate_expression_matrix} -eq 1 ]]; then
  gg_step_start "${task}"
  seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file "${og_id}.trait_matrix_input.fasta"

  python "${gg_support_dir}/get_trait_matrix.py" \
    --dir_trait "${dir_sp_expression}" \
    --seqfile "${og_id}.trait_matrix_input.fasta" \
    --ncpu "${GG_TASK_CPUS}" \
    --outfile expression_matrix.tsv
  rm -f -- "${og_id}.trait_matrix_input.fasta"
  if [[ -s expression_matrix.tsv ]]; then
    mv_out expression_matrix.tsv "${file_og_expression}"
  fi
  gg_artifact_record "${expression_provenance_args[@]}"
else
  gg_step_skip "${task}"
fi

task="Promoter fasta generation"
disable_if_no_input_file "run_extract_promoter_fasta" "${file_og_gff_info}"
promoter_fasta_needs_update=0
promoter_fasta_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.promoter_fasta.json"
  --step "promoter_fasta"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "gff_info=${file_og_gff_info}"
  --input "species_genome=${dir_sp_genome}"
  --output "promoter_fasta=${file_og_promoter_fasta}"
  --parameter "promoter_bp=${promoter_bp}"
)
gg_artifact_prepare_stage promoter_fasta_needs_update run_extract_promoter_fasta "${promoter_fasta_provenance_args[@]}" || exit $?
if [[ ${promoter_fasta_needs_update} -eq 1 && ${run_extract_promoter_fasta} -eq 1 ]]; then
  gg_step_start "${task}"

  python "${gg_support_dir}/get_promoter_fasta.py" \
    --dir_genome "${dir_sp_genome}" \
    --geneinfo_tsv "${file_og_gff_info}" \
    --seqkit_exe "seqkit" \
    --outfile "${og_id}.promoter.tmp.fasta" \
    --promoter_bp "${promoter_bp}" \
    --ncpu "${GG_TASK_CPUS}"
  if [[ -s "${og_id}.promoter.tmp.fasta" ]]; then
    seqkit seq --threads "${GG_TASK_CPUS}" "${og_id}.promoter.tmp.fasta" --out-file "${og_id}.promoter.out.fa.gz"
    mv_out "${og_id}.promoter.out.fa.gz" "${file_og_promoter_fasta}"
    rm -f -- "${og_id}.promoter.tmp.fasta"
  fi
  gg_artifact_record "${promoter_fasta_provenance_args[@]}"
else
  gg_step_skip "${task}"
fi

task="fimo"
jaspar_path=""
fimo_manifest="${dir_output_active}/artifact_provenance/${og_id}.fimo.json"
if [[ ${run_fimo} -eq 1 || -s "${file_og_fimo}" || -s "${fimo_manifest}" ]]; then
  if ! jaspar_path=$(ensure_jaspar_file "${gg_workspace_dir}" "${jaspar_file}") || [[ -z "${jaspar_path}" ]]; then
    echo "Failed to prepare JASPAR motif file (${jaspar_file}). Exiting."
    exit 1
  fi
fi
disable_if_no_input_file "run_fimo" "${file_og_promoter_fasta}" "${jaspar_path}"
fimo_needs_update=0
fimo_provenance_args=(
  --manifest "${fimo_manifest}"
  --step "fimo"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "promoter_fasta=${file_og_promoter_fasta}"
  --input "jaspar=${jaspar_path}"
  --output "fimo=${file_og_fimo}"
  --parameter "jaspar_selector=${jaspar_file}"
)
gg_artifact_prepare_stage fimo_needs_update run_fimo "${fimo_provenance_args[@]}" || exit $?
if [[ ${fimo_needs_update} -eq 1 && ${run_fimo} -eq 1 ]]; then
  gg_step_start "${task}"
  seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_promoter_fasta}" --out-file "${og_id}.fimo.input.fasta"

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
  gg_artifact_record "${fimo_provenance_args[@]}"
else
  gg_step_skip "${task}"
fi

task="Tree pruning"
disable_if_no_input_file "run_tree_pruning" "${file_og_expression}" "${file_og_untrimmed_aln_analysis}" "${file_og_trimmed_aln_analysis}" "${file_og_unrooted_tree_analysis}" "${file_og_rooted_tree_analysis}"
tree_pruning_needs_update=0
tree_pruning_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.tree_pruning.json"
  --step "tree_pruning"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "expression=${file_og_expression}"
  --input "untrimmed_alignment=${file_og_untrimmed_aln_analysis}"
  --input "trimmed_alignment=${file_og_trimmed_aln_analysis}"
  --input "unrooted_tree=${file_og_unrooted_tree_analysis}"
  --input "rooted_tree=${file_og_rooted_tree_analysis}"
  --output "untrimmed_alignment_pruned=${file_og_untrimmed_aln_pruned}"
  --output "trimmed_alignment_pruned=${file_og_trimmed_aln_pruned}"
  --output "unrooted_tree_pruned=${file_og_unrooted_tree_pruned}"
  --output "rooted_tree_pruned=${file_og_rooted_tree_pruned}"
  --optional-output "dated_tree_pruned=${file_og_dated_tree_pruned}"
  --parameter "dated_tree_present=$([[ -s "${file_og_dated_tree_analysis}" ]] && echo 1 || echo 0)"
)
if [[ -s "${file_og_dated_tree_analysis}" ]]; then
  tree_pruning_provenance_args+=(--input "dated_tree=${file_og_dated_tree_analysis}")
fi
gg_artifact_prepare_stage tree_pruning_needs_update run_tree_pruning "${tree_pruning_provenance_args[@]}" || exit $?
if [[ ${tree_pruning_needs_update} -eq 1 && ${run_tree_pruning} -eq 1 ]]; then
  gg_step_start "${task}"
  rm -f -- \
    "${file_og_untrimmed_aln_pruned}" \
    "${file_og_trimmed_aln_pruned}" \
    "${file_og_unrooted_tree_pruned}" \
    "${file_og_rooted_tree_pruned}" \
    "${file_og_dated_tree_pruned}"

  cut -f 1 "${file_og_expression}" | tail -n +2 > target_genes.txt

  if [[ -s "${file_og_untrimmed_aln_analysis}" ]]; then
    seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_untrimmed_aln_analysis}" --out-file "${og_id}.untrimmed.input.fasta"
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
    seqkit seq --threads "${GG_TASK_CPUS}" "${og_id}.untrimmed.pruned.tmp.fasta" --out-file "${og_id}.untrimmed.pruned.out.fa.gz"
    mv_out "${og_id}.untrimmed.pruned.out.fa.gz" "${file_og_untrimmed_aln_pruned}"
    rm -f -- "${og_id}.untrimmed.pruned.tmp.fasta"
  fi

  if [[ -s "${file_og_trimmed_aln_analysis}" ]]; then
    seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file "${og_id}.trimmed.input.fasta"
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
    seqkit seq --threads "${GG_TASK_CPUS}" "${og_id}.trimmed.pruned.tmp.fasta" --out-file "${og_id}.trimmed.pruned.out.fa.gz"
    mv_out "${og_id}.trimmed.pruned.out.fa.gz" "${file_og_trimmed_aln_pruned}"
    rm -f -- "${og_id}.trimmed.pruned.tmp.fasta"
  fi

  mapfile -t prune_genes < <(sed -e '/^[[:space:]]*$/d' target_genes.txt)
  prune_pattern=""
  if [[ ${#prune_genes[@]} -gt 0 ]]; then
    prune_pattern=$(
      printf '%s\n' "${prune_genes[@]}" |
        sed -e 's/[][(){}.^$+*?|\\-]/\\&/g' |
        paste -sd'|' -
    )
  fi

  if [[ -s "${file_og_unrooted_tree_analysis}" ]]; then
    if [[ ${#prune_genes[@]} -eq 0 ]]; then
      cat "${file_og_unrooted_tree_analysis}"
    else
      nwkit prune \
        --infile "${file_og_unrooted_tree_analysis}" \
        --pattern "^(${prune_pattern})$" \
        --invert-match yes
    fi |
      nwkit drop --target root --length yes |
      nwkit sanitize --remove-singleton yes --resolve-polytomy no \
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
        --invert-match yes \
        --preserve-properties yes
    fi |
      nwkit drop --target root --length yes --preserve-properties yes |
      nwkit sanitize --remove-singleton yes --resolve-polytomy no --preserve-properties yes \
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
        --invert-match yes
    fi |
      nwkit drop --target root --length yes |
      nwkit sanitize --remove-singleton yes --resolve-polytomy no \
        > "${og_id}.dated.pruned.tmp.nwk"
    mv_out "${og_id}.dated.pruned.tmp.nwk" "${file_og_dated_tree_pruned}"
  fi
  gg_artifact_record "${tree_pruning_provenance_args[@]}"
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
if [[ -s "${file_og_expression}" && ${run_l1ou} -eq 1 ]]; then
  # This block should be run after tree pruning.
  num_gene_trait=$(($(wc -l < "${file_og_expression}") - 1)) # -1 for header
  num_gene_tree=$(gg_count_fasta_records "${file_og_trimmed_aln_analysis}")
  if [[ ${num_gene_trait} -eq ${num_gene_tree} ]]; then
    echo "num_gene_trait (${num_gene_trait}) and num_gene_tree (${num_gene_tree}) matched."
  else
    echo "num_gene_trait (${num_gene_trait}) and num_gene_tree (${num_gene_tree}) did not match."
    if [[ ${run_tree_pruning} -ne 1 && ${run_l1ou} -eq 1 ]]; then
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
    "${file_og_l1ou_fit_rdata}"
    "${file_og_l1ou_fit_conv_rdata}"
    "${file_og_l1ou_fit_tree}"
    "${file_og_l1ou_fit_regime}"
    "${file_og_l1ou_fit_leaf}"
    "${file_og_l1ou_fit_plot}"
    "${file_og_iqtree_anc}"
    "${file_og_csubst_b}"
    "${file_og_csubst_cb_stats}"
    "${file_og_csubst_scan}"
    "${file_og_csubst_scan_units}"
    "${file_og_csubst_scan_foreground_branch}"
    "${file_og_csubst_scan_plot}"
    "${file_og_csubst_scan_log}"
    "${file_og_gene_pgls}"
    "${file_og_gene_pgls_plot}"
    "${file_og_species_nwkit_pgls}"
    "${file_og_species_rphylopars_pgls}"
    "${file_og_pgls_comparison}"
    "${file_og_pgls_method_status}"
    "${file_og_pgls_method_audit}"
    "${file_og_species_expression_summary}"
    "${file_og_species_expression_audit}"
    "${file_og_species_response_tip_summary}"
    "${file_og_species_response_sampling_covariance}"
    "${file_og_species_predictor_tip_summary}"
    "${file_og_species_predictor_sampling_covariance}"
    "${file_og_rsc_status}"
    "${file_og_rsc_pgls}"
    "${file_og_rsc_reconciliation}"
    "${file_og_rsc_gene_contrasts}"
    "${file_og_rsc_species_contrasts}"
    "${file_og_rsc_response_sampling_covariance}"
    "${file_og_rsc_response_tip_summary}"
    "${file_og_rsc_predictor_sampling_covariance}"
    "${file_og_rsc_predictor_tip_summary}"
    "${file_og_rsc_random_effects}"
    "${file_og_rsc_sensitivity}"
    "${file_og_rsc_trait_origins}"
    "${file_og_rsc_audit}"
    "${file_og_rsc_log}"
    "${file_og_expression_trait_pgls_log}"
  )
  for ((i = 2; i <= csubst_max_arity; i++)); do
    varname="file_og_csubst_cb_${i}"
    if [[ -n "${!varname:-}" ]]; then
      file_to_remove+=("${!varname}")
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
      if [[ $((${num_gene_pruned} * 2)) -lt ${num_stat_branch_row} ]]; then
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
# shellcheck shell=bash
# Sourced by gg_gene_evolution_core.sh.

task="Parameter estimation for mapdNdS"
disable_if_no_input_file "run_mapdnds_parameter_estimation" "${file_og_rooted_tree_analysis}" "${file_og_trimmed_aln_analysis}"
mapdnds_parameter_needs_update=0
mapdnds_parameter_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.mapdnds_parameter.json"
  --step "mapdnds_parameter"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "rooted_tree=${file_og_rooted_tree_analysis}"
  --input "trimmed_alignment=${file_og_trimmed_aln_analysis}"
  --output "mapdnds_parameter=${file_og_mapdnds_parameter}"
  --parameter "model=GY+F3X4+G4"
  --parameter "genetic_code=${genetic_code}"
)
gg_artifact_prepare_stage mapdnds_parameter_needs_update run_mapdnds_parameter_estimation "${mapdnds_parameter_provenance_args[@]}" || exit $?
if [[ ${mapdnds_parameter_needs_update} -eq 1 && ${run_mapdnds_parameter_estimation} -eq 1 ]]; then
  gg_step_start "${task}"

  nwkit drop --target intnode --support yes --name yes \
    --infile "${file_og_rooted_tree_analysis}" |
    nwkit sanitize --remove-singleton yes --resolve-polytomy no \
      > mapdnds_input.nwk

  seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file ./mapdnds_input.fasta

  # F3X4+G4 shouldn not be changed otherwise iqtree2mapnh.py has to be updated.
  build_iqtree_mem_args
  iqtree \
    -s mapdnds_input.fasta \
    -m "GY+F3X4+G4" \
    -te mapdnds_input.nwk \
    -T AUTO \
    --threads-max "${GG_TASK_CPUS}" \
    --seqtype "CODON${genetic_code}" \
    --prefix "${og_id}.iqtree2mapdNdS" \
    "${IQTREE_MEM_ARGS[@]}" \
    --ancestral \
    --seed 12345 \
    --redo

  python "${gg_support_dir}/iqtree2mapnh.py" \
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
    rm -f -- "${og_id}.mapdnds_parameter.zip"
    zip -r "${og_id}.mapdnds_parameter.zip" "${og_id}.mapdnds_parameter"
    python "${gg_support_dir}/atomic_zip_publish.py" \
      --source "${og_id}.mapdnds_parameter.zip" \
      --destination "${file_og_mapdnds_parameter}" \
      --expected-prefix "${og_id}.mapdnds_parameter" \
      --remove-source
  else
    echo "iqtree2mapnh.params was not generated."
  fi
  gg_artifact_record "${mapdnds_parameter_provenance_args[@]}"
else
  gg_step_skip "${task}"
fi

task="mapdNdS main run"
disable_if_no_input_file "run_mapdnds" "${file_og_mapdnds_parameter}" "${file_og_trimmed_aln_analysis}"
mapdnds_needs_update=0
mapdnds_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.mapdnds.json"
  --step "mapdnds"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "mapdnds_parameter=${file_og_mapdnds_parameter}"
  --input "trimmed_alignment=${file_og_trimmed_aln_analysis}"
  --output "dN_tree=${file_og_mapdnds_dn}"
  --output "dS_tree=${file_og_mapdnds_ds}"
  --parameter "genetic_code=${genetic_code}"
)
gg_artifact_prepare_stage mapdnds_needs_update run_mapdnds "${mapdnds_provenance_args[@]}" || exit $?
if [[ ${mapdnds_needs_update} -eq 1 && ${run_mapdnds} -eq 1 ]]; then
  gg_step_start "${task}"

  gg_extract_expected_zip_prefix \
    "${file_og_mapdnds_parameter}" \
    "${og_id}.mapdnds_parameter"
  cd "${dir_tmp}/${og_id}.mapdnds_parameter" || exit 1
  seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file ./mapdnds_input.fasta
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
  cd "${dir_tmp}" || exit 1
  gg_artifact_record "${mapdnds_provenance_args[@]}"
else
  gg_step_skip "${task}"
fi

task="CodeML two-ratio model"
disable_if_no_input_file "run_codeml_two_ratio" "${file_og_rooted_tree_analysis}" "${file_og_trimmed_aln_analysis}" "${file_sp_trait}"
codeml_two_ratio_needs_update=0
codeml_two_ratio_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.codeml_two_ratio.json"
  --step "codeml_two_ratio"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "rooted_tree=${file_og_rooted_tree_analysis}"
  --input "trimmed_alignment=${file_og_trimmed_aln_analysis}"
  --input "foreground_trait=${file_sp_trait}"
  --output "codeml_two_ratio=${file_og_codeml_two_ratio}"
  --parameter "genetic_code=${genetic_code}"
  --parameter "model=two_ratio"
)
gg_artifact_prepare_stage codeml_two_ratio_needs_update run_codeml_two_ratio "${codeml_two_ratio_provenance_args[@]}" || exit $?
if [[ ${codeml_two_ratio_needs_update} -eq 1 && ${run_codeml_two_ratio} -eq 1 ]]; then
  gg_step_start "${task}"

  binarize_species_trait "${file_sp_trait}" species_trait_binary.tsv
  write_species_trait_foreground_regex_table species_trait_binary.tsv foreground.tsv
  IFS=$'\t' read -r -a colname_array < foreground.tsv

  if [[ ${#colname_array[@]} -le 1 ]]; then
    echo "No trait columns were detected in ${file_sp_trait}. Skipping ${task}."
  else
    for ((i = 1; i < ${#colname_array[@]}; i++)); do
      trait="${colname_array[$i]}"
      echo "Processing trait: ${trait}"
      awk -F'\t' -v trait_col="$((i + 1))" 'NR>1 && $trait_col == 1 { print $1 }' foreground.tsv > "foreground_${trait}.txt"

      target_spnode=$(paste -sd'|' "foreground_${trait}.txt") # Bar-separated list of target species nodes
      echo "Regular expression for CodeML foreground node search: ${target_spnode}"

      nwkit drop \
        --infile "${file_og_rooted_tree_analysis}" \
        --target "intnode" \
        --support "yes" \
        --name "yes" |
        nwkit mark \
          --pattern "${target_spnode}" \
          --insert-text "#1" \
          --insert-separator "" \
          --target "mrca" \
          --target-only-clade "yes" \
          --outfile "codeml_input_${trait}.nwk"
      echo "CodeML input tree: $(< "codeml_input_${trait}.nwk")"

      bash "${gg_support_dir}/shorten_fasta_newick_names.sh" \
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
          "${gg_support_dir}/codeml/codeml_two_ratio.ctl.template" \
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
    for ((i = 1; i < ${#colname_array[@]}; i++)); do
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
  gg_artifact_record "${codeml_two_ratio_provenance_args[@]}"
else
  gg_step_skip "${task}"
fi

task="HyPhy dN-dS estimation"
disable_if_no_input_file "run_hyphy_dnds" "${file_og_rooted_tree_analysis}" "${file_og_trimmed_aln_analysis}"
hyphy_dnds_needs_update=0
hyphy_dnds_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.hyphy_dnds.json"
  --step "hyphy_dnds"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "rooted_tree=${file_og_rooted_tree_analysis}"
  --input "trimmed_alignment=${file_og_trimmed_aln_analysis}"
  --output "hyphy_dnds=${file_og_hyphy_dnds}"
  --parameter "genetic_code=${genetic_code}"
  --parameter "frequencies=CF3x4"
  --parameter "type=local"
)
gg_artifact_prepare_stage hyphy_dnds_needs_update run_hyphy_dnds "${hyphy_dnds_provenance_args[@]}" || exit $?
if [[ ${hyphy_dnds_needs_update} -eq 1 && ${run_hyphy_dnds} -eq 1 ]]; then
  gg_step_start "${task}"

  nwkit drop --target intnode --support yes --name yes \
    --infile "${file_og_rooted_tree_analysis}" |
    nwkit sanitize --remove-singleton yes --resolve-polytomy no \
      > "hyphy_input.nwk"

  seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file "hyphy_input.fasta"

  hyphy_genetic_code=$(get_hyphy_genetic_code "${genetic_code}")

  hyphy "${gg_support_dir}/hyphy/FitMG94.bf" \
    --alignment "hyphy_input.fasta" \
    --tree "hyphy_input.nwk" \
    --code "${hyphy_genetic_code}" \
    --frequencies "CF3x4" \
    --type "local" \
    --lrt "No" \
    --rooted "Yes" \
    --CPU "${GG_TASK_CPUS}"

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
  gg_artifact_record "${hyphy_dnds_provenance_args[@]}"
else
  gg_step_skip "${task}"
fi

# HyPhy RELAX


task="HyPhy RELAX"
disable_if_no_input_file "run_hyphy_relax" "${file_og_rooted_tree_analysis}" "${file_og_trimmed_aln_analysis}" "${file_sp_trait}"
hyphy_relax_needs_update=0
hyphy_relax_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.hyphy_relax.json"
  --step "hyphy_relax"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "rooted_tree=${file_og_rooted_tree_analysis}"
  --input "trimmed_alignment=${file_og_trimmed_aln_analysis}"
  --input "foreground_trait=${file_sp_trait}"
  --output "hyphy_relax=${file_og_hyphy_relax}"
  --parameter "genetic_code=${genetic_code}"
  --parameter "foreground=1"
  --parameter "mode=classic_minimal"
)
gg_artifact_prepare_stage hyphy_relax_needs_update run_hyphy_relax "${hyphy_relax_provenance_args[@]}" || exit $?
if [[ ${hyphy_relax_needs_update} -eq 1 && ${run_hyphy_relax} -eq 1 ]]; then
  gg_step_start "${task}"
  run_hyphy_relax_for_all_traits 1 "${file_og_hyphy_relax}"
  if [[ -s "${file_og_hyphy_relax}" ]]; then
    echo "The task has completed successfully: ${task}"
  else
    echo "The task has failed: ${task}"
  fi
  gg_artifact_record "${hyphy_relax_provenance_args[@]}"
else
  gg_step_skip "${task}"
fi

task="HyPhy RELAX with reversed foreground/background"
disable_if_no_input_file "run_hyphy_relax_reversed" "${file_og_rooted_tree_analysis}" "${file_og_trimmed_aln_analysis}" "${file_sp_trait}"
hyphy_relax_reversed_needs_update=0
hyphy_relax_reversed_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.hyphy_relax_reversed.json"
  --step "hyphy_relax_reversed"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "rooted_tree=${file_og_rooted_tree_analysis}"
  --input "trimmed_alignment=${file_og_trimmed_aln_analysis}"
  --input "foreground_trait=${file_sp_trait}"
  --output "hyphy_relax_reversed=${file_og_hyphy_relax_reversed}"
  --parameter "genetic_code=${genetic_code}"
  --parameter "foreground=0"
  --parameter "mode=classic_minimal"
)
gg_artifact_prepare_stage hyphy_relax_reversed_needs_update run_hyphy_relax_reversed "${hyphy_relax_reversed_provenance_args[@]}" || exit $?
if [[ ${hyphy_relax_reversed_needs_update} -eq 1 && ${run_hyphy_relax_reversed} -eq 1 ]]; then
  gg_step_start "${task}"
  run_hyphy_relax_for_all_traits 0 "${file_og_hyphy_relax_reversed}"
  if [[ -s "${file_og_hyphy_relax_reversed}" ]]; then
    echo "The task has completed successfully: ${task}"
  else
    echo "The task has failed: ${task}"
  fi
  gg_artifact_record "${hyphy_relax_reversed_provenance_args[@]}"
else
  gg_step_skip "${task}"
fi

task="Stochastic character mapping of intron evolution"
disable_if_no_input_file "run_scm_intron" "${file_og_gff_info}" "${file_og_dated_tree_analysis}"
scm_intron_needs_update=0
scm_intron_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.scm_intron.json"
  --step "scm_intron"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "gff_info=${file_og_gff_info}"
  --input "dated_tree=${file_og_dated_tree_analysis}"
  --output "scm_summary=${file_og_scm_intron_summary}"
  --optional-output "scm_plot=${file_og_scm_intron_plot}"
  --parameter "intron_gain_rate=${intron_gain_rate}"
  --parameter "retrotransposition_rate=${retrotransposition_rate}"
  --parameter "nrep=1000"
)
gg_artifact_prepare_stage scm_intron_needs_update run_scm_intron "${scm_intron_provenance_args[@]}" || exit $?
if [[ ${scm_intron_needs_update} -eq 1 && ${run_scm_intron} -eq 1 ]]; then
  gg_step_start "${task}"
  rm -f -- "${file_og_scm_intron_summary}" "${file_og_scm_intron_plot}"
  rm -f -- intron_evolution_summary.tsv intron_evolution_plot.pdf

  Rscript "${gg_support_dir}/scm_intron_evolution.r" \
    --tree_file="${file_og_dated_tree_analysis}" \
    --trait_file="${file_og_gff_info}" \
    --intron_gain_rate="${intron_gain_rate}" \
    --retrotransposition_rate="${retrotransposition_rate}" \
    --nrep=1000 \
    --nslots="${GG_TASK_CPUS}"

  cp_out intron_evolution_summary.tsv "${file_og_scm_intron_summary}"
  if [[ -e intron_evolution_plot.pdf ]]; then
    cp_out intron_evolution_plot.pdf "${file_og_scm_intron_plot}"
  fi
  gg_artifact_record "${scm_intron_provenance_args[@]}"
else
  gg_step_skip "${task}"
fi

task="l1ou"
disable_if_no_input_file "run_l1ou" "${file_og_trimmed_aln_analysis}" "${file_og_expression}" "${file_og_dated_tree_analysis}"
l1ou_needs_update=0
l1ou_provenance_args=(
  --manifest "${dir_output_active}/artifact_provenance/${og_id}.l1ou.json"
  --step "l1ou"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "trimmed_alignment=${file_og_trimmed_aln_analysis}"
  --input "expression=${file_og_expression}"
  --input "dated_tree=${file_og_dated_tree_analysis}"
  --output "fit=${file_og_l1ou_fit_rdata}"
  --output "tree=${file_og_l1ou_fit_tree}"
  --output "regime=${file_og_l1ou_fit_regime}"
  --output "leaf=${file_og_l1ou_fit_leaf}"
  --output "plot=${file_og_l1ou_fit_plot}"
  --parameter "criterion=${l1ou_criterion}"
  --parameter "alpha_upper=${l1ou_alpha_upper}"
  --parameter "convergence=${l1ou_convergence}"
  --parameter "nbootstrap=${l1ou_nbootstrap}"
  --parameter "large_tree_num_gene=${large_tree_num_gene}"
  --parameter "large_tree_max_nshift=${large_tree_max_nshift}"
)
if [[ ${l1ou_convergence} -eq 1 ]]; then
  l1ou_provenance_args+=(--output "convergent_fit=${file_og_l1ou_fit_conv_rdata}")
fi
if [[ ${l1ou_nbootstrap} -gt 0 ]]; then
  l1ou_provenance_args+=(--output "bootstrap=${file_og_l1ou_bootstrap}")
fi
gg_artifact_prepare_stage l1ou_needs_update run_l1ou "${l1ou_provenance_args[@]}" || exit $?
if [[ ${l1ou_needs_update} -eq 1 && ${run_l1ou} -eq 1 ]]; then
  gg_step_start "${task}"

  num_gene=$(gg_count_fasta_records "${file_og_trimmed_aln_analysis}")
  if [[ ${num_gene} -ge ${large_tree_num_gene} ]]; then
    max_nshift=${large_tree_max_nshift}
  else
    max_nshift=0
  fi

  fit_ind_file=''
  if [[ ${l1ou_use_fit_file} -eq 1 && -s "${file_og_l1ou_fit_rdata}" ]]; then
    fit_ind_file=${file_og_l1ou_fit_rdata}
  fi

  l1ou_cmd=(
    Rscript "${gg_support_dir}/detect_OU_shift_kfl1ou.r"
    --max_nshift="${max_nshift}"
    --tree_file="${file_og_dated_tree_analysis}"
    --trait_file="${file_og_expression}"
    --nslots="${GG_TASK_CPUS}"
    --criterion="${l1ou_criterion}"
    --nbootstrap="${l1ou_nbootstrap}"
    --fit_ind_file="${fit_ind_file}"
    --fit_conv_file=''
    --alpha_upper="${l1ou_alpha_upper}"
    --detect_convergence="${l1ou_convergence}"
    --replicate_sep="_"
  )
  "${l1ou_cmd[@]}"

  mv_out fit_ind.RData "${file_og_l1ou_fit_rdata}"
  mv_out l1ou_tree.tsv "${file_og_l1ou_fit_tree}"
  mv_out l1ou_regime.tsv "${file_og_l1ou_fit_regime}"
  mv_out l1ou_leaf.tsv "${file_og_l1ou_fit_leaf}"
  mv_out l1ou_plot.pdf "${file_og_l1ou_fit_plot}"
  if [[ ${l1ou_nbootstrap} -gt 0 ]]; then
    cp_out l1ou_bootstrap.tsv "${file_og_l1ou_bootstrap}"
  fi
  if [[ ${l1ou_convergence} -eq 1 ]]; then
    cp_out fit_conv.RData "${file_og_l1ou_fit_conv_rdata}"
  fi
  gg_artifact_record "${l1ou_provenance_args[@]}"

else
  gg_step_skip "${task}"
fi
# shellcheck shell=bash
# Sourced by gg_gene_evolution_core.sh.

task="Expression-trait phylogenetic regression"
rsc_gene_tree="${file_og_dated_tree_analysis}"
if [[ "${rsc_gene_branch_length}" == "unit" && ! -s "${rsc_gene_tree}" ]]; then
  rsc_gene_tree="${file_og_rooted_tree_analysis}"
fi
rsc_reconciliation_tree="${file_og_rooted_tree_analysis}"
if [[ ! -s "${rsc_reconciliation_tree}" ]]; then
  rsc_reconciliation_tree="${rsc_gene_tree}"
fi
rsc_effective_event_source="${rsc_event_source}"
if [[ "${rsc_effective_event_source}" == "auto" ]]; then
  rsc_effective_event_source="lca"
  if [[ ${run_expression_trait_pgls} -eq 1 || -s "${file_og_rsc_status}" || -s "${file_og_expression_trait_pgls_provenance}" ]] &&
    [[ -s "${rsc_reconciliation_tree}" ]] &&
    python "${gg_support_dir}/reconciled_speciation_contrast.py" has-nhx \
      --tree "${rsc_reconciliation_tree}"
  then
    rsc_effective_event_source="nhx"
  fi
fi

if [[ ${run_expression_trait_pgls} -eq 1 ]]; then
  if ! command -v nwkit >/dev/null 2>&1; then
    echo "run_expression_trait_pgls=1 requires nwkit with the pgls and reconcile commands." >&2
    exit 1
  fi
  if [[ ${pgls_run_species_rphylopars} -eq 1 ]]; then
    if ! command -v Rscript >/dev/null 2>&1; then
      echo "pgls_methods includes species-rphylopars, but Rscript is unavailable." >&2
      exit 1
    fi
  fi
  if [[ -s "${file_og_expression}" ]]; then
    if [[ ! -s "${file_sp_trait}" ]]; then
      echo "run_expression_trait_pgls=1 requires the species-trait table: ${file_sp_trait}" >&2
      exit 1
    fi
    if [[ ! -s "${species_tree_pruned}" ]]; then
      echo "run_expression_trait_pgls=1 requires the pruned species tree: ${species_tree_pruned}" >&2
      exit 1
    fi
    if [[ ${pgls_run_rsc} -eq 1 && ! -s "${rsc_gene_tree}" ]]; then
      echo "The RSC method requires a rooted gene tree: ${rsc_gene_tree}" >&2
      if [[ "${rsc_gene_branch_length}" == "original" ]]; then
        echo "rsc_gene_branch_length=original requires run_tree_dating=1 or an existing dated gene tree." >&2
      fi
      exit 1
    fi
    if [[ ! -s "${rsc_reconciliation_tree}" ]]; then
      echo "Expression-trait PGLS requires a rooted reconciliation tree to map genes to species: ${rsc_reconciliation_tree}" >&2
      exit 1
    fi
    if [[ "${rsc_species_branch_length}" == "original" && "${species_tree_basename}" != "dated_species_tree" ]]; then
      echo "rsc_species_branch_length=original requires output/species_tree/species_tree_summary/dated_species_tree.nwk." >&2
      echo "Use rsc_species_branch_length=unit only when a dated species tree is unavailable." >&2
      exit 1
    fi
    if [[ -n "${rsc_expression_sample_metadata}" && ! -s "${rsc_expression_sample_metadata}" ]]; then
      echo "RSC expression-sample metadata was not found: ${rsc_expression_sample_metadata}" >&2
      exit 1
    fi
    if [[ -n "${species_paralog_sampling_covariance}" && ! -s "${species_paralog_sampling_covariance}" ]]; then
      echo "Species-PGLS paralog sampling covariance was not found: ${species_paralog_sampling_covariance}" >&2
      exit 1
    fi
    if [[ "${rsc_effective_event_source}" == "nhx" ]] &&
      ! python "${gg_support_dir}/reconciled_speciation_contrast.py" has-nhx \
        --tree "${rsc_reconciliation_tree}"
    then
      echo "rsc_event_source=nhx requires GeneRax-style NHX D annotations: ${rsc_reconciliation_tree}" >&2
      exit 1
    fi
  fi
fi

rsc_nwkit_identity=""
if [[ -r /opt/pg/logs/source_revisions.tsv ]]; then
  rsc_nwkit_identity=$(awk -F '\t' '$1 == "nwkit" { print $2; exit }' /opt/pg/logs/source_revisions.tsv)
fi
if [[ -z "${rsc_nwkit_identity}" ]] && command -v nwkit >/dev/null 2>&1; then
  rsc_nwkit_identity=$(nwkit --version 2>&1 | tail -n 1 || true)
fi
rsc_nwkit_identity="${rsc_nwkit_identity:-unavailable}"
rsc_rphylopars_identity="not-requested"
if [[ ${pgls_run_species_rphylopars} -eq 1 ]]; then
  rsc_rphylopars_identity=$(Rscript -e 'cat(as.character(utils::packageVersion("Rphylopars")))' 2>/dev/null || true)
  rsc_rphylopars_identity="${rsc_rphylopars_identity:-unavailable}"
fi

rsc_needs_update=0
rsc_provenance_args=(
  --manifest "${file_og_expression_trait_pgls_provenance}"
  --step "expression_trait_pgls"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "adapter=${gg_support_dir}/reconciled_speciation_contrast.py"
  --input "species_comparator=${gg_support_dir}/species_tree_pgls.py"
  --input "rphylopars_adapter=${gg_support_dir}/species_tree_rphylopars.R"
  --output "status=${file_og_rsc_status}"
  --output "pgls=${file_og_rsc_pgls}"
  --output "reconciliation=${file_og_rsc_reconciliation}"
  --output "gene_contrasts=${file_og_rsc_gene_contrasts}"
  --output "species_contrasts=${file_og_rsc_species_contrasts}"
  --output "response_sampling_covariance=${file_og_rsc_response_sampling_covariance}"
  --output "response_tip_summary=${file_og_rsc_response_tip_summary}"
  --output "predictor_sampling_covariance=${file_og_rsc_predictor_sampling_covariance}"
  --output "predictor_tip_summary=${file_og_rsc_predictor_tip_summary}"
  --output "random_effects=${file_og_rsc_random_effects}"
  --output "sensitivity=${file_og_rsc_sensitivity}"
  --output "trait_origins=${file_og_rsc_trait_origins}"
  --output "audit=${file_og_rsc_audit}"
  --output "log=${file_og_rsc_log}"
  --output "species_nwkit=${file_og_species_nwkit_pgls}"
  --output "species_rphylopars=${file_og_species_rphylopars_pgls}"
  --output "comparison=${file_og_pgls_comparison}"
  --output "method_status=${file_og_pgls_method_status}"
  --output "method_audit=${file_og_pgls_method_audit}"
  --output "species_expression_summary=${file_og_species_expression_summary}"
  --output "species_expression_audit=${file_og_species_expression_audit}"
  --output "species_response_tip_summary=${file_og_species_response_tip_summary}"
  --output "species_response_sampling_covariance=${file_og_species_response_sampling_covariance}"
  --output "species_predictor_tip_summary=${file_og_species_predictor_tip_summary}"
  --output "species_predictor_sampling_covariance=${file_og_species_predictor_sampling_covariance}"
  --output "combined_log=${file_og_expression_trait_pgls_log}"
  --parameter "methods=${pgls_methods}"
  --parameter "species_expression_aggregation=${species_expression_aggregation}"
  --parameter "species_paralog_missing=${species_paralog_missing}"
  --parameter "rphylopars_sampling_covariance=${rphylopars_sampling_covariance}"
  --parameter "responses=${rsc_responses}"
  --parameter "predictors=${rsc_predictors}"
  --parameter "predictor_mode=${rsc_predictor_mode}"
  --parameter "requested_event_source=${rsc_event_source}"
  --parameter "effective_event_source=${rsc_effective_event_source}"
  --parameter "speciation_coverage=${rsc_speciation_coverage}"
  --parameter "event_weighting=${rsc_event_weighting}"
  --parameter "model=${rsc_model}"
  --parameter "gene_branch_length=${rsc_gene_branch_length}"
  --parameter "gene_evolution_model=${rsc_gene_evolution_model}"
  --parameter "gene_evolution_parameter=${rsc_gene_evolution_parameter}"
  --parameter "species_branch_length=${rsc_species_branch_length}"
  --parameter "species_evolution_model=${rsc_species_evolution_model}"
  --parameter "species_evolution_parameter=${rsc_species_evolution_parameter}"
  --parameter "inference=${rsc_inference}"
  --parameter "bootstrap_replicates=${rsc_bootstrap_replicates}"
  --parameter "seed=${rsc_seed}"
  --parameter "confidence_level=${rsc_confidence_level}"
  --parameter "reml=${rsc_reml}"
  --parameter "min_species_events=${rsc_min_species_events}"
  --parameter "unmatched=${rsc_unmatched}"
  --parameter "replicate_separator=${rsc_replicate_separator}"
  --parameter "within_variance=${rsc_within_variance}"
  --parameter "technical_aggregation=${rsc_technical_aggregation}"
  --parameter "predictor_biological_id=${rsc_predictor_biological_id}"
  --parameter "predictor_technical_id=${rsc_predictor_technical_id}"
  --parameter "predictor_batch=${rsc_predictor_batch}"
  --parameter "predictor_within_variance=${rsc_predictor_within_variance}"
  --parameter "predictor_technical_aggregation=${rsc_predictor_technical_aggregation}"
  --parameter "predictor_standard_error_columns=${rsc_predictor_standard_error_columns}"
  --parameter "predictor_sample_size_columns=${rsc_predictor_sample_size_columns}"
  --parameter "categorical_predictors=${rsc_categorical_predictors}"
  --parameter "ordered_predictors=${rsc_ordered_predictors}"
  --parameter "factor_reference=${rsc_factor_reference}"
  --parameter "factor_coding=${rsc_factor_coding}"
  --parameter "categorical_replicate_policy=${rsc_categorical_replicate_policy}"
  --parameter "event_random_effect=${rsc_event_random_effect}"
  --parameter "lineage_random_slope=${rsc_lineage_random_slope}"
  --parameter "lineage_inference=${rsc_lineage_inference}"
  --parameter "lineage_leave_one_out=${rsc_lineage_leave_one_out}"
  --parameter "categorical_origin_diagnostics=${rsc_categorical_origin_diagnostics}"
  --parameter "origin_map_replicates=${rsc_origin_map_replicates}"
  --parameter "origin_map_threads=${rsc_origin_map_threads}"
  --parameter "origin_min_posterior=${rsc_origin_min_posterior}"
  --parameter "origin_leave_one_out=${rsc_origin_leave_one_out}"
  --parameter "allow_large_dense=${rsc_allow_large_dense}"
  --parameter "species_parser=${species_label_parser}"
  --parameter "species_regex=${species_label_regex}"
  --parameter "nwkit_identity=${rsc_nwkit_identity}"
  --parameter "rphylopars_identity=${rsc_rphylopars_identity}"
)
if [[ -s "${file_og_expression}" ]]; then
  rsc_provenance_args+=(
    --input "reconciliation_tree=${rsc_reconciliation_tree}"
    --input "species_tree=${species_tree_pruned}"
    --input "expression=${file_og_expression}"
    --input "species_traits=${file_sp_trait}"
  )
  if [[ ${pgls_run_rsc} -eq 1 ]]; then
    rsc_provenance_args+=(--input "gene_tree=${rsc_gene_tree}")
  fi
else
  rsc_provenance_args+=(--parameter "expression_input=unavailable")
fi
if [[ -s "${file_og_expression}" && -n "${rsc_expression_sample_metadata}" ]]; then
  rsc_provenance_args+=(--input "expression_sample_metadata=${rsc_expression_sample_metadata}")
fi
if [[ -s "${file_og_expression}" && -n "${species_paralog_sampling_covariance}" ]]; then
  rsc_provenance_args+=(--input "paralog_sampling_covariance=${species_paralog_sampling_covariance}")
fi
if [[ -s "${file_og_expression}" && -n "${species_label_map_tsv}" ]]; then
  rsc_provenance_args+=(--input "species_map=${species_label_map_tsv}")
fi
gg_artifact_prepare_stage rsc_needs_update run_expression_trait_pgls "${rsc_provenance_args[@]}" || exit $?
if [[ ${rsc_needs_update} -eq 1 && ${run_expression_trait_pgls} -eq 1 ]]; then
  gg_step_start "${task}"
  rsc_work_dir="${og_id}.rsc"
  rsc_prepared_expression="${rsc_work_dir}/expression.tsv"
  rsc_prepared_traits="${rsc_work_dir}/species_traits.tsv"
  rsc_analysis_plan="${rsc_work_dir}/analysis_plan.tsv"
  rsc_metadata="${rsc_work_dir}/metadata.tsv"
  rsc_preflight_reconciliation="${rsc_work_dir}/preflight.reconciliation.tsv"
  rsc_preflight_audit="${rsc_work_dir}/preflight.audit.jsonl"
  rsc_preflight_metadata="${rsc_work_dir}/preflight.metadata.tsv"
  rsc_bundle_list="${rsc_work_dir}/bundle_list.tsv"
  rsc_combined_prefix="${rsc_work_dir}/${og_id}_rsc"
  rsc_combined_status="${rsc_work_dir}/${og_id}_rsc.status.tsv"
  rsc_combined_audit="${rsc_work_dir}/${og_id}_rsc.audit.jsonl"
  rsc_log_tmp="${rsc_work_dir}/${og_id}_rsc.log"
  pgls_log_tmp="${rsc_work_dir}/${og_id}_expression_trait_pgls.log"
  species_nwkit_tmp="${rsc_work_dir}/${og_id}_species_nwkit.pgls.tsv"
  species_rphylopars_tmp="${rsc_work_dir}/${og_id}_species_rphylopars.pgls.tsv"
  pgls_comparison_tmp="${rsc_work_dir}/${og_id}_pgls.comparison.tsv"
  pgls_method_status_tmp="${rsc_work_dir}/${og_id}_pgls.method-status.tsv"
  pgls_method_audit_tmp="${rsc_work_dir}/${og_id}_pgls.method-audit.jsonl"
  species_expression_summary_tmp="${rsc_work_dir}/${og_id}_species_expression.tsv"
  species_expression_audit_tmp="${rsc_work_dir}/${og_id}_species_expression.audit.tsv"
  species_response_tip_summary_tmp="${rsc_work_dir}/${og_id}_species.response-tip-summary.tsv"
  species_response_sampling_covariance_tmp="${rsc_work_dir}/${og_id}_species.response-sampling-covariance.tsv"
  species_predictor_tip_summary_tmp="${rsc_work_dir}/${og_id}_species.predictor-tip-summary.tsv"
  species_predictor_sampling_covariance_tmp="${rsc_work_dir}/${og_id}_species.predictor-sampling-covariance.tsv"
  rm -rf -- "${rsc_work_dir}"
  mkdir -p "${rsc_work_dir}/bundles"
  : > "${rsc_log_tmp}"
  : > "${pgls_log_tmp}"

  rsc_preparation_status="not_estimable"
  rsc_preparation_reason="expression_input_unavailable"
  rsc_effective_responses=""
  rsc_response_biological_id=""
  rsc_response_technical_id=""
  rsc_response_batch=""
  rsc_response_standard_error_columns=""
  rsc_response_sample_size_columns=""
  if [[ -s "${file_og_expression}" ]]; then
    rsc_prepare_args=(
      prepare
      --expression "${file_og_expression}"
      --species-traits "${file_sp_trait}"
      --responses "${rsc_responses}"
      --predictors "${rsc_predictors}"
      --predictor-mode "${rsc_predictor_mode}"
      --replicate-separator "${rsc_replicate_separator}"
      --within-variance "${rsc_within_variance}"
      --predictor-biological-id "${rsc_predictor_biological_id}"
      --predictor-technical-id "${rsc_predictor_technical_id}"
      --predictor-batch "${rsc_predictor_batch}"
      --predictor-within-variance "${rsc_predictor_within_variance}"
      --predictor-standard-error-columns "${rsc_predictor_standard_error_columns}"
      --predictor-sample-size-columns "${rsc_predictor_sample_size_columns}"
      --categorical-predictors "${rsc_categorical_predictors}"
      --ordered-predictors "${rsc_ordered_predictors}"
      --factor-reference "${rsc_factor_reference}"
      --expression-output "${rsc_prepared_expression}"
      --species-traits-output "${rsc_prepared_traits}"
      --analysis-plan-output "${rsc_analysis_plan}"
      --metadata-output "${rsc_metadata}"
    )
    if [[ -n "${rsc_expression_sample_metadata}" ]]; then
      rsc_prepare_args+=(--sample-metadata "${rsc_expression_sample_metadata}")
    fi
    python "${gg_support_dir}/reconciled_speciation_contrast.py" "${rsc_prepare_args[@]}" \
      2>&1 | tee -a "${rsc_log_tmp}"

    rsc_preparation_status=$(rsc_read_metadata_value "${rsc_metadata}" status)
    rsc_preparation_reason=$(rsc_read_metadata_value "${rsc_metadata}" reason)
    rsc_effective_responses=$(rsc_read_metadata_value "${rsc_metadata}" responses)
    rsc_response_biological_id=$(rsc_read_metadata_value "${rsc_metadata}" response_biological_id)
    rsc_response_technical_id=$(rsc_read_metadata_value "${rsc_metadata}" response_technical_id)
    rsc_response_batch=$(rsc_read_metadata_value "${rsc_metadata}" response_batch)
    rsc_response_standard_error_columns=$(rsc_read_metadata_value "${rsc_metadata}" standard_error_columns)
    rsc_response_sample_size_columns=$(rsc_read_metadata_value "${rsc_metadata}" sample_size_columns)
    rsc_response_sampling_uncertainty=$(rsc_read_metadata_value "${rsc_metadata}" response_sampling_uncertainty)
    rsc_predictor_sampling_uncertainty=$(rsc_read_metadata_value "${rsc_metadata}" predictor_sampling_uncertainty)
    if [[ ${pgls_run_rsc} -eq 1 && "${rsc_model}" == "legacy" && ( "${rsc_response_sampling_uncertainty}" == "yes" || "${rsc_predictor_sampling_uncertainty}" == "yes" ) ]]; then
      echo "rsc_model=legacy cannot use response or predictor sampling uncertainty after input preparation." >&2
      exit 1
    fi
  else
    echo "Expression-trait PGLS not estimable for ${og_id}: ${rsc_preparation_reason}" | tee -a "${rsc_log_tmp}"
  fi

  rsc_not_estimable_reason="${rsc_preparation_reason}"

  if [[ "${rsc_preparation_status}" == "ready" ]]; then
    rsc_reconcile_args=(
      reconcile
      --infile "${rsc_reconciliation_tree}"
      --species-tree "${species_tree_pruned}"
      --tree-id "${og_id}"
      --event-source "${rsc_effective_event_source}"
      --species-parser "${species_label_parser}"
      --unmatched "${rsc_unmatched}"
      --outfile "${rsc_preflight_reconciliation}"
      --audit "${rsc_preflight_audit}"
    )
    if [[ -n "${species_label_regex}" ]]; then
      rsc_reconcile_args+=(--species-regex "${species_label_regex}")
    fi
    if [[ -n "${species_label_map_tsv}" ]]; then
      rsc_reconcile_args+=(--species-map-tsv "${species_label_map_tsv}")
    fi
    nwkit "${rsc_reconcile_args[@]}" 2>&1 | tee -a "${rsc_log_tmp}"
    if [[ ${pgls_run_rsc} -eq 1 ]]; then
      python "${gg_support_dir}/reconciled_speciation_contrast.py" inspect-reconciliation \
        --reconciliation "${rsc_preflight_reconciliation}" \
        --coverage "${rsc_speciation_coverage}" \
        --min-species-events "${rsc_min_species_events}" \
        --output "${rsc_preflight_metadata}" \
        2>&1 | tee -a "${rsc_log_tmp}"
      rsc_preflight_status=$(rsc_read_metadata_value "${rsc_preflight_metadata}" status)
      rsc_preflight_reason=$(rsc_read_metadata_value "${rsc_preflight_metadata}" reason)
      if [[ -n "${rsc_preflight_reason}" ]]; then
        if [[ -n "${rsc_not_estimable_reason}" ]]; then
          rsc_not_estimable_reason+=";${rsc_preflight_reason}"
        else
          rsc_not_estimable_reason="${rsc_preflight_reason}"
        fi
      fi
    else
      rsc_preflight_status="not_requested"
      rsc_not_estimable_reason="method_not_requested"
    fi
  else
    rsc_preflight_status="not_estimable"
  fi

  if [[ ${pgls_run_rsc} -eq 0 || "${rsc_preflight_status}" != "ready" ]]; then
    if [[ -z "${rsc_not_estimable_reason}" ]]; then
      rsc_not_estimable_reason="rsc_preflight_not_estimable"
    fi
    rsc_empty_args=(
      empty-bundle
      --output-prefix "${rsc_combined_prefix}"
      --status-output "${rsc_combined_status}"
      --audit-output "${rsc_combined_audit}"
      --tree-id "${og_id}"
      --reason "${rsc_not_estimable_reason}"
    )
    if [[ -s "${rsc_preflight_reconciliation}" ]]; then
      rsc_empty_args+=(--reconciliation "${rsc_preflight_reconciliation}")
    fi
    if [[ -s "${rsc_preflight_audit}" ]]; then
      rsc_empty_args+=(--source-audit "${rsc_preflight_audit}")
    fi
    python "${gg_support_dir}/reconciled_speciation_contrast.py" "${rsc_empty_args[@]}" \
      2>&1 | tee -a "${rsc_log_tmp}"
    echo "RSC not estimable for ${og_id}: ${rsc_not_estimable_reason}" | tee -a "${rsc_log_tmp}"
  else
    printf 'analysis_id\tprefix\taudit\n' > "${rsc_bundle_list}"
    while IFS=$'\t' read -r \
      rsc_analysis_id \
      rsc_analysis_predictors \
      rsc_analysis_categorical \
      rsc_analysis_ordered \
      rsc_analysis_reference \
      rsc_analysis_predictor_se \
      rsc_analysis_predictor_n \
      rsc_origin_applicable
    do
      [[ -n "${rsc_analysis_id}" ]] || continue
      rsc_analysis_prefix="${rsc_work_dir}/bundles/${rsc_analysis_id}"
      rsc_analysis_audit="${rsc_analysis_prefix}.audit.jsonl"
      rsc_pgls_args=(
        pgls
        --gene-tree "${rsc_gene_tree}"
        --reconciliation-tree "${rsc_reconciliation_tree}"
        --species-tree "${species_tree_pruned}"
        --expression "${rsc_prepared_expression}"
        --species-traits "${rsc_prepared_traits}"
        --responses "${rsc_effective_responses}"
        --predictors "${rsc_analysis_predictors}"
        --tree-id "${og_id}"
        --out-prefix "${rsc_analysis_prefix}"
        --audit "${rsc_analysis_audit}"
        --event-source "${rsc_effective_event_source}"
        --species-parser "${species_label_parser}"
        --unmatched "${rsc_unmatched}"
        --gene-branch-length "${rsc_gene_branch_length}"
        --gene-evolution-model "${rsc_gene_evolution_model}"
        --species-branch-length "${rsc_species_branch_length}"
        --species-evolution-model "${rsc_species_evolution_model}"
        --event-weighting "${rsc_event_weighting}"
        --speciation-coverage "${rsc_speciation_coverage}"
        --confidence-level "${rsc_confidence_level}"
        --model "${rsc_model}"
        --inference "${rsc_inference}"
        --bootstrap-replicates "${rsc_bootstrap_replicates}"
        --seed "${rsc_seed}"
        --reml "${rsc_reml}"
        --factor-coding "${rsc_factor_coding}"
        --categorical-replicate-policy "${rsc_categorical_replicate_policy}"
        --event-random-effect "${rsc_event_random_effect}"
        --lineage-random-slope "${rsc_lineage_random_slope}"
        --lineage-inference "${rsc_lineage_inference}"
        --lineage-leave-one-out "${rsc_lineage_leave_one_out}"
        --allow-large-dense "${rsc_allow_large_dense}"
      )
      if [[ ! "${rsc_gene_evolution_model}" =~ ^(brownian|independent)$ ]]; then
        rsc_pgls_args+=(--gene-evolution-parameter "${rsc_gene_evolution_parameter}")
      fi
      if [[ ! "${rsc_species_evolution_model}" =~ ^(brownian|independent)$ ]]; then
        rsc_pgls_args+=(--species-evolution-parameter "${rsc_species_evolution_parameter}")
      fi
      if [[ -n "${species_label_regex}" ]]; then
        rsc_pgls_args+=(--species-regex "${species_label_regex}")
      fi
      if [[ -n "${species_label_map_tsv}" ]]; then
        rsc_pgls_args+=(--species-map-tsv "${species_label_map_tsv}")
      fi
      if [[ -n "${rsc_response_biological_id}" ]]; then
        rsc_pgls_args+=(
          --biological-id "${rsc_response_biological_id}"
          --technical-aggregation "${rsc_technical_aggregation}"
          --within-variance "${rsc_within_variance}"
        )
      elif [[ "${rsc_within_variance}" == "known-se" ]]; then
        rsc_pgls_args+=(--within-variance known-se)
      fi
      if [[ -n "${rsc_response_technical_id}" ]]; then
        rsc_pgls_args+=(--technical-id "${rsc_response_technical_id}")
      fi
      if [[ -n "${rsc_response_batch}" ]]; then
        rsc_pgls_args+=(--batch "${rsc_response_batch}")
      fi
      if [[ -n "${rsc_response_standard_error_columns}" ]]; then
        rsc_pgls_args+=(--standard-error-columns "${rsc_response_standard_error_columns}")
      fi
      if [[ -n "${rsc_response_sample_size_columns}" ]]; then
        rsc_pgls_args+=(--sample-size-columns "${rsc_response_sample_size_columns}")
      fi
      if [[ -n "${rsc_predictor_biological_id}" ]]; then
        rsc_pgls_args+=(
          --predictor-biological-id "${rsc_predictor_biological_id}"
          --predictor-technical-aggregation "${rsc_predictor_technical_aggregation}"
          --predictor-within-variance "${rsc_predictor_within_variance}"
        )
      elif [[ "${rsc_predictor_within_variance}" == "known-se" ]]; then
        rsc_pgls_args+=(--predictor-within-variance known-se)
      fi
      if [[ -n "${rsc_predictor_technical_id}" ]]; then
        rsc_pgls_args+=(--predictor-technical-id "${rsc_predictor_technical_id}")
      fi
      if [[ -n "${rsc_predictor_batch}" ]]; then
        rsc_pgls_args+=(--predictor-batch "${rsc_predictor_batch}")
      fi
      if [[ "${rsc_analysis_predictor_se}" != "." ]]; then
        rsc_pgls_args+=(--predictor-standard-error-columns "${rsc_analysis_predictor_se}")
      fi
      if [[ "${rsc_analysis_predictor_n}" != "." ]]; then
        rsc_pgls_args+=(--predictor-sample-size-columns "${rsc_analysis_predictor_n}")
      fi
      if [[ "${rsc_analysis_categorical}" != "." ]]; then
        rsc_pgls_args+=(--categorical-predictors "${rsc_analysis_categorical}")
      fi
      if [[ "${rsc_analysis_ordered}" != "." ]]; then
        rsc_pgls_args+=(--ordered-predictors "${rsc_analysis_ordered}")
      fi
      if [[ "${rsc_analysis_reference}" != "." ]]; then
        rsc_pgls_args+=(--factor-reference "${rsc_analysis_reference}")
      fi
      if [[ "${rsc_origin_applicable}" == "yes" && "${rsc_categorical_origin_diagnostics}" == "stochastic-map" ]]; then
        rsc_pgls_args+=(
          --categorical-origin-diagnostics stochastic-map
          --origin-map-replicates "${rsc_origin_map_replicates}"
          --origin-map-threads "${rsc_origin_map_threads}"
          --origin-min-posterior "${rsc_origin_min_posterior}"
          --origin-leave-one-out "${rsc_origin_leave_one_out}"
        )
      else
        rsc_pgls_args+=(--categorical-origin-diagnostics none)
      fi
      if nwkit "${rsc_pgls_args[@]}" 2>&1 | tee -a "${rsc_log_tmp}"; then
        rsc_analysis_result="ready"
      else
        rsc_nwkit_status=${PIPESTATUS[0]}
        rsc_error_metadata="${rsc_analysis_prefix}.error.tsv"
        python "${gg_support_dir}/reconciled_speciation_contrast.py" inspect-audit-error \
          --audit "${rsc_analysis_audit}" \
          --output "${rsc_error_metadata}"
        rsc_analysis_result=$(rsc_read_metadata_value "${rsc_error_metadata}" status)
        rsc_analysis_reason=$(rsc_read_metadata_value "${rsc_error_metadata}" reason)
        if [[ "${rsc_analysis_result}" != "not_estimable" ]]; then
          echo "Fatal NWKIT failure in RSC analysis ${rsc_analysis_id}: ${rsc_analysis_reason}" >&2
          exit "${rsc_nwkit_status}"
        fi
        echo "RSC analysis ${rsc_analysis_id} is not estimable: ${rsc_analysis_reason}" | tee -a "${rsc_log_tmp}"
        python "${gg_support_dir}/reconciled_speciation_contrast.py" empty-bundle \
          --output-prefix "${rsc_analysis_prefix}" \
          --status-output "${rsc_analysis_prefix}.status.tsv" \
          --audit-output "${rsc_analysis_prefix}.empty.audit.jsonl" \
          --tree-id "${og_id}" \
          --reason "analysis_not_estimable:${rsc_analysis_reason}"
      fi
      printf '%s\t%s\t%s\n' \
        "${rsc_analysis_id}" \
        "${rsc_analysis_prefix}" \
        "${rsc_analysis_audit}" \
        >> "${rsc_bundle_list}"
    done < <(tail -n +2 "${rsc_analysis_plan}")

    python "${gg_support_dir}/reconciled_speciation_contrast.py" aggregate \
      --bundle-list "${rsc_bundle_list}" \
      --output-prefix "${rsc_combined_prefix}" \
      --status-output "${rsc_combined_status}" \
      --audit-output "${rsc_combined_audit}" \
      --source-audit "${rsc_preflight_audit}" \
      --tree-id "${og_id}" \
      --reason "${rsc_preparation_reason}" \
      2>&1 | tee -a "${rsc_log_tmp}"
  fi

  species_pgls_args=(
    --methods "${pgls_methods}"
    --tree-id "${og_id}"
    --species-tree "${species_tree_pruned}"
    --reconciliation "${rsc_preflight_reconciliation}"
    --expression "${rsc_prepared_expression}"
    --species-traits "${rsc_prepared_traits}"
    --analysis-plan "${rsc_analysis_plan}"
    --metadata "${rsc_metadata}"
    --aggregation "${species_expression_aggregation}"
    --expression-value-type "${exp_value_type}"
    --paralog-missing "${species_paralog_missing}"
    --within-variance "${rsc_within_variance}"
    --technical-aggregation "${rsc_technical_aggregation}"
    --predictor-biological-id "${rsc_predictor_biological_id}"
    --predictor-technical-id "${rsc_predictor_technical_id}"
    --predictor-batch "${rsc_predictor_batch}"
    --predictor-within-variance "${rsc_predictor_within_variance}"
    --predictor-technical-aggregation "${rsc_predictor_technical_aggregation}"
    --categorical-replicate-policy "${rsc_categorical_replicate_policy}"
    --factor-coding "${rsc_factor_coding}"
    --branch-length "${rsc_species_branch_length}"
    --response-evolution-model "${rsc_gene_evolution_model}"
    --response-evolution-parameter "${rsc_gene_evolution_parameter}"
    --predictor-evolution-model "${rsc_species_evolution_model}"
    --predictor-evolution-parameter "${rsc_species_evolution_parameter}"
    --predictor-branch-length "${rsc_species_branch_length}"
    --inference "${rsc_inference}"
    --bootstrap-replicates "${rsc_bootstrap_replicates}"
    --seed "${rsc_seed}"
    --confidence-level "${rsc_confidence_level}"
    --reml "${rsc_reml}"
    --allow-large-dense "${rsc_allow_large_dense}"
    --rphylopars-sampling-covariance "${rphylopars_sampling_covariance}"
    --rphylopars-script "${gg_support_dir}/species_tree_rphylopars.R"
    --rsc-results "${rsc_combined_prefix}.pgls.tsv"
    --rsc-status "${rsc_combined_status}"
    --native-out "${species_nwkit_tmp}"
    --rphylopars-out "${species_rphylopars_tmp}"
    --comparison-out "${pgls_comparison_tmp}"
    --status-out "${pgls_method_status_tmp}"
    --audit-out "${pgls_method_audit_tmp}"
    --expression-summary-out "${species_expression_summary_tmp}"
    --expression-audit-out "${species_expression_audit_tmp}"
    --response-tip-summary-out "${species_response_tip_summary_tmp}"
    --response-sampling-covariance-out "${species_response_sampling_covariance_tmp}"
    --predictor-tip-summary-out "${species_predictor_tip_summary_tmp}"
    --predictor-sampling-covariance-out "${species_predictor_sampling_covariance_tmp}"
  )
  if [[ -n "${species_paralog_sampling_covariance}" ]]; then
    species_pgls_args+=(--paralog-sampling-covariance "${species_paralog_sampling_covariance}")
  fi
  if [[ "${rsc_preparation_status}" != "ready" || ! -s "${rsc_preflight_reconciliation}" ]]; then
    species_pgls_args+=(--empty-reason "${rsc_not_estimable_reason:-expression_trait_inputs_not_estimable}")
  fi
  python "${gg_support_dir}/species_tree_pgls.py" "${species_pgls_args[@]}" \
    2>&1 | tee -a "${pgls_log_tmp}"

  mv_out_bundle \
    "${rsc_combined_status}" "${file_og_rsc_status}" \
    "${rsc_combined_prefix}.pgls.tsv" "${file_og_rsc_pgls}" \
    "${rsc_combined_prefix}.reconciliation.tsv" "${file_og_rsc_reconciliation}" \
    "${rsc_combined_prefix}.gene-contrasts.tsv" "${file_og_rsc_gene_contrasts}" \
    "${rsc_combined_prefix}.species-contrasts.tsv" "${file_og_rsc_species_contrasts}" \
    "${rsc_combined_prefix}.response-sampling-covariance.tsv" "${file_og_rsc_response_sampling_covariance}" \
    "${rsc_combined_prefix}.response-tip-summary.tsv" "${file_og_rsc_response_tip_summary}" \
    "${rsc_combined_prefix}.predictor-sampling-covariance.tsv" "${file_og_rsc_predictor_sampling_covariance}" \
    "${rsc_combined_prefix}.predictor-tip-summary.tsv" "${file_og_rsc_predictor_tip_summary}" \
    "${rsc_combined_prefix}.random-effects.tsv" "${file_og_rsc_random_effects}" \
    "${rsc_combined_prefix}.sensitivity.tsv" "${file_og_rsc_sensitivity}" \
    "${rsc_combined_prefix}.trait-origins.tsv" "${file_og_rsc_trait_origins}" \
    "${rsc_combined_audit}" "${file_og_rsc_audit}" \
    "${rsc_log_tmp}" "${file_og_rsc_log}" \
    "${species_nwkit_tmp}" "${file_og_species_nwkit_pgls}" \
    "${species_rphylopars_tmp}" "${file_og_species_rphylopars_pgls}" \
    "${pgls_comparison_tmp}" "${file_og_pgls_comparison}" \
    "${pgls_method_status_tmp}" "${file_og_pgls_method_status}" \
    "${pgls_method_audit_tmp}" "${file_og_pgls_method_audit}" \
    "${species_expression_summary_tmp}" "${file_og_species_expression_summary}" \
    "${species_expression_audit_tmp}" "${file_og_species_expression_audit}" \
    "${species_response_tip_summary_tmp}" "${file_og_species_response_tip_summary}" \
    "${species_response_sampling_covariance_tmp}" "${file_og_species_response_sampling_covariance}" \
    "${species_predictor_tip_summary_tmp}" "${file_og_species_predictor_tip_summary}" \
    "${species_predictor_sampling_covariance_tmp}" "${file_og_species_predictor_sampling_covariance}" \
    "${pgls_log_tmp}" "${file_og_expression_trait_pgls_log}"
  nwkit_version_text=$(nwkit --version 2>&1 | tail -n 1 || true)
  gg_artifact_record \
    "${rsc_provenance_args[@]}" \
    --diagnostic "nwkit_version=${nwkit_version_text:-unknown}" \
    --diagnostic "rphylopars_version=recorded_in_method_status"
  rm -rf -- "${rsc_work_dir}"
else
  gg_step_skip "${task}"
fi

task="IQ-TREE ancestral codon sequence reconstruction for CSUBST"
disable_if_no_input_file "run_iqtree_anc" "${file_og_trimmed_aln_analysis}" "${file_og_rooted_tree_analysis}"
iqtree_anc_needs_update=0
iqtree_anc_provenance_args=(
  --manifest "${file_og_iqtree_anc_provenance}"
  --step "iqtree_anc"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "trimmed_alignment=${file_og_trimmed_aln_analysis}"
  --input "rooted_tree=${file_og_rooted_tree_analysis}"
  --output "iqtree_anc=${file_og_iqtree_anc}"
  --parameter "codon_model=${codon_model}"
  --parameter "genetic_code=${genetic_code}"
)
gg_artifact_prepare_stage iqtree_anc_needs_update run_iqtree_anc "${iqtree_anc_provenance_args[@]}" || exit $?
if [[ ${iqtree_anc_needs_update} -eq 1 && ${run_iqtree_anc} -eq 1 ]]; then
  gg_step_start "${task}"

  shopt -s nullglob
  csubst_cleanup_paths=(csubst.* csubst_search csubst_scan)
  shopt -u nullglob
  if [[ ${#csubst_cleanup_paths[@]} -gt 0 ]]; then
    rm -rf -- "${csubst_cleanup_paths[@]}"
  fi
  seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file "tmp.csubst.fasta"

  nwkit drop --target intnode --support yes --name no --infile "${file_og_rooted_tree_analysis}" |
    nwkit intersection --seqin tmp.csubst.fasta --seqout csubst.fasta |
    nwkit sanitize --remove-singleton yes --resolve-polytomy no \
      > csubst.nwk

  rm -f -- tmp.csubst.fasta

  build_iqtree_mem_args
  iqtree \
    -s csubst.fasta \
    -te csubst.nwk \
    -m "${codon_model}" \
    -T AUTO \
    --threads-max "${GG_TASK_CPUS}" \
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
    csubst_outputs=(csubst.*)
    shopt -u nullglob
    if [[ ${#csubst_outputs[@]} -eq 0 ]]; then
      echo "Expected csubst output files were not found."
      exit 1
    fi
    mv_out "${csubst_outputs[@]}" "${og_id}.iqtree.anc"
    rm -f -- "${og_id}.iqtree.anc.zip"
    zip -rq "${og_id}.iqtree.anc.zip" "${og_id}.iqtree.anc"
    python "${gg_support_dir}/atomic_zip_publish.py" \
      --source "${og_id}.iqtree.anc.zip" \
      --destination "${file_og_iqtree_anc}" \
      --expected-prefix "${og_id}.iqtree.anc" \
      --remove-source
    rm -rf -- "${og_id}.iqtree.anc"
    iqtree_version_text=$(iqtree --version 2>&1 || true)
    iqtree_version_text=${iqtree_version_text%%$'\n'*}
    gg_artifact_record \
      "${iqtree_anc_provenance_args[@]}" \
      --diagnostic "genegalleon_version=${GG_VERSION:-${SINGULARITYENV_GG_VERSION:-${APPTAINERENV_GG_VERSION:-unknown}}}" \
      --diagnostic "container_image=${gg_container_image_path:-unknown}" \
      --diagnostic "iqtree_version=${iqtree_version_text:-unknown}"
    iqtree_anc_needs_update=0
  fi
else
  gg_step_skip "${task}"
fi

if [[ ${iqtree_anc_needs_update} -eq 1 ]] && { [[ ${run_csubst} -eq 1 ]] || [[ ${run_csubst_scan} -eq 1 ]]; }; then
  echo "Refusing to consume stale or unverified IQ-TREE ancestral output: ${file_og_iqtree_anc}" >&2
  echo "Set run_iqtree_anc=1 so it can be regenerated from the declared alignment and rooted tree." >&2
  exit 1
fi

task="CSUBST"
disable_if_no_input_file "run_csubst" "${file_og_iqtree_anc}"
csubst_needs_update=0
csubst_provenance_args=(
  --manifest "${file_og_csubst_provenance}"
  --step "csubst"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "iqtree_anc=${file_og_iqtree_anc}"
  --output "csubst_b=${file_og_csubst_b}"
  --output "csubst_cb_stats=${file_og_csubst_cb_stats}"
  --parameter "codon_model=${codon_model}"
  --parameter "genetic_code=${genetic_code}"
  --parameter "max_arity=${csubst_max_arity}"
  --parameter "exhaustive_until=${csubst_exhaustive_until}"
  --parameter "cutoff_stat=${csubst_cutoff_stat}"
  --parameter "max_combination=${csubst_max_combination}"
  --parameter "fg_exclude_wg=${csubst_fg_exclude_wg}"
  --parameter "fg_stem_only=${csubst_fg_stem_only}"
  --parameter "nonsyn_recode=${csubst_nonsyn_recode}"
)
if [[ ${csubst_max_arity} -ge 2 ]]; then
  for ((i = 2; i <= csubst_max_arity; i++)); do
    csubst_cb_varname="file_og_csubst_cb_${i}"
    csubst_provenance_args+=(
      --optional-output "csubst_cb_${i}=${!csubst_cb_varname}"
    )
  done
fi
if [[ -s "${file_sp_trait}" ]]; then
  csubst_provenance_args+=(--input "foreground_trait=${file_sp_trait}")
  csubst_provenance_args+=(--parameter "foreground_present=1")
else
  csubst_provenance_args+=(--parameter "foreground_present=0")
fi
gg_artifact_prepare_stage csubst_needs_update run_csubst "${csubst_provenance_args[@]}" || exit $?
if [[ ${csubst_needs_update} -eq 1 && ${run_csubst} -eq 1 ]]; then
  gg_step_start "${task}"

  # csubst_cb_K files are managed optional outputs: a missing arity can be a
  # valid result, but a file from an older generation must never survive and
  # be imported into the gene-family database as if it were current.
  shopt -s nullglob
  old_csubst_cb_files=(
    "${dir_output_active}"/csubst_cb_[0-9]*/"${og_id}"_csubst_cb_*.tsv
    "${dir_output_active}"/csubst_cb_[0-9]*/"${og_id}".csubst_cb_*.tsv
  )
  shopt -u nullglob
  for old_csubst_cb_file in "${old_csubst_cb_files[@]}"; do
    old_csubst_cb_logical_path=${old_csubst_cb_file#"${dir_output_active}"/}
    if ! python "${gene_family_store_script}" delete \
      --root "${dir_output_active}" \
      --path "${old_csubst_cb_logical_path}" \
      --family-id "${og_id}"; then
      echo "Failed to remove stale managed CSUBST output: ${old_csubst_cb_file}" >&2
      exit 1
    fi
  done

  if [[ -s "${file_sp_trait}" ]]; then
    echo "CSUBST foreground specification file: ${file_sp_trait}"
    first_trait_header=""
    IFS= read -r first_trait_header < "${file_sp_trait}" || true
    if [[ "${first_trait_header}" == *" "* ]]; then
      echo "Column names should not contain spaces: ${file_sp_trait}"
      echo "Exiting."
      exit 1
    fi
    write_csubst_foreground_table "${file_sp_trait}" "foreground.tsv"
    foreground_params=(--foreground foreground.tsv --fg_format 2)
  else
    echo 'Foreground specification file was not found. CSUBST will run without it.'
    foreground_params=()
  fi
  shopt -s nullglob
  csubst_cleanup_paths=(csubst.*)
  shopt -u nullglob
  if [[ ${#csubst_cleanup_paths[@]} -gt 0 ]]; then
    rm -rf -- "${csubst_cleanup_paths[@]}"
  fi
  if [[ -z "${og_id:-}" ]]; then
    echo "og_id is empty. Refusing to run cleanup pattern."
    exit 1
  fi
  shopt -s nullglob
  og_cleanup_paths=("${og_id}".*)
  shopt -u nullglob
  if [[ ${#og_cleanup_paths[@]} -gt 0 ]]; then
    rm -rf -- "${og_cleanup_paths[@]}"
  fi
  gg_extract_expected_zip_prefix \
    "${file_og_iqtree_anc}" \
    "${og_id}.iqtree.anc"
  csubst_input_base="./${og_id}.iqtree.anc/csubst"
  csubst_search_dir="csubst_search"

  csubst search \
    --genetic_code "${genetic_code}" \
    --alignment_file "${csubst_input_base}.fasta" \
    --rooted_tree_file "${csubst_input_base}.nwk" \
    --iqtree_treefile "${csubst_input_base}.treefile" \
    --iqtree_state "${csubst_input_base}.state" \
    --iqtree_rate "${csubst_input_base}.rate" \
    --iqtree_iqtree "${csubst_input_base}.iqtree" \
    --iqtree_log "${csubst_input_base}.log" \
    --iqtree_model "${codon_model}" \
    --max_arity "${csubst_max_arity}" \
    --exhaustive_until "${csubst_exhaustive_until}" \
    --cutoff_stat "${csubst_cutoff_stat}" \
    --max_combination "${csubst_max_combination}" \
    --fg_exclude_wg "${csubst_fg_exclude_wg}" \
    --fg_stem_only "${csubst_fg_stem_only}" \
    --nonsyn_recode "${csubst_nonsyn_recode}" \
    --expectation_method "codon_model" \
    --threads "${GG_TASK_CPUS}" \
    "${foreground_params[@]}"

  csubst_b_src=""
  csubst_cb_stats_src=""
  for candidate in "${csubst_search_dir}/csubst_b.tsv" "csubst_b.tsv"; do
    if [[ -s "${candidate}" ]]; then
      csubst_b_src="${candidate}"
      break
    fi
  done
  for candidate in "${csubst_search_dir}/csubst_cb_stats.tsv" "csubst_cb_stats.tsv"; do
    if [[ -s "${candidate}" ]]; then
      csubst_cb_stats_src="${candidate}"
      break
    fi
  done

  if [[ -n "${csubst_b_src}" && -n "${csubst_cb_stats_src}" ]]; then
    echo "CSUBST was successful."
    mv_out "${csubst_b_src}" "${file_og_csubst_b}"
    mv_out "${csubst_cb_stats_src}" "${file_og_csubst_cb_stats}"
    if [[ ${csubst_max_arity} -ge 2 ]]; then
      for ((i = 2; i <= csubst_max_arity; i++)); do
        csubst_cb_src=""
        for candidate in "${csubst_search_dir}/csubst_cb_${i}.tsv" "csubst_cb_${i}.tsv"; do
          if [[ -e "${candidate}" ]]; then
            csubst_cb_src="${candidate}"
            break
          fi
        done
        if [[ -n "${csubst_cb_src}" ]]; then
          my_csubst_file=file_og_csubst_cb_${i}
          mv_out "${csubst_cb_src}" "${!my_csubst_file}"
        fi
      done
    fi
    csubst_version_text=$(csubst --version 2>&1 || true)
    csubst_version_text=${csubst_version_text%%$'\n'*}
    gg_artifact_record \
      "${csubst_provenance_args[@]}" \
      --diagnostic "genegalleon_version=${GG_VERSION:-${SINGULARITYENV_GG_VERSION:-${APPTAINERENV_GG_VERSION:-unknown}}}" \
      --diagnostic "container_image=${gg_container_image_path:-unknown}" \
      --diagnostic "csubst_version=${csubst_version_text:-unknown}"
    csubst_needs_update=0
  else
    echo "CSUBST failed."
  fi
else
  gg_step_skip "${task}"
fi

if [[ ${csubst_needs_update} -eq 1 && ${run_csubst} -ne 1 && ${run_summary} -eq 1 && -s "${file_og_csubst_b}" ]]; then
  echo "Refusing to include stale or unverified CSUBST output in stat_branch: ${file_og_csubst_b}" >&2
  echo "Set run_csubst=1 so it can be regenerated from the declared IQ-TREE archive and parameters." >&2
  exit 1
fi

task="CSUBST scan"
disable_if_no_input_file "run_csubst_scan" "${file_og_iqtree_anc}" "${file_sp_trait}"
csubst_scan_needs_update=0
csubst_scan_provenance_args=(
  --manifest "${file_og_csubst_scan_provenance}"
  --step "csubst_scan"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --input "iqtree_anc=${file_og_iqtree_anc}"
  --input "foreground_trait=${file_sp_trait}"
  --output "csubst_scan=${file_og_csubst_scan}"
  --output "csubst_scan_units=${file_og_csubst_scan_units}"
  --parameter "codon_model=${codon_model}"
  --parameter "genetic_code=${genetic_code}"
  --parameter "scan_unit_mode=${csubst_scan_unit_mode}"
  --parameter "scan_match=${csubst_scan_match}"
  --parameter "scan_min_event_pp=${csubst_scan_min_event_pp}"
  --parameter "scan_min_support=${csubst_scan_min_support}"
  --parameter "scan_rate_event_mode=${csubst_scan_rate_event_mode}"
  --parameter "scan_rate_length=${csubst_scan_rate_length}"
  --parameter "scan_rate_exposure=${csubst_scan_rate_exposure}"
  --parameter "scan_other_scope=${csubst_scan_other_scope}"
  --parameter "scan_pvalue_calibration=${csubst_scan_pvalue_calibration}"
  --parameter "scan_n_permutations=${csubst_scan_n_permutations}"
  --parameter "scan_site_plot=${csubst_scan_site_plot}"
  --parameter "tree_site_plot_format=${csubst_scan_tree_site_plot_format}"
  --parameter "tree_site_plot_max_sites=${csubst_scan_tree_site_plot_max_sites}"
  --parameter "nonsyn_recode=${csubst_nonsyn_recode}"
)
if [[ "${csubst_scan_site_plot}" == "yes" ]]; then
  csubst_scan_provenance_args+=(--output "csubst_scan_plot=${file_og_csubst_scan_plot}")
fi
gg_artifact_prepare_stage csubst_scan_needs_update run_csubst_scan "${csubst_scan_provenance_args[@]}" || exit $?
if [[ ${csubst_scan_needs_update} -eq 1 && ${run_csubst_scan} -eq 1 ]]; then
  gg_step_start "${task}"

  echo "CSUBST scan foreground specification file: ${file_sp_trait}"
  first_trait_header=""
  IFS= read -r first_trait_header < "${file_sp_trait}" || true
  if [[ "${first_trait_header}" == *" "* ]]; then
    echo "Column names should not contain spaces: ${file_sp_trait}"
    echo "Exiting."
    exit 1
  fi
  write_csubst_foreground_table "${file_sp_trait}" "foreground.tsv"

  if ! csubst scan -h >/dev/null 2>&1; then
    echo "csubst scan is unavailable in this runtime. Rebuild the GeneGalleon container with a csubst version that provides the scan subcommand."
    exit 1
  fi

  shopt -s nullglob
  csubst_cleanup_paths=(csubst.* csubst_scan)
  shopt -u nullglob
  if [[ ${#csubst_cleanup_paths[@]} -gt 0 ]]; then
    rm -rf -- "${csubst_cleanup_paths[@]}"
  fi
  if [[ -e "${og_id}.iqtree.anc" ]]; then
    rm -rf -- "${og_id}.iqtree.anc"
  fi
  gg_extract_expected_zip_prefix \
    "${file_og_iqtree_anc}" \
    "${og_id}.iqtree.anc"
  csubst_input_base="./${og_id}.iqtree.anc/csubst"
  csubst_scan_dir="csubst_scan"

  csubst scan \
    --genetic_code "${genetic_code}" \
    --alignment_file "${csubst_input_base}.fasta" \
    --rooted_tree_file "${csubst_input_base}.nwk" \
    --iqtree_treefile "${csubst_input_base}.treefile" \
    --iqtree_state "${csubst_input_base}.state" \
    --iqtree_rate "${csubst_input_base}.rate" \
    --iqtree_iqtree "${csubst_input_base}.iqtree" \
    --iqtree_log "${csubst_input_base}.log" \
    --iqtree_model "${codon_model}" \
    --foreground foreground.tsv \
    --fg_format 2 \
    --scan_unit_mode "${csubst_scan_unit_mode}" \
    --scan_match "${csubst_scan_match}" \
    --scan_min_event_pp "${csubst_scan_min_event_pp}" \
    --scan_min_support "${csubst_scan_min_support}" \
    --scan_rate_event_mode "${csubst_scan_rate_event_mode}" \
    --scan_rate_length "${csubst_scan_rate_length}" \
    --scan_rate_exposure "${csubst_scan_rate_exposure}" \
    --scan_other_scope "${csubst_scan_other_scope}" \
    --scan_pvalue_calibration "${csubst_scan_pvalue_calibration}" \
    --scan_n_permutations "${csubst_scan_n_permutations}" \
    --scan_site_plot "${csubst_scan_site_plot}" \
    --tree_site_plot_format "${csubst_scan_tree_site_plot_format}" \
    --tree_site_plot_max_sites "${csubst_scan_tree_site_plot_max_sites}" \
    --nonsyn_recode "${csubst_nonsyn_recode}" \
    --threads "${GG_TASK_CPUS}" \
    --outdir "${csubst_scan_dir}" \
    --output_prefix csubst

  if [[ -s "${csubst_scan_dir}/csubst_scan.tsv" && -s "${csubst_scan_dir}/csubst_scan_units.tsv" ]]; then
    echo "CSUBST scan was successful."
    mv_out "${csubst_scan_dir}/csubst_scan.tsv" "${file_og_csubst_scan}"
    mv_out "${csubst_scan_dir}/csubst_scan_units.tsv" "${file_og_csubst_scan_units}"
    shopt -s nullglob
    csubst_scan_foreground_branch_files=("${csubst_scan_dir}"/csubst_foreground_branch*.txt)
    shopt -u nullglob
    for csubst_scan_foreground_branch_file in "${csubst_scan_foreground_branch_files[@]}"; do
      csubst_scan_foreground_branch_name=$(basename "${csubst_scan_foreground_branch_file}")
      csubst_scan_foreground_branch_suffix="${csubst_scan_foreground_branch_name#csubst_foreground_branch}"
      mv_out \
        "${csubst_scan_foreground_branch_file}" \
        "${dir_output_active}/csubst_scan_foreground_branch/${og_id}_csubst_foreground_branch${csubst_scan_foreground_branch_suffix}"
    done
    if [[ -e "${csubst_scan_dir}/csubst_scan.tree_site.${csubst_scan_tree_site_plot_format}" ]]; then
      mv_out "${csubst_scan_dir}/csubst_scan.tree_site.${csubst_scan_tree_site_plot_format}" "${file_og_csubst_scan_plot}"
    fi
    if [[ -e "${csubst_scan_dir}/csubst.log" ]]; then
      mv_out "${csubst_scan_dir}/csubst.log" "${file_og_csubst_scan_log}"
    fi
    csubst_version_text=$(csubst --version 2>&1 || true)
    csubst_version_text=${csubst_version_text%%$'\n'*}
    gg_artifact_record \
      "${csubst_scan_provenance_args[@]}" \
      --diagnostic "genegalleon_version=${GG_VERSION:-${SINGULARITYENV_GG_VERSION:-${APPTAINERENV_GG_VERSION:-unknown}}}" \
      --diagnostic "container_image=${gg_container_image_path:-unknown}" \
      --diagnostic "csubst_version=${csubst_version_text:-unknown}"
    csubst_scan_needs_update=0
  else
    echo "CSUBST scan failed."
  fi
  rm -rf -- "${og_id}.iqtree.anc" "${csubst_scan_dir}"
else
  gg_step_skip "${task}"
fi
# shellcheck shell=bash
# Sourced by gg_gene_evolution_core.sh.

task="summary statistics"
summary_input_files=(
  "${species_tree_pruned}"
  "${file_og_primary_fasta}"
  "${file_og_trimmed_aln_analysis}"
  "${file_og_unrooted_tree_analysis}"
  "${file_og_rooted_tree_analysis}"
  "${file_og_rooted_log}"
  "${file_og_notung_reconcil}"
  "${file_og_dated_tree_analysis}"
  "${file_og_dated_tree_log}"
  "${file_og_generax_nhx}"
  "${file_og_hyphy_dnds}"
  "${file_og_hyphy_relax}"
  "${file_og_hyphy_relax_reversed}"
  "${file_og_l1ou_fit_tree}"
  "${file_og_l1ou_fit_regime}"
  "${file_og_l1ou_fit_leaf}"
  "${file_og_expression}"
  "${file_og_mapdnds_dn}"
  "${file_og_mapdnds_ds}"
  "${file_og_codeml_two_ratio}"
  "${file_og_gff_info}"
  "${file_og_fimo}"
  "${file_og_promoter_fasta}"
  "${file_og_scm_intron_summary}"
  "${file_og_csubst_b}"
  "${file_og_gene_pgls}"
  "${file_og_pgls_comparison}"
  "${file_og_pgls_method_status}"
  "${file_og_rsc_status}"
  "${file_og_rsc_pgls}"
  "${file_og_rpsblast}"
  "${file_og_uniprot_annotation}"
  "${file_og_cdskit_localize}"
  "${file_og_synteny}"
)
task="Synteny neighborhood grouping"
if [[ ${treevis_synteny} -eq 1 ]] && { [[ ${run_summary} -eq 1 ]] || [[ ${run_tree_plot} -eq 1 ]]; }; then
  synteny_source_dir="${dir_sp_cds}"
  synteny_sequence_mode="cds"
  if [[ "${input_sequence_mode}" == "protein" ]] && species_protein_input_has_files; then
    synteny_source_dir="${dir_sp_protein_input}"
    synteny_sequence_mode="protein"
  fi
  run_synteny_generation=1
  synteny_needs_update=0
  synteny_provenance_args=(
    --manifest "${dir_output_active}/artifact_provenance/${og_id}.synteny.json"
    --step "synteny"
    --family-id "${og_id}"
    --logical-root "${dir_output_active}"
    --workspace-root "${gg_workspace_dir}"
    --input "primary_fasta=${file_og_primary_fasta}"
    --optional-output "synteny=${file_og_synteny}"
    --parameter "input_sequence_mode=${synteny_sequence_mode}"
    --parameter "window=${treevis_synteny_window}"
    --parameter "query_blast_evalue=${query_blast_evalue}"
    --parameter "auto_evalue_cutoffs=${query_blast_auto_evalue_maxlen_cutoffs}"
    --parameter "genetic_code=${genetic_code}"
    --parameter "sequence_dir_present=$([[ -d "${synteny_source_dir}" ]] && echo 1 || echo 0)"
    --parameter "gff_dir_present=$([[ -d "${dir_sp_gff}" ]] && echo 1 || echo 0)"
  )
  if [[ -d "${synteny_source_dir}" ]]; then
    synteny_provenance_args+=(--input "species_sequences=${synteny_source_dir}")
  fi
  if [[ -d "${dir_sp_gff}" ]]; then
    synteny_provenance_args+=(--input "species_gff=${dir_sp_gff}")
  fi
  gg_artifact_prepare_stage synteny_needs_update run_synteny_generation "${synteny_provenance_args[@]}" || exit $?
  if [[ ${synteny_needs_update} -eq 1 ]]; then
    gg_step_start "${task}"
    rm -f -- "${file_og_synteny}"
    if [[ ! -d "${synteny_source_dir}" ]]; then
      echo "Sequence directory not found. Skipping synteny panel input generation: ${synteny_source_dir}"
    elif [[ ! -d "${dir_sp_gff}" ]]; then
      echo "species_gff directory not found. Skipping synteny panel input generation: ${dir_sp_gff}"
    elif [[ ! -s "${file_og_primary_fasta}" ]]; then
      echo "Focal sequence fasta file not found. Skipping synteny panel input generation: ${file_og_primary_fasta}"
    else
      synteny_evalue="${query_blast_evalue}"
      if [[ "$(printf '%s' "${query_blast_evalue}" | tr '[:upper:]' '[:lower:]')" == "auto" ]]; then
        synteny_evalue_query_fasta="${og_id}.synteny.evalue.aa.tmp.fasta"
        if ! prepare_synteny_evalue_fasta "${synteny_evalue_query_fasta}"; then
          rm -f -- "${synteny_evalue_query_fasta}"
          echo "Failed to prepare amino-acid FASTA for synteny query_blast_evalue=auto: ${file_og_primary_fasta}"
          exit 1
        fi
        if ! resolve_query_blast_evalue "${query_blast_evalue}" "${synteny_evalue_query_fasta}" "${query_blast_auto_evalue_maxlen_cutoffs}"; then
          rm -f -- "${synteny_evalue_query_fasta}"
          exit 1
        fi
        synteny_evalue="${effective_query_blast_evalue}"
        echo "synteny auto E-value: query_count=${query_blast_query_num_seqs} min_aa_len=${query_blast_query_min_aa_len} avg_aa_len=${query_blast_query_avg_aa_len} max_aa_len=${query_blast_query_max_aa_len}"
        echo "synteny auto E-value: cutoffs=${query_blast_auto_evalue_maxlen_cutoffs} effective_synteny_evalue=${synteny_evalue}"
        rm -f -- "${synteny_evalue_query_fasta}"
      fi
      python "${gg_support_dir}/synteny_neighbors.py" \
        --focal_cds_fasta "${file_og_primary_fasta}" \
        --dir_sp_cds "${synteny_source_dir}" \
        --dir_sp_gff "${dir_sp_gff}" \
        --cache_dir "${gg_workspace_output_dir}/species_gff_info" \
        --lock_dir "${file_og_parameters_dir}/synteny_locks" \
        --gff2genestat_script "${gg_support_dir}/gff2genestat.py" \
        --input_sequence_mode "${synteny_sequence_mode}" \
        --window "${treevis_synteny_window}" \
        --evalue "${synteny_evalue}" \
        --genetic_code "${genetic_code}" \
        --threads "${GG_TASK_CPUS}" \
        --outfile "${file_og_synteny}"
      if [[ ! -s "${file_og_synteny}" ]]; then
        echo "No synteny links were generated: ${file_og_synteny}"
      fi
    fi
    gg_artifact_record "${synteny_provenance_args[@]}"
  else
    gg_step_skip "${task}"
  fi
else
  gg_step_skip "${task}"
fi
task="summary statistics"
summary_needs_update=0
summary_outputs_changed=0
summary_provenance_args=(
  --manifest "${file_og_summary_provenance}"
  --step "summary_statistics"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --output "stat_branch=${file_og_stat_branch}"
  --output "stat_tree=${file_og_stat_tree}"
  --parameter "clade_ortholog_prefix=${treevis_clade_ortholog_prefix}"
  --parameter "mode_gene_evolution=${mode_gene_evolution}"
  --parameter "treevis_query_marker=${treevis_query_marker}"
  --parameter "query_blast_coverage=${query_blast_coverage}"
)
for summary_input_index in "${!summary_input_files[@]}"; do
  summary_input=${summary_input_files[${summary_input_index}]}
  if [[ -e "${summary_input}" ]]; then
    summary_provenance_args+=(--input "input_${summary_input_index}=${summary_input}")
  fi
done
if [[ "${mode_gene_evolution}" == "query2family" ]]; then
  for query_summary_input_spec in \
    "query_gene=${file_query_gene:-}" \
    "query_aa_fasta=${file_og_query_aa_fasta}" \
    "query_blast=${file_og_query_blast}"
  do
    query_summary_input_path=${query_summary_input_spec#*=}
    if [[ -e "${query_summary_input_path}" ]]; then
      summary_provenance_args+=(--input "${query_summary_input_spec}")
    fi
  done
fi
gg_artifact_prepare_stage summary_needs_update run_summary "${summary_provenance_args[@]}" || exit $?
disable_if_no_input_file "run_summary" "${file_og_rooted_tree_analysis}"
if [[ ${summary_needs_update} -eq 1 && ${run_summary} -eq 1 ]]; then
  gg_step_start "${task}"

  if [[ -s "${file_og_notung_reconcil}" ]]; then
    gg_extract_expected_zip_prefix \
      "${file_og_notung_reconcil}" \
      "${og_id}.notung_reconcile"
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
  if [[ -s "${file_og_primary_fasta}" ]]; then
    summary_unaligned_fasta="${og_id}.summary.unaligned.fasta"
    seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_primary_fasta}" --out-file "${summary_unaligned_fasta}"
    summary_tmp_files+=("${summary_unaligned_fasta}")
  fi
  if [[ -s "${file_og_trimmed_aln_analysis}" ]]; then
    summary_trimmed_fasta="${og_id}.summary.trimmed.fasta"
    seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_trimmed_aln_analysis}" --out-file "${summary_trimmed_fasta}"
    summary_tmp_files+=("${summary_trimmed_fasta}")
  fi
  if [[ -s "${file_og_promoter_fasta}" ]]; then
    summary_promoter_fasta="${og_id}.summary.promoter.fasta"
    seqkit seq --threads "${GG_TASK_CPUS}" "${file_og_promoter_fasta}" --out-file "${summary_promoter_fasta}"
    summary_tmp_files+=("${summary_promoter_fasta}")
  fi

  python "${gg_support_dir}/orthogroup_statistics.py" \
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
    --pgls_comparison "${file_og_pgls_comparison}" \
    --pgls_method_status "${file_og_pgls_method_status}" \
    --rsc_status "${file_og_rsc_status}" \
    --rsc_pgls "${file_og_rsc_pgls}" \
    --rpsblast "${file_og_rpsblast}" \
    --uniprot "${file_og_uniprot_annotation}" \
    --cdskit_localize "${file_og_cdskit_localize}" \
    --synteny "${file_og_synteny}" \
    --ncpu "${GG_TASK_CPUS}" \
    --clade_ortholog_prefix "${treevis_clade_ortholog_prefix}"
  if [[ ${#summary_tmp_files[@]} -gt 0 ]]; then
    rm -f -- "${summary_tmp_files[@]}"
  fi

  #--csubst_cb_stats ${file_og_csubst_cb_stats} \ # Does not support --arity 3 or larger

  cp_out orthogroup.branch.tsv "${file_og_stat_branch}"
  cp_out orthogroup.tree.tsv "${file_og_stat_tree}"
  summary_outputs_changed=1

else
  gg_step_skip "${task}"
fi

if [[ ${summary_needs_update} -eq 1 && ${run_summary} -ne 1 ]] && { [[ ${run_tree_plot} -eq 1 ]] || [[ "${mode_gene_evolution}" == "query2family" && ${treevis_query_marker} -eq 1 ]]; }; then
  echo "Refusing to consume stale or unverified summary outputs: ${file_og_stat_branch}, ${file_og_stat_tree}" >&2
  echo "Set run_summary=1 so they can be regenerated from the declared inputs and parameters." >&2
  exit 1
fi

task="Query marker annotation"
if [[ "${mode_gene_evolution}" == "query2family" && ${treevis_query_marker} -eq 1 && -s "${file_og_stat_branch}" ]]; then
  query_marker_needs_update=0
  if ! awk -F $'\t' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "query_marker") found = 1; exit(found ? 0 : 1) }' "${file_og_stat_branch}"; then
    query_marker_needs_update=1
  fi
  if [[ ${query_marker_needs_update} -eq 1 ]]; then
    gg_step_start "${task}"
    python "${gg_support_dir}/annotate_stat_branch_query_markers.py" \
      --stat_branch "${file_og_stat_branch}" \
      --query_gene "${file_query_gene}" \
      --query_aa_fasta "${file_og_query_aa_fasta}" \
      --query_blast "${file_og_query_blast}" \
      --min_query_blast_coverage "${query_blast_coverage}" \
      --outfile "orthogroup.branch.query_marker.tsv"
    cp_out "orthogroup.branch.query_marker.tsv" "${file_og_stat_branch}"
    summary_needs_update=1
    summary_outputs_changed=1
  else
    gg_step_skip "${task}"
  fi
else
  gg_step_skip "${task}"
fi

if [[ ${summary_outputs_changed} -eq 1 ]]; then
  gg_artifact_record \
    "${summary_provenance_args[@]}" \
    --diagnostic "genegalleon_version=${GG_VERSION:-${SINGULARITYENV_GG_VERSION:-${APPTAINERENV_GG_VERSION:-unknown}}}" \
    --diagnostic "container_image=${gg_container_image_path:-unknown}"
fi

if [[ -s "${file_og_iqtree_anc}" && -s "${file_og_stat_branch}" ]]; then
  python "${gg_support_dir}/validate_csubst_branch_identity.py" \
    --iqtree-anc "${file_og_iqtree_anc}" \
    --stat-branch "${file_og_stat_branch}"
fi

task="stat_branch2tree_plot"
disable_if_no_input_file "run_tree_plot" "${file_og_stat_branch}" "${file_og_stat_tree}"
tree_plot_needs_update=0
tree_plot_provenance_args=(
  --manifest "${file_og_tree_plot_provenance}"
  --step "tree_plot"
  --family-id "${og_id}"
  --logical-root "${dir_output_active}"
  --workspace-root "${gg_workspace_dir}"
  --output "tree_plot=${file_og_tree_plot}"
  --parameter "branch_length=${treevis_branch_length}"
  --parameter "support_value=${treevis_support_value}"
  --parameter "branch_color=${treevis_branch_color}"
  --parameter "heatmap_transform=${treevis_heatmap_transform}"
  --parameter "max_intergenic_dist=${treevis_max_intergenic_dist}"
  --parameter "synteny_window=${treevis_synteny_window}"
  --parameter "query_marker=${treevis_query_marker}"
  --parameter "retrotransposition_delta_intron=${treevis_retrotransposition_delta_intron}"
  --parameter "clade_ortholog=${treevis_clade_ortholog}"
  --parameter "event_method=${treevis_event_method}"
  --parameter "pie_chart_value_transformation=${treevis_pie_chart_value_transformation}"
  --parameter "long_branch_display=${treevis_long_branch_display}"
  --parameter "long_branch_ref_quantile=${treevis_long_branch_ref_quantile}"
  --parameter "long_branch_detect_ratio=${treevis_long_branch_detect_ratio}"
  --parameter "long_branch_cap_ratio=${treevis_long_branch_cap_ratio}"
  --parameter "long_branch_tail_shrink=${treevis_long_branch_tail_shrink}"
  --parameter "long_branch_max_fraction=${treevis_long_branch_max_fraction}"
  --parameter "csubst_max_arity=${csubst_max_arity}"
  --parameter "csubst_cutoff_stat=${csubst_cutoff_stat}"
  --parameter "promoter_bp=${promoter_bp}"
  --parameter "fimo_qvalue=${fimo_qvalue}"
  --parameter "species_label_parser=${species_label_parser}"
)
tree_plot_input_files=(
  "${file_og_stat_branch}"
  "${file_og_stat_tree}"
  "${file_og_synteny}"
  "${file_og_trimmed_aln_analysis}"
  "${file_og_untrimmed_aln_analysis}"
  "${file_og_clipkit}"
  "${file_og_orthogroup_extraction_fasta}"
  "${file_og_maxalign}"
  "${file_og_mafft}"
  "${file_og_primary_fasta}"
  "${file_og_rpsblast}"
  "${file_og_meme}"
  "${file_og_dated_tree}"
)
for ((i = 2; i <= csubst_max_arity; i++)); do
  csubst_cb_varname="file_og_csubst_cb_${i}"
  if [[ -n "${!csubst_cb_varname:-}" ]]; then
    tree_plot_input_files+=("${!csubst_cb_varname}")
  fi
done
for tree_plot_input_index in "${!tree_plot_input_files[@]}"; do
  tree_plot_input=${tree_plot_input_files[${tree_plot_input_index}]}
  if [[ -e "${tree_plot_input}" ]]; then
    tree_plot_provenance_args+=(--input "input_${tree_plot_input_index}=${tree_plot_input}")
  fi
done
gg_artifact_prepare_stage tree_plot_needs_update run_tree_plot "${tree_plot_provenance_args[@]}" || exit $?
if [[ ${run_tree_plot} -eq 1 ]]; then
  if ! Rscript -e "if (!requireNamespace('ggimage', quietly=TRUE)) quit(status=1)" > /dev/null 2>&1; then
    echo "ggimage package is unavailable. Disabling run_tree_plot."
    run_tree_plot=0
  fi
fi
if [[ ${tree_plot_needs_update} -eq 1 && ${run_tree_plot} -eq 1 ]]; then
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
      "${file_og_primary_fasta}"; do
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
      "${file_og_primary_fasta}"; do
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

  tree_plot_panel_args=(
    "--panel1=tree,${treevis_branch_length},${treevis_support_value},${treevis_branch_color},L"
    "--panel2=heatmap,${treevis_heatmap_transform},abs,_,expression_"
    "--panel3=pointplot,no,rel,_,expression_"
    "--panel4=cluster_membership,${treevis_max_intergenic_dist}"
    "--panel5=synteny,${file_og_synteny},${treevis_synteny_window}"
    "--panel6=tiplabel"
  )
  panel_index=7
  if [[ "${mode_gene_evolution}" == "query2family" && ${treevis_query_marker} -eq 1 ]]; then
    tree_plot_panel_args+=("--panel${panel_index}=categorical,query_marker,Query,-")
    panel_index=$((panel_index + 1))
  fi
  tree_plot_panel_args+=(
    "--panel${panel_index}=signal_peptide"
  )
  panel_index=$((panel_index + 1))
  tree_plot_panel_args+=(
    "--panel${panel_index}=transmembrane_domain"
  )
  panel_index=$((panel_index + 1))
  tree_plot_panel_args+=(
    "--panel${panel_index}=intron_number"
  )
  panel_index=$((panel_index + 1))
  tree_plot_panel_args+=(
    "--panel${panel_index}=domain,${file_og_rpsblast}"
  )
  panel_index=$((panel_index + 1))
  tree_plot_panel_args+=(
    "--panel${panel_index}=alignment,${panel11_trimmed_aln},${panel11_untrimmed_aln}"
  )
  panel_index=$((panel_index + 1))
  tree_plot_panel_args+=(
    "--panel${panel_index}=fimo,${promoter_bp},${fimo_qvalue}"
  )
  panel_index=$((panel_index + 1))
  tree_plot_panel_args+=(
    "--panel${panel_index}=meme,${file_og_meme}"
  )
  panel_index=$((panel_index + 1))
  tree_plot_panel_args+=(
    "--panel${panel_index}=ortholog,${ortholog_prefix},${file_og_dated_tree}"
  )

  TREEVIS_SPECIES_PARSER="${species_label_parser}" \
  Rscript "${gg_support_dir}/stat_branch2tree_plot.r" \
    --stat_branch="${file_og_stat_branch}" \
    --max_delta_intron_present="${treevis_retrotransposition_delta_intron}" \
    --width="7.2" \
    --rel_widths="" \
    "${tree_plot_panel_args[@]}" \
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
  gg_artifact_record \
    "${tree_plot_provenance_args[@]}" \
    --diagnostic "genegalleon_version=${GG_VERSION:-${SINGULARITYENV_GG_VERSION:-${APPTAINERENV_GG_VERSION:-unknown}}}" \
    --diagnostic "container_image=${gg_container_image_path:-unknown}"
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
  (
    if ! gg_shared_lock_acquire "${lock_file}" "parameter artifact copy (${file_to})"; then
      exit 1
    fi
    gg_shared_lock_start_heartbeat "${lock_file}"
    heartbeat_pid=${GG_SHARED_LOCK_HEARTBEAT_PID:-}
    cleanup_parameter_copy_lock() {
      gg_shared_lock_stop_heartbeat "${heartbeat_pid}"
      gg_shared_lock_release "${lock_file}"
    }
    trap cleanup_parameter_copy_lock EXIT
    filesize_from=$(stat -c%s "${file_from}")
    filesize_to=0
    if [[ -s "${file_to}" ]]; then
      filesize_to=$(stat -c%s "${file_to}")
    fi
    if [[ ${filesize_from} -ne ${filesize_to} ]]; then
      echo "Storing important files for record: ${file_to}"
      cp_out "${file_from}" "${file_to}"
    fi
  ) || exit 1
done

cd "${gg_workspace_dir}" || exit 1
remove_empty_subdirs "${dir_output_active}"

gene_family_outputs_complete=0
if [[ -s "${file_og_stat_branch}" && -s "${file_og_stat_tree}" && -s "${file_og_tree_plot}" ]]; then
  gene_family_outputs_complete=1
fi
if [[ "${gene_family_output_storage}" == "zip" && ${gg_debug_mode:-0} -eq 0 ]]; then
  if ! finalize_gene_family_run_success; then
    echo "Failed to record completed gene-family state for ${og_id}." >&2
    exit 1
  fi
fi
if [[ ${gene_family_run_succeeded:-0} -eq 1 && -n "${gene_family_materialization_receipt:-}" ]]; then
  rm -f -- "${gene_family_materialization_receipt}"
fi
if [[ "${gene_family_output_storage}" == "zip" && ${gg_debug_mode:-0} -eq 0 ]]; then
  if ! python "${gene_family_store_script}" archive-completed \
    "${gene_family_store_context_args[@]}" \
    "${gene_family_archive_write_args[@]}" \
    --min-files "${gene_family_zip_min_batch_files}" \
    --nonblocking
  then
    echo "Warning: Failed to archive completed gene-family outputs; live files were preserved." >&2
  fi
fi

if [[ ${gene_family_outputs_complete} -eq 1 && ${gg_debug_mode:-0} -eq 0 ]]; then
  echo "Output files detected."
  echo "${file_og_stat_branch}"
  echo "${file_og_stat_tree}"
  echo "${file_og_tree_plot}"
  echo "$(date): Exiting Singularity environment"
  exit 8
elif [[ ${gene_family_outputs_complete} -eq 1 && ${gg_debug_mode:-0} -eq 1 ]]; then
  echo "Output files detected & debug mode."
else
  echo "Output files not found."
fi

###################
echo "$(date): Exiting Singularity environment"
