# Named stage functions for gg_genome_evolution_core.sh.
# This file is sourced by workflow/core/gg_genome_evolution_core.sh.

resolve_busco_lineage_for_species_set() {
  if [[ -n "${busco_lineage_resolved}" ]]; then
    return 0
  fi
  if [[ $# -eq 0 ]]; then
    echo "No species names were provided for BUSCO lineage resolution." >&2
    return 1
  fi
  if ! busco_lineage_resolved=$(gg_resolve_busco_lineage "${gg_workspace_dir}" "${busco_lineage}" "$@"); then
    echo "Failed to resolve BUSCO lineage from request: ${busco_lineage}" >&2
    return 1
  fi
  echo "Resolved BUSCO lineage for species set (${#} species): ${busco_lineage_resolved}"
}

trim_ascii_whitespace() {
  local s=$1
  s="${s#"${s%%[![:space:]]*}"}"
  s="${s%"${s##*[![:space:]]}"}"
  printf '%s' "${s}"
}

species_tree_rooting_looks_like_label_list() {
  local token=""
  if [[ $# -eq 0 ]]; then
    return 1
  fi
  for token in "$@"; do
    if [[ -z "${token}" ]]; then
      continue
    fi
    if [[ ! "${token}" =~ ^[[:alnum:].-]+_[[:alnum:]_.-]+$ ]]; then
      return 1
    fi
  done
  return 0
}

species_genetic_code_table_path() {
  echo "${gg_workspace_input_dir}/species_genetic_code/species_genetic_code.tsv"
}

species_protein_input_dir_path() {
  echo "${gg_workspace_input_dir}/species_protein"
}

prepare_species_genetic_code_table() {
  local cds_dir=$1
  local default_code=$2
  local outfile=$3
  local input_table
  input_table=$(species_genetic_code_table_path)
  python - "${cds_dir}" "${default_code}" "${outfile}" "${input_table}" <<'PY'
import csv
import gzip
import os
import sys
from pathlib import Path


def species_from_filename(path: Path) -> str:
    name = path.name
    suffixes = (
        ".derived.cds.fa.gz",
        "_derived.cds.fa.gz",
        ".cds.all.fa.gz",
        "_cds.all.fa.gz",
        ".cds.fa.gz",
        "_cds.fa.gz",
        ".cds.fna.gz",
        "_cds.fna.gz",
        ".genome.fa.gz",
        "_genome.fa.gz",
        ".genomic.fna.gz",
        "_genomic.fna.gz",
        ".pep.fa.gz",
        "_pep.fa.gz",
        ".protein.fa.gz",
        "_protein.fa.gz",
        ".fa.gz",
        ".fas.gz",
        ".fasta.gz",
        ".fna.gz",
        ".fa",
        ".fas",
        ".fasta",
        ".fna",
    )
    lowered = name.lower()
    for suffix in sorted(suffixes, key=len, reverse=True):
        if lowered.endswith(suffix):
            name = name[: -len(suffix)]
            break
    parts = [part for part in name.split("_") if part]
    if len(parts) < 2:
        return path.stem
    def token_key(token):
        lowered_token = token.rstrip(".").lower()
        return {
            "spp": "sp",
            "ssp": "subsp",
            "subspecies": "subsp",
            "variety": "var",
            "form": "forma",
        }.get(lowered_token, lowered_token)
    def prefix_part(token):
        key = token_key(token)
        if key in {
            "sp", "cf", "aff", "nr", "subsp", "var", "forma", "f",
            "strain", "substrain", "serovar", "serotype", "serogroup",
            "pathovar", "pv", "biovar", "biotype", "chemovar", "morphovar",
            "cultivar", "cv", "isolate", "group", "subgroup", "complex",
            "clade", "lineage", "section", "series", "ecotype", "breed",
        }:
            return token
        return token.split(".", 1)[0]
    second = token_key(parts[1])
    third = token_key(parts[2]) if len(parts) >= 3 else ""
    if second == "sp":
        count = 3 if len(parts) >= 3 else 2
    elif second in {"cf", "aff", "nr"}:
        count = 3 if len(parts) >= 3 else 2
    elif third in {"cf", "aff", "nr"}:
        count = 3
    elif third in {
        "subsp", "var", "forma", "f", "strain", "substrain", "serovar",
        "serotype", "serogroup", "pathovar", "pv", "biovar", "biotype",
        "chemovar", "morphovar", "cultivar", "cv", "isolate", "group",
        "subgroup", "complex", "clade", "lineage", "section", "series",
        "ecotype", "breed",
    }:
        count = 4 if len(parts) >= 4 else 3
    else:
        count = 2
    return "_".join(prefix_part(part) for part in parts[:count] if prefix_part(part))


cds_dir = Path(sys.argv[1])
default_code = sys.argv[2].strip()
outfile = Path(sys.argv[3])
input_table = Path(sys.argv[4])
fasta_suffixes = (".fa", ".fa.gz", ".fas", ".fas.gz", ".fasta", ".fasta.gz", ".fna", ".fna.gz")

species = sorted(
    species_from_filename(path)
    for path in cds_dir.iterdir()
    if path.is_file() and not path.name.startswith(".") and path.name.endswith(fasta_suffixes)
)
if not species:
    sys.stderr.write(f"No species CDS files were found in: {cds_dir}\n")
    sys.exit(1)

overrides = {}
if input_table.exists():
    with input_table.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(
            (line for line in handle if line.strip() and not line.lstrip().startswith("#")),
            delimiter="\t",
        )
        if reader.fieldnames is None or "species" not in reader.fieldnames or "genetic_code" not in reader.fieldnames:
            sys.stderr.write(
                f"{input_table} must be a tab-delimited file with at least species and genetic_code columns.\n"
            )
            sys.exit(1)
        for row in reader:
            sp = (row.get("species") or "").strip().replace(" ", "_")
            code = (row.get("genetic_code") or "").strip()
            if not sp or not code:
                sys.stderr.write(f"Empty species or genetic_code entry was detected in {input_table}.\n")
                sys.exit(1)
            try:
                code_int = int(code)
            except ValueError:
                sys.stderr.write(f"Invalid genetic_code for {sp}: {code}\n")
                sys.exit(1)
            if code_int <= 0:
                sys.stderr.write(f"genetic_code must be a positive integer for {sp}: {code}\n")
                sys.exit(1)
            if sp in overrides:
                sys.stderr.write(f"Duplicate species entry in {input_table}: {sp}\n")
                sys.exit(1)
            overrides[sp] = str(code_int)

outfile.parent.mkdir(parents=True, exist_ok=True)
with outfile.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
    writer.writerow(["species", "genetic_code", "source"])
    for sp in species:
        if sp in overrides:
            writer.writerow([sp, overrides[sp], "species_genetic_code.tsv"])
        else:
            writer.writerow([sp, default_code, "default"])

unknown = sorted(sp for sp in overrides if sp not in set(species))
if unknown:
    sys.stderr.write(
        "Warning: species_genetic_code.tsv entries were ignored because no matching species_cds file was found: "
        + ", ".join(unknown)
        + "\n"
    )
PY
}

lookup_species_genetic_code() {
  local species_name=$1
  local table_path=$2
  local code=""
  if [[ -s "${table_path}" ]]; then
    code=$(awk -F'\t' -v sp="${species_name}" 'NR>1 && $1==sp {print $2; exit}' "${table_path}")
  fi
  if [[ -z "${code}" ]]; then
    code="${genetic_code}"
  fi
  printf '%s' "${code}"
}

parse_species_tree_rooting() {
  local raw_config=$1
  local -n method_ref=$2
  local -n value_ref=$3
  local -a fields=()
  local -a normalized_fields=()
  local idx
  local method=""

  IFS=',' read -r -a fields <<< "${raw_config}"
  if [[ ${#fields[@]} -eq 0 ]]; then
    echo "species_tree_rooting is empty."
    echo 'Expected one of: "outgroup,GENUS_SPECIES[,GENUS_SPECIES...]", "midpoint", "mad", "mv", "taxonomy[,ncbi[,opentree,timetree...]]".'
    return 1
  fi

  for idx in "${!fields[@]}"; do
    normalized_fields[idx]=$(trim_ascii_whitespace "${fields[idx]}")
  done
  method=$(printf '%s' "${normalized_fields[0]}" | tr '[:upper:]' '[:lower:]')
  if [[ -z "${method}" ]]; then
    echo "species_tree_rooting is missing a rooting method: ${raw_config}"
    return 1
  fi

  case "${method}" in
    outgroup)
      if [[ ${#normalized_fields[@]} -lt 2 ]]; then
        echo "species_tree_rooting=${raw_config} is invalid."
        echo 'Outgroup rooting requires at least one species label after "outgroup,".'
        return 1
      fi
      value_ref=""
      for ((idx = 1; idx < ${#normalized_fields[@]}; idx++)); do
        if [[ -z "${normalized_fields[idx]}" ]]; then
          continue
        fi
        if [[ -n "${value_ref}" ]]; then
          value_ref+=","
        fi
        value_ref+="${normalized_fields[idx]}"
      done
      if [[ -z "${value_ref}" ]]; then
        echo "species_tree_rooting=${raw_config} is invalid."
        echo 'Outgroup rooting requires one or more non-empty species labels after "outgroup,".'
        return 1
      fi
      ;;
    midpoint | mad | mv)
      if [[ ${#normalized_fields[@]} -gt 1 ]]; then
        echo "species_tree_rooting=${raw_config} is invalid."
        echo "Method ${method} does not accept additional comma-separated fields."
        return 1
      fi
      value_ref=""
      ;;
    md)
      if [[ ${#normalized_fields[@]} -gt 1 ]]; then
        echo "species_tree_rooting=${raw_config} is invalid."
        echo 'Method md does not accept additional comma-separated fields.'
        return 1
      fi
      method="mv"
      value_ref=""
      ;;
    taxonomy)
      value_ref=""
      for ((idx = 1; idx < ${#normalized_fields[@]}; idx++)); do
        if [[ -z "${normalized_fields[idx]}" ]]; then
          continue
        fi
        if [[ -n "${value_ref}" ]]; then
          value_ref+=","
        fi
        value_ref+="${normalized_fields[idx]}"
      done
      ;;
    *)
      if species_tree_rooting_looks_like_label_list "${normalized_fields[@]}"; then
        method="outgroup"
        value_ref=""
        for ((idx = 0; idx < ${#normalized_fields[@]}; idx++)); do
          if [[ -z "${normalized_fields[idx]}" ]]; then
            continue
          fi
          if [[ -n "${value_ref}" ]]; then
            value_ref+=","
          fi
          value_ref+="${normalized_fields[idx]}"
        done
        echo "species_tree_rooting=${raw_config} uses legacy outgroup-label syntax; interpreting it as outgroup,${value_ref}."
      else
        echo "Invalid species_tree_rooting: ${raw_config}"
        echo 'species_tree_rooting must be one of "outgroup,GENUS_SPECIES[,GENUS_SPECIES...]", "midpoint", "mad", "mv", or "taxonomy[,ncbi[,opentree,timetree...]]".'
        return 1
      fi
      ;;
  esac

  method_ref="${method}"
  return 0
}

parse_mcmctree_constraint_record() {
  local raw_record=$1
  local __out_name=$2
  local -n __out_ref="${__out_name}"
  local idx

  IFS=',' read -r -a __out_ref <<< "${raw_record}"
  if [[ ${#__out_ref[@]} -ne 4 ]]; then
    echo "Invalid mcmctree divergence-time constraint: ${raw_record}"
    echo "Expected exactly 4 comma-separated fields: left_species,right_species,lower_bound,upper_bound"
    return 1
  fi

  for idx in 0 1 2 3; do
    __out_ref[${idx}]=$(trim_ascii_whitespace "${__out_ref[${idx}]}")
  done

  if [[ -z "${__out_ref[0]}" || -z "${__out_ref[1]}" || -z "${__out_ref[2]}" || -z "${__out_ref[3]}" ]]; then
    echo "Invalid mcmctree divergence-time constraint (empty field): ${raw_record}"
    echo "Expected non-empty fields: left_species,right_species,lower_bound,upper_bound"
    return 1
  fi

  return 0
}

clear_directory_contents_safe() {
  local target_dir=$1
  local target_real
  local entries=()

  if [[ -z "${target_dir}" || ! -d "${target_dir}" ]]; then
    echo "Refusing to clear missing directory: ${target_dir}" >&2
    return 1
  fi

  target_real=$(cd "${target_dir}" && pwd -P)
  if [[ -z "${target_real}" || "${target_real}" == "/" ]]; then
    echo "Refusing to clear unsafe directory: ${target_dir}" >&2
    return 1
  fi

  shopt -s dotglob nullglob
  entries=("${target_real}"/*)
  shopt -u dotglob nullglob
  if [[ ${#entries[@]} -gt 0 ]]; then
    rm -rf -- "${entries[@]}"
  fi
}

sync_genome_busco_summary_table_from_shared() {
  if [[ "${file_genome_busco_summary_table}" == "${file_species_busco_summary_table}" ]]; then
    return 0
  fi
  if [[ ! -s "${file_species_busco_summary_table}" ]]; then
    return 1
  fi
  ensure_parent_dir "${file_genome_busco_summary_table}"
  if [[ -s "${file_genome_busco_summary_table}" ]] && cmp -s "${file_species_busco_summary_table}" "${file_genome_busco_summary_table}"; then
    return 0
  fi
  cp_out "${file_species_busco_summary_table}" "${file_genome_busco_summary_table}"
}

run_shared_species_busco_stage() {
  local task="BUSCO analysis of species-wise input files"
  local source_species_input_dir=""
  local -a source_species_input_fasta=()
  local -a species_input_fasta=()
  local -a input_species_set=()
  local -a busco_output_files=()
  local seq_full seq_file sp_ub file_sp_busco_full file_sp_busco_short
  local input_species busco_file busco_base busco_species busco_species_found
  local dir_busco_db="" dir_busco_lineage=""
  local full_found=0 short_found=0
  local missing_busco_outputs=0

  if [[ ${shared_species_busco_stage_done} -eq 1 ]]; then
    return 0
  fi
  shared_species_busco_stage_done=1

  if [[ ${run_species_busco} -ne 1 ]]; then
    gg_step_skip "${task}"
    return 0
  fi

  source_species_input_dir=$(effective_species_input_source_dir_path)
  ensure_dir "${dir_species_busco_full}"
  ensure_dir "${dir_species_busco_short}"
  mapfile -t source_species_input_fasta < <(gg_find_fasta_files "${source_species_input_dir}" 1)
  echo "Number of ${input_sequence_mode} files for BUSCO: ${#source_species_input_fasta[@]}"
  if [[ ${#source_species_input_fasta[@]} -eq 0 ]]; then
    echo "No ${input_sequence_mode} file found. Exiting."
    exit 1
  fi
  mapfile -t input_species_set < <(gg_species_names_from_fasta_dir "${source_species_input_dir}")
  mapfile -t busco_output_files < <(
    find "${dir_species_busco_full}" "${dir_species_busco_short}" -maxdepth 1 -type f \
      \( -name "*busco.full.tsv" -o -name "*busco.short.txt" \) \
      2> /dev/null | sort
  )
  for busco_file in "${busco_output_files[@]}"; do
    busco_base=$(basename "${busco_file}")
    busco_species=$(gg_species_name_from_path_or_dot "${busco_base}")
    busco_species_found=0
    for input_species in "${input_species_set[@]}"; do
      if [[ "${input_species}" == "${busco_species}" ]]; then
        busco_species_found=1
        break
      fi
    done
    if [[ ${busco_species_found} -eq 0 ]]; then
      echo "Removing stale BUSCO output for species not in current input: ${busco_file}"
      rm -f -- "${busco_file}"
    fi
  done
  for input_species in "${input_species_set[@]}"; do
    if busco_output_exists_for_species "${dir_species_busco_full}" "${input_species}" "*busco.full.tsv"; then
      full_found=1
    else
      full_found=0
    fi
    if busco_output_exists_for_species "${dir_species_busco_short}" "${input_species}" "*busco.short.txt"; then
      short_found=1
    else
      short_found=0
    fi
    if [[ ${full_found} -ne 1 || ${short_found} -ne 1 ]]; then
      missing_busco_outputs=1
      break
    fi
  done
  if [[ ${missing_busco_outputs} -ne 1 ]]; then
    for seq_full in "${source_species_input_fasta[@]}"; do
      echo "Skipped BUSCO: $(basename "${seq_full}")"
    done
    echo "$(date): End: ${task}"
    return 0
  fi
  if ! resolve_busco_lineage_for_species_set "${input_species_set[@]}"; then
    exit 1
  fi

  prepare_species_tree_input_dir
  mapfile -t species_input_fasta < <(gg_find_fasta_files "${species_tree_input_dir}" 1)
  for seq_full in "${species_input_fasta[@]}"; do
    seq_file=$(basename "${seq_full}")
    sp_ub=$(gg_species_name_from_path_or_dot "${seq_file}")
    file_sp_busco_full="${dir_species_busco_full}/${sp_ub}.busco.full.tsv"
    file_sp_busco_short="${dir_species_busco_short}/${sp_ub}.busco.short.txt"
    if busco_output_exists_for_species "${dir_species_busco_full}" "${sp_ub}" "*busco.full.tsv"; then
      full_found=1
    else
      full_found=0
    fi
    if busco_output_exists_for_species "${dir_species_busco_short}" "${sp_ub}" "*busco.short.txt"; then
      short_found=1
    else
      short_found=0
    fi
    if [[ ${full_found} -ne 1 || ${short_found} -ne 1 ]]; then
      gg_step_start "${task}: ${seq_file}"
      remove_busco_outputs_for_species "${dir_species_busco_full}" "${sp_ub}" "*busco.full.tsv"
      remove_busco_outputs_for_species "${dir_species_busco_short}" "${sp_ub}" "*busco.short.txt"
      seqkit seq --threads "${GG_TASK_CPUS}" "${species_tree_input_dir}/${seq_file}" --out-file "tmp.busco_input.fasta"

      if ! dir_busco_db=$(ensure_busco_download_path "${gg_workspace_dir}" "${busco_lineage_resolved}"); then
        echo "Failed to prepare BUSCO dataset: ${busco_lineage_resolved}"
        exit 1
      fi
      dir_busco_lineage="${dir_busco_db}/lineages/${busco_lineage_resolved}"

      gg_run_busco_with_metaeuk_modified_fas_compat \
        --in "tmp.busco_input.fasta" \
        --mode "${species_tree_busco_mode}" \
        --out "busco_tmp" \
        --cpu "${GG_TASK_CPUS}" \
        --force \
        --evalue 1e-03 \
        --limit 20 \
        --lineage_dataset "${dir_busco_lineage}" \
        --download_path "${dir_busco_db}" \
        --offline

      if copy_busco_tables "./busco_tmp" "${busco_lineage_resolved}" "${file_sp_busco_full}" "${file_sp_busco_short}"; then
        rm -rf -- "./busco_tmp"
      else
        echo "Failed to locate normalized BUSCO outputs for ${sp_ub}. Exiting."
        exit 1
      fi
    else
      echo "Skipped BUSCO: ${seq_file}"
    fi
  done
  echo "$(date): End: ${task}"
}

run_shared_busco_summary_stage() {
  local task="Collecting IDs of common BUSCO genes"
  local num_busco_ids=0
  local source_species_input_dir=""

  if [[ ${shared_busco_summary_stage_done} -eq 1 ]]; then
    sync_genome_busco_summary_table_from_shared || true
    return 0
  fi
  shared_busco_summary_stage_done=1

  if [[ ${run_build_species_busco_summary} -ne 1 ]]; then
    gg_step_skip "${task}"
    return 0
  fi

  source_species_input_dir=$(effective_species_input_source_dir_path)
  normalize_busco_table_naming "${dir_species_busco_full}" "${dir_species_busco_short}"
  if ! is_species_set_identical "${source_species_input_dir}" "${dir_species_busco_full}"; then
    echo "Exiting due to species-set mismatch between ${source_species_input_dir} and ${dir_species_busco_full}"
    exit 1
  fi
  if [[ ! -s "${file_species_busco_summary_table}" ]]; then
    gg_step_start "${task}"
    ensure_parent_dir "${file_species_busco_summary_table}"

    python "${gg_support_dir}/collect_common_BUSCO_genes.py" \
      --busco_outdir "${dir_species_busco_full}" \
      --ncpu "${GG_TASK_CPUS}" \
      --outfile "tmp.busco_summary_table.tsv"
    mv_out "tmp.busco_summary_table.tsv" "${file_species_busco_summary_table}"

    num_busco_ids=$(get_busco_summary_gene_count "${file_species_busco_summary_table}")
    echo "Number of BUSCO genes: ${num_busco_ids}"
    echo "$(date): End: ${task}"
  else
    gg_step_skip "${task}"
  fi
  sync_genome_busco_summary_table_from_shared || true
}

run_shared_species_omark_stage() {
  local task="OMArk analysis of species-wise protein input files"
  local source_species_input_dir=""
  local -a source_species_input_fasta=()
  local -a species_input_fasta=()
  local -a input_species_set=()
  local -a existing_species_dirs=()
  local protein_full protein_file sp_ub omamer_out omark_outdir sum_file omamer_query
  local input_species existing_dir existing_species existing_found
  local omark_db_file=""
  local missing_omark_outputs=0

  if [[ ${shared_species_omark_stage_done} -eq 1 ]]; then
    return 0
  fi
  shared_species_omark_stage_done=1

  if [[ ${run_species_omark} -ne 1 ]]; then
    gg_step_skip "${task}"
    return 0
  fi

  ensure_dir "${dir_species_omamer}"
  ensure_dir "${dir_species_omark}"
  source_species_input_dir=$(effective_species_input_source_dir_path)
  mapfile -t source_species_input_fasta < <(gg_find_fasta_files "${source_species_input_dir}" 1)
  echo "Number of protein files for OMArk: ${#source_species_input_fasta[@]}"
  if [[ ${#source_species_input_fasta[@]} -eq 0 ]]; then
    echo "No protein file found for OMArk. Exiting."
    exit 1
  fi

  mapfile -t input_species_set < <(gg_species_names_from_fasta_dir "${source_species_input_dir}")

  mapfile -t existing_species_dirs < <(find "${dir_species_omark}" -maxdepth 1 -mindepth 1 -type d ! -name '.*' | sort)
  for existing_dir in "${existing_species_dirs[@]}"; do
    existing_species=$(basename "${existing_dir}")
    existing_found=0
    for input_species in "${input_species_set[@]}"; do
      if [[ "${input_species}" == "${existing_species}" ]]; then
        existing_found=1
        break
      fi
    done
    if [[ ${existing_found} -eq 0 ]]; then
      echo "Removing stale OMArk output for species not in current input: ${existing_dir}"
      rm -rf -- "${existing_dir}"
      rm -f -- "${dir_species_omamer}/${existing_species}.omamer"
    fi
  done
  for input_species in "${input_species_set[@]}"; do
    omamer_out="${dir_species_omamer}/${input_species}.omamer"
    sum_file="${dir_species_omark}/${input_species}/${input_species}.sum"
    if [[ ! -s "${omamer_out}" || ! -s "${sum_file}" ]]; then
      missing_omark_outputs=1
      break
    fi
  done
  if [[ ${missing_omark_outputs} -ne 1 ]]; then
    for protein_full in "${source_species_input_fasta[@]}"; do
      echo "Skipped OMArk: $(basename "${protein_full}")"
    done
    echo "$(date): End: ${task}"
    return 0
  fi

  if ! omark_db_file=$(ensure_omark_database "${gg_workspace_dir}" "${omark_db_path}"); then
    echo "Failed to prepare OMArk database: ${omark_db_path}"
    exit 1
  fi
  omark_db_resolved="${omark_db_file}"
  prepare_species_protein_tmp
  mapfile -t species_input_fasta < <(gg_find_fasta_files "${dir_sp_protein}" 1)

  for protein_full in "${species_input_fasta[@]}"; do
    protein_file=$(basename "${protein_full}")
    sp_ub=$(gg_species_name_from_path_or_dot "${protein_file}")
    omamer_out="${dir_species_omamer}/${sp_ub}.omamer"
    omark_outdir="${dir_species_omark}/${sp_ub}"
    sum_file="${omark_outdir}/${sp_ub}.sum"

    if [[ -s "${omamer_out}" && -s "${sum_file}" ]]; then
      echo "Skipped OMArk: ${protein_file}"
      continue
    fi

    gg_step_start "${task}: ${protein_file}"
    ensure_dir "${omark_outdir}"
    if [[ ! -s "${omamer_out}" ]]; then
      omamer_query="${omark_outdir}/${sp_ub}.query.fa"
      stage_species_protein_fasta "${protein_full}" "${omamer_query}"
      omamer search \
        --db "${omark_db_file}" \
        --query "${omamer_query}" \
        --out "${omamer_out}"
    fi
    if [[ ! -s "${omamer_out}" ]]; then
      echo "OMAmer search output was not created: ${omamer_out}"
      exit 1
    fi

    if [[ ! -s "${sum_file}" ]]; then
      omark \
        -f "${omamer_out}" \
        -d "${omark_db_file}" \
        -o "${omark_outdir}"
    fi
    if [[ ! -s "${sum_file}" ]]; then
      echo "OMArk summary output was not created: ${sum_file}"
      exit 1
    fi
  done
  echo "$(date): End: ${task}"
}

run_shared_omark_summary_stage() {
  local task="Summarizing OMArk species quality results"

  if [[ ${shared_omark_summary_stage_done} -eq 1 ]]; then
    return 0
  fi
  shared_omark_summary_stage_done=1

  if [[ ${run_build_species_omark_summary} -ne 1 ]]; then
    gg_step_skip "${task}"
    return 0
  fi

  if [[ ! -s "${file_species_omark_summary_table}" ]]; then
    gg_step_start "${task}"
    ensure_parent_dir "${file_species_omark_summary_table}"

    python "${gg_support_dir}/summarize_omark.py" \
      --omark_outdir "${dir_species_omark}" \
      --outfile "tmp.omark_summary.tsv"
    mv_out "tmp.omark_summary.tsv" "${file_species_omark_summary_table}"
    echo "$(date): End: ${task}"
  else
    gg_step_skip "${task}"
  fi
}

species_protein_input_has_files() {
  local protein_files=()
  mapfile -t protein_files < <(gg_find_fasta_files "${dir_sp_protein_input}" 1)
  [[ ${#protein_files[@]} -gt 0 ]]
}

species_cds_input_has_files() {
  local cds_files=()
  mapfile -t cds_files < <(gg_find_fasta_files "${dir_sp_cds}" 1)
  [[ ${#cds_files[@]} -gt 0 ]]
}

effective_species_input_source_dir_path() {
  if [[ "${input_sequence_mode}" == "protein" ]] && species_protein_input_has_files; then
    echo "${dir_sp_protein_input}"
  else
    echo "${dir_sp_cds}"
  fi
}

compute_shared_protein_input_signature() {
  local input_file
  local -a input_files=()
  local -a stat_lines=()
  local metadata_source="species_cds"
  if [[ -s "${file_species_genetic_code}" ]]; then
    input_files+=( "${file_species_genetic_code}" )
  fi
  if [[ "${input_sequence_mode}" == "protein" ]] && species_protein_input_has_files; then
    metadata_source="species_protein"
    while IFS= read -r input_file; do
      input_files+=( "${input_file}" )
    done < <(gg_find_fasta_files "${dir_sp_protein_input}" 1)
  else
    while IFS= read -r input_file; do
      input_files+=( "${input_file}" )
    done < <(gg_find_fasta_files "${dir_sp_cds}" 1)
  fi

  if [[ ${#input_files[@]} -gt 0 ]]; then
    if stat --version > /dev/null 2>&1; then
      mapfile -t stat_lines < <(stat -c '%n:%s:%Y' "${input_files[@]}")
    else
      mapfile -t stat_lines < <(stat -f '%N:%z:%m' "${input_files[@]}")
    fi
  fi

  {
    printf 'input_sequence_mode=%s\n' "${input_sequence_mode}"
    printf 'genetic_code=%s\n' "${genetic_code}"
    printf 'metadata_source=%s\n' "${metadata_source}"
    printf 'species_tree_busco_mode=%s\n' "${species_tree_busco_mode}"
    printf '%s\n' "${stat_lines[@]}"
  } | cksum | awk '{print $1}'
}

print_effective_genome_evolution_config_summary() {
  local species_protein_input_available=0
  local species_cds_input_available=0
  local species_genetic_code_table_present=0

  if species_protein_input_has_files; then
    species_protein_input_available=1
  fi
  if species_cds_input_has_files; then
    species_cds_input_available=1
  fi
  if [[ -s "${file_species_genetic_code}" ]]; then
    species_genetic_code_table_present=1
  fi

  gg_print_registered_config_summary \
    "gg_genome_evolution_entrypoint.sh" \
    "effective config summary (gg_genome_evolution_core.sh)" \
    genetic_code \
    busco_lineage \
    annotation_species \
    annotation_species_resolved \
    omark_db_path \
    omark_db_resolved \
    species_tree_rooting_method \
    species_tree_rooting_value \
    species_tree_busco_mode \
    species_tree_sequence_label \
    species_protein_input_available \
    species_cds_input_available \
    species_genetic_code_table_present
}

refresh_dir_for_shared_protein_input_signature() {
  local target_dir=$1
  local description=$2
  local signature=$3
  local stamp_file="${target_dir}/.shared_protein_input_signature"
  local previous_signature=""

  ensure_dir "${target_dir}"
  if [[ -s "${stamp_file}" ]]; then
    previous_signature=$(< "${stamp_file}")
  fi
  if [[ -n "${previous_signature}" && "${previous_signature}" != "${signature}" ]]; then
    echo "Shared protein input signature changed for ${description}. Clearing derived outputs in ${target_dir}"
    if ! clear_directory_contents_safe "${target_dir}"; then
      echo "Failed to clear ${description} directory after input signature change: ${target_dir}"
      exit 1
    fi
  fi
  ensure_dir "${target_dir}"
  printf '%s\n' "${signature}" > "${stamp_file}"
}

species_tree_summary_generation_requested() {
  [[ ${run_concat_iqtree_protein} -eq 1 ||
    ${run_concat_iqtree_dna} -eq 1 ||
    ${run_astral_pep} -eq 1 ||
    ${run_astral_dna} -eq 1 ||
    ${run_convert_tree_format} -eq 1 ]]
}

refresh_species_tree_for_shared_protein_input_signature() {
  local signature=$1
  local stamp_file="${dir_species_tree}/.shared_protein_input_signature"
  local previous_signature=""

  ensure_dir "${dir_species_tree}"
  if [[ -s "${stamp_file}" ]]; then
    previous_signature=$(< "${stamp_file}")
  fi
  if [[ -n "${previous_signature}" && "${previous_signature}" != "${signature}" ]]; then
    if [[ ${species_tree_requested_for_orthofinder} -eq 0 ]]; then
      echo "Shared protein input signature changed for species_tree, but species-tree generation flags are disabled."
      echo "Keeping existing species_tree outputs for reuse: ${dir_species_tree}"
      echo "The species_tree signature stamp will be updated next time species-tree generation is enabled."
      return 0
    fi
    echo "Shared protein input signature changed for species_tree. Clearing derived outputs in ${dir_species_tree}"
    if ! clear_directory_contents_safe "${dir_species_tree}"; then
      echo "Failed to clear species_tree directory after input signature change: ${dir_species_tree}"
      exit 1
    fi
  fi
  ensure_dir "${dir_species_tree}"
  printf '%s\n' "${signature}" > "${stamp_file}"
}

cleanup_species_protein_tmp() {
  local cleanup_target
  for cleanup_target in "${dir_sp_protein}" "${dir_sp_protein}_orthofinder" "${dir_sp_protein}_core" "${dir_sp_protein}_additional"; do
    if [[ -d "${cleanup_target}" ]]; then
      echo "Removing temporary species_protein directory: ${cleanup_target}"
      rm -rf -- "${cleanup_target}"
    fi
  done
  species_protein_ready=0
}

stage_species_protein_fasta() {
  local source_fasta=$1
  local staged_fasta=$2
  ensure_parent_dir "${staged_fasta}"
  seqkit seq --threads "${GG_TASK_CPUS}" "${source_fasta}" --out-file "${staged_fasta}"
  if [[ ! -s "${staged_fasta}" ]]; then
    echo "Temporary species_protein FASTA was not created: ${staged_fasta}"
    exit 1
  fi
}

prepare_species_protein_orthofinder_dir() {
  local source_dir=$1
  local target_dir=$2
  local source_fasta source_name source_stem staged_fasta
  if [[ -e "${target_dir}" ]]; then
    rm -rf -- "${target_dir}"
  fi
  ensure_dir "${target_dir}"
  while IFS= read -r source_fasta; do
    source_name=$(basename "${source_fasta}")
    source_stem=${source_name%.gz}
    source_stem=${source_stem%.*}
    staged_fasta="${target_dir}/${source_stem}.fa"
    stage_species_protein_fasta "${source_fasta}" "${staged_fasta}"
  done < <(gg_find_fasta_files "${source_dir}" 1)
}

prepare_species_protein_tmp() {
  local cds_path cds sp_ub translated_file species_code protein_path
  local -a cds_files=()
  local -a protein_files=()

  if [[ ${species_protein_ready} -eq 1 ]]; then
    return 0
  fi

  gg_step_start "Prepare temporary species_protein FASTA files"
  cleanup_species_protein_tmp
  ensure_dir "${dir_sp_protein}"

  if [[ "${input_sequence_mode}" == "protein" ]] && species_protein_input_has_files; then
    check_species_protein_dir "${dir_sp_protein_input}"
    check_if_species_files_unique "${dir_sp_protein_input}"
    mapfile -t protein_files < <(gg_find_fasta_files "${dir_sp_protein_input}" 1)
    for protein_path in "${protein_files[@]}"; do
      sp_ub=$(gg_species_name_from_path "$(basename "${protein_path}")")
      translated_file="${sp_ub}.fa.gz"
      echo "Copying protein FASTA: $(basename "${protein_path}") -> ${translated_file}"
      seqkit seq --threads "${GG_TASK_CPUS}" "${protein_path}" --out-file "${dir_sp_protein}/${translated_file}"
    done
    if [[ -s "${file_species_genetic_code}" ]]; then
      echo "species_genetic_code.tsv is ignored because species_protein inputs are provided: ${file_species_genetic_code}"
    fi
    species_protein_source="species_protein"
    species_protein_ready=1
    return 0
  fi
  if [[ "${input_sequence_mode}" != "protein" ]] && species_protein_input_has_files; then
    echo "Ignoring species_protein inputs in cds mode: ${dir_sp_protein_input}"
    echo "cds mode always generates temporary species_protein FASTA files from species_cds."
  fi

  if ! species_cds_input_has_files; then
    echo "protein mode requires either species_protein inputs or species_cds inputs for fallback translation."
    echo "Checked: ${dir_sp_protein_input} and ${dir_sp_cds}"
    exit 1
  fi
  if [[ ${run_cds_translation} -ne 1 ]]; then
    echo "run_cds_translation must be 1 when species proteins need to be generated from species_cds."
    exit 1
  fi

  check_species_cds "${gg_workspace_dir}"
  check_if_species_files_unique "${dir_sp_cds}"
  prepare_species_genetic_code_table "${dir_sp_cds}" "${genetic_code}" "${file_species_genetic_code_resolved}"
  mapfile -t cds_files < <(gg_find_fasta_files "${dir_sp_cds}" 1)
  for cds_path in "${cds_files[@]}"; do
    cds=$(basename "${cds_path}")
    if zgrep -q -e "^>.*[[:blank:]]" "${cds_path}"; then
      echo "Space (\" \") is detected. Please remove all annotation info after spaces in sequence names. Exiting: ${cds}"
      exit 1
    fi
    sp_ub=$(gg_species_name_from_path "${cds}")
    translated_file="${sp_ub}.fa.gz"
    species_code=$(lookup_species_genetic_code "${sp_ub}" "${file_species_genetic_code_resolved}")
    echo "Translation started: ${cds} (genetic_code=${species_code}) -> ${translated_file}"

    seqkit seq --remove-gaps --threads "${GG_TASK_CPUS}" "${cds_path}" |
      gg_prepare_cds_fasta_stream "${GG_TASK_CPUS}" "${species_code}" |
      seqkit translate --transl-table "${species_code}" --threads "${GG_TASK_CPUS}" |
      sed -e "s/^>${sp_ub}[-_\.]/>/" -e "s/^>/>${sp_ub}_/" |
      sed -e '/^1 1$/d' -e 's/_frame=1[[:space:]]*//' |
      seqkit seq --threads "${GG_TASK_CPUS}" --out-file "${dir_sp_protein}/${translated_file}"
  done
  species_protein_source="species_cds"
  species_protein_ready=1
}

prepare_species_tree_input_dir() {
  if [[ "${input_sequence_mode}" == "protein" ]]; then
    prepare_species_protein_tmp
    species_tree_input_dir="${dir_sp_protein}"
  else
    species_tree_input_dir="${dir_sp_cds}"
  fi
}

root_species_tree() {
  local infile=$1
  local outfile=$2
  local tree_description=$3
  local root_log="${dir_tmp}/tmp.nwkit.root.$$.log"
  local outgroup_label_list=()
  local root_method="${species_tree_rooting_method}"
  local root_value="${species_tree_rooting_value}"
  local missing_outgroup=0
  local root_exit_code=0
  local -a nwkit_root_args=()
  ensure_parent_dir "${outfile}"
  rm -f -- "${outfile}" "${root_log}"

  nwkit_root_args=(--method "${root_method}" --infile "${infile}" --outfile "${outfile}")
  if [[ "${root_method}" == "outgroup" ]]; then
    mapfile -t outgroup_label_list < <(printf '%s' "${root_value}" | tr ',' '\n' | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' -e '/^$/d')
    for outgroup_label in "${outgroup_label_list[@]}"; do
      if ! grep -q -F -- "${outgroup_label}" "${infile}"; then
        missing_outgroup=1
        break
      fi
    done
    if [[ ${missing_outgroup} -eq 1 ]]; then
      echo "Error: Outgroup labels (${root_value}) are not present in ${tree_description}."
      return 1
    fi
    nwkit_root_args+=(--outgroup "${root_value}")
  elif [[ "${root_method}" == "taxonomy" ]]; then
    ensure_dir "${dir_nwkit_download_dir}"
    nwkit_root_args+=(--download_dir "${dir_nwkit_download_dir}")
    nwkit_root_args+=(--species-parser "${species_label_parser}")
    if [[ -n "${species_label_regex}" ]]; then
      nwkit_root_args+=(--species-regex "${species_label_regex}")
    fi
    if [[ -n "${species_label_map_tsv}" ]]; then
      nwkit_root_args+=(--species-map-tsv "${species_label_map_tsv}")
    fi
    if [[ -n "${root_value}" ]]; then
      nwkit_root_args+=(--taxonomy_source "${root_value}")
    fi
  fi

  if nwkit root "${nwkit_root_args[@]}" 2> "${root_log}"; then
    root_exit_code=0
  else
    root_exit_code=$?
  fi

  if [[ ${root_exit_code} -eq 0 && -s "${outfile}" ]]; then
    rm -f -- "${root_log}"
    return 0
  fi

  echo "Error: Failed to root ${tree_description}."
  [[ -s "${root_log}" ]] && cat "${root_log}"
  rm -f -- "${root_log}"
  return 1
}

restore_rooted_tree_internal_support() {
  local rooted_tree=$1
  local support_source_tree=$2
  local tree_description=$3
  local tmp_with_support="${dir_tmp}/tmp.nwkit.transfer_support.$$.nwk"
  local transfer_log="${dir_tmp}/tmp.nwkit.transfer_support.$$.log"
  local transfer_exit_code=0

  if [[ ! -s "${rooted_tree}" ]]; then
    echo "Warning: Cannot restore ${tree_description} support values. Missing rooted tree: ${rooted_tree}"
    return 1
  fi
  if [[ ! -s "${support_source_tree}" ]]; then
    echo "Warning: Cannot restore ${tree_description} support values. Missing support source tree: ${support_source_tree}"
    return 1
  fi

  rm -f -- "${tmp_with_support}" "${transfer_log}"
  if nwkit transfer \
    --infile "${rooted_tree}" \
    --infile2 "${support_source_tree}" \
    --target intnode \
    --support yes \
    --name no \
    --length no \
    --outfile "${tmp_with_support}" \
    2> "${transfer_log}"; then
    transfer_exit_code=0
  else
    transfer_exit_code=$?
  fi

  if [[ ${transfer_exit_code} -eq 0 && -s "${tmp_with_support}" ]]; then
    echo "Restored internal support values for ${tree_description} with nwkit transfer."
    cp_out "${tmp_with_support}" "${rooted_tree}"
    rm -f -- "${tmp_with_support}" "${transfer_log}"
    return 0
  fi

  echo "Warning: Failed to restore internal support values for ${tree_description}. Keeping rooted tree as-is."
  [[ -s "${transfer_log}" ]] && cat "${transfer_log}"
  rm -f -- "${tmp_with_support}" "${transfer_log}"
  return 1
}

save_astral_label_tree() {
  local infile=$1
  local outfile=$2
  local tree_description=$3

  if [[ "${species_tree_rooting_method}" == "mad" || "${species_tree_rooting_method}" == "mv" ]]; then
    echo "Saving unrooted ${tree_description}; ${species_tree_rooting_method} rooting will be applied after IQ-TREE branch-length optimization."
    cp_out "${infile}" "${outfile}"
  else
    root_species_tree "${infile}" "${outfile}" "${tree_description}"
  fi
}

astral_label_tree_rooting_is_deferred() {
  [[ "${species_tree_rooting_method}" == "mad" || "${species_tree_rooting_method}" == "mv" ]]
}

transfer_deferred_astral_label_tree_roots() {
  local rooted_reference_tree=$1
  local tree_dir=$2
  local file_prefix=$3
  local tree_description=$4
  local label_tree=""
  local tmp_transfer_tree=""
  local -a label_trees=()

  if ! astral_label_tree_rooting_is_deferred; then
    return 0
  fi
  if [[ ! -s "${rooted_reference_tree}" ]]; then
    echo "Warning: Cannot transfer ${species_tree_rooting_method} root to ${tree_description}. Missing rooted reference tree: ${rooted_reference_tree}"
    return 1
  fi

  shopt -s nullglob
  label_trees=("${tree_dir}/${file_prefix}".*.nwk)
  shopt -u nullglob
  for label_tree in "${label_trees[@]}"; do
    if [[ "$(basename "${label_tree}")" == "${file_prefix}.optimized.nwk" ]]; then
      continue
    fi
    tmp_transfer_tree="${tree_dir}/tmp.$(basename "${label_tree}").transfer.$$"
    rm -f -- "${tmp_transfer_tree}"
    if nwkit root \
      --method transfer \
      --infile "${label_tree}" \
      --infile2 "${rooted_reference_tree}" \
      --outfile "${tmp_transfer_tree}"; then
      mv_out "${tmp_transfer_tree}" "${label_tree}"
      Rscript "${gg_support_dir}/nwk2pdf.r" --underbar2space=yes --italic=yes --infile="${label_tree}"
    else
      echo "Warning: Failed to transfer ${species_tree_rooting_method} root to ${label_tree}."
      rm -f -- "${tmp_transfer_tree}"
    fi
  done
}

normalize_iq2mc_constraint_tree() {
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
    rm -f -- "${tmpfile}"
    return 1
  fi

  mv_out "${tmpfile}" "${infile}"
  return 0
}

resolve_mcmctree_time_scale_factor() {
  local scale_source=""
  if [[ -n "${mcmctree_time_scale_factor_cache}" ]]; then
    printf '%s\n' "${mcmctree_time_scale_factor_cache}"
    return 0
  fi
  if [[ -s "${file_constrained_tree}" ]]; then
    scale_source="${file_constrained_tree}"
  elif [[ -s "${file_iq2mc_rooted_tree:-}" ]]; then
    scale_source="${file_iq2mc_rooted_tree}"
  fi
  if [[ -z "${scale_source}" ]]; then
    mcmctree_time_scale_factor_cache="1"
    printf '%s\n' "${mcmctree_time_scale_factor_cache}"
    return 0
  fi
  mcmctree_time_scale_factor_cache=$(python "${gg_support_dir}/mcmctree_time_scale.py" \
    factor \
    --infile "${scale_source}" \
    --target-max "10")
  if [[ -z "${mcmctree_time_scale_factor_cache}" ]]; then
    echo "Error: Failed to resolve MCMCTree time scale factor." >&2
    return 1
  fi
  printf '%s\n' "${mcmctree_time_scale_factor_cache}"
}

scale_mcmctree_calibrations_file() {
  local infile=$1
  local outfile=$2
  local scale_factor=$3
  local direction=$4
  python "${gg_support_dir}/mcmctree_time_scale.py" \
    scale-calibrations \
    --infile "${infile}" \
    --outfile "${outfile}" \
    --scale "${scale_factor}" \
    --direction "${direction}"
}

scale_mcmctree_ctl_rootage_file() {
  local infile=$1
  local outfile=$2
  local scale_factor=$3
  local direction=$4
  python "${gg_support_dir}/mcmctree_time_scale.py" \
    scale-ctl-rootage \
    --infile "${infile}" \
    --outfile "${outfile}" \
    --scale "${scale_factor}" \
    --direction "${direction}"
}

extract_scaled_mcmctree_figtree() {
  local infile=$1
  local outfile=$2
  local scale_factor=$3
  python "${gg_support_dir}/mcmctree_time_scale.py" \
    extract-figtree \
    --infile "${infile}" \
    --outfile "${outfile}" \
    --scale "${scale_factor}" \
    --direction "up"
}

mcmctree_requires_bdparas_flag() {
  local probe_dir
  local probe_stdout
  local probe_stderr
  local rc
  probe_dir=$(mktemp -d "${dir_tmp}/tmp.mcmctree.bdparas_probe.XXXXXX")
  probe_stdout="${probe_dir}/stdout.txt"
  probe_stderr="${probe_dir}/stderr.txt"
  printf '%s\n' "3 1" "((a,b)'B(0.1,0.2,0.025,0.025)',c);" > "${probe_dir}/tree.nwk"
  cat > "${probe_dir}/dummy.phy" << 'EOF'
3 4
a    ACGT
b    ACGT
c    ACGT
EOF
  cat > "${probe_dir}/mcmctree.ctl" << 'EOF'
seed = 1
seqfile = dummy.phy
treefile = tree.nwk
outfile = mcmctree.out
mcmcfile = mcmc.txt
ndata = 1
seqtype = 0
usedata = 0
clock = 2
RootAge = <0.3
model = 4
alpha = 0.5
ncatG = 5
cleandata = 1
BDparas = 1 1 0.5
kappa_gamma = 6 2
alpha_gamma = 1 1
rgene_gamma = 2 20 1
sigma2_gamma = 1 10 1
burnin = 1
sampfreq = 1
nsample = 1
EOF
  if (
    cd "${probe_dir}"
    mcmctree mcmctree.ctl > "${probe_stdout}" 2> "${probe_stderr}"
  ); then
    rc=0
  else
    rc=$?
  fi
  if [[ ${rc} -ne 0 ]] && grep -q "BDparas: expect flag" "${probe_stderr}"; then
    rm -rf -- "${probe_dir}"
    return 0
  fi
  rm -rf -- "${probe_dir}"
  return 1
}

normalize_mcmctree_ctl_for_installed_paml() {
  local ctl_file=$1
  if [[ ! -s "${ctl_file}" ]]; then
    return 0
  fi
  if ! mcmctree_requires_bdparas_flag; then
    return 0
  fi
  python - "${ctl_file}" << 'PY'
from __future__ import annotations

import re
import sys
from pathlib import Path

path = Path(sys.argv[1])
lines = path.read_text(encoding="utf-8").splitlines()
updated = []
pattern = re.compile(
    r"^(\s*BDparas\s*=\s*)"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)"
    r"(\s+[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)"
    r"(\s+[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)"
    r"(\s*(?:#.*)?)$"
)
for line in lines:
    match = pattern.match(line)
    if match:
        line = "".join(match.groups()[:-1]) + " M" + match.group(5)
    updated.append(line)
path.write_text("\n".join(updated) + "\n", encoding="utf-8")
PY
}

iq2mc_option_supported() {
  local candidate=$1
  local resolved_candidate
  local probe_prefix
  local probe_aln
  local probe_tree
  local probe_rc=1
  local probe_hessian=""
  resolved_candidate=$(command -v "${candidate}" 2> /dev/null || true)

  # Avoid pipefail+SIGPIPE false negatives from `... | grep -q`.
  if grep -q -- "--mcmc-bds" < <("${candidate}" -h 2>&1 || true); then
    return 0
  fi
  if command -v strings > /dev/null 2>&1 && [[ -n "${resolved_candidate}" ]]; then
    if grep -q -- "--mcmc-bds" < <(strings "${resolved_candidate}" 2> /dev/null || true); then
      return 0
    fi
  fi

  # Fallback: run a tiny dating job to detect support reliably.
  probe_prefix="${dir_tmp}/tmp.iq2mc_probe.$$"
  probe_aln="${probe_prefix}.fa"
  probe_tree="${probe_prefix}.nwk"
  probe_hessian="${probe_prefix}.mcmctree.hessian"
  cat > "${probe_aln}" << 'EOF'
>a
ACGT
>b
ACGT
>c
ACGT
EOF
  cat > "${probe_tree}" << 'EOF'
((a:0.1,b:0.1):0.1,c:0.1);
EOF
  if [[ -s "${probe_aln}" && -s "${probe_tree}" ]]; then
    if "${candidate}" \
      -s "${probe_aln}" \
      -te "${probe_tree}" \
      -m JC \
      -n 0 \
      -T 1 \
      --dating mcmctree \
      --mcmc-bds 1,1,0.5 \
      --prefix "${probe_prefix}" \
      --redo > /dev/null 2>&1; then
      probe_rc=0
    else
      probe_rc=$?
    fi
  fi
  if [[ ${probe_rc} -eq 0 && -s "${probe_hessian}" ]]; then
    rm -f -- "${probe_aln}" "${probe_tree}" "${probe_prefix}".* 2> /dev/null || true
    return 0
  fi
  rm -f -- "${probe_aln}" "${probe_tree}" "${probe_prefix}".* 2> /dev/null || true

  return 1
}

resolve_iq2mc_binary() {
  local candidate
  for candidate in iqtree3 iqtree; do
    if ! command -v "${candidate}" > /dev/null 2>&1; then
      continue
    fi
    if iq2mc_option_supported "${candidate}"; then
      echo "${candidate}"
      return 0
    fi
  done
  return 1
}

count_newick_tips() {
  local tree_file=$1
  python - "${tree_file}" << 'PY'
import sys
from Bio import Phylo
try:
    tree = Phylo.read(sys.argv[1], "newick")
    print(len(tree.get_terminals()))
except Exception:
    print(0)
PY
}

build_astral_input() {
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
    rm -f -- "${merged_file}"
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
    rm -f -- "${merged_file}"
    return 1
  fi
  return 0
}

optimize_astral_tree_branch_lengths() {
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

  rm -f -- "${tmp_topology}" "${iqtree_prefix}".* "${tmp_concat_alignment}" "${tmp_rooted_optimized_tree}" "${tmp_optimized_tree_with_support}"
  if ! python "${gg_support_dir}/write_newick_topology.py" \
    --infile "${astral_support_tree}" \
    --outfile "${tmp_topology}" \
    --parser 0; then
    echo "Warning: Failed to sanitize ASTRAL topology for IQ-TREE branch-length optimization."
    rm -f -- "${tmp_topology}" "${iqtree_prefix}".* "${tmp_concat_alignment}" "${tmp_rooted_optimized_tree}" "${tmp_optimized_tree_with_support}"
    return 1
  fi
  seqkit seq --threads 1 "${concat_alignment}" --out-file "${tmp_concat_alignment}"
  if [[ ! -s "${tmp_concat_alignment}" ]]; then
    echo "Warning: Failed to prepare alignment for IQ-TREE branch-length optimization: ${concat_alignment}"
    rm -f -- "${tmp_topology}" "${iqtree_prefix}".* "${tmp_concat_alignment}" "${tmp_rooted_optimized_tree}" "${tmp_optimized_tree_with_support}"
    return 1
  fi

  iqtree \
    -s "${tmp_concat_alignment}" \
    -te "${tmp_topology}" \
    -m "${model}" \
    -n 0 \
    -T "${GG_TASK_CPUS}" \
    "${iqtree_full_mem_args[@]}" \
    --prefix "${iqtree_prefix}" \
    --seed 12345 \
    --redo

  if [[ ! -s "${iqtree_prefix}.treefile" ]]; then
    echo "Warning: IQ-TREE optimization did not produce a tree for ${tag}."
    rm -f -- "${tmp_topology}" "${iqtree_prefix}".* "${tmp_concat_alignment}" "${tmp_rooted_optimized_tree}" "${tmp_optimized_tree_with_support}"
    return 1
  fi

  local root_exit_code=0
  if root_species_tree \
    "${iqtree_prefix}.treefile" \
    "${tmp_rooted_optimized_tree}" \
    "optimized ASTRAL ${tag} tree"; then
    root_exit_code=0
  else
    root_exit_code=$?
  fi
  if [[ ${root_exit_code} -ne 0 || ! -s "${tmp_rooted_optimized_tree}" ]]; then
    echo "Warning: Failed to save optimized ASTRAL ${tag} tree."
    rm -f -- "${tmp_topology}" "${iqtree_prefix}".* "${tmp_concat_alignment}" "${tmp_rooted_optimized_tree}" "${tmp_optimized_tree_with_support}"
    return 1
  fi

  local transfer_exit_code=0
  if nwkit transfer \
    --infile "${tmp_rooted_optimized_tree}" \
    --infile2 "${astral_support_tree}" \
    --target intnode \
    --support yes \
    --name no \
    --length no \
    --outfile "${tmp_optimized_tree_with_support}"; then
    transfer_exit_code=0
  else
    transfer_exit_code=$?
  fi
  local copy_exit_code=0
  if [[ ${transfer_exit_code} -eq 0 && -s "${tmp_optimized_tree_with_support}" ]]; then
    if cp_out "${tmp_optimized_tree_with_support}" "${optimized_outfile}"; then
      copy_exit_code=0
    else
      copy_exit_code=$?
    fi
  else
    echo "Warning: Failed to transfer ASTRAL support values to optimized ${tag} tree. Keeping branch-length-optimized tree without support transfer."
    if cp_out "${tmp_rooted_optimized_tree}" "${optimized_outfile}"; then
      copy_exit_code=0
    else
      copy_exit_code=$?
    fi
  fi
  rm -f -- "${tmp_topology}" "${iqtree_prefix}".* "${tmp_concat_alignment}" "${tmp_rooted_optimized_tree}" "${tmp_optimized_tree_with_support}"

  if [[ ${copy_exit_code} -ne 0 || ! -s "${optimized_outfile}" ]]; then
    echo "Warning: Failed to write optimized ASTRAL ${tag} tree."
    return 1
  fi
  return 0
}

orthofinder_output_directory_cleanup() {
  local target_dir=$1
  local _threads=${2:-1}
  if [[ -d "${target_dir}" ]]; then
    remove_empty_subdirs "${target_dir}"
  fi
}

detect_orthofinder_version() {
  local version_output version
  version_output=$(orthofinder -v 2>&1 || true)
  version_output=$(printf '%s\n' "${version_output}" | sed -E $'s/\x1B\\[[0-9;?]*[ -/]*[@-~]//g')
  version=$(printf '%s\n' "${version_output}" | awk '
    match($0, /[Oo]rtho[Ff]inder:?v?[[:space:]]*[0-9]+([.][0-9]+)*/) {
      version = substr($0, RSTART, RLENGTH)
      sub(/^.*[Oo]rtho[Ff]inder:?v?[[:space:]]*/, "", version)
      print version
      exit
    }
    match($0, /[Vv]ersion[[:space:]:]*v?[0-9]+([.][0-9]+)*/) {
      version = substr($0, RSTART, RLENGTH)
      sub(/^.*[Vv]ersion[[:space:]:]*v?/, "", version)
      print version
      exit
    }
  ')
  printf '%s\n' "${version}"
}

orthofinder_supports_root_hog_equivalent() {
  local version=$1
  local major minor
  if [[ ! "${version}" =~ ^([0-9]+)(\.([0-9]+))? ]]; then
    return 1
  fi

  major=${BASH_REMATCH[1]}
  minor=${BASH_REMATCH[3]:-0}
  [[ ${major} -gt 3 || ( ${major} -eq 3 && ${minor} -ge 1 ) ]]
}

copy_root_hog_equivalent_from_orthogroups() {
  local source_dir=$1
  local target_dir=$2
  local orthofinder_version=$3
  local source_og="${source_dir}/Orthogroups/Orthogroups.tsv"
  local source_genecount="${source_dir}/Orthogroups/Orthogroups.GeneCount.tsv"
  local target_readme="${target_dir}/README.txt"

  if [[ ! -s "${source_og}" || ! -s "${source_genecount}" ]]; then
    echo "OrthoFinder v3.1+ root-HOG-equivalent source files were not found."
    echo "Expected: ${source_og}"
    echo "Expected: ${source_genecount}"
    return 1
  fi

  ensure_dir "${target_dir}"
  echo "Treating OrthoFinder Orthogroups.tsv as the root-level HOG equivalent for OrthoFinder ${orthofinder_version}."
  cp_out "${source_og}" "${target_dir}/Orthogroups.tsv" || return 1
  cp_out "${source_genecount}" "${target_dir}/Orthogroups.GeneCount.tsv" || return 1
  {
    printf 'The files in this directory were created by gg_genome_evolution_core.sh.\n'
    printf 'Detected OrthoFinder version: %s\n' "${orthofinder_version}"
    printf 'For OrthoFinder version >= 3.1, Orthogroups/Orthogroups.tsv is treated as the root-level HOG equivalent.\n'
    printf 'OrthoFinder v3.1+ writes root-level HOG content into Orthogroups/Orthogroups.tsv and may omit Phylogenetic_Hierarchical_Orthogroups/N0.tsv from the final output.\n'
    printf 'Source Orthogroups.tsv: %s\n' "${source_og}"
    printf 'Source Orthogroups.GeneCount.tsv: %s\n' "${source_genecount}"
  } > "${target_readme}" || return 1
}

busco_notung() {
  infile=$1
  indir=$2
  outdir=$3
  busco_id=${infile%%.*}
  outfile="${outdir}/${busco_id}.busco.notung.root.zip"
  if [[ -s "${outfile}" ]]; then
    return 0
  fi
  if [[ -e "./${busco_id}.notung.root" ]]; then
    rm -rf -- "./${busco_id}.notung.root"
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
    rm -rf -- "./${busco_id}.notung.root"
  fi
}

busco_species_tree_assisted_gene_tree_rooting() {
  infile=$1
  indir=$2
  intreedir=$3
  outdir_txt=$4
  outdir_nwk=$5
  busco_id=${infile%%.*}
  intree="${intreedir}/${busco_id}.busco.nwk"
  outfile_txt="${outdir_txt}/${busco_id}.busco.root.txt"
  outfile_nwk="${outdir_nwk}/${busco_id}.busco.root.nwk"
  if [[ -s "${outfile_txt}" && -s "${outfile_nwk}" ]]; then
    return 0
  fi
  echo "Start NOTUNG root: ${busco_id}"
  if [[ -e "./${busco_id}.notung.root" ]]; then
    rm -rf -- "./${busco_id}.notung.root"
  fi
  cp_out "${indir}"/"${infile}" .
  unzip -q "${infile}"

  Rscript "${gg_support_dir}/species_tree_guided_gene_tree_rooting.r" \
    "--notung_root_zip=./${infile}" \
    "--in_tree=${intree}" \
    "--out_tree=${busco_id}.root.nwk" \
    "--species_parser=${species_label_parser}" \
    "--ncpu=${GG_TASK_CPUS}" \
    2>&1 | tee "${busco_id}.root.txt"

  if [[ -s "${busco_id}.root.nwk" ]]; then
    mv_out "${busco_id}".root.txt "${outfile_txt}"
    mv_out "${busco_id}".root.nwk "${outfile_nwk}"
    rm -f -- "${infile}"
    rm -rf -- "${busco_id}.notung.root"
  fi
}

busco_grampa() {
  indir=$1
  outdir=$2
  outfile=$3
  ensure_dir "${outdir}"
  local nwk_files=()
  shopt -s nullglob
  nwk_files=("${indir}"/*.nwk)
  shopt -u nullglob
  if [[ ${#nwk_files[@]} -eq 0 ]]; then
    echo "Skipping Grampa because no rooted gene trees were found in: ${indir}"
    return 0
  fi
  if [[ -e "./grampa_out" ]]; then
    rm -rf -- "./grampa_out"
  fi
  nwkit drop --infile "${file_dated_species_tree}" --target 'intnode' --name 'yes' |
    sed -e "s/_/-/g" \
      > "grampa_input_species_tree.nwk"

  : > "grampa_input_gene_trees.nwk"
  : > "busco_genetree_filenames.txt"
  for nwk_file in "${nwk_files[@]}"; do
    transformed_tree=$(
      sed -E "s/([(,])([^_)(,:]+)_([^_)(,:]+)_([^)(,:]+)([)(,:])/\1\4|||\2\-\3\5/g" "${nwk_file}" |
        sed -e "s/_/-/g" -e "s/|||/_/g" |
        sed -E 's/\)([^():;,]+):/\):/g; s/\)([^():;,]+);/\);/g' |
        tr -d '\r\n'
    )
    if [[ "${transformed_tree}" == *"("* && "${transformed_tree}" == *")"* && "${transformed_tree}" == *";"* ]]; then
      printf "%s\n" "${transformed_tree}" >> "grampa_input_gene_trees.nwk"
      printf "%s\n" "${nwk_file##*/}" >> "busco_genetree_filenames.txt"
    fi
  done

  valid_tree_count=$(awk 'NF>0{n++} END{print n+0}' "grampa_input_gene_trees.nwk")
  if [[ ${valid_tree_count} -eq 0 ]]; then
    echo "Skipping Grampa because no valid rooted gene trees were prepared from: ${indir}"
    rm -f -- "grampa_input_species_tree.nwk" "grampa_input_gene_trees.nwk" "busco_genetree_filenames.txt"
    return 0
  fi
  echo "Number of rooted gene trees passed to Grampa: ${valid_tree_count}"

  grampa_args=(
    -s "grampa_input_species_tree.nwk"
    -g "grampa_input_gene_trees.nwk"
    -o "./grampa_out"
    -p "${GG_TASK_CPUS}"
    -v -1
    --maps
  )
  if [[ -n "${grampa_h1}" ]]; then
    local grampa_h1_normalized=${grampa_h1//_/-}
    grampa_h1_normalized=${grampa_h1_normalized//[[:space:]]/-}
    grampa_args+=(-h1 "${grampa_h1_normalized}")
  fi

  grampa.py "${grampa_args[@]}"

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

  python "${gg_support_dir}/parse_grampa.py" \
    --grampa_det "${grampa_det_file}" \
    --grampa_out "${grampa_out_file}" \
    --gene_trees "./grampa_input_gene_trees.nwk" \
    --species_tree "./grampa_input_species_tree.nwk" \
    --ncpu "${GG_TASK_CPUS}" \
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
    rm -rf -- "./grampa_out"
  else
    echo "Grampa output files are missing. Ending Grampa."
  fi
}
