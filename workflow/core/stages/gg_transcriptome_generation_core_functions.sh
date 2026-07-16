# Named stage functions for gg_transcriptome_generation_core.sh.
# This file is sourced by workflow/core/gg_transcriptome_generation_core.sh.

resolve_busco_lineage_for_current_species() {
  if [[ -n "${busco_lineage_resolved}" ]]; then
    return 0
  fi
  if ! busco_lineage_resolved=$(gg_resolve_busco_lineage "${gg_workspace_dir}" "${busco_lineage}" "${sp_ub}"); then
    echo "Failed to resolve BUSCO lineage for species ${sp_ub} from request: ${busco_lineage}" >&2
    return 1
  fi
  echo "Resolved BUSCO lineage for ${sp_ub}: ${busco_lineage_resolved}"
}

capture_busco_failure_context() {
  local stage_key=$1
  local input_fasta=$2
  local stderr_log=$3
  local repro_dir=""

  repro_dir="${dir_transcriptome_assembly_output}/busco_repro/${sp_ub}/${stage_key}"
  echo "Capturing BUSCO repro artifacts to: ${repro_dir}"
  capture_busco_repro_artifacts \
    "${repro_dir}" \
    "${input_fasta}" \
    "./busco_tmp" \
    "${busco_lineage_resolved}" \
    "${stage_key}" \
    "${stderr_log}"
}

cleanup_busco_stage_temp_artifacts() {
  local input_fasta=$1
  rm -rf -- "./busco_tmp"
  rm -f -- "${input_fasta}"
  rm -f -- "./busco_tmp.stderr.log"
}

run_busco_with_capture() {
  local stage_key=$1
  local input_fasta=$2
  local stderr_log="./busco_tmp.stderr.log"

  if [[ -e "./busco_tmp" ]]; then
    rm -rf -- "./busco_tmp"
  fi
  if [[ -e "${stderr_log}" ]]; then
    rm -f -- "${stderr_log}"
  fi

  if ! gg_run_busco_with_metaeuk_modified_fas_compat \
    --in "${input_fasta}" \
    --mode transcriptome \
    --out "busco_tmp" \
    --cpu "${GG_TASK_CPUS}" \
    --force \
    --evalue 1e-03 \
    --limit 20 \
    --lineage_dataset "${dir_busco_lineage}" \
    --download_path "${dir_busco_db}" \
    --offline \
    2> >(tee "${stderr_log}" >&2); then
    if gg_busco_stderr_matches_known_metaeuk_modified_fas_bug "${stderr_log}"; then
      echo "BUSCO hit the known MetaEuk transcriptome bug for ${stage_key}. Capturing repro artifacts and continuing without BUSCO outputs."
      capture_busco_failure_context "${stage_key}" "${input_fasta}" "${stderr_log}" || return 1
      return 10
    fi
    echo "BUSCO failed for ${stage_key}. Capturing repro artifacts."
    capture_busco_failure_context "${stage_key}" "${input_fasta}" "${stderr_log}" || return 1
    return 1
  fi

  return 0
}

invalidate_cached_query_table_if_prefix_mismatch() {
  local table_file=$1
  local expected_prefix=$2
  local label=$3
  local header_lines=${4:-0}
  local first_query=""
  local first_query_species=""
  local expected_species=""
  local stale_file=""

  if [[ ! -s "${table_file}" ]]; then
    return 0
  fi
  first_query=$(awk -F '\t' -v skip="${header_lines}" 'NR > skip && $1 != "" { print $1; exit }' "${table_file}")
  if [[ -z "${first_query}" ]]; then
    return 0
  fi
  expected_species=${expected_prefix%_}
  first_query_species=$(gg_species_name_from_path_or_dot "${first_query}")
  if [[ "${first_query_species}" != "${expected_species}" ]]; then
    stale_file="${table_file}.stale.$(date +%Y%m%d%H%M%S)"
    mv -f -- "${table_file}" "${stale_file}"
    echo "Cached ${label} is inconsistent with the current species prefix ${expected_prefix}."
    echo "First query ID: ${first_query}"
    echo "Archived stale file to: ${stale_file}"
  fi
}

get_total_fastq_len_from_files() {
  if [[ $# -eq 0 ]]; then
    echo 0
    return 0
  fi

  local sum_len
  sum_len=$(seqkit stats --tabular "$@" \
    | awk -F '\t' '
      NR == 1 {
        for (i = 1; i <= NF; i++) {
          if ($i == "sum_len") {
            col = i
          }
        }
        next
      }
      NR > 1 && col > 0 {
        gsub(/,/, "", $col)
        sum += $col
      }
      END {
        printf "%.0f\n", sum + 0
      }
    ')
  echo "${sum_len}"
}

fastq_num_seqs_from_file() {
  local fastq_path=$1
  local num_seqs=""

  if ! num_seqs=$(seqkit stats --tabular "${fastq_path}" \
    | awk -F '\t' '
      NR == 1 {
        for (i = 1; i <= NF; i++) {
          if ($i == "num_seqs") {
            col = i
          }
        }
        next
      }
      NR > 1 && col > 0 {
        gsub(/,/, "", $col)
        print $col
        exit
      }
    '); then
    echo "Failed to count FASTQ reads with seqkit stats: ${fastq_path}" >&2
    return 1
  fi

  if [[ -z "${num_seqs}" ]]; then
    echo "Failed to parse num_seqs from seqkit stats output: ${fastq_path}" >&2
    return 1
  fi

  printf '%s\n' "${num_seqs}"
}

filter_valid_paired_fastq_files() {
  local report_file=$1
  local left_file=""
  local right_file=""
  local expected_left=""
  local left_count=""
  local right_count=""
  local skipped_count=0
  local -a valid_left=()
  local -a valid_right=()

  printf 'left_fastq\tright_fastq\tleft_reads\tright_reads\treason\n' > "${report_file}"

  for left_file in "${files_left[@]}"; do
    right_file="${left_file%_1.amalgkit.fastq.gz}_2.amalgkit.fastq.gz"
    if [[ ! -f "${right_file}" ]]; then
      printf '%s\t%s\tNA\tNA\tmissing_right\n' "${left_file}" "${right_file}" >> "${report_file}"
      skipped_count=$((skipped_count + 1))
      continue
    fi

    if ! left_count=$(fastq_num_seqs_from_file "${left_file}"); then
      printf '%s\t%s\tNA\tNA\tleft_count_failed\n' "${left_file}" "${right_file}" >> "${report_file}"
      skipped_count=$((skipped_count + 1))
      continue
    fi
    if ! right_count=$(fastq_num_seqs_from_file "${right_file}"); then
      printf '%s\t%s\t%s\tNA\tright_count_failed\n' "${left_file}" "${right_file}" "${left_count}" >> "${report_file}"
      skipped_count=$((skipped_count + 1))
      continue
    fi
    if [[ "${left_count}" != "${right_count}" ]]; then
      printf '%s\t%s\t%s\t%s\tread_count_mismatch\n' "${left_file}" "${right_file}" "${left_count}" "${right_count}" >> "${report_file}"
      skipped_count=$((skipped_count + 1))
      continue
    fi

    valid_left+=("${left_file}")
    valid_right+=("${right_file}")
  done

  for right_file in "${files_right[@]}"; do
    expected_left="${right_file%_2.amalgkit.fastq.gz}_1.amalgkit.fastq.gz"
    if [[ ! -f "${expected_left}" ]]; then
      printf '%s\t%s\tNA\tNA\tmissing_left\n' "${expected_left}" "${right_file}" >> "${report_file}"
      skipped_count=$((skipped_count + 1))
    fi
  done

  files_left=("${valid_left[@]}")
  files_right=("${valid_right[@]}")

  if [[ ${skipped_count} -gt 0 ]]; then
    echo "Excluded ${skipped_count} invalid or orphan paired-end FASTQ file set(s) before transcriptome assembly."
    echo "Paired-end FASTQ validation report: ${report_file}"
  else
    echo "Paired-end FASTQ validation passed for ${#files_left[@]} pair(s)."
  fi

  if [[ ${#files_left[@]} -eq 0 ]]; then
    echo "No valid paired-end FASTQ pairs remain for transcriptome assembly. Exiting."
    return 1
  fi

  return 0
}

csv_join_from_array() {
  local out=""
  local item=""
  for item in "$@"; do
    if [[ -n "${out}" ]]; then
      out="${out},"
    fi
    out="${out}${item}"
  done
  printf '%s\n' "${out}"
}

normalize_amalgkit_download_limit_value() {
  local raw_value=$1
  local option_label=$2
  local normalized_value=""

  normalized_value=$(printf '%s' "${raw_value}" | tr '[:upper:]' '[:lower:]' | tr -d '[:space:]')
  if [[ -z "${normalized_value}" || "${normalized_value}" == "auto" ]]; then
    printf 'auto\n'
    return 0
  fi
  if [[ ! "${normalized_value}" =~ ^[0-9]+$ ]]; then
    echo "Invalid ${option_label}: ${raw_value}. Use a non-negative integer or auto. Exiting." >&2
    exit 1
  fi
  if (( 10#${normalized_value} == 0 )); then
    printf '0\n'
    return 0
  fi
  printf '%d\n' "$((10#${normalized_value}))"
}

amalgkit_command_supports_option() {
  local subcommand=$1
  local option_name=$2
  local help_text=""
  help_text=$(amalgkit "${subcommand}" --help 2>&1 || true)
  grep -Fq -- "${option_name}" <<< "${help_text}"
}

require_amalgkit_supported_options() {
  local subcommand=$1
  shift
  local option_name=""

  for option_name in "$@"; do
    if ! amalgkit_command_supports_option "${subcommand}" "${option_name}"; then
      echo "Installed amalgkit ${subcommand} does not support ${option_name}. Update amalgkit to a build with shared download throttling. Exiting." >&2
      exit 1
    fi
  done
}

build_entrez_or_search_string_from_file() {
  local input_file=$1
  local joined_terms=""

  joined_terms=$(awk '
    {
      line = $0
      gsub(/\r/, "", line)
      gsub(/^[[:space:]]+/, "", line)
      gsub(/[[:space:]]+$/, "", line)
      if (line != "") {
        if (!first) {
          printf " OR "
        }
        printf "%s", line
        first = 0
      }
    }
    END {
      if (first == 0) {
        printf "\n"
      } else {
        exit 1
      }
    }
  ' "${input_file}") || return 1

  printf '(%s)\n' "${joined_terms}"
}

metadata_table_has_data_rows() {
  local table_file=$1

  if [[ ! -s "${table_file}" ]]; then
    return 1
  fi
  [[ $(wc -l < "${table_file}") -gt 1 ]]
}

extract_sraid_metadata_rows_for_species() {
  local metadata_source=$1
  local species_name=$2
  local output_file=$3

  {
    head -n 1 "${metadata_source}"
    grep -F -- "${species_name}" "${metadata_source}" || true
  } | sed -e "s/\t\t\tno\t/\tyes\tyes\tno\t/g" > "${output_file}"
}

extract_requested_accessions_missing_from_metadata() {
  local metadata_table=$1
  local accession_file=$2

  python "${gg_support_dir}/amalgkit_metadata_accessions.py" \
    missing \
    "${metadata_table}" \
    "${accession_file}"
}

extract_transcriptomic_rows_for_requested_accessions() {
  local metadata_source=$1
  local accession_file=$2
  local output_file=$3
  local raw_output_file="${output_file}.raw.$$"

  python "${gg_support_dir}/amalgkit_metadata_accessions.py" \
    extract-transcriptomic \
    "${metadata_source}" \
    "${accession_file}" \
    "${raw_output_file}"

  mv_out < <(sed -e "s/\t\t\tno\t/\tyes\tyes\tno\t/g" "${raw_output_file}") "${output_file}"
  rm -f -- "${raw_output_file}"
}

merge_metadata_tables_by_run() {
  local primary_table=$1
  local extra_table=$2
  local output_file=$3

  python "${gg_support_dir}/amalgkit_metadata_accessions.py" \
    merge-by-run \
    "${primary_table}" \
    "${extra_table}" \
    "${output_file}"
}

repair_private_fastq_metadata_scientific_names() {
  local metadata_file=$1
  local species_label=$2
  local support_dir=$3
  local output_file="${metadata_file}.resolved.$$"

  python - "${metadata_file}" "${species_label}" "${support_dir}" "${output_file}" <<'PY'
import csv
import sys
from pathlib import Path

PLACEHOLDER = "Please add in format: Genus species"

metadata_path = Path(sys.argv[1])
species_label = str(sys.argv[2] or "").strip()
support_dir = Path(sys.argv[3])
output_path = Path(sys.argv[4])

if not metadata_path.exists():
    raise SystemExit("Metadata file was not found: {}".format(metadata_path))
if species_label == "":
    raise SystemExit("Species label is empty; cannot repair private fastq metadata scientific_name values.")
if not support_dir.exists():
    raise SystemExit("Support directory was not found: {}".format(support_dir))

sys.path.insert(0, str(support_dir))
from species_labeling import scientific_name_from_label

fallback_name = scientific_name_from_label(species_label)
fallback_name = str(fallback_name or "").strip()
if fallback_name == "":
    fallback_name = species_label.replace("_", " ").strip()
if fallback_name == "" or fallback_name.lower() == PLACEHOLDER.lower():
    raise SystemExit(
        "Could not derive a fallback scientific_name from species label: {}".format(species_label)
    )

with metadata_path.open("rt", encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    fieldnames = reader.fieldnames
    if not fieldnames:
        raise SystemExit("Metadata table is missing a header: {}".format(metadata_path))
    if "scientific_name" not in fieldnames:
        raise SystemExit("Metadata table is missing required 'scientific_name' column: {}".format(metadata_path))
    rows = list(reader)

replaced_rows = 0
for row in rows:
    scientific_name = str(row.get("scientific_name", "") or "").strip()
    if scientific_name == "" or scientific_name.lower() == PLACEHOLDER.lower():
        row["scientific_name"] = fallback_name
        replaced_rows += 1

output_path.parent.mkdir(parents=True, exist_ok=True)
with output_path.open("wt", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()
    for row in rows:
        writer.writerow({name: row.get(name, "") or "" for name in fieldnames})

if replaced_rows > 0:
    print(
        "Filled {} private fastq metadata row(s) with fallback scientific_name: {}".format(
            replaced_rows,
            fallback_name,
        )
    )
else:
    print("Private fastq metadata already had non-placeholder scientific_name values.")
PY

  mv_out "${output_file}" "${metadata_file}"
}

stage_quant_reference_fasta_aliases() {
  local metadata_file=$1
  local reference_fasta=$2
  local output_dir=$3
  local canonical_prefix=$4
  local support_dir=${5:-${gg_support_dir:-}}

  python - "${metadata_file}" "${reference_fasta}" "${output_dir}" "${canonical_prefix}" "${support_dir}" <<'PY'
import csv
import os
import re
import sys
from pathlib import Path

metadata_path = Path(sys.argv[1])
reference_path = Path(sys.argv[2])
output_dir = Path(sys.argv[3])
canonical_prefix = sys.argv[4]
support_dir = Path(sys.argv[5])

if not reference_path.exists():
    raise SystemExit("Reference fasta file was not found: {}".format(reference_path))
if output_dir.exists() and (not output_dir.is_dir()):
    raise SystemExit("Quant fasta staging path is not a directory: {}".format(output_dir))
if not support_dir.exists():
    raise SystemExit("Support directory was not found: {}".format(support_dir))

sys.path.insert(0, str(support_dir))
from species_labeling import base_species_label, species_label_from_taxonomic_text

def normalize_prefix(raw_value):
    text = str(raw_value or "").strip()
    if not text:
        return ""
    species_label = species_label_from_taxonomic_text(text)
    if species_label != "":
        return species_label
    text = re.sub(r"\s+", "_", text)
    text = text.replace("/", "_").replace("\\", "_")
    return text

canonical_prefix_normalized = normalize_prefix(canonical_prefix)
canonical_base_prefix = base_species_label(canonical_prefix_normalized)
if canonical_base_prefix == "":
    canonical_base_prefix = canonical_prefix_normalized
canonical_is_species_level = canonical_prefix_normalized == canonical_base_prefix
metadata_prefixes = []
seen_prefixes = set()

def add_metadata_prefix(raw_value):
    prefix = normalize_prefix(raw_value)
    if (not prefix) or (prefix in seen_prefixes):
        return
    seen_prefixes.add(prefix)
    metadata_prefixes.append(prefix)

if metadata_path.exists() and metadata_path.stat().st_size > 0:
    with metadata_path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        candidate_fields = ["scientific_name"] if "scientific_name" in fieldnames else []
        for row in reader:
            for fieldname in candidate_fields:
                add_metadata_prefix(row.get(fieldname, ""))

prefixes_to_stage = metadata_prefixes if metadata_prefixes else [canonical_prefix_normalized]
invalid_prefixes = []
for prefix in prefixes_to_stage:
    prefix_base = base_species_label(prefix)
    if prefix_base == "":
        prefix_base = prefix
    is_exact = prefix == canonical_prefix_normalized
    is_same_base = canonical_is_species_level and prefix_base == canonical_prefix_normalized
    if not (is_exact or is_same_base):
        invalid_prefixes.append(prefix)

if invalid_prefixes:
    raise SystemExit(
        "Metadata contains scientific_name prefixes outside GeneGalleon species_key '{}': {} (all prefixes: {})".format(
            canonical_prefix_normalized,
            ", ".join(invalid_prefixes),
            ", ".join(prefixes_to_stage),
        )
    )

if any(prefix == "" for prefix in prefixes_to_stage):
    raise SystemExit("Could not determine quant reference prefix for species_key: {}".format(canonical_prefix))

output_dir.mkdir(parents=True, exist_ok=True)
for prefix in prefixes_to_stage:
    alias_path = output_dir / "{}_for_kallisto_index.fasta".format(prefix)
    if alias_path.exists() or alias_path.is_symlink():
        if alias_path.resolve() != reference_path.resolve():
            raise SystemExit("Quant reference path already exists for a different target: {}".format(alias_path))
    else:
        os.symlink(reference_path, alias_path)
    print(alias_path.name)
PY
}

stage_amalgkit_merge_metadata_for_species() {
  local metadata_file=$1
  local output_file=$2
  local canonical_prefix=$3
  local support_dir=${4:-${gg_support_dir:-}}

  python - "${metadata_file}" "${output_file}" "${canonical_prefix}" "${support_dir}" <<'PY'
import csv
import re
import sys
from pathlib import Path

metadata_path = Path(sys.argv[1])
output_path = Path(sys.argv[2])
canonical_prefix = sys.argv[3]
support_dir = Path(sys.argv[4])

if not metadata_path.exists():
    raise SystemExit("Metadata file was not found: {}".format(metadata_path))
if not support_dir.exists():
    raise SystemExit("Support directory was not found: {}".format(support_dir))

sys.path.insert(0, str(support_dir))
from species_labeling import base_species_label, scientific_name_from_label, species_label_from_taxonomic_text

def normalize_prefix(raw_value):
    text = str(raw_value or "").strip()
    if not text:
        return ""
    species_label = species_label_from_taxonomic_text(text)
    if species_label != "":
        return species_label
    text = re.sub(r"\s+", "_", text)
    text = text.replace("/", "_").replace("\\", "_")
    return text

canonical_prefix_normalized = normalize_prefix(canonical_prefix)
canonical_base_prefix = base_species_label(canonical_prefix_normalized)
if canonical_base_prefix == "":
    canonical_base_prefix = canonical_prefix_normalized
canonical_is_species_level = canonical_prefix_normalized == canonical_base_prefix
canonical_scientific_name = scientific_name_from_label(canonical_prefix_normalized)
if canonical_scientific_name == "":
    canonical_scientific_name = canonical_prefix_normalized.replace("_", " ")

with metadata_path.open("rt", encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    fieldnames = reader.fieldnames
    if not fieldnames:
        raise SystemExit("Metadata table is missing a header: {}".format(metadata_path))
    if "scientific_name" not in fieldnames:
        raise SystemExit("Metadata table is missing required 'scientific_name' column: {}".format(metadata_path))
    rows = list(reader)

invalid_prefixes = []
seen_invalid_prefixes = set()
for row in rows:
    prefix = normalize_prefix(row.get("scientific_name", ""))
    prefix_base = base_species_label(prefix)
    if prefix_base == "":
        prefix_base = prefix
    is_exact = prefix == canonical_prefix_normalized
    is_same_base = canonical_is_species_level and prefix_base == canonical_prefix_normalized
    if not (is_exact or is_same_base):
        if prefix not in seen_invalid_prefixes:
            seen_invalid_prefixes.add(prefix)
            invalid_prefixes.append(prefix)
        continue
    row["scientific_name"] = canonical_scientific_name

if invalid_prefixes:
    raise SystemExit(
        "Metadata contains scientific_name prefixes outside GeneGalleon species_key '{}': {}".format(
            canonical_prefix_normalized,
            ", ".join(invalid_prefixes),
        )
    )

output_path.parent.mkdir(parents=True, exist_ok=True)
with output_path.open("wt", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=fieldnames,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()
    for row in rows:
        writer.writerow({name: row.get(name, "") or "" for name in fieldnames})
PY
}

resolve_amalgkit_merge_output_prefix() {
  local metadata_file=$1
  local merge_dir=$2
  local canonical_prefix=$3

  python - "${metadata_file}" "${merge_dir}" "${canonical_prefix}" <<'PY'
import csv
import re
import sys
from pathlib import Path

metadata_path = Path(sys.argv[1])
merge_dir = Path(sys.argv[2])
canonical_prefix = sys.argv[3]

def normalize_prefix(raw_value):
    text = str(raw_value or "").strip()
    if not text:
        return ""
    text = re.sub(r"\s+", "_", text)
    text = text.replace("/", "_").replace("\\", "_")
    return text

canonical_prefix_normalized = normalize_prefix(canonical_prefix)
metadata_prefixes = []
seen_prefixes = set()

def add_metadata_prefix(raw_value):
    prefix = normalize_prefix(raw_value)
    if (not prefix) or (prefix in seen_prefixes):
        return
    seen_prefixes.add(prefix)
    metadata_prefixes.append(prefix)

if metadata_path.exists() and metadata_path.stat().st_size > 0:
    with metadata_path.open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        candidate_fields = ["scientific_name"] if "scientific_name" in fieldnames else []
        for row in reader:
            for fieldname in candidate_fields:
                add_metadata_prefix(row.get(fieldname, ""))

if len(metadata_prefixes) > 1:
    raise SystemExit(
        "Metadata contains multiple scientific_name prefixes for one GeneGalleon species_key '{}': {}".format(
            canonical_prefix_normalized,
            ", ".join(metadata_prefixes),
        )
    )

prefix = metadata_prefixes[0] if metadata_prefixes else canonical_prefix_normalized
candidate = merge_dir / prefix / f"{prefix}_eff_length.tsv"
if candidate.is_file() and candidate.stat().st_size > 0:
    print(prefix)
    raise SystemExit(0)

raise SystemExit(1)
PY
}

run_amalgkit_metadata_query() {
  local search_string=$1
  local metadata_cmd=()

  echo "Entrez search string: ${search_string}"
  echo "amalgkit metadata NCBI metadata max concurrency: ${amalgkit_ncbi_metadata_max_concurrency}"
  require_amalgkit_supported_options "metadata" "--download_lock_dir" "--ncbi_metadata_max_concurrency"

  metadata_cmd=(
    amalgkit metadata
    --out_dir "./"
    --download_dir "${dir_amalgkit_download_dir}"
    --download_lock_dir "${dir_amalgkit_download_lock_dir}"
    --search_string "${search_string}"
    --ncbi_metadata_max_concurrency "${amalgkit_ncbi_metadata_max_concurrency}"
  )
  "${metadata_cmd[@]}"
}

cleanup_partial_getfastq_outputs() {
  rm -rf -- "${dir_tmp}/getfastq"
  rm -rf -- "${dir_amalgkit_getfastq_sp}"
  ensure_dir "${dir_amalgkit_getfastq_sp}"
}

clear_getfastq_safely_removed_markers() {
  local safely_removed_files=()
  if [[ -e "${file_amalgkit_getfastq_legacy_safely_removed_flag}" ]]; then
    rm -f -- "${file_amalgkit_getfastq_legacy_safely_removed_flag}"
  fi
  if [[ -d "${dir_amalgkit_getfastq_sp}" ]]; then
    mapfile -t safely_removed_files < <(find "${dir_amalgkit_getfastq_sp}" -type f -name "*.safely_removed" | sort)
    if [[ ${#safely_removed_files[@]} -gt 0 ]]; then
      rm -f -- "${safely_removed_files[@]}"
    fi
  fi
}

safe_delete_getfastq_fastq_files() {
  local fastq_files=()
  local fastq_file
  if [[ -e "${file_amalgkit_getfastq_legacy_safely_removed_flag}" ]]; then
    rm -f -- "${file_amalgkit_getfastq_legacy_safely_removed_flag}"
  fi
  if [[ ! -d "${dir_amalgkit_getfastq_sp}" ]]; then
    return 0
  fi
  mapfile -t fastq_files < <(find "${dir_amalgkit_getfastq_sp}" -type f -name "*.amalgkit.fastq.gz" | sort)
  if [[ ${#fastq_files[@]} -eq 0 ]]; then
    echo "No amalgkit getfastq FASTQ files were found to safe-delete in: ${dir_amalgkit_getfastq_sp}"
    return 0
  fi
  for fastq_file in "${fastq_files[@]}"; do
    rm -f -- "${fastq_file}"
    printf '%s\n' "This fastq file was safely removed after downstream completion in gg_transcriptome_generation." > "${fastq_file}.safely_removed"
  done
}

amalgkit_getfastq_log_has_fatal_message() {
  local log_file=$1
  [[ -s "${log_file}" ]] || return 1
  grep -Eq '^ERROR: ' "${log_file}"
}

run_amalgkit_getfastq_attempt() {
  local rrna_filter_value=$1
  local attempt_label=$2
  local log_file="${dir_tmp}/amalgkit_getfastq.${attempt_label}.log"
  local getfastq_outputs=()
  local getfastq_cmd=()

  rm -f -- "${log_file}"
  echo "Running amalgkit getfastq attempt '${attempt_label}' with --rrna_filter ${rrna_filter_value}"
  echo "amalgkit getfastq NCBI metadata max concurrency: ${amalgkit_ncbi_metadata_max_concurrency}"
  echo "amalgkit getfastq NCBI download max concurrency: ${amalgkit_ncbi_download_max_concurrency}"
  echo "amalgkit getfastq AWS download max concurrency: ${amalgkit_aws_download_max_concurrency}"
  echo "amalgkit getfastq GCP download max concurrency: ${amalgkit_gcp_download_max_concurrency}"
  require_amalgkit_supported_options "getfastq" \
    "--download_lock_dir" \
    "--ncbi_metadata_max_concurrency" \
    "--ncbi_download_max_concurrency" \
    "--aws_download_max_concurrency" \
    "--gcp_download_max_concurrency"

  getfastq_cmd=(
    amalgkit getfastq
    --out_dir "${dir_tmp}"
    --download_dir "${dir_amalgkit_download_dir}"
    --download_lock_dir "${dir_amalgkit_download_lock_dir}"
    --metadata "${file_amalgkit_metadata}"
    --threads "${GG_TASK_CPUS}"
    --ncbi_metadata_max_concurrency "${amalgkit_ncbi_metadata_max_concurrency}"
    --ncbi_download_max_concurrency "${amalgkit_ncbi_download_max_concurrency}"
    --aws_download_max_concurrency "${amalgkit_aws_download_max_concurrency}"
    --gcp_download_max_concurrency "${amalgkit_gcp_download_max_concurrency}"
    --rrna_filter "${rrna_filter_value}"
    --contam_filter "${amalgkit_contam_filter}"
    --contam_filter_rank "${contamination_removal_rank_for_amalgkit}"
    --contam_filter_db "${dir_mmseqs2_db}/UniRef90_DB"
    --remove_sra yes
    --remove_tmp yes
    --dump_print yes
    --read_name 'trinity'
    --aws yes
    --ncbi yes
    --redo no
  )

  if "${getfastq_cmd[@]}" > "${log_file}" 2>&1; then
    if amalgkit_getfastq_log_has_fatal_message "${log_file}"; then
      echo "Detected fatal message in amalgkit getfastq log despite a zero exit code: ${log_file}"
      grep -E '^ERROR: ' "${log_file}" || true
      cleanup_partial_getfastq_outputs
      return 2
    fi
    echo "amalgkit getfastq safely finished."
    echo "amalgkit getfastq log: ${log_file}"
    shopt -s nullglob
    getfastq_outputs=("${dir_tmp}"/getfastq/*)
    shopt -u nullglob
    if [[ ${#getfastq_outputs[@]} -eq 0 ]]; then
      echo "amalgkit getfastq finished but no output files were found under: ${dir_tmp}/getfastq"
      cleanup_partial_getfastq_outputs
      return 1
    fi
    mv_out_replace_dir "${dir_tmp}/getfastq" "${dir_amalgkit_getfastq_sp}"
    rm -rf -- "${dir_tmp}/getfastq"
    return 0
  fi

  local status_amalgkit_attempt=$?
  echo "amalgkit getfastq exit code: ${status_amalgkit_attempt}"
  echo "amalgkit getfastq log: ${log_file}"
  if amalgkit_getfastq_log_has_fatal_message "${log_file}"; then
    echo "Detected fatal message in amalgkit getfastq log."
    grep -E '^ERROR: ' "${log_file}" || true
    cleanup_partial_getfastq_outputs
    return 2
  fi
  tail -n 50 "${log_file}" || true
  cleanup_partial_getfastq_outputs
  return 1
}

detect_transcriptome_read_technology_from_metadata() {
  local metadata_tsv="$1"
  local classification_tsv="$2"
  local summary_sh="$3"

  python "${gg_support_dir}/detect_amalgkit_read_technology.py" \
    --metadata "${metadata_tsv}" \
    --classification-out "${classification_tsv}" \
    --summary-sh "${summary_sh}"

  # shellcheck disable=SC1090
  source "${summary_sh}"
}

load_classified_getfastq_files() {
  local classification_tsv="$1"
  local run=""
  local read_class=""
  local _platform_family=""
  local _long_read=""
  local _ont_direct_rna=""
  local _lib_layout=""
  local _instrument_value=""

  classified_short_single_fastq_files=()
  classified_short_left_fastq_files=()
  classified_short_right_fastq_files=()
  classified_long_fastq_files=()
  classified_pacbio_fastq_files=()
  classified_ont_fastq_files=()

  while IFS=$'\t' read -r run read_class _platform_family _long_read _ont_direct_rna _lib_layout _instrument_value; do
    local run_matches=()
    local match=""
    local single_fastq=""
    local left_fastq=""
    local right_fastq=""

    if [[ "${run}" == "run" || -z "${run}" ]]; then
      continue
    fi

    while IFS= read -r match; do
      run_matches+=( "${match}" )
    done < <(
      find "${dir_amalgkit_getfastq_sp}" -type f \
        \( -name "${run}.amalgkit.fastq.gz" -o -name "${run}_1.amalgkit.fastq.gz" -o -name "${run}_2.amalgkit.fastq.gz" \) \
        | sort
    )

    if [[ ${#run_matches[@]} -eq 0 ]]; then
      continue
    fi

    for match in "${run_matches[@]}"; do
      case "$(basename "${match}")" in
        "${run}.amalgkit.fastq.gz")
          single_fastq="${match}"
          ;;
        "${run}_1.amalgkit.fastq.gz")
          left_fastq="${match}"
          ;;
        "${run}_2.amalgkit.fastq.gz")
          right_fastq="${match}"
          ;;
      esac
    done

    case "${read_class}" in
      short_read)
        if [[ -n "${left_fastq}" && -n "${right_fastq}" ]]; then
          classified_short_left_fastq_files+=( "${left_fastq}" )
          classified_short_right_fastq_files+=( "${right_fastq}" )
        elif [[ -n "${single_fastq}" ]]; then
          classified_short_single_fastq_files+=( "${single_fastq}" )
        fi
        ;;
      pacbio|ont_cdna|ont_direct_rna|long_read_unknown)
        if [[ -n "${single_fastq}" ]]; then
          classified_long_fastq_files+=( "${single_fastq}" )
          if [[ "${read_class}" == "pacbio" ]]; then
            classified_pacbio_fastq_files+=( "${single_fastq}" )
          elif [[ "${read_class}" == "ont_cdna" || "${read_class}" == "ont_direct_rna" ]]; then
            classified_ont_fastq_files+=( "${single_fastq}" )
          fi
        else
          echo "Warning: long-read run ${run} did not produce a single-end FASTQ. Passing all matching files through RNA-Bloom2/Corset."
          for match in "${run_matches[@]}"; do
            classified_long_fastq_files+=( "${match}" )
            if [[ "${read_class}" == "pacbio" ]]; then
              classified_pacbio_fastq_files+=( "${match}" )
            elif [[ "${read_class}" == "ont_cdna" || "${read_class}" == "ont_direct_rna" ]]; then
              classified_ont_fastq_files+=( "${match}" )
            fi
          done
        fi
        ;;
    esac
  done < "${classification_tsv}"
}

configure_transcriptome_runtime_from_detected_metadata() {
  case "${requested_assembly_method}" in
    auto)
      if [[ ${detected_has_long_reads} -eq 1 ]]; then
        effective_assembly_method="rna-bloom2"
      else
        effective_assembly_method="rnaspades"
      fi
      ;;
    trinity|rnaspades)
      effective_assembly_method="${requested_assembly_method}"
      ;;
    rna-bloom2|rnabloom2)
      effective_assembly_method="rna-bloom2"
      ;;
    *)
      echo "Invalid value for 'assembly_method'. Please specify one of auto, Trinity, rnaSPAdes, or RNA-Bloom2."
      exit 1
      ;;
  esac

  echo "Detected metadata instrument column: ${detected_metadata_instrument_field}"
  echo "Detected transcriptome input class: ${detected_input_class}"
  echo "Detected run counts: total=${detected_metadata_run_count} short=${detected_short_read_run_count} pacbio=${detected_pacbio_run_count} ont_cdna=${detected_ont_cdna_run_count} ont_direct_rna=${detected_ont_direct_rna_run_count} long_read_unknown=${detected_long_read_unknown_run_count}"

  if [[ ${detected_has_pacbio} -eq 1 && ${detected_has_ont} -eq 1 ]]; then
    echo "Mixed PacBio and ONT long-read runs were detected in metadata: ${file_amalgkit_metadata}. Split these technologies into separate transcriptome-generation tasks. Exiting."
    exit 1
  fi
  if [[ ${detected_has_ont_cdna} -eq 1 && ${detected_has_ont_direct_rna} -eq 1 ]]; then
    echo "Mixed ONT cDNA and direct-RNA runs were detected in metadata: ${file_amalgkit_metadata}. Split these protocols into separate transcriptome-generation tasks. Exiting."
    exit 1
  fi
  if [[ ${detected_long_read_unknown_run_count} -gt 0 && ( ${detected_has_pacbio} -eq 1 || ${detected_has_ont} -eq 1 ) ]]; then
    echo "Metadata contains long-read runs whose platform could not be resolved to PacBio or ONT alongside explicitly classified long-read runs: ${file_amalgkit_metadata}. Split or relabel these runs before continuing. Exiting."
    exit 1
  fi

  if [[ ${detected_has_long_reads} -eq 1 ]]; then
    if [[ "${effective_assembly_method}" != "rna-bloom2" ]]; then
      echo "Long-read platforms were detected from metadata, but assembly_method=${assembly_method:-auto} resolves to ${effective_assembly_method}."
      echo "Use assembly_method=auto or assembly_method=RNA-Bloom2 for PacBio/ONT transcriptome assembly. Exiting."
      exit 1
    fi
    if [[ ${detected_has_short_reads} -eq 1 ]]; then
      echo "Mixed short-read and long-read runs were detected. RNA-Bloom2 assembly will use the long-read FASTQ files only."
    fi
    if [[ ${detected_long_read_unknown_run_count} -gt 0 ]]; then
      echo "Long-read runs with unresolved PacBio/ONT platform were inferred from metadata length heuristics."
      if [[ ${run_amalgkit_quant} -eq 1 && "${amalgkit_oarfish_seq_tech}" == "auto" ]]; then
        echo "amalgkit quant cannot auto-resolve oarfish sequencing technology for these runs. Set amalgkit_oarfish_seq_tech explicitly or disable run_amalgkit_quant. Exiting."
        exit 1
      fi
    fi
    if [[ ${run_amalgkit_quant} -eq 1 ]]; then
      echo "Detected long-read platforms in metadata. run_amalgkit_quant remains enabled; amalgkit quant will use quant_backend=${amalgkit_quant_backend}."
    fi
    if [[ ${run_amalgkit_merge} -eq 1 ]]; then
      echo "Detected long-read platforms in metadata. run_amalgkit_merge remains enabled because amalgkit merge accepts normalized abundance tables from long-read quant."
    fi
  else
    if [[ "${effective_assembly_method}" == "rna-bloom2" ]]; then
      echo "assembly_method=RNA-Bloom2 is currently reserved for long-read metadata-detected inputs. The metadata for ${sp_ub} did not indicate PacBio/ONT runs. Exiting."
      exit 1
    fi
  fi

  echo "Effective assembly method: ${effective_assembly_method}"
}

download_public_original_fastqs_for_metadata() {
  local metadata_tsv="$1"
  local output_dir="$2"
  python - "${metadata_tsv}" "${output_dir}" <<'PY'
import csv
import gzip
import shutil
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

metadata_path = Path(sys.argv[1])
output_root = Path(sys.argv[2])


def fetch_bytes(url: str) -> bytes:
    last_exc = None
    for attempt in range(1, 6):
        try:
            with urllib.request.urlopen(url, timeout=120) as response:
                return response.read()
        except Exception as exc:  # pragma: no cover - exercised via shell integration
            last_exc = exc
            if attempt == 5:
                raise
            time.sleep(2)
    raise last_exc  # pragma: no cover


def fetch_text(url: str) -> str:
    return fetch_bytes(url).decode("utf-8", "replace")


def sort_key(item):
    name = item[0].lower()
    if "forward" in name or "_1" in name or "r1" in name:
        return (0, name)
    if "reverse" in name or "_2" in name or "r2" in name:
        return (1, name)
    return (2, name)


with metadata_path.open("rt", encoding="utf-8", newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    fieldnames = list(reader.fieldnames or [])
    if "run" not in fieldnames:
        raise SystemExit("Metadata is missing required 'run' column: {}".format(metadata_path))
    runs = []
    seen = set()
    for row in reader:
        run = str(row.get("run", "") or "").strip()
        if not run or run in seen:
            continue
        seen.add(run)
        runs.append(run)

if not runs:
    raise SystemExit("No run accessions were found in metadata: {}".format(metadata_path))

output_root.mkdir(parents=True, exist_ok=True)

for run in runs:
    xml_url = "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/run_new?acc={}".format(
        urllib.parse.quote(run, safe="")
    )
    root = ET.fromstring(fetch_text(xml_url))
    fastq_files = []
    for node in root.iter("SRAFile"):
        if node.attrib.get("semantic_name") != "fastq":
            continue
        if node.attrib.get("supertype") != "Original":
            continue
        url = str(node.attrib.get("url", "") or "").strip()
        filename = str(node.attrib.get("filename", "") or "").strip()
        if not url.startswith("https://"):
            for alt in node.findall("Alternatives"):
                alt_url = str(alt.attrib.get("url", "") or "").strip()
                if alt_url.startswith("https://"):
                    url = alt_url
                    break
        if url.startswith("https://"):
            fastq_files.append((filename, url))

    if not fastq_files:
        raise SystemExit("No public original FASTQ URLs were found for run: {}".format(run))
    if len(fastq_files) > 2:
        raise SystemExit("Unexpected number of original FASTQ files for run {}: {}".format(run, len(fastq_files)))

    run_dir = output_root / run
    if run_dir.exists():
        shutil.rmtree(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    fastq_files.sort(key=sort_key)

    for idx, (filename, url) in enumerate(fastq_files, start=1):
        if len(fastq_files) == 1:
            dest = run_dir / "{}.amalgkit.fastq.gz".format(run)
        else:
            dest = run_dir / "{}_{}.amalgkit.fastq.gz".format(run, idx)
        payload = fetch_bytes(url)
        if payload[:2] == b"\x1f\x8b":
            dest.write_bytes(payload)
        else:
            with gzip.open(dest, "wb") as handle_out:
                handle_out.write(payload)
        print("Recovered original FASTQ for {}: {} -> {}".format(run, filename or url, dest))
PY
}

run_amalgkit_getfastq_or_fallback() {
  local status_amalgkit=0
  local fatal_retry_suffix=""

  if run_amalgkit_getfastq_attempt "${amalgkit_rrna_filter}" "initial"; then
    return 0
  else
    status_amalgkit=$?
  fi

  if [[ ${status_amalgkit} -eq 2 ]]; then
    if [[ "${amalgkit_rrna_filter}" != "no" ]]; then
      echo "Retrying amalgkit getfastq once with --rrna_filter no because the previous attempt logged a fatal error."
      fatal_retry_suffix=" after retrying with --rrna_filter no"
      if run_amalgkit_getfastq_attempt "no" "retry_rrna_filter_no"; then
        return 0
      else
        status_amalgkit=$?
      fi
    fi
    echo "amalgkit getfastq encountered a fatal error${fatal_retry_suffix}. Exiting without fallback download so partial outputs do not reach downstream steps."
    cleanup_partial_getfastq_outputs
    return 1
  fi

  echo "amalgkit getfastq did not safely finish. Attempting fallback download of public original FASTQ files."
  cleanup_partial_getfastq_outputs
  if download_public_original_fastqs_for_metadata "${file_amalgkit_metadata}" "${dir_amalgkit_getfastq_sp}"; then
    echo "Fallback download of public original FASTQ files succeeded."
    return 0
  fi
  echo "Fallback direct FASTQ recovery also failed. Exiting."
  return 1
}
