# shellcheck shell=bash
# Sourced by workflow/support/gg_util.sh.

exit_if_running_qstat() {
	local flag_start=0
	local qstat_output=""
	if [[ "${exit_if_running:-0}" -eq 1 ]]; then
		flag_start=1
	fi
	if ! command -v qstat >/dev/null 2>&1; then
		echo "qstat was not found. exit_if_running is deactivated."
		return 0
	fi
	if ! qstat > /dev/null 2>&1; then
		echo "qstat command failed. exit_if_running is deactivated."
		return 0
	fi
	echo "qstat was found."
	if [[ ${flag_start} -ne 1 ]]; then
		return 0
	fi

	qstat_output=$(qstat 2>/dev/null || true)
	if [[ -z "${qstat_output}" ]]; then
		echo "qstat output was empty. Skipping duplicate job check."
		return 0
	fi

	local job_name
	job_name=$(
		printf '%s\n' "${qstat_output}" \
		| awk -v job_id="${GG_JOB_ID:-}" 'NR>2 && $1==job_id && !seen {print $3; seen=1}'
	)
	if [[ -z "${job_name}" ]]; then
		echo "Could not parse job name from qstat for GG_JOB_ID=${GG_JOB_ID:-NA}. Skipping duplicate job check."
		return 0
	fi

	local running_id=()
	local task_token
	while IFS= read -r task_token; do
		running_id+=( "${task_token}" )
	done < <(
		printf '%s\n' "${qstat_output}" \
		| awk \
		-v target_job_name="${job_name}" \
		-v this_job_id="${GG_JOB_ID:-}" \
		-v this_task_id="${GG_ARRAY_TASK_ID:-1}" \
		'
		NR<=2 {next}
		$0 ~ / dr / {next}
		$3 != target_job_name {next}
		{
			task_id = (NF >= 11 ? $11 : "")
			if (task_id == "" || task_id ~ /:/) next
			if ($1 == this_job_id && task_id == this_task_id) next
			print task_id
		}
		' \
		| sort -u
	)
	echo "running_id: ${running_id[*]}"
	local flag=1
	local r
	for r in "${running_id[@]}"; do
		if [[ "${r}" =~ ^[0-9]+$ ]] && [[ "${r}" -eq "${GG_ARRAY_TASK_ID:-1}" ]]; then
			flag=0
		fi
	done
	if [[ ${flag} -eq 0 ]]; then
		echo "GG_ARRAY_TASK_ID=${GG_ARRAY_TASK_ID:-1} is running already. Exiting."
		exit 0
	fi
}

wait_until_jobn_le() {
  # https://qiita.com/jokester/items/34ce222e1702c6e120eb
  local max_jobn=$1
  while [[ "$(jobs | wc -l)" -gt "$max_jobn" ]] ; do
    sleep 1
  done
}

wait_for_background_jobs() {
  # `wait` without pid can miss failures when a non-final background job exits non-zero.
  local pids=()
  while IFS= read -r pid; do
    if [[ -n "${pid}" ]]; then
      pids+=( "${pid}" )
    fi
  done < <(jobs -pr)
  if [[ ${#pids[@]} -eq 0 ]]; then
    return 0
  fi

  local pid
  local status=0
  for pid in "${pids[@]}"; do
    if ! wait "${pid}"; then
      status=1
    fi
  done
  return "${status}"
}

check_species_cds_dir() {
  local dir_sp_cds=$1
  local species_cds_fasta=()
  local fasta_path
  while IFS= read -r fasta_path; do
    species_cds_fasta+=( "${fasta_path}" )
  done < <(gg_find_fasta_files "${dir_sp_cds}" 1)
  if [[ ${#species_cds_fasta[@]} -eq 0 ]]; then
    echo "No species_cds fasta files were found in: ${dir_sp_cds}"
    exit 1
  fi
  local error_log
  error_log=$(gg_mktemp)  # Temporary file to collect errors
  echo "$(date): Started validating the format of all species_cds files using ${GG_TASK_CPUS} processes."

  function check_single_species_cds () {
    local spfasta=$1
    local sp_ub
    local first_header
    local first_header_no_gt
    local first_header_sp
    local spfasta_startswith
    sp_ub=$(gg_species_name_from_path_or_dot "${spfasta}")
    local seq_names
    seq_names=$(seqkit seq "${spfasta}" | awk '/^>/ {print}')
    first_header=${seq_names%%$'\n'*}
    first_header=${first_header%%[[:space:]]*}
    first_header_no_gt=${first_header#>}
    first_header_sp=$(gg_species_name_from_path_or_dot "${first_header_no_gt}")
    spfasta_startswith=">${first_header_sp}"

    if [[ ${spfasta_startswith} != ">${sp_ub}" ]]; then
      echo "Sequence names start with ${spfasta_startswith} but this is not consistent with the species name (${sp_ub}) parsed from the file name of ${spfasta}" >> "${error_log}"
    fi

    local num_all_seq
    local num_uniq_seq
    num_all_seq=$(printf '%s\n' "${seq_names}" | grep -c '^>')
    num_uniq_seq=$(printf '%s\n' "${seq_names}" | sort -u | grep -c '^>')
    if [[ ${num_all_seq} -ne ${num_uniq_seq} ]]; then
      echo "Sequence names are not unique. # all seqs = ${num_all_seq} and # unique seqs = ${num_uniq_seq}" >> "${error_log}"
    fi

    for prohibited_character in "%" "/" "+" ":" ";" "&" "^" "$" "#" "@" "!" "~" "=" "\'" "\"" "\`" "*" "(" ")" "{" "}" "[" "]" "|" "?" " " "\t"; do
      if [[ "${seq_names}" == *"${prohibited_character}"* ]]; then
        echo "Sequence names contain '${prohibited_character}': ${spfasta}" >> "${error_log}"
      fi
    done
  }

  export -f check_single_species_cds
  export -f gg_species_name_from_path_or_dot _gg_strip_species_terminal_suffixes _gg_species_prefix_token_count _gg_species_rank_token_key _gg_species_is_rank_or_qualifier_token _gg_species_label_prefix_part
  export error_log
  if command -v parallel >/dev/null 2>&1; then
    parallel --jobs "${GG_TASK_CPUS}" check_single_species_cds ::: "${species_cds_fasta[@]}"
  else
    for spfasta in "${species_cds_fasta[@]}"; do
      check_single_species_cds "${spfasta}"
    done
  fi

  if [[ -s "${error_log}" ]]; then
    cat "${error_log}"
    rm -f -- "${error_log}"  # Clean up the temporary file
    echo "Exiting due to errors."
    exit 1
  else
    rm -f -- "${error_log}"  # Clean up the temporary file
    echo "$(date): All per-species CDS files are valid."
  fi
}

check_species_cds() {
  local gg_workspace_dir=$1
  check_species_cds_dir "$(workspace_input_root "${gg_workspace_dir}")/species_cds"
}

check_species_protein_dir() {
  local dir_sp_protein=$1
  local species_protein_fasta=()
  local fasta_path
  while IFS= read -r fasta_path; do
    species_protein_fasta+=( "${fasta_path}" )
  done < <(gg_find_fasta_files "${dir_sp_protein}" 1)
  if [[ ${#species_protein_fasta[@]} -eq 0 ]]; then
    echo "No species_protein fasta files were found in: ${dir_sp_protein}"
    exit 1
  fi
  local error_log
  error_log=$(gg_mktemp)
  echo "$(date): Started validating the format of all species_protein files using ${GG_TASK_CPUS} processes."

  function check_single_species_protein () {
    local spfasta=$1
    local sp_ub
    local first_header
    local first_header_no_gt
    local first_header_sp
    local spfasta_startswith
    local seq_names
    sp_ub=$(gg_species_name_from_path_or_dot "${spfasta}")
    seq_names=$(seqkit seq "${spfasta}" | awk '/^>/ {print}')
    first_header=${seq_names%%$'\n'*}
    first_header=${first_header%%[[:space:]]*}
    first_header_no_gt=${first_header#>}
    first_header_sp=$(gg_species_name_from_path_or_dot "${first_header_no_gt}")
    spfasta_startswith=">${first_header_sp}"

    if [[ ${spfasta_startswith} != ">${sp_ub}" ]]; then
      echo "Sequence names start with ${spfasta_startswith} but this is not consistent with the species name (${sp_ub}) parsed from the file name of ${spfasta}" >> "${error_log}"
    fi

    local num_all_seq
    local num_uniq_seq
    num_all_seq=$(printf '%s\n' "${seq_names}" | grep -c '^>')
    num_uniq_seq=$(printf '%s\n' "${seq_names}" | sort -u | grep -c '^>')
    if [[ ${num_all_seq} -ne ${num_uniq_seq} ]]; then
      echo "Sequence names are not unique. # all seqs = ${num_all_seq} and # unique seqs = ${num_uniq_seq}" >> "${error_log}"
    fi

    for prohibited_character in "%" "/" "+" ":" ";" "&" "^" "$" "#" "@" "!" "~" "=" "\'" "\"" "\`" "*" "(" ")" "{" "}" "[" "]" "|" "?" " " "\t"; do
      if [[ "${seq_names}" == *"${prohibited_character}"* ]]; then
        echo "Sequence names contain '${prohibited_character}': ${spfasta}" >> "${error_log}"
      fi
    done
  }

  export -f check_single_species_protein
  export -f gg_species_name_from_path_or_dot _gg_strip_species_terminal_suffixes _gg_species_prefix_token_count _gg_species_rank_token_key _gg_species_is_rank_or_qualifier_token _gg_species_label_prefix_part
  export error_log
  if command -v parallel >/dev/null 2>&1; then
    parallel --jobs "${GG_TASK_CPUS}" check_single_species_protein ::: "${species_protein_fasta[@]}"
  else
    for spfasta in "${species_protein_fasta[@]}"; do
      check_single_species_protein "${spfasta}"
    done
  fi

  if [[ -s "${error_log}" ]]; then
    cat "${error_log}"
    rm -f -- "${error_log}"
    echo "Exiting due to errors."
    exit 1
  else
    rm -f -- "${error_log}"
    echo "$(date): All per-species protein files are valid."
  fi
}

check_species_protein() {
  local gg_workspace_dir=$1
  check_species_protein_dir "$(workspace_input_root "${gg_workspace_dir}")/species_protein"
}

is_output_older_than_inputs() {
  local input_file_variable_regex=$1
  local output_file=$2
  local return_flag=0
  local i=0
  if [[ ! -s "${output_file}" ]]; then
    echo "Output file not found. Will be generated: ${output_file}"
    return_flag=1
	else
	  echo "Checking whether there are any input files that are newer than the output file: ${output_file}"
	  local infiles=()
	  local input_var_name
	  while IFS= read -r input_var_name; do
	    infiles+=( "${input_var_name}" )
	  done < <(compgen -A variable | grep -E -- "${input_file_variable_regex}" || true)
	  for file_var in "${infiles[@]}"; do
	    local file_path="${!file_var}"
      if [[ -e "${file_path}" ]]; then
        i=$((i+1))
        if [[ "${file_path}" -nt "${output_file}" ]]; then
          echo "Output file will be renewed. Detected a new input file: ${file_path}"
          return_flag=1
        fi
      fi
    done
  fi
  if [[ ${return_flag} -eq 0 ]]; then
    echo "All examined input files (${i}) are older than the output file: ${output_file}"
  fi
  return ${return_flag}
}

stop_if_species_not_found_in() {
  local dir=$1
  local species_name=$2
  species_name=$(echo "${species_name}" | sed -e "s/[[:space:]]/_/g")
  local files=()
  local file
  local file_species=""
  shopt -s nullglob dotglob
  for file in "${dir}"/*; do
    [[ -f "${file}" ]] || continue
    [[ "$(basename "${file}")" == .* ]] && continue
    files+=( "${file}" )
  done
  shopt -u nullglob dotglob
  local species_found_flag=0
  for file in "${files[@]}"; do
    file_species=$(gg_species_name_from_path_or_dot "$(basename "${file}")")
    if [[ "${file_species}" == "${species_name}" ]]; then
      species_found_flag=1
      break
    fi
  done
  if [[ ${species_found_flag} -eq 0 ]]; then
    echo "Exiting. The input file for ${species_name} not found in: ${dir}"
    exit 1
  else
    echo "A file for ${species_name} found in: ${dir}"
  fi
}

print_softmasked_percentage() {
  local fasta_path=$1
  local num_total_bp
  local num_masked_bp
  num_total_bp=$(grep -v '^>' "${fasta_path}" | tr -d '\n' | wc -c | sed -e "s/[[:space:]]//g")
  num_masked_bp=$(grep -v '^>' "${fasta_path}" | tr -d '\n' | tr -d -c 'atgc' | wc -c | sed -e "s/[[:space:]]//g")
  if [[ "${num_total_bp}" -eq 0 ]]; then
    echo "0.0% masked (0/0 bp)"
    return 0
  fi
  python -c 'import sys; num=int(sys.argv[1]); den=int(sys.argv[2]); print("{:,.1f}% masked ({:,}/{:,} bp)".format(num/den*100, num, den))' "${num_masked_bp}" "${num_total_bp}"
}

gg_normalize_input_sequence_mode() {
  local input_mode="${1:-cds}"
  input_mode=$(printf '%s' "${input_mode}" | tr '[:upper:]' '[:lower:]')
  if [[ "${input_mode}" != "cds" && "${input_mode}" != "protein" ]]; then
    echo "Invalid input_sequence_mode: ${input_mode}" >&2
    echo 'input_sequence_mode must be either "cds" or "protein".' >&2
    return 1
  fi
  printf '%s\n' "${input_mode}"
}

gg_species_genetic_code_table_path() {
  local input_root=${1:-}
  echo "${input_root}/species_genetic_code/species_genetic_code.tsv"
}

gg_species_protein_input_dir_path() {
  local input_root=${1:-}
  echo "${input_root}/species_protein"
}

gg_prepare_species_genetic_code_table() {
  local cds_dir=$1
  local default_code=$2
  local outfile=$3
  local input_table=${4:-}
  python - "${cds_dir}" "${default_code}" "${outfile}" "${input_table}" <<'PY'
import csv
import sys
from pathlib import Path


def species_from_filename(path: Path) -> str:
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
    name = path.name
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
input_table_raw = sys.argv[4].strip()
input_table = Path(input_table_raw) if input_table_raw else None
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
if input_table is not None and input_table.exists():
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

gg_lookup_species_genetic_code() {
  local species_name=$1
  local table_path=$2
  local default_code=${3:-1}
  local code=""
  if [[ -s "${table_path}" ]]; then
    code=$(awk -F'\t' -v sp="${species_name}" 'NR>1 && $1==sp {print $2; exit}' "${table_path}")
  fi
  if [[ -z "${code}" ]]; then
    code="${default_code}"
  fi
  printf '%s' "${code}"
}
