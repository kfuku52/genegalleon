#!/usr/bin/env bash

unset_singularity_envs() {
	unset SINGULARITY_BIND
	unset SINGULARITY_BINDPATH
	unset SINGULARITYENV_SGE_TASK_ID
	unset SINGULARITYENV_NSLOTS
	unset SINGULARITYENV_JOB_ID
	unset SINGULARITYENV_MEM_PER_SLOT
	unset SINGULARITYENV_MEM_PER_HOST
}

variable_SGEnizer() {
	local echo_header='variable_SGEnizer: '
	echo ${echo_header}'SLURM environmental variables are converted to the SGE style.'
	if [[ -z "${NSLOTS:-}" ]]; then
		echo ${echo_header}'${NSLOTS} is undefined or empty.'
		if [[ -n "${SLURM_CPUS_PER_TASK:-}" ]]; then
			echo ${echo_header}'NSLOTS=${SLURM_CPUS_PER_TASK}'
			NSLOTS=${SLURM_CPUS_PER_TASK}
		elif [[ -n "${PBS_NODEFILE:-}" && -f "${PBS_NODEFILE}" ]]; then
			PBS_NSLOTS=$(wc -l < "${PBS_NODEFILE}")
			echo ${echo_header}'NSLOTS=${PBS_NSLOTS}'
			NSLOTS=${PBS_NSLOTS}
		else
			echo ${echo_header}'both ${NSLOTS}, ${PBS_NODEFILE}, and ${SLURM_NPROCS} were undefined. NSLOTS=1'
			NSLOTS=1
		fi
	else
		echo ${echo_header}'${NSLOTS} is defined.'
	fi
	if [[ -z "${JOB_ID+x}" ]]; then
		echo ${echo_header}'${JOB_ID} is undefined.'
		if [[ -n "${SLURM_JOB_ID+x}" ]]; then
			echo ${echo_header}'JOB_ID=${SLURM_JOB_ID}'
			JOB_ID=${SLURM_JOB_ID}
		elif [[ -n "${PBS_JOBID+x}" ]]; then
			echo ${echo_header}'JOB_ID=${PBS_JOBID}'
			JOB_ID=${PBS_JOBID}
		else
			echo ${echo_header}'both ${JOB_ID} and ${SLURM_JOB_ID} were undefined. JOB_ID=1'
			JOB_ID=1
		fi
	else
		echo ${echo_header}'${JOB_ID} is defined.'
	fi
	if [[ -z "${SGE_TASK_ID+x}" ]]; then
		echo ${echo_header}'${SGE_TASK_ID} is undefined.'
		if [[ -n "${SLURM_ARRAY_TASK_ID+x}" ]]; then
			echo ${echo_header}'SGE_TASK_ID=${SLURM_ARRAY_TASK_ID}'
			SGE_TASK_ID=${SLURM_ARRAY_TASK_ID}
		elif [[ -n "${PBS_ARRAY_INDEX+x}" ]]; then
			echo ${echo_header}'SGE_TASK_ID=${PBS_ARRAY_INDEX}'
			SGE_TASK_ID=${PBS_ARRAY_INDEX}
		else
			echo ${echo_header}'both ${SGE_TASK_ID} and ${SLURM_ARRAY_TASK_ID} were undefined. SGE_TASK_ID=1'
			SGE_TASK_ID=1
		fi
	else
		echo ${echo_header}'${SGE_TASK_ID} is defined.'
	fi
	if type qstat >/dev/null 2>&1; then
		MEM_PER_SLOT=$(qstat -f -j "${JOB_ID}" 2>/dev/null | grep ",mem_req=" | sed -e "s/.*mem_req=//" -e "s/G,.*//" | sed -e "s/G.*//")
	fi
	if [[ -z "${MEM_PER_SLOT:-}" ]]; then
		echo ${echo_header}'${MEM_PER_SLOT} is undefined. '
		if [[ -n "${SLURM_MEM_PER_CPU+x}" ]]; then
			echo ${echo_header}'MEM_PER_SLOT=$((${SLURM_MEM_PER_CPU}/1024))'
			MEM_PER_SLOT=$((${SLURM_MEM_PER_CPU}/1024))
		else
			echo ${echo_header}'both ${MEM_PER_SLOT} and ${SLURM_MEM_PER_CPU} were undefined. MEM_PER_SLOT=3'
			MEM_PER_SLOT=3
		fi
	fi
	MEM_PER_HOST=$((${MEM_PER_SLOT}*${NSLOTS}))
	echo ${echo_header}"NSLOTS=${NSLOTS}"
	echo ${echo_header}"JOB_ID=${JOB_ID}"
	echo ${echo_header}"SGE_TASK_ID=${SGE_TASK_ID}"
	echo ${echo_header}"MEM_PER_SLOT=${MEM_PER_SLOT}"
	echo ${echo_header}"MEM_PER_HOST=${MEM_PER_HOST}"
	echo ""
  export NSLOTS
  export JOB_ID
  export SGE_TASK_ID
  export MEM_PER_SLOT
  export MEM_PER_HOST
}

set_singularityenv() {
  local echo_header="set_singularityenv: "
  echo ${echo_header}"original: dir_pg = ${dir_pg}"
  echo ${echo_header}"original: dir_script = ${dir_script}"
  echo ${echo_header}"original: gg_image = ${gg_image}"
  if [[ $(uname -s) != 'Darwin' ]]; then
    echo "OS is $(uname -s). Getting the original path of symlink."
    dir_pg=$(cd "$(dirname "${dir_pg}")" && pwd -P)/$(basename "${dir_pg}")
    dir_script=$(cd "$(dirname "${dir_script}")" && pwd -P)/$(basename "${dir_script}")
    gg_image=$(cd "$(dirname "${gg_image}")" && pwd -P)/$(basename "${gg_image}")
    echo ${echo_header}"formatted: dir_pg = ${dir_pg}"
    echo ${echo_header}"formatted: dir_script = ${dir_script}"
    echo ${echo_header}"formatted: gg_image = ${gg_image}"
	else
		echo "OS is $(uname -s). Symlink PATHs won't be updated."
	fi
	SINGULARITY_BIND="${dir_pg}:/workspace,${SINGULARITY_BIND}"
	SINGULARITY_BIND="${dir_script}:/script,${SINGULARITY_BIND}"
	SINGULARITY_BINDPATH="${dir_pg}:/workspace,${SINGULARITY_BINDPATH}" # For Julia SLURM with Singularity 2.5
	SINGULARITY_BINDPATH="${dir_script}:/script,${SINGULARITY_BINDPATH}" # For Julia SLURM with Singularity 2.5
	export SINGULARITY_BIND
	export SINGULARITY_BINDPATH # For Julia SLURM with Singularity 2.5
	export SINGULARITYENV_SGE_TASK_ID=${SGE_TASK_ID}
	export SINGULARITYENV_NSLOTS=${NSLOTS}
	export SINGULARITYENV_JOB_ID=${JOB_ID}
	export SINGULARITYENV_MEM_PER_SLOT=${MEM_PER_SLOT}
	export SINGULARITYENV_MEM_PER_HOST=${MEM_PER_HOST}
  if [[ -n ${delete_preexisting_tmp_dir} ]]; then
    export SINGULARITYENV_delete_preexisting_tmp_dir=${delete_preexisting_tmp_dir}
  else
    export SINGULARITYENV_delete_preexisting_tmp_dir=0
  fi
	if [[ -n ${delete_tmp_dir} ]]; then
		export SINGULARITYENV_delete_tmp_dir=${delete_tmp_dir}
	else
		export SINGULARITYENV_delete_tmp_dir=0
	fi
	set | grep "^SINGULARITY"
	echo ""
}

ensure_dir() {
	local d=$1
	if [[ -z "${d}" ]]; then
		return 0
	fi
	if [[ ! -d "${d}" ]]; then
		mkdir -p "${d}"
	fi
}

ensure_parent_dir() {
	local p=$1
	if [[ -z "${p}" ]]; then
		return 0
	fi
	ensure_dir "$(dirname -- "${p}")"
}

gg_species_name_from_path() {
  local path=$1
  local basename_path
  basename_path=$(basename "${path}")
  echo "${basename_path}" | sed -e "s/_/|/" -e "s/_.*//" -e "s/|/_/"
}

gg_species_name_from_path_or_dot() {
  local path=$1
  local basename_path
  basename_path=$(basename "${path}")
  echo "${basename_path}" | sed -e "s/_/|/" -e "s/[_\\.].*//" -e "s/|/_/"
}

gg_fasta_relabel_headers_to_species() {
  awk '
    /^>/ {
      header = substr($0, 2)
      sub("^.*/", "", header)
      n = split(header, a, "_")
      if (n >= 2) {
        print ">" a[1] "_" a[2]
      } else {
        print ">" header
      }
      next
    }
    { print }
  '
}

gg_busco_gene_tokens() {
  local mode=$1
  shift
  local token
  local split_token

  case "${mode}" in
    exclude_duplicated)
      for token in "$@"; do
        if [[ "${token}" == "-" || "${token}" == *","* ]]; then
          continue
        fi
        echo "${token}"
      done
      ;;
    split_duplicated)
      for token in "$@"; do
        if [[ "${token}" == "-" ]]; then
          continue
        fi
        for split_token in ${token//,/ }; do
          if [[ -n "${split_token}" ]]; then
            echo "${split_token}"
          fi
        done
      done
      ;;
    *)
      echo "gg_busco_gene_tokens: invalid mode: ${mode}" >&2
      return 1
      ;;
  esac
}

gg_seqkit_grep_by_patterns_from_infile_list() {
  local threads=$1
  local infile_list=$2
  shift 2
  local patterns=( "$@" )
  seqkit grep --threads "${threads}" --pattern-file <(printf '%s\n' "${patterns[@]}") --infile-list "${infile_list}"
}

gg_prepare_cds_fasta_stream() {
  local threads=${1:-1}
  local codontable=${2:-}
  seqkit replace --pattern X --replacement N --by-seq --ignore-case --threads "${threads}" \
  | seqkit replace --pattern " .*" --replacement "" --ignore-case --threads "${threads}"
  if [[ -n "${codontable}" ]]; then
    cdskit pad --codontable "${codontable}"
  else
    cdskit pad
  fi
}

workspace_input_root() {
  local dir_pg=$1
  if [[ -d "${dir_pg}/input" || -d "${dir_pg}/output" ]]; then
    echo "${dir_pg}/input"
  else
    echo "${dir_pg}"
  fi
}

workspace_output_root() {
  local dir_pg=$1
  if [[ -d "${dir_pg}/input" || -d "${dir_pg}/output" ]]; then
    echo "${dir_pg}/output"
  else
    echo "${dir_pg}"
  fi
}

workspace_db_root() {
  local dir_pg=$1
  local preferred_root="${dir_pg}/db"
  local legacy_split_root="${dir_pg}/output/db"
  if [[ -d "${preferred_root}" || ! -d "${legacy_split_root}" ]]; then
    echo "${preferred_root}"
  else
    echo "${legacy_split_root}"
  fi
}

workspace_taxonomy_root() {
  local dir_pg=$1
  local dir_db
  dir_db=$(workspace_db_root "${dir_pg}")
  echo "${dir_db}/ete_taxonomy"
}

workspace_pfam_root() {
  local dir_pg=$1
  local dir_db
  dir_db=$(workspace_db_root "${dir_pg}")
  echo "${dir_db}/pfam"
}

workspace_pfam_le_dir() {
  local dir_pg=$1
  local dir_pfam
  dir_pfam=$(workspace_pfam_root "${dir_pg}")
  echo "${dir_pfam}/Pfam_LE"
}

workspace_taxonomy_dbfile() {
  local dir_pg=$1
  local dir_taxonomy
  dir_taxonomy=$(workspace_taxonomy_root "${dir_pg}")
  echo "${dir_taxonomy}/taxa.sqlite"
}

gg_set_taxonomy_cache_env() {
  local dir_pg=$1
  local dir_taxonomy
  dir_taxonomy=$(workspace_taxonomy_root "${dir_pg}")
  ensure_dir "${dir_taxonomy}"
  ensure_dir "${dir_taxonomy}/ete"
  export ETE_DATA_HOME="${dir_taxonomy}"
  export ETE_CONFIG_HOME="${dir_taxonomy}"
  export XDG_DATA_HOME="${dir_taxonomy}"
  export XDG_CONFIG_HOME="${dir_taxonomy}"
  export GG_TAXONOMY_DBFILE="${dir_taxonomy}/taxa.sqlite"
}

_ensure_ete_taxonomy_db_locked() {
  local db_file=$1
  local py_exec=""
  local candidate

  for candidate in python python3 /opt/conda/bin/python /usr/bin/python3; do
    if [[ -x "${candidate}" ]]; then
      py_exec="${candidate}"
      break
    fi
    if command -v "${candidate}" >/dev/null 2>&1; then
      py_exec="$(command -v "${candidate}")"
      break
    fi
  done

  if [[ -z "${py_exec}" ]]; then
    echo "python/python3 command was not found. Cannot prepare ETE taxonomy DB." >&2
    return 1
  fi

  ETE_DATA_HOME="$(dirname "${db_file}")" \
  ETE_CONFIG_HOME="$(dirname "${db_file}")" \
  XDG_DATA_HOME="$(dirname "${db_file}")" \
  XDG_CONFIG_HOME="$(dirname "${db_file}")" \
  GG_TAXONOMY_DBFILE="${db_file}" \
  "${py_exec}" - <<'PY'
import importlib
import os

db_file = os.environ.get("GG_TAXONOMY_DBFILE", "").strip()
if not db_file:
    raise SystemExit("GG_TAXONOMY_DBFILE is empty.")

os.makedirs(os.path.dirname(db_file), exist_ok=True)

candidate_cache_dirs = []
for key in ("ETE_DATA_HOME", "ETE_CONFIG_HOME", "XDG_DATA_HOME", "XDG_CONFIG_HOME"):
    value = os.environ.get(key, "").strip()
    if value:
        candidate_cache_dirs.append(value)
        candidate_cache_dirs.append(os.path.join(value, "ete"))

home = os.path.expanduser("~")
candidate_cache_dirs.append(os.path.join(home, ".local", "share", "ete"))
candidate_cache_dirs.append(os.path.join(home, ".config", "ete"))

for cache_dir in dict.fromkeys(candidate_cache_dirs):
    os.makedirs(cache_dir, exist_ok=True)

def ensure_with_ete4():
    ncbiquery = importlib.import_module("ete4.ncbi_taxonomy.ncbiquery")
    if os.path.exists(db_file) and ncbiquery.is_taxadb_up_to_date(db_file):
        return "ete4:up_to_date"
    NCBITaxa = importlib.import_module("ete4").NCBITaxa
    NCBITaxa(dbfile=db_file, update=True)
    return "ete4:updated"

try:
    ensure_with_ete4()
except Exception as exc:  # pragma: no cover - runtime dependency handling
    raise SystemExit(f"Failed to prepare ETE taxonomy DB with ete4: {exc}")

if not os.path.exists(db_file) or os.path.getsize(db_file) == 0:
    raise SystemExit(f"ETE taxonomy DB was not generated: {db_file}")
PY

  if [[ $? -ne 0 || ! -s "${db_file}" ]]; then
    echo "Failed to prepare ETE taxonomy DB: ${db_file}" >&2
    return 1
  fi
}

ensure_ete_taxonomy_db() {
  local dir_pg=$1
  local dir_db
  local dir_taxonomy
  local db_file
  local lock_file
  local has_flock=1

  dir_db=$(workspace_db_root "${dir_pg}")
  dir_taxonomy=$(workspace_taxonomy_root "${dir_pg}")
  db_file=$(workspace_taxonomy_dbfile "${dir_pg}")
  lock_file="${dir_db}/locks/ete_taxonomy.lock"

  ensure_dir "${dir_db}"
  ensure_dir "${dir_taxonomy}"
  ensure_dir "${dir_taxonomy}/ete"
  ensure_dir "$(dirname "${lock_file}")"
  export ETE_DATA_HOME="${dir_taxonomy}"
  export ETE_CONFIG_HOME="${dir_taxonomy}"
  export XDG_DATA_HOME="${dir_taxonomy}"
  export XDG_CONFIG_HOME="${dir_taxonomy}"
  export GG_TAXONOMY_DBFILE="${db_file}"

  if ! command -v flock >/dev/null 2>&1; then
    has_flock=0
    echo "flock command was not found. Proceeding without lock file synchronization: ETE taxonomy DB" >&2
  fi

  if [[ ${has_flock} -eq 1 ]]; then
    exec 8> "${lock_file}"
    flock 8
  fi

  _ensure_ete_taxonomy_db_locked "${db_file}"
  local ensure_exit_code=$?

  if [[ ${has_flock} -eq 1 ]]; then
    flock -u 8
    exec 8>&-
  fi
  return ${ensure_exit_code}
}

gg_initialize_data_layout() {
  local dir_pg=$1
  local dir_input
  local dir_output
  local dir_db
  dir_input=$(workspace_input_root "${dir_pg}")
  dir_output=$(workspace_output_root "${dir_pg}")
  dir_db=$(workspace_db_root "${dir_pg}")
  if [[ "${dir_input}" != "${dir_pg}" || "${dir_output}" != "${dir_pg}" ]]; then
    ensure_dir "${dir_input}"
    ensure_dir "${dir_output}"
  fi
  ensure_dir "${dir_db}"
  gg_set_taxonomy_cache_env "${dir_pg}"
}

cp_out() {
	if [[ $# -lt 2 ]]; then
		echo "cp_out: at least 2 arguments are required."
		return 1
	fi
	local dest="${!#}"
	ensure_parent_dir "${dest}"
	cp "$@"
}

mv_out() {
	if [[ $# -lt 2 ]]; then
		echo "mv_out: at least 2 arguments are required."
		return 1
	fi
	local dest="${!#}"
	ensure_parent_dir "${dest}"
	mv "$@"
}

remove_empty_subdirs() {
  local dir_main=$1
  local echo_header="remove_empty_subdirs: "
	if [[ ! -d "${dir_main}" ]]; then
		echo "${echo_header}directory not found: ${dir_main}"
		echo ""
		return
	fi
	local sub_directories=()
	mapfile -d '' -t sub_directories < <(find "${dir_main}" -mindepth 1 -maxdepth 1 -type d -print0)
	for d in "${sub_directories[@]}"; do
		if [[ -z "$(find "${d}" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]]; then
			echo "${echo_header}deleting ${d}"
			rm -r "${d}"
		fi
	done
	echo ""
}

set_singularity_command() {
  local echo_header="set_singularity_command: "
  echo ${echo_header}"hostname = $(hostname)"
  if [[ $(hostname) == at* || $(hostname) == m* || $(hostname) == igt* || $(hostname) == it* ]]; then
    echo ${echo_header}"NIG Supercomputer: loading singularity 3"
    #module load singularity/3.2.0
    if [[ -e /var/spool/uge ]]; then
      export SINGULARITY_BIND="/var/spool/uge:/var/spool/uge,${SINGULARITY_BIND}"
    fi
    if [[ -e /var/spool/age ]]; then
      export SINGULARITY_BIND="/var/spool/age:/var/spool/age,${SINGULARITY_BIND}"
    fi
    export SINGULARITY_BIND="/opt/pkg:/opt/pkg,${SINGULARITY_BIND}"
    export SINGULARITY_BIND="/home/geadmin/UGER/uger/spool:/home/geadmin/UGER/uger/spool,${SINGULARITY_BIND}"
    #export SINGULARITY_BIND="/usr/local/seq/blast/ncbi:/usr/local/seq/blast/ncbi,${SINGULARITY_BIND}"
    #singularity_command="singularity shell -B /var/spool/uge:/var/spool/uge -B /opt/pkg:/opt/pkg"
    singularity_command="singularity shell"
  elif [[ $(hostname) == *nhr.fau.de ]]; then
  	echo ${echo_header}"NHR FAU Fritz"
  	singularity_command="singularity shell --contain"
  else
    echo ${echo_header}"No special process is activated in this host."
    singularity_command="singularity shell"
  fi
  echo ${echo_header}'${singularity_command}' = \"${singularity_command}\"
  echo ""
}

exit_if_running_qstat() {
	local flag_start=0
	if [[ "${exit_if_running:-0}" -eq 1 ]]; then
		flag_start=1
	fi
	if qstat > /dev/null 2>&1; then
		echo "qstat was found."
		flag_start=$((flag_start*1))
	else
		echo "qstat was not found. exit_if_running is deactivated."
		flag_start=0
	fi
	if [[ ${flag_start} -eq 1 ]]; then
		local job_name
		job_name=$(qstat | grep "${JOB_ID}" | head -n 1 | sed -e "s/^[[:space:]]*//" | cut -d " " -f 3)
		local running_id=()
		mapfile -t running_id < <(qstat | grep "${job_name}" | grep -v -e " dr " -e "${JOB_ID} .* ${SGE_TASK_ID}$" | sed -e "s/ qw / qw dummy/" -e "s/\s\+/\t/g" -e "1,2d" | cut -f11 | grep -v ":" | sort -u)
		echo "running_id: ${running_id[*]}"
		local flag=1
		for r in "${running_id[@]}"; do
			if [[ "${r}" =~ ^[0-9]+$ ]] && [[ "${r}" -eq "${SGE_TASK_ID}" ]]; then
				flag=0
			fi
		done
		if [[ ${flag} -eq 0 ]]; then
			echo "SGE_TASK_ID=${SGE_TASK_ID} is running already. Exiting."
			exit 0
		fi
	fi
}

wait_until_jobn_le() {
  # https://qiita.com/jokester/items/34ce222e1702c6e120eb
  local max_jobn=$1
  while [[ "$(jobs | wc -l)" -gt "$max_jobn" ]] ; do
    sleep 1
  done
}

check_species_cds() {
  local dir_pg=$1
  local dir_sp_cds="$(workspace_input_root "${dir_pg}")/species_cds"
  local species_cds_fasta=()
  mapfile -t species_cds_fasta < <(find "${dir_sp_cds}" -maxdepth 1 -type f | grep $(get_fasta_extensions_for_grep))
  if [[ ${#species_cds_fasta[@]} -eq 0 ]]; then
    echo "No species_cds fasta files were found in: ${dir_sp_cds}"
    exit 1
  fi
  local error_log=$(mktemp)  # Temporary file to collect errors
  echo "$(date): Started validating the format of all species_cds files using ${NSLOTS} processes."

  function check_single_species_cds () {
    local spfasta=$1
    local sp_ub
    local first_header
    local first_header_no_gt
    local first_header_sp
    local spfasta_startswith
    sp_ub=$(basename "${spfasta}" | sed -e "s/_/|/" -e "s/_.*//" -e "s/|/_/")
    local seq_names=$(seqkit seq "${spfasta}" | grep -e "^>")
    first_header=$(printf '%s\n' "${seq_names}" | head -n 1 | sed -e "s/[[:space:]].*//")
    first_header_no_gt=${first_header#>}
    first_header_sp=$(printf '%s' "${first_header_no_gt}" | sed -e "s/_/|/" -e "s/_.*//" -e "s/|/_/")
    spfasta_startswith=">${first_header_sp}"

    if [[ ${spfasta_startswith} != ">${sp_ub}" ]]; then
      echo "Sequence names start with ${spfasta_startswith} but this is not consistent with the species name (${sp_ub}) parsed from the file name of ${spfasta}" >> "${error_log}"
    fi

    local num_all_seq=$(seqkit seq "${spfasta}" | grep -e "^>" | wc -l)
    local num_uniq_seq=$(seqkit seq "${spfasta}" | grep -e "^>" | sort -u | wc -l)
    if [[ ${num_all_seq} -ne ${num_uniq_seq} ]]; then
      echo "Sequence names are not unique. # all seqs = ${num_all_seq} and # unique seqs = ${num_uniq_seq}" >> "${error_log}"
    fi

    for prohibited_character in "%" "/" "+" ":" ";" "&" "^" "$" "#" "@" "!" "~" "=" "\'" "\"" "\`" "*" "(" ")" "{" "}" "[" "]" "|" "?" " " "\t"; do
      printf '%s\n' "${seq_names}" | grep -e "^>" | grep -q -F -- "${prohibited_character}"
      local exit_code=$?
      if [[ ${exit_code} -eq 0 ]]; then
        echo "Sequence names contain '${prohibited_character}': ${spfasta}" >> "${error_log}"
      fi
    done
  }

  export -f check_single_species_cds
  export error_log
  parallel --jobs "${NSLOTS}" check_single_species_cds ::: "${species_cds_fasta[@]}"

  if [[ -s "${error_log}" ]]; then
    cat "${error_log}"
    rm -f "${error_log}"  # Clean up the temporary file
    echo "Exiting due to errors."
    exit 1
  else
    rm -f "${error_log}"  # Clean up the temporary file
    echo "$(date): All per-species CDS files are valid."
  fi
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
    mapfile -t infiles < <(set | grep "${input_file_variable_regex}" | sed -e "s/=.*//")
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
  shopt -s nullglob
  files=( "${dir}"/* )
  shopt -u nullglob
  local species_found_flag=0
  for file in "${files[@]}"; do
    if [[ "$(basename "${file}")" == *"${species_name}"* ]]; then
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
  local num_total_bp=$(grep -v '^>' "${fasta_path}" | wc -c | sed -e "s/[[:space:]]//g")
  local num_masked_bp=$(grep -v '^>' "${fasta_path}" | tr -d -c 'atgc'  | wc -c | sed -e "s/[[:space:]]//g")
  python -c 'import sys; num=int(sys.argv[1]); den=int(sys.argv[2]); print("{:,.1f}% masked ({:,}/{:,} bp)".format(num/den*100, num, den))' ${num_masked_bp} ${num_total_bp}
}

disable_if_no_input_file() {
  local run_variable_txt=$1
  shift
  local run_variable_value="${!run_variable_txt:-0}"
  local input_files=( "$@" )
  if [[ ${run_variable_value} -eq 0 ]]; then
    return
  fi
  local flag_deactivate=0
  for input_file in "${input_files[@]}"; do
    if [[ ! -s "${input_file}" ]]; then
      echo "Required input file undetected: ${input_file}"
      flag_deactivate=1
    fi
  done
  if [[ ${flag_deactivate} -eq 1 ]]; then
    printf -v "${run_variable_txt}" '%s' 0
    echo "Disabled ${run_variable_txt}"
  fi
}

check_if_species_files_unique() {
  local species_dir=$1
  local files=()
  shopt -s nullglob
  files=( "${species_dir}"/* )
  shopt -u nullglob
  local num_sp=${#files[@]}
  local sp_names=()
  local file_base
  for file in "${files[@]}"; do
    file_base=$(basename "${file}")
    sp_names+=( "$(gg_species_name_from_path "${file_base}")" )
  done
  local num_sp_uniq
  num_sp_uniq=$(printf '%s\n' "${sp_names[@]}" | sort -u | wc -l)
  echo "Number of species files and its scientific name unique: ${num_sp} and ${num_sp_uniq}"
  if [[ ${num_sp} -ne ${num_sp_uniq} ]]; then
    echo "Exiting. Species files are not unique in: ${species_dir}"
    echo "Species names: ${sp_names[@]}"
    exit 1
  else
    echo "Species files are unique in: ${species_dir}"
    echo "Species names: ${sp_names[@]}"
  fi
}

is_fastq_requiring_downstream_analysis_done() {
  local run_assembly_local=${run_assembly:-0}
  local run_rrna_local=${run_rRNA_contamination_report:-0}
  local run_merge_local=${run_amalgkit_merge:-0}
  local out=0

  if [[ ( ( ${run_assembly_local} -eq 1 && -s ${file_isoform} ) || ${run_assembly_local} -eq 0 ) && \
    ( ( ${run_rrna_local} -eq 1 && -s ${file_rRNA_contamination_report} ) || ${run_rrna_local} -eq 0 ) && \
    ( ( ${run_merge_local} -eq 1 && -s ${file_amalgkit_merge_count} ) || ${run_merge_local} -eq 0 ) ]]; then
    out=1
  fi

  if [[ ${run_assembly_local} -eq 0 && ${run_rrna_local} -eq 0 && ${run_merge_local} -eq 0 ]]; then
    out=0
  fi
  echo ${out}
}

get_total_fastq_len() {
  local input_dir=$1
  local regex=$2
  local files=()
  mapfile -d '' -t files < <(find "${input_dir}" -name "${regex}" -print0)
  if [[ ${#files[@]} -eq 0 ]]; then
    echo 0
    return 0
  fi

  local sum_len
  sum_len=$(seqkit stats --tabular "${files[@]}" \
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

gg_prepare_cmd_runtime() {
  local dir_pg_local=$1
  local conda_env=${2:-}
  local set_unlimited_stack=${3:-1}
  local print_start_message=${4:-1}

  gg_initialize_data_layout "${dir_pg_local}"
  dir_pg_input=$(workspace_input_root "${dir_pg_local}")
  dir_pg_output=$(workspace_output_root "${dir_pg_local}")
  dir_pg_db=$(workspace_db_root "${dir_pg_local}")
  export dir_pg_input dir_pg_output dir_pg_db

  if [[ -n "${conda_env}" ]]; then
    conda activate "${conda_env}"
  fi
  if [[ "${set_unlimited_stack}" -eq 1 ]]; then
    ulimit -s unlimited
  fi
  if [[ "${print_start_message}" -eq 1 ]]; then
    print_pg_container_starting_message
  fi
}

gg_source_home_bashrc() {
  if [[ -f /home/.bashrc ]]; then
    # shellcheck disable=SC1091
    source "/home/.bashrc"
  fi
}

gg_test_r_packages() {
  local rpackage
  for rpackage in "$@"; do
    echo "Testing: ${rpackage}" | tee >(cat >&2)
    R -q -e "suppressPackageStartupMessages(library(${rpackage}, quietly=TRUE))" > /dev/null
  done
}

gg_test_python_packages() {
  local pypackage
  for pypackage in "$@"; do
    echo "Testing: ${pypackage}" | tee >(cat >&2)
    python -c "import ${pypackage}" > /dev/null
  done
}

gg_test_shell_commands() {
  local command_text
  for command_text in "$@"; do
    echo "Testing: ${command_text}"
    eval "${command_text}" > /dev/null
  done
}

gg_print_spacer() {
  echo ""
  echo ""
}

gg_print_section() {
  local title=$1
  echo "### ${title} ###"
}

gg_list_or_not_found() {
  local target_path=$1
  if [[ -e "${target_path}" ]]; then
    ls -la "${target_path}"
  else
    echo "Not found: ${target_path}"
  fi
}

gg_print_labelled_path_listing() {
  local label=$1
  local target_path=$2
  echo "${label}"
  gg_list_or_not_found "${target_path}"
  gg_print_spacer
}

enable_all_run_flags_for_debug_mode() {
  local message=${1:-"gg debug mode: All run_* variables are forced to set 1."}
  local debug_mode=${gg_debug_mode:-0}
  local flags=()
  local flag

  if [[ "${debug_mode}" -ne 1 ]]; then
    return 1
  fi

  echo "${message}"
  mapfile -t flags < <(set | sed -n -E 's/^(run_[A-Za-z0-9_]+)=.*/\1/p')
  for flag in "${flags[@]}"; do
    printf -v "${flag}" '%s' 1
    export "${flag}"
  done
  return 0
}

ensure_scheduler_defaults() {
  local echo_header="ensure_scheduler_defaults: "
  if [[ -z "${NSLOTS:-}" ]]; then
    echo "${echo_header}NSLOTS is undefined or empty. NSLOTS=1"
    NSLOTS=1
  fi
  if [[ -z "${JOB_ID:-}" ]]; then
    echo "${echo_header}JOB_ID is undefined or empty. JOB_ID=1"
    JOB_ID=1
  fi
  if [[ -z "${SGE_TASK_ID:-}" ]]; then
    echo "${echo_header}SGE_TASK_ID is undefined or empty. SGE_TASK_ID=1"
    SGE_TASK_ID=1
  fi
  if [[ -z "${MEM_PER_SLOT:-}" ]]; then
    echo "${echo_header}MEM_PER_SLOT is undefined or empty. MEM_PER_SLOT=3"
    MEM_PER_SLOT=3
  fi
  if [[ -z "${MEM_PER_HOST:-}" ]]; then
    MEM_PER_HOST=$((MEM_PER_SLOT * NSLOTS))
    echo "${echo_header}MEM_PER_HOST is undefined or empty. MEM_PER_HOST=${MEM_PER_HOST}"
  fi
  export NSLOTS JOB_ID SGE_TASK_ID MEM_PER_SLOT MEM_PER_HOST
}

print_pg_container_starting_message() {
  local dir_pg_input_resolved
  local dir_pg_output_resolved
  local dir_pg_db_resolved
  ensure_scheduler_defaults
  dir_pg_input_resolved=$(workspace_input_root "${dir_pg}")
  dir_pg_output_resolved=$(workspace_output_root "${dir_pg}")
  dir_pg_db_resolved=$(workspace_db_root "${dir_pg}")
  echo "$(date): Starting gg Singularity/Apptainer environment"
  echo "pwd: $(pwd)"
  echo "dir_pg: ${dir_pg}"
  echo "dir_pg_input: ${dir_pg_input_resolved}"
  echo "dir_pg_output: ${dir_pg_output_resolved}"
  echo "dir_pg_db: ${dir_pg_db_resolved}"
  echo "python: $(which python)"
  echo "NSLOTS: ${NSLOTS}"
  echo "MEM_PER_SLOT: ${MEM_PER_SLOT}"
  echo "MEM_PER_HOST: ${MEM_PER_HOST}"
  echo "JOB_ID: ${JOB_ID}"
  echo "SGE_TASK_ID: ${SGE_TASK_ID}"
  echo "ulimit -Hn: $(ulimit -Hn)"
  echo "ulimit -Sn: $(ulimit -Sn)"
}

recreate_dir() {
  local dir=$1
  if [[ -e "${dir}" ]]; then
    echo "Deleting: ${dir}"
    rm -r "${dir}"
  fi
  echo "Creating: ${dir}"
  mkdir -p "${dir}"
}

gg_array_download_once() {
  local lock_file=$1
  local done_file=$2
  local description=$3
  shift 3
  local task_id=${SGE_TASK_ID:-1}
  local has_flock=1

  mkdir -p "$(dirname "${lock_file}")"

  if [[ -s "${done_file}" ]]; then
    return 0
  fi
  if ! command -v flock >/dev/null 2>&1; then
    has_flock=0
    echo "flock command was not found. Proceeding without lock file synchronization: ${description}" >&2
  fi
  if [[ ${has_flock} -eq 1 ]]; then
    exec 9> "${lock_file}"
    flock 9

    if [[ -s "${done_file}" ]]; then
      flock -u 9
      exec 9>&-
      return 0
    fi

    echo "SGE_TASK_ID=${task_id}: starting download: ${description}" >&2
    "$@"
    local download_exit_code=$?
    if [[ ${download_exit_code} -ne 0 ]]; then
      echo "SGE_TASK_ID=${task_id}: download failed: ${description}" >&2
      flock -u 9
      exec 9>&-
      return ${download_exit_code}
    fi
    if [[ ! -s "${done_file}" ]]; then
      echo "Downloaded DB file not found after synchronization: ${done_file}" >&2
      flock -u 9
      exec 9>&-
      return 1
    fi

    echo "SGE_TASK_ID=${task_id}: download completed: ${description}" >&2
    flock -u 9
    exec 9>&-
    return 0
  fi

  echo "SGE_TASK_ID=${task_id}: starting download (no flock): ${description}" >&2
  "$@"
  local download_exit_code=$?
  if [[ ${download_exit_code} -ne 0 ]]; then
    echo "SGE_TASK_ID=${task_id}: download failed (no flock): ${description}" >&2
    return ${download_exit_code}
  fi
  if [[ ! -s "${done_file}" ]]; then
    echo "Downloaded DB file not found after synchronization: ${done_file}" >&2
    return 1
  fi

  echo "SGE_TASK_ID=${task_id}: download completed (no flock): ${description}" >&2
}

gg_task_token_contains_sge_id() {
  local task_token=$1
  local target_task_id=$2
  local token_parts=()
  IFS=',' read -r -a token_parts <<< "${task_token}"
  local part
  for part in "${token_parts[@]}"; do
    if [[ "${part}" == "${target_task_id}" ]]; then
      return 0
    fi
    if [[ "${part}" =~ ^([0-9]+)-([0-9]+)(:([0-9]+))?$ ]]; then
      local start_id=${BASH_REMATCH[1]}
      local end_id=${BASH_REMATCH[2]}
      local step_id=${BASH_REMATCH[4]:-1}
      if (( target_task_id >= start_id && target_task_id <= end_id )); then
        if (( (target_task_id - start_id) % step_id == 0 )); then
          return 0
        fi
      fi
    fi
  done
  return 1
}

gg_is_sge_task_present_in_qstat() {
  local job_id=$1
  local target_task_id=$2
  if ! command -v qstat >/dev/null 2>&1; then
    return 2
  fi
  if [[ -z "${job_id}" ]]; then
    return 2
  fi
  local qstat_task_tokens=()
  mapfile -t qstat_task_tokens < <(qstat 2>/dev/null | awk -v jid="${job_id}" 'NR>2 && $1==jid {print $NF}')
  if [[ ${#qstat_task_tokens[@]} -eq 0 ]]; then
    return 1
  fi
  local token
  for token in "${qstat_task_tokens[@]}"; do
    if gg_task_token_contains_sge_id "${token}" "${target_task_id}"; then
      return 0
    fi
  done
  return 1
}

gg_get_slurm_array_master_job_id() {
  if [[ -n "${SLURM_ARRAY_JOB_ID:-}" ]]; then
    echo "${SLURM_ARRAY_JOB_ID}"
    return 0
  fi
  local id_candidate="${SLURM_JOB_ID:-${JOB_ID:-}}"
  if [[ -z "${id_candidate}" ]]; then
    return 0
  fi
  if [[ "${id_candidate}" =~ ^([0-9]+)_([0-9]+)$ ]]; then
    echo "${BASH_REMATCH[1]}"
  else
    echo "${id_candidate}"
  fi
}

gg_task_token_contains_slurm_id() {
  local job_token=$1
  local array_job_id=$2
  local target_task_id=$3
  local suffix
  local expr
  local expr_parts=()
  local part

  if [[ "${job_token}" == "${array_job_id}_${target_task_id}" ]]; then
    return 0
  fi
  if [[ "${job_token}" != "${array_job_id}_"* ]]; then
    return 1
  fi

  suffix="${job_token#${array_job_id}_}"
  if [[ "${suffix}" == \[*\] ]]; then
    expr="${suffix#[}"
    expr="${expr%]}"
  else
    expr="${suffix}"
  fi

  IFS=',' read -r -a expr_parts <<< "${expr}"
  for part in "${expr_parts[@]}"; do
    part="${part%%%*}" # remove SLURM throttling suffix, e.g. 1-100%10
    if [[ "${part}" == "${target_task_id}" ]]; then
      return 0
    fi
    if [[ "${part}" =~ ^([0-9]+)-([0-9]+)(:([0-9]+))?$ ]]; then
      local start_id=${BASH_REMATCH[1]}
      local end_id=${BASH_REMATCH[2]}
      local step_id=${BASH_REMATCH[4]:-1}
      if (( target_task_id >= start_id && target_task_id <= end_id )); then
        if (( (target_task_id - start_id) % step_id == 0 )); then
          return 0
        fi
      fi
    fi
  done
  return 1
}

gg_is_sge_task_present_in_squeue() {
  local array_job_id=$1
  local target_task_id=$2
  if ! command -v squeue >/dev/null 2>&1; then
    return 2
  fi
  if [[ -z "${array_job_id}" ]]; then
    return 2
  fi

  # Fast path: query the specific array task directly.
  if [[ -n "$(squeue -h -j "${array_job_id}_${target_task_id}" -o "%i" 2>/dev/null)" ]]; then
    return 0
  fi

  local slurm_job_tokens=()
  mapfile -t slurm_job_tokens < <(squeue -h -j "${array_job_id}" -o "%i" 2>/dev/null)
  if [[ ${#slurm_job_tokens[@]} -eq 0 ]]; then
    return 1
  fi
  local token
  for token in "${slurm_job_tokens[@]}"; do
    if gg_task_token_contains_slurm_id "${token}" "${array_job_id}" "${target_task_id}"; then
      return 0
    fi
  done
  return 1
}

_download_uniprot_sprot_to_prefix() {
  local output_prefix=$1
  local output_dir
  output_dir=$(dirname "${output_prefix}")
  local tmp_dir
  tmp_dir=$(mktemp -d "${output_dir}/tmp.uniprot_sprot.XXXXXX")
  local pep_tmp="${tmp_dir}/uniprot_sprot.pep"
  local dmnd_tmp_prefix="${tmp_dir}/uniprot_sprot"
  local uniprot_url="${GG_UNIPROT_SPROT_URL:-https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz}"

  curl -fsSL "${uniprot_url}" | gzip -dc > "${pep_tmp}"
  if [[ ! -s "${pep_tmp}" ]]; then
    echo "Failed to download UniProt Swiss-Prot FASTA from: ${uniprot_url}" >&2
    rm -rf "${tmp_dir}"
    return 1
  fi

  diamond makedb --in "${pep_tmp}" --db "${dmnd_tmp_prefix}"
  if [[ ! -s "${dmnd_tmp_prefix}.dmnd" ]]; then
    echo "Failed to build DIAMOND DB from UniProt Swiss-Prot." >&2
    rm -rf "${tmp_dir}"
    return 1
  fi

  mv "${pep_tmp}" "${output_prefix}.pep"
  mv "${dmnd_tmp_prefix}.dmnd" "${output_prefix}.dmnd"
  rm -rf "${tmp_dir}"
}

ensure_uniprot_sprot_db() {
  local dir_pg=$1
  local sys_prefix="/usr/local/db/uniprot_sprot"
  local runtime_root="$(workspace_db_root "${dir_pg}")"
  local runtime_prefix="${runtime_root}/uniprot_sprot"
  local lock_file="${runtime_root}/locks/uniprot_sprot.lock"

  if [[ -s "${sys_prefix}.pep" && -s "${sys_prefix}.dmnd" ]]; then
    echo "${sys_prefix}"
    return 0
  fi
  if [[ -s "${runtime_prefix}.pep" && -s "${runtime_prefix}.dmnd" ]]; then
    echo "${runtime_prefix}"
    return 0
  fi

  mkdir -p "${runtime_root}"
  gg_array_download_once "${lock_file}" "${runtime_prefix}.dmnd" "UniProt Swiss-Prot (FASTA + DIAMOND)" \
    _download_uniprot_sprot_to_prefix "${runtime_prefix}" || return 1

  if [[ ! -s "${runtime_prefix}.pep" || ! -s "${runtime_prefix}.dmnd" ]]; then
    echo "UniProt DB download/build failed." >&2
    return 1
  fi
  echo "${runtime_prefix}"
}

_download_pfam_le_to_dir() {
  local output_dir=$1
  local parent_dir
  parent_dir=$(dirname "${output_dir}")
  local tmp_dir
  tmp_dir=$(mktemp -d "${parent_dir}/tmp.pfam_le.XXXXXX")
  local url="${GG_PFAM_LE_URL:-https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/little_endian/Pfam_LE.tar.gz}"
  local archive_path="${tmp_dir}/Pfam_LE.tar.gz"
  local staged_dir="${tmp_dir}/Pfam_LE"
  local pfam_files=()
  local f

  curl -fsSL "${url}" -o "${archive_path}"
  tar -xzf "${archive_path}" -C "${tmp_dir}"
  mkdir -p "${staged_dir}"
  mapfile -t pfam_files < <(find "${tmp_dir}" -type f \( -name "Pfam.*" -o -name "Pfam.pal" \))
  if [[ ${#pfam_files[@]} -eq 0 ]]; then
    echo "Pfam_LE archive did not contain Pfam DB files." >&2
    rm -rf "${tmp_dir}"
    return 1
  fi
  for f in "${pfam_files[@]}"; do
    cp "${f}" "${staged_dir}/"
  done

  # Old format: Pfam.loo/Pfam.rps, New split format: Pfam.pal + Pfam.00.* volumes.
  if [[ ! -s "${staged_dir}/Pfam.pal" && ( ! -s "${staged_dir}/Pfam.loo" || ! -s "${staged_dir}/Pfam.rps" ) ]]; then
    echo "Pfam_LE archive did not contain expected RPS-BLAST index files." >&2
    rm -rf "${tmp_dir}"
    return 1
  fi

  rm -rf "${output_dir}.tmp"
  mv "${staged_dir}" "${output_dir}.tmp"
  rm -rf "${output_dir}"
  mv "${output_dir}.tmp" "${output_dir}"
  rm -rf "${tmp_dir}"
}

pfam_le_db_is_ready() {
  local db_dir=$1
  [[ -s "${db_dir}/Pfam.pal" || ( -s "${db_dir}/Pfam.loo" && -s "${db_dir}/Pfam.rps" ) ]]
}

ensure_pfam_le_db() {
  local dir_pg=$1
  local sys_dir="/usr/local/db/Pfam_LE"
  local runtime_root="$(workspace_db_root "${dir_pg}")"
  local runtime_parent
  runtime_parent=$(workspace_pfam_root "${dir_pg}")
  local runtime_dir
  runtime_dir=$(workspace_pfam_le_dir "${dir_pg}")
  local legacy_runtime_dir="${runtime_root}/Pfam_LE"
  local lock_file="${runtime_root}/locks/pfam_le.lock"

  if pfam_le_db_is_ready "${sys_dir}"; then
    echo "${sys_dir}"
    return 0
  fi

  # Backward compatibility: migrate old workspace layout db/Pfam_LE -> db/pfam/Pfam_LE.
  if [[ ! -e "${runtime_dir}" ]] && pfam_le_db_is_ready "${legacy_runtime_dir}"; then
    mkdir -p "${runtime_parent}"
    if mv "${legacy_runtime_dir}" "${runtime_dir}" 2>/dev/null; then
      echo "Migrated legacy Pfam_LE directory to: ${runtime_dir}" >&2
    fi
  fi

  if pfam_le_db_is_ready "${runtime_dir}"; then
    echo "${runtime_dir}"
    return 0
  fi
  if pfam_le_db_is_ready "${legacy_runtime_dir}"; then
    echo "${legacy_runtime_dir}"
    return 0
  fi

  mkdir -p "${runtime_root}" "${runtime_parent}"
  gg_array_download_once "${lock_file}" "${runtime_dir}/Pfam.pal" "Pfam_LE RPS-BLAST DB" \
    _download_pfam_le_to_dir "${runtime_dir}" || return 1

  if ! pfam_le_db_is_ready "${runtime_dir}"; then
    echo "Pfam_LE DB download failed." >&2
    return 1
  fi
  echo "${runtime_dir}"
}

# Backward compatibility alias used by older workflow scripts.
ensure_pfam_domain_db() {
  ensure_pfam_le_db "$@"
}

_download_jaspar_meme_to_file() {
  local output_file=$1
  local output_dir
  output_dir=$(dirname "${output_file}")
  local filename
  filename=$(basename "${output_file}")
  local year
  year=$(echo "${filename}" | sed -n 's/^JASPAR\([0-9][0-9][0-9][0-9]\)_.*/\1/p')
  local tmp_file="${output_file}.tmp"
  local downloaded=0
  local url
  local urls=(
    "https://jaspar.elixir.no/download/data/${year}/CORE/${filename}"
    "https://jaspar.genereg.net/download/data/${year}/CORE/${filename}"
    "https://jaspar${year}.elixir.no/download/data/${year}/CORE/${filename}"
  )

  if [[ -z "${year}" ]]; then
    echo "Could not parse JASPAR year from file name: ${filename}" >&2
    return 1
  fi

  mkdir -p "${output_dir}"
  for url in "${urls[@]}"; do
    if curl -fsSL "${url}" -o "${tmp_file}"; then
      if [[ -s "${tmp_file}" ]]; then
        downloaded=1
        break
      fi
    fi
  done

  if [[ ${downloaded} -ne 1 ]]; then
    echo "Failed to download JASPAR file: ${filename}" >&2
    rm -f "${tmp_file}"
    return 1
  fi

  mv "${tmp_file}" "${output_file}"
}

ensure_jaspar_file() {
  local dir_pg=$1
  local jaspar_filename=$2
  local sys_file="/usr/local/db/jaspar/${jaspar_filename}"
  local runtime_root="$(workspace_db_root "${dir_pg}")"
  local runtime_file="${runtime_root}/jaspar/${jaspar_filename}"
  local lock_basename
  lock_basename=$(echo "${jaspar_filename}" | tr '/ ' '__')
  local lock_file="${runtime_root}/locks/jaspar_${lock_basename}.lock"

  if [[ -s "${sys_file}" ]]; then
    echo "${sys_file}"
    return 0
  fi
  if [[ -s "${runtime_file}" ]]; then
    echo "${runtime_file}"
    return 0
  fi

  mkdir -p "${runtime_root}"
  gg_array_download_once "${lock_file}" "${runtime_file}" "JASPAR motif file (${jaspar_filename})" \
    _download_jaspar_meme_to_file "${runtime_file}" || return 1

  if [[ ! -s "${runtime_file}" ]]; then
    echo "JASPAR download failed: ${jaspar_filename}" >&2
    return 1
  fi
  echo "${runtime_file}"
}

_download_silva_rrna_ref_to_file() {
  local output_file=$1
  local output_dir
  output_dir=$(dirname "${output_file}")
  local tmp_file="${output_file}.tmp"
  local url="${GG_SILVA_RRNA_URL:-https://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz}"

  mkdir -p "${output_dir}"
  curl -fsSL "${url}" -o "${tmp_file}"
  if [[ ! -s "${tmp_file}" ]]; then
    echo "Failed to download SILVA rRNA reference: ${url}" >&2
    rm -f "${tmp_file}"
    return 1
  fi
  gzip -t "${tmp_file}"
  mv "${tmp_file}" "${output_file}"
}

ensure_silva_rrna_ref_db() {
  local dir_pg=$1
  local sys_file="/usr/local/db/silva/rRNA_ref.fa.gz"
  local runtime_root="$(workspace_db_root "${dir_pg}")"
  local runtime_file="${runtime_root}/silva/rRNA_ref.fa.gz"
  local lock_file="${runtime_root}/locks/silva_rrna_ref.lock"

  if [[ -s "${sys_file}" ]]; then
    echo "${sys_file}"
    return 0
  fi
  if [[ -s "${runtime_file}" ]]; then
    echo "${runtime_file}"
    return 0
  fi

  mkdir -p "${runtime_root}"
  gg_array_download_once "${lock_file}" "${runtime_file}" "SILVA rRNA reference FASTA" \
    _download_silva_rrna_ref_to_file "${runtime_file}" || return 1

  if [[ ! -s "${runtime_file}" ]]; then
    echo "SILVA rRNA reference download failed." >&2
    return 1
  fi
  echo "${runtime_file}"
}

_download_mmseqs_uniref90_db() {
  local db_dir=$1
  local nthreads=$2
  local uniref_db="UniRef90"
  mkdir -p "${db_dir}"
  mmseqs databases "${uniref_db}" "${uniref_db}_DB" "${db_dir}" --threads "${nthreads}" || return 1
  if [[ ! -s "${db_dir}/UniRef90_DB" || ! -s "${db_dir}/UniRef90_DB.dbtype" ]]; then
    return 1
  fi
  touch "${db_dir}/UniRef90_DB.ready"
}

ensure_mmseqs_uniref90_db() {
  local db_dir=$1
  local nthreads=${2:-1}
  local db_file="${db_dir}/UniRef90_DB"
  local dbtype_file="${db_file}.dbtype"
  local done_file="${db_dir}/UniRef90_DB.ready"
  local lock_file="${db_dir}/locks/uniref90.lock"

  if [[ -s "${db_file}" && -s "${dbtype_file}" && -s "${done_file}" ]]; then
    return 0
  fi

  mkdir -p "${db_dir}"
  if [[ -s "${db_file}" && -s "${dbtype_file}" && ! -s "${done_file}" ]]; then
    echo "MMseqs2 UniRef90 DB file exists without ready marker. Reusing existing DB and creating ready marker." >&2
    touch "${done_file}"
    return 0
  fi
  gg_array_download_once "${lock_file}" "${done_file}" "MMseqs2 UniRef90 taxonomy DB" \
    _download_mmseqs_uniref90_db "${db_dir}" "${nthreads}" || return 1

  if [[ ! -s "${db_file}" || ! -s "${dbtype_file}" || ! -s "${done_file}" ]]; then
    echo "MMseqs2 UniRef90 DB download failed." >&2
    return 1
  fi
}

get_hyphy_genetic_code() {
  local ncbi_genetic_code=$1
  if [[ ${ncbi_genetic_code} -eq 1 ]]; then
    echo "Universal"
  elif [[ ${ncbi_genetic_code} -eq 2 ]]; then
    echo "Vertebrate-mtDNA"
  elif [[ ${ncbi_genetic_code} -eq 3 ]]; then
    echo "Yeast-mtDNA"
  elif [[ ${ncbi_genetic_code} -eq 4 ]]; then
    echo "Mold-Protozoan-mtDNA"
  elif [[ ${ncbi_genetic_code} -eq 5 ]]; then
    echo "Invertebrate-mtDNA"
  elif [[ ${ncbi_genetic_code} -eq 6 ]]; then
    echo "Ciliate-Nuclear"
  elif [[ ${ncbi_genetic_code} -eq 9 ]]; then
    echo "Echinoderm-mtDNA"
  elif [[ ${ncbi_genetic_code} -eq 10 ]]; then
    echo "Euplotid-Nuclear"
  elif [[ ${ncbi_genetic_code} -eq 12 ]]; then
    echo "Alt-Yeast-Nuclear"
  elif [[ ${ncbi_genetic_code} -eq 13 ]]; then
    echo "Ascidian-mtDNA"
  elif [[ ${ncbi_genetic_code} -eq 14 ]]; then
    echo "Flatworm-mtDNA"
  elif [[ ${ncbi_genetic_code} -eq 15 ]]; then
    echo "Blepharisma-Nuclear"
  elif [[ ${ncbi_genetic_code} -eq 16 ]]; then
    echo "Chlorophycean-mtDNA"
  elif [[ ${ncbi_genetic_code} -eq 21 ]]; then
    echo "Trematode-mtDNA"
  elif [[ ${ncbi_genetic_code} -eq 22 ]]; then
    echo "Scenedesmus-obliquus-mtDNA"
  elif [[ ${ncbi_genetic_code} -eq 23 ]]; then
    echo "Thraustochytrium-mtDNA"
  elif [[ ${ncbi_genetic_code} -eq 24 ]]; then
    echo "Pterobranchia-mtDNA"
  elif [[ ${ncbi_genetic_code} -eq 25 ]]; then
    echo "SR1-and-Gracilibacteria"
  elif [[ ${ncbi_genetic_code} -eq 26 ]]; then
    echo "Pachysolen-Nuclear"
  elif [[ ${ncbi_genetic_code} -eq 29 ]]; then
    echo "Mesodinium-Nuclear"
  elif [[ ${ncbi_genetic_code} -eq 30 ]]; then
    echo "Peritrich-Nuclear"
  elif [[ ${ncbi_genetic_code} -eq 33 ]]; then
    echo "Cephalodiscidae-mtDNA"
  else
    echo "Unknown"
    >&2 echo "This NCBI genetic code cannot be used in HyPhy: ${ncbi_genetic_code}"
  fi
}

get_fasta_extensions_for_grep() {
  printf -- "-e %s " ".fa$" ".fa.gz$" ".fas$" ".fas.gz$" ".fasta$" ".fasta.gz$" ".fna$" ".fna.gz$"
}

is_species_set_identical() {
  local return_flag=0
  local dir1=$1
  local dir2=$2
  local files1=()
  local files2=()
  shopt -s nullglob
  files1=( "${dir1}"/* )
  files2=( "${dir2}"/* )
  shopt -u nullglob
  local num_sp1=${#files1[@]}
  local num_sp2=${#files2[@]}
  if [[ ${num_sp1} -ne ${num_sp2} ]]; then
    echo "Number of species in ${dir1} and ${dir2} are different."
    return_flag=1
  fi
  local sp1=()
  local sp2=()
  local f
  for f in "${files1[@]}"; do
    sp1+=( "$(gg_species_name_from_path_or_dot "${f}")" )
  done
  for f in "${files2[@]}"; do
    sp2+=( "$(gg_species_name_from_path_or_dot "${f}")" )
  done
  mapfile -t sp1 < <(printf '%s\n' "${sp1[@]}" | sort -u)
  mapfile -t sp2 < <(printf '%s\n' "${sp2[@]}" | sort -u)
  if [[ "${sp1[*]}" != "${sp2[*]}" ]]; then
    echo "Species names are different between ${dir1} and ${dir2}"
    echo "dir1 (${dir1}): ${sp1[*]}"
    echo "dir2 (${dir2}): ${sp2[*]}"
    return_flag=1
  fi
  return "${return_flag}"
}
