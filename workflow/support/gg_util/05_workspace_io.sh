# shellcheck shell=bash
# Workspace layout, portable file movement, and input validation helpers.
# This file is sourced by workflow/support/gg_util.sh.

workspace_input_root() {
  local gg_workspace_dir=$1
  echo "${gg_workspace_dir}/input"
}

workspace_output_root() {
  local gg_workspace_dir=$1
  echo "${gg_workspace_dir}/output"
}

workspace_downloads_root() {
  local gg_workspace_dir=$1
  echo "${gg_workspace_dir}/downloads"
}

gg_resolve_workspace_layout() {
  local gg_workspace_dir=$1
  : "${gg_workspace_dir:?workspace directory is required}"
  echo "split"
}

workspace_taxonomy_root() {
  local gg_workspace_dir=$1
  local dir_db
  dir_db=$(workspace_downloads_root "${gg_workspace_dir}")
  echo "${dir_db}/ete_taxonomy"
}

workspace_pfam_root() {
  local gg_workspace_dir=$1
  local dir_db
  dir_db=$(workspace_downloads_root "${gg_workspace_dir}")
  echo "${dir_db}/pfam"
}

workspace_omark_root() {
  local gg_workspace_dir=$1
  local dir_db
  dir_db=$(workspace_downloads_root "${gg_workspace_dir}")
  echo "${dir_db}/omark"
}

workspace_omark_dbfile() {
  local gg_workspace_dir=$1
  local dir_omark
  dir_omark=$(workspace_omark_root "${gg_workspace_dir}")
  echo "${dir_omark}/LUCA.h5"
}

workspace_pfam_le_dir() {
  local gg_workspace_dir=$1
  local dir_pfam
  dir_pfam=$(workspace_pfam_root "${gg_workspace_dir}")
  echo "${dir_pfam}/Pfam_LE"
}

workspace_taxonomy_dbfile() {
  local gg_workspace_dir=$1
  local dir_taxonomy
  dir_taxonomy=$(workspace_taxonomy_root "${gg_workspace_dir}")
  echo "${dir_taxonomy}/taxa.sqlite"
}

workspace_taxonomy_taxdumpfile() {
  local gg_workspace_dir=$1
  local dir_taxonomy
  dir_taxonomy=$(workspace_taxonomy_root "${gg_workspace_dir}")
  echo "${dir_taxonomy}/taxdump.tar.gz"
}

gg_set_taxonomy_cache_env() {
  local gg_workspace_dir=$1
  local dir_taxonomy
  dir_taxonomy=$(workspace_taxonomy_root "${gg_workspace_dir}")
  ensure_dir "${dir_taxonomy}"
  ensure_dir "${dir_taxonomy}/ete"
  ensure_dir "${dir_taxonomy}/ete4"
  export ETE_DATA_HOME="${dir_taxonomy}"
  export ETE_CONFIG_HOME="${dir_taxonomy}"
  export XDG_DATA_HOME="${dir_taxonomy}"
  export XDG_CONFIG_HOME="${dir_taxonomy}"
  export GG_TAXONOMY_DBFILE="${dir_taxonomy}/taxa.sqlite"
  export GG_TAXONOMY_TAXDUMPFILE="${dir_taxonomy}/taxdump.tar.gz"
}

_ensure_ete_taxonomy_db_locked() {
  local db_file=$1
  local taxdump_file=$2
  local py_exec=""
  py_exec=$(gg_find_python_exec || true)

  if [[ -z "${py_exec}" ]]; then
    echo "python/python3 command was not found. Cannot prepare ETE taxonomy DB." >&2
    return 1
  fi

  if ! ETE_DATA_HOME="$(dirname "${db_file}")" \
    ETE_CONFIG_HOME="$(dirname "${db_file}")" \
    XDG_DATA_HOME="$(dirname "${db_file}")" \
    XDG_CONFIG_HOME="$(dirname "${db_file}")" \
    GG_TAXONOMY_DBFILE="${db_file}" \
    GG_TAXONOMY_TAXDUMPFILE="${taxdump_file}" \
    "${py_exec}" - <<'PY'
import importlib
import os
import urllib.request

db_file = os.environ.get("GG_TAXONOMY_DBFILE", "").strip()
if not db_file:
    raise SystemExit("GG_TAXONOMY_DBFILE is empty.")
taxdump_file = os.environ.get("GG_TAXONOMY_TAXDUMPFILE", "").strip()
if not taxdump_file:
    raise SystemExit("GG_TAXONOMY_TAXDUMPFILE is empty.")

os.makedirs(os.path.dirname(db_file), exist_ok=True)
os.makedirs(os.path.dirname(taxdump_file), exist_ok=True)

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
    os.makedirs(os.path.join(cache_dir, "ete4"), exist_ok=True)

def ensure_taxdump_file():
    if os.path.exists(taxdump_file):
        return taxdump_file
    tmp_path = taxdump_file + ".tmp"
    if os.path.lexists(tmp_path):
        os.remove(tmp_path)
    urllib.request.urlretrieve(
        "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
        tmp_path,
    )
    os.replace(tmp_path, taxdump_file)
    return taxdump_file

def ensure_with_ete4():
    ncbiquery = importlib.import_module("ete4.ncbi_taxonomy.ncbiquery")
    if os.path.exists(db_file) and ncbiquery.is_taxadb_up_to_date(db_file):
        return "ete4:up_to_date"
    NCBITaxa = importlib.import_module("ete4").NCBITaxa
    NCBITaxa(dbfile=db_file, taxdump_file=ensure_taxdump_file(), update=True)
    return "ete4:updated"

try:
    ensure_with_ete4()
except Exception as exc:  # pragma: no cover - runtime dependency handling
    raise SystemExit(f"Failed to prepare ETE taxonomy DB with ete4: {exc}")

if not os.path.exists(db_file) or os.path.getsize(db_file) == 0:
    raise SystemExit(f"ETE taxonomy DB was not generated: {db_file}")
PY
  then
    echo "Failed to prepare ETE taxonomy DB: ${db_file}" >&2
    return 1
  fi
  if [[ ! -s "${db_file}" ]]; then
    echo "Failed to prepare ETE taxonomy DB: ${db_file}" >&2
    return 1
  fi
}

ensure_ete_taxonomy_db() {
  local gg_workspace_dir=$1
  local dir_db
  local dir_taxonomy
  local db_file
  local taxdump_file
  local lock_file

  dir_db=$(workspace_downloads_root "${gg_workspace_dir}")
  dir_taxonomy=$(workspace_taxonomy_root "${gg_workspace_dir}")
  db_file=$(workspace_taxonomy_dbfile "${gg_workspace_dir}")
  taxdump_file=$(workspace_taxonomy_taxdumpfile "${gg_workspace_dir}")
  lock_file="${dir_db}/locks/ete_taxonomy.lock"

  ensure_dir "${dir_db}"
  ensure_dir "${dir_taxonomy}"
  ensure_dir "${dir_taxonomy}/ete"
  ensure_dir "${dir_taxonomy}/ete4"
  ensure_dir "$(dirname "${lock_file}")"
  export ETE_DATA_HOME="${dir_taxonomy}"
  export ETE_CONFIG_HOME="${dir_taxonomy}"
  export XDG_DATA_HOME="${dir_taxonomy}"
  export XDG_CONFIG_HOME="${dir_taxonomy}"
  export GG_TAXONOMY_DBFILE="${db_file}"
  export GG_TAXONOMY_TAXDUMPFILE="${taxdump_file}"

  gg_array_download_once "${lock_file}" "${db_file}" "ETE taxonomy DB" \
    _ensure_ete_taxonomy_db_locked "${db_file}" "${taxdump_file}" || return 1
  if [[ ! -s "${db_file}" ]]; then
    echo "Failed to prepare ETE taxonomy DB: ${db_file}" >&2
    return 1
  fi
  return 0
}

_download_omark_database_locked() {
  local target_file=$1
  local ready_file=$2
  local url=$3
  local tmp_file="${target_file}.part"

  ensure_parent_dir "${target_file}"
  rm -f -- "${ready_file}"
  if [[ -f "${tmp_file}" && ! -s "${tmp_file}" ]]; then
    rm -f -- "${tmp_file}"
  fi

  echo "Starting OMArk database download: ${url}" >&2
  if ! curl -fL --retry 5 --retry-all-errors --retry-delay 2 --continue-at - "${url}" -o "${tmp_file}"; then
    echo "OMArk database download failed: ${url}" >&2
    return 1
  fi
  mv -- "${tmp_file}" "${target_file}" || return 1
  if [[ ! -s "${target_file}" ]]; then
    echo "Downloaded OMArk database is empty: ${target_file}" >&2
    return 1
  fi
  gg_write_ready_marker "${ready_file}"
  echo "OMArk database download has been finished: ${target_file}" >&2
}

ensure_omark_database() {
  local gg_workspace_dir=$1
  local requested=${2:-auto}
  local normalized_requested=""
  local requested_lc=""
  local target_file=""
  local ready_file=""
  local lock_file=""
  local resolved_url=""

  normalized_requested=$(printf '%s' "${requested}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
  requested_lc=$(printf '%s' "${normalized_requested}" | tr '[:upper:]' '[:lower:]')
  if [[ -n "${normalized_requested}" && "${requested_lc}" != "auto" ]]; then
    target_file="${normalized_requested}"
  else
    target_file=$(workspace_omark_dbfile "${gg_workspace_dir}")
  fi

  if [[ -n "${normalized_requested}" && "${requested_lc}" != "auto" && -s "${target_file}" ]]; then
    printf '%s\n' "${target_file}"
    return 0
  fi

  ready_file="${target_file}.download.ready"
  lock_file="$(dirname "${target_file}")/.LUCA.h5.download.lock"
  ensure_parent_dir "${target_file}"

  if [[ -s "${target_file}" && -s "${ready_file}" ]]; then
    printf '%s\n' "${target_file}"
    return 0
  fi

  if ! resolved_url=$(gg_resolve_omark_database_url); then
    echo "Failed to resolve OMArk database download URL." >&2
    return 1
  fi

  gg_array_download_once "${lock_file}" "${ready_file}" "OMArk OMAmer database (LUCA.h5)" \
    _download_omark_database_locked "${target_file}" "${ready_file}" "${resolved_url}" || return 1

  if [[ ! -s "${target_file}" ]]; then
    echo "Failed to prepare OMArk database: ${target_file}" >&2
    return 1
  fi

  printf '%s\n' "${target_file}"
}

gg_initialize_data_layout() {
  local gg_workspace_dir=$1
  local dir_input
  local dir_output
  local dir_db
  dir_input=$(workspace_input_root "${gg_workspace_dir}")
  dir_output=$(workspace_output_root "${gg_workspace_dir}")
  dir_db=$(workspace_downloads_root "${gg_workspace_dir}")
  ensure_dir "${dir_input}"
  ensure_dir "${dir_output}"
  ensure_dir "${dir_db}"
  gg_set_taxonomy_cache_env "${gg_workspace_dir}"
}

gg_copy_file_portable() {
	if [[ $# -ne 2 ]]; then
		echo "gg_copy_file_portable: exactly 2 arguments are required."
		return 1
	fi
	local source_path=$1
	local dest_path=$2
	local python_exec
	python_exec=$(gg_find_python_exec) || return 1
	"${python_exec}" - "${source_path}" "${dest_path}" <<'PY'
import os
import shutil
import sys

source_path, dest_path = sys.argv[1:3]
if dest_path.endswith(os.sep) or os.path.isdir(dest_path):
    os.makedirs(dest_path, exist_ok=True)
    dest_path = os.path.join(dest_path, os.path.basename(source_path))
else:
    parent = os.path.dirname(dest_path)
    if parent:
        os.makedirs(parent, exist_ok=True)
shutil.copy2(source_path, dest_path)
PY
}

gg_copy_dir_portable() {
	if [[ $# -ne 2 ]]; then
		echo "gg_copy_dir_portable: exactly 2 arguments are required."
		return 1
	fi
	local source_dir=$1
	local dest_dir=$2
	local python_exec
	python_exec=$(gg_find_python_exec) || return 1
	"${python_exec}" - "${source_dir}" "${dest_dir}" <<'PY'
import os
import shutil
import sys

source_dir, dest_dir = sys.argv[1:3]
parent = os.path.dirname(dest_dir)
if parent:
    os.makedirs(parent, exist_ok=True)
shutil.copytree(source_dir, dest_dir, dirs_exist_ok=True)
PY
}

_gg_atomic_stream_to_path() (
	local destination=$1
	local parent
	local basename
	local token
	local staged
	parent=$(dirname -- "${destination}")
	basename=$(basename -- "${destination}")
	ensure_dir "${parent}" || return 1
	token="${BASHPID:-$$}.${RANDOM}"
	staged="${parent}/.${basename}.gg-stream.${token}"
	if [[ -e "${staged}" || -L "${staged}" ]]; then
		echo "Atomic stream staging path already exists: ${staged}" >&2
		return 1
	fi
	trap 'rm -f -- "${staged}"' EXIT HUP INT TERM
	cat > "${staged}" || return 1
	mv_out_bundle "${staged}" "${destination}" || return 1
	trap - EXIT HUP INT TERM
)

_gg_atomic_copy_one() (
	local source=$1
	local destination=$2
	local parent
	local basename
	local token
	local staged
	parent=$(dirname -- "${destination}")
	basename=$(basename -- "${destination}")
	ensure_dir "${parent}" || return 1
	token="${BASHPID:-$$}.${RANDOM}"
	staged="${parent}/.${basename}.gg-copy.${token}"
	if [[ -e "${staged}" || -L "${staged}" ]]; then
		echo "Atomic copy staging path already exists: ${staged}" >&2
		return 1
	fi
	trap 'rm -rf -- "${staged}"' EXIT HUP INT TERM
	cp -- "${source}" "${staged}" || return 1
	mv_out_bundle "${staged}" "${destination}" || return 1
	trap - EXIT HUP INT TERM
)

cp_out() {
	if [[ $# -eq 1 ]]; then
		if [[ -p /dev/stdin ]]; then
			_gg_atomic_stream_to_path "$1"
			return $?
		fi
		echo "cp_out: at least 2 arguments are required unless stdin is piped."
		return 1
	fi
	if [[ $# -eq 2 && "$1" == "-" ]]; then
		_gg_atomic_stream_to_path "$2"
		return $?
	fi
	if [[ $# -lt 2 ]]; then
		echo "cp_out: at least 2 arguments are required."
		return 1
	fi
	local destination_argument="${!#}"
	local source_count=$(( $# - 1 ))
	local source_index
	local source
	local destination
	if [[ ${source_count} -gt 1 || "${destination_argument}" == */ || -d "${destination_argument}" ]]; then
		ensure_dir "${destination_argument%/}" || return 1
		for ((source_index = 1; source_index <= source_count; source_index++)); do
			source=${!source_index}
			destination="${destination_argument%/}/$(basename -- "${source}")"
			_gg_atomic_copy_one "${source}" "${destination}" || return 1
		done
		return 0
	fi
	_gg_atomic_copy_one "$1" "${destination_argument}"
}

mv_out() {
	if [[ $# -eq 1 ]]; then
		if [[ -p /dev/stdin ]]; then
			_gg_atomic_stream_to_path "$1"
			return $?
		fi
		echo "mv_out: at least 2 arguments are required unless stdin is piped."
		return 1
	fi
	if [[ $# -eq 2 && "$1" == "-" ]]; then
		_gg_atomic_stream_to_path "$2"
		return $?
	fi
	if [[ $# -lt 2 ]]; then
		echo "mv_out: at least 2 arguments are required."
		return 1
	fi
	local dest="${!#}"
	if [[ $# -gt 2 || "${dest}" == */ ]]; then
		ensure_dir "${dest%/}"
	else
		ensure_parent_dir "${dest}"
	fi
	mv -- "$@"
}

_gg_publish_lock_acquire() {
	local lock_dir=$1
	local description=$2
	local poll_seconds
	local timeout_seconds
	local stale_seconds
	local wait_started
	local wait_logged=0
	poll_seconds=$(gg_lock_poll_seconds)
	timeout_seconds=$(gg_lock_acquire_timeout_seconds)
	stale_seconds=$(gg_lock_stale_seconds)
	wait_started=$(date +%s)
	ensure_parent_dir "${lock_dir}" || return 1
	while true; do
		if mkdir -- "${lock_dir}" 2>/dev/null; then
			return 0
		fi
		local now_epoch
		local lock_mtime
		now_epoch=$(date +%s)
		lock_mtime=$(gg_stat_mtime_epoch "${lock_dir}")
		if [[ "${lock_mtime}" =~ ^[0-9]+$ ]] && (( now_epoch - lock_mtime >= stale_seconds )); then
			if rmdir -- "${lock_dir}" 2>/dev/null; then
				echo "Recovered stale publication lock: ${description}" >&2
				continue
			fi
		fi
		if [[ ${wait_logged} -eq 0 ]]; then
			echo "Waiting for publication lock: ${description}" >&2
			wait_logged=1
		fi
		if (( now_epoch - wait_started >= timeout_seconds )); then
			echo "Timed out waiting for publication lock: ${description}" >&2
			return 1
		fi
		sleep "${poll_seconds}"
	done
}

_gg_publish_lock_release() {
	local lock_dir=$1
	rmdir -- "${lock_dir}"
}

mv_out_bundle() (
	# Publish source/destination pairs as one recoverable transaction.  Every
	# source is first staged on its destination filesystem.  Existing outputs
	# are retained until all staging succeeds, and any ordinary error or caught
	# signal restores the complete previous bundle.
	if [[ $# -eq 0 || $(( $# % 2 )) -ne 0 ]]; then
		echo "mv_out_bundle: source/destination argument pairs are required."
		return 1
	fi
	local -a sources=()
	local -a destinations=()
	local -a canonical_sources=()
	local -a canonical_destinations=()
	local -a staged=()
	local -a backups=()
	local -a had_destination=()
	local -a stage_attempted=()
	local -a publish_attempted=()
	local -a bundle_lock_paths=()
	local -a bundle_acquired_paths=()
	local -a bundle_lock_heartbeat_pids=()
	local argument_index pair_index previous_index source destination parent basename token canonical_source canonical_destination lock_path swap_path
	local pair_count=$(( $# / 2 ))
	local LC_ALL=C
	for ((argument_index = 1; argument_index <= $#; argument_index += 2)); do
		source=${!argument_index}
		pair_index=$(( argument_index + 1 ))
		destination=${!pair_index}
		if [[ ! -e "${source}" && ! -L "${source}" ]]; then
			echo "mv_out_bundle: source not found: ${source}"
			return 1
		fi
		if [[ -z "${destination}" || "${destination}" == "/" || "${destination}" == "." || "${destination}" == ".." ]]; then
			echo "mv_out_bundle: unsafe destination: ${destination}"
			return 1
		fi
		parent=$(dirname -- "${destination}")
		ensure_dir "${parent}" || return 1
		canonical_source=$(gg_resolve_physical_path "${source}") || return 1
		canonical_destination=$(gg_resolve_physical_path "${destination}") || return 1
		if [[ "${canonical_source}" == "${canonical_destination}" || ( -e "${destination}" && "${source}" -ef "${destination}" ) ]]; then
			echo "mv_out_bundle: source and destination must differ: ${source}"
			return 1
		fi
		for ((previous_index = 0; previous_index < ${#destinations[@]}; previous_index++)); do
			if [[ "${canonical_destinations[previous_index]}" == "${canonical_destination}" || ( -e "${destinations[previous_index]}" && -e "${destination}" && "${destinations[previous_index]}" -ef "${destination}" ) ]]; then
				echo "mv_out_bundle: duplicate destination: ${destination}"
				return 1
			fi
			if [[ "${canonical_sources[previous_index]}" == "${canonical_source}" || "${sources[previous_index]}" -ef "${source}" ]]; then
				echo "mv_out_bundle: duplicate source: ${source}"
				return 1
			fi
		done
		sources+=("${source}")
		destinations+=("${destination}")
		canonical_sources+=("${canonical_source}")
		canonical_destinations+=("${canonical_destination}")
	done
	for ((pair_index = 0; pair_index < pair_count; pair_index++)); do
		for ((previous_index = 0; previous_index < pair_count; previous_index++)); do
			if [[ "${canonical_sources[pair_index]}" == "${canonical_destinations[previous_index]}" || ( -e "${destinations[previous_index]}" && "${sources[pair_index]}" -ef "${destinations[previous_index]}" ) ]]; then
				echo "mv_out_bundle: sources and destinations must not overlap: ${sources[pair_index]}"
				return 1
			fi
		done
	done
	token="${BASHPID:-$$}.${RANDOM}"
	for ((pair_index = 0; pair_index < pair_count; pair_index++)); do
		destination=${destinations[pair_index]}
		parent=$(dirname -- "${destination}")
		basename=$(basename -- "${destination}")
		ensure_dir "${parent}" || return 1
		staged+=("${parent}/.${basename}.gg-stage.${token}.${pair_index}")
		backups+=("${parent}/.${basename}.gg-backup.${token}.${pair_index}")
		had_destination+=("no")
		stage_attempted+=("no")
		publish_attempted+=("no")
		bundle_lock_paths+=("${parent}/.${basename}.gg-bundle.lock")
		if [[ -e "${staged[pair_index]}" || -L "${staged[pair_index]}" || -e "${backups[pair_index]}" || -L "${backups[pair_index]}" ]]; then
			echo "mv_out_bundle: transaction path collision for destination: ${destination}"
			return 1
		fi
	done
	# Acquire every destination lock in a deterministic order so overlapping
	# bundles cannot race even when their first destinations differ.
	for ((pair_index = 0; pair_index < pair_count; pair_index++)); do
		for ((previous_index = pair_index + 1; previous_index < pair_count; previous_index++)); do
			if [[ "${bundle_lock_paths[pair_index]}" > "${bundle_lock_paths[previous_index]}" ]]; then
				swap_path=${bundle_lock_paths[pair_index]}
				bundle_lock_paths[pair_index]="${bundle_lock_paths[previous_index]}"
				bundle_lock_paths[previous_index]="${swap_path}"
			fi
		done
	done
	local committed=0
	_mv_out_bundle_release_lock() {
		local release_index
		for ((release_index = ${#bundle_acquired_paths[@]} - 1; release_index >= 0; release_index--)); do
			gg_shared_lock_stop_heartbeat "${bundle_lock_heartbeat_pids[release_index]:-}"
			_gg_publish_lock_release "${bundle_acquired_paths[release_index]}" || true
		done
		bundle_acquired_paths=()
		bundle_lock_heartbeat_pids=()
	}
	_mv_out_bundle_rollback() {
		local rollback_index
		if [[ ${committed} -eq 1 ]]; then
			return
		fi
		# Complete one deterministic rollback even if a scheduler repeats the
		# termination signal while the EXIT trap is restoring files.
		trap '' HUP INT TERM
		for ((rollback_index = pair_count - 1; rollback_index >= 0; rollback_index--)); do
			if [[ "${publish_attempted[rollback_index]}" == "yes" && ! -e "${staged[rollback_index]}" && ! -L "${staged[rollback_index]}" && ( -e "${destinations[rollback_index]}" || -L "${destinations[rollback_index]}" ) ]]; then
				if [[ ! -e "${sources[rollback_index]}" && ! -L "${sources[rollback_index]}" ]]; then
					mv -- "${destinations[rollback_index]}" "${sources[rollback_index]}" || true
				else
					mv -- "${destinations[rollback_index]}" "${staged[rollback_index]}" || true
				fi
			fi
		done
		for ((rollback_index = pair_count - 1; rollback_index >= 0; rollback_index--)); do
			if [[ "${had_destination[rollback_index]}" == "yes" && ( -e "${backups[rollback_index]}" || -L "${backups[rollback_index]}" ) && ! -e "${destinations[rollback_index]}" && ! -L "${destinations[rollback_index]}" ]]; then
				mv -- "${backups[rollback_index]}" "${destinations[rollback_index]}" || true
			fi
		done
		for ((rollback_index = pair_count - 1; rollback_index >= 0; rollback_index--)); do
			if [[ "${stage_attempted[rollback_index]}" == "yes" && ( -e "${staged[rollback_index]}" || -L "${staged[rollback_index]}" ) ]]; then
				if [[ ! -e "${sources[rollback_index]}" && ! -L "${sources[rollback_index]}" ]]; then
					mv -- "${staged[rollback_index]}" "${sources[rollback_index]}" || true
				else
					rm -rf -- "${staged[rollback_index]}" || true
				fi
			fi
		done
	}
	_mv_out_bundle_cleanup() {
		_mv_out_bundle_rollback
		_mv_out_bundle_release_lock
	}
	trap '_mv_out_bundle_cleanup' EXIT
	trap 'exit 130' HUP INT TERM

	for ((pair_index = 0; pair_index < pair_count; pair_index++)); do
		stage_attempted[pair_index]="yes"
		mv -- "${sources[pair_index]}" "${staged[pair_index]}" || return 1
	done
	for lock_path in "${bundle_lock_paths[@]}"; do
		if ! _gg_publish_lock_acquire "${lock_path}" "output bundle publication (${lock_path})"; then
			echo "mv_out_bundle: failed to acquire publication lock: ${lock_path}" >&2
			return 1
		fi
		bundle_acquired_paths+=("${lock_path}")
		gg_shared_lock_start_heartbeat "${lock_path}"
		bundle_lock_heartbeat_pids+=("${GG_SHARED_LOCK_HEARTBEAT_PID:-}")
	done
	for ((pair_index = 0; pair_index < pair_count; pair_index++)); do
		if [[ -e "${destinations[pair_index]}" || -L "${destinations[pair_index]}" ]]; then
			had_destination[pair_index]="yes"
			mv -- "${destinations[pair_index]}" "${backups[pair_index]}" || return 1
		fi
	done
	for ((pair_index = 0; pair_index < pair_count; pair_index++)); do
		publish_attempted[pair_index]="yes"
		mv -- "${staged[pair_index]}" "${destinations[pair_index]}" || return 1
	done
	committed=1
	for ((pair_index = 0; pair_index < pair_count; pair_index++)); do
		if [[ "${had_destination[pair_index]}" == "yes" ]]; then
			rm -rf -- "${backups[pair_index]}" || return 1
		fi
	done
	_mv_out_bundle_release_lock
	trap - EXIT HUP INT TERM
)

mv_out_replace_dir() {
	if [[ $# -ne 2 ]]; then
		echo "mv_out_replace_dir: exactly 2 arguments are required."
		return 1
	fi
	local staged_dir=$1
	local dest_dir=$2
	if [[ ! -d "${staged_dir}" ]]; then
		echo "mv_out_replace_dir: source directory not found: ${staged_dir}"
		return 1
	fi
	mv_out_bundle "${staged_dir}" "${dest_dir}"
}

resolve_rnaspades_transcript_fasta() {
	if [[ $# -ne 1 ]]; then
		echo "resolve_rnaspades_transcript_fasta: exactly 1 argument is required."
		return 1
	fi
	local output_dir=$1
	local candidate=""
	local candidate_files=(
		"transcripts.fasta"
		"soft_filtered_transcripts.fasta"
		"hard_filtered_transcripts.fasta"
	)
	for candidate in "${candidate_files[@]}"; do
		if [[ -s "${output_dir%/}/${candidate}" ]]; then
			printf '%s\n' "${output_dir%/}/${candidate}"
			return 0
		fi
	done
	return 1
}

capture_busco_repro_artifacts() {
	if [[ $# -lt 3 || $# -gt 6 ]]; then
		echo "capture_busco_repro_artifacts: expected 3 to 6 arguments: repro_dir input_fasta busco_tmp_dir [lineage] [stage_key] [stderr_log]"
		return 1
	fi
	local repro_dir=$1
	local input_fasta=$2
	local busco_tmp_dir=$3
	local lineage=${4:-}
	local stage_key=${5:-}
	local stderr_log=${6:-}

	recreate_dir "${repro_dir}" || return 1
	if [[ -s "${input_fasta}" ]]; then
		gg_copy_file_portable "${input_fasta}" "${repro_dir}/" || return 1
	fi
	if [[ -d "${busco_tmp_dir}" ]]; then
		gg_copy_dir_portable "${busco_tmp_dir}" "${repro_dir}/busco_tmp" || return 1
	fi
	if [[ -s "${stderr_log}" ]]; then
		gg_copy_file_portable "${stderr_log}" "${repro_dir}/busco.stderr.log" || return 1
	fi
	{
		printf 'stage_key\t%s\n' "${stage_key}"
		printf 'lineage\t%s\n' "${lineage}"
		printf 'input_fasta\t%s\n' "${input_fasta}"
		printf 'busco_tmp_dir\t%s\n' "${busco_tmp_dir}"
		printf 'stderr_log\t%s\n' "${stderr_log}"
	} > "${repro_dir}/capture_info.tsv"
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
	local d
	while IFS= read -r -d '' d; do
		sub_directories+=( "${d}" )
	done < <(find "${dir_main}" -mindepth 1 -maxdepth 1 -type d -print0)
	for d in "${sub_directories[@]}"; do
		if [[ -z "$(find "${d}" -mindepth 1 -maxdepth 1 -print -quit 2>/dev/null)" ]]; then
			echo "${echo_header}deleting ${d}"
				rm -rf -- "${d}"
		fi
	done
	echo ""
}

set_singularity_command() {
  local echo_header="set_singularity_command: "
  local requested_runtime=""
  local runtime_bin
  local command_display=""
  requested_runtime="$(gg_requested_container_runtime)" || return 1
  if [[ "${requested_runtime}" == "docker" || "${requested_runtime}" == "auto" ]]; then
    gg_auto_enable_docker_runtime_if_available || true
  fi
  if ! runtime_bin=$(gg_detect_container_runtime_binary); then
    if [[ "${GG_CONTAINER_RUNTIME:-auto}" == "docker" ]]; then
      echo "${echo_header}Docker-backed runtime is unavailable."
    else
      echo "${echo_header}Neither singularity nor apptainer was found on PATH."
    fi
    return 1
  fi
  echo ${echo_header}"hostname = $(hostname)"
  echo ${echo_header}"container runtime = ${runtime_bin}"
  if declare -F gg_site_container_shell_command >/dev/null 2>&1; then
    gg_site_container_shell_command "${runtime_bin}" singularity_command || return 1
  else
    echo ${echo_header}"No site adapter was loaded. Using default exec command."
    singularity_command=( "${runtime_bin}" exec )
  fi
  command_display=$(gg_container_shell_command_display || true)
  echo ${echo_header}'${singularity_command[@]}' = \"${command_display}\"
  echo ""
}

gg_entrypoint_prepare_container_runtime() {
	local call_exit_if_running=${1:-0}
	gg_scheduler_runtime_prelude
	unset_singularity_envs
	gg_normalize_scheduler_env
	if [[ "${call_exit_if_running}" -eq 1 ]]; then
		exit_if_running_qstat
	fi
	if ! set_singularity_command; then
		return 1
	fi
	return 0
}

gg_entrypoint_activate_container_runtime() {
	set_singularityenv
	gg_print_scheduler_runtime_summary
	gg_entrypoint_print_version_summary
}

gg_entrypoint_enter_workspace() {
	mkdir -p "${gg_workspace_dir}"
	cd "${gg_workspace_dir}" || exit 1
}
