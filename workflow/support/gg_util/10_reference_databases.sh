# shellcheck shell=bash
# Sourced by workflow/support/gg_util.sh.

ensure_pfam_le_db() {
  local gg_workspace_dir=$1
  local sys_dir="/usr/local/db/Pfam_LE"
  local runtime_root
  runtime_root=$(workspace_downloads_root "${gg_workspace_dir}")
  local runtime_parent
  runtime_parent=$(workspace_pfam_root "${gg_workspace_dir}")
  local runtime_dir
  runtime_dir=$(workspace_pfam_le_dir "${gg_workspace_dir}")
  local ready_file="${runtime_dir}/.pfam_le.ready"
  local lock_file="${runtime_root}/locks/pfam_le.lock"

  if pfam_le_db_is_ready "${sys_dir}"; then
    echo "${sys_dir}"
    return 0
  fi

  if pfam_le_db_is_ready "${runtime_dir}"; then
    gg_write_ready_marker "${ready_file}"
    echo "${runtime_dir}"
    return 0
  fi

  mkdir -p "${runtime_root}" "${runtime_parent}"
  gg_array_download_once "${lock_file}" "${ready_file}" "Pfam_LE RPS-BLAST DB" \
    _download_pfam_le_to_dir "${runtime_dir}" || return 1

  if ! pfam_le_db_is_ready "${runtime_dir}"; then
    echo "Pfam_LE DB download failed." >&2
    return 1
  fi
  echo "${runtime_dir}"
}

_download_jaspar_meme_to_file() {
  local output_file=$1
  local output_dir
  output_dir=$(dirname "${output_file}")
  local filename
  filename=$(basename "${output_file}")
  local tmp_file="${output_file}.tmp"
  local downloaded=0
  local url
  local urls=()
  local candidate_url

  while IFS= read -r candidate_url; do
    urls+=( "${candidate_url}" )
  done < <(_jaspar_meme_url_candidates "${filename}")
  if [[ ${#urls[@]} -eq 0 ]]; then
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
    rm -f -- "${tmp_file}"
    return 1
  fi

  mv -- "${tmp_file}" "${output_file}"
}

_jaspar_year_from_filename() {
  local filename=$1
  echo "${filename}" | sed -n 's/^JASPAR\([0-9][0-9][0-9][0-9]\)_.*/\1/p'
}

_jaspar_meme_url_candidates() {
  local filename=$1
  local year
  year=$(_jaspar_year_from_filename "${filename}")
  if [[ -z "${year}" ]]; then
    return 1
  fi
  echo "https://jaspar.elixir.no/download/data/${year}/CORE/${filename}"
  echo "https://jaspar.genereg.net/download/data/${year}/CORE/${filename}"
  echo "https://jaspar${year}.elixir.no/download/data/${year}/CORE/${filename}"
}

_jaspar_is_plants_core_meme_filename() {
  local filename=$1
  [[ "${filename}" == JASPAR[0-9][0-9][0-9][0-9]_CORE_plants_non-redundant_pfms_meme.txt ]]
}

_jaspar_is_latest_selector() {
  local selector=${1:-}
  [[ -z "${selector}" || "${selector}" == "latest" || "${selector}" == "auto" ]]
}

_jaspar_remote_meme_file_exists() {
  local filename=$1
  local url
  while IFS= read -r url; do
    if curl -fsSI --max-time 20 "${url}" >/dev/null 2>&1; then
      return 0
    fi
    if curl -fsSL --max-time 20 --range 0-0 "${url}" -o /dev/null >/dev/null 2>&1; then
      return 0
    fi
  done < <(_jaspar_meme_url_candidates "${filename}")
  return 1
}

_jaspar_find_latest_meme_filename_remote() {
  local skip_remote_lookup=${GG_JASPAR_SKIP_REMOTE_LOOKUP:-0}
  local max_year=${GG_JASPAR_MAX_YEAR:-$(date +%Y)}
  local min_year=${GG_JASPAR_MIN_YEAR:-2000}
  local year
  local candidate_filename

  if [[ "${skip_remote_lookup}" == "1" ]]; then
    return 1
  fi
  if [[ ! "${max_year}" =~ ^[0-9]+$ ]]; then
    max_year=$(date +%Y)
  fi
  if [[ ! "${min_year}" =~ ^[0-9]+$ ]]; then
    min_year=2000
  fi
  if (( max_year < min_year )); then
    min_year=${max_year}
  fi

  for ((year=max_year; year>=min_year; year--)); do
    candidate_filename="JASPAR${year}_CORE_plants_non-redundant_pfms_meme.txt"
    if _jaspar_remote_meme_file_exists "${candidate_filename}"; then
      echo "${candidate_filename}"
      return 0
    fi
  done
  return 1
}

_jaspar_find_latest_meme_filename_local() {
  local sys_dir=$1
  local runtime_dir=$2
  local candidate_files=()
  local file

  shopt -s nullglob
  for file in "${sys_dir}"/JASPAR*_CORE_plants_non-redundant_pfms_meme.txt "${runtime_dir}"/JASPAR*_CORE_plants_non-redundant_pfms_meme.txt; do
    if [[ -s "${file}" ]]; then
      candidate_files+=( "$(basename "${file}")" )
    fi
  done
  shopt -u nullglob

  if [[ ${#candidate_files[@]} -eq 0 ]]; then
    return 1
  fi

  local latest_candidate=""
  while IFS= read -r candidate; do
    if [[ -n "${candidate}" ]]; then
      latest_candidate="${candidate}"
      break
    fi
  done < <(
    printf '%s\n' "${candidate_files[@]}" \
      | sed -n 's/^\(JASPAR[0-9][0-9][0-9][0-9]_CORE_plants_non-redundant_pfms_meme\.txt\)$/\1/p' \
      | sort -u \
      | sort -r
  )
  if [[ -z "${latest_candidate}" ]]; then
    return 1
  fi
  printf '%s\n' "${latest_candidate}"
}

_ensure_jaspar_file_named() {
  local gg_workspace_dir=$1
  local jaspar_filename=$2
  local sys_file="/usr/local/db/jaspar/${jaspar_filename}"
  local runtime_root
  runtime_root=$(workspace_downloads_root "${gg_workspace_dir}")
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

_resolve_latest_jaspar_path_from_marker() {
  local latest_marker=$1
  local sys_dir=$2
  local runtime_dir=$3
  local resolved_filename=""

  if [[ ! -s "${latest_marker}" ]]; then
    return 1
  fi
  read -r resolved_filename < "${latest_marker}"
  if ! _jaspar_is_plants_core_meme_filename "${resolved_filename}"; then
    return 1
  fi
  if [[ -s "${sys_dir}/${resolved_filename}" ]]; then
    printf '%s\n' "${sys_dir}/${resolved_filename}"
    return 0
  fi
  if [[ -s "${runtime_dir}/${resolved_filename}" ]]; then
    printf '%s\n' "${runtime_dir}/${resolved_filename}"
    return 0
  fi
  return 1
}

_prepare_latest_jaspar_file_locked() {
  local gg_workspace_dir=$1
  local latest_marker=$2
  local sys_dir=$3
  local runtime_dir=$4
  local resolved_filename=""
  local resolved_path=""

  if resolved_path=$(_resolve_latest_jaspar_path_from_marker "${latest_marker}" "${sys_dir}" "${runtime_dir}"); then
    printf '%s\n' "${resolved_path}"
    return 0
  fi
  if resolved_filename=$(_jaspar_find_latest_meme_filename_remote); then
    :
  else
    resolved_filename=""
  fi
  if [[ -z "${resolved_filename}" ]]; then
    if resolved_filename=$(_jaspar_find_latest_meme_filename_local "${sys_dir}" "${runtime_dir}"); then
      :
    else
      resolved_filename=""
    fi
  fi
  if [[ -z "${resolved_filename}" ]]; then
    echo "Could not resolve latest JASPAR plants CORE MEME motif file." >&2
    return 1
  fi
  if ! resolved_path=$(_ensure_jaspar_file_named "${gg_workspace_dir}" "${resolved_filename}"); then
    return 1
  fi
  printf '%s\n' "${resolved_filename}" > "${latest_marker}.tmp"
  mv -- "${latest_marker}.tmp" "${latest_marker}"
  printf '%s\n' "${resolved_path}"
}

ensure_latest_jaspar_file() {
  local gg_workspace_dir=$1
  local runtime_root
  runtime_root=$(workspace_downloads_root "${gg_workspace_dir}")
  local runtime_dir="${runtime_root}/jaspar"
  local sys_dir="/usr/local/db/jaspar"
  local lock_file="${runtime_root}/locks/jaspar_latest.lock"
  local latest_marker="${runtime_dir}/latest_core_plants_non-redundant_pfms_meme.filename"
  local resolved_path=""
  local heartbeat_pid=""
  local ensure_exit_code=0

  mkdir -p "${runtime_root}" "${runtime_dir}" "$(dirname "${lock_file}")"

  if resolved_path=$(_resolve_latest_jaspar_path_from_marker "${latest_marker}" "${sys_dir}" "${runtime_dir}"); then
    echo "${resolved_path}"
    return 0
  fi

  if ! gg_shared_lock_acquire "${lock_file}" "latest JASPAR motif file"; then
    return 1
  fi
  gg_shared_lock_start_heartbeat "${lock_file}"
  heartbeat_pid=${GG_SHARED_LOCK_HEARTBEAT_PID:-}

  if resolved_path=$(_prepare_latest_jaspar_file_locked "${gg_workspace_dir}" "${latest_marker}" "${sys_dir}" "${runtime_dir}"); then
    ensure_exit_code=0
  else
    ensure_exit_code=$?
  fi

  gg_shared_lock_stop_heartbeat "${heartbeat_pid}"
  gg_shared_lock_release "${lock_file}"

  if [[ ${ensure_exit_code} -ne 0 ]]; then
    return ${ensure_exit_code}
  fi
  if [[ -z "${resolved_path}" || ! -s "${resolved_path}" ]]; then
    echo "Latest JASPAR motif file was not prepared correctly: ${resolved_path}" >&2
    return 1
  fi
  echo "${resolved_path}"
}

ensure_jaspar_file() {
  local gg_workspace_dir=$1
  local jaspar_filename=${2:-latest}

  if _jaspar_is_latest_selector "${jaspar_filename}"; then
    ensure_latest_jaspar_file "${gg_workspace_dir}"
    return $?
  fi
  _ensure_jaspar_file_named "${gg_workspace_dir}" "${jaspar_filename}"
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
    rm -f -- "${tmp_file}"
    return 1
  fi
  gzip -t "${tmp_file}"
  mv -- "${tmp_file}" "${output_file}"
}

ensure_silva_rrna_ref_db() {
  local gg_workspace_dir=$1
  local sys_file="/usr/local/db/silva/rRNA_ref.fa.gz"
  local runtime_root
  runtime_root=$(workspace_downloads_root "${gg_workspace_dir}")
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
  local output_db="${db_dir}/${uniref_db}_DB"
  local attempt
  local max_attempts=3
  mkdir -p "${db_dir}"
  for attempt in 1 2 3; do
    echo "Preparing MMseqs2 UniRef90 taxonomy DB in: ${db_dir} (attempt ${attempt}/${max_attempts})" >&2
    if mmseqs databases "${uniref_db}" "${output_db}" "${db_dir}" --threads "${nthreads}"; then
      break
    fi
    echo "MMseqs2 UniRef90 taxonomy DB preparation failed in: ${db_dir} (attempt ${attempt}/${max_attempts})" >&2
    if command -v df >/dev/null 2>&1; then
      df -h "${db_dir}" >&2 || true
    fi
    if [[ ${attempt} -lt ${max_attempts} ]]; then
      sleep 5
    else
      return 1
    fi
  done
  if [[ ! -s "${db_dir}/UniRef90_DB" || ! -s "${db_dir}/UniRef90_DB.dbtype" ]]; then
    return 1
  fi
  gg_write_ready_marker "${db_dir}/UniRef90_DB.ready"
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
    gg_write_ready_marker "${done_file}"
    return 0
  fi
  gg_array_download_once "${lock_file}" "${done_file}" "MMseqs2 UniRef90 taxonomy DB" \
    _download_mmseqs_uniref90_db "${db_dir}" "${nthreads}" || return 1

  if [[ ! -s "${db_file}" || ! -s "${dbtype_file}" || ! -s "${done_file}" ]]; then
    echo "MMseqs2 UniRef90 DB download failed." >&2
    return 1
  fi
}
