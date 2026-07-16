# shellcheck shell=bash
# Reference database download and preparation helpers.
# This file is sourced by workflow/support/gg_util.sh.

_download_uniprot_sprot_to_prefix() {
  local output_prefix=$1
  local output_dir
  output_dir=$(dirname "${output_prefix}")
  local tmp_dir
  tmp_dir=$(mktemp -d "${output_dir}/tmp.uniprot_sprot.XXXXXX")
  local pep_tmp="${tmp_dir}/uniprot_sprot.pep"
  local dmnd_tmp_prefix="${tmp_dir}/uniprot_sprot"
  local uniprot_url="${GG_UNIPROT_SPROT_URL:-https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz}"

  if ! curl -fsSL "${uniprot_url}" | gzip -dc > "${pep_tmp}"; then
    echo "Failed to download/decompress UniProt Swiss-Prot FASTA from: ${uniprot_url}" >&2
    rm -rf -- "${tmp_dir}"
    return 1
  fi
  if [[ ! -s "${pep_tmp}" ]]; then
    echo "Failed to download UniProt Swiss-Prot FASTA from: ${uniprot_url}" >&2
    rm -rf -- "${tmp_dir}"
    return 1
  fi

  if ! diamond makedb --in "${pep_tmp}" --db "${dmnd_tmp_prefix}"; then
    echo "Failed to run DIAMOND makedb for UniProt Swiss-Prot." >&2
    rm -rf -- "${tmp_dir}"
    return 1
  fi
  if [[ ! -s "${dmnd_tmp_prefix}.dmnd" ]]; then
    echo "Failed to build DIAMOND DB from UniProt Swiss-Prot." >&2
    rm -rf -- "${tmp_dir}"
    return 1
  fi

  mv -- "${pep_tmp}" "${output_prefix}.pep"
  mv -- "${dmnd_tmp_prefix}.dmnd" "${output_prefix}.dmnd"
  rm -rf -- "${tmp_dir}"
}

_download_uniprot_sprot_pep_to_file() {
  local output_file=$1
  local output_dir
  output_dir=$(dirname "${output_file}")
  local tmp_dir
  tmp_dir=$(mktemp -d "${output_dir}/tmp.uniprot_sprot_pep.XXXXXX")
  local pep_tmp="${tmp_dir}/uniprot_sprot.pep"
  local uniprot_url="${GG_UNIPROT_SPROT_URL:-https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz}"

  if ! curl -fsSL "${uniprot_url}" | gzip -dc > "${pep_tmp}"; then
    echo "Failed to download/decompress UniProt Swiss-Prot FASTA from: ${uniprot_url}" >&2
    rm -rf -- "${tmp_dir}"
    return 1
  fi
  if [[ ! -s "${pep_tmp}" ]]; then
    echo "Failed to download UniProt Swiss-Prot FASTA from: ${uniprot_url}" >&2
    rm -rf -- "${tmp_dir}"
    return 1
  fi

  mv -- "${pep_tmp}" "${output_file}"
  rm -rf -- "${tmp_dir}"
}

_download_uniprot_sprot_dat_to_file() {
  local output_file=$1
  local output_dir
  output_dir=$(dirname "${output_file}")
  local tmp_dir
  tmp_dir=$(mktemp -d "${output_dir}/tmp.uniprot_sprot_dat.XXXXXX")
  local dat_tmp="${tmp_dir}/uniprot_sprot.dat.gz"
  local uniprot_dat_url="${GG_UNIPROT_SPROT_DAT_URL:-https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz}"

  if ! curl -fsSL "${uniprot_dat_url}" -o "${dat_tmp}"; then
    echo "Failed to download UniProt Swiss-Prot DAT from: ${uniprot_dat_url}" >&2
    rm -rf -- "${tmp_dir}"
    return 1
  fi
  if [[ ! -s "${dat_tmp}" ]]; then
    echo "Downloaded UniProt Swiss-Prot DAT is empty: ${uniprot_dat_url}" >&2
    rm -rf -- "${tmp_dir}"
    return 1
  fi

  mv -- "${dat_tmp}" "${output_file}"
  rm -rf -- "${tmp_dir}"
}

_build_uniprot_sprot_metadata_from_dat() {
  local dat_path=$1
  local out_path=$2
  local py_exec=""
  py_exec=$(gg_find_python_exec || true)
  if [[ -z "${py_exec}" ]]; then
    echo "python/python3 command was not found. Cannot build UniProt metadata." >&2
    return 1
  fi
  local util_dir="${GG_UTIL_SUPPORT_DIR:-}"
  if [[ -z "${util_dir}" ]]; then
    util_dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)"
  fi
  local script_path="${util_dir}/build_uniprot_metadata_from_dat.py"
  if [[ ! -s "${script_path}" ]]; then
    echo "UniProt metadata builder script was not found: ${script_path}" >&2
    return 1
  fi
  local out_tmp="${out_path}.tmp.$$"
  if [[ "${out_path}" == *.gz ]]; then
    out_tmp="${out_tmp}.gz"
  fi
  if ! "${py_exec}" "${script_path}" --uniprot_dat "${dat_path}" --outfile "${out_tmp}"; then
    rm -f -- "${out_tmp}"
    return 1
  fi
  if [[ ! -s "${out_tmp}" ]]; then
    echo "UniProt metadata builder produced an empty file: ${out_tmp}" >&2
    rm -f -- "${out_tmp}"
    return 1
  fi
  mv -- "${out_tmp}" "${out_path}"
}

_prepare_uniprot_sprot_meta_to_file() {
  local output_meta=$1
  local sys_prefix=$2
  local runtime_prefix=$3

  if [[ -s "${output_meta}" ]]; then
    return 0
  fi
  if [[ -s "${sys_prefix}.meta.tsv.gz" ]]; then
    cp -- "${sys_prefix}.meta.tsv.gz" "${output_meta}"
    return 0
  fi
  if [[ -s "${runtime_prefix}.meta.tsv.gz" && "${runtime_prefix}.meta.tsv.gz" != "${output_meta}" ]]; then
    cp -- "${runtime_prefix}.meta.tsv.gz" "${output_meta}"
    return 0
  fi

  local tmp_dir=""
  tmp_dir=$(mktemp -d "$(dirname "${output_meta}")/tmp.uniprot_sprot_meta.XXXXXX")
  local dat_source=""

  if [[ -s "${runtime_prefix}.dat.gz" ]]; then
    dat_source="${runtime_prefix}.dat.gz"
  elif [[ -s "${sys_prefix}.dat.gz" ]]; then
    dat_source="${sys_prefix}.dat.gz"
  else
    dat_source="${tmp_dir}/uniprot_sprot.dat.gz"
    if ! _download_uniprot_sprot_dat_to_file "${dat_source}"; then
      rm -rf -- "${tmp_dir}"
      return 1
    fi
    if [[ ! -s "${runtime_prefix}.dat.gz" ]]; then
      cp -- "${dat_source}" "${runtime_prefix}.dat.gz" || true
    fi
  fi

  if ! _build_uniprot_sprot_metadata_from_dat "${dat_source}" "${output_meta}"; then
    rm -rf -- "${tmp_dir}"
    return 1
  fi
  rm -rf -- "${tmp_dir}"
}

ensure_uniprot_sprot_metadata_tsv() {
  local gg_workspace_dir=$1
  local db_prefix="${2:-}"
  local sys_prefix="/usr/local/db/uniprot_sprot"
  local runtime_root
  runtime_root="$(workspace_downloads_root "${gg_workspace_dir}")"
  local runtime_prefix="${runtime_root}/uniprot_sprot/uniprot_sprot"
  local runtime_meta="${runtime_prefix}.meta.tsv.gz"
  local lock_file="${runtime_root}/locks/uniprot_sprot.meta.lock"
  local target_meta="${runtime_meta}"

  if [[ -n "${db_prefix}" && -s "${db_prefix}.meta.tsv.gz" ]]; then
    echo "${db_prefix}.meta.tsv.gz"
    return 0
  fi
  if [[ -n "${db_prefix}" && "${db_prefix}" == "${sys_prefix}" && -s "${sys_prefix}.meta.tsv.gz" ]]; then
    echo "${sys_prefix}.meta.tsv.gz"
    return 0
  fi
  if [[ -n "${db_prefix}" && "${db_prefix}" == "${runtime_prefix}" ]]; then
    target_meta="${runtime_meta}"
  fi
  if [[ -s "${target_meta}" ]]; then
    echo "${target_meta}"
    return 0
  fi

  mkdir -p "${runtime_root}" "$(dirname "${runtime_meta}")"
  gg_array_download_once "${lock_file}" "${target_meta}" "UniProt Swiss-Prot metadata TSV" \
    _prepare_uniprot_sprot_meta_to_file "${target_meta}" "${sys_prefix}" "${runtime_prefix}" || return 1

  if [[ ! -s "${target_meta}" ]]; then
    echo "UniProt Swiss-Prot metadata generation failed: ${target_meta}" >&2
    return 1
  fi
  echo "${target_meta}"
}

_build_mmseqs_db_from_fasta() {
  local fasta_file=$1
  local db_prefix=$2
  local done_file=$3

  if ! mmseqs createdb "${fasta_file}" "${db_prefix}"; then
    echo "Failed to build MMseqs2 DB from FASTA: ${fasta_file}" >&2
    return 1
  fi
  if [[ ! -s "${db_prefix}" || ! -s "${db_prefix}.dbtype" ]]; then
    echo "MMseqs2 DB files were not generated: ${db_prefix}" >&2
    return 1
  fi
  gg_write_ready_marker "${done_file}"
}

_build_blast_db_from_fasta() {
  local fasta_file=$1
  local db_prefix=$2
  local done_file=$3

  if ! makeblastdb -dbtype prot -in "${fasta_file}" -out "${db_prefix}"; then
    echo "Failed to build BLASTP DB from FASTA: ${fasta_file}" >&2
    return 1
  fi
  if [[ ! -s "${db_prefix}.pin" || ! -s "${db_prefix}.phr" || ! -s "${db_prefix}.psq" ]]; then
    echo "BLASTP DB files were not generated: ${db_prefix}" >&2
    return 1
  fi
  gg_write_ready_marker "${done_file}"
}

validate_uniprot_sprot_db_prefix() {
  local db_prefix=$1
  local db_kind=$2

  if [[ -z "${db_prefix}" ]]; then
    echo "UniProt Swiss-Prot ${db_kind} DB prefix is empty." >&2
    return 1
  fi
  if [[ "${db_prefix}" == *$'\n'* || "${db_prefix}" == *$'\r'* || "${db_prefix}" == *$'\t'* ]]; then
    printf 'UniProt Swiss-Prot %s DB prefix is malformed: %q\n' "${db_kind}" "${db_prefix}" >&2
    return 1
  fi

  case "${db_kind}" in
    blastp)
      if [[ ! -s "${db_prefix}.pin" || ! -s "${db_prefix}.phr" || ! -s "${db_prefix}.psq" ]]; then
        echo "UniProt Swiss-Prot BLASTP DB files were not found for prefix: ${db_prefix}" >&2
        return 1
      fi
      ;;
    mmseqs2)
      if [[ ! -s "${db_prefix}.mmseqs" || ! -s "${db_prefix}.mmseqs.dbtype" ]]; then
        echo "UniProt Swiss-Prot MMseqs2 DB files were not found for prefix: ${db_prefix}" >&2
        return 1
      fi
      ;;
    *)
      echo "Unsupported UniProt Swiss-Prot DB kind: ${db_kind}" >&2
      return 1
      ;;
  esac
}

ensure_uniprot_sprot_db() {
  local gg_workspace_dir=$1
  local sys_prefix="/usr/local/db/uniprot_sprot"
  local runtime_root
  runtime_root=$(workspace_downloads_root "${gg_workspace_dir}")
  local runtime_dir="${runtime_root}/uniprot_sprot"
  local runtime_prefix="${runtime_dir}/uniprot_sprot"
  local lock_file="${runtime_root}/locks/uniprot_sprot.lock"

  if [[ -s "${sys_prefix}.pep" && -s "${sys_prefix}.dmnd" ]]; then
    if ! ensure_uniprot_sprot_metadata_tsv "${gg_workspace_dir}" "${sys_prefix}" >/dev/null; then
      echo "Warning: UniProt Swiss-Prot metadata TSV is unavailable for prefix: ${sys_prefix}" >&2
    fi
    echo "${sys_prefix}"
    return 0
  fi
  if [[ -s "${runtime_prefix}.pep" && -s "${runtime_prefix}.dmnd" ]]; then
    if ! ensure_uniprot_sprot_metadata_tsv "${gg_workspace_dir}" "${runtime_prefix}" >/dev/null; then
      echo "Warning: UniProt Swiss-Prot metadata TSV is unavailable for prefix: ${runtime_prefix}" >&2
    fi
    echo "${runtime_prefix}"
    return 0
  fi

  mkdir -p "${runtime_root}" "${runtime_dir}"
  gg_array_download_once "${lock_file}" "${runtime_prefix}.dmnd" "UniProt Swiss-Prot (FASTA + DIAMOND)" \
    _download_uniprot_sprot_to_prefix "${runtime_prefix}" || return 1

  if [[ ! -s "${runtime_prefix}.pep" || ! -s "${runtime_prefix}.dmnd" ]]; then
    echo "UniProt DB download/build failed." >&2
    return 1
  fi
  if ! ensure_uniprot_sprot_metadata_tsv "${gg_workspace_dir}" "${runtime_prefix}" >/dev/null; then
    echo "Warning: UniProt Swiss-Prot metadata TSV is unavailable for prefix: ${runtime_prefix}" >&2
  fi
  echo "${runtime_prefix}"
}

ensure_uniprot_sprot_mmseqs_db() {
  local gg_workspace_dir=$1
  local sys_prefix="/usr/local/db/uniprot_sprot"
  local runtime_root
  runtime_root=$(workspace_downloads_root "${gg_workspace_dir}")
  local runtime_dir="${runtime_root}/uniprot_sprot"
  local runtime_prefix="${runtime_dir}/uniprot_sprot"
  local runtime_pep="${runtime_prefix}.pep"
  local runtime_mmseqs_db="${runtime_prefix}.mmseqs"
  local runtime_mmseqs_dbtype="${runtime_mmseqs_db}.dbtype"
  local runtime_ready="${runtime_mmseqs_db}.ready"
  local lock_file_pep="${runtime_root}/locks/uniprot_sprot.pep.lock"
  local lock_file_mmseqs="${runtime_root}/locks/uniprot_sprot.mmseqs.lock"

  if [[ -s "${sys_prefix}.pep" && -s "${sys_prefix}.mmseqs" && -s "${sys_prefix}.mmseqs.dbtype" ]]; then
    if ! ensure_uniprot_sprot_metadata_tsv "${gg_workspace_dir}" "${sys_prefix}" >/dev/null; then
      echo "Warning: UniProt Swiss-Prot metadata TSV is unavailable for prefix: ${sys_prefix}" >&2
    fi
    echo "${sys_prefix}"
    return 0
  fi

  mkdir -p "${runtime_root}" "${runtime_dir}"

  if [[ ! -s "${runtime_pep}" ]]; then
    if [[ -s "${sys_prefix}.pep" ]]; then
      gg_array_download_once "${lock_file_pep}" "${runtime_pep}" "UniProt Swiss-Prot FASTA" \
        cp -- "${sys_prefix}.pep" "${runtime_pep}" || return 1
    else
      gg_array_download_once "${lock_file_pep}" "${runtime_pep}" "UniProt Swiss-Prot FASTA" \
        _download_uniprot_sprot_pep_to_file "${runtime_pep}" || return 1
    fi
  fi
  if [[ ! -s "${runtime_pep}" ]]; then
    echo "UniProt Swiss-Prot FASTA was not found after synchronization: ${runtime_pep}" >&2
    return 1
  fi

  if [[ -s "${runtime_mmseqs_db}" && -s "${runtime_mmseqs_dbtype}" && -s "${runtime_ready}" ]]; then
    echo "${runtime_prefix}"
    return 0
  fi
  if [[ -s "${runtime_mmseqs_db}" && -s "${runtime_mmseqs_dbtype}" && ! -s "${runtime_ready}" ]]; then
    echo "MMseqs2 UniProt Swiss-Prot DB exists without ready marker. Reusing existing DB and creating ready marker." >&2
    gg_write_ready_marker "${runtime_ready}"
    echo "${runtime_prefix}"
    return 0
  fi

  gg_array_download_once "${lock_file_mmseqs}" "${runtime_ready}" "UniProt Swiss-Prot MMseqs2 DB" \
    _build_mmseqs_db_from_fasta "${runtime_pep}" "${runtime_mmseqs_db}" "${runtime_ready}" || return 1

  if [[ ! -s "${runtime_mmseqs_db}" || ! -s "${runtime_mmseqs_dbtype}" || ! -s "${runtime_ready}" ]]; then
    echo "MMseqs2 UniProt Swiss-Prot DB download/build failed." >&2
    return 1
  fi
  if ! ensure_uniprot_sprot_metadata_tsv "${gg_workspace_dir}" "${runtime_prefix}" >/dev/null; then
    echo "Warning: UniProt Swiss-Prot metadata TSV is unavailable for prefix: ${runtime_prefix}" >&2
  fi
  echo "${runtime_prefix}"
}

ensure_uniprot_sprot_blast_db() {
  local gg_workspace_dir=$1
  local sys_prefix="/usr/local/db/uniprot_sprot"
  local runtime_root
  runtime_root=$(workspace_downloads_root "${gg_workspace_dir}")
  local runtime_dir="${runtime_root}/uniprot_sprot"
  local runtime_prefix="${runtime_dir}/uniprot_sprot"
  local runtime_pep="${runtime_prefix}.pep"
  local runtime_blast_ready="${runtime_prefix}.blast.ready"
  local lock_file_pep="${runtime_root}/locks/uniprot_sprot.pep.lock"
  local lock_file_blast="${runtime_root}/locks/uniprot_sprot.blast.lock"

  if [[ -s "${sys_prefix}.pep" && -s "${sys_prefix}.pin" && -s "${sys_prefix}.phr" && -s "${sys_prefix}.psq" ]]; then
    if ! ensure_uniprot_sprot_metadata_tsv "${gg_workspace_dir}" "${sys_prefix}" >/dev/null; then
      echo "Warning: UniProt Swiss-Prot metadata TSV is unavailable for prefix: ${sys_prefix}" >&2
    fi
    echo "${sys_prefix}"
    return 0
  fi
  if [[ -s "${runtime_pep}" && -s "${runtime_prefix}.pin" && -s "${runtime_prefix}.phr" && -s "${runtime_prefix}.psq" ]]; then
    if ! ensure_uniprot_sprot_metadata_tsv "${gg_workspace_dir}" "${runtime_prefix}" >/dev/null; then
      echo "Warning: UniProt Swiss-Prot metadata TSV is unavailable for prefix: ${runtime_prefix}" >&2
    fi
    echo "${runtime_prefix}"
    return 0
  fi

  mkdir -p "${runtime_root}" "${runtime_dir}"

  if [[ ! -s "${runtime_pep}" ]]; then
    if [[ -s "${sys_prefix}.pep" ]]; then
      gg_array_download_once "${lock_file_pep}" "${runtime_pep}" "UniProt Swiss-Prot FASTA" \
        cp -- "${sys_prefix}.pep" "${runtime_pep}" || return 1
    else
      gg_array_download_once "${lock_file_pep}" "${runtime_pep}" "UniProt Swiss-Prot FASTA" \
        _download_uniprot_sprot_pep_to_file "${runtime_pep}" || return 1
    fi
  fi
  if [[ ! -s "${runtime_pep}" ]]; then
    echo "UniProt Swiss-Prot FASTA was not found after synchronization: ${runtime_pep}" >&2
    return 1
  fi

  if [[ -s "${runtime_prefix}.pin" && -s "${runtime_prefix}.phr" && -s "${runtime_prefix}.psq" && -s "${runtime_blast_ready}" ]]; then
    echo "${runtime_prefix}"
    return 0
  fi
  if [[ -s "${runtime_prefix}.pin" && -s "${runtime_prefix}.phr" && -s "${runtime_prefix}.psq" && ! -s "${runtime_blast_ready}" ]]; then
    echo "BLASTP UniProt Swiss-Prot DB exists without ready marker. Reusing existing DB and creating ready marker." >&2
    gg_write_ready_marker "${runtime_blast_ready}"
    echo "${runtime_prefix}"
    return 0
  fi

  gg_array_download_once "${lock_file_blast}" "${runtime_blast_ready}" "UniProt Swiss-Prot BLASTP DB" \
    _build_blast_db_from_fasta "${runtime_pep}" "${runtime_prefix}" "${runtime_blast_ready}" || return 1

  if [[ ! -s "${runtime_prefix}.pin" || ! -s "${runtime_prefix}.phr" || ! -s "${runtime_prefix}.psq" || ! -s "${runtime_blast_ready}" ]]; then
    echo "BLASTP UniProt Swiss-Prot DB download/build failed." >&2
    return 1
  fi
  if ! ensure_uniprot_sprot_metadata_tsv "${gg_workspace_dir}" "${runtime_prefix}" >/dev/null; then
    echo "Warning: UniProt Swiss-Prot metadata TSV is unavailable for prefix: ${runtime_prefix}" >&2
  fi
  echo "${runtime_prefix}"
}

_download_pfam_le_to_dir() {
  local output_dir=$1
  local ready_file="${output_dir}/.pfam_le.ready"
  local parent_dir
  parent_dir=$(dirname "${output_dir}")
  local tmp_dir
  tmp_dir=$(mktemp -d "${parent_dir}/tmp.pfam_le.XXXXXX")
  local url="${GG_PFAM_LE_URL:-https://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/little_endian/Pfam_LE.tar.gz}"
  local archive_path="${tmp_dir}/Pfam_LE.tar.gz"
  local staged_dir="${tmp_dir}/Pfam_LE"
  local pfam_files=()
  local f

  if ! curl -fsSL "${url}" -o "${archive_path}"; then
    echo "Failed to download Pfam_LE archive: ${url}" >&2
    rm -rf -- "${tmp_dir}"
    return 1
  fi
  if ! tar -xzf "${archive_path}" -C "${tmp_dir}"; then
    echo "Failed to extract Pfam_LE archive: ${archive_path}" >&2
    rm -rf -- "${tmp_dir}"
    return 1
  fi
  if ! mkdir -p "${staged_dir}"; then
    echo "Failed to create staging directory: ${staged_dir}" >&2
    rm -rf -- "${tmp_dir}"
    return 1
  fi
  local pfam_path
  while IFS= read -r pfam_path; do
    pfam_files+=( "${pfam_path}" )
  done < <(find "${tmp_dir}" -type f \( -name "Pfam.*" -o -name "Pfam.pal" \))
  if [[ ${#pfam_files[@]} -eq 0 ]]; then
    echo "Pfam_LE archive did not contain Pfam DB files." >&2
    rm -rf -- "${tmp_dir}"
    return 1
  fi
  for f in "${pfam_files[@]}"; do
    cp -- "${f}" "${staged_dir}/"
  done

  # Old format: Pfam.loo/Pfam.rps, New split format: Pfam.pal + Pfam.00.* volumes.
  if [[ ! -s "${staged_dir}/Pfam.pal" && ( ! -s "${staged_dir}/Pfam.loo" || ! -s "${staged_dir}/Pfam.rps" ) ]]; then
    echo "Pfam_LE archive did not contain expected RPS-BLAST index files." >&2
    rm -rf -- "${tmp_dir}"
    return 1
  fi

  if [[ -z "${output_dir}" || "${output_dir}" == "/" ]]; then
    echo "Unsafe Pfam output directory: '${output_dir}'" >&2
    rm -rf -- "${tmp_dir}"
    return 1
  fi
  rm -rf -- "${output_dir}.tmp"
  mv -- "${staged_dir}" "${output_dir}.tmp"
  rm -rf -- "${output_dir}"
  mv -- "${output_dir}.tmp" "${output_dir}"
  gg_write_ready_marker "${ready_file}"
  rm -rf -- "${tmp_dir}"
}

pfam_le_db_is_ready() {
  local db_dir=$1
  [[ -s "${db_dir}/Pfam.pal" || ( -s "${db_dir}/Pfam.loo" && -s "${db_dir}/Pfam.rps" ) ]]
}
