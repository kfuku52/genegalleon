#!/usr/bin/env bash
set -euo pipefail

copy_busco_tables() {
  local busco_root_dir="$1"
  local lineage="$2"
  local file_full="$3"
  local file_short="$4"
  local run_dir="${busco_root_dir}/run_${lineage}"
  local full_src=""
  local short_src=""
  local tmp_full_src=""
  local -a short_candidates=()

  if [[ -s "${run_dir}/full_table.tsv" ]]; then
    full_src="${run_dir}/full_table.tsv"
  elif [[ -s "${run_dir}/full_table.tsv.gz" ]]; then
    tmp_full_src=$(mktemp "${TMPDIR:-/tmp}/gg_busco.full_table.XXXXXX.tsv")
    gzip -cd "${run_dir}/full_table.tsv.gz" > "${tmp_full_src}"
    full_src="${tmp_full_src}"
  fi

  if [[ -s "${run_dir}/short_summary.txt" ]]; then
    short_src="${run_dir}/short_summary.txt"
  else
    while IFS= read -r short_candidate; do
      [[ -n "${short_candidate}" ]] || continue
      short_candidates+=( "${short_candidate}" )
    done < <(find "${run_dir}" -maxdepth 1 -type f -name "short_summary*.txt" | sort)
    if [[ ${#short_candidates[@]} -gt 0 ]]; then
      short_src="${short_candidates[0]}"
    fi
  fi

  if [[ -z "${full_src}" || -z "${short_src}" ]]; then
    echo "BUSCO outputs were not found under ${run_dir}."
    echo "Expected full_table.tsv(.gz) and short_summary*.txt."
    rm -f -- "${tmp_full_src}"
    return 1
  fi

  cp_out "${full_src}" "${file_full}"
  cp_out "${short_src}" "${file_short}"
  rm -f -- "${tmp_full_src}"
}

get_busco_summary_gene_count() {
  local summary_table=$1
  local line_count=0
  if [[ ! -s "${summary_table}" ]]; then
    echo 0
    return 0
  fi
  line_count=$(wc -l < "${summary_table}")
  if [[ ${line_count} -le 1 ]]; then
    echo 0
  else
    echo $((line_count - 1))
  fi
}

normalize_busco_table_naming() {
  local full_dir=$1
  local short_dir=$2
  if [[ -n "${full_dir}" ]]; then
    ensure_dir "${full_dir}"
  fi
  if [[ -n "${short_dir}" ]]; then
    ensure_dir "${short_dir}"
  fi
}

busco_output_exists_for_species() {
  local search_dir=$1
  local species_name=$2
  local name_glob=$3
  local busco_file busco_base busco_species

  if [[ -z "${search_dir}" || ! -d "${search_dir}" ]]; then
    return 1
  fi

  while IFS= read -r busco_file; do
    busco_base=$(basename "${busco_file}")
    busco_species=$(gg_species_name_from_path_or_dot "${busco_base}")
    if [[ "${busco_species}" == "${species_name}" ]]; then
      return 0
    fi
  done < <(find "${search_dir}" -maxdepth 1 -type f ! -name '.*' -name "${name_glob}" 2> /dev/null | sort)

  return 1
}

remove_busco_outputs_for_species() {
  local search_dir=$1
  local species_name=$2
  local name_glob=$3
  local busco_file busco_base busco_species

  if [[ -z "${search_dir}" || ! -d "${search_dir}" ]]; then
    return 0
  fi

  while IFS= read -r busco_file; do
    busco_base=$(basename "${busco_file}")
    busco_species=$(gg_species_name_from_path_or_dot "${busco_base}")
    if [[ "${busco_species}" == "${species_name}" ]]; then
      rm -f -- "${busco_file}"
    fi
  done < <(find "${search_dir}" -maxdepth 1 -type f ! -name '.*' -name "${name_glob}" 2> /dev/null | sort)
}
