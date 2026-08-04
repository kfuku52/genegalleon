# shellcheck shell=bash
# Workflow execution, environment, version, and scheduler helpers.
# This file is sourced by workflow/support/gg_util.sh.

disable_flag_with_reason() {
  local variable_name=$1
  local reason=$2
  local variable_value="${!variable_name:-0}"
  if [[ "${variable_value}" == "0" ]]; then
    return 0
  fi
  printf -v "${variable_name}" '%s' 0
  echo "Disabled ${variable_name}: ${reason}"
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
  local file
  shopt -s nullglob dotglob
  for file in "${species_dir}"/*; do
    [[ -f "${file}" ]] || continue
    [[ "$(basename "${file}")" == .* ]] && continue
    files+=( "${file}" )
  done
  shopt -u nullglob dotglob
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
    printf 'Species names: %s\n' "${sp_names[*]}"
    exit 1
  else
    echo "Species files are unique in: ${species_dir}"
    printf 'Species names: %s\n' "${sp_names[*]}"
  fi
}

is_fastq_requiring_downstream_analysis_done() {
  local run_assembly_local=${run_assembly:-0}
  local run_merge_local=${run_amalgkit_merge:-0}
  local out=0

  if [[ ( ( ${run_assembly_local} -eq 1 && -s "${file_isoform}" ) || ${run_assembly_local} -eq 0 ) && \
    ( ( ${run_merge_local} -eq 1 && -s "${file_amalgkit_merge_count}" ) || ${run_merge_local} -eq 0 ) ]]; then
    out=1
  fi

  if [[ ${run_assembly_local} -eq 0 && ${run_merge_local} -eq 0 ]]; then
    out=0
  fi
  echo "${out}"
}

get_total_fastq_len() {
  local input_dir=$1
  local regex=$2
  local files=()
  local f
  while IFS= read -r -d '' f; do
    files+=( "${f}" )
  done < <(find "${input_dir}" -type f ! -name '.*' -name "${regex}" -print0)
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
  local gg_workspace_dir_local=$1
  local conda_env=${2:-}
  local set_unlimited_stack=${3:-1}
  local print_start_message=${4:-1}

	gg_initialize_data_layout "${gg_workspace_dir_local}"
	gg_workspace_layout_resolved=$(gg_resolve_workspace_layout "${gg_workspace_dir_local}")
	gg_workspace_input_dir=$(workspace_input_root "${gg_workspace_dir_local}")
	gg_workspace_output_dir=$(workspace_output_root "${gg_workspace_dir_local}")
	gg_workspace_downloads_dir=$(workspace_downloads_root "${gg_workspace_dir_local}")
	export gg_workspace_layout_resolved gg_workspace_input_dir gg_workspace_output_dir gg_workspace_downloads_dir

	if [[ -n "${conda_env}" ]]; then
		gg_activate_conda_env "${conda_env}"
	fi
	if [[ "${set_unlimited_stack}" -eq 1 ]]; then
		ulimit -s unlimited 2>/dev/null || true
	fi
	if [[ "${print_start_message}" -eq 1 ]]; then
		print_gg_container_starting_message
	fi
}

gg_initialize_conda_shell() {
	# Conda shell functions are local to the current Bash process. Keep the
	# initialization marker shell-local too, and do not trust a marker inherited
	# from a parent process unless the function itself is also present.
	export -n GG_CONDA_SHELL_INITIALIZED 2>/dev/null || true
	if declare -F conda >/dev/null 2>&1; then
		GG_CONDA_SHELL_INITIALIZED=1
		return 0
	fi
	GG_CONDA_SHELL_INITIALIZED=0
	if command -v micromamba >/dev/null 2>&1; then
		export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-/opt/conda}"
		eval "$(micromamba shell hook --shell bash)"
		if ! declare -F conda >/dev/null 2>&1; then
			conda() {
				micromamba "$@"
			}
		fi
	elif [[ -f /opt/conda/etc/profile.d/conda.sh ]]; then
		# shellcheck disable=SC1091
		source /opt/conda/etc/profile.d/conda.sh
	elif command -v conda >/dev/null 2>&1; then
		eval "$(conda shell.bash hook)"
	else
		return 1
	fi
	if ! declare -F conda >/dev/null 2>&1; then
		return 1
	fi
	GG_CONDA_SHELL_INITIALIZED=1
	return 0
}

gg_activate_conda_env() {
	local conda_env=${1:-}
	local had_nounset=0
	local activate_status=0
	if [[ -z "${conda_env}" ]]; then
		return 0
	fi
	if ! gg_initialize_conda_shell; then
		echo "gg_activate_conda_env: failed to initialize conda shell support." >&2
		return 1
	fi
	if ! command -v conda >/dev/null 2>&1; then
		echo "gg_activate_conda_env: conda command is unavailable after initialization." >&2
		return 1
	fi
	if [[ $- == *u* ]]; then
		had_nounset=1
		set +u
	fi
	conda activate "${conda_env}"
	activate_status=$?
	if [[ ${had_nounset} -eq 1 ]]; then
		set -u
	fi
	return "${activate_status}"
}

gg_deactivate_conda_env() {
	local had_nounset=0
	if [[ "${GG_CONDA_SHELL_INITIALIZED:-0}" -ne 1 ]]; then
		return 0
	fi
	if declare -F conda >/dev/null 2>&1; then
		if [[ $- == *u* ]]; then
			had_nounset=1
			set +u
		fi
		conda deactivate >/dev/null 2>&1 || true
		if [[ ${had_nounset} -eq 1 ]]; then
			set -u
		fi
	fi
}

gg_test_r_packages() {
  local rpackage
  for rpackage in "$@"; do
    echo "Testing: ${rpackage}"
    echo "Testing: ${rpackage}" >&2
    R -q -e "suppressPackageStartupMessages(library(${rpackage}, quietly=TRUE))" > /dev/null
  done
}

gg_test_python_packages() {
  local pypackage
  for pypackage in "$@"; do
    echo "Testing: ${pypackage}"
    echo "Testing: ${pypackage}" >&2
    python -c 'import importlib, sys; importlib.import_module(sys.argv[1])' "${pypackage}" > /dev/null
  done
}

gg_test_shell_commands() {
  local command_text
  local -a command_parts=()
  for command_text in "$@"; do
    echo "Testing: ${command_text}"
    command_parts=()
    read -r -a command_parts <<< "${command_text}"
    if [[ ${#command_parts[@]} -eq 0 ]]; then
      continue
    fi
    "${command_parts[@]}" > /dev/null
  done
}

gg_step_start() {
  local task_name=$1
  echo "$(date): Start: ${task_name}"
  echo "$(date): Start: ${task_name}" >&2
}

gg_step_skip() {
  local task_name=$1
  echo "$(date): Skipped: ${task_name}"
}

gg_count_fasta_records() {
  local infile=$1
  if [[ ! -s "${infile}" ]]; then
    echo 0
    return 0
  fi
  seqkit seq --threads 1 "${infile}" | awk '/^>/{n++} END{print n+0}'
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

gg_read_repo_version() {
  local version_file=${1:-}
  local version="unknown"

  if [[ -n "${version_file}" && -s "${version_file}" ]]; then
    version="$(head -n 1 "${version_file}" | tr -d '\r')"
  fi
  echo "${version}"
}

gg_trim_text() {
  printf '%s' "${1:-}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//'
}

gg_unquote_text() {
  local value=""

  value="$(gg_trim_text "${1:-}")"
  case "${value}" in
    \"*\")
      value="${value#\"}"
      value="${value%\"}"
      ;;
    \'*\')
      value="${value#\'}"
      value="${value%\'}"
      ;;
  esac
  printf '%s\n' "${value}"
}

gg_extract_inspect_value() {
  local inspect_text=${1:-}
  local key=${2:-}
  local line=""

  if [[ -z "${inspect_text}" || -z "${key}" ]]; then
    return 1
  fi

  while IFS= read -r line; do
    line="${line%$'\r'}"
    line="$(gg_trim_text "${line}")"
    case "${line}" in
      "${key}:"*)
        gg_unquote_text "${line#*:}"
        return 0
        ;;
      "${key}="*)
        gg_unquote_text "${line#*=}"
        return 0
        ;;
      "${key} "*)
        gg_unquote_text "${line#"${key}"}"
        return 0
        ;;
    esac
  done <<EOF
${inspect_text}
EOF

  return 1
}

gg_guess_container_tag_from_ref() {
  local image_ref=${1:-}
  local last_component=""

  if [[ -z "${image_ref}" ]]; then
    return 1
  fi
  last_component="${image_ref##*/}"
  if [[ "${last_component}" == *@* ]]; then
    return 1
  fi
  if [[ "${last_component}" != *:* ]]; then
    return 1
  fi
  printf '%s\n' "${last_component##*:}"
}

gg_print_version_summary() {
  local section_title=${1:-}
  local container_runtime_bin=${2:-}
  local inspect_text=${3:-}
  local version_file=""
  local gg_version=""
  local container_image_path=""
  local container_version_display=""
  local image_version=""
  local image_version_source="unknown"
  local image_revision=""
  local image_ref=""
  local image_tag=""

  if [[ -n "${gg_workflow_dir:-}" ]]; then
    version_file="${gg_workflow_dir}/../VERSION"
  fi
  gg_version="${SINGULARITYENV_GG_VERSION:-${APPTAINERENV_GG_VERSION:-}}"
  if [[ -z "${gg_version}" ]]; then
    gg_version="$(gg_read_repo_version "${version_file}")"
  fi

  container_image_path="${gg_container_image_path:-unknown}"
  if [[ -z "${inspect_text}" && -n "${container_runtime_bin}" && -n "${gg_container_image_path:-}" ]]; then
    if command -v "${container_runtime_bin}" >/dev/null 2>&1; then
      inspect_text=$("${container_runtime_bin}" inspect "${gg_container_image_path}" 2>/dev/null || true)
    fi
  fi

  image_version="$(gg_extract_inspect_value "${inspect_text}" "org.opencontainers.image.version" || true)"
  if [[ -n "${image_version}" ]]; then
    image_version_source="label"
  fi
  image_revision="$(gg_extract_inspect_value "${inspect_text}" "org.opencontainers.image.revision" || true)"
  image_ref="$(gg_extract_inspect_value "${inspect_text}" "io.genegalleon.local_image_ref" || true)"
  if [[ -z "${image_ref}" ]]; then
    image_ref="$(gg_extract_inspect_value "${inspect_text}" "docker-image" || true)"
  fi
  image_tag="$(gg_extract_inspect_value "${inspect_text}" "io.genegalleon.local_image_tag" || true)"
  if [[ -z "${image_tag}" && -n "${image_ref}" ]]; then
    image_tag="$(gg_guess_container_tag_from_ref "${image_ref}" || true)"
  fi
  if [[ -z "${image_version}" && -n "${image_tag}" ]]; then
    image_version="${image_tag}"
    image_version_source="tag"
  fi
  if [[ -z "${image_version}" && -n "${image_revision}" ]]; then
    image_version="revision:${image_revision}"
    image_version_source="revision"
  fi
  if [[ -z "${image_version}" ]]; then
    image_version="unknown"
  fi
  container_version_display="${image_version}"
  if [[ -n "${image_tag}" && "${image_tag}" != "${image_version}" ]]; then
    container_version_display="${container_version_display} (${image_tag})"
  fi

  if [[ -n "${section_title}" ]]; then
    gg_print_section "${section_title}"
  fi
  echo "genegalleon version: ${gg_version}"
  echo "container version: ${container_version_display}"
  if [[ "${image_version_source}" == "label" && "${gg_version}" != "unknown" && "${image_version}" != "unknown" && "${gg_version}" != "${image_version}" ]]; then
    echo "WARNING: genegalleon version (${gg_version}) does not match container version (${image_version})."
  fi
  echo "container path: ${container_image_path}"
  if [[ -n "${image_ref}" ]]; then
    echo "container image ref: ${image_ref}"
  fi
  if [[ -n "${image_revision}" ]]; then
    echo "container revision: ${image_revision}"
  fi
  if [[ -n "${section_title}" ]]; then
    gg_print_spacer
  fi
  return 0
}

gg_entrypoint_print_version_summary() {
  local container_runtime_bin=""

  container_runtime_bin="$(gg_container_shell_command_runtime_bin || true)"
  gg_print_version_summary "genegalleon startup versions" "${container_runtime_bin}" || true
  return 0
}
