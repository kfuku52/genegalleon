# shellcheck shell=bash
# Container, scheduler, configuration, and temporary-file runtime helpers.
# This file is sourced by workflow/support/gg_util.sh.

unset_singularity_envs() {
	unset SINGULARITY_BIND
	unset SINGULARITY_BINDPATH
	unset APPTAINER_BIND
	unset APPTAINER_BINDPATH
	unset GG_CONTAINER_BIND_MOUNTS
	unset SINGULARITYENV_GG_ARRAY_TASK_ID
	unset SINGULARITYENV_GG_TASK_CPUS
	unset SINGULARITYENV_GG_JOB_ID
	unset SINGULARITYENV_GG_MEM_PER_CPU_GB
	unset SINGULARITYENV_GG_MEM_TOTAL_GB
	unset SINGULARITYENV_GG_MEM_TOOL_RESERVE_GB
	unset SINGULARITYENV_GG_MEM_TOOL_GB
	unset SINGULARITYENV_OMP_NUM_THREADS
	unset SINGULARITYENV_OPENBLAS_NUM_THREADS
	unset SINGULARITYENV_MKL_NUM_THREADS
	unset SINGULARITYENV_NUMEXPR_NUM_THREADS
	unset SINGULARITYENV_OMP_THREAD_LIMIT
	unset SINGULARITYENV_KMP_ALL_THREADS
	unset SINGULARITYENV_KMP_DEVICE_THREAD_LIMIT
	unset SINGULARITYENV_KMP_TEAMS_THREAD_LIMIT
	unset SINGULARITYENV_SGE_TASK_ID
	unset SINGULARITYENV_NSLOTS
	unset SINGULARITYENV_JOB_ID
	unset SINGULARITYENV_MEM_PER_SLOT
	unset SINGULARITYENV_MEM_PER_HOST
	unset APPTAINERENV_GG_ARRAY_TASK_ID
	unset APPTAINERENV_GG_TASK_CPUS
	unset APPTAINERENV_GG_JOB_ID
	unset APPTAINERENV_GG_MEM_PER_CPU_GB
	unset APPTAINERENV_GG_MEM_TOTAL_GB
	unset APPTAINERENV_GG_MEM_TOOL_RESERVE_GB
	unset APPTAINERENV_GG_MEM_TOOL_GB
	unset APPTAINERENV_OMP_NUM_THREADS
	unset APPTAINERENV_OPENBLAS_NUM_THREADS
	unset APPTAINERENV_MKL_NUM_THREADS
	unset APPTAINERENV_NUMEXPR_NUM_THREADS
	unset APPTAINERENV_OMP_THREAD_LIMIT
	unset APPTAINERENV_KMP_ALL_THREADS
	unset APPTAINERENV_KMP_DEVICE_THREAD_LIMIT
	unset APPTAINERENV_KMP_TEAMS_THREAD_LIMIT
	unset APPTAINERENV_SGE_TASK_ID
	unset APPTAINERENV_NSLOTS
	unset APPTAINERENV_JOB_ID
	unset APPTAINERENV_MEM_PER_SLOT
	unset APPTAINERENV_MEM_PER_HOST
}

gg_scheduler_runtime_prelude() {
	if declare -F gg_site_scheduler_prelude >/dev/null 2>&1; then
		gg_site_scheduler_prelude || return 1
	elif [[ -n "${PBS_O_WORKDIR:-}" ]]; then
		cd "${PBS_O_WORKDIR}" || return 1
	fi
	ulimit -s unlimited 2>/dev/null || true
}

gg_resolve_physical_path() {
	local path=${1:-}
	local parent=""
	local resolved_path=""

	if [[ -z "${path}" ]]; then
		return 1
	fi
	if command -v readlink >/dev/null 2>&1; then
		if resolved_path=$(readlink -f -- "${path}" 2>/dev/null); then
			printf '%s\n' "${resolved_path}"
			return 0
		fi
	fi
	parent=$(cd "$(dirname "${path}")" && pwd -P) || return 1
	resolved_path="${parent}/$(basename -- "${path}")"
	printf '%s\n' "${resolved_path}"
}

gg_normalize_contamination_removal_rank_for_amalgkit() {
	local rank=${1:-}

	rank=$(printf '%s' "${rank}" | tr '[:upper:]' '[:lower:]')
	case "${rank}" in
		""|"superkingdom")
			printf '%s\n' "domain"
			;;
		*)
			printf '%s\n' "${rank}"
			;;
	esac
}

gg_normalize_contamination_removal_rank_for_remove_contaminated_sequences() {
	local rank=${1:-}

	rank=$(printf '%s' "${rank}" | tr '[:upper:]' '[:lower:]')
	case "${rank}" in
		""|"domain")
			printf '%s\n' "superkingdom"
			;;
		*)
			printf '%s\n' "${rank}"
			;;
	esac
}

gg_csv_prepend() {
	local item=$1
	local existing=${2:-}
	if [[ -z "${existing}" ]]; then
		printf '%s' "${item}"
	else
		printf '%s,%s' "${item}" "${existing}"
	fi
}

gg_csv_append() {
	local existing=${1:-}
	local item=${2:-}
	if [[ -z "${item}" ]]; then
		printf '%s' "${existing}"
	elif [[ -z "${existing}" ]]; then
		printf '%s' "${item}"
	else
		printf '%s,%s' "${existing}" "${item}"
	fi
}

gg_trim_ascii_whitespace() {
	local value=${1:-}
	value="${value#"${value%%[![:space:]]*}"}"
	value="${value%"${value##*[![:space:]]}"}"
	printf '%s' "${value}"
}

gg_is_valid_shell_var_name() {
	local var_name=${1:-}
	[[ "${var_name}" =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]]
}

gg_normalize_config_var_name() {
	local var_name=${1:-}
	var_name=$(gg_trim_ascii_whitespace "${var_name}")
	if [[ -z "${var_name}" || "${var_name}" == \#* ]]; then
		return 1
	fi
	printf '%s\n' "${var_name}"
}

gg_should_mask_config_var_name() {
	local var_name=${1:-}
	local var_name_upper=""

	var_name_upper=$(printf '%s' "${var_name}" | tr '[:lower:]' '[:upper:]')
	case "${var_name_upper}" in
		*TOKEN*|*SECRET*|*PASSWORD*|*PASSWD*|*API_KEY*|*HTTP_HEADER*|*AUTH_HEADER*|*PRIVATE_KEY*|*ACCESS_KEY*)
			return 0
			;;
	esac
	return 1
}

gg_print_named_config_summary() {
	local title=${1:-}
	shift || true
	local raw_var_name=""
	local var_name=""
	local value=""

	if [[ -z "${title}" ]]; then
		title="config summary"
	fi

	echo "### ${title} ###"
	for raw_var_name in "$@"; do
		var_name=$(gg_normalize_config_var_name "${raw_var_name}" || true)
		if [[ -z "${var_name}" ]]; then
			continue
		fi
		if ! gg_is_valid_shell_var_name "${var_name}"; then
			echo "${var_name}=<invalid-var-name>"
			continue
		fi
		if [[ -z "${!var_name+x}" ]]; then
			echo "${var_name}=<unset>"
			continue
		fi
		if gg_should_mask_config_var_name "${var_name}"; then
			echo "${var_name}=<masked>"
			continue
		fi
		value=${!var_name}
		if [[ -z "${value}" ]]; then
			value="<empty>"
		fi
		printf '%s=%s\n' "${var_name}" "${value}"
	done
	echo ""
}

gg_print_registered_config_summary() {
	local entrypoint_name=${1:-}
	local title=${2:-}
	shift 2 || true
	local raw_var_name=""
	local var_name=""
	local -a summary_vars=()

	if [[ -z "${entrypoint_name}" ]]; then
		return 0
	fi
	if ! declare -F gg_print_entrypoint_config_vars >/dev/null 2>&1; then
		echo "gg_print_registered_config_summary: config var registry helper is unavailable." >&2
		return 0
	fi

	while IFS= read -r raw_var_name; do
		var_name=$(gg_normalize_config_var_name "${raw_var_name}" || true)
		if [[ -z "${var_name}" ]]; then
			continue
		fi
		summary_vars+=( "${var_name}" )
	done < <(gg_print_entrypoint_config_vars "${entrypoint_name}")

	while [[ $# -gt 0 ]]; do
		var_name=$(gg_normalize_config_var_name "$1" || true)
		if [[ -n "${var_name}" ]]; then
			summary_vars+=( "${var_name}" )
		fi
		shift
	done

	if [[ ${#summary_vars[@]} -eq 0 ]]; then
		return 0
	fi
	gg_print_named_config_summary "${title}" "${summary_vars[@]}"
}

gg_container_mount_destination() {
	local mount_spec=${1:-}
	local remainder
	local destination

	if [[ -z "${mount_spec}" ]]; then
		return 1
	fi
	remainder=${mount_spec#*:}
	if [[ "${remainder}" == "${mount_spec}" ]]; then
		return 1
	fi
	destination=${remainder%%:*}
	if [[ -z "${destination}" ]]; then
		return 1
	fi
	printf '%s\n' "${destination}"
}

gg_container_bind_destination_exists_in_csv() {
	local destination=${1:-}
	local csv_mounts=${2:-}
	local mount_entry
	local entry_destination

	if [[ -z "${destination}" || -z "${csv_mounts}" ]]; then
		return 1
	fi

	IFS=',' read -r -a mount_entries <<< "${csv_mounts}"
	for mount_entry in "${mount_entries[@]}"; do
		[[ -n "${mount_entry}" ]] || continue
		entry_destination=$(gg_container_mount_destination "${mount_entry}" || true)
		if [[ -n "${entry_destination}" && "${entry_destination}" == "${destination}" ]]; then
			return 0
		fi
	done
	return 1
}

gg_container_bind_csv_normalize() {
	local csv_mounts
	local mount_entry
	local entry_destination
	local normalized=""

	for csv_mounts in "$@"; do
		[[ -n "${csv_mounts}" ]] || continue
		IFS=',' read -r -a mount_entries <<< "${csv_mounts}"
		for mount_entry in "${mount_entries[@]}"; do
			[[ -n "${mount_entry}" ]] || continue
			entry_destination=$(gg_container_mount_destination "${mount_entry}" || true)
			if [[ -z "${entry_destination}" ]]; then
				continue
			fi
			if gg_container_bind_destination_exists_in_csv "${entry_destination}" "${normalized}"; then
				continue
			fi
			normalized=$(gg_csv_append "${normalized}" "${mount_entry}")
		done
	done

	printf '%s\n' "${normalized}"
}

gg_container_shell_command_is_set() {
	if ! declare -p singularity_command >/dev/null 2>&1; then
		return 1
	fi
	case "$(declare -p singularity_command 2>/dev/null)" in
		declare\ -a*)
			[[ ${#singularity_command[@]} -gt 0 ]]
			;;
		*)
			[[ -n "${singularity_command:-}" ]]
			;;
	esac
}

gg_container_shell_command_runtime_bin() {
	if ! gg_container_shell_command_is_set; then
		return 1
	fi
	case "$(declare -p singularity_command 2>/dev/null)" in
		declare\ -a*)
			printf '%s\n' "${singularity_command[0]}"
			;;
		*)
			printf '%s\n' "${singularity_command%% *}"
			;;
	esac
}

gg_container_shell_command_subcommand() {
	if ! gg_container_shell_command_is_set; then
		return 1
	fi
	case "$(declare -p singularity_command 2>/dev/null)" in
		declare\ -a*)
			if [[ ${#singularity_command[@]} -lt 2 ]]; then
				return 1
			fi
			printf '%s\n' "${singularity_command[1]}"
			;;
		*)
			local runtime_bin=""
			local subcommand=""
			read -r runtime_bin subcommand _ <<< "${singularity_command}"
			if [[ -z "${subcommand}" ]]; then
				return 1
			fi
			printf '%s\n' "${subcommand}"
			;;
	esac
}

gg_container_shell_command_display() {
	if ! gg_container_shell_command_is_set; then
		return 1
	fi
	case "$(declare -p singularity_command 2>/dev/null)" in
		declare\ -a*)
			local display=""
			local arg
			local quoted_arg=""
			for arg in "${singularity_command[@]}"; do
				printf -v quoted_arg '%q' "${arg}"
				if [[ -n "${display}" ]]; then
					display="${display} "
				fi
				display="${display}${quoted_arg}"
			done
			printf '%s\n' "${display}"
			;;
		*)
			printf '%s\n' "${singularity_command}"
			;;
	esac
}

gg_run_container_shell_script() {
	local image_path=$1
	local script_path=$2
	local subcommand=""

	if ! gg_container_shell_command_is_set; then
		echo "gg_run_container_shell_script: container shell command is not initialized." >&2
		return 1
	fi
	if [[ ! -s "${script_path}" ]]; then
		echo "gg_run_container_shell_script: script not found: ${script_path}" >&2
		return 1
	fi
	subcommand=$(gg_container_shell_command_subcommand || true)
	case "$(declare -p singularity_command 2>/dev/null)" in
		declare\ -a*)
			if [[ "${subcommand}" == "exec" ]]; then
				"${singularity_command[@]}" "${image_path}" bash -s -- < "${script_path}"
			else
				"${singularity_command[@]}" "${image_path}" < "${script_path}"
			fi
			;;
		*)
			# Backward compatibility for external site adapters that still set a string command.
			if [[ "${subcommand}" == "exec" ]]; then
				${singularity_command} "${image_path}" bash -s -- < "${script_path}"
			else
				${singularity_command} "${image_path}" < "${script_path}"
			fi
			;;
	esac
}

gg_sanitize_runtime_path_component() {
	local value=${1:-}
	value="${value//[^[:alnum:]._-]/_}"
	if [[ -z "${value}" ]]; then
		value="unknown"
	fi
	printf '%s\n' "${value}"
}

gg_entrypoint_runtime_snapshot_dir() {
	local entrypoint_name=${1:-}
	local entrypoint_stem=""
	local job_id=""
	local task_id=""
	local output_root=""
	local scheduler_kind=""

	entrypoint_stem="$(basename "${entrypoint_name}")"
	entrypoint_stem="${entrypoint_stem%.sh}"
	job_id="${GG_JOB_ID:-${SLURM_JOB_ID:-${PBS_JOBID:-${JOB_ID:-local_$$}}}}"
	task_id="${GG_ARRAY_TASK_ID:-${SLURM_ARRAY_TASK_ID:-${PBS_ARRAY_INDEX:-${PBS_ARRAYID:-${SGE_TASK_ID:-1}}}}}"
	scheduler_kind="${GG_SCHEDULER_KIND:-$(gg_detect_scheduler_kind)}"
	if [[ "${scheduler_kind}" == "local" ]]; then
		job_id="${job_id}.${BASHPID:-$$}"
	fi

	output_root="${gg_workspace_output_dir:-}"
	if [[ -z "${output_root}" && -n "${gg_workspace_dir:-}" ]]; then
		output_root="$(workspace_output_root "${gg_workspace_dir}")"
	fi
	if [[ -z "${output_root}" ]]; then
		output_root="${PWD}/workspace/output"
	fi

	printf '%s\n' \
		"${output_root}/runtime/$(gg_sanitize_runtime_path_component "${entrypoint_stem}")/$(gg_sanitize_runtime_path_component "${job_id}")_$(gg_sanitize_runtime_path_component "${task_id}")"
}

gg_prepare_entrypoint_runtime_snapshot() {
	local entrypoint_name=${1:-}
	local script_path=${2:-}
	local runtime_dir=""
	local snapshot_script=""
	local entrypoint_path=""
	local metadata_path=""

	if [[ -z "${entrypoint_name}" ]]; then
		echo "gg_prepare_entrypoint_runtime_snapshot: entrypoint name is required." >&2
		return 1
	fi
	if [[ ! -s "${script_path}" ]]; then
		echo "gg_prepare_entrypoint_runtime_snapshot: script not found: ${script_path}" >&2
		return 1
	fi

	runtime_dir="$(gg_entrypoint_runtime_snapshot_dir "${entrypoint_name}")"
	mkdir -p -- "${runtime_dir}"

	snapshot_script="${runtime_dir}/$(basename "${script_path}")"
	cp -- "${script_path}" "${snapshot_script}"
	chmod --reference="${script_path}" "${snapshot_script}" 2>/dev/null || chmod a+r "${snapshot_script}" 2>/dev/null || true

	if [[ -n "${gg_workflow_dir:-}" ]]; then
		entrypoint_path="${gg_workflow_dir}/$(basename "${entrypoint_name}")"
		if [[ -s "${entrypoint_path}" ]]; then
			cp -- "${entrypoint_path}" "${runtime_dir}/$(basename "${entrypoint_path}")"
		fi
		for metadata_path in "${gg_workflow_dir}/gg_common_params.sh" "${gg_workflow_dir}/../VERSION"; do
			if [[ -s "${metadata_path}" ]]; then
				cp -- "${metadata_path}" "${runtime_dir}/$(basename "${metadata_path}")"
			fi
		done
	fi

	echo "Runtime snapshot directory: ${runtime_dir}" >&2
	echo "Runtime snapshot core script: ${snapshot_script}" >&2
	printf '%s\n' "${snapshot_script}"
}

gg_detect_active_container_runtime() {
	local runtime_bin=""

	if runtime_bin=$(gg_container_shell_command_runtime_bin 2>/dev/null); then
		:
	fi
	if [[ -z "${runtime_bin}" ]]; then
		runtime_bin=$(gg_detect_container_runtime_binary || true)
	fi
	printf '%s\n' "${runtime_bin}"
}

# Optional runtime override is used by callers outside the ShellCheck driver set.
# shellcheck disable=SC2120
gg_sync_container_bind_envs() {
	local runtime_bin=${1:-}
	local bind_mounts=""

	if [[ -z "${runtime_bin}" ]]; then
		runtime_bin=$(gg_detect_active_container_runtime)
	fi
	bind_mounts=$(gg_container_bind_csv_normalize \
		"${GG_CONTAINER_BIND_MOUNTS:-}" \
		"${SINGULARITY_BIND:-}" \
		"${SINGULARITY_BINDPATH:-}" \
		"${APPTAINER_BIND:-}" \
		"${APPTAINER_BINDPATH:-}")

	unset SINGULARITY_BIND
	unset SINGULARITY_BINDPATH
	unset APPTAINER_BIND
	unset APPTAINER_BINDPATH

	GG_CONTAINER_BIND_MOUNTS="${bind_mounts}"
	export GG_CONTAINER_BIND_MOUNTS

	if [[ -z "${bind_mounts}" ]]; then
		return 0
	fi

	case "${runtime_bin}" in
		apptainer)
			APPTAINER_BINDPATH="${bind_mounts}"
			export APPTAINER_BINDPATH
			;;
		*)
			SINGULARITY_BINDPATH="${bind_mounts}"
			export SINGULARITY_BINDPATH
			;;
	esac
}

gg_container_bind_destination_exists() {
	local mount_spec=${1:-}
	local destination

	destination=$(gg_container_mount_destination "${mount_spec}" || true)
	if [[ -z "${destination}" ]]; then
		return 1
	fi
	if gg_container_bind_destination_exists_in_csv "${destination}" "${GG_CONTAINER_BIND_MOUNTS:-}"; then
		return 0
	fi
	if gg_container_bind_destination_exists_in_csv "${destination}" "${SINGULARITY_BIND:-}"; then
		return 0
	fi
	if gg_container_bind_destination_exists_in_csv "${destination}" "${SINGULARITY_BINDPATH:-}"; then
		return 0
	fi
	if gg_container_bind_destination_exists_in_csv "${destination}" "${APPTAINER_BIND:-}"; then
		return 0
	fi
	if gg_container_bind_destination_exists_in_csv "${destination}" "${APPTAINER_BINDPATH:-}"; then
		return 0
	fi
	return 1
}

gg_add_container_bind_mount() {
	local mount_spec=$1
	if gg_container_bind_destination_exists "${mount_spec}"; then
		return 0
	fi
	GG_CONTAINER_BIND_MOUNTS=$(gg_csv_prepend "${mount_spec}" "${GG_CONTAINER_BIND_MOUNTS:-}")
	export GG_CONTAINER_BIND_MOUNTS
	gg_sync_container_bind_envs
}

gg_export_var_to_container_env_if_set() {
	local var_name=$1
	var_name=$(gg_normalize_config_var_name "${var_name}" || true)
	if [[ -z "${var_name}" ]]; then
		return 0
	fi
	if ! gg_is_valid_shell_var_name "${var_name}"; then
		echo "gg_export_var_to_container_env_if_set: ignoring invalid variable name: ${var_name}" >&2
		return 0
	fi
	if [[ -n "${!var_name+x}" ]]; then
		export "${var_name?}"
		export "SINGULARITYENV_${var_name}=${!var_name}"
		export "APPTAINERENV_${var_name}=${!var_name}"
	fi
}

gg_forward_env_vars_with_prefix_to_container_env() {
	local prefix=${1:-}
	local var_name

	if [[ -z "${prefix}" ]]; then
		return 0
	fi

	while IFS= read -r var_name; do
		[[ -n "${var_name}" ]] || continue
		gg_export_var_to_container_env_if_set "${var_name}"
	done < <(compgen -A variable "${prefix}" | LC_ALL=C sort || true)
}

gg_requested_container_runtime() {
	local requested_runtime="${GG_CONTAINER_RUNTIME:-auto}"

	requested_runtime=$(printf '%s' "${requested_runtime}" | tr '[:upper:]' '[:lower:]')
	case "${requested_runtime}" in
		""|auto)
			printf '%s\n' "auto"
			;;
		docker)
			printf '%s\n' "docker"
			;;
		*)
			echo "Unsupported GG_CONTAINER_RUNTIME=${GG_CONTAINER_RUNTIME:-}. Use auto or docker." >&2
			return 1
			;;
	esac
}
