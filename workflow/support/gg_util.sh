#!/usr/bin/env bash

gg_util_dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd 2>/dev/null || true)"
if [[ -n "${gg_util_dir}" && -s "${gg_util_dir}/gg_site_runtime.sh" ]]; then
	# shellcheck disable=SC1090
	source "${gg_util_dir}/gg_site_runtime.sh"
fi
if [[ -n "${gg_util_dir}" && -s "${gg_util_dir}/gg_entrypoint_config_vars.sh" ]]; then
	# shellcheck disable=SC1090
	source "${gg_util_dir}/gg_entrypoint_config_vars.sh"
fi
unset gg_util_dir

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
	resolved_path=$(cd "$(dirname "${path}")" && pwd -P)/$(basename "${path}")
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

gg_detect_active_container_runtime() {
	local runtime_bin=""

	if [[ -n "${singularity_command:-}" ]]; then
		runtime_bin="${singularity_command%% *}"
	fi
	if [[ -z "${runtime_bin}" ]]; then
		runtime_bin=$(gg_detect_container_runtime_binary || true)
	fi
	printf '%s\n' "${runtime_bin}"
}

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
	if [[ -z "${var_name}" ]]; then
		return 0
	fi
	if [[ -n "${!var_name+x}" ]]; then
		export "${var_name}"
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

gg_apply_named_env_overrides() {
	local var_name
	local env_name
	local env_value

	if [[ $(($# % 2)) -ne 0 ]]; then
		echo "gg_apply_named_env_overrides: expected var/env pairs." >&2
		return 1
	fi

	while [[ $# -gt 0 ]]; do
		var_name=$1
		env_name=$2
		shift 2
		env_value="${!env_name:-}"
		if [[ -n "${env_value}" ]]; then
			printf -v "${var_name}" '%s' "${env_value}"
		fi
	done
}

forward_config_vars_to_container_env() {
	local job_script=$1
	shift || true
	local entrypoint_name
	local var_name
	entrypoint_name="$(basename "${job_script}")"
	if ! declare -F gg_print_entrypoint_config_vars >/dev/null 2>&1; then
		echo "forward_config_vars_to_container_env: config var registry helper is unavailable." >&2
		return 1
	fi
	while IFS= read -r var_name; do
		if [[ -z "${var_name}" ]]; then
			continue
		fi
		gg_export_var_to_container_env_if_set "${var_name}"
	done < <(gg_print_entrypoint_config_vars "${entrypoint_name}")

	# gg_debug_mode is read by gg_*_core.sh scripts but defined outside the config block.
	gg_export_var_to_container_env_if_set "gg_debug_mode"

	while [[ $# -gt 0 ]]; do
		gg_export_var_to_container_env_if_set "$1"
		shift
	done
}

gg_detect_container_runtime_binary() {
	if command -v singularity >/dev/null 2>&1; then
		echo singularity
		return 0
	fi
	if command -v apptainer >/dev/null 2>&1; then
		echo apptainer
		return 0
	fi
	return 1
}

gg_print_container_env_summary() {
	local echo_header="set_singularityenv: "
	local forwarded_env_count=0
	forwarded_env_count=$(env | awk -F= '/^(SINGULARITYENV_|APPTAINERENV_)/{n++} END{print n+0}')
	echo "${echo_header}GG_CONTAINER_BIND_MOUNTS=${GG_CONTAINER_BIND_MOUNTS:-}"
	echo "${echo_header}SINGULARITY_BIND=${SINGULARITY_BIND:-}"
	echo "${echo_header}SINGULARITY_BINDPATH=${SINGULARITY_BINDPATH:-}"
	echo "${echo_header}APPTAINER_BIND=${APPTAINER_BIND:-}"
	echo "${echo_header}APPTAINER_BINDPATH=${APPTAINER_BINDPATH:-}"
	echo "${echo_header}forwarded_container_env_vars=${forwarded_env_count}"
	echo ""
}

gg_detect_scheduler_kind() {
	local normalized_job_id="${GG_JOB_ID:-${JOB_ID:-}}"
	if [[ -n "${SLURM_JOB_ID:-}" ]]; then
		echo slurm
		return 0
	fi
	if [[ -n "${PBS_JOBID:-}" || -n "${PBS_O_WORKDIR:-}" || -n "${PBS_NODEFILE:-}" ]]; then
		echo pbs
		return 0
	fi
	if [[ -n "${SGE_ROOT:-}" || -n "${PE_HOSTFILE:-}" || -n "${QUEUE:-}" || -n "${ARC:-}" ]]; then
		echo uge
		return 0
	fi
	if [[ -n "${normalized_job_id}" && "${normalized_job_id}" != "1" ]]; then
		echo uge
		return 0
	fi
	echo local
	return 0
}

gg_sync_legacy_scheduler_aliases() {
	NSLOTS=${GG_TASK_CPUS:-1}
	JOB_ID=${GG_JOB_ID:-1}
	SGE_TASK_ID=${GG_ARRAY_TASK_ID:-1}
	MEM_PER_SLOT=${GG_MEM_PER_CPU_GB:-3}
	MEM_PER_HOST=${GG_MEM_TOTAL_GB:-$((MEM_PER_SLOT * NSLOTS))}
	export GG_TASK_CPUS GG_JOB_ID GG_ARRAY_TASK_ID GG_MEM_PER_CPU_GB GG_MEM_TOTAL_GB
	export NSLOTS JOB_ID SGE_TASK_ID MEM_PER_SLOT MEM_PER_HOST
}

gg_normalize_scheduler_env() {
	local echo_header='gg_normalize_scheduler_env: '
	local scheduler_kind
	local pbs_slots=""
	scheduler_kind=$(gg_detect_scheduler_kind)
	GG_SCHEDULER_KIND=${scheduler_kind}
	echo ${echo_header}'Scheduler metadata is normalized to GG_* variables.'
	if [[ -z "${GG_TASK_CPUS:-}" ]]; then
		if [[ -n "${NSLOTS:-}" ]]; then
			echo ${echo_header}'GG_TASK_CPUS=${NSLOTS} (from legacy NSLOTS)'
			GG_TASK_CPUS=${NSLOTS}
		else
			echo ${echo_header}'${GG_TASK_CPUS} is undefined or empty.'
		fi
	fi
	if [[ -z "${GG_TASK_CPUS:-}" ]]; then
		if [[ -n "${SLURM_CPUS_PER_TASK:-}" ]]; then
			echo ${echo_header}'GG_TASK_CPUS=${SLURM_CPUS_PER_TASK}'
			GG_TASK_CPUS=${SLURM_CPUS_PER_TASK}
		elif [[ -n "${PBS_NODEFILE:-}" && -f "${PBS_NODEFILE}" ]]; then
			pbs_slots=$(wc -l < "${PBS_NODEFILE}")
			echo ${echo_header}'GG_TASK_CPUS=${pbs_slots}'
			GG_TASK_CPUS=${pbs_slots}
		else
			echo ${echo_header}'No scheduler CPU metadata was detected. GG_TASK_CPUS=1'
			GG_TASK_CPUS=1
		fi
	fi
	if [[ -z "${GG_JOB_ID:-}" ]]; then
		if [[ -n "${JOB_ID:-}" ]]; then
			echo ${echo_header}'GG_JOB_ID=${JOB_ID} (from legacy JOB_ID)'
			GG_JOB_ID=${JOB_ID}
		else
			echo ${echo_header}'${GG_JOB_ID} is undefined.'
		fi
	fi
	if [[ -z "${GG_JOB_ID:-}" ]]; then
		if [[ -n "${SLURM_JOB_ID:-}" ]]; then
			echo ${echo_header}'GG_JOB_ID=${SLURM_JOB_ID}'
			GG_JOB_ID=${SLURM_JOB_ID}
		elif [[ -n "${PBS_JOBID:-}" ]]; then
			echo ${echo_header}'GG_JOB_ID=${PBS_JOBID}'
			GG_JOB_ID=${PBS_JOBID}
		else
			echo ${echo_header}'No scheduler job ID was detected. GG_JOB_ID=1'
			GG_JOB_ID=1
		fi
	fi
	if [[ -z "${GG_ARRAY_TASK_ID:-}" ]]; then
		if [[ -n "${SGE_TASK_ID:-}" ]]; then
			echo ${echo_header}'GG_ARRAY_TASK_ID=${SGE_TASK_ID} (from legacy SGE_TASK_ID)'
			GG_ARRAY_TASK_ID=${SGE_TASK_ID}
		else
			echo ${echo_header}'${GG_ARRAY_TASK_ID} is undefined.'
		fi
	fi
	if [[ -z "${GG_ARRAY_TASK_ID:-}" ]]; then
		if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
			echo ${echo_header}'GG_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}'
			GG_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}
		elif [[ -n "${PBS_ARRAY_INDEX:-}" ]]; then
			echo ${echo_header}'GG_ARRAY_TASK_ID=${PBS_ARRAY_INDEX}'
			GG_ARRAY_TASK_ID=${PBS_ARRAY_INDEX}
		elif [[ -n "${PBS_ARRAYID:-}" ]]; then
			echo ${echo_header}'GG_ARRAY_TASK_ID=${PBS_ARRAYID}'
			GG_ARRAY_TASK_ID=${PBS_ARRAYID}
		else
			echo ${echo_header}'No scheduler array ID was detected. GG_ARRAY_TASK_ID=1'
			GG_ARRAY_TASK_ID=1
		fi
	fi
	if [[ -z "${GG_MEM_PER_CPU_GB:-}" && -n "${MEM_PER_SLOT:-}" ]]; then
		echo ${echo_header}'GG_MEM_PER_CPU_GB='"${MEM_PER_SLOT}"' (from legacy MEM_PER_SLOT)'
		GG_MEM_PER_CPU_GB=${MEM_PER_SLOT}
	fi
	if [[ -z "${GG_MEM_PER_CPU_GB:-}" ]] && type qstat >/dev/null 2>&1; then
		GG_MEM_PER_CPU_GB=$(
			{ qstat -f -j "${GG_JOB_ID}" 2>/dev/null || true; } \
			| awk -F'mem_req=|G' '/mem_req=[0-9]+G/ { if (!seen) { print $2; seen=1 } }'
		)
	fi
	if [[ -z "${GG_MEM_PER_CPU_GB:-}" ]]; then
		echo ${echo_header}'${GG_MEM_PER_CPU_GB} is undefined.'
		if [[ -n "${SLURM_MEM_PER_CPU:-}" ]]; then
			echo ${echo_header}'GG_MEM_PER_CPU_GB=$((${SLURM_MEM_PER_CPU}/1024))'
			GG_MEM_PER_CPU_GB=$((${SLURM_MEM_PER_CPU}/1024))
		else
			echo ${echo_header}'No scheduler memory-per-CPU metadata was detected. GG_MEM_PER_CPU_GB=3'
			GG_MEM_PER_CPU_GB=3
		fi
	fi
	if [[ -z "${GG_MEM_TOTAL_GB:-}" && -n "${MEM_PER_HOST:-}" ]]; then
		echo ${echo_header}'GG_MEM_TOTAL_GB='"${MEM_PER_HOST}"' (from legacy MEM_PER_HOST)'
		GG_MEM_TOTAL_GB=${MEM_PER_HOST}
	fi
	if [[ -z "${GG_MEM_TOTAL_GB:-}" ]]; then
		GG_MEM_TOTAL_GB=$((${GG_MEM_PER_CPU_GB}*${GG_TASK_CPUS}))
	fi
	gg_sync_legacy_scheduler_aliases
	echo ${echo_header}"GG_TASK_CPUS=${GG_TASK_CPUS}"
	echo ${echo_header}"GG_JOB_ID=${GG_JOB_ID}"
	echo ${echo_header}"GG_ARRAY_TASK_ID=${GG_ARRAY_TASK_ID}"
	echo ${echo_header}"GG_MEM_PER_CPU_GB=${GG_MEM_PER_CPU_GB}"
	echo ${echo_header}"GG_MEM_TOTAL_GB=${GG_MEM_TOTAL_GB}"
	echo ""
	export GG_SCHEDULER_KIND
}

variable_SGEnizer() {
	gg_normalize_scheduler_env "$@"
}

gg_print_scheduler_runtime_summary() {
	local echo_header='scheduler_runtime_summary: '
	local scheduler='local'
	local pbs_slots='NA'

	if [[ -n "${GG_SCHEDULER_KIND:-}" ]]; then
		scheduler="${GG_SCHEDULER_KIND}"
	else
		scheduler=$(gg_detect_scheduler_kind)
	fi
	if [[ -n "${PBS_NODEFILE:-}" && -f "${PBS_NODEFILE}" ]]; then
		pbs_slots=$(wc -l < "${PBS_NODEFILE}")
	fi

	echo "${echo_header}scheduler=${scheduler}"
	echo "${echo_header}requested.slurm cpus_per_task=${SLURM_CPUS_PER_TASK:-NA} mem_per_cpu_mb=${SLURM_MEM_PER_CPU:-NA} array_task_id=${SLURM_ARRAY_TASK_ID:-NA} job_id=${SLURM_JOB_ID:-NA}"
	echo "${echo_header}requested.pbs nodefile_slots=${pbs_slots} array_index=${PBS_ARRAY_INDEX:-NA} array_id=${PBS_ARRAYID:-NA} job_id=${PBS_JOBID:-NA}"
	echo "${echo_header}legacy_aliases NSLOTS=${NSLOTS:-NA} SGE_TASK_ID=${SGE_TASK_ID:-NA} JOB_ID=${JOB_ID:-NA} MEM_PER_SLOT=${MEM_PER_SLOT:-NA} MEM_PER_HOST=${MEM_PER_HOST:-NA}"
	echo "${echo_header}detected GG_TASK_CPUS=${GG_TASK_CPUS:-NA} GG_MEM_PER_CPU_GB=${GG_MEM_PER_CPU_GB:-NA} GG_MEM_TOTAL_GB=${GG_MEM_TOTAL_GB:-NA} GG_JOB_ID=${GG_JOB_ID:-NA} GG_ARRAY_TASK_ID=${GG_ARRAY_TASK_ID:-NA}"
	echo "${echo_header}forwarded SINGULARITYENV_GG_TASK_CPUS=${SINGULARITYENV_GG_TASK_CPUS:-unset} SINGULARITYENV_GG_MEM_PER_CPU_GB=${SINGULARITYENV_GG_MEM_PER_CPU_GB:-unset} SINGULARITYENV_GG_ARRAY_TASK_ID=${SINGULARITYENV_GG_ARRAY_TASK_ID:-unset}"
	echo "${echo_header}forwarded APPTAINERENV_GG_TASK_CPUS=${APPTAINERENV_GG_TASK_CPUS:-unset} APPTAINERENV_GG_MEM_PER_CPU_GB=${APPTAINERENV_GG_MEM_PER_CPU_GB:-unset} APPTAINERENV_GG_ARRAY_TASK_ID=${APPTAINERENV_GG_ARRAY_TASK_ID:-unset}"
	echo ""
}

set_singularityenv() {
  local echo_header="set_singularityenv: "
  local resolved_workspace_dir="${gg_workspace_dir}"
  local resolved_workflow_dir="${gg_workflow_dir}"
  local resolved_container_image_path="${gg_container_image_path}"
  local resolved_workspace_layout=""
  echo ${echo_header}"original: gg_workspace_dir = ${gg_workspace_dir}"
  echo ${echo_header}"original: gg_workflow_dir = ${gg_workflow_dir}"
  echo ${echo_header}"original: gg_container_image_path = ${gg_container_image_path}"
  if [[ $(uname -s) != 'Darwin' ]]; then
    echo "OS is $(uname -s). Getting the original path of symlink."
    resolved_workspace_dir=$(gg_resolve_physical_path "${gg_workspace_dir}")
    resolved_workflow_dir=$(gg_resolve_physical_path "${gg_workflow_dir}")
    resolved_container_image_path=$(gg_resolve_physical_path "${gg_container_image_path}")
    echo ${echo_header}"formatted: gg_workspace_dir = ${resolved_workspace_dir}"
    echo ${echo_header}"formatted: gg_workflow_dir = ${resolved_workflow_dir}"
    echo ${echo_header}"formatted: gg_container_image_path = ${resolved_container_image_path}"
	else
		echo "OS is $(uname -s). Symlink PATHs won't be updated."
	fi
	gg_workspace_dir="${resolved_workspace_dir}"
	gg_workflow_dir="${resolved_workflow_dir}"
	gg_container_image_path="${resolved_container_image_path}"
  export gg_workspace_dir
  export gg_workflow_dir
  export gg_container_image_path
	resolved_workspace_layout=$(gg_resolve_workspace_layout "${gg_workspace_dir}")
	gg_add_container_bind_mount "${resolved_workspace_dir}:/workspace"
	gg_add_container_bind_mount "${resolved_workflow_dir}:/script"
	export SINGULARITYENV_GG_ARRAY_TASK_ID=${GG_ARRAY_TASK_ID:-1}
	export APPTAINERENV_GG_ARRAY_TASK_ID=${GG_ARRAY_TASK_ID:-1}
	export SINGULARITYENV_GG_TASK_CPUS=${GG_TASK_CPUS:-1}
	export APPTAINERENV_GG_TASK_CPUS=${GG_TASK_CPUS:-1}
	export SINGULARITYENV_GG_JOB_ID=${GG_JOB_ID:-1}
	export APPTAINERENV_GG_JOB_ID=${GG_JOB_ID:-1}
	export SINGULARITYENV_GG_MEM_PER_CPU_GB=${GG_MEM_PER_CPU_GB:-3}
	export APPTAINERENV_GG_MEM_PER_CPU_GB=${GG_MEM_PER_CPU_GB:-3}
	export SINGULARITYENV_GG_MEM_TOTAL_GB=${GG_MEM_TOTAL_GB:-3}
	export APPTAINERENV_GG_MEM_TOTAL_GB=${GG_MEM_TOTAL_GB:-3}
	export SINGULARITYENV_SGE_TASK_ID=${SGE_TASK_ID:-1}
	export APPTAINERENV_SGE_TASK_ID=${SGE_TASK_ID:-1}
	export SINGULARITYENV_NSLOTS=${NSLOTS:-1}
	export APPTAINERENV_NSLOTS=${NSLOTS:-1}
	export SINGULARITYENV_JOB_ID=${JOB_ID:-1}
	export APPTAINERENV_JOB_ID=${JOB_ID:-1}
	export SINGULARITYENV_MEM_PER_SLOT=${MEM_PER_SLOT:-3}
	export APPTAINERENV_MEM_PER_SLOT=${MEM_PER_SLOT:-3}
	export SINGULARITYENV_MEM_PER_HOST=${MEM_PER_HOST:-3}
	export APPTAINERENV_MEM_PER_HOST=${MEM_PER_HOST:-3}
  local gg_common_var_name
  for gg_common_var_name in ${!GG_COMMON_@}; do
    if [[ -n "${!gg_common_var_name+x}" ]]; then
      export "SINGULARITYENV_${gg_common_var_name}=${!gg_common_var_name}"
      export "APPTAINERENV_${gg_common_var_name}=${!gg_common_var_name}"
    fi
  done
  if [[ -n "${delete_preexisting_tmp_dir:-}" ]]; then
    export SINGULARITYENV_delete_preexisting_tmp_dir=${delete_preexisting_tmp_dir}
    export APPTAINERENV_delete_preexisting_tmp_dir=${delete_preexisting_tmp_dir}
  else
    export SINGULARITYENV_delete_preexisting_tmp_dir=0
    export APPTAINERENV_delete_preexisting_tmp_dir=0
  fi
	if [[ -n "${delete_tmp_dir:-}" ]]; then
		export SINGULARITYENV_delete_tmp_dir=${delete_tmp_dir}
		export APPTAINERENV_delete_tmp_dir=${delete_tmp_dir}
	else
		export SINGULARITYENV_delete_tmp_dir=1
		export APPTAINERENV_delete_tmp_dir=1
	fi
	gg_print_container_env_summary
}

ensure_dir() {
	local d=$1
	if [[ -z "${d}" ]]; then
		return 0
	fi
	if [[ -d "${d}" ]]; then
		return 0
	fi

	local attempt
	local mkdir_err
	mkdir_err=$(mktemp)
	for attempt in 1 2 3 4 5; do
		if mkdir -p "${d}" 2>"${mkdir_err}"; then
			rm -f -- "${mkdir_err}"
			return 0
		fi
		if grep -qi "Resource deadlock avoided" "${mkdir_err}"; then
			sleep 1
			continue
		fi
		cat "${mkdir_err}" >&2
		rm -f -- "${mkdir_err}"
		return 1
	done

	cat "${mkdir_err}" >&2
	rm -f -- "${mkdir_err}"
	return 1
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

gg_annotation_species_priority() {
  cat <<'EOF'
Arabidopsis_thaliana
Oryza_sativa
Homo_sapiens
Mus_musculus
Danio_rerio
Drosophila_melanogaster
Caenorhabditis_elegans
Saccharomyces_cerevisiae
Schizosaccharomyces_pombe
Escherichia_coli
EOF
}

gg_normalize_annotation_species() {
  local species=${1:-}
  species=$(printf '%s' "${species}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' -e 's/[[:space:]]\+/_/g' -e 's/_$//')
  printf '%s\n' "${species}"
}

gg_normalize_busco_lineage_request() {
  local requested=${1:-}
  requested=$(printf '%s' "${requested}" | sed -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//')
  printf '%s\n' "${requested}"
}

gg_finalize_auto_busco_lineage_name() {
  local lineage=${1:-}
  local odb_version=${2:-12}
  if [[ -z "${lineage}" ]]; then
    return 1
  fi
  if [[ "${lineage}" == *_odb* ]]; then
    printf '%s\n' "${lineage}"
    return 0
  fi
  printf '%s_odb%s\n' "${lineage}" "${odb_version}"
}

gg_species_names_from_fasta_dir() {
  local search_dir=${1:-}
  local file
  local species_name=""

  while IFS= read -r file; do
    species_name=$(gg_species_name_from_path_or_dot "${file}")
    [[ -n "${species_name}" ]] || continue
    printf '%s\n' "${species_name}"
  done < <(gg_find_fasta_files "${search_dir}" 1) | sort -u
}

gg_species_names_from_annotation_dir() {
  local search_dir=${1:-}
  local file_base
  local species_name=""

  while IFS= read -r file_base; do
    species_name=$(gg_species_name_from_path_or_dot "${file_base}")
    [[ -n "${species_name}" ]] || continue
    printf '%s\n' "${species_name}"
  done < <(gg_find_file_basenames "${search_dir}" "*_annotation.tsv" 1) | sort -u
}

gg_resolve_annotation_species() {
  local requested=${1:-auto}
  shift || true
  local normalized_requested=""
  local requested_lc=""
  local candidates=()
  local candidate=""
  local available=()
  local preferred=""

  normalized_requested=$(gg_normalize_annotation_species "${requested}")
  requested_lc=$(printf '%s' "${normalized_requested}" | tr '[:upper:]' '[:lower:]')
  if [[ -n "${normalized_requested}" && "${requested_lc}" != "auto" ]]; then
    printf '%s\n' "${normalized_requested}"
    return 0
  fi

  for candidate in "$@"; do
    candidate=$(gg_normalize_annotation_species "${candidate}")
    [[ -n "${candidate}" ]] || continue
    candidates+=( "${candidate}" )
  done
  if [[ ${#candidates[@]} -eq 0 ]]; then
    return 1
  fi

  while IFS= read -r candidate; do
    [[ -n "${candidate}" ]] || continue
    available+=( "${candidate}" )
  done < <(printf '%s\n' "${candidates[@]}" | sed -e '/^[[:space:]]*$/d' | sort -u)
  while IFS= read -r preferred; do
    [[ -n "${preferred}" ]] || continue
    for candidate in "${available[@]}"; do
      if [[ "${candidate}" == "${preferred}" ]]; then
        printf '%s\n' "${candidate}"
        return 0
      fi
    done
  done < <(gg_annotation_species_priority)

  printf '%s\n' "${available[0]}"
  return 0
}

gg_find_python_exec() {
  local candidate

  for candidate in python python3 /opt/conda/bin/python /usr/bin/python3; do
    if [[ -x "${candidate}" ]]; then
      printf '%s\n' "${candidate}"
      return 0
    fi
    if command -v "${candidate}" >/dev/null 2>&1; then
      command -v "${candidate}"
      return 0
    fi
  done
  return 1
}

workspace_busco_placement_root() {
  local gg_workspace_dir=$1
  local dir_db
  dir_db=$(workspace_downloads_root "${gg_workspace_dir}")
  echo "${dir_db}/busco_downloads/placement_files"
}

gg_latest_busco_mapping_odb_version_from_dir() {
  local mapping_dir=${1:-}
  local py_exec=""

  if [[ -z "${mapping_dir}" || ! -d "${mapping_dir}" ]]; then
    return 1
  fi
  py_exec=$(gg_find_python_exec || true)
  if [[ -z "${py_exec}" ]]; then
    return 1
  fi

  GG_BUSCO_MAPPING_DIR="${mapping_dir}" \
  "${py_exec}" - <<'PY'
import os
import re
import sys
from pathlib import Path

mapping_dir = Path(os.environ.get("GG_BUSCO_MAPPING_DIR", "").strip())
if not mapping_dir:
    raise SystemExit(1)

pattern = re.compile(r"mapping_taxids-busco_dataset_name\.(archaea|bacteria|eukaryota)_odb(\d+)\..*\.txt$")
versions = set()
for mapping_file in mapping_dir.glob("mapping_taxids-busco_dataset_name.*_odb*.txt"):
    match = pattern.fullmatch(mapping_file.name)
    if not match:
        continue
    _, version = match.groups()
    versions.add(int(version))

if not versions:
    raise SystemExit(1)
print(max(versions))
PY
}

gg_fetch_latest_busco_mapping_odb_version() {
  local py_exec=""

  py_exec=$(gg_find_python_exec || true)
  if [[ -z "${py_exec}" ]]; then
    return 1
  fi

  "${py_exec}" - <<'PY'
import re
import urllib.request

base_url = "https://busco-data.ezlab.org/v5/data/placement_files/"
html = urllib.request.urlopen(base_url, timeout=120).read().decode("utf-8", "replace")
pattern = re.compile(
    r"mapping_taxids-busco_dataset_name\.(archaea|bacteria|eukaryota)_odb(\d+)\.\d{4}-\d{2}-\d{2}\.txt\.tar\.gz"
)
required_domains = {"archaea", "bacteria", "eukaryota"}
versions = {}
for domain, version in pattern.findall(html):
    versions.setdefault(int(version), set()).add(domain)

eligible = [version for version, domains in versions.items() if required_domains.issubset(domains)]
if not eligible:
    raise SystemExit("No common BUSCO ODB placement mapping version was found across archaea, bacteria, and eukaryota.")
print(max(eligible))
PY
}

_download_busco_dataset_mapping_files_locked() {
  local mapping_dir=$1
  local stamp_file=$2
  local odb_version=${3:-}
  local py_exec=""

  py_exec=$(gg_find_python_exec || true)
  if [[ -z "${py_exec}" ]]; then
    echo "python/python3 command was not found. Cannot prepare BUSCO placement mappings." >&2
    return 1
  fi

  if ! GG_BUSCO_MAPPING_DIR="${mapping_dir}" \
    GG_BUSCO_MAPPING_STAMP="${stamp_file}" \
    GG_BUSCO_MAPPING_ODB_VERSION="${odb_version}" \
    "${py_exec}" - <<'PY'
import io
import os
import re
import tarfile
import urllib.request
from pathlib import Path

base_url = "https://busco-data.ezlab.org/v5/data/placement_files/"
mapping_dir = Path(os.environ.get("GG_BUSCO_MAPPING_DIR", "").strip())
stamp_file = Path(os.environ.get("GG_BUSCO_MAPPING_STAMP", "").strip())
odb_version = os.environ.get("GG_BUSCO_MAPPING_ODB_VERSION", "").strip()
if not mapping_dir:
    raise SystemExit("GG_BUSCO_MAPPING_DIR is empty.")
if not stamp_file:
    raise SystemExit("GG_BUSCO_MAPPING_STAMP is empty.")
if not odb_version:
    raise SystemExit("GG_BUSCO_MAPPING_ODB_VERSION is empty.")

mapping_dir.mkdir(parents=True, exist_ok=True)
html = urllib.request.urlopen(base_url, timeout=120).read().decode("utf-8", "replace")
selected = []
for domain in ("archaea", "bacteria", "eukaryota"):
    pattern = re.compile(
        rf'mapping_taxids-busco_dataset_name\.{domain}_odb{re.escape(odb_version)}\.[0-9]{{4}}-[0-9]{{2}}-[0-9]{{2}}\.txt\.tar\.gz'
    )
    matches = sorted(set(pattern.findall(html)))
    if not matches:
        raise SystemExit(f"BUSCO placement mapping for domain/version not found: {domain}, odb{odb_version}")
    archive_name = matches[-1]
    archive_url = base_url + archive_name
    archive_bytes = urllib.request.urlopen(archive_url, timeout=120).read()
    with tarfile.open(fileobj=io.BytesIO(archive_bytes), mode="r:gz") as tar:
        members = [member for member in tar.getmembers() if member.isfile() and member.name.endswith(".txt")]
        if not members:
            raise SystemExit(f"No mapping text file found inside archive: {archive_name}")
        member = members[0]
        extracted = tar.extractfile(member)
        if extracted is None:
            raise SystemExit(f"Failed to extract mapping text from archive: {archive_name}")
        out_path = mapping_dir / Path(member.name).name
        out_path.write_bytes(extracted.read())
        selected.append((domain, archive_name, out_path.name))

stamp_lines = [f"odb_version\t{odb_version}", "domain\tarchive\tmapping_file"]
stamp_lines.extend("\t".join(item) for item in selected)
stamp_file.write_text("\n".join(stamp_lines) + "\n", encoding="utf-8")
PY
  then
    echo "Failed to prepare BUSCO placement mappings in: ${mapping_dir}" >&2
    return 1
  fi

  if [[ ! -s "${stamp_file}" ]]; then
    echo "BUSCO placement mapping stamp file is missing: ${stamp_file}" >&2
    return 1
  fi
}

ensure_busco_dataset_mapping_files() {
  local gg_workspace_dir=$1
  local dir_db
  local mapping_dir
  local lock_file
  local stamp_file
  local remote_odb_version=""
  local odb_version=""

  dir_db="$(workspace_downloads_root "${gg_workspace_dir}")/busco_downloads"
  mapping_dir=$(workspace_busco_placement_root "${gg_workspace_dir}")
  lock_file="${dir_db}/locks/busco_dataset_mappings.lock"

  ensure_dir "${dir_db}"
  ensure_dir "${dir_db}/locks"
  ensure_dir "${mapping_dir}"

  if remote_odb_version=$(gg_fetch_latest_busco_mapping_odb_version 2>/dev/null); then
    odb_version="${remote_odb_version}"
  else
    odb_version=$(gg_latest_busco_mapping_odb_version_from_dir "${mapping_dir}" || true)
  fi
  if [[ -z "${odb_version}" ]]; then
    echo "Failed to determine a BUSCO placement mapping ODB version." >&2
    return 1
  fi
  stamp_file="${mapping_dir}/mapping_taxids-busco_dataset_name.odb${odb_version}.ready.tsv"

  if [[ -n "${remote_odb_version}" ]]; then
    gg_array_download_once "${lock_file}" "${stamp_file}" "BUSCO taxid-to-dataset mapping files (odb${odb_version})" \
      _download_busco_dataset_mapping_files_locked "${mapping_dir}" "${stamp_file}" "${odb_version}" || return 1
  elif [[ ! -s "${stamp_file}" ]]; then
    printf 'odb_version\t%s\n' "${odb_version}" > "${stamp_file}"
  fi

  if [[ ! -s "${stamp_file}" ]]; then
    echo "Failed to prepare BUSCO placement mapping stamp file: ${stamp_file}" >&2
    return 1
  fi

  echo "${mapping_dir}"
}

gg_resolve_busco_lineage_from_lineages() {
  local requested=${1:-auto}
  local mapping_dir=${2:-}
  shift 2 || true
  local normalized_requested=""
  local requested_lc=""
  local py_exec=""

  normalized_requested=$(gg_normalize_busco_lineage_request "${requested}")
  requested_lc=$(printf '%s' "${normalized_requested}" | tr '[:upper:]' '[:lower:]')
  if [[ -n "${normalized_requested}" && "${requested_lc}" != "auto" ]]; then
    printf '%s\n' "${normalized_requested}"
    return 0
  fi
  if [[ -z "${mapping_dir}" ]]; then
    echo "gg_resolve_busco_lineage_from_lineages: mapping_dir is empty." >&2
    return 1
  fi
  if [[ $# -eq 0 ]]; then
    echo "gg_resolve_busco_lineage_from_lineages: no lineage taxid strings were provided." >&2
    return 1
  fi

  py_exec=$(gg_find_python_exec || true)
  if [[ -z "${py_exec}" ]]; then
    echo "python/python3 command was not found. Cannot resolve BUSCO lineage from taxid lineages." >&2
    return 1
  fi

  GG_BUSCO_MAPPING_DIR="${mapping_dir}" \
  "${py_exec}" - "$@" <<'PY'
import os
import re
import sys
from pathlib import Path

mapping_dir = Path(os.environ.get("GG_BUSCO_MAPPING_DIR", "").strip())
if not mapping_dir:
    raise SystemExit("GG_BUSCO_MAPPING_DIR is empty.")

mapping = {}
pattern = re.compile(r"mapping_taxids-busco_dataset_name\.(archaea|bacteria|eukaryota)_odb(\d+)\..*\.txt$")
mapping_files_by_version = {}
for mapping_file in mapping_dir.glob("mapping_taxids-busco_dataset_name.*_odb*.txt"):
    match = pattern.fullmatch(mapping_file.name)
    if not match:
        continue
    _, version = match.groups()
    mapping_files_by_version.setdefault(int(version), []).append(mapping_file)
eligible_versions = sorted(mapping_files_by_version.keys(), reverse=True)
if not eligible_versions:
    raise SystemExit(f"No BUSCO mapping files were found in {mapping_dir}")

lineages = []
for raw_lineage in sys.argv[1:]:
    tokens = [token.strip() for token in str(raw_lineage).split(",") if token.strip()]
    if not tokens:
        raise SystemExit("Encountered an empty lineage taxid string.")
    lineages.append([int(token) for token in tokens])

common_taxids = set(lineages[0])
for lineage in lineages[1:]:
    common_taxids &= set(lineage)

for selected_version in eligible_versions:
    mapping = {}
    for mapping_file in sorted(mapping_files_by_version[selected_version]):
        with mapping_file.open(encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if not line:
                    continue
                parts = line.split("\t", 1)
                if len(parts) != 2:
                    raise SystemExit(f"Malformed BUSCO mapping line in {mapping_file}: {raw_line.rstrip()}")
                taxid_str, dataset = parts
                mapping[int(taxid_str)] = dataset.strip()
    for taxid in reversed(lineages[0]):
        if taxid in common_taxids and taxid in mapping:
            print(mapping[taxid])
            raise SystemExit(0)

raise SystemExit(
    "No BUSCO dataset mapping matched the deepest common taxid across the provided species lineages in any available ODB version."
)
PY
}

gg_busco_taxid_lineages_for_species() {
  local gg_workspace_dir=$1
  shift || true
  local py_exec=""
  local db_file=""
  local species_name=""
  local normalized_species=()

  if [[ $# -eq 0 ]]; then
    echo "gg_busco_taxid_lineages_for_species: no species names were provided." >&2
    return 1
  fi

  if ! ensure_ete_taxonomy_db "${gg_workspace_dir}" >&2; then
    echo "Failed to prepare ETE taxonomy DB for BUSCO lineage resolution." >&2
    return 1
  fi
  db_file=$(workspace_taxonomy_dbfile "${gg_workspace_dir}")
  py_exec=$(gg_find_python_exec || true)
  if [[ -z "${py_exec}" ]]; then
    echo "python/python3 command was not found. Cannot resolve species taxid lineages." >&2
    return 1
  fi

  for species_name in "$@"; do
    species_name=$(gg_normalize_annotation_species "${species_name}")
    [[ -n "${species_name}" ]] || continue
    normalized_species+=( "${species_name}" )
  done
  if [[ ${#normalized_species[@]} -eq 0 ]]; then
    echo "gg_busco_taxid_lineages_for_species: all provided species names were empty after normalization." >&2
    return 1
  fi

  GG_TAXONOMY_DBFILE="${db_file}" \
  "${py_exec}" - "${normalized_species[@]}" <<'PY'
import os
import sys

db_file = os.environ.get("GG_TAXONOMY_DBFILE", "").strip()
if not db_file:
    raise SystemExit("GG_TAXONOMY_DBFILE is empty.")

from ete4 import NCBITaxa  # pragma: no cover - runtime dependency

ncbi = NCBITaxa(dbfile=db_file)

def resolve_species_taxid(species_name: str) -> int:
    candidates = []
    normalized = species_name.replace("_", " ").strip()
    if normalized:
        candidates.append(normalized)
    tokens = normalized.split()
    if len(tokens) >= 2:
        genus_species = " ".join(tokens[:2])
        if genus_species not in candidates:
            candidates.append(genus_species)
    for candidate in candidates:
        translated = ncbi.get_name_translator([candidate])
        if translated:
            return next(iter(translated.values()))[0]
    raise SystemExit(f"Species name was not found in the NCBI taxonomy database: {species_name}")

for species_name in sys.argv[1:]:
    taxid = resolve_species_taxid(species_name)
    lineage = ncbi.get_lineage(taxid)
    if not lineage:
        raise SystemExit(f"NCBI taxonomy lineage was empty for species: {species_name}")
    print(",".join(str(taxid_item) for taxid_item in lineage))
PY
}

gg_resolve_busco_lineage() {
  local gg_workspace_dir=$1
  local requested=${2:-auto}
  shift 2 || true
  local normalized_requested=""
  local requested_lc=""
  local mapping_dir=""
  local lineage_output=""
  local lineage_item=""
  local lineages=()

  normalized_requested=$(gg_normalize_busco_lineage_request "${requested}")
  requested_lc=$(printf '%s' "${normalized_requested}" | tr '[:upper:]' '[:lower:]')
  if [[ -n "${normalized_requested}" && "${requested_lc}" != "auto" ]]; then
    printf '%s\n' "${normalized_requested}"
    return 0
  fi
  if [[ $# -eq 0 ]]; then
    echo "gg_resolve_busco_lineage: no species names were provided for auto resolution." >&2
    return 1
  fi

  if ! mapping_dir=$(ensure_busco_dataset_mapping_files "${gg_workspace_dir}"); then
    echo "Failed to prepare BUSCO placement mapping files for auto lineage resolution." >&2
    return 1
  fi
  if ! lineage_output=$(gg_busco_taxid_lineages_for_species "${gg_workspace_dir}" "$@"); then
    echo "Failed to resolve NCBI taxid lineages for BUSCO auto resolution." >&2
    return 1
  fi
  while IFS= read -r lineage_item; do
    [[ -n "${lineage_item}" ]] || continue
    lineages+=( "${lineage_item}" )
  done <<< "${lineage_output}"
  if [[ ${#lineages[@]} -eq 0 ]]; then
    echo "No lineage taxid strings were produced for BUSCO auto resolution." >&2
    return 1
  fi

  local auto_lineage=""
  if ! auto_lineage=$(gg_resolve_busco_lineage_from_lineages "auto" "${mapping_dir}" "${lineages[@]}"); then
    return 1
  fi
  gg_finalize_auto_busco_lineage_name "${auto_lineage}"
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
  if [[ -n "${codontable}" ]]; then
    seqkit replace --pattern X --replacement N --by-seq --ignore-case --threads "${threads}" \
    | seqkit replace --pattern " .*" --replacement "" --ignore-case --threads "${threads}" \
    | cdskit pad --codontable "${codontable}"
  else
    seqkit replace --pattern X --replacement N --by-seq --ignore-case --threads "${threads}" \
    | seqkit replace --pattern " .*" --replacement "" --ignore-case --threads "${threads}" \
    | cdskit pad
  fi
}

_download_busco_lineage_to_runtime() {
  local busco_lineage=$1
  local runtime_busco_db=$2
  local runtime_busco_lineage=$3
  if [[ -e "${runtime_busco_lineage}" ]]; then
    return 0
  fi
  echo "Starting BUSCO dataset download: ${busco_lineage}" >&2
  if ! busco --download "${busco_lineage}" >&2; then
    echo "BUSCO dataset download failed: ${busco_lineage}" >&2
    return 1
  fi
  if [[ -d busco_downloads ]]; then
    find busco_downloads -mindepth 1 -maxdepth 1 -exec mv -f -- {} "${runtime_busco_db}"/ \;
    rm -rf -- busco_downloads
  fi
  if [[ ! -e "${runtime_busco_lineage}" ]]; then
    echo "BUSCO lineage dataset is still missing after download: ${runtime_busco_lineage}" >&2
    return 1
  fi
  echo "BUSCO dataset download has been finished: ${busco_lineage}" >&2
}

ensure_busco_download_path() {
  local gg_workspace_dir=$1
  local busco_lineage=$2
  local sys_busco_db="/usr/local/db/busco_downloads"
  local sys_busco_lineage="${sys_busco_db}/lineages/${busco_lineage}"
  local runtime_busco_db
  local runtime_busco_lineage
  local lock_file

  if [[ -z "${busco_lineage}" ]]; then
    echo "ensure_busco_download_path: busco_lineage is empty." >&2
    return 1
  fi

  if [[ -e "${sys_busco_lineage}" ]]; then
    echo "${sys_busco_db}"
    return 0
  fi

  runtime_busco_db="$(workspace_downloads_root "${gg_workspace_dir}")/busco_downloads"
  runtime_busco_lineage="${runtime_busco_db}/lineages/${busco_lineage}"
  lock_file="${runtime_busco_db}/locks/busco_downloads.lock"
  ensure_dir "${runtime_busco_db}"
  ensure_dir "${runtime_busco_db}/locks"

  gg_array_download_once "${lock_file}" "${runtime_busco_lineage}" "BUSCO dataset download (${busco_lineage})" \
    _download_busco_lineage_to_runtime "${busco_lineage}" "${runtime_busco_db}" "${runtime_busco_lineage}" || return 1

  if [[ ! -e "${runtime_busco_lineage}" ]]; then
    echo "Failed to prepare BUSCO lineage dataset: ${busco_lineage}" >&2
    return 1
  fi

  echo "${runtime_busco_db}"
}

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
}

_ensure_ete_taxonomy_db_locked() {
  local db_file=$1
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
    os.makedirs(os.path.join(cache_dir, "ete4"), exist_ok=True)

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
  local lock_file

  dir_db=$(workspace_downloads_root "${gg_workspace_dir}")
  dir_taxonomy=$(workspace_taxonomy_root "${gg_workspace_dir}")
  db_file=$(workspace_taxonomy_dbfile "${gg_workspace_dir}")
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

  gg_array_download_once "${lock_file}" "${db_file}" "ETE taxonomy DB" \
    _ensure_ete_taxonomy_db_locked "${db_file}" || return 1
  if [[ ! -s "${db_file}" ]]; then
    echo "Failed to prepare ETE taxonomy DB: ${db_file}" >&2
    return 1
  fi
  return 0
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

cp_out() {
	if [[ $# -eq 1 ]]; then
		if [[ -p /dev/stdin ]]; then
			ensure_parent_dir "$1"
			cat > "$1"
			return $?
		fi
		echo "cp_out: at least 2 arguments are required unless stdin is piped."
		return 1
	fi
	if [[ $# -eq 2 && "$1" == "-" ]]; then
		ensure_parent_dir "$2"
		cat > "$2"
		return $?
	fi
	if [[ $# -lt 2 ]]; then
		echo "cp_out: at least 2 arguments are required."
		return 1
	fi
	local dest="${!#}"
	if [[ $# -gt 2 || "${dest}" == */ ]]; then
		ensure_dir "${dest%/}"
	else
		ensure_parent_dir "${dest}"
	fi
	cp -- "$@"
}

mv_out() {
	if [[ $# -eq 1 ]]; then
		if [[ -p /dev/stdin ]]; then
			ensure_parent_dir "$1"
			cat > "$1"
			return $?
		fi
		echo "mv_out: at least 2 arguments are required unless stdin is piped."
		return 1
	fi
	if [[ $# -eq 2 && "$1" == "-" ]]; then
		ensure_parent_dir "$2"
		cat > "$2"
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
  local runtime_bin
  if ! runtime_bin=$(gg_detect_container_runtime_binary); then
    echo "${echo_header}Neither singularity nor apptainer was found on PATH."
    return 1
  fi
  echo ${echo_header}"hostname = $(hostname)"
  echo ${echo_header}"container runtime = ${runtime_bin}"
  if declare -F gg_site_container_shell_command >/dev/null 2>&1; then
    gg_site_container_shell_command "${runtime_bin}" singularity_command || return 1
  else
    echo ${echo_header}"No site adapter was loaded. Using default shell."
    singularity_command="${runtime_bin} shell"
  fi
  echo ${echo_header}'${singularity_command}' = \"${singularity_command}\"
  echo ""
}

gg_entrypoint_prepare_container_runtime() {
	local call_exit_if_running=${1:-0}
	gg_scheduler_runtime_prelude
	unset_singularity_envs
	if [[ "${call_exit_if_running}" -eq 1 ]]; then
		exit_if_running_qstat
	fi
	if ! set_singularity_command; then
		return 1
	fi
	gg_normalize_scheduler_env
	return 0
}

gg_entrypoint_activate_container_runtime() {
	set_singularityenv
	gg_print_scheduler_runtime_summary
}

gg_entrypoint_enter_workspace() {
	mkdir -p "${gg_workspace_dir}"
	cd "${gg_workspace_dir}"
}

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

check_species_cds() {
  local gg_workspace_dir=$1
  local dir_sp_cds="$(workspace_input_root "${gg_workspace_dir}")/species_cds"
  local species_cds_fasta=()
  local fasta_path
  while IFS= read -r fasta_path; do
    species_cds_fasta+=( "${fasta_path}" )
  done < <(gg_find_fasta_files "${dir_sp_cds}" 1)
  if [[ ${#species_cds_fasta[@]} -eq 0 ]]; then
    echo "No species_cds fasta files were found in: ${dir_sp_cds}"
    exit 1
  fi
  local error_log=$(mktemp)  # Temporary file to collect errors
  echo "$(date): Started validating the format of all species_cds files using ${GG_TASK_CPUS} processes."

  function check_single_species_cds () {
    local spfasta=$1
    local sp_ub
    local first_header
    local first_header_no_gt
    local first_header_sp
    local spfasta_startswith
    sp_ub=$(basename "${spfasta}" | sed -e "s/_/|/" -e "s/_.*//" -e "s/|/_/")
    local seq_names
    seq_names=$(seqkit seq "${spfasta}" | awk '/^>/ {print}')
    first_header=${seq_names%%$'\n'*}
    first_header=${first_header%%[[:space:]]*}
    first_header_no_gt=${first_header#>}
    first_header_sp=$(printf '%s' "${first_header_no_gt}" | sed -e "s/_/|/" -e "s/_.*//" -e "s/|/_/")
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
  export error_log
  parallel --jobs "${GG_TASK_CPUS}" check_single_species_cds ::: "${species_cds_fasta[@]}"

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
  shopt -s nullglob dotglob
  for file in "${dir}"/*; do
    [[ -f "${file}" ]] || continue
    [[ "$(basename "${file}")" == .* ]] && continue
    files+=( "${file}" )
  done
  shopt -u nullglob dotglob
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
    echo "Species names: ${sp_names[@]}"
    exit 1
  else
    echo "Species files are unique in: ${species_dir}"
    echo "Species names: ${sp_names[@]}"
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
		ulimit -s unlimited
	fi
	if [[ "${print_start_message}" -eq 1 ]]; then
		print_gg_container_starting_message
	fi
}

gg_initialize_conda_shell() {
	if [[ "${GG_CONDA_SHELL_INITIALIZED:-0}" -eq 1 ]]; then
		return 0
	fi
	if command -v micromamba >/dev/null 2>&1; then
		export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-/opt/conda}"
		eval "$(micromamba shell hook --shell bash)"
		if ! declare -F conda >/dev/null 2>&1 && ! command -v conda >/dev/null 2>&1; then
			conda() {
				micromamba "$@"
			}
		fi
		GG_CONDA_SHELL_INITIALIZED=1
		export GG_CONDA_SHELL_INITIALIZED
		return 0
	fi
	if [[ -f /opt/conda/etc/profile.d/conda.sh ]]; then
		# shellcheck disable=SC1091
		source /opt/conda/etc/profile.d/conda.sh
		GG_CONDA_SHELL_INITIALIZED=1
		export GG_CONDA_SHELL_INITIALIZED
		return 0
	fi
	if command -v conda >/dev/null 2>&1; then
		eval "$(conda shell.bash hook)"
		GG_CONDA_SHELL_INITIALIZED=1
		export GG_CONDA_SHELL_INITIALIZED
		return 0
	fi
	return 1
}

gg_activate_conda_env() {
	local conda_env=${1:-}
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
	conda activate "${conda_env}"
}

gg_deactivate_conda_env() {
	if [[ "${GG_CONDA_SHELL_INITIALIZED:-0}" -ne 1 ]]; then
		return 0
	fi
	if declare -F conda >/dev/null 2>&1; then
		conda deactivate >/dev/null 2>&1 || true
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

gg_trigger_versions_dump() {
  local trigger_name=${1:-unknown_job}
  local versions_script
  local version_file
  local gg_version
  local dir_output
  local versions_dir
  local container_key_seed
  local container_key_hash
  local inspect_snapshot
  local image_file_hash
  local versions_script_hash=""
  local log_file
  local tmp_log_file
  local lock_file
  local failed_log_file
  local had_flock=0
  local singularity_bin
  local container_runtime_bin
  local versions_exit_code=0
  local block_exit_code=0
  local had_errexit=0

  if [[ -z "${gg_workflow_dir:-}" || -z "${gg_support_dir:-}" || -z "${gg_workspace_dir:-}" || -z "${gg_container_image_path:-}" ]]; then
    echo "gg_trigger_versions_dump: gg_workflow_dir/gg_support_dir/gg_workspace_dir/gg_container_image_path are required." >&2
    return 1
  fi

  version_file="${gg_workflow_dir}/../VERSION"
  gg_version="${SINGULARITYENV_GG_VERSION:-}"
  if [[ -z "${gg_version}" ]]; then
    gg_version="$(gg_read_repo_version "${version_file}")"
  fi
  export SINGULARITYENV_GG_VERSION="${gg_version}"
  export APPTAINERENV_GG_VERSION="${gg_version}"

  versions_script="${gg_support_dir}/gg_versions.sh"
  if [[ ! -s "${versions_script}" ]]; then
    echo "gg_trigger_versions_dump: versions script not found: ${versions_script}" >&2
    return 1
  fi

  if [[ -z "${singularity_command:-}" ]]; then
    set_singularity_command
  fi
  container_runtime_bin="${singularity_command%% *}"
  if [[ -z "${container_runtime_bin}" ]]; then
    if ! container_runtime_bin=$(gg_detect_container_runtime_binary); then
      echo "gg_trigger_versions_dump: container runtime command not found. Skipping version dump." >&2
      return 0
    fi
  fi
  if ! command -v "${container_runtime_bin}" >/dev/null 2>&1; then
    echo "gg_trigger_versions_dump: container runtime command not found (${container_runtime_bin}). Skipping version dump." >&2
    return 0
  fi

  if ! gg_container_bind_destination_exists "${gg_workspace_dir}:/workspace" \
    || ! gg_container_bind_destination_exists "${gg_workflow_dir}:/script"; then
    gg_normalize_scheduler_env
    set_singularityenv
  fi

  dir_output=$(workspace_output_root "${gg_workspace_dir}")
  versions_dir="${dir_output}/versions"
  ensure_dir "${versions_dir}"

  container_key_seed="gg_container_image_path=${gg_container_image_path};runtime=${container_runtime_bin};gg_version=${gg_version}"
  if [[ -s "${versions_script}" ]]; then
    if command -v sha256sum >/dev/null 2>&1; then
      versions_script_hash=$(sha256sum "${versions_script}" | awk '{print $1}')
      container_key_seed="${container_key_seed};versions_script_sha256=${versions_script_hash}"
    elif command -v shasum >/dev/null 2>&1; then
      versions_script_hash=$(shasum -a 256 "${versions_script}" | awk '{print $1}')
      container_key_seed="${container_key_seed};versions_script_sha256=${versions_script_hash}"
    else
      versions_script_hash=$(cksum "${versions_script}" | awk '{print $1 "-" $2}')
      container_key_seed="${container_key_seed};versions_script_cksum=${versions_script_hash}"
    fi
  fi
  if [[ -s "${gg_container_image_path}" ]]; then
    if command -v sha256sum >/dev/null 2>&1; then
      image_file_hash=$(sha256sum "${gg_container_image_path}" | awk '{print $1}')
      container_key_seed="${container_key_seed};image_sha256=${image_file_hash}"
    elif command -v shasum >/dev/null 2>&1; then
      image_file_hash=$(shasum -a 256 "${gg_container_image_path}" | awk '{print $1}')
      container_key_seed="${container_key_seed};image_sha256=${image_file_hash}"
    else
      image_file_hash=$(cksum "${gg_container_image_path}" | awk '{print $1 "-" $2}')
      container_key_seed="${container_key_seed};image_cksum=${image_file_hash}"
    fi
  fi
  inspect_snapshot=""
  if inspect_snapshot=$("${container_runtime_bin}" inspect "${gg_container_image_path}" 2>/dev/null); then
    if [[ -n "${inspect_snapshot}" ]]; then
      container_key_seed="${container_key_seed};inspect=${inspect_snapshot}"
    fi
  fi
  if command -v sha256sum >/dev/null 2>&1; then
    container_key_hash=$(printf '%s' "${container_key_seed}" | sha256sum | awk '{print $1}')
  elif command -v shasum >/dev/null 2>&1; then
    container_key_hash=$(printf '%s' "${container_key_seed}" | shasum -a 256 | awk '{print $1}')
  else
    container_key_hash=$(printf '%s' "${container_key_seed}" | cksum | awk '{print $1}')
  fi
  if [[ -z "${container_key_hash}" ]]; then
    container_key_hash=$(echo "${gg_container_image_path}" | tr '[:space:]/:' '_' | tr -cd '[:alnum:]_.-')
    if [[ -z "${container_key_hash}" ]]; then
      container_key_hash="unknown_container"
    fi
  fi
  log_file="${versions_dir}/container.${container_key_hash}.versions.log"
  tmp_log_file="${log_file}.tmp.$$"
  lock_file="${log_file}.lock"

  if command -v flock >/dev/null 2>&1; then
    exec 9>"${lock_file}"
    flock 9
    had_flock=1
  fi
  if [[ -s "${log_file}" ]]; then
    echo "gg_trigger_versions_dump: skipped existing ${log_file}"
    if [[ ${had_flock} -eq 1 ]]; then
      flock -u 9
      exec 9>&-
    fi
    return 0
  fi

  if [[ $- == *e* ]]; then
    had_errexit=1
    set +e
  fi
  {
    echo "### genegalleon version ###"
    echo "${gg_version}"
    echo ""
    echo "$(date): Triggered gg_versions by ${trigger_name}"
    ${singularity_command} "${gg_container_image_path}" < "${versions_script}" || {
      cmd_rc=$?
      if [[ ${versions_exit_code} -eq 0 ]]; then
        versions_exit_code=${cmd_rc}
      fi
    }
    echo ""
    echo "### gg_container ###"
    singularity_bin="$(command -v "${container_runtime_bin}")"
    if [[ "${singularity_bin}" == *"/gg_wrapper_bin/"* ]]; then
      echo "Skipping container inspect/version under Docker-backed singularity shim: ${singularity_bin}"
    else
      "${container_runtime_bin}" inspect "${gg_container_image_path}" || {
        cmd_rc=$?
        if [[ ${versions_exit_code} -eq 0 ]]; then
          versions_exit_code=${cmd_rc}
        fi
      }
      echo ""
      echo "### ${container_runtime_bin} version ###"
      "${container_runtime_bin}" version || {
        cmd_rc=$?
        if [[ ${versions_exit_code} -eq 0 ]]; then
          versions_exit_code=${cmd_rc}
        fi
      }
    fi
    echo ""
    echo "### Host OS info ###"
    if [[ -f /etc/os-release ]]; then
      cat /etc/os-release
    else
      echo "/etc/os-release was not found. Falling back to uname output."
      uname -a
    fi
    echo ""
    echo "$(date): gg_versions trigger completed"
  } > "${tmp_log_file}" 2>&1
  block_exit_code=$?
  if [[ ${block_exit_code} -ne 0 && ${versions_exit_code} -eq 0 ]]; then
    versions_exit_code=${block_exit_code}
  fi
  if [[ ${had_errexit} -eq 1 ]]; then
    set -e
  fi

  if [[ ${versions_exit_code} -ne 0 ]]; then
    failed_log_file="${versions_dir}/container.${container_key_hash}.versions.failed.$(date '+%Y%m%d_%H%M%S').log"
    if [[ -s "${tmp_log_file}" ]]; then
      mv_out "${tmp_log_file}" "${failed_log_file}"
    else
      : > "${failed_log_file}"
    fi
    if [[ ${had_flock} -eq 1 ]]; then
      flock -u 9
      exec 9>&-
    fi
    echo "gg_trigger_versions_dump: failed (exit=${versions_exit_code}). Log: ${failed_log_file}" >&2
    return "${versions_exit_code}"
  fi

  if [[ -s "${tmp_log_file}" ]]; then
    mv_out "${tmp_log_file}" "${log_file}"
  fi
  if [[ ${had_flock} -eq 1 ]]; then
    flock -u 9
    exec 9>&-
  fi
  echo "gg_trigger_versions_dump: wrote ${log_file}"
  return 0
}

gg_require_versions_dump() {
  local trigger_name=${1:-unknown_job}
  local versions_exit_code=0
  local had_errexit=0

  if [[ $- == *e* ]]; then
    had_errexit=1
    set +e
  fi
  gg_trigger_versions_dump "${trigger_name}"
  versions_exit_code=$?
  if [[ ${had_errexit} -eq 1 ]]; then
    set -e
  fi
  if [[ ${versions_exit_code} -ne 0 ]]; then
    echo "gg_require_versions_dump: gg_versions trigger failed for ${trigger_name}." >&2
    return "${versions_exit_code}"
  fi
  return 0
}

enable_all_run_flags_for_debug_mode() {
  local message=${1:-"gg debug mode: All run_* variables are forced to set 1."}
  local debug_mode=${gg_debug_mode:-0}
  local flags=()
  local flag

  if [[ "${debug_mode}" -ne 1 ]]; then
    return 0
  fi

  echo "${message}"
  while IFS= read -r flag; do
    flags+=( "${flag}" )
  done < <(set | sed -n -E 's/^(run_[A-Za-z0-9_]+)=.*/\1/p')
  for flag in "${flags[@]}"; do
    printf -v "${flag}" '%s' 1
    export "${flag}"
  done
  return 0
}

ensure_gg_scheduler_defaults() {
  local echo_header="ensure_gg_scheduler_defaults: "
  if [[ -z "${GG_TASK_CPUS:-}" ]]; then
    echo "${echo_header}GG_TASK_CPUS is undefined or empty. GG_TASK_CPUS=1"
    GG_TASK_CPUS=1
  fi
  if [[ -z "${GG_JOB_ID:-}" ]]; then
    echo "${echo_header}GG_JOB_ID is undefined or empty. GG_JOB_ID=1"
    GG_JOB_ID=1
  fi
  if [[ -z "${GG_ARRAY_TASK_ID:-}" ]]; then
    echo "${echo_header}GG_ARRAY_TASK_ID is undefined or empty. GG_ARRAY_TASK_ID=1"
    GG_ARRAY_TASK_ID=1
  fi
  if [[ -z "${GG_MEM_PER_CPU_GB:-}" ]]; then
    echo "${echo_header}GG_MEM_PER_CPU_GB is undefined or empty. GG_MEM_PER_CPU_GB=3"
    GG_MEM_PER_CPU_GB=3
  fi
  if [[ -z "${GG_MEM_TOTAL_GB:-}" ]]; then
    GG_MEM_TOTAL_GB=$((GG_MEM_PER_CPU_GB * GG_TASK_CPUS))
    echo "${echo_header}GG_MEM_TOTAL_GB is undefined or empty. GG_MEM_TOTAL_GB=${GG_MEM_TOTAL_GB}"
  fi
  gg_sync_legacy_scheduler_aliases
}

ensure_scheduler_defaults() {
  ensure_gg_scheduler_defaults "$@"
}

print_gg_container_starting_message() {
  local gg_workspace_input_dir_resolved
  local gg_workspace_output_dir_resolved
  local gg_workspace_downloads_dir_resolved
  local gg_workspace_layout_local
  ensure_gg_scheduler_defaults
  gg_workspace_layout_local=$(gg_resolve_workspace_layout "${gg_workspace_dir}")
  gg_workspace_input_dir_resolved=$(workspace_input_root "${gg_workspace_dir}")
  gg_workspace_output_dir_resolved=$(workspace_output_root "${gg_workspace_dir}")
  gg_workspace_downloads_dir_resolved=$(workspace_downloads_root "${gg_workspace_dir}")
  echo "$(date): Starting gg Singularity/Apptainer environment"
  echo "pwd: $(pwd)"
  echo "gg_workspace_dir: ${gg_workspace_dir}"
  echo "gg_workspace_layout: ${gg_workspace_layout_local}"
  echo "gg_workspace_input_dir: ${gg_workspace_input_dir_resolved}"
  echo "gg_workspace_output_dir: ${gg_workspace_output_dir_resolved}"
  echo "gg_workspace_downloads_dir: ${gg_workspace_downloads_dir_resolved}"
  echo "python: $(command -v python)"
  echo "GG_TASK_CPUS: ${GG_TASK_CPUS}"
  echo "GG_MEM_PER_CPU_GB: ${GG_MEM_PER_CPU_GB}"
  echo "GG_MEM_TOTAL_GB: ${GG_MEM_TOTAL_GB}"
  echo "GG_JOB_ID: ${GG_JOB_ID}"
  echo "GG_ARRAY_TASK_ID: ${GG_ARRAY_TASK_ID}"
  echo "ulimit -Hn: $(ulimit -Hn)"
  echo "ulimit -Sn: $(ulimit -Sn)"
}

recreate_dir() {
  local dir=$1
  if [[ -z "${dir}" || "${dir}" == "/" ]]; then
    echo "Refusing to recreate unsafe directory path: ${dir}"
    return 1
  fi
  if [[ -e "${dir}" ]]; then
    echo "Deleting: ${dir}"
    rm -rf -- "${dir}"
  fi
  echo "Creating: ${dir}"
  mkdir -p "${dir}"
}

gg_lock_stale_seconds() {
  local stale_seconds="${GG_LOCK_STALE_SECONDS:-86400}"
  if [[ ! "${stale_seconds}" =~ ^[0-9]+$ ]]; then
    stale_seconds=86400
  fi
  if (( stale_seconds < 1 )); then
    stale_seconds=1
  fi
  echo "${stale_seconds}"
}

gg_stat_mtime_epoch() {
  local target=$1
  if [[ ! -e "${target}" ]]; then
    echo ""
    return 0
  fi
  if stat --version >/dev/null 2>&1; then
    stat -c '%Y' "${target}" 2>/dev/null || true
  else
    stat -f '%m' "${target}" 2>/dev/null || true
  fi
}

gg_lock_pid_is_alive() {
  local pid=$1
  if [[ ! "${pid}" =~ ^[0-9]+$ ]]; then
    return 1
  fi
  kill -0 "${pid}" 2>/dev/null
}

gg_lock_marker_path() {
  local lock_file=$1
  echo "${lock_file}.owner"
}

gg_maybe_recover_stale_lock_marker() {
  local lock_file=$1
  local description=$2
  local marker_file
  marker_file=$(gg_lock_marker_path "${lock_file}")
  if [[ ! -e "${marker_file}" ]]; then
    return 0
  fi

  local marker_pid=""
  local marker_start_ts=""
  local marker_mtime=""
  local stale_seconds
  stale_seconds=$(gg_lock_stale_seconds)
  local now_epoch
  now_epoch=$(date +%s)

  if [[ -s "${marker_file}" ]]; then
    IFS=$'\t' read -r marker_pid marker_start_ts _ < "${marker_file}" || true
  fi
  if [[ ! "${marker_start_ts}" =~ ^[0-9]+$ ]]; then
    marker_mtime=$(gg_stat_mtime_epoch "${marker_file}")
    if [[ "${marker_mtime}" =~ ^[0-9]+$ ]]; then
      marker_start_ts=${marker_mtime}
    else
      marker_start_ts=${now_epoch}
    fi
  fi

  local marker_age=$(( now_epoch - marker_start_ts ))
  if (( marker_age < 0 )); then
    marker_age=0
  fi
  local stale_reason=""
  if [[ -n "${marker_pid}" ]]; then
    if gg_lock_pid_is_alive "${marker_pid}"; then
      stale_reason=""
    else
      stale_reason="owner_not_running"
    fi
  else
    if (( marker_age >= stale_seconds )); then
      stale_reason="owner_unknown_timeout"
    fi
  fi
  if [[ -n "${stale_reason}" ]]; then
    rm -f -- "${marker_file}"
    echo "Recovered stale lock marker: ${description} (${stale_reason}, age=${marker_age}s, marker=${marker_file})" >&2
  fi
}

gg_write_lock_marker() {
  local lock_file=$1
  local description=$2
  local marker_file
  marker_file=$(gg_lock_marker_path "${lock_file}")
  printf '%s\t%s\t%s\n' "$$" "$(date +%s)" "${description}" > "${marker_file}"
}

gg_remove_lock_marker() {
  local lock_file=$1
  local marker_file
  marker_file=$(gg_lock_marker_path "${lock_file}")
  rm -f -- "${marker_file}"
}

gg_safe_remove_lock_dir() {
  local lock_dir=$1
  if [[ ! -e "${lock_dir}" ]]; then
    return 0
  fi
  if [[ -z "${lock_dir}" || "${lock_dir}" == "/" || "${lock_dir}" == "." ]]; then
    return 1
  fi
  if [[ "${lock_dir}" != *.dlock ]]; then
    return 1
  fi
  rm -rf -- "${lock_dir}"
}

gg_acquire_mkdir_lock() {
  local lock_dir=$1
  local description=$2
  local owner_file="${lock_dir}/owner"
  local stale_seconds
  stale_seconds=$(gg_lock_stale_seconds)

  while true; do
    if mkdir "${lock_dir}" 2>/dev/null; then
      printf '%s\t%s\t%s\n' "$$" "$(date +%s)" "${description}" > "${owner_file}"
      return 0
    fi

    local marker_pid=""
    local marker_start_ts=""
    local marker_mtime=""
    local now_epoch
    now_epoch=$(date +%s)

    if [[ -s "${owner_file}" ]]; then
      IFS=$'\t' read -r marker_pid marker_start_ts _ < "${owner_file}" || true
    fi
    if [[ ! "${marker_start_ts}" =~ ^[0-9]+$ ]]; then
      marker_mtime=$(gg_stat_mtime_epoch "${lock_dir}")
      if [[ "${marker_mtime}" =~ ^[0-9]+$ ]]; then
        marker_start_ts=${marker_mtime}
      else
        marker_start_ts=${now_epoch}
      fi
    fi

    local marker_age=$(( now_epoch - marker_start_ts ))
    if (( marker_age < 0 )); then
      marker_age=0
    fi
    local stale_reason=""
    if [[ -n "${marker_pid}" ]]; then
      if gg_lock_pid_is_alive "${marker_pid}"; then
        stale_reason=""
      else
        stale_reason="owner_not_running"
      fi
    else
      if (( marker_age >= stale_seconds )); then
        stale_reason="owner_unknown_timeout"
      fi
    fi

    if [[ -n "${stale_reason}" ]]; then
      if gg_safe_remove_lock_dir "${lock_dir}"; then
        echo "Recovered stale lock directory: ${description} (${stale_reason}, age=${marker_age}s, lock=${lock_dir})" >&2
        continue
      fi
    fi
    sleep 1
  done
}

gg_release_mkdir_lock() {
  local lock_dir=$1
  gg_safe_remove_lock_dir "${lock_dir}" || true
}

gg_array_download_once() {
  local lock_file=$1
  local done_file=$2
  local description=$3
  shift 3
  local task_id=${GG_ARRAY_TASK_ID:-1}
  local has_flock=1
  local lock_dir="${lock_file}.dlock"

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
    gg_maybe_recover_stale_lock_marker "${lock_file}" "${description}"
    gg_write_lock_marker "${lock_file}" "${description}"

    if [[ -s "${done_file}" ]]; then
      gg_remove_lock_marker "${lock_file}"
      flock -u 9
      exec 9>&-
      return 0
    fi

    echo "GG_ARRAY_TASK_ID=${task_id}: starting download: ${description}" >&2
    "$@"
    local download_exit_code=$?
    if [[ ${download_exit_code} -ne 0 ]]; then
      echo "GG_ARRAY_TASK_ID=${task_id}: download failed: ${description}" >&2
      gg_remove_lock_marker "${lock_file}"
      flock -u 9
      exec 9>&-
      return ${download_exit_code}
    fi
    if [[ ! -s "${done_file}" ]]; then
      echo "Downloaded DB file not found after synchronization: ${done_file}" >&2
      gg_remove_lock_marker "${lock_file}"
      flock -u 9
      exec 9>&-
      return 1
    fi

    echo "GG_ARRAY_TASK_ID=${task_id}: download completed: ${description}" >&2
    gg_remove_lock_marker "${lock_file}"
    flock -u 9
    exec 9>&-
    return 0
  fi

  gg_acquire_mkdir_lock "${lock_dir}" "${description}" || return 1
  if [[ -s "${done_file}" ]]; then
    gg_release_mkdir_lock "${lock_dir}"
    return 0
  fi
  echo "GG_ARRAY_TASK_ID=${task_id}: starting download (mkdir lock fallback): ${description}" >&2
  "$@"
  local download_exit_code=$?
  if [[ ${download_exit_code} -ne 0 ]]; then
    echo "GG_ARRAY_TASK_ID=${task_id}: download failed (mkdir lock fallback): ${description}" >&2
    gg_release_mkdir_lock "${lock_dir}"
    return ${download_exit_code}
  fi
  if [[ ! -s "${done_file}" ]]; then
    echo "Downloaded DB file not found after synchronization: ${done_file}" >&2
    gg_release_mkdir_lock "${lock_dir}"
    return 1
  fi

  echo "GG_ARRAY_TASK_ID=${task_id}: download completed (mkdir lock fallback): ${description}" >&2
  gg_release_mkdir_lock "${lock_dir}"
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
  local qstat_token
  while IFS= read -r qstat_token; do
    qstat_task_tokens+=( "${qstat_token}" )
  done < <(qstat 2>/dev/null | awk -v jid="${job_id}" 'NR>2 && $1==jid {print $NF}')
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
  local id_candidate="${SLURM_JOB_ID:-${GG_JOB_ID:-}}"
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
  local slurm_token
  while IFS= read -r slurm_token; do
    slurm_job_tokens+=( "${slurm_token}" )
  done < <(squeue -h -j "${array_job_id}" -o "%i" 2>/dev/null)
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
  touch "${done_file}"
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
  touch "${done_file}"
}

ensure_uniprot_sprot_db() {
  local gg_workspace_dir=$1
  local sys_prefix="/usr/local/db/uniprot_sprot"
  local runtime_root="$(workspace_downloads_root "${gg_workspace_dir}")"
  local runtime_dir="${runtime_root}/uniprot_sprot"
  local runtime_prefix="${runtime_dir}/uniprot_sprot"
  local lock_file="${runtime_root}/locks/uniprot_sprot.lock"

  if [[ -s "${sys_prefix}.pep" && -s "${sys_prefix}.dmnd" ]]; then
    echo "${sys_prefix}"
    return 0
  fi
  if [[ -s "${runtime_prefix}.pep" && -s "${runtime_prefix}.dmnd" ]]; then
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
  echo "${runtime_prefix}"
}

ensure_uniprot_sprot_mmseqs_db() {
  local gg_workspace_dir=$1
  local sys_prefix="/usr/local/db/uniprot_sprot"
  local runtime_root="$(workspace_downloads_root "${gg_workspace_dir}")"
  local runtime_dir="${runtime_root}/uniprot_sprot"
  local runtime_prefix="${runtime_dir}/uniprot_sprot"
  local runtime_pep="${runtime_prefix}.pep"
  local runtime_mmseqs_db="${runtime_prefix}.mmseqs"
  local runtime_mmseqs_dbtype="${runtime_mmseqs_db}.dbtype"
  local runtime_ready="${runtime_mmseqs_db}.ready"
  local lock_file_pep="${runtime_root}/locks/uniprot_sprot.pep.lock"
  local lock_file_mmseqs="${runtime_root}/locks/uniprot_sprot.mmseqs.lock"

  if [[ -s "${sys_prefix}.pep" && -s "${sys_prefix}.mmseqs" && -s "${sys_prefix}.mmseqs.dbtype" ]]; then
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
    touch "${runtime_ready}"
    echo "${runtime_prefix}"
    return 0
  fi

  gg_array_download_once "${lock_file_mmseqs}" "${runtime_ready}" "UniProt Swiss-Prot MMseqs2 DB" \
    _build_mmseqs_db_from_fasta "${runtime_pep}" "${runtime_mmseqs_db}" "${runtime_ready}" || return 1

  if [[ ! -s "${runtime_mmseqs_db}" || ! -s "${runtime_mmseqs_dbtype}" || ! -s "${runtime_ready}" ]]; then
    echo "MMseqs2 UniProt Swiss-Prot DB download/build failed." >&2
    return 1
  fi
  echo "${runtime_prefix}"
}

ensure_uniprot_sprot_blast_db() {
  local gg_workspace_dir=$1
  local sys_prefix="/usr/local/db/uniprot_sprot"
  local runtime_root="$(workspace_downloads_root "${gg_workspace_dir}")"
  local runtime_dir="${runtime_root}/uniprot_sprot"
  local runtime_prefix="${runtime_dir}/uniprot_sprot"
  local runtime_pep="${runtime_prefix}.pep"
  local runtime_blast_ready="${runtime_prefix}.blast.ready"
  local lock_file_pep="${runtime_root}/locks/uniprot_sprot.pep.lock"
  local lock_file_blast="${runtime_root}/locks/uniprot_sprot.blast.lock"

  if [[ -s "${sys_prefix}.pep" && -s "${sys_prefix}.pin" && -s "${sys_prefix}.phr" && -s "${sys_prefix}.psq" ]]; then
    echo "${sys_prefix}"
    return 0
  fi
  if [[ -s "${runtime_pep}" && -s "${runtime_prefix}.pin" && -s "${runtime_prefix}.phr" && -s "${runtime_prefix}.psq" ]]; then
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
    touch "${runtime_blast_ready}"
    echo "${runtime_prefix}"
    return 0
  fi

  gg_array_download_once "${lock_file_blast}" "${runtime_blast_ready}" "UniProt Swiss-Prot BLASTP DB" \
    _build_blast_db_from_fasta "${runtime_pep}" "${runtime_prefix}" "${runtime_blast_ready}" || return 1

  if [[ ! -s "${runtime_prefix}.pin" || ! -s "${runtime_prefix}.phr" || ! -s "${runtime_prefix}.psq" || ! -s "${runtime_blast_ready}" ]]; then
    echo "BLASTP UniProt Swiss-Prot DB download/build failed." >&2
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
  rm -rf -- "${tmp_dir}"
}

pfam_le_db_is_ready() {
  local db_dir=$1
  [[ -s "${db_dir}/Pfam.pal" || ( -s "${db_dir}/Pfam.loo" && -s "${db_dir}/Pfam.rps" ) ]]
}

ensure_pfam_le_db() {
  local gg_workspace_dir=$1
  local sys_dir="/usr/local/db/Pfam_LE"
  local runtime_root="$(workspace_downloads_root "${gg_workspace_dir}")"
  local runtime_parent
  runtime_parent=$(workspace_pfam_root "${gg_workspace_dir}")
  local runtime_dir
  runtime_dir=$(workspace_pfam_le_dir "${gg_workspace_dir}")
  local lock_file="${runtime_root}/locks/pfam_le.lock"

  if pfam_le_db_is_ready "${sys_dir}"; then
    echo "${sys_dir}"
    return 0
  fi

  if pfam_le_db_is_ready "${runtime_dir}"; then
    echo "${runtime_dir}"
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
  local runtime_root="$(workspace_downloads_root "${gg_workspace_dir}")"
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

ensure_latest_jaspar_file() {
  local gg_workspace_dir=$1
  local runtime_root="$(workspace_downloads_root "${gg_workspace_dir}")"
  local runtime_dir="${runtime_root}/jaspar"
  local sys_dir="/usr/local/db/jaspar"
  local lock_file="${runtime_root}/locks/jaspar_latest.lock"
  local lock_dir="${lock_file}.dlock"
  local latest_marker="${runtime_dir}/latest_core_plants_non-redundant_pfms_meme.filename"
  local has_flock=1
  local resolved_filename=""
  local resolved_path=""
  local ensure_exit_code=0

  mkdir -p "${runtime_root}" "${runtime_dir}" "$(dirname "${lock_file}")"

  if ! command -v flock >/dev/null 2>&1; then
    has_flock=0
    echo "flock command was not found. Using mkdir lock fallback synchronization: latest JASPAR motif file" >&2
  fi
  if [[ ${has_flock} -eq 1 ]]; then
    exec 8> "${lock_file}"
    flock 8
    gg_maybe_recover_stale_lock_marker "${lock_file}" "latest JASPAR motif file"
    gg_write_lock_marker "${lock_file}" "latest JASPAR motif file"
  else
    gg_acquire_mkdir_lock "${lock_dir}" "latest JASPAR motif file" || return 1
  fi

  if [[ -s "${latest_marker}" ]]; then
    read -r resolved_filename < "${latest_marker}"
    if ! _jaspar_is_plants_core_meme_filename "${resolved_filename}"; then
      resolved_filename=""
    fi
    if [[ -n "${resolved_filename}" ]]; then
      if [[ -s "${sys_dir}/${resolved_filename}" ]]; then
        resolved_path="${sys_dir}/${resolved_filename}"
      elif [[ -s "${runtime_dir}/${resolved_filename}" ]]; then
        resolved_path="${runtime_dir}/${resolved_filename}"
      fi
    fi
  fi

  if [[ -z "${resolved_path}" ]]; then
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
      ensure_exit_code=1
    else
      if resolved_path=$(_ensure_jaspar_file_named "${gg_workspace_dir}" "${resolved_filename}"); then
        ensure_exit_code=0
      else
        ensure_exit_code=$?
      fi
      if [[ ${ensure_exit_code} -eq 0 ]]; then
        printf '%s\n' "${resolved_filename}" > "${latest_marker}.tmp"
        mv -- "${latest_marker}.tmp" "${latest_marker}"
      fi
    fi
  fi

  if [[ ${has_flock} -eq 1 ]]; then
    gg_remove_lock_marker "${lock_file}"
    flock -u 8
    exec 8>&-
  else
    gg_release_mkdir_lock "${lock_dir}"
  fi

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
  local runtime_root="$(workspace_downloads_root "${gg_workspace_dir}")"
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

gg_find_fasta_files() {
  local search_dir=$1
  local maxdepth=${2:-1}
  if [[ -z "${search_dir}" || ! -d "${search_dir}" ]]; then
    return 0
  fi
  find "${search_dir}" -maxdepth "${maxdepth}" -type f ! -name '.*' \
  \( -name "*.fa" -o -name "*.fa.gz" -o -name "*.fas" -o -name "*.fas.gz" -o -name "*.fasta" -o -name "*.fasta.gz" -o -name "*.fna" -o -name "*.fna.gz" \) \
  | sort
}

gg_find_file_basenames() {
  local search_dir=$1
  local name_pattern=${2:-*}
  local maxdepth=${3:-1}
  if [[ -z "${search_dir}" || ! -d "${search_dir}" ]]; then
    return 0
  fi
  find "${search_dir}" -maxdepth "${maxdepth}" -type f ! -name '.*' -name "${name_pattern}" \
  | awk -F'/' '{print $NF}' \
  | sort
}

is_species_set_identical() {
  local return_flag=0
  local dir1=$1
  local dir2=$2
  if [[ ! -d "${dir1}" || ! -d "${dir2}" ]]; then
    echo "Directory not found for species-set comparison: ${dir1} / ${dir2}"
    return 1
  fi
  local files1=()
  local files2=()
  local f
  shopt -s nullglob dotglob
  for f in "${dir1}"/*; do
    [[ -f "${f}" ]] || continue
    [[ "$(basename "${f}")" == .* ]] && continue
    files1+=( "${f}" )
  done
  for f in "${dir2}"/*; do
    [[ -f "${f}" ]] || continue
    [[ "$(basename "${f}")" == .* ]] && continue
    files2+=( "${f}" )
  done
  shopt -u nullglob dotglob
  local num_sp1=${#files1[@]}
  local num_sp2=${#files2[@]}
  if [[ ${num_sp1} -ne ${num_sp2} ]]; then
    echo "Number of species in ${dir1} and ${dir2} are different."
    return_flag=1
  fi
  local sp1=()
  local sp2=()
  for f in "${files1[@]}"; do
    sp1+=( "$(gg_species_name_from_path_or_dot "${f}")" )
  done
  for f in "${files2[@]}"; do
    sp2+=( "$(gg_species_name_from_path_or_dot "${f}")" )
  done
  local sp1_unique=()
  local sp2_unique=()
  local sp_name
  while IFS= read -r sp_name; do
    sp1_unique+=( "${sp_name}" )
  done < <(printf '%s\n' "${sp1[@]}" | sort -u)
  while IFS= read -r sp_name; do
    sp2_unique+=( "${sp_name}" )
  done < <(printf '%s\n' "${sp2[@]}" | sort -u)
  sp1=( "${sp1_unique[@]}" )
  sp2=( "${sp2_unique[@]}" )
  if [[ "${sp1[*]}" != "${sp2[*]}" ]]; then
    echo "Species names are different between ${dir1} and ${dir2}"
    echo "dir1 (${dir1}): ${sp1[*]}"
    echo "dir2 (${dir2}): ${sp2[*]}"
    return_flag=1
  fi
  return "${return_flag}"
}
