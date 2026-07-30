# shellcheck shell=bash
# Sourced by workflow/support/gg_util.sh.

gg_docker_singularity_shim_source_path() {
	local support_dir="${gg_support_dir:-}"

	if [[ -z "${support_dir}" ]]; then
		support_dir="${GG_UTIL_SUPPORT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)}"
	fi
	printf '%s\n' "${support_dir}/gg_wrapper_bin/singularity"
}

gg_docker_singularity_shim_path() {
	local wrapper_bin="${GG_WRAPPER_BIN:-/tmp/gg_wrapper_bin}"
	local source_shim=""
	local target_shim=""

	source_shim="$(gg_docker_singularity_shim_source_path)"
	if [[ ! -x "${source_shim}" ]]; then
		echo "Docker runtime shim source is missing or not executable: ${source_shim}" >&2
		return 1
	fi

	target_shim="${wrapper_bin}/singularity"
	mkdir -p "${wrapper_bin}"
	if ! ln -sf "${source_shim}" "${target_shim}" 2>/dev/null; then
		cp -- "${source_shim}" "${target_shim}"
		chmod +x "${target_shim}"
	fi

	printf '%s\n' "${target_shim}"
}

gg_default_docker_image_candidates() {
	printf '%s\n' "ghcr.io/kfuku52/genegalleon:latest"
	printf '%s\n' "local/genegalleon:dev"
}

gg_docker_image_exists() {
	local image_ref=${1:-}

	if [[ -z "${image_ref}" ]]; then
		return 1
	fi
	if ! command -v docker >/dev/null 2>&1; then
		return 1
	fi
	docker image inspect "${image_ref}" >/dev/null 2>&1
}

gg_resolve_docker_image_ref() {
	local image_ref="${GG_CONTAINER_DOCKER_IMAGE:-${GG_WRAPPER_IMAGE:-}}"
	local candidate=""

	if [[ -n "${image_ref}" ]]; then
		printf '%s\n' "${image_ref}"
		return 0
	fi

	while IFS= read -r candidate; do
		[[ -n "${candidate}" ]] || continue
		if gg_docker_image_exists "${candidate}"; then
			printf '%s\n' "${candidate}"
			return 0
		fi
	done < <(gg_default_docker_image_candidates)

	return 1
}

gg_auto_enable_docker_runtime_if_available() {
	local resolved_image_ref=""

	if [[ -z "${GG_CONTAINER_DOCKER_IMAGE:-}" && -z "${GG_WRAPPER_IMAGE:-}" ]]; then
		if [[ -n "${gg_container_image_path:-}" && -s "${gg_container_image_path}" ]]; then
			return 1
		fi
	fi
	if ! resolved_image_ref="$(gg_resolve_docker_image_ref)"; then
		return 1
	fi
	export GG_CONTAINER_DOCKER_IMAGE="${resolved_image_ref}"
	export GG_WRAPPER_IMAGE="${resolved_image_ref}"
	return 0
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

gg_entrypoint_env_override_prefix() {
	local entrypoint_name=${1:-}
	local entrypoint_base=""
	local entrypoint_stem=""

	entrypoint_base="$(basename "${entrypoint_name}")"
	case "${entrypoint_base}" in
		gg_input_generation_entrypoint.sh)
			printf '%s\n' "GG_INPUT_"
			return 0
			;;
		gg_transcriptome_generation_entrypoint.sh)
			printf '%s\n' "GG_TRANSCRIPTOME_"
			return 0
			;;
		gg_genome_annotation_entrypoint.sh)
			printf '%s\n' "GG_GENOME_ANNOTATION_"
			return 0
			;;
		gg_genome_evolution_entrypoint.sh)
			printf '%s\n' "GG_GENOME_EVOLUTION_"
			return 0
			;;
		gg_gene_evolution_entrypoint.sh)
			printf '%s\n' "GG_GENE_EVOLUTION_"
			return 0
			;;
		gg_gene_summary_entrypoint.sh)
			printf '%s\n' "GG_GENE_SUMMARY_"
			return 0
			;;
		gg_progress_summary_entrypoint.sh)
			printf '%s\n' "GG_PROGRESS_SUMMARY_"
			return 0
			;;
	esac

	entrypoint_stem="${entrypoint_base%.sh}"
	entrypoint_stem="${entrypoint_stem#gg_}"
	entrypoint_stem="${entrypoint_stem%_entrypoint}"
	entrypoint_stem="$(printf '%s' "${entrypoint_stem}" | tr '[:lower:]' '[:upper:]')"
	if [[ -z "${entrypoint_stem}" ]]; then
		return 1
	fi
	printf 'GG_%s_\n' "${entrypoint_stem}"
}

gg_config_var_env_override_name() {
	local env_prefix=${1:-}
	local var_name=${2:-}
	local upper_var_name=""

	if [[ -z "${env_prefix}" || -z "${var_name}" ]]; then
		return 1
	fi
	upper_var_name="$(printf '%s' "${var_name}" | tr '[:lower:]' '[:upper:]')"
	printf '%s%s\n' "${env_prefix}" "${upper_var_name}"
}

gg_apply_env_override_to_config_var() {
	local env_prefix=${1:-}
	local raw_var_name=${2:-}
	local var_name=""
	local env_name=""
	local env_value=""

	var_name=$(gg_normalize_config_var_name "${raw_var_name}" || true)
	if [[ -z "${var_name}" ]]; then
		return 0
	fi
	env_name=$(gg_config_var_env_override_name "${env_prefix}" "${var_name}") || return 1
	if [[ -n "${!env_name+x}" ]]; then
		env_value="${!env_name}"
		printf -v "${var_name}" '%s' "${env_value}"
	fi
}

gg_apply_registered_env_overrides() {
	local job_script=$1
	shift || true
	local entrypoint_name=""
	local env_prefix=""
	local raw_var_name=""

	entrypoint_name="$(basename "${job_script}")"
	env_prefix="$(gg_entrypoint_env_override_prefix "${entrypoint_name}")" || {
		echo "gg_apply_registered_env_overrides: failed to derive environment prefix for ${entrypoint_name}" >&2
		return 1
	}
	if ! declare -F gg_print_entrypoint_config_vars >/dev/null 2>&1; then
		echo "gg_apply_registered_env_overrides: config var registry helper is unavailable." >&2
		return 1
	fi

	while IFS= read -r raw_var_name; do
		gg_apply_env_override_to_config_var "${env_prefix}" "${raw_var_name}"
	done < <(gg_print_entrypoint_config_vars "${entrypoint_name}")

	while [[ $# -gt 0 ]]; do
		gg_apply_env_override_to_config_var "${env_prefix}" "$1"
		shift
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
		var_name=$(gg_normalize_config_var_name "${var_name}" || true)
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
	local requested_runtime
	local shim_path

	requested_runtime="$(gg_requested_container_runtime)" || return 1
	if [[ "${requested_runtime}" == "docker" ]]; then
		shim_path="$(gg_docker_singularity_shim_path)"
		if [[ ! -x "${shim_path}" ]]; then
			echo "Docker runtime shim not found or not executable: ${shim_path}" >&2
			return 1
		fi
		echo "${shim_path}"
		return 0
	fi
	if [[ "${requested_runtime}" == "auto" ]] && gg_auto_enable_docker_runtime_if_available; then
		shim_path="$(gg_docker_singularity_shim_path)"
		if [[ ! -x "${shim_path}" ]]; then
			echo "Docker runtime shim not found or not executable: ${shim_path}" >&2
			return 1
		fi
		echo "${shim_path}"
		return 0
	fi
	if declare -F gg_detect_site_profile >/dev/null 2>&1 \
		&& [[ "$(gg_detect_site_profile)" == "shirokane" ]] \
		&& command -v apptainer >/dev/null 2>&1; then
		echo apptainer
		return 0
	fi
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
	if [[ -n "${GG_MEM_TOOL_RESERVE_GB+x}" ]]; then
		export GG_MEM_TOOL_RESERVE_GB
	fi
	if [[ -n "${GG_MEM_TOOL_GB+x}" ]]; then
		export GG_MEM_TOOL_GB
	fi
	export NSLOTS JOB_ID SGE_TASK_ID MEM_PER_SLOT MEM_PER_HOST
}

gg_is_nonnegative_integer() {
	[[ "${1:-}" =~ ^[0-9]+$ ]]
}

gg_default_tool_memory_reserve_gb() {
	local total_gb=${1:-1}
	local reserve_gb=0
	if ! gg_is_nonnegative_integer "${total_gb}" || [[ "${total_gb}" -le 1 ]]; then
		echo 0
		return 0
	fi
	reserve_gb=$((total_gb / 8))
	if [[ "${reserve_gb}" -lt 1 ]]; then
		reserve_gb=1
	fi
	if [[ "${reserve_gb}" -gt 4 ]]; then
		reserve_gb=4
	fi
	echo "${reserve_gb}"
}

gg_normalize_memory_budget() {
	local echo_header='gg_normalize_memory_budget: '
	local max_tool_memory_gb=""

	if ! gg_is_nonnegative_integer "${GG_MEM_TOTAL_GB:-}" || [[ "${GG_MEM_TOTAL_GB:-0}" -lt 1 ]]; then
		echo "${echo_header}GG_MEM_TOTAL_GB is invalid or less than 1 (${GG_MEM_TOTAL_GB:-unset}). GG_MEM_TOTAL_GB=1"
		GG_MEM_TOTAL_GB=1
	fi

	if [[ -z "${GG_MEM_TOOL_RESERVE_GB:-}" ]]; then
		GG_MEM_TOOL_RESERVE_GB=$(gg_default_tool_memory_reserve_gb "${GG_MEM_TOTAL_GB}")
		echo "${echo_header}GG_MEM_TOOL_RESERVE_GB=${GG_MEM_TOOL_RESERVE_GB} (derived from GG_MEM_TOTAL_GB=${GG_MEM_TOTAL_GB})"
	elif ! gg_is_nonnegative_integer "${GG_MEM_TOOL_RESERVE_GB}"; then
		echo "${echo_header}GG_MEM_TOOL_RESERVE_GB is invalid (${GG_MEM_TOOL_RESERVE_GB}). GG_MEM_TOOL_RESERVE_GB=0"
		GG_MEM_TOOL_RESERVE_GB=0
	fi

	max_tool_memory_gb=$((GG_MEM_TOTAL_GB - GG_MEM_TOOL_RESERVE_GB))
	if [[ "${max_tool_memory_gb}" -lt 1 ]]; then
		max_tool_memory_gb=1
	fi

	if [[ -z "${GG_MEM_TOOL_GB:-}" ]]; then
		GG_MEM_TOOL_GB=${max_tool_memory_gb}
		echo "${echo_header}GG_MEM_TOOL_GB=${GG_MEM_TOOL_GB} (GG_MEM_TOTAL_GB - GG_MEM_TOOL_RESERVE_GB, floored at 1G)"
	elif ! gg_is_nonnegative_integer "${GG_MEM_TOOL_GB}" || [[ "${GG_MEM_TOOL_GB}" -lt 1 ]]; then
		echo "${echo_header}GG_MEM_TOOL_GB is invalid or less than 1 (${GG_MEM_TOOL_GB}). GG_MEM_TOOL_GB=1"
		GG_MEM_TOOL_GB=1
	elif [[ "${GG_MEM_TOOL_GB}" -gt "${max_tool_memory_gb}" ]]; then
		echo "${echo_header}GG_MEM_TOOL_GB=${GG_MEM_TOOL_GB} exceeds the reserved budget. GG_MEM_TOOL_GB=${max_tool_memory_gb}"
		GG_MEM_TOOL_GB=${max_tool_memory_gb}
	fi
}

gg_memory_fraction_gb() {
	local total_gb=${1:-1}
	local numerator=${2:-1}
	local denominator=${3:-1}
	local value_gb=1

	if ! gg_is_nonnegative_integer "${total_gb}" || [[ "${total_gb}" -lt 1 ]]; then
		total_gb=1
	fi
	if ! gg_is_nonnegative_integer "${numerator}"; then
		numerator=1
	fi
	if ! gg_is_nonnegative_integer "${denominator}" || [[ "${denominator}" -lt 1 ]]; then
		denominator=1
	fi

	value_gb=$((total_gb * numerator / denominator))
	if [[ "${value_gb}" -lt 1 ]]; then
		value_gb=1
	fi
	echo "${value_gb}"
}

gg_sge_memory_value_to_gb() {
	local raw_value=${1:-}
	[[ "${raw_value}" =~ ^([0-9]+([.][0-9]+)?)([bBkKmMgGtTpP]?)$ ]] || return 1
	local number="${BASH_REMATCH[1]}"
	local unit
	unit=$(printf '%s' "${BASH_REMATCH[3]}" | tr '[:lower:]' '[:upper:]')

	awk -v number="${number}" -v unit="${unit}" '
		BEGIN {
			factor = 1 / 1073741824
			if (unit == "K") factor = 1 / 1048576
			else if (unit == "M") factor = 1 / 1024
			else if (unit == "G") factor = 1
			else if (unit == "T") factor = 1024
			else if (unit == "P") factor = 1048576
			value = number * factor
			rounded = int(value)
			if (rounded < value) rounded++
			if (rounded < 1 && value > 0) rounded = 1
			print rounded
		}
	'
}

gg_sge_requested_mem_per_slot_gb() {
	local job_id=${1:-}
	local qstat_output=""
	local requested_value=""

	if [[ -z "${job_id}" ]] || ! command -v qstat >/dev/null 2>&1; then
		return 1
	fi
	qstat_output=$(qstat -j "${job_id}" 2>/dev/null || true)
	if [[ -z "${qstat_output}" ]]; then
		qstat_output=$(qstat -f -j "${job_id}" 2>/dev/null || true)
	fi
	requested_value=$(
		printf '%s\n' "${qstat_output}" \
			| awk '
				match($0, /s_vmem=[0-9]+([.][0-9]+)?[bBkKmMgGtTpP]?/) {
					print substr($0, RSTART + 7, RLENGTH - 7)
					exit
				}
			'
	)
	if [[ -z "${requested_value}" ]]; then
		return 1
	fi
	gg_sge_memory_value_to_gb "${requested_value}"
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
		if [[ "${SGE_TASK_ID:-}" =~ ^[1-9][0-9]*$ ]]; then
			echo ${echo_header}'GG_ARRAY_TASK_ID=${SGE_TASK_ID} (from legacy SGE_TASK_ID)'
			GG_ARRAY_TASK_ID=${SGE_TASK_ID}
		else
			if [[ -n "${SGE_TASK_ID:-}" ]]; then
				echo ${echo_header}'Ignoring non-numeric SGE_TASK_ID='"${SGE_TASK_ID}"'.'
			else
				echo ${echo_header}'${GG_ARRAY_TASK_ID} is undefined.'
			fi
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
		GG_MEM_PER_CPU_GB=$(gg_sge_requested_mem_per_slot_gb "${GG_JOB_ID}" || true)
		if [[ -n "${GG_MEM_PER_CPU_GB}" ]]; then
			echo ${echo_header}'GG_MEM_PER_CPU_GB='"${GG_MEM_PER_CPU_GB}"' (from AGE s_vmem)'
		fi
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
	gg_normalize_memory_budget
	gg_sync_legacy_scheduler_aliases
	echo ${echo_header}"GG_TASK_CPUS=${GG_TASK_CPUS}"
	echo ${echo_header}"GG_JOB_ID=${GG_JOB_ID}"
	echo ${echo_header}"GG_ARRAY_TASK_ID=${GG_ARRAY_TASK_ID}"
	echo ${echo_header}"GG_MEM_PER_CPU_GB=${GG_MEM_PER_CPU_GB}"
	echo ${echo_header}"GG_MEM_TOTAL_GB=${GG_MEM_TOTAL_GB}"
	echo ${echo_header}"GG_MEM_TOOL_RESERVE_GB=${GG_MEM_TOOL_RESERVE_GB}"
	echo ${echo_header}"GG_MEM_TOOL_GB=${GG_MEM_TOOL_GB}"
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
	echo "${echo_header}detected GG_TASK_CPUS=${GG_TASK_CPUS:-NA} GG_MEM_PER_CPU_GB=${GG_MEM_PER_CPU_GB:-NA} GG_MEM_TOTAL_GB=${GG_MEM_TOTAL_GB:-NA} GG_MEM_TOOL_RESERVE_GB=${GG_MEM_TOOL_RESERVE_GB:-NA} GG_MEM_TOOL_GB=${GG_MEM_TOOL_GB:-NA} GG_JOB_ID=${GG_JOB_ID:-NA} GG_ARRAY_TASK_ID=${GG_ARRAY_TASK_ID:-NA}"
	echo "${echo_header}forwarded SINGULARITYENV_GG_TASK_CPUS=${SINGULARITYENV_GG_TASK_CPUS:-unset} SINGULARITYENV_GG_MEM_PER_CPU_GB=${SINGULARITYENV_GG_MEM_PER_CPU_GB:-unset} SINGULARITYENV_GG_ARRAY_TASK_ID=${SINGULARITYENV_GG_ARRAY_TASK_ID:-unset}"
	echo "${echo_header}forwarded APPTAINERENV_GG_TASK_CPUS=${APPTAINERENV_GG_TASK_CPUS:-unset} APPTAINERENV_GG_MEM_PER_CPU_GB=${APPTAINERENV_GG_MEM_PER_CPU_GB:-unset} APPTAINERENV_GG_ARRAY_TASK_ID=${APPTAINERENV_GG_ARRAY_TASK_ID:-unset}"
	echo "${echo_header}forwarded SINGULARITYENV_OMP_NUM_THREADS=${SINGULARITYENV_OMP_NUM_THREADS:-unset} SINGULARITYENV_OPENBLAS_NUM_THREADS=${SINGULARITYENV_OPENBLAS_NUM_THREADS:-unset} SINGULARITYENV_MKL_NUM_THREADS=${SINGULARITYENV_MKL_NUM_THREADS:-unset} SINGULARITYENV_NUMEXPR_NUM_THREADS=${SINGULARITYENV_NUMEXPR_NUM_THREADS:-unset}"
	echo "${echo_header}forwarded APPTAINERENV_OMP_NUM_THREADS=${APPTAINERENV_OMP_NUM_THREADS:-unset} APPTAINERENV_OPENBLAS_NUM_THREADS=${APPTAINERENV_OPENBLAS_NUM_THREADS:-unset} APPTAINERENV_MKL_NUM_THREADS=${APPTAINERENV_MKL_NUM_THREADS:-unset} APPTAINERENV_NUMEXPR_NUM_THREADS=${APPTAINERENV_NUMEXPR_NUM_THREADS:-unset}"
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
	export SINGULARITYENV_OMP_NUM_THREADS=${GG_TASK_CPUS:-1}
	export APPTAINERENV_OMP_NUM_THREADS=${GG_TASK_CPUS:-1}
	export SINGULARITYENV_OPENBLAS_NUM_THREADS=${GG_TASK_CPUS:-1}
	export APPTAINERENV_OPENBLAS_NUM_THREADS=${GG_TASK_CPUS:-1}
	export SINGULARITYENV_MKL_NUM_THREADS=${GG_TASK_CPUS:-1}
	export APPTAINERENV_MKL_NUM_THREADS=${GG_TASK_CPUS:-1}
	export SINGULARITYENV_NUMEXPR_NUM_THREADS=${GG_TASK_CPUS:-1}
	export APPTAINERENV_NUMEXPR_NUM_THREADS=${GG_TASK_CPUS:-1}
	unset SINGULARITYENV_OMP_THREAD_LIMIT
	unset APPTAINERENV_OMP_THREAD_LIMIT
	unset SINGULARITYENV_KMP_ALL_THREADS
	unset APPTAINERENV_KMP_ALL_THREADS
	unset SINGULARITYENV_KMP_DEVICE_THREAD_LIMIT
	unset APPTAINERENV_KMP_DEVICE_THREAD_LIMIT
	unset SINGULARITYENV_KMP_TEAMS_THREAD_LIMIT
	unset APPTAINERENV_KMP_TEAMS_THREAD_LIMIT
	export SINGULARITYENV_GG_JOB_ID=${GG_JOB_ID:-1}
	export APPTAINERENV_GG_JOB_ID=${GG_JOB_ID:-1}
	export SINGULARITYENV_GG_MEM_PER_CPU_GB=${GG_MEM_PER_CPU_GB:-3}
	export APPTAINERENV_GG_MEM_PER_CPU_GB=${GG_MEM_PER_CPU_GB:-3}
	export SINGULARITYENV_GG_MEM_TOTAL_GB=${GG_MEM_TOTAL_GB:-3}
	export APPTAINERENV_GG_MEM_TOTAL_GB=${GG_MEM_TOTAL_GB:-3}
	export SINGULARITYENV_GG_MEM_TOOL_RESERVE_GB=${GG_MEM_TOOL_RESERVE_GB:-0}
	export APPTAINERENV_GG_MEM_TOOL_RESERVE_GB=${GG_MEM_TOOL_RESERVE_GB:-0}
	export SINGULARITYENV_GG_MEM_TOOL_GB=${GG_MEM_TOOL_GB:-${GG_MEM_TOTAL_GB:-3}}
	export APPTAINERENV_GG_MEM_TOOL_GB=${GG_MEM_TOOL_GB:-${GG_MEM_TOTAL_GB:-3}}
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
	export SINGULARITYENV_PYTHONPYCACHEPREFIX=/tmp/genegalleon_pycache
	export APPTAINERENV_PYTHONPYCACHEPREFIX=/tmp/genegalleon_pycache
	export SINGULARITYENV_PYTHONNOUSERSITE=1
	export APPTAINERENV_PYTHONNOUSERSITE=1
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
