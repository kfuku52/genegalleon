#!/usr/bin/env bash
set -euo pipefail

# Run every *_entrypoint.sh in dependency-aware debug order.
#
# Required chain requested by project workflow:
#   gg_genome_evolution_entrypoint.sh
#     -> gg_gene_evolution_entrypoint.sh (mode_orthogroup=1, mode_query2family=0)
#     -> gg_gene_evolution_entrypoint.sh (mode_orthogroup=0, mode_query2family=1)
#     -> gg_gene_database_entrypoint.sh
#     -> gg_gene_convergence_entrypoint.sh
#
# Other entrypoints are executed around this chain so all entrypoints are covered.

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/.." && pwd)"

wrapper_bin="${GG_WRAPPER_BIN:-/tmp/gg_wrapper_bin}"
if [[ -x "${wrapper_bin}/singularity" ]]; then
  run_path="${wrapper_bin}:${PATH}"
else
  run_path="${PATH}"
fi

logs_dir="${GG_DEBUG_LOG_DIR:-${repo_root}/workspace/output/debug_entrypoint_logs}"
timeout_sec="${GG_ENTRYPOINT_TIMEOUT_SEC:-1800}"   # default 30 min per entrypoint; set 0 to disable
dry_run="${GG_ENTRYPOINT_DRY_RUN:-0}"          # 1 = print commands only
only_steps_csv="${GG_ENTRYPOINT_ONLY_STEPS:-}" # optional comma-separated step IDs
benchmark_enabled="${GG_BENCHMARK:-1}"         # 1 = collect runtime metrics in summary.tsv
benchmark_raw="${GG_BENCHMARK_RAW:-1}"         # 1 = keep raw /usr/bin/time outputs
run_id="${GG_RUN_ID:-$(date -u +"%Y%m%dT%H%M%SZ")}"

nslots="${NSLOTS:-2}"
job_id="${JOB_ID:-1}"
sge_task_id="${SGE_TASK_ID:-1}"
mem_per_slot="${MEM_PER_SLOT:-4096}"
gg_wrapper_image="${GG_WRAPPER_IMAGE:-local/genegalleon:dev}"
pymol_headless="${PYMOL_HEADLESS:-1}"
qt_qpa_platform="${QT_QPA_PLATFORM:-offscreen}"

mkdir -p "${logs_dir}"
summary_tsv="${logs_dir}/summary.tsv"
# columns:
# step_id, entrypoint, status, exit_code, log_file,
# wall_sec, user_sec, sys_sec, max_rss_kb, ncpu, timeout_sec, start_iso, end_iso, run_id
: > "${summary_tsv}"
time_raw_dir="${logs_dir}/time_raw"
if [[ "${benchmark_raw}" -eq 1 ]]; then
  mkdir -p "${time_raw_dir}"
fi

append_summary_row() {
  local step_id="$1"
  local ep_name="$2"
  local result_status="$3"
  local rc="$4"
  local log_file="$5"
  local wall_sec="$6"
  local user_sec="$7"
  local sys_sec="$8"
  local max_rss_kb="$9"
  local start_iso="${10}"
  local end_iso="${11}"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${step_id}" \
    "${ep_name}" \
    "${result_status}" \
    "${rc}" \
    "${log_file}" \
    "${wall_sec}" \
    "${user_sec}" \
    "${sys_sec}" \
    "${max_rss_kb}" \
    "${nslots}" \
    "${timeout_sec}" \
    "${start_iso}" \
    "${end_iso}" \
    "${run_id}" >> "${summary_tsv}"
}

step_passed() {
  local step_id="$1"
  awk -F '\t' -v step="${step_id}" '
    $1 == step { status = $3 }
    END { exit(status == "PASS" ? 0 : 1) }
  ' "${summary_tsv}"
}

prereq_fail_reason=""
check_step_prerequisites() {
  local step_id="$1"
  local file_orthogroup_genecount_selected="${repo_root}/workspace/output/orthofinder/Orthogroups_filtered/Orthogroups.GeneCount.selected.tsv"
  local dir_orthogroup="${repo_root}/workspace/output/orthogroup"
  local file_orthogroup_db="${dir_orthogroup}/gg_orthogroup.db"

  prereq_fail_reason=""
  case "${step_id}" in
    gg_gene_evolution_mode_orthogroup)
      if [[ ! -s "${file_orthogroup_genecount_selected}" ]]; then
        prereq_fail_reason="missing required file: ${file_orthogroup_genecount_selected}"
        return 1
      fi
      ;;
    gg_gene_database)
      if [[ ! -d "${dir_orthogroup}/stat_tree" || ! -d "${dir_orthogroup}/stat_branch" ]]; then
        prereq_fail_reason="missing required directories: ${dir_orthogroup}/stat_tree or ${dir_orthogroup}/stat_branch"
        return 1
      fi
      ;;
    gg_gene_convergence)
      if [[ ! -s "${file_orthogroup_db}" ]]; then
        prereq_fail_reason="missing required database: ${file_orthogroup_db}"
        return 1
      fi
      ;;
  esac
  return 0
}

detect_time_mode() {
  local probe_file="${logs_dir}/.time_probe.$$"
  rm -f -- "${probe_file}"
  if /usr/bin/time -l -p /usr/bin/true >/dev/null 2> "${probe_file}"; then
    rm -f -- "${probe_file}"
    echo "bsd_full"
    return 0
  fi
  if /usr/bin/time -p /usr/bin/true >/dev/null 2> "${probe_file}"; then
    rm -f -- "${probe_file}"
    echo "bsd_basic"
    return 0
  fi
  if /usr/bin/time -v /usr/bin/true >/dev/null 2> "${probe_file}"; then
    rm -f -- "${probe_file}"
    echo "gnu"
    return 0
  fi
  rm -f -- "${probe_file}"
  echo "none"
}

extract_time_metrics() {
  local time_file="$1"
  local mode="$2"
  local user_sec="NA"
  local sys_sec="NA"
  local max_rss_kb="NA"
  local value=""

  if [[ ! -s "${time_file}" ]]; then
    printf '%s\t%s\t%s\n' "${user_sec}" "${sys_sec}" "${max_rss_kb}"
    return 0
  fi

  if [[ "${mode}" == "bsd_full" || "${mode}" == "bsd_basic" ]]; then
    value=$(awk '$1=="user"{print $2; exit}' "${time_file}" || true)
    if [[ -n "${value}" ]]; then user_sec="${value}"; fi
    value=$(awk '$1=="sys"{print $2; exit}' "${time_file}" || true)
    if [[ -n "${value}" ]]; then sys_sec="${value}"; fi
    if [[ "${mode}" == "bsd_full" ]]; then
      value=$(awk '/maximum resident set size/{print $1; exit}' "${time_file}" || true)
      if [[ -n "${value}" ]]; then
        if [[ "${value}" =~ ^[0-9]+$ ]]; then
          # BSD time may report bytes on some platforms; normalize to kbytes for cross-platform comparability.
          if [[ "${value}" -gt 1048576 ]]; then
            max_rss_kb=$(( (value + 1023) / 1024 ))
          else
            max_rss_kb="${value}"
          fi
        else
          max_rss_kb="${value}"
        fi
      fi
    fi
  elif [[ "${mode}" == "gnu" ]]; then
    value=$(awk -F': *' '/^User time \(seconds\):/{print $2; exit}' "${time_file}" || true)
    if [[ -n "${value}" ]]; then user_sec="${value}"; fi
    value=$(awk -F': *' '/^System time \(seconds\):/{print $2; exit}' "${time_file}" || true)
    if [[ -n "${value}" ]]; then sys_sec="${value}"; fi
    value=$(awk -F': *' '/^Maximum resident set size \(kbytes\):/{print $2; exit}' "${time_file}" || true)
    if [[ -n "${value}" ]]; then max_rss_kb="${value}"; fi
  fi

  printf '%s\t%s\t%s\n' "${user_sec}" "${sys_sec}" "${max_rss_kb}"
}

time_mode="none"
if [[ "${benchmark_enabled}" -eq 1 ]]; then
  time_mode="$(detect_time_mode)"
fi

# format: step_id|entrypoint_script|extra_envs
ordered_steps=(
  "gg_input_generation|gg_input_generation_entrypoint.sh|"
  "gg_transcriptome_generation|gg_transcriptome_generation_entrypoint.sh|"
  "gg_genome_annotation|gg_genome_annotation_entrypoint.sh|"
  "gg_genome_evolution|gg_genome_evolution_entrypoint.sh|"
  "gg_gene_evolution_mode_orthogroup|gg_gene_evolution_entrypoint.sh|mode_orthogroup=1 mode_query2family=0 run_hyphy_relax=0 run_hyphy_relax_reversed=0"
  "gg_gene_evolution_mode_query2family|gg_gene_evolution_entrypoint.sh|mode_orthogroup=0 mode_query2family=1 run_hyphy_relax=0 run_hyphy_relax_reversed=0"
  "gg_gene_database|gg_gene_database_entrypoint.sh|"
  "gg_gene_convergence|gg_gene_convergence_entrypoint.sh|"
  "gg_progress_summary|gg_progress_summary_entrypoint.sh|"
)

step_is_selected() {
  local step_id="$1"
  if [[ -z "${only_steps_csv}" ]]; then
    return 0
  fi
  local token
  IFS=',' read -r -a token <<< "${only_steps_csv}"
  local s
  for s in "${token[@]}"; do
    if [[ "${step_id}" == "${s}" ]]; then
      return 0
    fi
  done
  return 1
}

cleanup_docker_containers() {
  local ids
  ids="$(docker ps -q --filter "ancestor=${gg_wrapper_image}" 2>/dev/null || true)"
  if [[ -n "${ids}" ]]; then
    docker kill ${ids} >/dev/null 2>&1 || true
  fi
}

run_with_timeout() {
  local limit="$1"
  shift

  if [[ "${limit}" -le 0 ]]; then
    "$@"
    return $?
  fi

  local marker
  marker="$(mktemp "${TMPDIR:-/tmp}/gg_timeout.XXXXXX")"
  rm -f "${marker}"

  (
    "$@"
  ) &
  local pid=$!

  (
    sleep "${limit}"
    echo timeout > "${marker}"
    kill -TERM "-${pid}" 2>/dev/null || kill -TERM "${pid}" 2>/dev/null || true
    sleep 5
    kill -KILL "-${pid}" 2>/dev/null || kill -KILL "${pid}" 2>/dev/null || true
  ) &
  local watchdog=$!

  wait "${pid}"
  local rc=$?
  kill "${watchdog}" 2>/dev/null || true
  wait "${watchdog}" 2>/dev/null || true

  if [[ -s "${marker}" ]]; then
    rm -f "${marker}"
    return 124
  fi
  rm -f "${marker}"
  return "${rc}"
}

run_one_step() {
  local step_id="$1"
  local ep_name="$2"
  local extra_env_str="${3:-}"
  local log_file="${logs_dir}/${step_id}.log"
  local time_file=""
  local rc=0
  local result_status="FAIL"
  local wall_sec="NA"
  local user_sec="NA"
  local sys_sec="NA"
  local max_rss_kb="NA"
  local start_iso
  local end_iso
  local start_epoch
  local end_epoch
  local metrics
  local -a time_cmd=()
  local -a extra_env=()

  local -a ep_env=(
    "PATH=${run_path}"
    "GG_WRAPPER_IMAGE=${gg_wrapper_image}"
    "gg_debug_mode=1"
    "NSLOTS=${nslots}"
    "JOB_ID=${job_id}"
    "SGE_TASK_ID=${sge_task_id}"
    "MEM_PER_SLOT=${mem_per_slot}"
    "PYMOL_HEADLESS=${pymol_headless}"
    "QT_QPA_PLATFORM=${qt_qpa_platform}"
  )

  if [[ -n "${extra_env_str}" ]]; then
    # shellcheck disable=SC2206
    extra_env=( ${extra_env_str} )
    ep_env+=( "${extra_env[@]}" )
  fi

  echo "=== RUN ${step_id} (${ep_name}) ==="
  start_iso="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  start_epoch="$(date +%s)"

  if ! check_step_prerequisites "${step_id}"; then
    end_iso="${start_iso}"
    wall_sec=0
    result_status="SKIP_DEP"
    rc=125
    printf '%s\n' "${prereq_fail_reason}" > "${log_file}"
    append_summary_row \
      "${step_id}" \
      "${ep_name}" \
      "${result_status}" \
      "${rc}" \
      "${log_file}" \
      "${wall_sec}" \
      "${user_sec}" \
      "${sys_sec}" \
      "${max_rss_kb}" \
      "${start_iso}" \
      "${end_iso}"
    echo "RESULT ${step_id} ${result_status} rc=${rc} (${prereq_fail_reason})"
    cleanup_docker_containers
    return 0
  fi

  if [[ "${dry_run}" -eq 1 ]]; then
    echo "DRY-RUN: env ${ep_env[*]} bash ${script_dir}/${ep_name}"
    result_status="DRYRUN"
    rc=0
    wall_sec=0
  else
    if [[ "${benchmark_raw}" -eq 1 ]]; then
      time_file="${time_raw_dir}/${step_id}.time.txt"
    else
      time_file="${logs_dir}/.${step_id}.time.txt"
    fi
    rm -f -- "${time_file}"

    time_cmd=()
    if [[ "${benchmark_enabled}" -eq 1 && "${time_mode}" != "none" ]]; then
      time_cmd=( "/usr/bin/time" )
      if [[ "${time_mode}" == "bsd_full" ]]; then
        time_cmd+=( "-l" "-p" )
      elif [[ "${time_mode}" == "bsd_basic" ]]; then
        time_cmd+=( "-p" )
      elif [[ "${time_mode}" == "gnu" ]]; then
        time_cmd+=( "-v" )
      fi
      time_cmd+=( "-o" "${time_file}" )
    fi

    set +e
    if [[ ${#time_cmd[@]} -gt 0 ]]; then
      run_with_timeout "${timeout_sec}" "${time_cmd[@]}" env "${ep_env[@]}" bash "${script_dir}/${ep_name}" > "${log_file}" 2>&1
    else
      run_with_timeout "${timeout_sec}" env "${ep_env[@]}" bash "${script_dir}/${ep_name}" > "${log_file}" 2>&1
    fi
    rc=$?
    set -e

    if [[ "${rc}" -eq 0 ]]; then
      result_status="PASS"
    elif [[ "${rc}" -eq 124 ]]; then
      result_status="TIMEOUT"
    fi

    if [[ "${benchmark_enabled}" -eq 1 && "${time_mode}" != "none" ]]; then
      metrics="$(extract_time_metrics "${time_file}" "${time_mode}")"
      IFS=$'\t' read -r user_sec sys_sec max_rss_kb <<< "${metrics}"
    fi

    if [[ "${benchmark_raw}" -ne 1 ]]; then
      rm -f -- "${time_file}"
    fi
  fi
  end_iso="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  end_epoch="$(date +%s)"
  if [[ "${wall_sec}" == "NA" ]]; then
    wall_sec=$((end_epoch - start_epoch))
  fi

  append_summary_row \
    "${step_id}" \
    "${ep_name}" \
    "${result_status}" \
    "${rc}" \
    "${log_file}" \
    "${wall_sec}" \
    "${user_sec}" \
    "${sys_sec}" \
    "${max_rss_kb}" \
    "${start_iso}" \
    "${end_iso}"
  echo "RESULT ${step_id} ${result_status} rc=${rc}"
  cleanup_docker_containers
}

echo "Dependency-aware debug order:"
for i in "${!ordered_steps[@]}"; do
  IFS='|' read -r step_id ep_name extra_env_str <<< "${ordered_steps[$i]}"
  n=$((i + 1))
  if [[ -n "${extra_env_str}" ]]; then
    echo "  ${n}. ${step_id}: ${ep_name} (${extra_env_str})"
  else
    echo "  ${n}. ${step_id}: ${ep_name}"
  fi
done
if [[ -n "${only_steps_csv}" ]]; then
  echo "Selected steps only: ${only_steps_csv}"
fi
echo ""

for step in "${ordered_steps[@]}"; do
  IFS='|' read -r step_id ep_name extra_env_str <<< "${step}"
  if ! step_is_selected "${step_id}"; then
    continue
  fi
  run_one_step "${step_id}" "${ep_name}" "${extra_env_str}"
done

echo ""
echo "=== SUMMARY ==="
cat "${summary_tsv}"

if awk -F '\t' '($3!="PASS" && $3!="DRYRUN"){exit 1}' "${summary_tsv}"; then
  exit 0
else
  exit 1
fi
