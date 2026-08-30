#!/usr/bin/env bash

gg_advisory_shared_lock_acquire() {
  local lock_file=${1:-}
  if [[ -z "${lock_file}" ]]; then
    echo "Advisory shared lock path is empty." >&2
    return 1
  fi
  if [[ -n "${GG_ADVISORY_SHARED_LOCK_TOKEN:-}" ]]; then
    echo "A gene-family shared lock is already held." >&2
    return 1
  fi
  local py_exec helper_script
  py_exec=$(gg_shared_lock_python) || return 1
  helper_script=$(gg_shared_lock_helper_script) || return 1
  helper_script="$(dirname "${helper_script}")/shared_namespace_lock.py"
  GG_ADVISORY_SHARED_LOCK_TOKEN=$("${py_exec}" "${helper_script}" acquire-shared \
    "${lock_file}" --owner-pid "$$" --timeout "$(gg_lock_acquire_timeout_seconds)") || return 1
  GG_ADVISORY_SHARED_LOCK_PATH=${lock_file}
}

gg_advisory_shared_lock_release() {
  if [[ -z "${GG_ADVISORY_SHARED_LOCK_TOKEN:-}" ]]; then
    return 0
  fi
  local py_exec helper_script
  py_exec=$(gg_shared_lock_python) || return 1
  helper_script=$(gg_shared_lock_helper_script) || return 1
  helper_script="$(dirname "${helper_script}")/shared_namespace_lock.py"
  "${py_exec}" "${helper_script}" release-shared "${GG_ADVISORY_SHARED_LOCK_PATH}" \
    --token "${GG_ADVISORY_SHARED_LOCK_TOKEN}" || return 1
  GG_ADVISORY_SHARED_LOCK_TOKEN=""
  GG_ADVISORY_SHARED_LOCK_PATH=""
}

gg_lock_stale_seconds() {
  local stale_seconds="${GG_LOCK_STALE_SECONDS:-900}"
  if [[ ! "${stale_seconds}" =~ ^[0-9]+$ ]]; then
    stale_seconds=900
  fi
  if (( stale_seconds < 1 )); then
    stale_seconds=1
  fi
  echo "${stale_seconds}"
}

gg_lock_heartbeat_seconds() {
  local heartbeat_seconds="${GG_LOCK_HEARTBEAT_SECONDS:-60}"
  if [[ ! "${heartbeat_seconds}" =~ ^[0-9]+$ ]]; then
    heartbeat_seconds=60
  fi
  if (( heartbeat_seconds < 1 )); then
    heartbeat_seconds=1
  fi
  echo "${heartbeat_seconds}"
}

gg_lock_acquire_timeout_seconds() {
  local acquire_timeout_seconds="${GG_LOCK_ACQUIRE_TIMEOUT_SECONDS:-86400}"
  if [[ ! "${acquire_timeout_seconds}" =~ ^[0-9]+$ ]]; then
    acquire_timeout_seconds=86400
  fi
  if (( acquire_timeout_seconds < 1 )); then
    acquire_timeout_seconds=1
  fi
  echo "${acquire_timeout_seconds}"
}

gg_lock_poll_seconds() {
  local poll_seconds="${GG_LOCK_POLL_SECONDS:-5}"
  if [[ ! "${poll_seconds}" =~ ^[0-9]+$ ]]; then
    poll_seconds=5
  fi
  if (( poll_seconds < 1 )); then
    poll_seconds=1
  fi
  echo "${poll_seconds}"
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

gg_write_ready_marker() {
  local marker_path=$1
  mkdir -p "$(dirname "${marker_path}")"
  printf 'ready\n' > "${marker_path}"
}

gg_artifact_ready() {
  local artifact_path=$1
  if [[ -d "${artifact_path}" ]]; then
    return 0
  fi
  [[ -s "${artifact_path}" ]]
}

gg_array_expected_task_count() {
  local fallback_count=${1:-1}
  local first_task=""
  local last_task=""
  local task_step=""
  if [[ "${GG_ARRAY_TASK_COUNT:-}" =~ ^[1-9][0-9]*$ ]]; then
    echo "${GG_ARRAY_TASK_COUNT}"
    return 0
  fi
  if [[ "${SLURM_ARRAY_TASK_COUNT:-}" =~ ^[1-9][0-9]*$ ]]; then
    echo "${SLURM_ARRAY_TASK_COUNT}"
    return 0
  fi
  if [[ "${PBS_ARRAY_TASK_COUNT:-}" =~ ^[1-9][0-9]*$ ]]; then
    echo "${PBS_ARRAY_TASK_COUNT}"
    return 0
  fi
  first_task=${SGE_TASK_FIRST:-}
  last_task=${SGE_TASK_LAST:-}
  task_step=${SGE_TASK_STEPSIZE:-1}
  if [[ "${first_task}" =~ ^[1-9][0-9]*$ \
    && "${last_task}" =~ ^[1-9][0-9]*$ \
    && "${task_step}" =~ ^[1-9][0-9]*$ \
    && ${last_task} -ge ${first_task} ]]; then
    echo $(( ((last_task - first_task) / task_step) + 1 ))
    return 0
  fi
  if [[ -n "${PBS_ARRAY_INDEX:-${PBS_ARRAYID:-}}" \
    && "${fallback_count}" =~ ^[1-9][0-9]*$ ]]; then
    echo "${fallback_count}"
    return 0
  fi
  if [[ "${GG_SCHEDULER_KIND:-local}" != "local" \
    && "${fallback_count}" =~ ^[1-9][0-9]*$ ]]; then
    echo "${fallback_count}"
    return 0
  fi
  echo 1
}

gg_array_finalizer_claim() {
  local state_root=${1:-}
  local stage_name=${2:-}
  local fallback_count=${3:-1}
  local expected_count=""
  local run_id=""
  local task_id=""
  local run_dir=""
  local lock_file=""
  local ready_count=0
  if [[ -z "${state_root}" || -z "${stage_name}" ]]; then
    echo "gg_array_finalizer_claim requires STATE_ROOT and STAGE_NAME." >&2
    return 2
  fi
  if ! command -v flock >/dev/null 2>&1; then
    echo "Array finalizer coordination requires flock inside the runtime." >&2
    return 2
  fi
  expected_count=$(gg_array_expected_task_count "${fallback_count}")
  run_id=$(printf '%s' "${GG_JOB_ID:-local}" | sed 's/[^[:alnum:]._-]/_/g')
  task_id=$(printf '%s' "${GG_ARRAY_TASK_ID:-1}" | sed 's/[^[:alnum:]._-]/_/g')
  stage_name=$(printf '%s' "${stage_name}" | sed 's/[^[:alnum:]._-]/_/g')
  # Local, non-array invocations commonly reuse the synthetic job ID "1".
  # Give each process its own state directory so a completed prior invocation
  # cannot suppress the summary in a later invocation.
  if [[ "${GG_SCHEDULER_KIND:-local}" == "local" && ${expected_count} -eq 1 ]]; then
    run_id="${run_id}.$$"
  fi
  run_dir="${state_root%/}/${stage_name}/${run_id}"
  lock_file="${run_dir}.lock"
  mkdir -p "${run_dir}"
  exec {GG_ARRAY_FINALIZER_LOCK_FD}>"${lock_file}" || return 2
  if ! flock -x "${GG_ARRAY_FINALIZER_LOCK_FD}"; then
    echo "Failed to acquire the array finalizer lock: ${lock_file}. Verify that the shared workspace provides cross-node flock semantics." >&2
    exec {GG_ARRAY_FINALIZER_LOCK_FD}>&-
    GG_ARRAY_FINALIZER_LOCK_FD=""
    return 2
  fi
  printf 'ready\n' > "${run_dir}/task.${task_id}.ready"
  ready_count=$(find "${run_dir}" -maxdepth 1 -type f -name 'task.*.ready' -print | wc -l | awk '{print $1}')
  if [[ -e "${run_dir}/done" || -e "${run_dir}/claim" || ${ready_count} -lt ${expected_count} ]]; then
    flock -u "${GG_ARRAY_FINALIZER_LOCK_FD}"
    exec {GG_ARRAY_FINALIZER_LOCK_FD}>&-
    GG_ARRAY_FINALIZER_LOCK_FD=""
    echo "Array finalizer ${stage_name}: ready ${ready_count}/${expected_count}; this task will not finalize."
    return 1
  fi
  printf '%s\n' "${task_id}" > "${run_dir}/claim"
  flock -u "${GG_ARRAY_FINALIZER_LOCK_FD}"
  exec {GG_ARRAY_FINALIZER_LOCK_FD}>&-
  GG_ARRAY_FINALIZER_LOCK_FD=""
  GG_ARRAY_FINALIZER_RUN_DIR="${run_dir}"
  echo "Array finalizer ${stage_name}: ready ${ready_count}/${expected_count}; task ${task_id} claimed finalization."
}

gg_array_finalizer_complete() {
  local run_dir=${GG_ARRAY_FINALIZER_RUN_DIR:-}
  if [[ -z "${run_dir}" || ! -d "${run_dir}" ]]; then
    echo "No array finalizer claim is active." >&2
    return 1
  fi
  exec {GG_ARRAY_FINALIZER_LOCK_FD}>"${run_dir}.lock" || return 1
  flock -x "${GG_ARRAY_FINALIZER_LOCK_FD}" || return 1
  rm -f -- "${run_dir}/claim"
  printf 'done\n' > "${run_dir}/done"
  flock -u "${GG_ARRAY_FINALIZER_LOCK_FD}"
  exec {GG_ARRAY_FINALIZER_LOCK_FD}>&-
  GG_ARRAY_FINALIZER_LOCK_FD=""
  GG_ARRAY_FINALIZER_RUN_DIR=""
}

gg_array_finalizer_release() {
  local run_dir=${GG_ARRAY_FINALIZER_RUN_DIR:-}
  if [[ -z "${run_dir}" ]]; then
    return 0
  fi
  exec {GG_ARRAY_FINALIZER_LOCK_FD}>"${run_dir}.lock" || return 1
  flock -x "${GG_ARRAY_FINALIZER_LOCK_FD}" || return 1
  rm -f -- "${run_dir}/claim"
  flock -u "${GG_ARRAY_FINALIZER_LOCK_FD}"
  exec {GG_ARRAY_FINALIZER_LOCK_FD}>&-
  GG_ARRAY_FINALIZER_LOCK_FD=""
  GG_ARRAY_FINALIZER_RUN_DIR=""
}

gg_stage_transaction_lock_acquire() {
  local state_root=${1:-}
  local stage_name=${2:-}
  local safe_stage_name=""
  local lock_file=""
  if [[ -z "${state_root}" || "${state_root}" == "/" || "${state_root}" == "." \
    || "${state_root}" == ".." || -z "${stage_name}" ]]; then
    echo "gg_stage_transaction_lock_acquire requires a safe STATE_ROOT and STAGE_NAME." >&2
    return 2
  fi
  if [[ -n "${GG_STAGE_TRANSACTION_LOCK_FILE:-}" ]]; then
    echo "A stage transaction lock is already held: ${GG_STAGE_TRANSACTION_LOCK_FILE}" >&2
    return 2
  fi
  safe_stage_name=$(printf '%s' "${stage_name}" | sed 's/[^[:alnum:]._-]/_/g')
  if [[ -z "${safe_stage_name}" ]]; then
    echo "The stage transaction name is empty after normalization." >&2
    return 2
  fi
  lock_file="${state_root%/}/stage_transactions/${safe_stage_name}.lock"
  if ! gg_shared_lock_acquire "${lock_file}" "${safe_stage_name} stage transaction"; then
    return 1
  fi
  gg_shared_lock_start_heartbeat "${lock_file}"
  GG_STAGE_TRANSACTION_LOCK_FILE="${lock_file}"
  GG_STAGE_TRANSACTION_LOCK_HEARTBEAT_PID="${GG_SHARED_LOCK_HEARTBEAT_PID:-}"
  echo "Acquired stage transaction lock: ${safe_stage_name}" >&2
}

gg_stage_transaction_lock_release() {
  local lock_file=${GG_STAGE_TRANSACTION_LOCK_FILE:-}
  local heartbeat_pid=${GG_STAGE_TRANSACTION_LOCK_HEARTBEAT_PID:-}
  local release_status=0
  if [[ -z "${lock_file}" ]]; then
    return 0
  fi
  gg_shared_lock_stop_heartbeat "${heartbeat_pid}" || release_status=$?
  gg_shared_lock_release "${lock_file}" || release_status=$?
  GG_STAGE_TRANSACTION_LOCK_FILE=""
  GG_STAGE_TRANSACTION_LOCK_HEARTBEAT_PID=""
  return "${release_status}"
}

gg_shared_lock_helper_script() {
  local helper_dir="${gg_support_dir:-}"
  if [[ -z "${helper_dir}" ]]; then
    helper_dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd 2>/dev/null || true)"
  fi
  if [[ -z "${helper_dir}" || ! -s "${helper_dir}/shared_lock.py" ]]; then
    echo "shared_lock.py was not found relative to gg_util.sh" >&2
    return 1
  fi
  printf '%s\n' "${helper_dir}/shared_lock.py"
}

gg_shared_lock_python() {
  if declare -F gg_find_python_exec >/dev/null 2>&1; then
    gg_find_python_exec
    return $?
  fi

  # gg_shared_lock.sh is also sourced directly by lightweight utilities and
  # tests, without the rest of gg_util.sh. Keep lock coordination usable in
  # that mode instead of spinning forever when the shared helper cannot run.
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

gg_lock_hostname() {
  hostname 2>/dev/null || uname -n
}

gg_lock_boot_id() {
  if [[ -r /proc/sys/kernel/random/boot_id ]]; then
    tr -d '\r\n' < /proc/sys/kernel/random/boot_id
    return 0
  fi
  if command -v sysctl >/dev/null 2>&1; then
    sysctl -n kern.bootsessionuuid 2>/dev/null | tr -d '[:space:]'
    return 0
  fi
  echo ""
}

gg_shared_lock_read_metadata() {
  local lock_file=$1
  local py_exec
  local helper_script
  py_exec=$(gg_shared_lock_python) || return 1
  helper_script=$(gg_shared_lock_helper_script) || return 1
  "${py_exec}" "${helper_script}" read-metadata "${lock_file}"
}

gg_shared_lock_owner_summary() {
  local lock_file=$1
  local py_exec
  local helper_script
  py_exec=$(gg_shared_lock_python) || return 1
  helper_script=$(gg_shared_lock_helper_script) || return 1
  "${py_exec}" "${helper_script}" owner-summary "${lock_file}"
}

gg_shared_lock_try_create() {
  local lock_file=$1
  local owner_pid=${2:-$$}
  local py_exec
  local helper_script
  py_exec=$(gg_shared_lock_python) || return 1
  helper_script=$(gg_shared_lock_helper_script) || return 1
  "${py_exec}" "${helper_script}" try-create "${lock_file}" --pid "${owner_pid}"
}

gg_shared_lock_remove_if_unchanged() {
  local lock_file=$1
  local expected_device=$2
  local expected_inode=$3
  local py_exec
  local helper_script
  py_exec=$(gg_shared_lock_python) || return 1
  helper_script=$(gg_shared_lock_helper_script) || return 1
  "${py_exec}" "${helper_script}" remove-if-unchanged "${lock_file}" "${expected_device}" "${expected_inode}"
}

gg_shared_lock_reclaim_if_stale() {
  local lock_file=$1
  local description=$2
  if [[ ! -e "${lock_file}" ]]; then
    return 1
  fi
  local stale_seconds
  stale_seconds=$(gg_lock_stale_seconds)
  local py_exec
  local helper_script
  local stale_summary=""
  py_exec=$(gg_shared_lock_python) || return 1
  helper_script=$(gg_shared_lock_helper_script) || return 1
  if stale_summary=$("${py_exec}" "${helper_script}" reclaim-if-stale "${lock_file}" --stale-seconds "${stale_seconds}"); then
    echo "Recovered stale shared lock: ${description} (${stale_summary})" >&2
    return 0
  fi
  return 1
}

gg_shared_lock_start_heartbeat() {
  local lock_file=$1
  local interval_seconds
  interval_seconds=$(gg_lock_heartbeat_seconds)
  (
    local heartbeat_sleep_pid=""

    stop_heartbeat_process() {
      local child_pid
      trap '' HUP INT TERM
      # TERM may arrive after sleep was forked but before $! was assigned.
      # The shell job table already owns that child; reap it there so it cannot
      # retain the workflow's output pipes or inherited advisory-lock fd.
      while IFS= read -r child_pid; do
        kill "${child_pid}" 2>/dev/null || true
        wait "${child_pid}" 2>/dev/null || true
      done < <(jobs -pr)
      exit 0
    }

    trap 'stop_heartbeat_process' HUP INT TERM
    while [[ -e "${lock_file}" ]]; do
      sleep "${interval_seconds}" &
      heartbeat_sleep_pid=$!
      wait "${heartbeat_sleep_pid}" || exit 0
      heartbeat_sleep_pid=""
      if [[ -e "${lock_file}" ]]; then
        touch -c -- "${lock_file}" 2>/dev/null || true
      fi
    done
  ) &
  GG_SHARED_LOCK_HEARTBEAT_PID=$!
}

gg_shared_lock_stop_heartbeat() {
  local heartbeat_pid=${1:-}
  if [[ ! "${heartbeat_pid}" =~ ^[0-9]+$ ]]; then
    return 0
  fi
  kill "${heartbeat_pid}" 2>/dev/null || true
  wait "${heartbeat_pid}" 2>/dev/null || true
}

gg_shared_lock_release() {
  local lock_file=$1
  rm -f -- "${lock_file}"
}

gg_shared_lock_acquire() {
  local lock_file=$1
  local description=$2
  local task_id=${GG_ARRAY_TASK_ID:-1}
  local poll_seconds
  local timeout_seconds
  local wait_started
  local wait_logged=0
  poll_seconds=$(gg_lock_poll_seconds)
  timeout_seconds=$(gg_lock_acquire_timeout_seconds)
  wait_started=$(date +%s)
  mkdir -p "$(dirname "${lock_file}")"
  while true; do
    if gg_shared_lock_try_create "${lock_file}"; then
      return 0
    fi
    if gg_shared_lock_reclaim_if_stale "${lock_file}" "${description}"; then
      continue
    fi
    local owner_summary
    owner_summary=$(gg_shared_lock_owner_summary "${lock_file}")
    if [[ ${wait_logged} -eq 0 ]]; then
      echo "GG_ARRAY_TASK_ID=${task_id}: waiting for shared lock: ${description} (${owner_summary})" >&2
      wait_logged=1
    fi
    local now_epoch
    now_epoch=$(date +%s)
    if (( now_epoch - wait_started >= timeout_seconds )); then
      echo "GG_ARRAY_TASK_ID=${task_id}: timed out waiting for shared lock: ${description} (${owner_summary})" >&2
      return 1
    fi
    sleep "${poll_seconds}"
  done
}

gg_shared_semaphore_max_slots() {
  local max_slots="${1:-1}"
  if [[ ! "${max_slots}" =~ ^[0-9]+$ ]]; then
    max_slots=1
  fi
  if (( max_slots < 0 )); then
    max_slots=1
  fi
  echo "${max_slots}"
}

gg_shared_semaphore_acquire() {
  local semaphore_dir=$1
  local requested_slots=$2
  local description=$3
  local max_slots
  local task_id=${GG_ARRAY_TASK_ID:-1}
  local poll_seconds
  local timeout_seconds
  local wait_started
  local wait_logged=0
  local slot_idx=0
  local slot_lock=""

  max_slots=$(gg_shared_semaphore_max_slots "${requested_slots}")
  if (( max_slots < 1 )); then
    GG_SHARED_SEMAPHORE_SLOT_LOCK_FILE=""
    GG_SHARED_SEMAPHORE_SLOT_INDEX=""
    GG_SHARED_SEMAPHORE_MAX_SLOTS="${max_slots}"
    return 0
  fi

  poll_seconds=$(gg_lock_poll_seconds)
  timeout_seconds=$(gg_lock_acquire_timeout_seconds)
  wait_started=$(date +%s)
  mkdir -p "${semaphore_dir}"
  while true; do
    for (( slot_idx=1; slot_idx<=max_slots; slot_idx++ )); do
      slot_lock="${semaphore_dir}/slot.${slot_idx}.lock"
      if gg_shared_lock_try_create "${slot_lock}"; then
        GG_SHARED_SEMAPHORE_SLOT_LOCK_FILE="${slot_lock}"
        GG_SHARED_SEMAPHORE_SLOT_INDEX="${slot_idx}"
        GG_SHARED_SEMAPHORE_MAX_SLOTS="${max_slots}"
        return 0
      fi
      if gg_shared_lock_reclaim_if_stale "${slot_lock}" "${description} slot ${slot_idx}/${max_slots}"; then
        if gg_shared_lock_try_create "${slot_lock}"; then
          GG_SHARED_SEMAPHORE_SLOT_LOCK_FILE="${slot_lock}"
          GG_SHARED_SEMAPHORE_SLOT_INDEX="${slot_idx}"
          GG_SHARED_SEMAPHORE_MAX_SLOTS="${max_slots}"
          return 0
        fi
      fi
    done
    if [[ ${wait_logged} -eq 0 ]]; then
      echo "GG_ARRAY_TASK_ID=${task_id}: waiting for shared semaphore slot: ${description} (max_concurrent=${max_slots}; lock_dir=${semaphore_dir})" >&2
      wait_logged=1
    fi
    local now_epoch
    now_epoch=$(date +%s)
    if (( now_epoch - wait_started >= timeout_seconds )); then
      echo "GG_ARRAY_TASK_ID=${task_id}: timed out waiting for shared semaphore slot: ${description} (max_concurrent=${max_slots}; lock_dir=${semaphore_dir})" >&2
      return 1
    fi
    sleep "${poll_seconds}"
  done
}

gg_shared_semaphore_release() {
  local slot_lock="${1:-${GG_SHARED_SEMAPHORE_SLOT_LOCK_FILE:-}}"
  if [[ -z "${slot_lock}" ]]; then
    return 0
  fi
  gg_shared_lock_release "${slot_lock}"
  if [[ "${slot_lock}" == "${GG_SHARED_SEMAPHORE_SLOT_LOCK_FILE:-}" ]]; then
    GG_SHARED_SEMAPHORE_SLOT_LOCK_FILE=""
    GG_SHARED_SEMAPHORE_SLOT_INDEX=""
    GG_SHARED_SEMAPHORE_MAX_SLOTS=""
  fi
}

gg_run_with_shared_semaphore() {
  local semaphore_dir=$1
  local requested_slots=$2
  local description=$3
  shift 3
  local task_id=${GG_ARRAY_TASK_ID:-1}
  local max_slots
  local slot_lock=""
  local heartbeat_pid=""
  local command_exit_code=0
  local saved_exit_trap=""
  local saved_hup_trap=""
  local saved_int_trap=""
  local saved_term_trap=""
  local initial_slot_lock="${GG_SHARED_SEMAPHORE_SLOT_LOCK_FILE:-}"

  max_slots=$(gg_shared_semaphore_max_slots "${requested_slots}")
  if (( max_slots < 1 )); then
    if "$@"; then
      return 0
    fi
    return $?
  fi

  cleanup_shared_semaphore() {
    local acquired_slot_lock="${slot_lock}"
    if [[ -z "${acquired_slot_lock}" \
      && "${GG_SHARED_SEMAPHORE_SLOT_LOCK_FILE:-}" != "${initial_slot_lock}" ]]; then
      acquired_slot_lock="${GG_SHARED_SEMAPHORE_SLOT_LOCK_FILE:-}"
    fi
    gg_shared_lock_stop_heartbeat "${heartbeat_pid}"
    gg_shared_semaphore_release "${acquired_slot_lock}"
  }

  restore_shared_semaphore_traps() {
    trap - EXIT
    trap - HUP
    trap - INT
    trap - TERM
    if [[ -n "${saved_exit_trap}" ]]; then
      eval "${saved_exit_trap}"
    fi
    if [[ -n "${saved_hup_trap}" ]]; then
      eval "${saved_hup_trap}"
    fi
    if [[ -n "${saved_int_trap}" ]]; then
      eval "${saved_int_trap}"
    fi
    if [[ -n "${saved_term_trap}" ]]; then
      eval "${saved_term_trap}"
    fi
  }

  shared_semaphore_signal_handler() {
    local signal_name=$1
    cleanup_shared_semaphore
    restore_shared_semaphore_traps
    kill "-${signal_name}" "${BASHPID:-$$}"
  }

  saved_exit_trap=$(trap -p EXIT || true)
  saved_hup_trap=$(trap -p HUP || true)
  saved_int_trap=$(trap -p INT || true)
  saved_term_trap=$(trap -p TERM || true)
  trap 'cleanup_shared_semaphore; restore_shared_semaphore_traps' EXIT
  trap 'shared_semaphore_signal_handler HUP' HUP
  trap 'shared_semaphore_signal_handler INT' INT
  trap 'shared_semaphore_signal_handler TERM' TERM

  if ! gg_shared_semaphore_acquire "${semaphore_dir}" "${max_slots}" "${description}"; then
    restore_shared_semaphore_traps
    return 1
  fi
  slot_lock="${GG_SHARED_SEMAPHORE_SLOT_LOCK_FILE:-}"
  if [[ -z "${slot_lock}" ]]; then
    echo "Failed to resolve shared semaphore slot after acquire: ${description}" >&2
    cleanup_shared_semaphore
    restore_shared_semaphore_traps
    return 1
  fi

  gg_shared_lock_start_heartbeat "${slot_lock}"
  heartbeat_pid=${GG_SHARED_LOCK_HEARTBEAT_PID:-}
  echo "GG_ARRAY_TASK_ID=${task_id}: acquired shared semaphore slot ${GG_SHARED_SEMAPHORE_SLOT_INDEX:-?}/${GG_SHARED_SEMAPHORE_MAX_SLOTS:-${max_slots}}: ${description}" >&2
  if "$@"; then
    command_exit_code=0
  else
    command_exit_code=$?
  fi
  cleanup_shared_semaphore
  restore_shared_semaphore_traps
  return "${command_exit_code}"
}

gg_lock_pid_is_alive() {
  local pid=$1
  if [[ ! "${pid}" =~ ^[0-9]+$ ]]; then
    return 1
  fi
  kill -0 "${pid}" 2>/dev/null
}

gg_array_download_once() {
  local lock_file=$1
  local artifact_path=$2
  local description=$3
  shift 3
  local task_id=${GG_ARRAY_TASK_ID:-1}
  local heartbeat_pid=""
  local artifact_exit_code=0

  mkdir -p "$(dirname "${lock_file}")"

  if gg_artifact_ready "${artifact_path}"; then
    return 0
  fi
  if ! gg_shared_lock_acquire "${lock_file}" "${description}"; then
    return 1
  fi
  gg_shared_lock_start_heartbeat "${lock_file}"
  heartbeat_pid=${GG_SHARED_LOCK_HEARTBEAT_PID:-}

  if gg_artifact_ready "${artifact_path}"; then
    gg_shared_lock_stop_heartbeat "${heartbeat_pid}"
    gg_shared_lock_release "${lock_file}"
    return 0
  fi

  echo "GG_ARRAY_TASK_ID=${task_id}: starting shared artifact preparation: ${description}" >&2
  # Artifact builders may run tools that log to stdout; keep caller-captured
  # return values, such as DB prefixes, clean.
  "$@" >&2
  artifact_exit_code=$?
  if [[ ${artifact_exit_code} -ne 0 ]]; then
    echo "GG_ARRAY_TASK_ID=${task_id}: shared artifact preparation failed: ${description}" >&2
    gg_shared_lock_stop_heartbeat "${heartbeat_pid}"
    gg_shared_lock_release "${lock_file}"
    return "${artifact_exit_code}"
  fi
  if ! gg_artifact_ready "${artifact_path}"; then
    echo "Shared artifact was not ready after synchronization: ${artifact_path}" >&2
    gg_shared_lock_stop_heartbeat "${heartbeat_pid}"
    gg_shared_lock_release "${lock_file}"
    return 1
  fi

  echo "GG_ARRAY_TASK_ID=${task_id}: shared artifact ready: ${description}" >&2
  gg_shared_lock_stop_heartbeat "${heartbeat_pid}"
  gg_shared_lock_release "${lock_file}"
}
