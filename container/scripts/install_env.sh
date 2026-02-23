#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <env_name> <required_pkg_file> <optional_pkg_file> [bootstrap_pkg]"
  exit 1
fi

env_name=$1
required_file=$2
optional_file=$3
bootstrap_pkg=${4:-python=3.10}

max_attempts=${MICROMAMBA_MAX_ATTEMPTS:-10}
netcheck_attempts=${MICROMAMBA_NETCHECK_ATTEMPTS:-3}
use_defaults_channel=${USE_DEFAULTS_CHANNEL:-1}
strict_channel_priority=${STRICT_CHANNEL_PRIORITY:-1}
remote_connect_timeout_secs=${MICROMAMBA_REMOTE_CONNECT_TIMEOUT_SECS:-30}
remote_max_retries=${MICROMAMBA_REMOTE_MAX_RETRIES:-10}
remote_backoff_factor=${MICROMAMBA_REMOTE_BACKOFF_FACTOR:-2}
repodata_use_zst=${MICROMAMBA_REPODATA_USE_ZST:-0}
repodata_ttl=${MICROMAMBA_REPODATA_TTL:-3600}
retry_clean_mode=${MICROMAMBA_RETRY_CLEAN_MODE:-locks}
mamba_no_low_speed_limit=${MAMBA_NO_LOW_SPEED_LIMIT:-1}
retry_sleep_base_sec=${MICROMAMBA_RETRY_SLEEP_BASE_SEC:-20}
retry_sleep_max_sec=${MICROMAMBA_RETRY_SLEEP_MAX_SEC:-120}

channel_args=(-c conda-forge -c bioconda)
if [[ "${use_defaults_channel}" == "1" ]]; then
  channel_args+=(-c defaults)
fi

priority_args=()
if [[ "${strict_channel_priority}" == "1" ]]; then
  priority_args+=(--strict-channel-priority)
fi

run_micromamba_with_retry() {
  local label=$1
  shift
  local attempt rc sleep_secs tmp_log
  for ((attempt=1; attempt<=max_attempts; attempt++)); do
    echo "[install_env] ${label} (attempt ${attempt}/${max_attempts})"
    tmp_log=$(mktemp)
    set +e
    "$@" 2>&1 | tee "${tmp_log}"
    rc=${PIPESTATUS[0]}
    set -e
    rm -f "${tmp_log}"
    if (( rc == 0 )); then
      return 0
    fi
    if (( attempt == max_attempts )); then
      echo "[install_env] ${label} failed after ${max_attempts} attempt(s) (exit ${rc})"
      return "${rc}"
    fi
    echo "[install_env] ${label} failed (exit ${rc}); retry cleanup mode=${retry_clean_mode}"
    case "${retry_clean_mode}" in
      full)
        micromamba clean -a -y || true
        ;;
      index)
        micromamba clean -i -l -y || true
        ;;
      locks|*)
        micromamba clean -l -y || true
        ;;
    esac
    sleep_secs=$(( attempt * retry_sleep_base_sec ))
    if (( sleep_secs > retry_sleep_max_sec )); then
      sleep_secs=${retry_sleep_max_sec}
    fi
    echo "[install_env] waiting ${sleep_secs}s before retry"
    sleep "${sleep_secs}"
  done
}

set_micromamba_config_safe() {
  local key=$1
  local value=$2
  local out
  out=$(micromamba config set "${key}" "${value}" 2>&1 || true)
  if [[ "${out}" == *"Key is invalid"* ]]; then
    echo "[install_env] WARNING: micromamba config key unsupported: ${key}"
    return 1
  fi
  if [[ -n "${out}" ]]; then
    echo "[install_env] WARNING: micromamba config ${key} returned: ${out}"
    return 1
  fi
  if micromamba config list | awk -v key="${key}:" '$1 == key {found=1} END {exit(found ? 0 : 1)}'; then
    echo "[install_env] micromamba config ${key}=${value}"
  else
    echo "[install_env] WARNING: micromamba config ${key} did not persist"
    return 1
  fi
}

configure_micromamba_network() {
  local zst_value=true
  case "${repodata_use_zst,,}" in
    0|false|no|off)
      zst_value=false
      ;;
  esac
  export MAMBA_NO_LOW_SPEED_LIMIT="${mamba_no_low_speed_limit}"
  echo "[install_env] env MAMBA_NO_LOW_SPEED_LIMIT=${MAMBA_NO_LOW_SPEED_LIMIT}"
  set_micromamba_config_safe remote_connect_timeout_secs "${remote_connect_timeout_secs}"
  set_micromamba_config_safe remote_max_retries "${remote_max_retries}"
  set_micromamba_config_safe remote_backoff_factor "${remote_backoff_factor}"
  set_micromamba_config_safe local_repodata_ttl "${repodata_ttl}"
  set_micromamba_config_safe repodata_use_zst "${zst_value}"
}

preflight_conda_endpoints() {
  local url attempt
  local urls=(
    "https://conda.anaconda.org/conda-forge/noarch/repodata.json"
    "https://conda.anaconda.org/bioconda/noarch/repodata.json"
  )
  if [[ "${use_defaults_channel}" == "1" ]]; then
    urls+=("https://repo.anaconda.com/pkgs/main/noarch/repodata.json")
  fi

  for url in "${urls[@]}"; do
    for ((attempt=1; attempt<=netcheck_attempts; attempt++)); do
      if curl -fsSIL --retry 3 --retry-delay 2 --retry-max-time 45 --connect-timeout 10 --max-time 30 "${url}" >/dev/null; then
        break
      fi
      if (( attempt == netcheck_attempts )); then
        echo "[install_env] WARNING: network preflight failed for ${url}; continuing with retries enabled"
      else
        sleep $(( attempt * 5 ))
      fi
    done
  done
}

read_pkg_file() {
  local f=$1
  if [[ -f "${f}" ]]; then
    sed -e 's/[[:space:]]*#.*$//' -e '/^[[:space:]]*$/d' "${f}"
  fi
}

mapfile -t required_pkgs < <(read_pkg_file "${required_file}")
mapfile -t optional_pkgs < <(read_pkg_file "${optional_file}")

create_specs=()
if [[ -n "${bootstrap_pkg}" ]]; then
  create_specs+=("${bootstrap_pkg}")
fi
if [[ ${#required_pkgs[@]} -gt 0 ]]; then
  create_specs+=("${required_pkgs[@]}")
fi

configure_micromamba_network
preflight_conda_endpoints

env_prefix() {
  if [[ "${env_name}" == "base" ]]; then
    echo "${MAMBA_ROOT_PREFIX:-/opt/conda}"
  else
    echo "${MAMBA_ROOT_PREFIX:-/opt/conda}/envs/${env_name}"
  fi
}

target_prefix=$(env_prefix)
if [[ "${env_name}" == "base" || -f "${target_prefix}/conda-meta/history" ]]; then
  echo "[install_env] Updating env '${env_name}' with ${#required_pkgs[@]} required package(s)"
  if [[ ${#create_specs[@]} -gt 0 ]]; then
    install_cmd=(micromamba install -y --retry-clean-cache --repodata-ttl "${repodata_ttl}" "${priority_args[@]}" -n "${env_name}" "${channel_args[@]}" "${create_specs[@]}")
    run_micromamba_with_retry "install required specs into '${env_name}'" "${install_cmd[@]}"
  else
    echo "[install_env] No required package specs provided for '${env_name}'"
  fi
else
  echo "[install_env] Creating env '${env_name}' with ${#required_pkgs[@]} required package(s)"
  create_cmd=(micromamba create -y --retry-clean-cache --repodata-ttl "${repodata_ttl}" "${priority_args[@]}" -n "${env_name}" "${channel_args[@]}")
  if [[ ${#create_specs[@]} -gt 0 ]]; then
    create_cmd+=("${create_specs[@]}")
  fi
  run_micromamba_with_retry "create '${env_name}'" "${create_cmd[@]}"
fi

failed_optional="/opt/pg/logs/failed_optional_${env_name}.txt"
mkdir -p "$(dirname "${failed_optional}")"
: > "${failed_optional}"

if [[ ${#optional_pkgs[@]} -gt 0 ]]; then
  echo "[install_env] Installing ${#optional_pkgs[@]} optional package(s) for '${env_name}' in a single transaction"
  if ! run_micromamba_with_retry \
      "install optional batch for '${env_name}'" \
      micromamba install -y --retry-clean-cache --repodata-ttl "${repodata_ttl}" "${priority_args[@]}" -n "${env_name}" "${channel_args[@]}" "${optional_pkgs[@]}"; then
    echo "[install_env] Optional batch install failed for '${env_name}', retrying per package"
    for pkg in "${optional_pkgs[@]}"; do
      echo "  [optional] ${pkg}"
      if ! run_micromamba_with_retry \
          "install optional package '${pkg}' for '${env_name}'" \
          micromamba install -y --retry-clean-cache --repodata-ttl "${repodata_ttl}" "${priority_args[@]}" -n "${env_name}" "${channel_args[@]}" "${pkg}"; then
        echo "${pkg}" >> "${failed_optional}"
      fi
    done
  fi
else
  echo "[install_env] No optional packages listed for '${env_name}'"
fi

echo "[install_env] Completed '${env_name}'"
if [[ -s "${failed_optional}" ]]; then
  echo "[install_env] Optional package failures for '${env_name}':"
  cat "${failed_optional}"
fi
