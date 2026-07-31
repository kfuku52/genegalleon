#!/usr/bin/env bash

# SHIROKANE AGE resources for validating and installing a GeneGalleon SIF.
#$ -S /bin/bash
#$ -cwd
#$ -pe def_slot 1
#$ -l s_vmem=8G
#$ -l ljob

set -euo pipefail

submit_dir="${SGE_O_WORKDIR:-${PWD}}"
repo_root=""
for candidate_root in "${submit_dir}" "${submit_dir}/.."; do
  if [[ -s "${candidate_root}/workflow/support/gg_site_runtime.sh" ]]; then
    repo_root="$(cd "${candidate_root}" && pwd -P)"
    break
  fi
done
if [[ -z "${repo_root}" ]]; then
  echo "Could not locate the GeneGalleon repository from SGE_O_WORKDIR=${submit_dir}." >&2
  exit 1
fi

# shellcheck disable=SC1091
source "${repo_root}/workflow/support/gg_site_runtime.sh"

source_sif="${GG_SHIROKANE_PREBUILT_SIF:-}"
tag="${GG_SHIROKANE_SIF_TAG:-}"
expected_sha256="${GG_SHIROKANE_SIF_SHA256:-}"
expected_arch="${GG_SHIROKANE_EXPECTED_ARCH:-amd64}"
sif_dir="${GG_SHIROKANE_SIF_DIR:-${repo_root}/containers}"
stable_link="${GG_SHIROKANE_SIF_LINK:-${repo_root}/genegalleon.sif}"

if [[ -z "${source_sif}" ]]; then
  echo "GG_SHIROKANE_PREBUILT_SIF is required." >&2
  echo "Transfer a completed SIF to SHIROKANE before submitting this job." >&2
  exit 2
fi
if [[ -z "${tag}" ]]; then
  echo "GG_SHIROKANE_SIF_TAG is required." >&2
  echo "Use an immutable label such as release-v1.2.3 or sha256-<digest12>." >&2
  exit 2
fi
if [[ ! "${tag}" =~ ^[A-Za-z0-9][A-Za-z0-9_.-]*$ ]]; then
  echo "GG_SHIROKANE_SIF_TAG contains unsupported characters: ${tag}" >&2
  exit 2
fi
if [[ -z "${expected_sha256}" ]]; then
  echo "GG_SHIROKANE_SIF_SHA256 is required." >&2
  echo "Calculate the transferred SIF checksum before submitting this job." >&2
  exit 2
fi
if [[ ! "${expected_sha256}" =~ ^[A-Fa-f0-9]{64}$ ]]; then
  echo "GG_SHIROKANE_SIF_SHA256 must contain exactly 64 hexadecimal characters." >&2
  exit 2
fi
expected_sha256="$(printf '%s' "${expected_sha256}" | tr '[:upper:]' '[:lower:]')"
case "${expected_arch}" in
  amd64|x86_64)
    expected_arch="amd64"
    ;;
  *)
    echo "GG_SHIROKANE_EXPECTED_ARCH must be amd64 or x86_64: ${expected_arch}" >&2
    exit 2
    ;;
esac

case "${source_sif}" in
  /*)
    ;;
  *)
    source_sif="${repo_root}/${source_sif}"
    ;;
esac
if [[ ! -s "${source_sif}" ]]; then
  echo "Prebuilt SIF was not found or is empty: ${source_sif}" >&2
  exit 1
fi
source_sif="$(cd "$(dirname "${source_sif}")" && pwd -P)/$(basename "${source_sif}")"

mkdir -p "${sif_dir}"
sif_dir="$(cd "${sif_dir}" && pwd -P)"
sif_path="${sif_dir}/genegalleon-${tag}.sif"
partial_path="${sif_path}.partial.${JOB_ID:-$$}"
checksum_path="${sif_path}.sha256"
checksum_partial_path="${checksum_path}.partial.${JOB_ID:-$$}"
case "${stable_link}" in
  /*)
    ;;
  *)
    stable_link="${repo_root}/${stable_link}"
    ;;
esac
mkdir -p "$(dirname "${stable_link}")"
stable_link_dir="$(cd "$(dirname "${stable_link}")" && pwd -P)"
stable_link="${stable_link_dir}/$(basename "${stable_link}")"
stable_link_partial="${stable_link_dir}/.$(basename "${stable_link}").partial.${JOB_ID:-$$}"
lock_dir="${sif_dir}/.genegalleon-${tag}.install.lock"
created_partial=0
created_checksum_partial=0
created_link_partial=0
lock_acquired=0

cleanup() {
  local exit_code=$?
  if [[ "${created_partial}" -eq 1 && -f "${partial_path}" ]]; then
    rm -f -- "${partial_path}"
  fi
  if [[ "${created_checksum_partial}" -eq 1 && -f "${checksum_partial_path}" ]]; then
    rm -f -- "${checksum_partial_path}"
  fi
  if [[ "${created_link_partial}" -eq 1 && -L "${stable_link_partial}" ]]; then
    rm -f -- "${stable_link_partial}"
  fi
  if [[ "${lock_acquired}" -eq 1 ]]; then
    rm -f -- "${lock_dir}/owner"
    rmdir -- "${lock_dir}" 2>/dev/null || true
  fi
  return "${exit_code}"
}
trap cleanup EXIT

gg_shirokane_load_apptainer_module

verify_sha256() {
  local path=$1
  local actual_sha256
  actual_sha256="$(sha256sum "${path}" | awk '{print $1}')"
  if [[ "${actual_sha256}" != "${expected_sha256}" ]]; then
    echo "SHA-256 mismatch for ${path}" >&2
    echo "Expected: ${expected_sha256}" >&2
    echo "Actual:   ${actual_sha256}" >&2
    return 1
  fi
}

normalize_arch() {
  case "${1:-}" in
    amd64|x86_64)
      echo "amd64"
      ;;
    arm64|aarch64)
      echo "arm64"
      ;;
    *)
      printf '%s\n' "${1:-unknown}"
      ;;
  esac
}

validate_sif() {
  local path=$1
  local detected_arch=""

  apptainer inspect "${path}" >/dev/null
  verify_sha256 "${path}"
  if ! detected_arch="$(apptainer exec --cleanenv "${path}" uname -m)"; then
    echo "SIF execution check failed: ${path}" >&2
    return 1
  fi
  detected_arch="$(normalize_arch "${detected_arch}")"
  if [[ "${detected_arch}" != "${expected_arch}" ]]; then
    echo "SIF architecture mismatch for ${path}" >&2
    echo "Expected: ${expected_arch}" >&2
    echo "Actual:   ${detected_arch}" >&2
    return 1
  fi
}

write_checksum_sidecar() {
  local sif_basename

  sif_basename="$(basename "${sif_path}")"
  if [[ -e "${checksum_partial_path}" ]]; then
    echo "Refusing to overwrite an existing partial checksum: ${checksum_partial_path}" >&2
    return 1
  fi
  created_checksum_partial=1
  printf '%s  %s\n' "${expected_sha256}" "${sif_basename}" > "${checksum_partial_path}"
  mv -- "${checksum_partial_path}" "${checksum_path}"
  created_checksum_partial=0
}

update_stable_link() {
  local link_target="${sif_path}"

  if [[ -e "${stable_link}" && ! -L "${stable_link}" ]]; then
    echo "Refusing to replace a non-symlink container path: ${stable_link}" >&2
    return 1
  fi
  case "${sif_path}" in
    "${stable_link_dir}"/*)
      link_target="${sif_path#"${stable_link_dir}/"}"
      ;;
  esac
  if [[ -e "${stable_link_partial}" || -L "${stable_link_partial}" ]]; then
    echo "Refusing to overwrite an existing partial symlink: ${stable_link_partial}" >&2
    return 1
  fi
  created_link_partial=1
  ln -s -- "${link_target}" "${stable_link_partial}"
  if [[ "$(mv --help 2>&1 || true)" == *"--no-target-directory"* ]]; then
    mv -Tf -- "${stable_link_partial}" "${stable_link}"
  else
    mv -f -- "${stable_link_partial}" "${stable_link}"
  fi
  created_link_partial=0
}

echo "Validating prebuilt SIF: ${source_sif}"
validate_sif "${source_sif}"

if ! mkdir -- "${lock_dir}" 2>/dev/null; then
  echo "Another SIF installation is using this tag: ${tag}" >&2
  if [[ -s "${lock_dir}/owner" ]]; then
    echo "Lock owner: $(cat "${lock_dir}/owner")" >&2
  fi
  exit 1
fi
lock_acquired=1
printf 'job_id=%s host=%s pid=%s\n' "${JOB_ID:-local}" "$(hostname)" "$$" > "${lock_dir}/owner"

if [[ -s "${sif_path}" ]]; then
  echo "Reusing existing SIF: ${sif_path}"
  validate_sif "${sif_path}"
else
  if [[ -e "${partial_path}" ]]; then
    echo "Refusing to overwrite an existing partial SIF: ${partial_path}" >&2
    exit 1
  fi
  echo "Installing versioned SIF: ${sif_path}"
  created_partial=1
  cp -- "${source_sif}" "${partial_path}"
  validate_sif "${partial_path}"
  mv -- "${partial_path}" "${sif_path}"
  created_partial=0
fi

write_checksum_sidecar
update_stable_link

echo "GeneGalleon SIF ready: ${sif_path}"
echo "Default SIF link: ${stable_link}"
apptainer inspect "${sif_path}"
