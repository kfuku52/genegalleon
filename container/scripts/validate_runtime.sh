#!/usr/bin/env bash
set -euo pipefail

required_file=""
optional_file=""
report_file=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --required)
      required_file=$2
      shift 2
      ;;
    --optional)
      optional_file=$2
      shift 2
      ;;
    --report)
      report_file=$2
      shift 2
      ;;
    *)
      echo "Unknown argument: $1"
      exit 1
      ;;
  esac
done

if [[ -z "${required_file}" || -z "${optional_file}" || -z "${report_file}" ]]; then
  echo "Usage: $0 --required <tsv> --optional <tsv> --report <tsv>"
  exit 1
fi

check_one() {
  local tier=$1
  local env_name=$2
  local command_name=$3
  local found=0
  if micromamba run -n "${env_name}" bash -lc "command -v ${command_name} >/dev/null 2>&1"; then
    found=1
  fi
  printf '%s\t%s\t%s\t%s\n' "${tier}" "${env_name}" "${command_name}" "${found}" >> "${report_file}"
  if [[ "${tier}" == "required" && "${found}" -eq 0 ]]; then
    echo "[validate_runtime] missing required command: env=${env_name} command=${command_name}"
  fi
  return $((1-found))
}

read_tsv() {
  local file=$1
  [[ -f "${file}" ]] || return 0
  sed -e 's/[[:space:]]*#.*$//' -e '/^[[:space:]]*$/d' "${file}"
}

check_trinity_salmon_compatibility() {
  local env_name="base"
  local testdata_dir
  local tmp_dir
  local validation_log
  local version_info
  local passed=1

  if ! micromamba run -n "${env_name}" bash -lc "command -v Trinity >/dev/null 2>&1"; then
    return 0
  fi

  testdata_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/../testdata/salmon_trinity_smoke" && pwd)"
  tmp_dir="$(mktemp -d)"
  validation_log="${tmp_dir}/validation.log"

  version_info="$(micromamba run -n "${env_name}" salmon --version 2>&1 || true)"
  if [[ ! "${version_info}" =~ ^salmon[[:space:]]+1\. ]]; then
    echo "[validate_runtime] Trinity requires Salmon 1.x; found: ${version_info}"
    passed=0
  fi

  if [[ ${passed} -eq 1 ]] && ! micromamba run -n "${env_name}" \
      salmon --no-version-check index \
      -t "${testdata_dir}/transcripts.fa" \
      -i "${tmp_dir}/salmon.index" \
      -k 25 \
      -p 1 >"${validation_log}" 2>&1; then
    echo "[validate_runtime] Salmon 1.x failed to build the Trinity compatibility index."
    cat "${validation_log}"
    passed=0
  fi

  if [[ ${passed} -eq 1 ]] && ! micromamba run -n "${env_name}" \
      salmon --no-version-check quant \
      -i "${tmp_dir}/salmon.index" \
      -l U \
      -r "${testdata_dir}/reads.fa" \
      -o "${tmp_dir}/salmon.quant" \
      -p 1 \
      --minAssignedFrags 1 \
      --validateMappings >"${validation_log}" 2>&1; then
    echo "[validate_runtime] Salmon 1.x rejected the options required by Trinity."
    cat "${validation_log}"
    passed=0
  fi

  rm -rf -- "${tmp_dir}"
  printf '%s\t%s\t%s\t%s\n' \
    "required" "${env_name}" "Trinity/Salmon compatibility" "${passed}" >> "${report_file}"
  [[ ${passed} -eq 1 ]]
}

mkdir -p "$(dirname "${report_file}")"
printf 'tier\tenv\tcommand\tfound\n' > "${report_file}"

required_failed=0
while IFS=$'\t' read -r env_name command_name; do
  if ! check_one "required" "${env_name}" "${command_name}"; then
    required_failed=1
  fi
done < <(read_tsv "${required_file}")

if micromamba run -n base bash -lc "command -v csubst >/dev/null 2>&1"; then
  if micromamba run -n base csubst scan -h >/dev/null 2>&1; then
    printf '%s\t%s\t%s\t%s\n' "required" "base" "csubst scan" "1" >> "${report_file}"
  else
    printf '%s\t%s\t%s\t%s\n' "required" "base" "csubst scan" "0" >> "${report_file}"
    echo "[validate_runtime] csubst is installed but the scan subcommand is unavailable."
    required_failed=1
  fi
fi

while IFS=$'\t' read -r env_name command_name; do
  check_one "optional" "${env_name}" "${command_name}" || true
done < <(read_tsv "${optional_file}")

if ! check_trinity_salmon_compatibility; then
  required_failed=1
fi

if [[ ${required_failed} -ne 0 ]]; then
  echo "[validate_runtime] Required command checks failed. See: ${report_file}"
  exit 1
fi

echo "[validate_runtime] Required checks passed. Report: ${report_file}"
