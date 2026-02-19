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

mkdir -p "$(dirname "${report_file}")"
printf 'tier\tenv\tcommand\tfound\n' > "${report_file}"

required_failed=0
while IFS=$'\t' read -r env_name command_name; do
  if ! check_one "required" "${env_name}" "${command_name}"; then
    required_failed=1
  fi
done < <(read_tsv "${required_file}")

while IFS=$'\t' read -r env_name command_name; do
  check_one "optional" "${env_name}" "${command_name}" || true
done < <(read_tsv "${optional_file}")

if [[ ${required_failed} -ne 0 ]]; then
  echo "[validate_runtime] Required command checks failed. See: ${report_file}"
  exit 1
fi

echo "[validate_runtime] Required checks passed. Report: ${report_file}"
