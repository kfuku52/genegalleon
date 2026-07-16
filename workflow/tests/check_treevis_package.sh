#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
package_dir="${repo_root}/workflow/support/treevis"
work_dir="$(mktemp -d)"

cleanup() {
  rm -rf -- "${work_dir}"
}
trap cleanup EXIT

cd "${work_dir}"
R CMD build "${package_dir}"
package_tarball="$(find "${work_dir}" -maxdepth 1 -name 'genegalleon.treevis_*.tar.gz' -print -quit)"
if [[ -z "${package_tarball}" ]]; then
  echo "treevis package tarball was not created." >&2
  exit 1
fi
R CMD check --no-manual --no-build-vignettes "${package_tarball}"
check_log="$(find "${work_dir}" -maxdepth 2 -name '00check.log' -print -quit)"
if [[ -z "${check_log}" ]]; then
  echo "treevis R CMD check log was not created." >&2
  exit 1
fi
check_status="$(grep '^Status:' "${check_log}" || true)"
if [[ "${check_status}" != "Status: OK" ]]; then
  cat "${check_log}" >&2
  echo "treevis R CMD check reported warnings, notes, or errors." >&2
  exit 1
fi
