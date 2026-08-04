#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
sif_path="${GENEGALLEON_SIF:-${gg_container_image_path:-${repo_root}/genegalleon.sif}}"
site_runtime="${repo_root}/workflow/support/gg_site_runtime.sh"

usage() {
  cat <<'EOF'
Usage:
  bash workflow/tests/run_in_sif.sh <command> [args...]

Runs a validation command inside the GeneGalleon SIF runtime.

Environment overrides:
  GENEGALLEON_SIF=/path/to/genegalleon.sif
  gg_container_image_path=/path/to/genegalleon.sif
  GENEGALLEON_SIF_EXTRA_BINDS=/absolute/read-only/path

Multiple extra bind paths may be separated by newlines. Each source is mounted
at the same absolute path inside the container with read-only permissions.

Examples:
  bash workflow/tests/run_in_sif.sh python -m pytest -q workflow/tests/test_hgt_end_to_end.py
  bash workflow/tests/run_in_sif.sh Rscript workflow/tests/test_treevis_main.R
EOF
}

if [[ $# -eq 0 ]]; then
  usage >&2
  exit 2
fi

if [[ ! -f "${sif_path}" ]]; then
  cat >&2 <<EOF
GeneGalleon SIF was not found: ${sif_path}

Build or place the SIF at the repository root as ./genegalleon.sif, or set:
  GENEGALLEON_SIF=/path/to/genegalleon.sif
EOF
  exit 1
fi

if [[ ! -r "${site_runtime}" ]]; then
  echo "GeneGalleon site runtime helper was not found: ${site_runtime}" >&2
  exit 1
fi
# shellcheck disable=SC1090
source "${site_runtime}"
gg_prepend_container_runtime_path \
  "${GG_CONTAINER_RUNTIME_PACKAGE_ROOT:-/opt/pkg}" \
  "${GG_CONTAINER_RUNTIME_LEGACY_DIR:-/bio/package/singularity/singularity_3.0/bin}"

container_engine=""
if command -v apptainer >/dev/null 2>&1; then
  container_engine="$(command -v apptainer)"
elif command -v singularity >/dev/null 2>&1; then
  container_engine="$(command -v singularity)"
fi

if [[ -z "${container_engine}" ]]; then
  cat >&2 <<'EOF'
Neither apptainer nor singularity was found.

Install Apptainer/Singularity on this host before claiming SIF validation.
EOF
  exit 1
fi

bind_args=(--bind "${repo_root}:${repo_root}")
if [[ -n "${GENEGALLEON_SIF_EXTRA_BINDS:-}" ]]; then
  while IFS= read -r bind_path; do
    [[ -z "${bind_path}" ]] && continue
    if [[ "${bind_path}" != /* || "${bind_path}" == *:* || ! -e "${bind_path}" ]]; then
      echo "Invalid GENEGALLEON_SIF_EXTRA_BINDS path: ${bind_path}" >&2
      exit 1
    fi
    bind_args+=(--bind "${bind_path}:${bind_path}:ro")
  done <<< "${GENEGALLEON_SIF_EXTRA_BINDS}"
fi

exec "${container_engine}" exec \
  --cleanenv \
  "${bind_args[@]}" \
  --pwd "${repo_root}" \
  "${sif_path}" \
  "$@"
