#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
sif_path="${GENEGALLEON_SIF:-${gg_container_image_path:-${repo_root}/genegalleon.sif}}"

usage() {
  cat <<'EOF'
Usage:
  bash workflow/tests/run_in_sif.sh <command> [args...]

Runs a validation command inside the GeneGalleon SIF runtime.

Environment overrides:
  GENEGALLEON_SIF=/path/to/genegalleon.sif
  gg_container_image_path=/path/to/genegalleon.sif

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

exec "${container_engine}" exec \
  --cleanenv \
  --bind "${repo_root}:${repo_root}" \
  --pwd "${repo_root}" \
  "${sif_path}" \
  "$@"
