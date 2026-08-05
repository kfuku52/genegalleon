#!/usr/bin/env bash
set -euo pipefail

gg_workspace_storage_script_dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
# shellcheck disable=SC1091
source "${gg_workspace_storage_script_dir}/gg_common_params.sh"
export GG_COMMON_GENE_FAMILY_LARGE_ZIP_WARNING_BYTES
export GG_COMMON_GENE_FAMILY_FINAL_ZIP_MAX_BYTES
exec python3 "${gg_workspace_storage_script_dir}/support/workspace_output_storage.py" "$@"
