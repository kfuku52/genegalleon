#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${gg_workflow_dir:-}" ]]; then
  echo "gg_path_defaults.sh requires gg_workflow_dir to be set." >&2
  return 1 2>/dev/null || exit 1
fi

: "${gg_support_dir:=${gg_workflow_dir}/support}"
: "${gg_core_dir:=${gg_workflow_dir}/core}"
: "${gg_workspace_dir:=${gg_workflow_dir}/../workspace}"
: "${gg_container_image_path:=${gg_workflow_dir}/../genegalleon.sif}"
