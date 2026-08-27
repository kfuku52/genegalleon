#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../source_branches.env
source "${script_dir}/../source_branches.env"

format="env"
scope="all"

usage() {
  cat <<'EOF'
Usage: bash container/scripts/resolve_source_revisions.sh [--format env|tsv] [--scope all|owned]

Resolve the configured moving source branches in parallel. Explicit *_REPO_SHA
environment overrides are preserved for one-off reproduction/debugging builds.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --format)
      format="${2:-}"
      shift 2
      ;;
    --scope)
      scope="${2:-}"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      usage >&2
      exit 2
      ;;
  esac
done

case "${format}" in
  env|tsv) ;;
  *)
    echo "--format must be env or tsv" >&2
    exit 2
    ;;
esac
case "${scope}" in
  all|owned) ;;
  *)
    echo "--scope must be all or owned" >&2
    exit 2
    ;;
esac

BUSCO_REPO_URL=${BUSCO_REPO_URL:-https://gitlab.com/ezlab/busco.git}
PAML_REPO_URL=${PAML_REPO_URL:-https://github.com/iqtree/paml.git}
KFL1OU_REPO_URL=${KFL1OU_REPO_URL:-https://github.com/kfuku52/kfl1ou.git}
KFFRACTBIAS_REPO_URL=${KFFRACTBIAS_REPO_URL:-https://github.com/kfuku52/kfFractBias.git}
KFTOOLS_REPO_URL=${KFTOOLS_REPO_URL:-https://github.com/kfuku52/kftools.git}
RKFTOOLS_REPO_URL=${RKFTOOLS_REPO_URL:-https://github.com/kfuku52/rkftools.git}
RADTE_REPO_URL=${RADTE_REPO_URL:-https://github.com/kfuku52/RADTE.git}

KFU52_AMALGKIT_REPO_REF=${KFU52_AMALGKIT_REPO_REF:-${GG_SOURCE_AMALGKIT_REPO_REF}}
KFU52_CDSKIT_REPO_REF=${KFU52_CDSKIT_REPO_REF:-${GG_SOURCE_CDSKIT_REPO_REF}}
KFU52_CSUBST_REPO_REF=${KFU52_CSUBST_REPO_REF:-${GG_SOURCE_CSUBST_REPO_REF}}
KFU52_NWKIT_REPO_REF=${KFU52_NWKIT_REPO_REF:-${GG_SOURCE_NWKIT_REPO_REF}}
BUSCO_REPO_REF=${BUSCO_REPO_REF:-${GG_SOURCE_BUSCO_REPO_REF}}
PAML_REPO_REF=${PAML_REPO_REF:-${GG_SOURCE_PAML_REPO_REF}}
KFL1OU_REPO_REF=${KFL1OU_REPO_REF:-${GG_SOURCE_KFL1OU_REPO_REF}}
KFFRACTBIAS_REPO_REF=${KFFRACTBIAS_REPO_REF:-${GG_SOURCE_KFFRACTBIAS_REPO_REF}}
KFTOOLS_REPO_REF=${KFTOOLS_REPO_REF:-${GG_SOURCE_KFTOOLS_REPO_REF}}
RKFTOOLS_REPO_REF=${RKFTOOLS_REPO_REF:-${GG_SOURCE_RKFTOOLS_REPO_REF}}
RADTE_REPO_REF=${RADTE_REPO_REF:-${GG_SOURCE_RADTE_REPO_REF}}

records=(
  "amalgkit|KFU52_AMALGKIT_REPO_SHA|https://github.com/kfuku52/amalgkit.git|${KFU52_AMALGKIT_REPO_REF}|owned"
  "cdskit|KFU52_CDSKIT_REPO_SHA|https://github.com/kfuku52/cdskit.git|${KFU52_CDSKIT_REPO_REF}|owned"
  "csubst|KFU52_CSUBST_REPO_SHA|https://github.com/kfuku52/csubst.git|${KFU52_CSUBST_REPO_REF}|owned"
  "nwkit|KFU52_NWKIT_REPO_SHA|https://github.com/kfuku52/nwkit.git|${KFU52_NWKIT_REPO_REF}|owned"
  "BUSCO|BUSCO_REPO_SHA|${BUSCO_REPO_URL}|${BUSCO_REPO_REF}|third-party"
  "paml|PAML_REPO_SHA|${PAML_REPO_URL}|${PAML_REPO_REF}|third-party"
  "kfl1ou|KFL1OU_REPO_SHA|${KFL1OU_REPO_URL}|${KFL1OU_REPO_REF}|owned"
  "kfFractBias|KFFRACTBIAS_REPO_SHA|${KFFRACTBIAS_REPO_URL}|${KFFRACTBIAS_REPO_REF}|owned"
  "kftools|KFTOOLS_REPO_SHA|${KFTOOLS_REPO_URL}|${KFTOOLS_REPO_REF}|owned"
  "rkftools|RKFTOOLS_REPO_SHA|${RKFTOOLS_REPO_URL}|${RKFTOOLS_REPO_REF}|owned"
  "RADTE|RADTE_REPO_SHA|${RADTE_REPO_URL}|${RADTE_REPO_REF}|owned"
)

work_dir="$(mktemp -d)"
cleanup() {
  rm -rf -- "${work_dir}"
}
trap cleanup EXIT

pids=()
indexes=()
for index in "${!records[@]}"; do
  IFS='|' read -r source_name sha_var repo_url repo_ref ownership <<< "${records[$index]}"
  if [[ "${scope}" == "owned" && "${ownership}" != "owned" ]]; then
    continue
  fi
  (
    resolved_sha="${!sha_var:-}"
    if [[ -z "${resolved_sha}" ]]; then
      resolved_sha="$(bash "${script_dir}/resolve_git_branch_sha.sh" "${repo_url}" "${repo_ref}")"
    fi
    if [[ ! "${resolved_sha}" =~ ^[0-9a-f]{40}$ ]]; then
      echo "Invalid resolved revision for ${source_name}: ${resolved_sha}" >&2
      exit 1
    fi
    printf '%s\t%s\t%s\n' "${source_name}" "${sha_var}" "${resolved_sha}" > "${work_dir}/${index}"
  ) &
  pids+=("$!")
  indexes+=("${index}")
done

failed=0
for pid in "${pids[@]}"; do
  if ! wait "${pid}"; then
    failed=1
  fi
done
if [[ "${failed}" == "1" ]]; then
  echo "One or more upstream source revisions could not be resolved." >&2
  exit 1
fi

if [[ "${format}" == "tsv" ]]; then
  printf 'source\trevision\n'
fi
for index in "${indexes[@]}"; do
  IFS=$'\t' read -r source_name sha_var resolved_sha < "${work_dir}/${index}"
  if [[ "${format}" == "env" ]]; then
    printf '%s=%s\n' "${sha_var}" "${resolved_sha}"
  else
    printf '%s\t%s\n' "${source_name}" "${resolved_sha}"
  fi
done
