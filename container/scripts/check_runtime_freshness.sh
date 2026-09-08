#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
mode="${GG_RUNTIME_FRESHNESS:-daily}"
scope="${GG_RUNTIME_FRESHNESS_SCOPE:-owned}"
runtime_kind=""
runtime_ref=""
engine=""
expected_hash=""

usage() {
  cat <<'EOF'
Usage:
  bash container/scripts/check_runtime_freshness.sh --docker IMAGE
  bash container/scripts/check_runtime_freshness.sh --sif PATH --engine apptainer|singularity

GG_RUNTIME_FRESHNESS=daily  Resolve moving upstreams at most once per UTC day (default).
GG_RUNTIME_FRESHNESS=always Resolve moving upstreams on every check.
GG_RUNTIME_FRESHNESS=off    Skip the freshness check for intentional offline use.
GG_RUNTIME_FRESHNESS_SCOPE=owned checks repository-owned upstreams (default).
GG_RUNTIME_FRESHNESS_SCOPE=all also checks BUSCO and PAML branch tips.

--expected-hash SHA256 compares only the embedded runtime identity. CI uses
this after restoring a cache entry and does not resolve moving branches again.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --docker)
      runtime_kind="docker"
      runtime_ref="${2:-}"
      shift 2
      ;;
    --sif)
      runtime_kind="sif"
      runtime_ref="${2:-}"
      shift 2
      ;;
    --engine)
      engine="${2:-}"
      shift 2
      ;;
    --expected-hash)
      expected_hash="${2:-}"
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

case "${mode}" in
  daily|always|off) ;;
  *)
    echo "GG_RUNTIME_FRESHNESS must be one of: daily, always, off" >&2
    exit 2
    ;;
esac
case "${scope}" in
  owned|all) ;;
  *)
    echo "GG_RUNTIME_FRESHNESS_SCOPE must be owned or all" >&2
    exit 2
    ;;
esac
if [[ -z "${runtime_kind}" || -z "${runtime_ref}" ]]; then
  usage >&2
  exit 2
fi
if [[ -n "${expected_hash}" && ! "${expected_hash}" =~ ^[0-9a-f]{64}$ ]]; then
  echo "--expected-hash must be a 64-character lowercase SHA-256" >&2
  exit 2
fi
if [[ "${mode}" == "off" && -z "${expected_hash}" ]]; then
  exit 0
fi

runtime_hash=""
security_epoch=""
runtime_arch=""
build_target=""
if [[ "${runtime_kind}" == "docker" ]]; then
  metadata="$(docker image inspect --format '{{.Architecture}}|{{index .Config.Labels "io.genegalleon.runtime-input"}}|{{index .Config.Labels "io.genegalleon.security-refresh-epoch"}}|{{index .Config.Labels "io.genegalleon.build-target"}}' "${runtime_ref}")"
  IFS='|' read -r runtime_arch runtime_hash security_epoch build_target <<< "${metadata}"
else
  if [[ -z "${engine}" ]]; then
    echo "--engine is required with --sif" >&2
    exit 2
  fi
  inspect_json="$("${engine}" inspect --json "${runtime_ref}")"
  metadata="$(python3 -c '
import json, sys
payload = json.load(sys.stdin)
found = {}
def visit(value):
    if isinstance(value, dict):
        for key, child in value.items():
            if key in {"io.genegalleon.runtime-input", "io.genegalleon.security-refresh-epoch", "io.genegalleon.build-target"}:
                found[key] = child
            visit(child)
    elif isinstance(value, list):
        for child in value:
            visit(child)
visit(payload)
print("|".join(str(found.get(key, "")) for key in (
    "io.genegalleon.runtime-input", "io.genegalleon.security-refresh-epoch", "io.genegalleon.build-target")))
' <<< "${inspect_json}")"
  IFS='|' read -r runtime_hash security_epoch build_target <<< "${metadata}"
  case "$(uname -m)" in
    x86_64|amd64) runtime_arch="amd64" ;;
    aarch64|arm64) runtime_arch="arm64" ;;
    *)
      echo "Unsupported runtime architecture: $(uname -m)" >&2
      exit 1
      ;;
  esac
fi

if [[ ! "${runtime_hash}" =~ ^[0-9a-f]{64}$ || ! "${security_epoch}" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ \
   || ( "${build_target}" != "runtime" && "${build_target}" != "development" ) ]]; then
  cat >&2 <<EOF
GeneGalleon runtime freshness metadata is missing or invalid: ${runtime_ref}

Rebuild the runtime from this checkout, or set GG_RUNTIME_FRESHNESS=off only
for an intentional offline check with a known older image.
EOF
  exit 1
fi
if [[ -n "${expected_hash}" ]]; then
  if [[ "${runtime_hash}" != "${expected_hash}" ]]; then
    echo "GeneGalleon runtime identity mismatch: observed=${runtime_hash} expected=${expected_hash}" >&2
    exit 1
  fi
  echo "[runtime freshness] ${runtime_ref} has exact runtime input ${runtime_hash}."
  exit 0
fi

sha256_file() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}

sha256_stream() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum | awk '{print $1}'
  else
    shasum -a 256 | awk '{print $1}'
  fi
}

git_dir="$(git -C "${repo_root}" rev-parse --git-dir 2>/dev/null || true)"
if [[ -n "${git_dir}" && "${git_dir}" != /* ]]; then
  git_dir="${repo_root}/${git_dir}"
fi
private_cache_parent=""
if [[ -n "${GG_RUNTIME_FRESHNESS_CACHE_DIR:-}" ]]; then
  cache_root="${GG_RUNTIME_FRESHNESS_CACHE_DIR}"
elif [[ -n "${git_dir}" ]]; then
  cache_root="${git_dir}/gg-runtime-freshness"
else
  private_cache_parent="${TMPDIR:-/tmp}/genegalleon_$(id -u)"
  if [[ -L "${private_cache_parent}" || ( -e "${private_cache_parent}" && ( ! -d "${private_cache_parent}" || ! -O "${private_cache_parent}" ) ) ]]; then
    echo "Refusing unsafe runtime freshness cache parent: ${private_cache_parent}" >&2
    exit 1
  fi
  (umask 077; mkdir -p -- "${private_cache_parent}")
  if [[ -L "${private_cache_parent}" || ! -d "${private_cache_parent}" || ! -O "${private_cache_parent}" ]]; then
    echo "Runtime freshness cache parent is not an owned directory: ${private_cache_parent}" >&2
    exit 1
  fi
  chmod 700 "${private_cache_parent}"
  cache_root="${private_cache_parent}/gg-runtime-freshness"
fi
source_config_hash="$(sha256_file "${repo_root}/container/source_branches.env")"
resolver_hash="$(sha256_file "${script_dir}/resolve_source_revisions.sh")"
sha_variables=(
  KFU52_AMALGKIT_REPO_SHA
  KFU52_CDSKIT_REPO_SHA
  KFU52_CSUBST_REPO_SHA
  KFU52_NWKIT_REPO_SHA
  BUSCO_REPO_SHA
  PAML_REPO_SHA
  KFL1OU_REPO_SHA
  KFFRACTBIAS_REPO_SHA
  KFTOOLS_REPO_SHA
  RKFTOOLS_REPO_SHA
  RADTE_REPO_SHA
)
resolution_variables=(
  "${sha_variables[@]}"
  KFU52_AMALGKIT_REPO_REF
  KFU52_CDSKIT_REPO_REF
  KFU52_CSUBST_REPO_REF
  KFU52_NWKIT_REPO_REF
  BUSCO_REPO_URL
  BUSCO_REPO_REF
  PAML_REPO_URL
  PAML_REPO_REF
  KFL1OU_REPO_URL
  KFL1OU_REPO_REF
  KFFRACTBIAS_REPO_URL
  KFFRACTBIAS_REPO_REF
  KFTOOLS_REPO_URL
  KFTOOLS_REPO_REF
  RKFTOOLS_REPO_URL
  RKFTOOLS_REPO_REF
  RADTE_REPO_URL
  RADTE_REPO_REF
)
override_fingerprint="$(
  for variable in "${resolution_variables[@]}"; do
    if [[ "${scope}" == "owned" && ( "${variable}" == BUSCO_REPO_* || "${variable}" == PAML_REPO_* ) ]]; then
      continue
    fi
    printf '%s=%s\n' "${variable}" "${!variable:-}"
  done | sha256_stream
)"
cache_file="${cache_root}/$(date -u +%F)-${source_config_hash}-${resolver_hash}-${scope}-${override_fingerprint}.env"

mkdir -p "${cache_root}"
if [[ "${mode}" == "always" || ! -s "${cache_file}" ]]; then
  cache_tmp="${cache_file}.tmp.$$"
  cleanup_cache_tmp() {
    rm -f -- "${cache_tmp}"
  }
  trap cleanup_cache_tmp EXIT HUP INT TERM
  if ! bash "${script_dir}/resolve_source_revisions.sh" --format env --scope "${scope}" > "${cache_tmp}"; then
    echo "Runtime freshness check failed while resolving upstream sources." >&2
    exit 1
  fi
  mv -f -- "${cache_tmp}" "${cache_file}"
  cache_tmp=""
  trap - EXIT HUP INT TERM
fi

# The resolver has already captured any explicit overrides that belong to the
# selected scope. Clear inherited values so third-party revisions in owned mode
# always come from the runtime itself rather than an unrelated caller shell.
for variable in "${sha_variables[@]}"; do
  unset "${variable}"
done
while IFS='=' read -r variable value; do
  if [[ ! "${variable}" =~ ^[A-Z0-9_]+_REPO_SHA$ || ! "${value}" =~ ^[0-9a-f]{40}$ ]]; then
    echo "Invalid runtime freshness cache record: ${variable}=${value}" >&2
    exit 1
  fi
  printf -v "${variable}" '%s' "${value}"
  export "${variable?}"
done < "${cache_file}"

missing_manifest_revisions=0
for variable in "${sha_variables[@]}"; do
  if [[ -z "${!variable:-}" ]]; then
    missing_manifest_revisions=1
    break
  fi
done

if [[ "${missing_manifest_revisions}" == "1" ]]; then
  if [[ "${runtime_kind}" == "docker" ]]; then
    source_manifest="$(
      docker run --rm --entrypoint cat \
        "${runtime_ref}" \
        /opt/pg/logs/source_revisions.tsv
    )"
  else
    source_manifest="$(
      "${engine}" exec \
        "${runtime_ref}" \
        cat /opt/pg/logs/source_revisions.tsv
    )"
  fi
  while IFS=$'\t' read -r source_name revision; do
    case "${source_name}" in
      source) continue ;;
      amalgkit) variable=KFU52_AMALGKIT_REPO_SHA ;;
      cdskit) variable=KFU52_CDSKIT_REPO_SHA ;;
      csubst) variable=KFU52_CSUBST_REPO_SHA ;;
      nwkit) variable=KFU52_NWKIT_REPO_SHA ;;
      BUSCO) variable=BUSCO_REPO_SHA ;;
      paml) variable=PAML_REPO_SHA ;;
      kfl1ou) variable=KFL1OU_REPO_SHA ;;
      kfFractBias) variable=KFFRACTBIAS_REPO_SHA ;;
      kftools) variable=KFTOOLS_REPO_SHA ;;
      rkftools) variable=RKFTOOLS_REPO_SHA ;;
      RADTE) variable=RADTE_REPO_SHA ;;
      *)
        echo "Unknown source in runtime manifest: ${source_name}" >&2
        exit 1
        ;;
    esac
    if [[ ! "${revision}" =~ ^[0-9a-f]{40}$ ]]; then
      echo "Invalid revision in runtime manifest for ${source_name}: ${revision}" >&2
      exit 1
    fi
    if [[ -z "${!variable:-}" ]]; then
      printf -v "${variable}" '%s' "${revision}"
      export "${variable?}"
    fi
  done <<< "${source_manifest}"
fi

for variable in "${sha_variables[@]}"; do
  if [[ ! "${!variable:-}" =~ ^[0-9a-f]{40}$ ]]; then
    echo "Runtime freshness check lacks an exact revision for ${variable}." >&2
    exit 1
  fi
done

export GG_BUILD_PLATFORMS="linux/${runtime_arch}"
export GG_BUILD_VCS_REF="runtime-freshness"
export GG_BUILD_VERSION="runtime-freshness"
export GG_BUILD_TARGET="${build_target}"
expected_hash="$(bash "${script_dir}/compute_build_input_hash.sh" --runtime "${security_epoch}")"

if [[ "${runtime_hash}" != "${expected_hash}" ]]; then
  cat >&2 <<EOF
GeneGalleon runtime is stale for the current container inputs or upstream revisions.
  runtime:  ${runtime_ref}
  observed: ${runtime_hash}
  expected: ${expected_hash}

Rebuild with:
  BUILD_SIF=0 IMAGE_SOURCE=local IMAGE=local/genegalleon TAG=dev bash ./gg_container_build_entrypoint.sh

Set GG_RUNTIME_FRESHNESS=off only for an intentional offline check.
EOF
  exit 1
fi

echo "[runtime freshness] ${runtime_ref} matches the daily ${scope} upstream snapshot and current container inputs."
