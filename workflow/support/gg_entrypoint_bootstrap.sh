#!/usr/bin/env bash
set -euo pipefail

gg_search_workflow_dir_from_base() {
  local base_dir=${1:-}
  local search_dir=""
  local parent_dir=""

  [[ -n "${base_dir}" ]] || return 1
  search_dir="$(cd "$(dirname "${base_dir}")" 2>/dev/null && pwd -P || true)"
  if [[ -d "${base_dir}" ]]; then
    search_dir="$(cd "${base_dir}" 2>/dev/null && pwd -P || true)"
  fi

  while [[ -n "${search_dir}" ]]; do
    if [[ -f "${search_dir}/support/gg_util.sh" ]]; then
      (cd "${search_dir}" && pwd -P)
      return 0
    fi
    if [[ -f "${search_dir}/workflow/support/gg_util.sh" ]]; then
      (cd "${search_dir}/workflow" && pwd -P)
      return 0
    fi
    parent_dir="$(dirname "${search_dir}")"
    if [[ "${parent_dir}" == "${search_dir}" ]]; then
      break
    fi
    search_dir="${parent_dir}"
  done

  return 1
}

gg_validate_explicit_runtime_release() {
  local release_root=${1:-}
  local expected_workflow=${2:-}
  local bootstrap_path=${3:-${BASH_SOURCE[0]:-}}
  local snapshot_sha256=${KFAUTO_RUNTIME_RELEASE_SNAPSHOT_SHA256:-}

  [[ -n "${release_root}" && -n "${expected_workflow}" ]] || {
    echo "Explicit runtime release requires a release root and workflow identity." >&2
    return 1
  }
  [[ "${snapshot_sha256}" =~ ^[0-9a-f]{64}$ ]] || {
    echo "Explicit runtime release snapshot marker is missing or invalid." >&2
    return 1
  }

  python - "${release_root}" "${expected_workflow}" "${bootstrap_path}" "${snapshot_sha256}" <<'PY'
import hashlib
import json
import os
import pathlib
import stat
import sys

requested_release_root = pathlib.Path(sys.argv[1])
expected_workflow = sys.argv[2]
bootstrap_path = pathlib.Path(sys.argv[3])
expected_snapshot = sys.argv[4]


def fail(message):
    raise SystemExit(message)


def canonical_sha256(payload):
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


if requested_release_root.is_symlink() or not requested_release_root.is_dir():
    fail(
        "Explicit runtime release root is unavailable or symlinked: {}".format(
            requested_release_root
        )
    )
release_root = requested_release_root.resolve(strict=True)
workflow_dir = release_root / "workflow"
expected_bootstrap = workflow_dir / "support/gg_entrypoint_bootstrap.sh"
manifest_path = release_root / ".kfauto-release.json"
for required in (expected_bootstrap, manifest_path):
    if required.is_symlink() or not required.is_file():
        fail("Explicit runtime release contract file is unavailable or symlinked: {}".format(required))
if bootstrap_path.resolve(strict=True) != expected_bootstrap.resolve(strict=True):
    fail("Explicit runtime release bootstrap was shadowed: {}".format(bootstrap_path))

try:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
    fail("Explicit runtime release manifest is invalid: {}".format(exc))
if manifest.get("kind") != "immutable-genegalleon-runtime-release-v1":
    fail("Explicit runtime release manifest kind is invalid.")
if manifest.get("workflow") != expected_workflow:
    fail(
        "Explicit runtime release workflow changed: expected={} observed={}".format(
            expected_workflow, manifest.get("workflow")
        )
    )
manifest_release_root_text = manifest.get("release_root")
if not isinstance(manifest_release_root_text, str) or not manifest_release_root_text:
    fail("Explicit runtime release manifest root is invalid.")
manifest_release_root = pathlib.Path(manifest_release_root_text)
if not manifest_release_root.is_absolute():
    fail("Explicit runtime release manifest root is not absolute.")
try:
    manifest_physical_release_root = manifest_release_root.resolve(strict=True)
except OSError as exc:
    fail("Explicit runtime release manifest root is unavailable: {}".format(exc))
if manifest_physical_release_root != release_root:
    fail("Explicit runtime release manifest root changed.")
target_ref = str(manifest.get("target_ref", ""))
if len(target_ref) != 40 or any(char not in "0123456789abcdef" for char in target_ref):
    fail("Explicit runtime release target ref is invalid.")
contract_sha256 = str(manifest.get("contract_sha256", ""))
if len(contract_sha256) != 64 or any(char not in "0123456789abcdef" for char in contract_sha256):
    fail("Explicit runtime release contract hash is invalid.")
contract = {key: value for key, value in manifest.items() if key != "contract_sha256"}
if canonical_sha256(contract) != contract_sha256:
    fail("Explicit runtime release manifest contract hash changed.")

installed = manifest.get("installed_files")
if not isinstance(installed, list):
    fail("Explicit runtime release file inventory is unavailable.")
inventory = {}
for item in installed:
    if not isinstance(item, dict) or set(item) != {"path", "mode", "sha256"}:
        fail("Explicit runtime release file inventory is invalid.")
    relative = item["path"]
    mode = item["mode"]
    digest = item["sha256"]
    if not isinstance(relative, str) or not relative or "\x00" in relative:
        fail("Explicit runtime release file inventory contains an invalid path.")
    relative_path = pathlib.PurePosixPath(relative)
    if (
        relative_path.is_absolute()
        or relative_path.as_posix() != relative
        or any(part in {"", ".", ".."} for part in relative_path.parts)
    ):
        fail("Explicit runtime release file inventory contains an unsafe path: {}".format(relative))
    if not isinstance(mode, int) or isinstance(mode, bool) or not 0 <= mode <= 0o7777:
        fail("Explicit runtime release file inventory contains an invalid mode: {}".format(relative))
    if (
        not isinstance(digest, str)
        or len(digest) != 64
        or any(char not in "0123456789abcdef" for char in digest)
        or relative in inventory
    ):
        fail("Explicit runtime release file inventory contains an invalid entry.")
    inventory[relative] = {
        "mode": mode,
        "sha256": digest,
        "parts": relative_path.parts,
    }

snapshot_paths = {
    "VERSION",
    "workflow/{}_entrypoint.sh".format(expected_workflow),
    "workflow/core/{}_core.sh".format(expected_workflow),
    "workflow/support/artifact_provenance.py",
    "workflow/support/gg_entrypoint_bootstrap.sh",
    "workflow/support/gg_util.sh",
}
required_paths = snapshot_paths.union(
    {
        "workflow/gg_common_params.sh",
        "workflow/gg_path_defaults.sh",
        "workflow/support/gg_core_bootstrap.sh",
        "workflow/support/gg_entrypoint_config_vars.sh",
        "workflow/support/gg_shared_lock.sh",
        "workflow/support/gg_site_runtime.sh",
    }
)
required_paths.update(
    "workflow/support/gg_util/{:02d}_{}.sh".format(number, name)
    for number, name in enumerate(
        (
            "runtime_config",
            "container_scheduler",
            "species_helpers",
            "busco_runtime",
            "workspace_io",
            "workspace_validation",
            "execution_runtime",
            "execution_reporting",
            "sequence_databases",
            "reference_databases",
            "fasta",
        ),
        start=1,
    )
)
if not required_paths.issubset(inventory):
    fail("Explicit runtime release support closure is incomplete.")

verified_sha256 = {}
for relative in sorted(inventory):
    item = inventory[relative]
    path = release_root.joinpath(*item["parts"])
    current = release_root
    for part in item["parts"]:
        current /= part
        if current.is_symlink():
            fail("Explicit runtime release file is unavailable or symlinked: {}".format(relative))
    if not path.is_file():
        fail("Explicit runtime release file is unavailable or symlinked: {}".format(relative))
    try:
        resolved_path = path.resolve(strict=True)
        resolved_path.relative_to(release_root)
        payload = path.read_bytes()
        observed_mode = stat.S_IMODE(path.stat().st_mode)
    except (OSError, ValueError) as exc:
        fail("Explicit runtime release file is unsafe or unreadable: {} ({})".format(relative, exc))
    observed = hashlib.sha256(payload).hexdigest()
    if observed != item["sha256"] or observed_mode != item["mode"]:
        fail("Explicit runtime release file changed: {}".format(relative))
    verified_sha256[relative] = observed

file_sha256 = {relative: verified_sha256[relative] for relative in sorted(snapshot_paths)}

snapshot = {
    "schema_version": 1,
    "kind": "gridengine-runtime-release-snapshot-v1",
    "project_root": manifest.get("project_root"),
    "workflow": expected_workflow,
    "release_root": manifest_release_root_text,
    "target_ref": target_ref,
    "release_contract_sha256": contract_sha256,
    "file_sha256": file_sha256,
}
observed_snapshot = canonical_sha256(snapshot)
if observed_snapshot != expected_snapshot:
    fail(
        "Explicit runtime release snapshot hash changed: expected={} observed={}".format(
            expected_snapshot, observed_snapshot
        )
    )
print(workflow_dir)
PY
}

gg_find_workflow_dir() {
  local entrypoint_path=${1:-}
  local bootstrap_path=${2:-${BASH_SOURCE[0]:-}}
  local expected_workflow=${3:-}
  local base_dir=""

  if [[ -n "${KFAUTO_RUNTIME_RELEASE_ROOT:-}" ]]; then
    gg_validate_explicit_runtime_release \
      "${KFAUTO_RUNTIME_RELEASE_ROOT}" \
      "${expected_workflow}" \
      "${BASH_SOURCE[0]:-}"
    return $?
  fi

  for base_dir in "${entrypoint_path}" "${bootstrap_path}"; do
    if gg_search_workflow_dir_from_base "${base_dir}"; then
      return 0
    fi
  done

  for base_dir in "${SLURM_SUBMIT_DIR:-}" "${PBS_O_WORKDIR:-}" "${PWD:-}"; do
    [[ -n "${base_dir}" ]] || continue
    if gg_search_workflow_dir_from_base "${base_dir}"; then
      return 0
    fi
  done

  return 1
}

gg_set_workflow_dir() {
  local entrypoint_path=${1:-${BASH_SOURCE[1]:-}}
  local bootstrap_path=${2:-${BASH_SOURCE[0]:-}}
  local expected_workflow=${3:-}
  local resolved_workflow_dir=""

  resolved_workflow_dir="$(gg_find_workflow_dir "${entrypoint_path}" "${bootstrap_path}" "${expected_workflow}")" || {
    echo "Failed to resolve workflow directory from entrypoint_path=${entrypoint_path:-unset}" >&2
    return 1
  }

  gg_workflow_dir="${resolved_workflow_dir}"
}

gg_source_common_params_if_available() {
  local resolved_workflow_dir=${1:-${gg_workflow_dir:-}}

  if [[ -n "${resolved_workflow_dir}" && -s "${resolved_workflow_dir}/gg_common_params.sh" ]]; then
    # shellcheck disable=SC1090
    source "${resolved_workflow_dir}/gg_common_params.sh"
  fi
}

gg_configure_python_pycacheprefix() {
  local default_pycache_prefix=""

  if [[ -n "${PYTHONPYCACHEPREFIX:-}" ]]; then
    return 0
  fi

  default_pycache_prefix="${TMPDIR:-/tmp}/genegalleon_pycache"
  mkdir -p -- "${default_pycache_prefix}" 2>/dev/null || true
  export PYTHONPYCACHEPREFIX="${default_pycache_prefix}"
}

gg_entrypoint_initialize() {
  local entrypoint_path=${1:-${BASH_SOURCE[1]:-}}
  local load_common_params=${2:-0}
  local expected_workflow=${3:-}

  if ! gg_set_workflow_dir "${entrypoint_path}" "${BASH_SOURCE[0]:-}" "${expected_workflow}"; then
    return 1
  fi
  # shellcheck disable=SC1090
  source "${gg_workflow_dir}/gg_path_defaults.sh"
  if [[ "${load_common_params}" -eq 1 ]]; then
    gg_source_common_params_if_available "${gg_workflow_dir}"
  fi
  gg_configure_python_pycacheprefix
}
