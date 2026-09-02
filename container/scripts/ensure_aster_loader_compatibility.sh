#!/usr/bin/env bash
set -euo pipefail

prefix="${1:-${MAMBA_ROOT_PREFIX:-/opt/conda}}"
python_bin="${prefix}/bin/python"

if [[ ! -x "${python_bin}" ]]; then
  echo "ASTER dependency correction requires an executable Python at ${python_bin}." >&2
  exit 1
fi

"${python_bin}" - "${prefix}" <<'PY'
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path


prefix = Path(sys.argv[1]).resolve()
fix_url = "https://github.com/bioconda/bioconda-recipes/pull/68643"
expected_hooks = {
    prefix / "etc/conda/activate.d/aster_activate.sh": (
        "2bfbab3165b7d77e0cca87c633053b14c36aa01192d461c5f25dd193075ac3f1"
    ),
    prefix / "etc/conda/deactivate.d/aster_deactivate.sh": (
        "953b87a56054dac5b788ec65334af7380fb094d81975cecf7ee0717e24837926"
    ),
}
installed_hooks = {
    path
    for directory in ("activate.d", "deactivate.d")
    for path in (prefix / "etc/conda" / directory).glob("aster_*.sh")
}

if not installed_hooks:
    print("[aster compatibility] Hook-free ASTER package detected; no correction needed.")
    raise SystemExit(0)

if installed_hooks != set(expected_hooks):
    found = ", ".join(str(path) for path in sorted(installed_hooks))
    raise SystemExit(
        "Refusing to modify an unexpected set of ASTER activation hooks: " + found
    )

records = sorted((prefix / "conda-meta").glob("aster-*.json"))
if len(records) != 1:
    raise SystemExit(f"Expected one ASTER package record, found {len(records)}: {records}")
metadata = json.loads(records[0].read_text(encoding="utf-8"))
identity = (
    metadata.get("name"),
    metadata.get("version"),
    metadata.get("build_number"),
)
if identity != ("aster", "1.25", 0):
    raise SystemExit(
        "Refusing to modify ASTER hooks for an unrecognized package identity: "
        f"name={identity[0]!r} version={identity[1]!r} build_number={identity[2]!r}"
    )

for path, expected_digest in expected_hooks.items():
    if path.is_symlink() or not path.is_file():
        raise SystemExit(f"Refusing to modify a non-regular ASTER hook: {path}")
    actual_digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual_digest != expected_digest:
        raise SystemExit(
            f"Refusing to modify changed ASTER hook {path}: sha256={actual_digest}"
        )

for path in expected_hooks:
    path.unlink()

log_dir = prefix.parent / "pg/logs"
log_dir.mkdir(parents=True, exist_ok=True)
(log_dir / "aster_dependency_correction.txt").write_text(
    "status=removed legacy loader hooks\n"
    f"package=aster {metadata['version']} {metadata.get('build', '')}\n"
    f"upstream_fix={fix_url}\n"
    "remove_after=hook-free ASTER is published by Bioconda and a normal solve passes\n",
    encoding="utf-8",
)
print(
    "[aster compatibility] Removed the two verified ASTER 1.25 build-0 loader hooks; "
    f"upstream fix: {fix_url}"
)
PY
