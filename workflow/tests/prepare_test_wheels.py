#!/usr/bin/env python3
"""Resolve the moving test dependency only in temporary wheel-build inputs."""

import argparse
import hashlib
import json
import re
from pathlib import Path

TEST_DIR = Path(__file__).resolve().parent


def prepare(directory: Path, csubst_sha: str):
    if not re.fullmatch(r"[0-9a-f]{40}", csubst_sha):
        raise ValueError("csubst source must be a resolved 40-character commit SHA")
    requirements = TEST_DIR / "requirements.txt"
    constraints = TEST_DIR / "requirements.lock.txt"
    build_lines = []
    install_lines = []
    found = 0
    for line in requirements.read_text().splitlines():
        if line.startswith("csubst @ git+https://github.com/kfuku52/csubst.git@"):
            found += 1
            build_lines.append("csubst @ git+https://github.com/kfuku52/csubst.git@" + csubst_sha)
            install_lines.append("csubst")
        elif "git+" in line and not line.lstrip().startswith("#"):
            raise ValueError("Every VCS test requirement needs a resolved wheel-cache identity")
        else:
            build_lines.append(line)
            install_lines.append(line)
    if found != 1:
        raise ValueError("Expected one moving csubst test requirement")
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "build-requirements.txt").write_text("\n".join(build_lines) + "\n")
    (directory / "install-requirements.txt").write_text("\n".join(install_lines) + "\n")
    (directory / "source-identity.json").write_text(json.dumps({
        "csubst_sha": csubst_sha,
        "requirements_sha256": hashlib.sha256(requirements.read_bytes()).hexdigest(),
        "constraints_sha256": hashlib.sha256(constraints.read_bytes()).hexdigest(),
    }, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", required=True, type=Path)
    parser.add_argument("--csubst-sha", required=True)
    args = parser.parse_args()
    prepare(args.directory, args.csubst_sha)


if __name__ == "__main__":
    main()
