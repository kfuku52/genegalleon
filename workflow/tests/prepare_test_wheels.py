#!/usr/bin/env python3
"""Resolve the moving test dependencies only in temporary wheel-build inputs."""

import argparse
import hashlib
import json
import re
from pathlib import Path

TEST_DIR = Path(__file__).resolve().parent


def prepare(directory: Path, csubst_sha: str, nwkit_sha: str):
    sources = {"csubst": csubst_sha, "nwkit": nwkit_sha}
    for name, sha in sources.items():
        if not re.fullmatch(r"[0-9a-f]{40}", sha):
            raise ValueError(f"{name} source must be a resolved 40-character commit SHA")
    requirements = TEST_DIR / "requirements.txt"
    constraints = TEST_DIR / "requirements.lock.txt"
    build_lines = []
    install_lines = []
    found = dict.fromkeys(sources, 0)
    for line in requirements.read_text().splitlines():
        for name, sha in sources.items():
            prefix = f"{name} @ git+https://github.com/kfuku52/{name}.git@"
            if line.startswith(prefix):
                found[name] += 1
                build_lines.append(prefix + sha)
                install_lines.append(name)
                break
        else:
            if "git+" in line and not line.lstrip().startswith("#"):
                raise ValueError("Every VCS test requirement needs a resolved wheel-cache identity")
            build_lines.append(line)
            install_lines.append(line)
    for name, count in found.items():
        if count != 1:
            raise ValueError(f"Expected one moving {name} test requirement")
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "build-requirements.txt").write_text("\n".join(build_lines) + "\n")
    (directory / "install-requirements.txt").write_text("\n".join(install_lines) + "\n")
    (directory / "source-identity.json").write_text(json.dumps({
        **{f"{name}_sha": sha for name, sha in sources.items()},
        "requirements_sha256": hashlib.sha256(requirements.read_bytes()).hexdigest(),
        "constraints_sha256": hashlib.sha256(constraints.read_bytes()).hexdigest(),
    }, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", required=True, type=Path)
    parser.add_argument("--csubst-sha", required=True)
    parser.add_argument("--nwkit-sha", required=True)
    args = parser.parse_args()
    prepare(args.directory, args.csubst_sha, args.nwkit_sha)


if __name__ == "__main__":
    main()
