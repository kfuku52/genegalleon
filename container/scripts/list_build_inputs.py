#!/usr/bin/env python3
"""List repository files consumed by Docker COPY, excluding generated context files."""

import argparse
import json
import shlex
import shutil
from pathlib import Path

# These exclusions also appear at the end of .dockerignore. Keep the Docker
# context and its identity independent of Python imports and editor caches.
IGNORED_DIRECTORIES = {"__pycache__", ".pytest_cache", ".ruff_cache"}
IGNORED_SUFFIXES = (".pyc", ".pyo", ".swp", ".swo", "~")


def generated(path: Path) -> bool:
    return any(
        part in IGNORED_DIRECTORIES or part == ".DS_Store" or part.endswith(IGNORED_SUFFIXES)
        for part in path.parts
    )


def build_inputs(root: Path) -> list[str]:
    inputs = {".dockerignore", "container/Dockerfile"}
    dockerfile = (root / "container/Dockerfile").read_text(encoding="utf-8")
    for instruction in dockerfile.replace("\\\n", " ").splitlines():
        parts = instruction.strip().split(None, 1)
        if parts and parts[0].upper() == "ADD":
            raise ValueError("ADD inputs need an explicit build-input identity; use COPY for local files")
        if not parts or parts[0].upper() != "COPY":
            continue
        arguments = parts[1].strip()
        from_stage = False
        while arguments.startswith("--"):
            option, arguments = arguments.split(None, 1)
            from_stage = from_stage or option.startswith("--from=")
        if from_stage:
            continue
        sources = json.loads(arguments) if arguments.startswith("[") else shlex.split(arguments)
        for source in sources[:-1]:
            if Path(source).is_absolute() or ".." in Path(source).parts:
                raise ValueError(f"Docker COPY source must remain in the context: {source}")
            matches = sorted(root.glob(source))
            if not matches:
                raise ValueError(f"Docker COPY input is missing: {source}")
            for match in matches:
                if match.is_symlink():
                    raise ValueError(f"Runtime inputs must not use symlinked COPY sources: {source}")
                paths = match.rglob("*") if match.is_dir() else [match]
                for path in paths:
                    relative = path.relative_to(root)
                    if generated(relative):
                        continue
                    if path.is_symlink():
                        raise ValueError(f"Runtime inputs must be regular files, not symlinks: {relative}")
                    if path.is_file():
                        if any(character in str(relative) for character in "\n\r\t"):
                            raise ValueError(f"Unsupported control character in build input: {relative!s}")
                        inputs.add(relative.as_posix())
    return sorted(inputs)


def stage_native_context(root: Path, destination: Path) -> None:
    """Give the native SIF builder the same source files as the Docker context."""
    for relative in build_inputs(root):
        path = Path(relative)
        if relative in {".dockerignore", "container/Dockerfile"}:
            continue
        if path.is_relative_to("container"):
            target = destination / path.relative_to("container")
        elif path.is_relative_to("workflow/support/treevis"):
            target = destination / "treevis" / path.relative_to("workflow/support/treevis")
        else:
            raise ValueError(f"Docker input needs a native context mapping: {relative}")
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(root / path, target)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--stage-native-context", type=Path,
                        help="Copy the declared inputs into a native SIF build staging directory")
    args = parser.parse_args()
    if args.stage_native_context:
        stage_native_context(args.repo_root, args.stage_native_context)
        return
    for path in build_inputs(args.repo_root):
        print(path)


if __name__ == "__main__":
    main()
