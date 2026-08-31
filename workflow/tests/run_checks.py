#!/usr/bin/env python3
"""Run the same declared validation commands in dev, CI, and published runtimes."""

import argparse
import json
import os
import shlex
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
MANIFEST_PATH = Path(__file__).with_name("validation_manifest.json")
SUITES = ("smoke", "fast", "static", "integration-download", "integration-workflow", "runtime", "full", "r")


def load_manifest():
    manifest = json.loads(MANIFEST_PATH.read_text())
    if manifest.get("version") != 1:
        raise ValueError("Unsupported validation manifest version")
    return manifest


def commands_for(suite, workers, extra):
    manifest = load_manifest()
    commands = []
    if suite != "r":
        command = [sys.executable, "-m", "pytest", "-q"]
        if suite in {"fast", "integration-download", "integration-workflow", "full"}:
            command += ["-n", workers, "--dist", "load"]
        if suite != "full":
            # This option is defined in the test-directory conftest. Keeping
            # its value attached prevents pytest's initial parse from treating
            # it as a path before that conftest has been loaded via testpaths.
            command += [f"--gg-suite={suite}"]
        if suite in {"runtime", "full"}:
            command += ["--gg-strict-runtime"]
        # pytest's configured testpaths supplies the default; explicit file
        # paths and -k/-x/--lf retain their native meaning without a second root.
        commands.append(command + extra)
    if suite == "runtime":
        commands.append([sys.executable, "-m", "pytest", "-q", "--gg-strict-runtime"]
                        + manifest["runtime_extra_python"] + extra)
    if suite in {"runtime", "full", "r"}:
        commands.extend(manifest["r_commands"])
    return commands


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("suite", choices=SUITES)
    parser.add_argument("--workers", default=os.environ.get("GG_PYTEST_WORKERS", "2"))
    parser.add_argument("--list", action="store_true", help="Print commands without running them")
    args, extra = parser.parse_known_args()
    if args.workers != "auto" and (not args.workers.isdigit() or int(args.workers) < 1):
        parser.error("--workers must be a positive integer or auto")
    if extra[:1] == ["--"]:
        extra = extra[1:]
    if args.suite == "r" and extra:
        parser.error("pytest arguments apply only to Python suites")
    commands = commands_for(args.suite, args.workers, extra)
    if args.list:
        print(json.dumps(commands, indent=2))
        return 0
    env = os.environ.copy()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    if args.suite in {"runtime", "full", "r"}:
        env.update(load_manifest()["environment"])
    for command in commands:
        print("[validation] " + shlex.join(command), flush=True)
        result = subprocess.run(command, cwd=REPO_ROOT, env=env, check=False)
        if result.returncode:
            return result.returncode
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
