#!/usr/bin/env python3
"""Lint local composite steps with the same actionlint checks as workflow steps."""

import argparse
import json
import subprocess
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--actionlint", default="actionlint")
    parser.add_argument("--shellcheck", default="shellcheck")
    args = parser.parse_args()
    for path in sorted((REPO_ROOT / ".github/actions").glob("*/action.y*ml")):
        action = yaml.load(path.read_text(), Loader=yaml.BaseLoader)
        if action.get("runs", {}).get("using") != "composite":
            continue
        # actionlint accepts workflows, not action metadata. Give the actual
        # composite steps an equivalent typed inputs context without executing
        # them; this preserves expression and shell lint after extraction.
        inputs = {
            name: {"description": value.get("description", name), "type": "string",
                   "required": value.get("required") == "true", "default": value.get("default", "")}
            for name, value in action.get("inputs", {}).items()
        }
        workflow = {
            "name": action["name"],
            "on": {"workflow_dispatch": {"inputs": inputs}},
            "permissions": {"contents": "read"},
            "jobs": {"composite": {"runs-on": "ubuntu-24.04", "steps": action["runs"]["steps"]}},
        }
        result = subprocess.run(
            [args.actionlint, f"-shellcheck={args.shellcheck}", "-stdin-filename", str(path), "-"],
            input=json.dumps(workflow), text=True, cwd=REPO_ROOT, check=False,
        )
        if result.returncode:
            return result.returncode
        print(f"Composite action lint passed: {path.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
