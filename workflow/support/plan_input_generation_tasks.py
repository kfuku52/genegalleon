#!/usr/bin/env python3

import argparse
import json
from pathlib import Path
import sys

import format_species_inputs as fsi


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Discover gg_input_generation species tasks and write a reusable task-plan JSON."
    )
    parser.add_argument(
        "--provider",
        choices=("all",) + fsi.PROVIDERS,
        required=True,
        help="Input provider type. Use 'all' with --input-dir pointing to a provider-root directory.",
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Provider input directory to scan for species tasks.",
    )
    parser.add_argument(
        "--outfile",
        required=True,
        help="Output JSON path for the discovered task plan.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Exit with error if task discovery reports any species-level errors.",
    )
    return parser


def serialize_task(task):
    payload = {}
    for key, value in task.items():
        if isinstance(value, Path):
            payload[key] = str(value)
        elif value is None:
            payload[key] = ""
        else:
            payload[key] = value
    return payload


def main():
    parser = build_arg_parser()
    args = parser.parse_args()

    try:
      provider_inputs = fsi.resolve_provider_inputs(args)
    except ValueError as exc:
        parser.error(str(exc))

    all_tasks = []
    all_warnings = []
    all_errors = []
    resolved_inputs = []

    for provider, input_dir in provider_inputs:
        resolved_inputs.append({"provider": provider, "input_dir": str(input_dir)})
        if not input_dir.exists() or not input_dir.is_dir():
            message = "[{}] input directory not found: {}".format(provider, input_dir)
            if args.provider == "all":
                all_warnings.append(message)
            else:
                all_errors.append(message)
            continue
        tasks, warnings, errors = fsi.discover_tasks(provider, input_dir)
        all_tasks.extend(tasks)
        all_warnings.extend(warnings)
        all_errors.extend(errors)

    for warning in all_warnings:
        sys.stderr.write("Warning: {}\n".format(warning))
    for error in all_errors:
        sys.stderr.write("Error: {}\n".format(error))

    if args.strict and all_errors:
        return 1
    if not all_tasks:
        sys.stderr.write("No species tasks were discovered.\n")
        return 1

    outfile = Path(args.outfile).expanduser().resolve()
    outfile.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "version": 1,
        "provider": args.provider,
        "input_dir": str(Path(args.input_dir).expanduser().resolve()),
        "provider_inputs": resolved_inputs,
        "task_count": len(all_tasks),
        "species": [task["species_prefix"] for task in all_tasks],
        "tasks": [serialize_task(task) for task in all_tasks],
    }
    with open(outfile, "wt", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=True, indent=2, sort_keys=True)

    print("Discovered {} species tasks -> {}".format(len(all_tasks), outfile))
    return 0


if __name__ == "__main__":
    sys.exit(main())
