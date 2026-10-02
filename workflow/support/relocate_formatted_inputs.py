#!/usr/bin/env python3
"""Relocate generated output references before publishing staged species inputs."""

import argparse
import csv
import json
from pathlib import Path


def relocate_metadata(mappings, summary):
    mappings = [(Path(source).resolve(), Path(target).resolve()) for source, target in mappings]

    def relocate(value):
        if isinstance(value, dict):
            return {key: relocate(item) for key, item in value.items()}
        if isinstance(value, list):
            return [relocate(item) for item in value]
        if isinstance(value, str):
            for source, target in mappings:
                if value == str(source) or value.startswith(str(source) + "/"):
                    return str(target) + value[len(str(source)):]
        return value

    # Only generated audits and summary fields carry formatted-output paths.
    # Sequence/annotation bytes and their content/stat fingerprints stay intact.
    for source, _target in mappings:
        for path in source.rglob("*.json"):
            payload = json.loads(path.read_text(encoding="utf-8"))
            path.write_text(json.dumps(relocate(payload), ensure_ascii=True, indent=2,
                                       sort_keys=True) + "\n", encoding="utf-8")
    summary = Path(summary)
    with summary.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames
        rows = [relocate(row) for row in reader]
    with summary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mapping", nargs=2, action="append", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()
    relocate_metadata(args.mapping, args.summary)


if __name__ == "__main__":
    main()
