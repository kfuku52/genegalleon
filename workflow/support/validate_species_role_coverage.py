#!/usr/bin/env python3
"""Check declared per-species output coverage without equating optional roles."""

import argparse
import csv
from pathlib import Path

from validate_required_species_outputs import OUTPUT_COLUMNS, OUTPUT_LABELS

EXTENSIONS = {
    "cds": (".fa", ".fas", ".fasta", ".fna"),
    "gff": (".gff", ".gff3", ".gtf"),
    "genome": (".fa", ".fas", ".fasta", ".fna"),
}


def output_files(directory, role):
    if not directory.is_dir():
        return {}
    return {path.name: path for path in directory.iterdir()
            if not path.name.startswith(".") and path.is_file()
            and path.name.removesuffix(".gz").lower().endswith(EXTENSIONS[role])}


def validate(summary, directories, required=()):
    with summary.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"species_prefix", *OUTPUT_COLUMNS.values()}.issubset(reader.fieldnames or ()):
            raise ValueError("Species summary lacks declared output coverage columns")
        rows = list(reader)
    if not rows:
        raise ValueError("Species summary is empty")
    observed = {role: output_files(directory, role) for role, directory in directories.items()}
    declared = {role: set() for role in directories}
    species = set()
    for row in rows:
        name = (row.get("species_prefix") or "").strip()
        if not name or name in species:
            raise ValueError("Missing or duplicate species in output coverage: " + name)
        species.add(name)
        for role, column in OUTPUT_COLUMNS.items():
            raw = (row.get(column) or "").strip()
            supplied = (row.get(role + "_input_path") or "").strip()
            if not raw:
                if role in required or supplied:
                    raise ValueError(f"Missing declared/required {OUTPUT_LABELS[role]} output for {name}")
                continue
            if role == "gff" and not (row.get("cds_output_path") or "").strip():
                raise ValueError("GFF output has no paired CDS for " + name)
            path = Path(raw)
            candidate = observed[role].get(path.name)
            if (candidate is None or not path.is_file() or not candidate.samefile(path)
                    or path.stat().st_size == 0 or path.name in declared[role]):
                raise ValueError(f"Missing, duplicate or misplaced {OUTPUT_LABELS[role]} output for {name}: {path}")
            declared[role].add(path.name)
    for role in directories:
        extra = set(observed[role]) - declared[role]
        if extra:
            raise ValueError(f"Undeclared formatted {OUTPUT_LABELS[role]} outputs: " + ", ".join(sorted(extra)))
    return {"species": len(species), **{role: len(paths) for role, paths in declared.items()}}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species-summary", type=Path, required=True)
    for role in OUTPUT_COLUMNS:
        parser.add_argument("--species-" + role + "-dir", type=Path, required=True)
        parser.add_argument("--require-" + role, action="store_true")
    args = parser.parse_args(argv)
    try:
        result = validate(args.species_summary,
                          {role: getattr(args, "species_" + role + "_dir") for role in OUTPUT_COLUMNS},
                          [role for role in OUTPUT_COLUMNS if getattr(args, "require_" + role)])
    except (OSError, ValueError, csv.Error) as error:
        parser.error(str(error))
    print("Verified declared species roles: " + ", ".join(f"{key}={value}" for key, value in result.items()))


if __name__ == "__main__":
    main()
