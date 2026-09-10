#!/usr/bin/env python3
"""Extract BUSCO C (S+D) percentages and provenance from shared short summaries."""

import argparse
import csv
import math
import re
from pathlib import Path

from species_labeling import strip_species_label_terminal_suffixes

COLUMNS = ("species", "busco_complete_pct", "busco_duplicated_pct", "lineage", "mode", "busco_version", "source")


def read_short_summary(path):
    text = path.read_text(encoding="utf-8")
    values = re.findall(r"(?:^|[,\s])C:([0-9.]+)%\[S:([0-9.]+)%,D:([0-9.]+)%\]", text)
    if len(values) != 1:
        raise ValueError(f"Expected one BUSCO C/S/D summary in {path}")
    complete, single, duplicated = map(float, values[0])
    if not all(math.isfinite(x) and 0 <= x <= 100 for x in (complete, single, duplicated)) or duplicated > complete:
        raise ValueError(f"Invalid BUSCO percentages in {path}")
    # Allow the combined rounding uncertainty of the three displayed values.
    tolerance = sum(0.5 * 10 ** (-len(value.partition(".")[2])) for value in values[0])
    if abs(complete - single - duplicated) > tolerance + 1e-9:
        raise ValueError(f"Inconsistent BUSCO C versus S+D percentages in {path}")

    def field(pattern):
        match = re.search(pattern, text, re.IGNORECASE)
        return match.group(1).strip() if match else ""

    return {
        "species": strip_species_label_terminal_suffixes(path.name).replace(" ", "_"),
        "busco_complete_pct": complete,
        "busco_duplicated_pct": duplicated,
        "lineage": field(r"lineage dataset is:\s*(\S+)"),
        "mode": field(r"BUSCO was run in mode:\s*(\S+)"),
        "busco_version": field(r"BUSCO version is:\s*(\S+)"),
        "source": str(path.resolve()),
    }


def collect(directory):
    if directory.exists() and not directory.is_dir():
        raise ValueError(f"BUSCO summary directory is not a directory: {directory}")
    rows = []
    seen = set()
    # An absent shared BUSCO directory is explicitly represented by an empty table.
    for path in sorted(directory.glob("*busco.short.txt")):
        row = read_short_summary(path)
        if row["species"] in seen:
            raise ValueError(f"Ambiguous BUSCO summaries for {row['species']}")
        seen.add(row["species"])
        rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    rows = collect(args.directory)
    if any(args.output.resolve() == Path(row["source"]).resolve() or
           (args.output.exists() and args.output.samefile(row["source"])) for row in rows):
        raise ValueError("BUSCO metadata output must not replace a source summary")
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
