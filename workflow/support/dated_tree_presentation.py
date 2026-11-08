"""Publication panels for saved dated trees, parsed and validated by NWKIT.

Geological boundaries: ICS International Chronostratigraphic Chart 2026/06,
https://stratigraphy.org/ICSchart/ChronostratChart2026-06.pdf . Precambrian is
shown as a single broad interval. Backgrounds alternate between two light greys.
"""

import csv
import re
from pathlib import Path

STATUS = ("single", "duplicated", "fragmented", "missing")
STATUS_LABELS = ("Single-copy", "Duplicated", "Fragmented", "Missing")
STATUS_COLOURS = ("#000000", "#B22222", "#666666", "#CCCCCC")
GEOLOGICAL_DATA = Path(__file__).with_name("geological_periods.tsv")
GEOLOGICAL_SOURCE = "https://stratigraphy.org/ICSchart/ChronostratChart2026-06.pdf"


def species_key(value):
    return str(value).strip().replace(" ", "_")


def read_metadata_table(path, required=()):
    """Reject ambiguous TSVs before DictReader can silently discard columns."""
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if not fields or not all(fields) or len(fields) != len(set(fields)) or set(required) - set(fields):
            raise ValueError(f"Invalid metadata header in {path}: columns must be unique and include required fields.")
        rows = list(reader)
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError(f"Metadata row width does not match its header in {path}.")
    return rows


def read_busco(summary, species, results=None, prefix="busco_cds", *, input_paths=None):
    rows = read_metadata_table(summary, [prefix + "_" + status for status in (*STATUS, "total")])
    by_species = {}
    for row in rows:
        key = species_key(row.get("Species", row.get("species", "")))
        if not key or key in by_species:
            raise ValueError("Missing or duplicate species in BUSCO summary.")
        by_species[key] = row
    datasets, counts, sources = {}, {}, {}
    for name in species:
        key = species_key(name)
        if key not in by_species:
            raise ValueError(f"Missing BUSCO counts for {name}.")
        row = by_species[key]
        raw = [row.get(prefix + "_" + status, "") for status in STATUS]
        if not all(re.fullmatch(r"[0-9]+", value) for value in raw):
            raise ValueError(f"Invalid BUSCO counts for {name}.")
        values = list(map(int, raw))
        total = row.get(prefix + "_total", "")
        if not re.fullmatch(r"[0-9]+", total) or sum(values) != int(total) or int(total) <= 0:
            raise ValueError(f"BUSCO counts do not sum to the positive total for {name}.")
        counts[name] = values
        datasets[name] = row.get(prefix + "_lineage", "").strip()
        if datasets[name].upper() in {"NA", "N/A", "NAN", "NONE", "NULL"}:
            datasets[name] = ""
    if results is None and not all(datasets.values()):
        stem = "species_cds" if prefix == "busco_cds" else "species_genome"
        for parent in Path(summary).resolve().parents[:3]:
            for suffix in ("_busco_short", "_busco_full"):
                candidate = parent / (stem + suffix)
                if candidate.is_dir():
                    results = candidate
                    break
            if results is not None:
                break
    if results is not None:
        if not Path(results).is_dir():
            raise ValueError(f"BUSCO results directory does not exist: {results}.")
        paths = sorted(Path(results).glob("*busco.short.txt")) + sorted(Path(results).glob("*busco.full.tsv"))
        for path in paths:
            key = species_key(re.sub(r"\.busco\.(short\.txt|full\.tsv)$", "", path.name))
            matched = [name for name in species if species_key(name) == key]
            if not matched:
                continue
            if input_paths is not None:
                input_paths.append(path.resolve())
            with path.open(encoding="utf-8") as handle:
                header = "".join(line for _, line in zip(range(30), handle, strict=False))
            found = set(re.findall(r"lineage dataset is:[ \t]*(\S+)", header, flags=re.IGNORECASE))
            if len(found) > 1:
                raise ValueError(f"Conflicting BUSCO dataset metadata in {path}.")
            if found:
                name = matched[0]
                dataset = found.pop()
                if datasets[name] and datasets[name] != dataset:
                    raise ValueError(f"Conflicting BUSCO datasets for {name}.")
                datasets[name] = dataset
                sources[name] = str(path.resolve())
    known = {dataset for dataset in datasets.values() if dataset}
    if len(known) > 1:
        raise ValueError("Mixed BUSCO lineage datasets cannot share a completeness axis.")
    totals = {sum(values) for values in counts.values()}
    if len(totals) != 1:
        raise ValueError("BUSCO totals differ between species; a shared percentage axis would be misleading.")
    dataset = next(iter(known)) if known and all(datasets.values()) else None
    return counts, dataset, sources


from nwkit.dated_presentation import geological_intervals as geological_intervals
from nwkit.dated_presentation import interval_label as interval_label
from nwkit.dated_presentation import read_branch_annotations as read_branch_annotations


def render_dated_tree(infile, outfile, **kwargs):
    from nwkit.dated_presentation import render_dated_tree as render
    kwargs.setdefault("geological_data", GEOLOGICAL_DATA)
    return render(infile, outfile, panel_loader=read_busco, **kwargs)
