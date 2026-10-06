"""Explicit, frozen target-species overrides for rescue/refinement prediction."""
import csv
import math
from pathlib import Path

FIELDS = {"max_intron": int, "max_interval": int, "padding": int,
          "minimum_coverage": float, "minimum_identity": float,
          "min_support": int, "candidate_limit": int}


def read_profiles(path, species):
    if path is None:
        return {}
    result = {}
    with Path(path).open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = reader.fieldnames or []
        if (len(columns) != len(set(columns)) or "species" not in columns
                or set(columns) - {"species", *FIELDS}):
            raise ValueError("Species profiles require species and supported parameter columns")
        for row in reader:
            name = (row.get("species") or "").strip()
            if None in row or any(v is None for v in row.values()) or name not in species or name in result:
                raise ValueError("Unknown, duplicate or malformed species profile: " + name)
            values = {}
            for key, converter in FIELDS.items():
                if not (row.get(key) or "").strip():
                    continue
                value = converter(row[key])
                valid = (0 <= value <= 1 if key.startswith("minimum_") else value >= 0 if key == "padding" else value >= 1)
                if not math.isfinite(value) or not valid:
                    raise ValueError("Invalid species profile: " + name + " / " + key)
                values[key] = value
            result[name] = values
    return result


def parameters_for(request, species):
    return {**request["parameters"], **request.get("species_profiles", {}).get(species, {})}
