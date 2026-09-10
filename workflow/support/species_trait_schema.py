"""Bind declared trait types to a table and select numeric analysis inputs."""

import argparse
import csv
import hashlib
import io
import json
import math
import os
import re
import tempfile
from pathlib import Path

TYPES = {"numeric", "binary", "categorical", "text"}
MISSING = {"", "NA", "NaN", "nan"}


def schema_path(table):
    return Path(str(table) + ".schema.json")


def schema_payload(table_payload, trait_types):
    if any(kind not in TYPES for kind in trait_types.values()):
        raise ValueError("Unknown declared trait type.")
    return (json.dumps({"schema_version": 1, "table_sha256": hashlib.sha256(table_payload).hexdigest(),
                        "traits": trait_types}, ensure_ascii=False, indent=2) + "\n").encode("utf-8")


def read_table(table):
    payload = Path(table).read_bytes()
    rows = list(csv.reader(io.StringIO(payload.decode("utf-8")), delimiter="\t", strict=True))
    if not rows or len(rows[0]) < 2:
        raise ValueError("Trait table requires a species column and at least one trait column.")
    header = rows[0]
    if any(not name.strip() or name != name.strip() or any(c in name for c in "\t\r\n") for name in header) or len(header) != len(set(header)):
        raise ValueError("Trait headers must be unique, non-empty and unpadded.")
    if any(len(row) != len(header) for row in rows[1:]):
        raise ValueError("Trait table rows must have the same number of fields as the header.")
    return payload, header, rows[1:]


def select_traits(table, requested="all", eligible=None):
    payload, header, rows = read_table(table)
    available = header[1:]
    kinds = dict.fromkeys(available, "unspecified")
    sidecar = schema_path(table)
    if sidecar.exists() or sidecar.is_symlink():
        schema = json.loads(sidecar.read_text(encoding="utf-8"))
        if not isinstance(schema, dict) or schema.get("schema_version") != 1:
            raise ValueError("Unsupported trait schema.")
        if schema.get("table_sha256") != hashlib.sha256(payload).hexdigest():
            raise ValueError("Trait schema does not match the table content; regenerate both files.")
        declared = schema.get("traits")
        if (not isinstance(declared, dict) or set(declared) != set(available)
                or any(not isinstance(kind, str) or kind not in TYPES for kind in declared.values())):
            raise ValueError("Trait schema must declare exactly every trait column with a supported type.")
        kinds = declared
    automatic = requested.lower() == "all"
    chosen = available if automatic else re.findall(r"[^,\s]+", requested)
    if not chosen or len(set(chosen)) != len(chosen) or not set(chosen) <= set(available):
        raise ValueError("Trait selection must contain unique existing trait columns.")
    report = []
    for column, trait in enumerate(available, 1):
        kind = kinds[trait]
        selected = trait in chosen
        reason = "" if selected else "not_requested"
        if selected and eligible is not None and trait not in eligible:
            selected, reason = False, "observation_contract_excluded"
        if selected and kind in {"text", "categorical"}:
            if not automatic:
                raise ValueError(f"Trait {trait} is declared {kind}; explicitly encode it in a numeric/binary trait column first.")
            selected, reason = False, "declared_" + kind
        if selected:
            for row in rows:
                value = row[column]
                if value in MISSING:
                    continue
                try:
                    numeric = float(value)
                except ValueError:
                    raise ValueError(f"Trait {trait} contains an invalid numeric value: {value!r}. "
                                     "Text/category columns require a declared trait schema.") from None
                if not math.isfinite(numeric):
                    raise ValueError(f"Trait {trait} contains a non-finite value.")
        report.append({"trait": trait, "value_type": kind, "status": "selected" if selected else "excluded", "reason": reason})
    if not any(row["status"] == "selected" for row in report):
        raise ValueError("No numeric/binary traits selected for analysis.")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--table", required=True)
    parser.add_argument("--trait", default="all")
    parser.add_argument("--report", required=True)
    args = parser.parse_args()
    report = select_traits(args.table, args.trait)
    target = Path(args.report)
    inputs = [Path(args.table), schema_path(args.table), Path(__file__)]
    if target.is_symlink() or (target.exists() and not target.is_file()):
        raise ValueError("Trait selection report must be a regular file.")
    if any(target.resolve() == path.resolve() or
           (target.exists() and path.exists() and target.samefile(path)) for path in inputs):
        raise ValueError("Trait selection report must not replace an input.")
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="", dir=target.parent,
                                         prefix=".trait-selection-", delete=False) as handle:
            temporary = Path(handle.name)
            writer = csv.DictWriter(handle, fieldnames=["trait", "value_type", "status", "reason"], delimiter="\t")
            writer.writeheader()
            writer.writerows(report)
        os.replace(temporary, target)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
