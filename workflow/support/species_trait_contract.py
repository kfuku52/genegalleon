"""Carry observation provenance through species-trait selection.

Ordinary user trait tables remain usable without a sidecar. GBIF-derived
observations require their hash-bound contract and an explicit column selection.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import tempfile
from contextlib import contextmanager
from pathlib import Path

import pandas
from gbif_observations import COMPLETE_STATES, METRIC_DEFINITIONS, file_sha256
from species_trait_schema import schema_path, schema_payload, select_traits

CONTRACT_VERSION = 1


@contextmanager
def trait_output_transaction(outputs, inputs=()):
    """Publish a trait bundle with rollback and metadata-last visibility.

    This input-generation I/O contract uses only the standard library; ordinary
    trait generation does not require an unrelated regression package.
    """
    originals = [Path(path) for path in outputs]
    targets = {path: path.resolve() for path in originals}
    sources = [Path(path).resolve() for path in inputs]
    if len(set(targets.values())) != len(originals):
        raise ValueError("Trait output paths must be distinct")
    for index, (original, target) in enumerate(targets.items()):
        if original.is_symlink() or (target.exists() and not target.is_file()):
            raise ValueError(f"Trait output must be a regular file, not a symlink/directory: {original}")
        if any(target.exists() and other.exists() and os.path.samefile(target, other) for other in list(targets.values())[:index]):
            raise ValueError("Trait output paths must be distinct, including hard links")
        if any(target == source or (target.exists() and source.exists() and os.path.samefile(target, source)) for source in sources):
            raise ValueError(f"Trait output must not overwrite an input: {original}")
    directories, staged, backups, installed = [], {}, {}, []
    preserve_backups = False
    try:
        for index, (original, target) in enumerate(targets.items()):
            target.parent.mkdir(parents=True, exist_ok=True)
            directory = Path(tempfile.mkdtemp(prefix=".trait-bundle-", dir=target.parent))
            directories.append(directory)
            staged[original] = directory / f"{index}.new"
        yield staged
        if any(not path.is_file() for path in staged.values()):
            raise ValueError("Trait bundle publication is missing a staged output")
        try:
            for original, target in targets.items():
                if target.exists():
                    backup = staged[original].with_suffix(".backup")
                    shutil.copy2(target, backup)
                    backups[original] = backup
                    os.chmod(staged[original], target.stat().st_mode & 0o777)
                os.replace(staged[original], target)
                installed.append(original)
        except BaseException:
            for original in reversed(installed):
                try:
                    if original in backups:
                        os.replace(backups[original], targets[original])
                    else:
                        targets[original].unlink(missing_ok=True)
                except OSError:
                    preserve_backups = True
            if preserve_backups:
                print("Trait bundle rollback needs filesystem recovery; backups retained in: " + ", ".join(map(str, directories)), flush=True)
            raise
    finally:
        if not preserve_backups:
            for directory in directories:
                shutil.rmtree(directory)


def metadata_path(path: Path) -> Path:
    return Path(str(path) + ".metadata.json")


def sidecar_paths(path: Path) -> dict[str, Path]:
    return {"metadata": metadata_path(path),
            "gbif_quality": Path(str(path) + ".gbif-quality.tsv"),
            "gbif_observations": Path(str(path) + ".gbif-observations.tsv")}


def read_contract(path: Path) -> dict:
    sidecar = metadata_path(path)
    if not sidecar.exists():
        return {}
    metadata = json.loads(sidecar.read_text())
    if metadata.get("schema_version") != CONTRACT_VERSION:
        raise ValueError(f"Unsupported species-trait metadata schema: {sidecar}")
    if metadata.get("table_sha256") != file_sha256(path):
        raise ValueError(f"Species-trait table differs from its metadata: {path}; regenerate the bundle")
    if not isinstance(metadata.get("traits"), dict):
        raise ValueError(f"Species-trait metadata lacks trait roles: {sidecar}")
    return metadata


def select_species_traits(path: Path | str, selection: str = "all", auxiliary: tuple[str, ...] = (),
                          allow_empty: bool = False) -> tuple[pandas.DataFrame, dict]:
    path = Path(path)
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        for line, row in enumerate(reader, 2):
            if len(row) != len(header):
                raise ValueError(f"Species-trait row {line} has the wrong number of fields: {path}")
    if not header or len(header) != len(set(header)) or any(not name for name in header):
        raise ValueError(f"Species-trait header must contain distinct nonempty column names: {path}")
    frame = pandas.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if frame.shape[1] < 1:
        raise ValueError(f"Empty species-trait table: {path}")
    metadata = read_contract(path)
    definitions = metadata.get("traits", {})
    identifier = frame.columns[0]
    available = [name for name in frame.columns[1:] if name not in auxiliary]
    automatic = str(selection).strip().lower() in {"", "all"}
    protected = {name for name in available if name.startswith("gbif_") or definitions.get(name, {}).get("role", "trait") != "trait"}
    requested = [name for name in available if name not in protected] if automatic else list(dict.fromkeys(re.split(r"[,\s]+", selection.strip())))
    absent = sorted(set(requested + list(auxiliary)) - set(frame.columns[1:]))
    if absent:
        raise ValueError("Species-trait columns not found: " + ", ".join(absent))
    if not requested and not allow_empty:
        raise ValueError("No automatically eligible species traits. Select GBIF observed columns explicitly and retain their .metadata.json sidecar.")
    audit = {"source": str(path.resolve()), "table_sha256": file_sha256(path),
             "selection": "automatic" if automatic else "explicit", "selected": requested,
             "excluded_columns": sorted(protected) if automatic else [], "masked_species": {},
             "traits": {name: definitions.get(name, {"role": "trait", "source": "user"}) for name in requested}}
    gbif_selected = [name for name in requested if name.startswith("gbif_") or definitions.get(name, {}).get("source") == "gbif"]
    if gbif_selected:
        if not metadata:
            raise ValueError("GBIF observations require hash-bound species-trait metadata; legacy GBIF quality/range columns cannot be analyzed directly")
        if frame[identifier].duplicated().any():
            raise ValueError("GBIF species-level observations cannot have duplicated species rows")
        quality = metadata.get("gbif", {}).get("quality", [])
        by_species = {row["species"]: row for row in quality}
        if len(by_species) != len(quality):
            raise ValueError("Duplicate species in GBIF quality metadata")
        for name in gbif_selected:
            definition = definitions.get(name, {})
            if definition.get("role") != "observation" or definition.get("source_column") not in METRIC_DEFINITIONS:
                raise ValueError(f"{name} is a GBIF quality/legacy/undefined field, not an eligible observation metric")
            masked = []
            for index, species in frame[identifier].items():
                record = by_species.get(species)
                if record is None:
                    raise ValueError(f"GBIF quality metadata is missing species {species!r}")
                if record.get("status") not in COMPLETE_STATES or record.get("analysis_eligible") is not True:
                    frame.at[index, name] = ""
                    masked.append({"species": species, "status": record.get("status", "unknown"),
                                   "reason": record.get("ineligibility_reason", record.get("termination_reason", "unknown"))})
            audit["masked_species"][name] = masked
        audit["interpretation"] = "Association with retained GBIF observations; not evidence of an association with true biological ranges or of adaptation"
        audit["gbif"] = metadata["gbif"]
    if metadata:
        audit["metadata_sha256"] = file_sha256(metadata_path(path))
    if automatic and protected:
        print("Excluded observation/quality columns from automatic trait selection: " + ", ".join(sorted(protected)), flush=True)
    return frame.loc[:, list(dict.fromkeys([identifier, *requested, *auxiliary]))].copy(), audit



def select_analysis_traits(path, selection="all"):
    """Apply observation eligibility before numeric response type selection."""
    frame, audit = select_species_traits(path, selection)
    report = select_traits(path, "all" if audit["selection"] == "automatic" else selection, eligible=audit["selected"])
    selected = [row["trait"] for row in report if row["status"] == "selected"]
    audit["selected"] = selected
    audit["type_selection"] = report
    audit["trait_types"] = {row["trait"]: "numeric" if row["value_type"] == "unspecified" else row["value_type"]
                            for row in report if row["status"] == "selected"}
    audit["traits"] = {name: audit["traits"][name] for name in selected}
    audit["masked_species"] = {name: value for name, value in audit["masked_species"].items() if name in selected}
    if schema_path(path).exists():
        audit["schema_sha256"] = file_sha256(schema_path(path))
    return frame.loc[:, [frame.columns[0], *selected]].copy(), audit


def select_foreground_traits(path: Path | str, selection: str = "all", allow_empty: bool = False):
    frame, audit = select_species_traits(path, selection, allow_empty=allow_empty)
    if any(definition.get("role") == "observation" for definition in audit["traits"].values()):
        raise ValueError("GBIF observation metrics are not foreground definitions; supply independently defined biological foreground traits")
    return frame, audit


def write_derived_contract(output: Path, audit: dict) -> None:
    """Attach the selected source meaning to an already written derived table."""
    columns = pandas.read_csv(output, sep="\t", nrows=0).columns
    metadata = {"schema_version": CONTRACT_VERSION, "table_sha256": file_sha256(output),
                "traits": {name: definition for name, definition in audit["traits"].items() if name in columns}, "derivation": {key: value for key, value in audit.items() if key not in {"gbif", "traits"}}}
    if "gbif" in audit:
        metadata["gbif"] = audit["gbif"]
    metadata_path(output).write_text(json.dumps(metadata, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    if "trait_types" in audit:
        schema_path(output).write_bytes(schema_payload(output.read_bytes(), {name: kind for name, kind in audit["trait_types"].items() if name in columns}))


def trait_bundle_payloads(frame: pandas.DataFrame, path: Path, definitions: dict, gbif_bundle: dict | None) -> dict[Path, bytes]:
    sidecars = sidecar_paths(path)
    metadata = {"schema_version": CONTRACT_VERSION, "traits": definitions}
    quality_rows, observation_rows = [], []
    if gbif_bundle is not None:
        metadata["gbif"] = {key: value for key, value in gbif_bundle.items() if key not in {"analysis_rows", "observations"}}
        quality_rows = gbif_bundle["quality"]
        observation_rows = gbif_bundle["observations"]
    # JSON holds nested counts/reasons without lossy flattening. The TSV is a
    # companion for inspection, not an implicit extra set of numeric traits.
    quality_frame = pandas.DataFrame.from_records([
        {key: json.dumps(value, sort_keys=True, ensure_ascii=False) if isinstance(value, (list, dict)) else value
         for key, value in row.items()} for row in quality_rows
    ]) if quality_rows else pandas.DataFrame(columns=["species", "status", "analysis_eligible", "retained_count"])
    observation_frame = pandas.DataFrame.from_records(observation_rows) if observation_rows else pandas.DataFrame(columns=["species", *METRIC_DEFINITIONS])
    payloads = {path: frame.to_csv(sep="\t", index=False).encode("utf-8"),
                sidecars["gbif_quality"]: quality_frame.to_csv(sep="\t", index=False).encode("utf-8"),
                sidecars["gbif_observations"]: observation_frame.to_csv(sep="\t", index=False).encode("utf-8")}
    metadata["table_sha256"] = hashlib.sha256(payloads[path]).hexdigest()
    metadata["sidecars"] = {key: {"filename": item.name, "sha256": hashlib.sha256(payloads[item]).hexdigest()}
                            for key, item in sidecars.items() if key != "metadata"}
    payloads[sidecars["metadata"]] = (json.dumps(metadata, indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")
    return payloads


def write_trait_bundle(frame: pandas.DataFrame, path: Path, definitions: dict, gbif_bundle: dict | None,
                       inputs: list[Path]) -> None:
    payloads = trait_bundle_payloads(frame, path, definitions, gbif_bundle)
    with trait_output_transaction(list(payloads), inputs) as staged:
        for target, payload in payloads.items():
            staged[target].write_bytes(payload)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--select", default="all")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--allow-empty", action="store_true")
    parser.add_argument("--numeric-only", action="store_true")
    parser.add_argument("--selection-report", type=Path)
    args = parser.parse_args(argv)
    if args.selection_report and not args.numeric_only:
        parser.error("--selection-report requires --numeric-only")
    if args.numeric_only:
        frame, audit = select_analysis_traits(args.input, args.select)
    else:
        frame, audit = select_species_traits(args.input, args.select, allow_empty=args.allow_empty)
    outputs = [args.output] + ([args.report] if args.report else [])
    if args.numeric_only:
        outputs.append(schema_path(args.output))
    if args.selection_report:
        outputs.append(args.selection_report)
    outputs.append(metadata_path(args.output))
    with trait_output_transaction(outputs, [args.input, metadata_path(args.input), schema_path(args.input)]) as staged:
        frame.to_csv(staged[args.output], sep="\t", index=False)
        derived = {"schema_version": CONTRACT_VERSION, "table_sha256": file_sha256(staged[args.output]), "traits": audit["traits"], "derivation": audit}
        if "gbif" in audit:
            derived["gbif"] = audit["gbif"]
        Path(staged[metadata_path(args.output)]).write_text(json.dumps(derived, indent=2, ensure_ascii=False) + "\n")
        if args.numeric_only:
            staged[schema_path(args.output)].write_bytes(schema_payload(staged[args.output].read_bytes(), audit["trait_types"]))
        if args.selection_report:
            pandas.DataFrame(audit["type_selection"]).to_csv(staged[args.selection_report], sep="\t", index=False)
        if args.report:
            Path(staged[args.report]).write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
