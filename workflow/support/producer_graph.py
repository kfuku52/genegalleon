"""Read-only producer dependencies from current declarations and recorded contracts.

This module identifies possible regeneration chains. It does not execute a stage,
infer whether a stage is enabled from historical output, or authorize a restart.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import artifact_provenance as provenance


class ProducerPlanError(ValueError):
    """A declared producer graph is incomplete or ambiguous."""


def _pairs(values: list[str], option: str) -> list[tuple[str, str]]:
    return provenance.parse_path_pairs(values, option)


def _recorded_entries(recorded: dict, collection: str) -> dict[str, dict]:
    rows = recorded.get(collection, [])
    if not isinstance(rows, list) or any(not isinstance(row, dict) or not row.get("label") for row in rows):
        raise ProducerPlanError(f"invalid recorded {collection}")
    result = {row["label"]: row for row in rows}
    if len(result) != len(rows):
        raise ProducerPlanError(f"duplicate recorded {collection} label")
    return result


def _declarations(args: argparse.Namespace) -> tuple[list[tuple[str, str, str]], list[tuple[str, str, str]]]:
    if (args.input_gene_family_store or args.input_gene_family_subdir or
            args.input_gene_family_artifact or args.output_fasta_type or args.recover_output):
        raise ProducerPlanError("this producer query does not cover store inputs, FASTA policy, or recovery recipes")
    inputs = [(label, raw, "file_or_directory") for label, raw in _pairs(args.input, "--input")]
    inputs += [(label, raw, "logical_directory") for label, raw in
               _pairs(args.input_logical_directory, "--input-logical-directory")]
    outputs = [(label, raw, "required") for label, raw in _pairs(args.output, "--output")]
    outputs += [(label, raw, "logical_directory") for label, raw in
                _pairs(args.output_logical_directory, "--output-logical-directory")]
    outputs += [(label, raw, "optional") for label, raw in
                _pairs(args.optional_output, "--optional-output")]
    if len({label for label, _, _ in inputs}) != len(inputs) or len({label for label, _, _ in outputs}) != len(outputs):
        raise ProducerPlanError("duplicate declared input or output label")
    if not outputs:
        raise ProducerPlanError("a producer must declare an output")
    return inputs, outputs


def _check_recorded(args: argparse.Namespace, recorded: dict,
                    inputs: list[tuple[str, str, str]],
                    outputs: list[tuple[str, str, str]]) -> dict[str, dict]:
    if (recorded.get("schema_version") != provenance.SCHEMA_VERSION or
            recorded.get("family_id") != args.family_id or recorded.get("step") != args.step or
            recorded.get("parameters") != provenance.normalized_parameters(args.parameter)):
        raise ProducerPlanError("recorded producer identity or parameters differ")
    input_rows = _recorded_entries(recorded, "inputs")
    output_rows = _recorded_entries(recorded, "outputs")
    optional_rows = _recorded_entries(recorded, "optional_outputs")
    expected_inputs = {label for label, _, _ in inputs}
    expected_outputs = {label for label, _, kind in outputs if kind != "optional"}
    expected_optional = {label for label, _, kind in outputs if kind == "optional"}
    if (set(input_rows) != expected_inputs or set(output_rows) != expected_outputs or
            set(optional_rows) != expected_optional):
        raise ProducerPlanError("recorded producer input/output set differs")
    for label, raw, kind in inputs:
        row = input_rows[label]
        reference = provenance.path_reference(Path(raw), args.logical_root.absolute(), args.workspace_root.absolute())
        if any(row.get(key) != value for key, value in reference.items()):
            raise ProducerPlanError(f"recorded input path differs: {label}")
        if kind == "logical_directory" and row.get("artifact_type") != "logical_directory":
            raise ProducerPlanError(f"recorded logical input type differs: {label}")
        if kind == "file_or_directory" and row.get("artifact_type") not in {"file", "directory"}:
            raise ProducerPlanError(f"unsupported recorded input type: {label}")
        if not isinstance(row.get("sha256"), str) or len(row["sha256"]) != 64:
            raise ProducerPlanError(f"recorded input digest is missing: {label}")
    for label, raw, kind in outputs:
        row = (optional_rows if kind == "optional" else output_rows)[label]
        reference = provenance.path_reference(Path(raw), args.logical_root.absolute(), args.workspace_root.absolute())
        if any(row.get(key) != value for key, value in reference.items()):
            raise ProducerPlanError(f"recorded output path differs: {label}")
        if kind == "logical_directory" and row.get("artifact_type") != "logical_directory":
            raise ProducerPlanError(f"recorded logical output type differs: {label}")
    return input_rows


def inspect(plan: dict, parse_contract) -> dict:
    """Return bounded dependency evidence; all paths are inspected without writes."""
    if not isinstance(plan, dict) or plan.get("schema") != "genegalleon-producer-plan-v1":
        raise ProducerPlanError("unsupported producer plan schema")
    raw_contracts, targets = plan.get("contracts"), plan.get("targets")
    if (not isinstance(raw_contracts, list) or not 1 <= len(raw_contracts) <= 1024 or
            not isinstance(targets, list) or not 1 <= len(targets) <= 1024):
        raise ProducerPlanError("producer plan requires 1..1024 contracts and targets")
    producers: dict[Path, list[dict]] = {}
    identities = set()
    workspace_root = None
    for item in raw_contracts:
        if not isinstance(item, dict) or type(item.get("enabled")) is not bool:
            raise ProducerPlanError("each producer requires argv and an explicit enabled boolean")
        args = parse_contract(item.get("argv"))
        if (not args.workspace_root.is_absolute() or not args.logical_root.is_absolute() or
                not args.manifest.is_absolute() or ".." in args.workspace_root.parts or
                ".." in args.manifest.parts):
            raise ProducerPlanError("producer roots and manifests must be absolute without traversal")
        if args.stale_policy != "rebuild":
            raise ProducerPlanError("producer candidates require the declared rebuild policy")
        if workspace_root is None:
            workspace_root = args.workspace_root.absolute()
        elif args.workspace_root.absolute() != workspace_root:
            raise ProducerPlanError("producer contracts must share one workspace")
        identity = (args.family_id, args.step)
        if identity in identities:
            raise ProducerPlanError("duplicate family/step producer")
        identities.add(identity)
        inputs, outputs = _declarations(args)
        if any(not Path(raw).is_absolute() or ".." in Path(raw).parts
               for _label, raw, _kind in [*inputs, *outputs]):
            raise ProducerPlanError("producer inputs and outputs must be absolute without traversal")
        if not provenance.declared_path_exists(args.manifest):
            raise ProducerPlanError(f"producer has no recorded manifest: {args.manifest}")
        recorded = provenance.load_manifest(args.manifest)
        input_rows = _check_recorded(args, recorded, inputs, outputs)
        node = {"family_id": args.family_id, "step": args.step, "enabled": item["enabled"],
                "inputs": inputs, "input_rows": input_rows, "outputs": outputs,
                "manifest": str(args.manifest), "recorded": recorded}
        for _label, raw, kind in outputs:
            path = Path(raw).absolute()
            if ".." in Path(raw).parts:
                raise ProducerPlanError(f"producer output contains traversal: {raw}")
            if not path.is_relative_to(args.workspace_root.absolute()):
                raise ProducerPlanError(f"producer output is outside its workspace: {path}")
            producers.setdefault(path, []).append({**node, "output_kind": kind})

    def visit(path: Path, active: set[Path]) -> dict:
        rows = producers.get(path, [])
        if len(rows) != 1:
            return {"path": str(path), "state": "ambiguous_producer" if rows else "missing_raw_input"}
        row = rows[0]
        if not row["enabled"]:
            return {"path": str(path), "state": "producer_disabled", "step": row["step"]}
        if path in active:
            return {"path": str(path), "state": "producer_cycle", "step": row["step"]}
        if row["output_kind"] == "optional":
            optional = True
        else:
            optional = False
        dependencies = []
        for label, raw, kind in row["inputs"]:
            source = Path(raw).absolute()
            try:
                if kind == "logical_directory":
                    current_sha, _size, _members = provenance.raw_or_zip_directory_digest(source)
                    current_type = "logical_directory"
                else:
                    current_sha, _size, current_type = provenance.sha256_path(source)
            except FileNotFoundError:
                dependencies.append(visit(source, active | {path}))
                continue
            except (OSError, ValueError, provenance.ProvenanceError) as exc:
                dependencies.append({"path": str(source), "state": "unavailable", "detail": str(exc)})
                continue
            expected = row["input_rows"][label]
            valid = current_sha == expected["sha256"] and current_type == expected["artifact_type"]
            dependencies.append({"path": str(source), "state": "verified_input" if valid else "changed_input"})
        # An optional producer may validly finish without publishing its file;
        # it cannot satisfy a later stage's required input by itself.
        allowed = {"verified_input", "producer_candidate"}
        state = ("optional_producer_candidate" if optional else "producer_candidate") if all(
            child["state"] in allowed for child in dependencies) else "blocked"
        return {"path": str(path), "state": state, "family_id": row["family_id"],
                "step": row["step"], "manifest": row["manifest"], "dependencies": dependencies}

    inspected = []
    for raw in targets:
        if not isinstance(raw, str) or not Path(raw).is_absolute():
            raise ProducerPlanError("targets must be absolute paths")
        if ".." in Path(raw).parts:
            raise ProducerPlanError("target contains traversal")
        target = Path(raw).absolute()
        if not target.is_relative_to(workspace_root):
            raise ProducerPlanError("target is outside the declared workspace")
        if target.is_symlink():
            inspected.append({"path": str(target), "state": "unavailable", "detail": "symlinked target"})
        elif target.exists():
            inspected.append({"path": str(target), "state": "already_present"})
        else:
            inspected.append(visit(target, set()))
    return {"targets": inspected, "coverage": "declared-recorded-producers-only",
            "execution_authorized": False, "requires_runtime_revalidation": True}
