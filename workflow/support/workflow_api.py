#!/usr/bin/env python3
"""Read-only GeneGalleon API. JSON on stdout; no scheduler or repair actions."""
from __future__ import annotations

import argparse
import base64
import contextlib
import csv
import hashlib
import io
import json
import os
import re
import sys
import time
from pathlib import Path

# Disable bytecode and the optional persistent digest cache before importing
# workflow libraries: even a query against a read-only workspace must work.
sys.dont_write_bytecode = True
os.environ.pop("GG_CONTENT_DIGEST_CACHE", None)
import artifact_provenance as provenance
from gene_family_output_store import GeneFamilyOutputStore, archive_queue_status, read_only_observation
from workflow_observation import contract_arguments, read_regular_json, strict_json_loads

SCHEMA = "genegalleon-api-v1"
MAX_JSON_BYTES = 8 * 1024 * 1024
STATUS_PAGE_BYTES = 12 * 1024 * 1024
STATUS_PAGE_SIZE = 512


class InvalidCursor(ValueError):
    """Only this error permits a client to discard a cursor and start again."""


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def read_json(path):
    return read_regular_json(path, MAX_JSON_BYTES)


def legacy_baseline(payload):
    diagnostics = payload.get("diagnostics", {})
    return any(str(diagnostics.get(key, "")).startswith("adopted_legacy") for key in (
        "provenance_state", "optional_output_provenance_state"))


def inventory(directory, *, pattern=None, directories=False, limit=10000):
    """Bound allocation and propagate I/O errors that pathlib glob suppresses."""
    selected = []
    with os.scandir(directory) as entries:
        for entry in entries:
            if entry.name.startswith("."):
                continue
            if directories and not (entry.is_dir(follow_symlinks=False) or entry.is_symlink()):
                continue
            if pattern is not None and not (entry.name.startswith(pattern) and entry.name.endswith(".json")):
                continue
            selected.append(Path(entry.path))
            if len(selected) > limit:
                raise ValueError("observation inventory exceeds configured bound")
    return sorted(selected)


def read_run(path):
    if Path(path).parent.is_symlink():
        raise ValueError("symlinked attempt directory")
    record = read_json(path)
    if not isinstance(record, dict) or record.get("schema") != "genegalleon-observation-v1":
        raise ValueError("unsupported run observation")
    attempt = record.get("attempt_id")
    codes = record.get("accepted_exit_codes")
    if (not isinstance(attempt, str) or re.fullmatch(r"[0-9a-f]{32}", attempt) is None
            or attempt != Path(path).parent.name or not isinstance(record.get("workflow"), str)
            or type(record.get("started_at_ns")) is not int
            or record["started_at_ns"] < 0
            or not isinstance(codes, list) or 0 not in codes
            or any(type(code) is not int or not 0 <= code <= 255 for code in codes)):
        raise ValueError("invalid run identity or execution contract")
    if record.get("execution_state") == "started":
        if (record.get("exit_code") is not None or "finished_at_ns" in record
                or record.get("execution_accepted") not in (None, False)):
            raise ValueError("started observation has contradictory completion fields")
    elif record.get("execution_state") == "exited":
        code = record.get("exit_code")
        if (type(code) is not int or not 0 <= code <= 255 or type(record.get("finished_at_ns")) is not int
                or record["finished_at_ns"] < record["started_at_ns"]
                or record.get("execution_accepted") is not (code in codes)):
            raise ValueError("invalid final execution observation")
    else:
        raise ValueError("unknown execution state")
    return record


def read_contract(path, attempt_id):
    row = read_json(path)
    if (not isinstance(row, dict) or row.get("schema") != "genegalleon-contract-observation-v1"
            or row.get("attempt_id") != attempt_id
            or any(not isinstance(row.get(key), str) or not row[key] for key in ("family_id", "step"))
            or row.get("operation") not in ("record", "needs-run")
            or type(row.get("exit_code")) is not int or not 0 <= row["exit_code"] <= 255
            or type(row.get("observed_at_ns")) is not int or row["observed_at_ns"] < 0
            or not isinstance(row.get("argv"), list) or any(not isinstance(value, str) for value in row["argv"])
            or (row.get("manifest_sha256") is not None and
                (not isinstance(row["manifest_sha256"], str) or re.fullmatch(r"[0-9a-f]{64}", row["manifest_sha256"]) is None))):
        raise ValueError("invalid contract observation")
    if Path(path).name != "contract-" + digest([row["family_id"], row["step"]]) + ".json":
        raise ValueError("contract filename does not match family/step identity")
    return row


def envelope(command, **values):
    return {"schema": SCHEMA, "command": command, "observed_at_ns": time.time_ns(),
            "read_only": True, **values}


def capabilities(_args):
    return envelope("capabilities", capabilities={
        "status": "attempt-delta-v1", "preflight": "declared-artifact-contracts-v1",
        "status_pages": "attempt-pages-v1",
        "verify": "family-provenance-v1", "runtime": "registered-config-v1",
        "errors": "owned-boundary-codes-v1",
    }, provenance_schema_versions=[provenance.SCHEMA_VERSION],
        verify_workspace_relocation="explicit-recorded-workspace-root-v1",
        verify_profiles=["gene-evolution-terminal-v1"],
        preflight_stale_policy_override=True,
        requires_kfauto=False, observations_opt_in="GG_OBSERVABILITY=1")


def status(args):
    if args.page_size is not None:
        return status_page(args)
    if args.page_cursor or args.known_records_file:
        raise ValueError("paged status options require --page-size")
    root = args.directory.absolute()
    if not root.is_dir():
        raise ValueError(f"observation directory does not exist: {root}")
    identity = str(root.resolve())
    previous = {}
    since = read_json(args.since_file) if args.since_file else args.since
    if since:
        if not isinstance(since, str) or len(since) > MAX_JSON_BYTES:
            raise ValueError("cursor exceeds size bound")
        cursor = json.loads(base64.b64decode(since.encode(), validate=True))
        if cursor.get("schema") != SCHEMA or cursor.get("root") != identity:
            raise ValueError("cursor belongs to another API version or directory")
        previous = cursor["records"]
        if not isinstance(previous, dict) or any(not isinstance(k, str) or not isinstance(v, str)
                                                for k, v in previous.items()):
            raise ValueError("invalid cursor records")
    records = {}
    changed = []
    errors = []
    # Per-attempt files avoid a shared NFS SQLite writer. A cursor records content
    # fingerprints, not wall-clock order, so late writes and clock skew are safe.
    for count, attempt_dir in enumerate(inventory(root, directories=True, limit=args.max_records), 1):
        path = attempt_dir / "run.json"
        if count > args.max_records:
            raise ValueError("observation inventory exceeds --max-records; narrow the directory")
        key = path.parent.name
        try:
            if path.parent.is_symlink():
                raise ValueError("symlinked attempt directory")
            record = attempt_observation(path.parent)
            records[key] = digest(record)
            if previous.get(key) != records[key]:
                changed.append(record)
        except (OSError, ValueError, TypeError, AttributeError, KeyError) as exc:
            errors.append({"attempt_id": key, "error_code": "observation_unavailable", "detail": str(exc)})
    cursor = base64.b64encode(json.dumps({"schema": SCHEMA, "root": identity, "records": records},
                                        separators=(",", ":")).encode()).decode()
    if len(cursor) > MAX_JSON_BYTES:
        raise ValueError("next cursor exceeds size bound")
    return envelope("status", records=changed, removed=sorted(previous.keys() - records.keys()
                    - {item["attempt_id"] for item in errors}), errors=errors,
                    next_cursor=cursor if not errors else None, coverage="recorded-attempts-only",
                    complete=not errors, snapshot_consistency="per-file atomic; not a global transaction",
                    inventory_count=len(records), missing_history="legacy runs have no observations")


def attempt_observation(directory):
    """One bounded attempt, with the same receipts and errors as legacy status."""
    record = read_run(directory / "run.json")
    errors, contracts = [], []
    for path in inventory(directory, pattern="error-"):
        error = read_json(path)
        if (not isinstance(error, dict) or error.get("schema") != "genegalleon-observation-v1"
                or not isinstance(error.get("error_code"), str)):
            raise ValueError("invalid error observation")
        errors.append(error)
    for path in inventory(directory, pattern="contract-"):
        receipt = read_contract(path, record["attempt_id"])
        contracts.append({**{name: receipt[name] for name in (
            "family_id", "step", "operation", "exit_code", "observed_at_ns", "manifest_sha256")},
            "receipt_sha256": digest(receipt), "workspace_root": receipt.get("workspace_root")})
    record.update(error_evidence=errors, contract_evidence=contracts)
    if len(json.dumps(record, sort_keys=True, allow_nan=False).encode()) > MAX_JSON_BYTES:
        raise ValueError("one attempt observation exceeds size bound")
    return record


def status_page(args):
    """Bounded pages and per-record deltas; cursors never contain the full inventory."""
    if (not 1 <= args.page_size <= STATUS_PAGE_SIZE or args.since or args.since_file
            or not 1 <= args.max_records <= 100000):
        raise ValueError("invalid paged status bounds or incompatible legacy cursor")
    root = args.directory.absolute()
    physical = root.resolve(strict=True)
    st = physical.stat()
    identity = digest([str(physical), st.st_dev, st.st_ino])
    directories = inventory(root, directories=True, limit=args.max_records)
    names = [path.name for path in directories]
    if any(re.fullmatch(r"[0-9a-f]{32}", name) is None for name in names):
        raise ValueError("invalid attempt directory identity")
    membership = digest(names)
    after = ""
    if args.page_cursor:
        try:
            if len(args.page_cursor) > 4096:
                raise ValueError("oversized cursor")
            token = strict_json_loads(base64.b64decode(args.page_cursor, validate=True))
            if (not isinstance(token, dict) or token.get("schema") != "genegalleon-status-page-v1"
                    or token.get("root_identity") != identity or token.get("inventory_sha256") != membership
                    or token.get("page_size") != args.page_size
                    or not isinstance(token.get("after"), str) or token["after"] not in names):
                raise ValueError("cursor identity or inventory changed")
            after = token["after"]
        except (ValueError, TypeError, KeyError, RecursionError) as exc:
            raise InvalidCursor("status cursor is invalid or its inventory changed") from exc
    known = {}
    if args.known_records_file:
        try:
            value = read_json(args.known_records_file)
            if (not isinstance(value, dict) or value.get("schema") != "genegalleon-status-known-v1"
                    or value.get("root_identity") != identity or not isinstance(value.get("records"), dict)
                    or len(value["records"]) > STATUS_PAGE_SIZE
                    or any(re.fullmatch(r"[0-9a-f]{32}", key) is None or not isinstance(sha, str)
                           or re.fullmatch(r"[0-9a-f]{64}", sha) is None for key, sha in value["records"].items())):
                raise ValueError("invalid known-record cursor")
            known = value["records"]
        except (ValueError, TypeError, KeyError, RecursionError) as exc:
            raise InvalidCursor("known-record cursor is invalid or belongs to another directory") from exc
    selected = [path for path in directories if path.name > after][:args.page_size]
    records, fingerprints = [], []
    size = 0
    for path in selected:
        row = attempt_observation(path)
        sha = digest(row)
        changed = known.get(path.name) != sha
        row_size = len(json.dumps(row, sort_keys=True, allow_nan=False).encode()) if changed else 0
        if fingerprints and size + row_size + 256 > STATUS_PAGE_BYTES:
            break
        size += row_size + 256
        fingerprints.append([path.name, sha])
        if changed:
            records.append(row)
    last = fingerprints[-1][0] if fingerprints else after
    more = bool(names and last != names[-1])
    token = {"schema": "genegalleon-status-page-v1", "root_identity": identity,
             "inventory_sha256": membership, "page_size": args.page_size, "after": last}
    cursor = base64.b64encode(json.dumps(token, separators=(",", ":")).encode()).decode() if more else None
    # Membership must remain unchanged while this page's records were read.
    if [path.name for path in inventory(root, directories=True, limit=args.max_records)] != names:
        raise InvalidCursor("status inventory changed during the page read")
    return envelope("status", status_protocol="attempt-pages-v1", complete=True, errors=[],
                    root_identity=identity, inventory_sha256=membership, inventory_count=len(names),
                    page_after=after, page_last=last, fingerprints=fingerprints, records=records,
                    next_cursor=cursor, snapshot_complete=not more,
                    snapshot_consistency="stable membership; per-file atomic, not a global content transaction")


def relocate_contract(args, workspace_root):
    old = args.workspace_root.absolute()
    new = workspace_root.absolute()
    if not new.is_dir():
        raise ValueError("replacement workspace root must exist")

    def moved(path):
        path = Path(path)
        return str(new / path.relative_to(old)) if path.is_absolute() and path.is_relative_to(old) else str(path)

    for name in ("manifest", "logical_root", "workspace_root"):
        setattr(args, name, Path(moved(getattr(args, name))))
    for name in ("input", "input_gene_family_store", "input_gene_family_subdir", "input_gene_family_artifact",
                 "input_logical_directory", "output", "output_logical_directory", "optional_output", "recover_output"):
        values = []
        for raw in getattr(args, name):
            label, value = provenance.parse_key_value(raw, "--" + name.replace("_", "-"))
            components = value.split("::") if name in {"input_gene_family_subdir", "input_gene_family_artifact"} else [value]
            components[0] = moved(components[0])
            values.append(label + "=" + "::".join(components))
        setattr(args, name, values)


def preflight_step(argv, workspace_root=None, stale_policy=None):
    if not isinstance(argv, list) or not all(isinstance(value, str) for value in argv):
        raise ValueError("each contract must be a list of provenance needs-run arguments")
    if any(value.split("=", 1)[0] in {"--dry-run", "--help", "-h"} for value in argv):
        raise ValueError("API controls dry-run; help is not a contract argument")
    parser = provenance.build_parser()
    try:
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            args = parser.parse_args(["needs-run", *argv, "--dry-run"])
    except SystemExit as exc:
        raise ValueError("invalid provenance contract arguments") from exc
    if workspace_root:
        relocate_contract(args, workspace_root)
    if stale_policy is not None:
        args.stale_policy = stale_policy
    if workspace_root or stale_policy is not None:
        argv = contract_arguments(args)
    with provenance.logical_observation(args):
        return inspect_preflight_step(args, argv)


def inspect_preflight_step(args, argv):
    provenance.configure_digest_cache(None)
    before = None
    if provenance.declared_path_exists(args.manifest):
        before = provenance.load_manifest(args.manifest)
        if not isinstance(before, dict) or not isinstance(before.get("diagnostics", {}), dict):
            raise ValueError("invalid manifest object or diagnostics")
    capture = io.StringIO()
    with contextlib.redirect_stdout(capture), contextlib.redirect_stderr(capture):
        result = provenance.dispatch(["needs-run", *argv, "--dry-run"])
    current = None
    try:
        current = provenance.build_contract(args, include_diagnostics=False)
    except (OSError, ValueError, provenance.ProvenanceError):
        pass
    matched = before is not None and current is not None and (
        provenance.contract_comparison_payload(before) == provenance.contract_comparison_payload(current))
    adopted = before is not None and legacy_baseline(before)
    if result == provenance.CURRENT:
        state = "verified_current" if matched and not adopted else (
            "legacy_reusable" if before is None or adopted else "policy_reuse_or_adoption")
        action = "reuse"
    elif result == provenance.NEEDS_RUN:
        state, action = "needs_run", "run_or_restore"
    elif result == provenance.STALE_STOP:
        state, action = "blocked", "resolve_stale_contract"
    else:
        state, action = "unavailable", "inspect_inputs_or_contract"
    # Bind the preview to exactly the declarations and observed contents. This
    # digest is evidence, not permission to skip the runtime validator later.
    after = provenance.load_manifest(args.manifest) if provenance.declared_path_exists(args.manifest) else None
    if before != after:
        state, action = "unavailable", "retry_observation"
    code = {provenance.STALE_STOP: "artifact_stale", provenance.ERROR: "provenance_error"}.get(result)
    if result == provenance.STALE_STOP and before is not None and current is not None:
        code = provenance.contract_difference_code(before, current)
    return {"family_id": args.family_id, "step": args.step, "state": state,
            "action": action, "runtime_exit_code": result, "stale_policy": args.stale_policy,
            "error_code": code,
            "binding_sha256": digest({"argv": argv, "recorded": before, "current": current}),
            "historical_generation_verified": state == "verified_current",
            "diagnostics": capture.getvalue().splitlines()}


def preflight(args):
    if args.attempt:
        run = read_run(args.attempt / "run.json")
        observed = []
        for count, path in enumerate(inventory(args.attempt, pattern="contract-"), 1):
            if count > 10000:
                raise ValueError("contract inventory exceeds bound")
            row = read_contract(path, run["attempt_id"])
            observed.append(row)
        plan = {"schema": "genegalleon-preflight-plan-v1",
                "contracts": [row["argv"] for row in sorted(observed, key=lambda row: row["observed_at_ns"])]}
    else:
        plan = read_json(args.plan)
    if not isinstance(plan, dict) or plan.get("schema") != "genegalleon-preflight-plan-v1":
        raise ValueError("unsupported preflight plan schema")
    contracts = plan.get("contracts")
    if not isinstance(contracts, list) or not 1 <= len(contracts) <= 10000:
        raise ValueError("plan requires 1..10000 declared contracts")
    results = [preflight_step(argv, args.workspace_root, args.stale_policy) for argv in contracts]
    keys = [(row["family_id"], row["step"]) for row in results]
    if len(set(keys)) != len(keys):
        raise ValueError("duplicate family/step contracts are ambiguous")
    return envelope("preflight", contracts=results, plan_sha256=digest(plan),
                    blocked=any(row["state"] in {"blocked", "unavailable"} for row in results),
                    coverage="declared-contracts-only", execution_authorized=False,
                    requires_runtime_revalidation=True,
                    stale_policy_override=args.stale_policy,
                    workspace_root_override=str(args.workspace_root) if args.workspace_root else None)


def verify_terminal_profile(store, family):
    members = []
    try:
        for subdir, suffix in (("stat_branch", "_stat.branch.tsv"), ("stat_tree", "_stat.tree.tsv"),
                               ("tree_plot", "_tree_plot.pdf")):
            name = family + suffix
            artifact = store.artifact(subdir, name)
            if artifact is None:
                raise ValueError(f"missing terminal artifact: {subdir}/{name}")
            with store.open_binary(subdir, name) as handle:
                raw = handle.read(16 * 1024 * 1024 + 1)
            if not 0 < len(raw) <= 16 * 1024 * 1024:
                raise ValueError("terminal artifact exceeds verification bounds")
            if subdir == "tree_plot":
                if not raw.startswith(b"%PDF-") or b"%%EOF" not in raw[-1024:]:
                    raise ValueError("invalid terminal PDF framing")
            else:
                reader = csv.DictReader(io.StringIO(raw.decode("utf-8")), delimiter="\t")
                rows = list(reader)
                fields = reader.fieldnames
                if (not rows or not fields or len(set(fields)) != len(fields)
                        or any(None in row or None in row.values() for row in rows)):
                    raise ValueError("malformed terminal statistics")
                if subdir == "stat_branch":
                    ids = [row.get("branch_id", "") for row in rows]
                    if not {"branch_id", "node_name"}.issubset(fields) or not all(ids) or len(set(ids)) != len(ids):
                        raise ValueError("missing or duplicated terminal branch identities")
                elif len(rows) != 1:
                    raise ValueError("terminal tree statistics must have one row")
            members.append({"path": artifact.logical_path, "sha256": hashlib.sha256(raw).hexdigest(),
                            "size_bytes": len(raw), "mtime_ns": artifact.mtime_ns})
        return {"profile": "gene-evolution-terminal-v1", "state": "verified", "members": members}
    except (OSError, ValueError, csv.Error) as exc:
        return {"profile": "gene-evolution-terminal-v1", "state": "unverified", "detail": str(exc)}


def verify(args):
    root, workspace = args.root.absolute(), args.workspace_root.absolute()
    if not root.is_dir() or not workspace.is_dir():
        raise ValueError("logical root and workspace must exist")
    if args.recorded_workspace_root and not args.attempt:
        raise ValueError("recorded workspace mapping requires --attempt")
    if args.profile and not {"summary_statistics", "tree_plot"}.issubset(args.require_step):
        raise ValueError("terminal profile requires summary_statistics and tree_plot contracts")
    store = GeneFamilyOutputStore(root, family_filter=args.family_id)
    family_before = store.family_observation(args.family_id)
    rows = []
    provenance.configure_digest_cache(None)
    if len(set(args.require_step)) != len(args.require_step):
        raise ValueError("duplicate required steps")
    manifest_names = dict(provenance.parse_unique_pairs(args.manifest, "--manifest"))
    if manifest_names.keys() - set(args.require_step):
        raise ValueError("manifest overrides must name a required step")
    for step in args.require_step:
        name = manifest_names.get(step, f"{args.family_id}.{step}.json")
        artifact = store.artifact(provenance.MANIFEST_SUBDIR, name)
        if artifact is None:
            rows.append({"step": step, "state": "unverified", "error_code": "evidence_missing",
                         "detail": "no matching manifest; this does not imply failure or require a rerun"})
            continue
        payload = provenance.read_manifest_from_store(store, artifact.name)
        if payload.get("family_id") != args.family_id or payload.get("step") != step:
            raise ValueError("manifest identity does not match the requested family/step")
        if not payload.get("outputs") and not payload.get("optional_outputs"):
            outcome, reason = "invalid_manifest", "no declared outputs"
        else:
            outcome, reason = provenance.audit_manifest(payload, store, root, workspace)
        diagnostics = payload.get("diagnostics", {})
        if not isinstance(diagnostics, dict):
            raise ValueError("invalid manifest diagnostics")
        adopted = legacy_baseline(payload)
        rows.append({"step": payload.get("step"), "state": (
            "legacy_reusable" if outcome == "current" and adopted else
            "verified_current" if outcome == "current" else "unverified"),
            "error_code": None if outcome == "current" else outcome, "detail": reason,
            "manifest": artifact.logical_path, "manifest_sha256": digest(payload),
            "parameters": payload.get("parameters"), "producer_diagnostics": payload.get("diagnostics", {})})
    family_after = GeneFamilyOutputStore(root, family_filter=args.family_id).family_observation(args.family_id)
    analytical_state = None if family_after is None else family_after["status"]
    all_current = all(row["state"] == "verified_current" for row in rows)
    # A family marked running/failed must not inherit completion from older files.
    completion = "verified_declared_steps" if all_current and family_before == family_after and analytical_state not in {"running", "failed"} else "unverified"
    attempt_id = None
    if args.attempt:
        run = read_run(args.attempt / "run.json")
        attempt_id = run["attempt_id"]
        if run.get("execution_state") != "exited" or run.get("execution_accepted") is not True:
            completion = "unverified"
        for row in rows:
            path = args.attempt / ("contract-" + digest([args.family_id, row["step"]]) + ".json")
            receipt = read_contract(path, attempt_id) if path.exists() else {}
            valid = receipt.get("schema") == "genegalleon-contract-observation-v1" and (
                receipt.get("attempt_id") == attempt_id and receipt.get("family_id") == args.family_id
                and receipt.get("step") == row["step"] and receipt.get("manifest_sha256") is not None
                and receipt.get("manifest_sha256") == row.get("manifest_sha256")
                and (receipt.get("operation"), receipt.get("exit_code")) in {("record", 0), ("needs-run", 1)})
            if valid:
                try:
                    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                        declared = provenance.build_parser().parse_args(["needs-run", *receipt["argv"], "--dry-run"])
                    if args.recorded_workspace_root:
                        if declared.workspace_root.absolute() != args.recorded_workspace_root.absolute():
                            raise ValueError("receipt does not belong to the explicitly mapped workspace")
                        relocate_contract(declared, workspace)
                    valid = (declared.logical_root.resolve() == root.resolve()
                             and declared.workspace_root.resolve() == workspace.resolve()
                             and declared.manifest.absolute() == (root / row["manifest"]).absolute()
                             and declared.family_id == args.family_id and declared.step == row["step"]
                             and type(receipt.get("observed_at_ns")) is int
                             and run["started_at_ns"] <= receipt["observed_at_ns"] <= run.get("finished_at_ns", -1))
                    if valid:
                        valid = preflight_step(contract_arguments(declared))["state"] == "verified_current"
                except (SystemExit, KeyError, OSError, ValueError, TypeError):
                    valid = False
            row["attempt_bound"] = valid
            if not valid:
                completion = "unverified"
    terminal = verify_terminal_profile(store, args.family_id) if args.profile else None
    if terminal is not None and terminal["state"] != "verified":
        completion = "unverified"
    final_family = GeneFamilyOutputStore(root, family_filter=args.family_id).family_observation(args.family_id)
    if final_family != family_after:
        completion = "unverified"
    family_after = final_family
    analytical_state = None if family_after is None else family_after["status"]
    return envelope("verify", family_id=args.family_id, completion_state=completion,
                    analytical_state=analytical_state, family_observation=family_after, contracts=rows,
                    coverage="required-steps-only; no inference about disabled or unlisted steps",
                    queue=archive_queue_status(root) if args.include_queue else None,
                    verification_scope="recorded inputs/outputs; use preflight for proposed settings",
                    attempt_id=attempt_id, historical_attempt_binding="requires --attempt with matching contract receipts",
                    terminal_validation=terminal,
                    recorded_workspace_root=str(args.recorded_workspace_root) if args.recorded_workspace_root else None,
                    validator_sha256=hashlib.sha256(Path(provenance.__file__).read_bytes()).hexdigest())


def runtime(args):
    record = read_run(args.attempt / "run.json")
    return envelope("runtime", attempt_id=record["attempt_id"], runtime=record.get("runtime"),
                    coverage="recorded-attempt-only")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("capabilities").set_defaults(handler=capabilities)
    status_parser = commands.add_parser("status")
    status_parser.add_argument("--directory", type=Path, required=True)
    cursor_source = status_parser.add_mutually_exclusive_group()
    cursor_source.add_argument("--since", help="opaque next_cursor from the last complete response")
    cursor_source.add_argument("--since-file", type=Path, help="JSON string containing next_cursor, for large inventories")
    status_parser.add_argument("--max-records", type=int, default=100000)
    status_parser.add_argument("--page-size", type=int, help="bounded status pages, at most 512 attempts per page")
    status_parser.add_argument("--page-cursor", help="opaque continuation from a previous page")
    status_parser.add_argument("--known-records-file", type=Path, help="bounded known-record fingerprints for this page")
    status_parser.set_defaults(handler=status)
    preflight_parser = commands.add_parser("preflight")
    preflight_source = preflight_parser.add_mutually_exclusive_group(required=True)
    preflight_source.add_argument("--plan", type=Path)
    preflight_source.add_argument("--attempt", type=Path)
    preflight_parser.add_argument("--workspace-root", type=Path, help="explicitly relocate workspace-owned declared paths")
    preflight_parser.add_argument("--stale-policy", choices=("stop", "rebuild", "reuse"),
                                  help="explicitly preview a proposed policy instead of the recorded policy")
    preflight_parser.set_defaults(handler=preflight)
    verify_parser = commands.add_parser("verify")
    verify_parser.add_argument("--root", type=Path, required=True)
    verify_parser.add_argument("--workspace-root", type=Path, required=True)
    verify_parser.add_argument("--family-id", required=True)
    verify_parser.add_argument("--require-step", action="append", required=True)
    verify_parser.add_argument("--manifest", action="append", default=[], metavar="STEP=FILENAME",
                               help="override a required step's FAMILY.STEP.json name within artifact_provenance")
    verify_parser.add_argument("--include-queue", action="store_true", help="also inventory the root-wide archive queue")
    verify_parser.add_argument("--attempt", type=Path, help="require matching successful-attempt contract receipts")
    verify_parser.add_argument("--recorded-workspace-root", type=Path,
                               help="explicit original workspace mapping for --attempt; relocate only its declared paths")
    verify_parser.add_argument("--profile", choices=("gene-evolution-terminal-v1",),
                               help="also validate the owned terminal statistics and PDF structure")
    verify_parser.set_defaults(handler=verify)
    runtime_parser = commands.add_parser("runtime")
    runtime_parser.add_argument("--attempt", type=Path, required=True)
    runtime_parser.set_defaults(handler=runtime)
    args = parser.parse_args(argv)
    try:
        with read_only_observation():
            result = args.handler(args)
        encoded = json.dumps(result, sort_keys=True, allow_nan=False)
    except Exception as exc:
        print(json.dumps(envelope(args.command, complete=False,
                                  error_code="cursor_invalid" if isinstance(exc, InvalidCursor) else "query_unavailable",
                                  detail=str(exc))))
        return 2
    print(encoded)
    return 0


if __name__ == "__main__":
    sys.exit(main())
