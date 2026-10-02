"""Sealed native staging evidence and invocation-local fresh-read fences."""

import hashlib
import json
import re
from pathlib import Path
from urllib.parse import unquote, urlparse

from input_generation_array_state import _stat_identity

REUSE_FIELDS = ("reuse_staged_plan", "reuse_staged_plan_sha256", "reuse_staged_workspace",
                "reuse_staged_task_index", "reuse_staged_receipt_sha256")
ROLES = ("cds", "gff", "gbff", "genome")


def identity(path):
    path = Path(path)
    if path.is_symlink() or not path.is_file() or path.stat().st_size == 0:
        raise ValueError("Missing, empty or unsafe staging reuse input: " + str(path))
    return _stat_identity(path.stat())


class StagedProofReader:
    """Memoize small sealed documents, with change fences, within one invocation."""

    def __init__(self):
        self.documents = {}
        self.identities = {}

    def check(self):
        for path, before in self.identities.items():
            if identity(path) != before:
                raise ValueError("Staging reuse evidence changed: " + path)

    def read(self, path, expected):
        if not re.fullmatch(r"[0-9a-f]{64}", expected):
            raise ValueError("Invalid staging reuse evidence SHA256")
        key = (str(path), expected)
        if key not in self.documents:
            before = identity(path)
            if str(path) in self.identities and self.identities[str(path)] != before:
                raise ValueError("Staging reuse evidence changed: " + str(path))
            data = Path(path).read_bytes()
            if identity(path) != before or hashlib.sha256(data).hexdigest() != expected:
                raise ValueError("Staging reuse evidence SHA256 mismatch: " + str(path))
            self.documents[key] = json.loads(data)
            self.identities[str(path)] = before
        if identity(path) != self.identities[str(path)]:
            raise ValueError("Staging reuse evidence changed: " + str(path))
        return self.documents[key]

    def resolve(self, task):
        row = task["manifest_row"]
        supplied = [bool(row.get(key)) for key in REUSE_FIELDS]
        if not any(supplied):
            return None
        if not all(supplied) or row.get("bind_local_sources") != "1":
            raise ValueError("Staging reuse requires all proof fields and bind_local_sources=1")
        plan_path = Path(row["reuse_staged_plan"])
        workspace = Path(row["reuse_staged_workspace"])
        if not plan_path.is_absolute() or not workspace.is_absolute():
            raise ValueError("Staging reuse evidence paths must be absolute")
        old_plan = self.read(plan_path, row["reuse_staged_plan_sha256"])
        try:
            index = int(row["reuse_staged_task_index"])
            if index < 1 or index > old_plan["task_count"]:
                raise ValueError("Invalid staging reuse task index")
            original = old_plan["tasks"][index - 1]
        except (KeyError, IndexError, TypeError) as exc:
            raise ValueError("Invalid staging reuse task plan") from exc
        if old_plan.get("version") != 2 or old_plan.get("download_mode") != "staged":
            raise ValueError("Staging reuse requires a native staged version-2 plan")
        receipt_path = Path(str(plan_path) + ".tasks") / (str(index) + ".json")
        receipt = self.read(receipt_path, row["reuse_staged_receipt_sha256"])
        if receipt.get("plan_sha256") != row["reuse_staged_plan_sha256"] or receipt.get("task_index") != index:
            raise ValueError("Staging reuse receipt belongs to another plan/task")
        actual = receipt["task"]
        for key in ("provider", "species_key", "species_prefix"):
            if actual.get(key) != task[key] or original.get(key) != task[key]:
                raise ValueError("Staging reuse species/provider mismatch: " + key)
        hashes = {}
        for role in ROLES:
            url = row.get(role + "_url", "") or ""
            previous = actual.get(role + "_path") or ""
            if bool(url) != bool(previous) or row.get(role + "_archive_member"):
                raise ValueError("Staging reuse roles do not match the receipt: " + role)
            if not url:
                continue
            parsed = urlparse(url)
            source = Path(unquote(parsed.path))
            previous_path = Path(previous)
            if previous_path.is_relative_to("/workspace"):
                previous_path = workspace / previous_path.relative_to("/workspace")
            if (parsed.scheme != "file" or parsed.netloc or parsed.query or parsed.fragment
                    or not source.is_absolute() or source != previous_path
                    or source.resolve() != source):
                raise ValueError("Staging reuse source path does not match the receipt: " + role)
            identity(source)
            expected = actual["input_sha256"].get(previous, "")
            if not re.fullmatch(r"[0-9a-f]{64}", expected):
                raise ValueError("Staging reuse receipt has no valid source SHA256: " + role)
            if str(source) in hashes and hashes[str(source)] != expected:
                raise ValueError("Conflicting staging reuse source SHA256")
            hashes[str(source)] = expected
        if not hashes:
            raise ValueError("Staging reuse receipt has no input roles")
        return {"plan_sha256": row["reuse_staged_plan_sha256"],
                "receipt_sha256": row["reuse_staged_receipt_sha256"],
                "task_index": index, "input_sha256": hashes}


class FreshReadFence:
    """A full-content check stays valid only while every file identity is unchanged."""

    def __init__(self, paths):
        self.identities = {str(path): identity(path) for path in paths}
        self.hashes = None

    def check(self):
        for path, before in self.identities.items():
            if identity(path) != before:
                raise ValueError("Input changed during staging reuse: " + path)

    def certify(self, observed, expected):
        self.check()
        if observed != expected or set(expected) != set(self.identities):
            raise ValueError("Staging reuse input SHA256 mismatch")
        self.hashes = dict(observed)

    def verified(self, paths):
        self.check()
        if self.hashes is None or set(map(str, paths)) != set(self.hashes):
            raise ValueError("Staging reuse fresh-read scope does not match input roles")
        return dict(self.hashes)
