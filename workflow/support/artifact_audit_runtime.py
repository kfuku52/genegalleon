"""Run-scoped digest reuse and bounded progress for full provenance audits."""
from __future__ import annotations

import concurrent.futures
import hashlib
import json
import math
import stat
import sys
import threading
import time
from pathlib import Path

from workflow_observation import atomic_json

PROGRESS_SCHEMA = "genegalleon-step-progress-v1"
RESULT_SCHEMA = "genegalleon-audit-result-v1"


def source_identity(support):
    hashes = {}
    for path in sorted(Path(support).glob("*.py")):
        before = signature(path)
        if len(hashes) >= 256 or before[3] > 4 * 1024 * 1024:
            raise ValueError("Audit source closure exceeds bound")
        hashes["workflow/support/" + path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        if signature(path) != before:
            raise ValueError("Audit source changed while identifying its runtime")
    return hashlib.sha256(json.dumps(hashes, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def validate_observation(value, attempt_id, *, result=False, started_at_ns=0):
    """Validate bounded optional evidence without granting workflow completion."""
    if (not isinstance(value, dict) or value.get("schema") != (RESULT_SCHEMA if result else PROGRESS_SCHEMA)
            or value.get("step") != "artifact_audit" or value.get("attempt_id") != attempt_id
            or not isinstance(value.get("source_sha256"), str)
            or len(value["source_sha256"]) != 64 or any(c not in "0123456789abcdef" for c in value["source_sha256"])
            or len(json.dumps(value, allow_nan=False).encode()) > 65536):
        raise ValueError("Invalid artifact audit observation identity")
    def integer(name, maximum=10**30):
        number = value.get(name)
        if type(number) is not int or not 0 <= number <= maximum:
            raise ValueError("Invalid audit observation counter: " + name)
        return number
    def duration(name, nullable=False):
        number = value.get(name)
        if nullable and number is None:
            return
        if type(number) not in (int, float) or not math.isfinite(number) or number < 0:
            raise ValueError("Invalid audit observation duration: " + name)
    duration("elapsed_seconds")
    metrics = value.get("metrics")
    if (not isinstance(metrics, dict) or len(metrics) > 32
            or any(type(number) is not int or number < 0 for number in metrics.values())):
        raise ValueError("Invalid audit observation metrics")
    if result:
        integer("exit_code", 255)
        checked = integer("checked")
        if integer("finished_at_ns") < started_at_ns:
            raise ValueError("Audit result predates its attempt")
        counts = value.get("status_counts")
        if (not isinstance(counts, dict) or len(counts) > 128
                or any(not isinstance(key, str) or len(key) > 64 or type(count) is not int or count < 0
                       for key, count in counts.items()) or sum(counts.values()) != checked):
            raise ValueError("Invalid audit result status counts")
        for name in ("inventory_sha256", "report_sha256"):
            if (not isinstance(value.get(name), str) or len(value[name]) != 64
                    or any(c not in "0123456789abcdef" for c in value[name])):
                raise ValueError("Invalid audit result digest")
    else:
        completed = integer("completed")
        if value.get("total") is not None and integer("total") < completed:
            raise ValueError("Audit progress exceeds its total")
        if (value.get("phase") not in {"inventory", "manifest_hash", "legacy_inventory", "branch_identity",
                                      "source_revalidation", "report"}
                or value.get("state") not in {"running", "completed", "failed", "interrupted"}
                or value.get("eta_scope") != "current-phase-only" or not 1 <= integer("workers", 64)):
            raise ValueError("Invalid audit progress phase or state")
        start, advance, observed = (integer(name) for name in
                                    ("started_at_ns", "last_advanced_at_ns", "observed_at_ns"))
        if not started_at_ns <= start <= advance <= observed:
            raise ValueError("Invalid audit progress time order")
        duration("phase_elapsed_seconds")
        duration("eta_seconds", nullable=True)
    return value


def signature(path):
    value = Path(path).lstat()
    if stat.S_ISLNK(value.st_mode):
        raise ValueError(f"Symlinked audit source is unsupported: {path}")
    return (value.st_dev, value.st_ino, value.st_mode, value.st_size,
            value.st_mtime_ns, value.st_ctime_ns)


def directory_signature(path):
    root = Path(path)
    return tuple((str(child.relative_to(root)), signature(child)) for child in sorted(root.rglob("*")))


class AuditDigests:
    """Reuse one hash per source, rehash unique sources before publication.

    File/ZIP identities are checked on every reuse; directory membership and all
    source identities are rechecked before accepting the complete audit. No
    digest from an earlier audit is accepted instead of reading a file here.
    """
    def __init__(self, limit=250000):
        self.limit = limit
        self.entries = {}
        self.guards = {}
        self.sources = {}
        self.lock = threading.RLock()
        self.stripes = [threading.RLock() for _ in range(256)]
        self.bytes_hashed = 0
        self.cache_hits = 0
        self.unique_sources = 0
        self.revalidation_reads = 0

    def read(self, path, compute, *, member=None, directory=False, constraint=None):
        path = Path(path).absolute()
        key = str(path), member, directory, constraint
        with self.stripes[hash(key) % len(self.stripes)]:
            before = signature(path)
            with self.lock:
                previous = self.entries.get(key)
                if previous is not None and previous[0] == before:
                    self.cache_hits += 1
                    return previous[1]
            members = directory_signature(path) if directory else None
            value = compute()
            if signature(path) != before or directory and directory_signature(path) != members:
                raise ValueError(f"Audit source changed while hashing: {path}")
            with self.lock:
                self.bytes_hashed += int(value[1])
                self.unique_sources += 1
                if len(self.entries) < self.limit or key in self.entries:
                    self.entries[key] = before, value, members, compute
            return value

    def guard(self, key, identity, reader):
        with self.lock:
            if key in self.guards and self.guards[key][0] != identity:
                raise ValueError("Audit source changed during the audit")
            if len(self.guards) < self.limit or key in self.guards:
                self.guards[key] = identity, reader
            else:
                raise ValueError("Audit source guard inventory exceeds its bounded capacity")

    def resolve(self, key, reader):
        with self.stripes[hash(key) % len(self.stripes)]:
            with self.lock:
                if key in self.sources:
                    return self.sources[key]
            value = reader()
            with self.lock:
                if len(self.sources) < self.limit:
                    self.sources[key] = value
            return value

    def validate(self, workers=1, progress=None):
        with self.lock:
            entries = list(self.entries.items())
            guards = list(self.guards.values())
        def validate_entry(item):
            (path, _member, directory, _constraint), (before, value, members, compute) = item
            if signature(path) != before or directory and directory_signature(path) != members:
                raise ValueError(f"Audit source changed before publication: {path}")
            # Content is the final fence even if filesystem attributes remain
            # unchanged or cached after a rapid same-size write.
            actual = compute()
            with self.lock:
                self.bytes_hashed += int(actual[1])
                self.revalidation_reads += 1
            if actual != value or signature(path) != before or directory and directory_signature(path) != members:
                raise ValueError(f"Audit source content changed before publication: {path}")
        def validate_guard(item):
            identity, reader = item
            if reader() != identity:
                raise ValueError("Logical audit source changed before publication")
        if progress:
            progress.phase("source_revalidation", len(entries) + len(guards))
        with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
            for values, check in ((entries, validate_entry), (guards, validate_guard)):
                for offset in range(0, len(values), workers * 8):
                    for _ in executor.map(check, values[offset:offset + workers * 8]):
                        if progress:
                            progress.advance(metrics=self.metrics())

    def metrics(self):
        with self.lock:
            return {"digest_bytes_read": self.bytes_hashed, "digest_cache_hits": self.cache_hits,
                    "unique_digest_reads": self.unique_sources, "cache_entries": len(self.entries),
                    "digest_revalidation_reads": self.revalidation_reads}


class AuditProgress:
    """One small atomic snapshot; stdout remains reserved for the serve protocol."""
    def __init__(self, path, *, interval=10, attempt_dir=None, workers=1, source_sha256=None):
        if not math.isfinite(interval) or interval < 0:
            raise ValueError("Progress interval must be finite and nonnegative")
        self.path = Path(path)
        self.attempt_dir = Path(attempt_dir) if attempt_dir else None
        self.interval = interval
        self.begin = time.monotonic()
        self.phase_begin = self.begin
        self.lock = threading.RLock()
        self.stop = threading.Event()
        self.thread = None
        self.write_ns = 0
        self.state = {"schema": PROGRESS_SCHEMA, "step": "artifact_audit", "phase": "inventory",
                      "state": "running", "completed": 0, "total": None,
                      "attempt_id": self.attempt_dir.name if self.attempt_dir else None,
                      "started_at_ns": time.time_ns(), "last_advanced_at_ns": time.time_ns(),
                      "workers": workers, "source_sha256": source_sha256, "metrics": {}}

    def __enter__(self):
        if self.interval > 0:
            self.publish()
            self.thread = threading.Thread(target=self._loop, daemon=True)
            self.thread.start()
        return self

    def _loop(self):
        while not self.stop.wait(self.interval):
            self.publish()

    def phase(self, name, total=None):
        with self.lock:
            self.phase_begin = time.monotonic()
            self.state.update(phase=name, completed=0, total=total, last_advanced_at_ns=time.time_ns())
        if self.interval > 0:
            self.publish()

    def advance(self, count=1, metrics=None):
        with self.lock:
            self.state["completed"] += count
            self.state["last_advanced_at_ns"] = time.time_ns()
            if metrics is not None:
                self.state["metrics"] = metrics

    def payload(self):
        with self.lock:
            value = json.loads(json.dumps(self.state))
            elapsed = time.monotonic() - self.begin
            phase_elapsed = time.monotonic() - self.phase_begin
        total, completed = value["total"], value["completed"]
        value.update(observed_at_ns=time.time_ns(), elapsed_seconds=elapsed,
                     phase_elapsed_seconds=phase_elapsed, eta_seconds=None, eta_scope="current-phase-only")
        if total and completed >= 20 and phase_elapsed >= 10 and completed < total:
            value["eta_seconds"] = phase_elapsed / completed * (total - completed)
        return value

    def publish(self):
        value = self.payload()
        began = time.monotonic_ns()
        try:
            atomic_json(self.path, value)
            if self.attempt_dir:
                atomic_json(self.attempt_dir / "progress-artifact_audit.json", value)
        except (OSError, ValueError) as exc:
            print(f"Artifact audit progress unavailable: {exc}", file=sys.stderr)
        finally:
            with self.lock:
                self.write_ns += time.monotonic_ns() - began

    def finish(self, *, status, rows, report, inventory_sha256, metrics):
        report = Path(report)
        with self.lock:
            self.state.update(state="completed" if status == 0 else "failed", metrics=metrics)
        self.stop.set()
        if self.thread:
            self.thread.join()
        if self.interval > 0:
            self.publish()
        metrics = {**metrics, "progress_write_ns": self.write_ns}
        result = {"schema": RESULT_SCHEMA, "step": "artifact_audit",
                  "attempt_id": self.state["attempt_id"], "source_sha256": self.state["source_sha256"],
                  "exit_code": status, "checked": len(rows), "inventory_sha256": inventory_sha256,
                  "report_sha256": hashlib.sha256(report.read_bytes()).hexdigest(),
                  "status_counts": {state: sum(row["status"] == state for row in rows)
                                    for state in sorted({row["status"] for row in rows})},
                  "finished_at_ns": time.time_ns(), "elapsed_seconds": time.monotonic() - self.begin,
                  "metrics": metrics, "scope": "artifact-audit-only; not workflow completion"}
        atomic_json(report.with_suffix(report.suffix + ".result.json"), result)
        if self.attempt_dir:
            atomic_json(self.attempt_dir / "result-artifact_audit.json", result)

    def __exit__(self, exc_type, _exc, _traceback):
        self.stop.set()
        if self.thread:
            self.thread.join()
        if exc_type:
            with self.lock:
                self.state["state"] = "interrupted"
            if self.interval > 0:
                self.publish()
