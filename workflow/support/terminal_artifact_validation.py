"""Stream the native gene-family terminal checks without retaining whole tables."""
from __future__ import annotations

import csv
import hashlib

CHUNK_BYTES = 1024 * 1024
MAX_RECORD_BYTES = 16 * 1024 * 1024
MAX_ARTIFACT_BYTES = 4 * 1024 * 1024 * 1024


class StatisticsLines:
    """Hash exact UTF-8 bytes and bound a complete, possibly quoted CSV record."""

    def __init__(self, handle):
        self.handle = handle
        self.digest = hashlib.sha256()
        self.size = self.record_size = 0

    def __iter__(self):
        return self

    def __next__(self):
        raw = self.handle.readline(MAX_RECORD_BYTES - self.record_size + 1)
        if not raw:
            raise StopIteration
        self.size += len(raw)
        self.record_size += len(raw)
        if self.record_size > MAX_RECORD_BYTES or self.size > MAX_ARTIFACT_BYTES:
            raise ValueError("terminal statistics exceeds record or artifact bounds")
        self.digest.update(raw)
        return raw.decode("utf-8")


def verify_statistics(handle, branch):
    lines = StatisticsLines(handle)
    previous_limit = csv.field_size_limit()
    try:
        csv.field_size_limit(MAX_RECORD_BYTES)
        reader = csv.DictReader(lines, delimiter="\t", strict=True)
        fields = reader.fieldnames
        if not fields or len(set(fields)) != len(fields):
            raise ValueError("malformed terminal statistics")
        if branch and not {"branch_id", "node_name"}.issubset(fields):
            raise ValueError("missing or duplicated terminal branch identities")
        lines.record_size = 0
        count = 0
        identities = set()
        for row in reader:
            if None in row or None in row.values():
                raise ValueError("malformed terminal statistics")
            count += 1
            if branch:
                identity = row["branch_id"]
                if not identity or identity in identities:
                    raise ValueError("missing or duplicated terminal branch identities")
                identities.add(identity)
            elif count > 1:
                raise ValueError("terminal tree statistics must have one row")
            lines.record_size = 0
        if not count:
            raise ValueError("malformed terminal statistics")
        return lines.digest.hexdigest(), lines.size
    finally:
        csv.field_size_limit(previous_limit)


def verify_pdf(handle):
    digest = hashlib.sha256()
    size = 0
    prefix = tail = b""
    while raw := handle.read(CHUNK_BYTES):
        size += len(raw)
        if size > MAX_ARTIFACT_BYTES:
            raise ValueError("terminal artifact exceeds verification bounds")
        digest.update(raw)
        prefix = (prefix + raw[:5])[:5]
        tail = (tail + raw)[-1024:]
    if not prefix.startswith(b"%PDF-") or b"%%EOF" not in tail:
        raise ValueError("invalid terminal PDF framing")
    return digest.hexdigest(), size


def verify_terminal_profile(store, family):
    members = []
    try:
        for subdir, suffix in (("stat_branch", "_stat.branch.tsv"), ("stat_tree", "_stat.tree.tsv"),
                               ("tree_plot", "_tree_plot.pdf")):
            name = family + suffix
            artifact = store.source_artifact(subdir, name)
            if artifact is None:
                raise ValueError(f"missing terminal artifact: {subdir}/{name}")
            with store.open_binary(subdir, name) as handle:
                digest, size = (verify_pdf(handle) if subdir == "tree_plot"
                                else verify_statistics(handle, subdir == "stat_branch"))
            members.append({"path": artifact.logical_path, "sha256": digest,
                            "size_bytes": size, "mtime_ns": artifact.mtime_ns})
        return {"profile": "gene-evolution-terminal-v1", "state": "verified", "members": members}
    except (OSError, ValueError, csv.Error) as exc:
        return {"profile": "gene-evolution-terminal-v1", "state": "unverified", "detail": str(exc)}
