"""Bounded file-count benchmark; scientific jobs and real outputs are never used."""
import hashlib
import json
import resource
import statistics
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from workflow.support import coordination_lock_storage as storage
from workflow.support.shared_namespace_lock import namespace_lock


def trial(count):
    with tempfile.TemporaryDirectory() as temporary:
        base = Path(temporary)
        root = base / "output"
        root.mkdir()
        (root / ".gg_store").mkdir()
        science = root / "scientific.tsv"
        science.write_bytes(b"original scientific result\n")
        before_sha = hashlib.sha256(science.read_bytes()).hexdigest()
        for number in range(count):
            digest = hashlib.sha256(f"OG{number:07d}".encode()).hexdigest()
            for marker in (root / ".gg_run_locks" / f".family.{digest}.lock.guard",
                           root / ".gg_store_locks/family-locks-v2" / digest[:2] / f"{digest}.lock"):
                with namespace_lock(marker, exclusive=True):
                    pass
        journal = base / "journal"
        journal.mkdir(mode=0o700)
        before = sum(1 for _ in base.rglob("*"))
        expected = storage.snapshot(root)
        started = time.perf_counter()
        result = storage.retire(root, journal, expected, writers_stopped=True)
        elapsed = time.perf_counter() - started
        after = sum(1 for _ in base.rglob("*"))
        assert result["retired_entries"] == count * 6 and result["verified"]
        assert hashlib.sha256(science.read_bytes()).hexdigest() == before_sha
        return {"families": count, "entries_before": before, "entries_after": after,
                "removed_guard_entries": result["retired_entries"],
                "protected_science_sha256": before_sha, "wall_seconds": elapsed}


def main():
    trial(32)
    rows = [trial(1000) for _ in range(3)]
    assert len({(r["entries_before"], r["entries_after"], r["protected_science_sha256"])
                for r in rows}) == 1
    print(json.dumps({"trials": rows, "median_wall_seconds": statistics.median(
        r["wall_seconds"] for r in rows), "process_peak_rss_kib": resource.getrusage(
            resource.RUSAGE_SELF).ru_maxrss}, sort_keys=True))


if __name__ == "__main__":
    main()
