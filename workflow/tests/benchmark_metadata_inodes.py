"""Repeatable physical-entry benchmark with exact logical-body equivalence."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import statistics
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
import metadata_maintenance as maintenance
import metadata_store as metadata
import performance_log as performance


def entries(root):
    return sum(1 for _ in root.rglob("*"))


def digest(values):
    return hashlib.sha256(b"".join(sorted(values))).hexdigest()


def trial(base, count):
    observations = base / "observations"
    cache = base / "version_inventory"
    old = base / "baseline" / "versions" / "performance"
    new = base / "current" / "versions" / "performance"
    for root in (observations, cache, old, new):
        root.mkdir(parents=True)
    for number in range(count):
        name = f"{number:032x}"
        target = observations / name
        target.mkdir()
        value = {
            "schema": "genegalleon-observation-v1",
            "attempt_id": name,
            "execution_state": "exited",
            "finished_at_ns": 1,
        }
        (target / "run.json").write_text(json.dumps(value) + "\n")
        (cache / (f"{number:064x}" + ".log")).write_bytes(b"unchanged version body\n")
    metadata.pack(observations, writers_stopped=True, delete_originals=True)
    metadata.pack(cache, writers_stopped=True, delete_originals=True)
    for number in range(count):
        marker = cache / (f"{number:064x}" + ".log.lock")
        marker.write_bytes(b"")
        (cache / (marker.name + ".namespace-v1") / "readers").mkdir(parents=True)
    with metadata.read_snapshot():
        body_before = digest(metadata.archived_bytes(observations / f"{n:032x}/run.json") for n in range(count))
    cache_before = (cache / metadata.DATABASE).read_bytes()
    obs_before, cache_entries_before = entries(observations), entries(cache)
    frames = [(json.dumps({"id": n, "elapsed": n / 1000}, sort_keys=True) + "\n").encode() for n in range(count)]
    start = time.perf_counter()
    for pid, raw in enumerate(frames, 1):
        # Exact pre-change per-process append layout, using the same bytes.
        with (old / f"{pid}.jsonl").open("ab") as stream:
            stream.write(raw)
    baseline_seconds = time.perf_counter() - start
    start = time.perf_counter()
    for pid, raw in enumerate(frames, 1):
        performance.append(new, pid, raw)
    current_seconds = time.perf_counter() - start
    assert sorted(performance.records(old), key=lambda r: r["id"]) == sorted(
        performance.records(new), key=lambda r: r["id"]
    )
    live = entries(new)
    start = time.perf_counter()
    maintenance.maintain(
        observations,
        base / "observation-journal",
        kind="prune-observations",
        operation_id="benchmark-observations",
        writers_stopped=True,
    )
    maintenance.maintain(
        cache, base / "cache-journal", kind="retire-version-locks", operation_id="benchmark-cache", writers_stopped=True
    )
    metadata.pack(old, writers_stopped=True, delete_originals=True, operation_id="benchmark-performance")
    metadata.pack(new, writers_stopped=True, delete_originals=True, operation_id="benchmark-stripes")
    collection_seconds = time.perf_counter() - start
    with metadata.read_snapshot():
        assert digest(metadata.archived_bytes(observations / f"{n:032x}/run.json") for n in range(count)) == body_before
    assert (cache / metadata.DATABASE).read_bytes() == cache_before
    assert sorted(performance.records(old), key=lambda r: r["id"]) == [json.loads(raw) for raw in frames]
    assert sorted(performance.records(new), key=lambda r: r["id"]) == [json.loads(raw) for raw in frames]
    obs_after, cache_after = entries(observations), entries(cache)
    assert (
        maintenance.restore(observations, base / "observation-journal", writers_stopped=True)["restored_entries"]
        == count
    )
    assert maintenance.restore(cache, base / "cache-journal", writers_stopped=True)["restored_entries"] == count * 3
    # Count targets before restore; two journals add only six retained files.
    return {
        "records": count,
        "observations_before": obs_before,
        "observations_after": obs_after,
        "cache_before": cache_entries_before,
        "cache_after": cache_after,
        "performance_legacy_live": count,
        "performance_current_live": live,
        "performance_legacy_collected": entries(old),
        "performance_stripes_collected": entries(new),
        "baseline_append_seconds": baseline_seconds,
        "current_append_seconds": current_seconds,
        "collection_seconds": collection_seconds,
        "logical_bodies_equal": True,
        "restore_verified": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--count", type=int, default=2048)
    parser.add_argument("--trials", type=int, default=3)
    args = parser.parse_args()
    if not 64 <= args.count <= 10000 or not 1 <= args.trials <= 5:
        parser.error("bounded fixture requires 64..10000 records and 1..5 trials")
    values = []
    for _ in range(args.trials):
        with tempfile.TemporaryDirectory(prefix="gg-metadata-inodes-") as temporary:
            values.append(trial(Path(temporary), args.count))
    timing = {key: statistics.median(row[key] for row in values) for key in values[0] if key.endswith("seconds")}
    print(
        json.dumps(
            {"python": sys.version.split()[0], "platform": platform.platform(), "trials": values, "medians": timing},
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
