"""Reproducible file/body-equivalent collector benchmark, without scientific inputs."""

import argparse
import hashlib
import json
import platform
import resource
import sqlite3
import sys
import tempfile
import time
import zlib
from pathlib import Path

parser = argparse.ArgumentParser(description="Benchmark lossless metadata collection in owned temporary storage.")
parser.add_argument("--support", type=Path, default=Path(__file__).resolve().parents[1] / "support")
parser.add_argument("--work-dir", type=Path, required=True)
parser.add_argument("--attempts", type=int, default=1000)
parser.add_argument("--layout", choices=("observations", "version-cache"), default="observations")
parser.add_argument("--inventory-only", action="store_true")
parser.add_argument("--run-read-latency-ms", type=float, default=0)
parser.add_argument("--source-read-latency-ms", type=float, default=0)
parser.add_argument("--keep-originals", action="store_true")
args = parser.parse_args()
if not 0 <= args.run_read_latency_ms <= 100 or not 0 <= args.source_read_latency_ms <= 100:
    parser.error("choose simulated read latencies from 0 to 100 ms")
if not 1 <= args.attempts <= 10000 or not args.work_dir.is_dir():
    parser.error("choose 1..10000 attempts and an existing owned work directory")
support, base, count = args.support.resolve(), args.work_dir.resolve(), args.attempts
sys.path.insert(0, str(support))
import metadata_store as m  # noqa: E402

original_read = m._read_source
if args.run_read_latency_ms or args.source_read_latency_ms:
    def delayed_read(path, **kwargs):
        delay = args.source_read_latency_ms
        if Path(path).name == "run.json":
            delay += args.run_read_latency_ms
        if delay:
            time.sleep(delay / 1000)
        return original_read(path, **kwargs)
    m._read_source = delayed_read

with tempfile.TemporaryDirectory(prefix="pack-bench-", dir=base) as raw:
    work = Path(raw)
    root = work / ("observations" if args.layout == "observations" else "version_inventory")
    before = {}
    root.mkdir()
    for i in range(count if args.layout == "observations" else 0):
        aid = f"{i:032x}"
        folder = root / aid
        folder.mkdir(parents=True)
        values = {
            "run.json": {
                "attempt_id": aid,
                "workflow": "gg_gene_evolution",
                "execution_state": "exited",
                "finished_at_ns": 2,
                "started_at_ns": 1,
            }
        }
        values.update(
            {
                "contract-" + str(j) + ".json": {
                    "attempt_id": aid,
                    "family_id": "OG0000001",
                    "step": str(j),
                    "exit_code": 0,
                }
                for j in range(4)
            }
        )
        for name, value in values.items():
            body = json.dumps(value, sort_keys=True).encode()
            path = folder / name
            path.write_bytes(body)
            before[str(path.relative_to(root))] = hashlib.sha256(body).hexdigest()
    if args.layout == "version-cache":
        # Build the legacy flat payload/lock layout without touching its controls.
        before = {}
        for number in range(count):
            key = f"{number:064x}"
            body = f"tool version {number}".encode()
            for suffix, value in [
                ("log", body),
                ("json", json.dumps({"schema_version": 1, "key": key, "sha256": hashlib.sha256(body).hexdigest()}).encode()),
            ]:
                path = root / (key + "." + suffix)
                path.write_bytes(value)
                before[path.name] = hashlib.sha256(value).hexdigest()
            gate = root / (key + ".log.lock") / "gate"
            gate.mkdir(parents=True)
            (gate / "owner.json").write_text("{}")
    records_values_sha256 = None
    start = time.monotonic()
    if args.inventory_only:
        selected = m._inventory(root, 500000)
        elapsed = time.monotonic() - start
        assert selected == sorted(before) and not (root / m.DATABASE).exists()
        after = {key: hashlib.sha256((root / key).read_bytes()).hexdigest() for key in selected}
        assert before == after
    else:
        result = m.pack(root, scratch=work, writers_stopped=True, delete_originals=not args.keep_originals, operation_id="collector-benchmark")
        elapsed = time.monotonic() - start
        with sqlite3.connect(root / m.DATABASE) as db:
            after = {
                key: hashlib.sha256(zlib.decompress(body)).hexdigest()
                for key, body in db.execute("SELECT path,body FROM files")
            }
            records_values_sha256 = hashlib.sha256(
                json.dumps(db.execute("SELECT * FROM records ORDER BY path").fetchall()).encode()
            ).hexdigest()
        assert before == after and result["removed_files"] == (0 if args.keep_originals else len(before))
        if args.keep_originals:
            assert all(hashlib.sha256((root / key).read_bytes()).hexdigest() == sha for key, sha in before.items())
    print(
        json.dumps(
            {
                "seconds": round(elapsed, 6),
                "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                * (1 if sys.platform == "darwin" else 1024),
                "python": platform.python_version(),
                "platform": platform.platform(),
                "records_sha256": hashlib.sha256(json.dumps(after, sort_keys=True).encode()).hexdigest(),
                "files": len(after),
                "records_values_sha256": records_values_sha256,
                "inventory_only": args.inventory_only,
                "simulated_run_read_latency_ms": args.run_read_latency_ms,
                "simulated_source_read_latency_ms": args.source_read_latency_ms,
                "delete_originals": not args.keep_originals,
            }
        )
    )

    for connection in getattr(m._LOCAL, "connections", {}).values():
        connection.close()
