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
args = parser.parse_args()
if not 1 <= args.attempts <= 10000 or not args.work_dir.is_dir():
    parser.error("choose 1..10000 attempts and an existing owned work directory")
support, base, count = args.support.resolve(), args.work_dir.resolve(), args.attempts
sys.path.insert(0, str(support))
import metadata_store as m  # noqa: E402

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
    start = time.monotonic()
    result = m.pack(root, scratch=work, writers_stopped=True, delete_originals=True, operation_id="collector-benchmark")
    elapsed = time.monotonic() - start
    with sqlite3.connect(root / m.DATABASE) as db:
        after = {
            key: hashlib.sha256(zlib.decompress(body)).hexdigest()
            for key, body in db.execute("SELECT path,body FROM files")
        }
    assert before == after and result["removed_files"] == len(before)
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
            }
        )
    )

    for connection in getattr(m._LOCAL, "connections", {}).values():
        connection.close()
