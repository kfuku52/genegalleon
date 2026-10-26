"""Compare fresh native release verification over a many-file runtime closure."""

import argparse
import contextlib
import hashlib
import importlib.util
import io
import json
import os
import resource
import statistics
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
BOOTSTRAP = REPO / "workflow/support/gg_entrypoint_bootstrap.sh"


def body(path):
    return path.read_text().split("<<'PY'\n", 1)[1].split("\nPY\n", 1)[0]


def canonical(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def worker(script, root, snapshot, input_delay_ms):
    code = compile(body(script), str(script), "exec")
    counts = {"stat_calls": 0, "lstat_calls": 0, "hashed_bytes": 0}
    counter_lock = threading.Lock()
    original_hash = hashlib.sha256

    def delivered(data):
        if data and input_delay_ms:
            time.sleep(input_delay_ms / 1000)
        with counter_lock:
            counts["hashed_bytes"] += len(data)

    class Hasher:
        def __init__(self, data=b"", **kwargs):
            delivered(data)
            self.wrapped = original_hash(data, **kwargs)

        def update(self, data):
            delivered(data)
            self.wrapped.update(data)

        def hexdigest(self):
            return self.wrapped.hexdigest()

    def counted(function, key):
        def call(*args, **kwargs):
            with counter_lock:
                counts[key] += 1
            return function(*args, **kwargs)
        return call

    os.stat = counted(os.stat, "stat_calls")
    os.lstat = counted(os.lstat, "lstat_calls")
    # Python 3.9 pathlib captures these functions in its accessor at import.
    accessor = getattr(Path(), "_accessor", None)
    if accessor is not None:
        type(accessor).stat = staticmethod(os.stat)
        type(accessor).lstat = staticmethod(os.lstat)
    hashlib.sha256 = Hasher
    sys.argv = ["-", str(root), "gg_input_generation", str(root / "workflow/support/gg_entrypoint_bootstrap.sh"), snapshot]
    output = io.StringIO()
    started = time.perf_counter()
    with contextlib.redirect_stdout(output):
        exec(code, {"__name__": "__main__"})
    elapsed = time.perf_counter() - started
    print(json.dumps({**counts, "seconds": elapsed,
                      "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 if sys.platform == "darwin" else 1),
                      "output_sha256": original_hash(output.getvalue().encode()).hexdigest()}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, default=BOOTSTRAP)
    parser.add_argument("--files", type=int, default=2000)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--input-delay-ms", type=float, default=0,
                        help="Simulated latency at hash input delivery; not a filesystem measurement")
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--root", type=Path)
    parser.add_argument("--snapshot")
    args = parser.parse_args()
    if args.worker:
        worker(args.candidate, args.root, args.snapshot, args.input_delay_ms)
        return
    if args.files < 1 or args.repeats < 2 or not 0 <= args.input_delay_ms <= 100:
        parser.error("positive file count and at least two repetitions required")
    fixture_path = REPO / "workflow/tests/test_gg_util_paths.py"
    spec = importlib.util.spec_from_file_location("release_fixture", fixture_path)
    fixture = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(fixture)
    with tempfile.TemporaryDirectory(prefix="runtime-release-benchmark-") as directory:
        root = Path(directory) / "release"
        fixture._write_runtime_release_contract(root, "gg_input_generation")
        manifest_path = root / ".kfauto-release.json"
        manifest = json.loads(manifest_path.read_text())
        for index in range(args.files):
            relative = f"workflow/support/fixture/group{index // 100:03d}/file{index:05d}.txt"
            target = root / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(b"runtime integrity fixture\n" * 40)
            manifest["installed_files"].append({"path": relative, "mode": target.stat().st_mode & 0o7777,
                                                "sha256": hashlib.sha256(target.read_bytes()).hexdigest()})
        manifest["contract_sha256"] = canonical({k: v for k, v in manifest.items() if k != "contract_sha256"})
        manifest_path.write_text(json.dumps(manifest, sort_keys=True, indent=2) + "\n")
        required = {"VERSION", "workflow/gg_input_generation_entrypoint.sh", "workflow/core/gg_input_generation_core.sh",
                    "workflow/support/artifact_provenance.py", "workflow/support/gg_entrypoint_bootstrap.sh", "workflow/support/gg_util.sh"}
        inventory = {row["path"]: row["sha256"] for row in manifest["installed_files"]}
        snapshot = canonical({"schema_version": 1, "kind": "gridengine-runtime-release-snapshot-v1",
                              "project_root": manifest["project_root"], "workflow": manifest["workflow"],
                              "release_root": manifest["release_root"], "target_ref": manifest["target_ref"],
                              "release_contract_sha256": manifest["contract_sha256"],
                              "file_sha256": {key: inventory[key] for key in sorted(required)}})
        results = {"baseline": [], "candidate": []}
        for iteration in range(args.repeats + 1):
            order = ["baseline", "candidate"] if iteration % 2 == 0 else ["candidate", "baseline"]
            for name in order:
                path = args.baseline if name == "baseline" else args.candidate
                command = [sys.executable, "-B", __file__, "--worker", "--baseline", str(args.baseline),
                           "--candidate", str(path), "--root", str(root), "--snapshot", snapshot,
                           "--input-delay-ms", str(args.input_delay_ms)]
                result = subprocess.run(command, capture_output=True, text=True)
                if result.returncode:
                    raise RuntimeError(result.stderr)
                if iteration:
                    results[name].append(json.loads(result.stdout))
        outputs = {row["output_sha256"] for values in results.values() for row in values}
        hash_bytes = {row["hashed_bytes"] for values in results.values() for row in values}
        if len(outputs) != 1 or len(hash_bytes) != 1:
            raise RuntimeError("verification output or fresh hashed bytes differ")
        report = {"files": len(manifest["installed_files"]), "repeats": args.repeats,
                  "simulated_input_delay_ms": args.input_delay_ms,
                  "scope": "fresh runtime verification only; no scientific inputs or scheduler submission",
                  "equal_output_and_hashed_bytes": True, "results": results,
                  "median_seconds": {name: statistics.median(row["seconds"] for row in values)
                                     for name, values in results.items()}}
        print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
