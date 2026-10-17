"""Compare verified selection on the same synthetic publication.

Run through workflow/tests/run_in_runtime.sh with positional mode baseline or
selected. Both modes hash every publication member. Setup is excluded from
measurements; two warmups precede five fresh-process repetitions. This measures
targeted model inspection, not rescue throughput or cold operating-system I/O.
"""
import argparse
import hashlib
import json
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
import rescue_model_store as store  # noqa: E402


def digest_records(records):
    result = hashlib.sha256()
    count = 0
    for record in records:
        result.update(json.dumps(record, ensure_ascii=False, separators=(",", ":"),
                                 allow_nan=False).encode() + b"\n")
        count += 1
    return count, result.hexdigest()


def records(count):
    for i in range(count):
        body = i // 7
        sequence = "ATG" + "".join("ACGT"[int(c, 16) % 4] for c in
                                   hashlib.shake_256(str(body).encode()).hexdigest(1500)) + "TAA"
        yield {"model_id": f"model_{i}", "seqid": "chr1", "strand": "+",
               "cds": [[body * 4, body * 4 + len(sequence), 0]], "sequence": sequence,
               "status": "proposal", "support": [{"donor": f"species_{i % 23}",
               "query": f"gene_{i}"}], "identity": .9, "coverage": .95,
               "raw_prediction": {"seqid": "chr1", "strand": "+", "sequence": sequence,
                                  "cds": [[body * 4, body * 4 + len(sequence), 0]],
                                  "identity": .9, "coverage": .95}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("baseline", "selected"))
    parser.add_argument("--worker", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        mode, directory = args.mode, args.worker
        selected = {"model_0", "model_34999"}
        started = time.perf_counter()
        if mode == "baseline":
            iterator = (m for m in store.iter_models(directory) if m.get("model_id") in selected)
        else:
            iterator = store.iter_models(directory, model_ids=selected)
        count, sha = digest_records(iterator)
        print(json.dumps({"seconds": time.perf_counter() - started,
                          "maxrss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                          "count": count, "output_sha256": sha}))
        return
    mode = args.mode
    with tempfile.TemporaryDirectory(prefix="gg-selection-component-") as temporary:
        directory = Path(temporary)
        manifest = store.write_model_store(directory, records(35000), codec="gzip",
                                           shard_bytes=1024 * 1024)
        files = {p.relative_to(directory).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                 for p in directory.rglob("*") if p.is_file()}
        (directory / "receipt.json").write_text(json.dumps({"key": {"species": "synthetic"},
                                                           "files": files}))
        runs = []
        for index in range(7):
            run = json.loads(subprocess.check_output([sys.executable, str(Path(__file__).resolve()),
                                                     mode, "--worker", str(directory)], text=True))
            if index >= 2:
                runs.append(run)
        print(json.dumps({"mode": mode, "records": 35000, "bodies": manifest["counts"]["bodies"],
                          "input_bytes": sum(p.stat().st_size for p in directory.rglob("*") if p.is_file()),
                          "warmups": 2, "repeats": 5, "runs": runs,
                          "median_seconds": statistics.median(r["seconds"] for r in runs),
                          "median_maxrss_kib": statistics.median(r["maxrss_kib"] for r in runs),
                          "python": sys.version, "platform": platform.platform(),
                          "scope": "synthetic selection only; full member checksums remain enabled"}))


if __name__ == "__main__":
    main()
