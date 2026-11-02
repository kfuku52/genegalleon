#!/usr/bin/env python3
"""Compare exact-output metadata extraction in fresh, sequential Linux workers.

Pass a saved pre-change amalgkit_metadata_accessions.py as --baseline. Synthetic
TSV generation, module imports and output hashing are outside the timed region.
Workers share identical input bytes; their order alternates across repeats.
"""

import argparse
import csv
import hashlib
import importlib.util
import json
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

CURRENT = Path(__file__).resolve().parents[1] / "support/amalgkit_metadata_accessions.py"


def measure(module_path, root, output):
    spec = importlib.util.spec_from_file_location("measured_metadata", module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    start = time.perf_counter()
    module.extract_transcriptomic_rows(root / "metadata.tsv", root / "runs.txt", output)
    elapsed = time.perf_counter() - start
    peak = next(int(line.split()[1]) / 1024 for line in Path("/proc/self/status").read_text().splitlines()
                if line.startswith("VmHWM:"))
    with output.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    return {"wall_seconds": elapsed, "peak_mib": peak, "output_sha256": digest}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, default=CURRENT)
    parser.add_argument("--rows", type=int, default=40000)
    parser.add_argument("--columns", type=int, default=24)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--worker", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--output", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(measure(args.candidate, args.worker, args.output)))
        return
    if args.rows < 1 or args.columns < 4 or args.repeats < 2:
        parser.error("Use positive rows, at least four columns and at least two repeats")
    with tempfile.TemporaryDirectory(prefix="gg-metadata-benchmark-") as directory:
        root = Path(directory)
        with (root / "metadata.tsv").open("w", newline="") as stream, (root / "runs.txt").open("w") as runs:
            writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
            writer.writerow(["run", "scientific_name", "lib_source", "lib_strategy"] +
                            [f"field_{n}" for n in range(args.columns - 4)])
            for n in range(args.rows):
                writer.writerow([f"SRR{n}", "Plant one", "TRANSCRIPTOMIC", "RNA-Seq"] +
                                [f"value_{n}_{c}_" + "x" * 64 for c in range(args.columns - 4)])
                runs.write(f"SRR{n}\n")
        results = {"baseline": [], "candidate": []}
        for repeat in range(args.repeats):
            order = ("baseline", "candidate") if repeat % 2 == 0 else ("candidate", "baseline")
            for variant in order:
                result = subprocess.run([sys.executable, "-B", str(Path(__file__).resolve()),
                                         "--baseline", str(args.baseline), "--candidate", str(getattr(args, variant)),
                                         "--worker", str(root), "--output", str(root / "result.tsv")],
                                        text=True, capture_output=True, check=True)
                results[variant].append(json.loads(result.stdout))
        hashes = {r["output_sha256"] for values in results.values() for r in values}
        if len(hashes) != 1:
            raise RuntimeError("Metadata extraction outputs differ")
        summary = {k: {metric: statistics.median(r[metric] for r in values)
                       for metric in ("wall_seconds", "peak_mib")} for k, values in results.items()}
        print(json.dumps({"rows": args.rows, "columns": args.columns,
                          "input_bytes": (root / "metadata.tsv").stat().st_size,
                          "baseline_source_sha256": hashlib.sha256(args.baseline.read_bytes()).hexdigest(),
                          "candidate_source_sha256": hashlib.sha256(args.candidate.read_bytes()).hexdigest(),
                          "output_sha256": hashes.pop(), "medians": summary, "runs": results}, indent=2))


if __name__ == "__main__":
    main()
