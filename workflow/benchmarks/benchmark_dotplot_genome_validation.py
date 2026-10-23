#!/usr/bin/env python3
"""Compare complete genome-length parsing with bounded synthetic FASTA fixtures.

Fixtures are generated with <= 1 MiB chunks outside timers. Workers are independent
processes, and the first alternating pair is discarded. Wall and CPU time cover
the complete chromosome_lengths call. Linux VmHWM reports the worker's own peak
RSS, avoiding inherited parent resource.ru_maxrss artifacts. Exact typed output or
exception hashes must agree for every pair. Run through run_in_runtime.sh.
"""
import argparse
import gzip
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
from types import SimpleNamespace

PROFILES = {
    "raw64": (64, False, False, False),
    "raw80": (80, False, False, False),
    "rawlong": (0, False, False, False),
    "gz64": (64, True, False, False),
    "gz80": (80, True, False, False),
    "gzlong": (0, True, False, False),
    "mixed80": (80, False, True, False),
    "mixedgz80": (80, True, True, False),
    "late_invalid80": (80, False, True, True),
    "lategz_invalid80": (80, True, True, True),
    "small80": (80, False, True, False),
}


def fixture(path, profile, bases):
    width, compressed, mixed, invalid = PROFILES[profile]
    pattern = b"ACGTURYSWKMBDHVNacgturyswkmbdhvn" if mixed else b"ACGTN"
    if compressed:
        handle = gzip.open(path, "wb", compresslevel=1)
    else:
        handle = path.open("wb")
    with handle:
        handle.write(b">chr1 benchmark assembly\n")
        if width:
            line = (pattern * (width // len(pattern) + 1))[:width] + b"\n"
            block = line * max(1, 1048576 // len(line))
            block_bases = (len(block) // len(line)) * width
            while bases >= block_bases:
                handle.write(block)
                bases -= block_bases
            full, remainder = divmod(bases, width)
            handle.write(line * full)
            if remainder:
                handle.write(line[:remainder] + b"\n")
        else:
            block = (pattern * (1048576 // len(pattern) + 1))[:1048576]
            while bases:
                count = min(bases, len(block))
                handle.write(block[:count])
                bases -= count
            handle.write(b"\n")
        if invalid:
            handle.write(b"!\n")


def rss_bytes():
    status = Path("/proc/self/status")
    if status.exists():
        for line in status.read_text().splitlines():
            if line.startswith("VmHWM:"):
                return int(line.split()[1]) * 1024, "linux_proc_VmHWM"
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return raw * (1 if sys.platform == "darwin" else 1024), "resource_ru_maxrss"


def typed(value):
    if isinstance(value, dict):
        return ("dict", [(typed(key), typed(item)) for key, item in value.items()])
    if isinstance(value, (tuple, list)):
        return (type(value).__name__, [typed(item) for item in value])
    return (type(value).__name__, value)


def worker(args):
    sys.path.insert(0, str(args.fallback_support))
    sys.path.insert(0, str(args.support_root))
    from pairwise_synteny_dotplot import chromosome_lengths
    genes = [SimpleNamespace(seqid="chr1", end=args.bases)]
    start_wall, start_cpu = time.perf_counter(), time.process_time()
    try:
        result = chromosome_lengths({"kind": "genome", "path": str(args.fixture)}, genes, 1)
        outcome = ("ok", typed(result))
    except Exception as error:
        outcome = ("error", type(error).__module__, type(error).__qualname__, str(error))
    wall, cpu = time.perf_counter() - start_wall, time.process_time() - start_cpu
    peak_rss, rss_source = rss_bytes()
    assert (outcome[0] == "error") == ("invalid" in args.profile), outcome
    print(json.dumps({"wall_seconds": wall, "cpu_seconds": cpu,
                      "peak_rss_bytes": peak_rss, "rss_source": rss_source,
                      "outcome_sha256": hashlib.sha256(json.dumps(outcome).encode()).hexdigest()}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-support", type=Path, required=True)
    parser.add_argument("--candidate-support", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--bases", type=int, default=32 * 1048576)
    parser.add_argument("--trials", type=int, default=5)
    parser.add_argument("--profiles", nargs="+", choices=PROFILES, default=list(PROFILES))
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--support-root", type=Path)
    parser.add_argument("--fallback-support", type=Path)
    parser.add_argument("--fixture", type=Path)
    parser.add_argument("--profile", choices=PROFILES)
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    assert 1 <= args.bases <= 128 * 1048576 and 1 <= args.trials <= 7
    roots = {"baseline": args.baseline_support.resolve(), "candidate": args.candidate_support.resolve()}
    report = {"python": sys.version, "platform": platform.platform(), "arguments": sys.argv,
              "implementations": {key: {"path": str(root), "sha256": hashlib.sha256(
                  (root / "pairwise_synteny_dotplot.py").read_bytes()).hexdigest()} for key, root in roots.items()},
              "profiles": {}}
    with tempfile.TemporaryDirectory(prefix="gg-genome-validator-benchmark-") as directory:
        for profile in args.profiles:
            bases = min(args.bases, 4096) if profile == "small80" else args.bases
            compressed = PROFILES[profile][1]
            path = Path(directory) / (profile + (".fa.gz" if compressed else ".fa"))
            fixture(path, profile, bases)
            samples = {key: [] for key in roots}
            reference = None
            for trial in range(args.trials + 1):
                order = list(roots) if trial % 2 == 0 else list(reversed(roots))
                for implementation in order:
                    raw = subprocess.check_output([
                        sys.executable, str(Path(__file__).resolve()), "--worker",
                        "--baseline-support", str(roots["baseline"]),
                        "--candidate-support", str(roots["candidate"]),
                        "--support-root", str(roots[implementation]), "--fallback-support", str(roots["baseline"]),
                        "--fixture", str(path), "--profile", profile, "--bases", str(bases)], text=True)
                    row = json.loads(raw)
                    if reference is None:
                        reference = row["outcome_sha256"]
                    else:
                        assert reference == row["outcome_sha256"], (profile, implementation, trial)
                    if trial:
                        samples[implementation].append(row)
                    print(profile, implementation, "warmup" if not trial else trial,
                          row["wall_seconds"], row["cpu_seconds"], row["peak_rss_bytes"], flush=True)
            report["profiles"][profile] = {"bases": bases, "input_bytes": path.stat().st_size,
                "samples": samples, "medians": {name: {metric: statistics.median(row[metric] for row in rows)
                for metric in ("wall_seconds", "cpu_seconds", "peak_rss_bytes")} for name, rows in samples.items()}}
            print(profile, "medians", json.dumps(report["profiles"][profile]["medians"]), flush=True)
            path.unlink()
    report["parent_peak_rss_bytes"], report["parent_rss_source"] = rss_bytes()
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
