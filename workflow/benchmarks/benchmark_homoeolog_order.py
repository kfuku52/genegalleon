#!/usr/bin/env python3
"""Compare complete homoeolog grouping against a saved support-directory snapshot.

Synthetic chromosome/scaffold fixtures and imports are outside the timer. Workers
are isolated processes; a warmup pair precedes alternating paired samples. Timings
cover the complete homoeolog_order call, including anchor reading; small profiles
use batches and report time per call. Linux VmHWM measures each worker's own peak
RSS. Typed result hashes must match.
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
from types import SimpleNamespace

SHAPES = {"normal": (24, 24, "groups"), "dense": (8, 8, "dense"),
          "sparse": (500, 500, "groups"), "unplaced": (500, 500, "diagonal"),
          "sparse_large": (4000, 4000, "groups"),
          "dense_two_queries": (500, 2, "dense"),
          "dense_two_targets": (2, 500, "dense")}
BATCHES = {"normal": 100, "dense": 20, "sparse": 1, "unplaced": 2,
           "sparse_large": 1, "dense_two_queries": 1, "dense_two_targets": 1}


def fixture(root, profile):
    targets, queries, topology = SHAPES[profile]
    selected = [[f"t{i}" for i in range(targets)], [f"q{i}" for i in range(queries)]]
    genomes = [[], []]
    path = root / "anchors"
    with path.open("w") as handle:
        handle.write("###\n")
        for i, target in enumerate(selected[0]):
            partners = selected[1] if topology == "dense" else (
                [selected[1][i % queries]] if topology == "diagonal" else
                [selected[1][(i // 2 * 2) % queries], selected[1][(i // 2 * 2 + 1) % queries]])
            for j, query in enumerate(partners):
                # Distinct IDs retain unique support and vary exact harmonic scores.
                for k in range(1 + (i + j) % 4):
                    target_gene, query_gene = f"T{i}_{j}_{k}", f"Q{i}_{j}_{k}"
                    genomes[0].append(SimpleNamespace(gene_id=target_gene, seqid=target))
                    genomes[1].append(SimpleNamespace(gene_id=query_gene, seqid=query))
                    handle.write(f"{target_gene} {query_gene} 10\n")
                    if k == 0:
                        handle.write(f"{target_gene} {query_gene} 99\n")
    return selected, genomes, path


def typed(value):
    if isinstance(value, dict):
        return ("dict", [(typed(key), typed(item)) for key, item in value.items()])
    if isinstance(value, (tuple, list)):
        return (type(value).__name__, [typed(item) for item in value])
    return (type(value).__name__, value)


def rss_bytes():
    status = Path("/proc/self/status")
    if status.exists():
        for line in status.read_text().splitlines():
            if line.startswith("VmHWM:"):
                return int(line.split()[1]) * 1024, "linux_proc_VmHWM"
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak * (1 if sys.platform == "darwin" else 1024), "resource_ru_maxrss"


def worker(args):
    sys.path.insert(0, str(args.baseline_support))
    sys.path.insert(0, str(args.support_root))
    from pairwise_synteny_dotplot import homoeolog_order
    selected, genomes, path = fixture(args.fixture_root, args.profile)
    start_wall, start_cpu = time.perf_counter(), time.process_time()
    for _ in range(BATCHES[args.profile]):
        result = homoeolog_order(selected, genomes, path)
    wall = (time.perf_counter() - start_wall) / BATCHES[args.profile]
    cpu = (time.process_time() - start_cpu) / BATCHES[args.profile]
    peak, source = rss_bytes()
    print(json.dumps({"wall_seconds": wall, "cpu_seconds": cpu, "peak_rss_bytes": peak,
                      "rss_source": source, "anchor_bytes": path.stat().st_size,
                      "calls_per_sample": BATCHES[args.profile],
                      "result_sha256": hashlib.sha256(json.dumps(typed(result)).encode()).hexdigest()}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-support", type=Path, required=True)
    parser.add_argument("--candidate-support", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--trials", type=int, default=5)
    parser.add_argument("--profiles", nargs="+", choices=SHAPES, default=list(SHAPES))
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--support-root", type=Path)
    parser.add_argument("--fixture-root", type=Path)
    parser.add_argument("--profile", choices=SHAPES)
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    roots = {"baseline": args.baseline_support.resolve()}
    if args.candidate_support:
        roots["candidate"] = args.candidate_support.resolve()
    assert 1 <= args.trials <= 7
    report = {"python": sys.version, "platform": platform.platform(), "arguments": sys.argv,
              "implementations": {key: {"root": str(root), "sha256": hashlib.sha256(
                  (root / "pairwise_synteny_dotplot.py").read_bytes()).hexdigest()} for key, root in roots.items()},
              "profiles": {}}
    with tempfile.TemporaryDirectory(prefix="gg-homoeolog-order-") as temporary:
        for profile in args.profiles:
            samples = {name: [] for name in roots}
            reference = None
            for trial in range(args.trials + 1):
                order = list(roots) if trial % 2 == 0 else list(reversed(roots))
                for name in order:
                    fixture_root = Path(temporary) / f"{profile}-{name}-{trial}"
                    fixture_root.mkdir()
                    raw = subprocess.check_output([
                        sys.executable, str(Path(__file__).resolve()), "--worker",
                        "--baseline-support", str(roots["baseline"]), "--support-root", str(roots[name]),
                        "--fixture-root", str(fixture_root), "--profile", profile], text=True)
                    row = json.loads(raw)
                    if reference is None:
                        reference = row["result_sha256"]
                    else:
                        assert reference == row["result_sha256"], (profile, name, trial)
                    if trial:
                        samples[name].append(row)
                    print(profile, name, "warmup" if not trial else trial,
                          row["wall_seconds"], row["cpu_seconds"], row["peak_rss_bytes"], flush=True)
            report["profiles"][profile] = {"samples": samples, "medians": {name: {
                metric: statistics.median(row[metric] for row in rows)
                for metric in ("wall_seconds", "cpu_seconds", "peak_rss_bytes")}
                for name, rows in samples.items()}}
            print(profile, "medians", json.dumps(report["profiles"][profile]["medians"]), flush=True)
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
