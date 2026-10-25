#!/usr/bin/env python3
"""Compare family-scoped WGD evidence loading in alternating fresh processes.

Fixtures are written by the parent before any worker starts. Operation timing
excludes imports, fixture creation and output hashing; peak RSS includes imports
and evidence loading. The baseline executes the original classify() loading
block extracted from its AST, then selects the same family scope in the timer.
"""

import argparse
import ast
import hashlib
import importlib.util
import json
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

SUPPORT = Path(__file__).resolve().parents[1] / "support"
PROFILES = ("sparse", "dense")


def load_source(path):
    # Saved baseline files need the same installed/current support dependencies.
    sys.path.insert(0, str(SUPPORT))
    spec = importlib.util.spec_from_file_location("_wgd_memory_source", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    if hasattr(module, "_load_family_evidence"):
        return module._load_family_evidence, False

    source = ast.parse(path.read_text(encoding="utf-8"))
    classify = next(node for node in source.body if isinstance(node, ast.FunctionDef) and node.name == "classify")
    start = next(i for i, node in enumerate(classify.body)
                 if isinstance(node, ast.Assign)
                 and any(isinstance(target, ast.Name) and target.id == "position_rows" for target in node.targets))
    stop = next(i for i in range(start, len(classify.body))
                if isinstance(classify.body[i], ast.For)
                and any(isinstance(node, ast.Constant) and node.value == "anchor_evidence.tsv"
                        for node in ast.walk(classify.body[i].iter)))
    wrapper = ast.parse("def _benchmark_load(evidence, gene_ids):\n    return positions, events, anchors\n")
    wrapper.body[0].body = classify.body[start:stop + 1] + wrapper.body[0].body
    namespace = dict(vars(module))
    exec(compile(ast.fix_missing_locations(wrapper), str(path), "exec"), namespace)
    return namespace["_benchmark_load"], True


def write_fixture(root, positions, anchors, targets):
    with (root / "gene_positions.tsv").open("w", encoding="utf-8", newline="") as handle:
        handle.write("gene_id\tspecies\tlocus_id\tseqid\trank\tstart\tend\n")
        for i in range(positions):
            handle.write(f"gene_{i:07d}\tSpecies_{i // 10000:03d}\tlocus_{i:07d}\tchr{i // 1000}\t"
                         f"{i % 1000 + 1}\t{i * 100}\t{i * 100 + 80}\n")
    with (root / "wgd_events.tsv").open("w", encoding="utf-8", newline="") as handle:
        handle.write("species_event_id\tevent_support\tnum_branch_matched_synteny_species\n")
        for i in range(32):
            handle.write(f"event_{i:02d}\t{'WGD-supported' if i % 2 else 'unresolved'}\t{i % 3}\n")
    with (root / "anchor_evidence.tsv").open("w", encoding="utf-8", newline="") as handle:
        handle.write("species\tblock_id\tgene_a\tgene_b\tks\tks_status\tspecies_event_id\tplacement_status\n")
        for i in range(anchors):
            a = i % positions
            base = a // 10000 * 10000
            width = min(10000, positions - base)
            b = base + (a * 17 + 11) % width
            if i % 1000 == 0:
                a, b = i // 1000 % targets, (i // 1000 + 1) % targets
            # Keep duplicate pairs, missing Ks and unmatched events: all are
            # observable evidence, including for tandem conflict detection.
            ks, status, placement = (("0.5", "ok", "interval_supported") if i % 5
                                     else ("NA", "cds_unavailable", "missing_ks"))
            handle.write(f"Species_{a // 10000:03d}\tblock_{i // 64:07d}\tgene_{a:07d}\tgene_{b:07d}\t"
                         f"{ks}\t{status}\tevent_{a % 32:02d}\t{placement}\n")


def fingerprint(positions, events, anchors):
    proof = hashlib.sha256()
    encoder = json.JSONEncoder(ensure_ascii=False, separators=(",", ":"))

    def append(value):
        for chunk in encoder.iterencode(value):
            proof.update(chunk.encode("utf-8"))
        proof.update(b"\n")

    for name, rows in (("positions", positions), ("events", events)):
        append(name)
        for key, row in rows.items():
            append([key, row])
    append("anchors")
    for pair, rows in anchors.items():
        append(pair)
        for row in rows:
            append(row)
    return proof.hexdigest()


def file_metadata(path):
    with path.open("rb") as handle:
        return {"bytes": path.stat().st_size, "sha256": hashlib.file_digest(handle, "sha256").hexdigest()}


def worker(args):
    loader, needs_scope = load_source(args.script)
    count = args.targets if args.profile == "sparse" else args.positions
    gene_ids = {f"gene_{i:07d}" for i in range(count)}
    started, started_cpu = time.perf_counter(), time.process_time()
    positions, events, anchors = loader(args.evidence, gene_ids)
    if needs_scope:
        positions = {key: row for key, row in positions.items() if key in gene_ids}
        anchors = {pair: rows for pair, rows in anchors.items() if all(gene in gene_ids for gene in pair)}
    seconds, cpu_seconds = time.perf_counter() - started, time.process_time() - started_cpu
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == "darwin" else 1024)
    result = {"seconds": seconds, "cpu_seconds": cpu_seconds, "peak_rss_bytes": peak,
              "positions": len(positions), "events": len(events), "anchor_pairs": len(anchors),
              "anchor_rows": sum(map(len, anchors.values())),
              "output_sha256": fingerprint(positions, events, anchors)}
    print(json.dumps(result))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--implementation", type=Path, default=SUPPORT / "wgd_ssd.py")
    parser.add_argument("--baseline-only", action="store_true")
    parser.add_argument("--positions", type=int, default=200000)
    parser.add_argument("--anchors", type=int, default=500000)
    parser.add_argument("--targets", type=int, default=100)
    parser.add_argument("--profiles", nargs="+", choices=PROFILES, default=list(PROFILES))
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--script", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--evidence", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--profile", choices=PROFILES, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if not 1 <= args.targets <= args.positions or args.anchors < 1 or not 1 <= args.repeats <= 20:
        parser.error("Require 1 <= targets <= positions, positive anchors and repeats in 1..20")
    if args.worker:
        worker(args)
        return
    if args.baseline is None or args.output is None:
        parser.error("--baseline and --output are required")
    scripts = {"baseline": args.baseline.resolve()}
    if not args.baseline_only:
        scripts["implementation"] = args.implementation.resolve()
    results = {"python": sys.version, "platform": platform.platform(), "arguments": sys.argv,
               "warmups": 1, "repeats": args.repeats, "positions": args.positions,
               "anchors": args.anchors, "sparse_targets": args.targets,
               "cache": "warm OS caches; fresh process per sample; fixtures prepared by parent",
               "scripts": {name: {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
                           for name, path in scripts.items()}, "samples": {}}
    with tempfile.TemporaryDirectory(prefix="gg-wgd-evidence-memory-") as temporary:
        evidence = Path(temporary)
        write_fixture(evidence, args.positions, args.anchors, args.targets)
        results["inputs"] = {path.name: file_metadata(path) for path in sorted(evidence.iterdir())}
        for profile in args.profiles:
            samples = {name: [] for name in scripts}
            reference = None
            for trial in range(args.repeats + 1):
                order = list(scripts) if trial % 2 == 0 else list(reversed(scripts))
                for name in order:
                    command = [sys.executable, str(Path(__file__).resolve()), "--worker", "--script", str(scripts[name]),
                               "--evidence", str(evidence), "--profile", profile, "--positions", str(args.positions),
                               "--anchors", str(args.anchors), "--targets", str(args.targets)]
                    sample = json.loads(subprocess.check_output(command, text=True))
                    proof = {key: value for key, value in sample.items()
                             if key not in {"seconds", "cpu_seconds", "peak_rss_bytes"}}
                    if reference is None:
                        reference = proof
                    elif proof != reference:
                        raise RuntimeError(f"WGD selected evidence differs: {profile} {name} trial {trial}")
                    if trial:
                        samples[name].append(sample)
                    print(profile, name, "warmup" if trial == 0 else trial, sample, flush=True)
            results["samples"][profile] = samples
    results["medians"] = {profile: {name: {metric: statistics.median(row[metric] for row in samples)
                                                   for metric in ("seconds", "cpu_seconds", "peak_rss_bytes")}
                                      for name, samples in groups.items()}
                          for profile, groups in results["samples"].items()}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(results["medians"], indent=2))


if __name__ == "__main__":
    main()
