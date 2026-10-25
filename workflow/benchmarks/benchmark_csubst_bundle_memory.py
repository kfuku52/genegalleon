#!/usr/bin/env python3
"""Compare structural-tip export RAM using complete synthetic CSUBST bundles.

Run through workflow/tests/run_in_runtime.sh. The controller creates fixtures
before starting fresh measured workers; real finalize()/bundle verification is
used without mocks. Fit files are synthetic payloads, not predictor results.
Timing covers write_structural_tip_alignment(), including bundle verification.
Peak RSS covers the whole worker, including imports and result/input hashing.
"""

import argparse
import gzip
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
PROFILES = ("small", "tips", "ancestors")


def digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def load(path):
    spec = importlib.util.spec_from_file_location("measured_csubst_bundle", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_record(handle, identifier, sequence):
    handle.write(f">{identifier} synthetic description\n")
    for offset in range(0, len(sequence), 80):
        handle.write(sequence[offset:offset + 80] + "\n")


def prepare_fixture(root, profile, args):
    import numpy

    module = load(args.implementation)
    tips = 16 if profile == "small" else args.tips
    codons = 100 if profile == "small" else args.codons
    ancestors = (tips - 1 if args.ancestors is None else args.ancestors) if profile == "ancestors" else 0
    directory = root / module.STRUCTURAL_DIR
    (directory / "inspect").mkdir(parents=True)
    dna = "ATGGCTCAATGG---NNNACCGGT" * ((codons + 7) // 8)
    states = "ACDEFGHIKLMNPQRSTVWY-X" * ((codons + 20) // 21)
    with (directory / "csubst.fasta").open("w") as handle:
        for index in range(tips):
            shift = index % 8 * 3
            sequence = (dna[shift:] + dna[:shift])[:codons * 3]
            write_record(handle, f"tip{index:07d}", sequence)
    structural = directory / "inspect/csubst_alignment_3di.fa"
    with structural.open("w") as handle:
        for index in range(ancestors):
            write_record(handle, f"Node{index:07d}", states[:codons])
        # Structural order differs from CDS order, exercising tip selection.
        for index in reversed(range(tips)):
            shift = index % 21
            sequence = (states[shift:] + states[:shift])[:codons]
            write_record(handle, f"tip{index:07d}", sequence)
    tree = "(" + ",".join(f"tip{i:07d}:1" for i in range(tips)) + ");\n"
    for suffix in module.FIT_SUFFIXES:
        path = directory / f"csubst.{suffix}"
        if suffix == "fasta":
            continue
        if suffix == "ckp.gz":
            path.write_bytes(gzip.compress(b"synthetic checkpoint\n", mtime=0))
        else:
            path.write_text(tree if suffix in {"nwk", "treefile"} else f"synthetic {suffix}\n")
    numpy.savez(directory / "csubst_3di_state_cache.npz", synthetic=numpy.arange(20, dtype="uint8"))
    module.finalize(root, 1)
    alignment = directory / "csubst.fasta"
    if args.gzip_cds:
        alignment = directory / "csubst.fasta.gz"
        with (directory / "csubst.fasta").open("rb") as source, alignment.open("wb") as raw:
            with gzip.GzipFile(fileobj=raw, filename="", mode="wb", mtime=0) as target:
                while block := source.read(1024 * 1024):
                    target.write(block)
    files = {path.relative_to(root).as_posix(): {"sha256": digest(path), "bytes": path.stat().st_size}
             for path in sorted(root.rglob("*")) if path.is_file()}
    fixture = {"profile": profile, "tips": tips, "codons": codons, "ancestors": ancestors,
               "cds_sequence_bytes": tips * codons * 3,
               "alignment": alignment.relative_to(root).as_posix(), "files": files}
    metadata = root / "fixture.json"
    metadata.write_text(json.dumps(fixture, indent=2) + "\n")
    return metadata, fixture


def worker(args):
    module = load(args.implementation)
    fixture = json.loads(args.worker.read_text())
    root = args.worker.parent
    alignment = root / fixture["alignment"]
    started_wall, started_cpu = time.perf_counter(), time.process_time()
    module.write_structural_tip_alignment(alignment, args.worker_output)
    wall, cpu = time.perf_counter() - started_wall, time.process_time() - started_cpu
    # Hash outside the timer, with bounded buffers. The production call above
    # still performs its complete normal manifest/member verification.
    inputs = {name: digest(root / name) for name in fixture["files"]}
    if inputs != {name: item["sha256"] for name, item in fixture["files"].items()}:
        raise RuntimeError("Benchmark inputs changed")
    output_hash = digest(args.worker_output)
    peak_bytes = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == "darwin" else 1024)
    return {"seconds": wall, "cpu_seconds": cpu, "peak_rss_bytes": peak_bytes,
            "output_sha256": output_hash, "input_sha256": inputs}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--baseline-only", action="store_true")
    parser.add_argument("--implementation", type=Path, default=SUPPORT / "csubst_input_bundle.py")
    parser.add_argument("--profiles", choices=PROFILES, nargs="+", default=list(PROFILES))
    parser.add_argument("--tips", type=int, default=4096)
    parser.add_argument("--codons", type=int, default=5120)
    parser.add_argument("--ancestors", type=int, help="Defaults to tips minus one in the ancestors profile")
    parser.add_argument("--gzip-cds", action="store_true")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--worker", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--worker-output", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    args.implementation = args.implementation.resolve()
    if args.worker:
        print(json.dumps(worker(args)))
        return
    if min(args.tips, args.codons) < 1 or not 1 <= args.repeats <= 20:
        parser.error("--tips/--codons must be positive and --repeats must be 1..20")
    if args.ancestors is not None and args.ancestors < 0:
        parser.error("--ancestors must be non-negative")
    if args.output is None or args.output.exists():
        parser.error("--output must name a new file")
    if not args.baseline:
        parser.error("--baseline is required, including for --baseline-only")
    variants = {"baseline": args.baseline.resolve()}
    if not args.baseline_only:
        variants["implementation"] = args.implementation
    else:
        args.implementation = variants["baseline"]
    payload = {"python": sys.version, "platform": platform.platform(),
               "warmups": 1, "repeats": args.repeats, "gzip_cds": args.gzip_cds,
               "scripts": {name: {"path": str(path), "sha256": digest(path)}
                           for name, path in variants.items()},
               "scope": "Synthetic alignment export; complete production bundle verification included; "
                        "fixture creation/imports/result hashing excluded from time; whole-worker peak RSS; "
                        "warm OS cache, fresh process per sample; no predictor or whole-workflow speed claim.",
               "fixtures": {}, "warmup_samples": {}, "samples": {}, "medians": {}}
    with tempfile.TemporaryDirectory(prefix="gg-csubst-bundle-memory-") as temporary:
        for profile in dict.fromkeys(args.profiles):
            root = Path(temporary) / profile
            fixture_path, fixture = prepare_fixture(root, profile, args)
            payload["fixtures"][profile] = fixture
            samples = {name: [] for name in variants}
            warmups = {}
            expected_hash = None
            for trial in range(args.repeats + 1):
                names = list(variants)
                if trial % 2:
                    names.reverse()
                for name in names:
                    output = root / f"{name}-{trial}.fasta"
                    command = [sys.executable, str(Path(__file__).resolve()), "--implementation", str(variants[name]),
                               "--worker", str(fixture_path), "--worker-output", str(output)]
                    sample = json.loads(subprocess.check_output(command, text=True))
                    if expected_hash is None:
                        expected_hash = sample["output_sha256"]
                    if sample["output_sha256"] != expected_hash:
                        raise RuntimeError(f"Structural-tip output differs: {profile}, {name}, trial {trial}")
                    if trial:
                        samples[name].append(sample)
                    else:
                        warmups[name] = sample
                    output.unlink()
                    print(profile, name, "warmup" if not trial else trial,
                          sample["seconds"], sample["peak_rss_bytes"], flush=True)
            payload["samples"][profile] = samples
            payload["warmup_samples"][profile] = warmups
            payload["medians"][profile] = {
                name: {metric: statistics.median(sample[metric] for sample in group)
                       for metric in ("seconds", "cpu_seconds", "peak_rss_bytes")}
                for name, group in samples.items()}
    payload["equivalent"] = True
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload["medians"], indent=2))


if __name__ == "__main__":
    main()
