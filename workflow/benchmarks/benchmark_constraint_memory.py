#!/usr/bin/env python3
"""Compare constraint validation RAM with fixtures created outside fresh workers.

Timing covers the real run(args), including its full input hashing. Linux VmHWM
is captured before benchmark result hashing; imports remain in whole-worker RSS.
No IQ-TREE search is timed. Default large input is 128 taxa by 1,000,000 columns.
"""

import argparse
import contextlib
import gzip
import hashlib
import importlib.util
import io
import json
import platform
import random
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from types import SimpleNamespace

SUPPORT = Path(__file__).resolve().parents[1] / "support"


def digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def load(path):
    spec = importlib.util.spec_from_file_location("measured_constraint", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def peak_bytes():
    for line in Path("/proc/self/status").read_text().splitlines():
        if line.startswith("VmHWM:"):
            return int(line.split()[1]) * 1024
    raise RuntimeError("Linux VmHWM is required")


def prepare_profile(root, taxa, columns, transports):
    root.mkdir(parents=True)
    labels = [f"taxon{i:04d}" for i in range(taxa)]
    alignment = root / "alignment.fa"
    rng = random.Random(20261025)
    lines = ["".join(rng.choices("ACDEFGHIKLMNPQRSTVWY-", k=80)) for _ in range(4096)]
    with alignment.open("w") as handle:
        for index, label in enumerate(labels):
            handle.write(f">{label} retained description\n")
            for offset in range(0, columns, 80):
                handle.write(lines[(index * 97 + offset // 80) % len(lines)][:min(80, columns - offset)] + "\n")
    paths = {"plain": alignment}
    if "gzip" in transports:
        compressed = alignment.with_suffix(".fa.gz")
        with alignment.open("rb") as source, compressed.open("wb") as raw:
            with gzip.GzipFile(fileobj=raw, filename="", mode="wb", mtime=0) as target:
                while block := source.read(1024 * 1024):
                    target.write(block)
        paths["gzip"] = compressed
    selection = root / "taxa.txt"
    selected = labels[::2]
    selection.write_text("\n".join(reversed(selected)) + "\n")
    constraints = {}
    for mode, names in (("all", labels), ("subset", selected)):
        constraint = root / f"{mode}.nwk"
        constraint.write_text("((" + ",".join(names[:2]) + ")," + ",".join(names[2:]) + ");\n")
        constraints[mode] = constraint
    inputs = {path.name: {"sha256": digest(path), "bytes": path.stat().st_size}
              for path in (*paths.values(), selection, *constraints.values())}
    return paths, selection, constraints, inputs


def write_case(path, alignment, constraint, taxa_list, output, inputs, **extra):
    value = {"args": {"alignment": str(alignment), "constraint": str(constraint),
                      "taxa_list": str(taxa_list) if taxa_list else "", "scope": "complete",
                      "tree": str(constraint), "output": output},
             "input_sha256": inputs, **extra}
    path.write_text(json.dumps(value, indent=2) + "\n")
    return path


def prepare_probes(root):
    root.mkdir()
    good = b"".join(b">" + label + b"\nACD\n" for label in (b"A", b"B", b"C", b"D"))
    cases = [("empty", b"", "ValueError", "unique, nonempty"),
             ("late-duplicate", good + b">A\nACD\n", "ValueError", "unique, nonempty"),
             ("late-length", good + b">E\nACDE\n", "ValueError", "same positive length"),
             ("gzip-tail", gzip.compress(good + b">A\nACD\n>late\n" + b"ACD" * 4000 + b"\n", mtime=0)[:-8],
              "EOFError", "end-of-stream")]
    probes = {}
    for name, payload, error, message in cases:
        source = root / (name + (".fa.gz" if name == "gzip-tail" else ".fa"))
        source.write_bytes(payload)
        # Missing later inputs ensure alignment parse/axis errors take priority.
        for output in (False, True):
            case_name = name + ("-output" if output else "-validate")
            probes[case_name] = write_case(root / f"{case_name}.json", source, root / "missing.nwk",
                                           root / "missing.txt", output,
                                           {source.name: {"sha256": digest(source), "bytes": len(payload)}},
                                           error=error, message=message)
    return probes


def worker(args):
    module = load(args.implementation)
    case = json.loads(args.worker.read_text())
    call_args = SimpleNamespace(**case["args"])
    call_args.output = str(args.worker_output) if call_args.output else None
    if "error" in case:
        args.worker_output.write_text("previous output\n")
        try:
            module.run(call_args)
        except Exception as error:
            if type(error).__name__ != case["error"] or case["message"] not in str(error):
                raise
            if args.worker_output.read_text() != "previous output\n":
                raise RuntimeError("Invalid input replaced the existing output") from error
            return {"exception_type": type(error).__name__, "message": str(error), "output_preserved": True,
                    "input_sha256": case["input_sha256"]}
        raise RuntimeError("Invalid benchmark input was accepted")
    stdout = io.StringIO()
    started_wall, started_cpu = time.perf_counter(), time.process_time()
    with contextlib.redirect_stdout(stdout):
        module.run(call_args)
    wall, cpu = time.perf_counter() - started_wall, time.process_time() - started_cpu
    peak = peak_bytes()
    report = json.loads(stdout.getvalue())
    for role in ("alignment", "constraint", "taxa_list"):
        source = getattr(call_args, role)
        if source and report[role + "_sha256"] != case["input_sha256"][Path(source).name]["sha256"]:
            raise RuntimeError(f"Input identity mismatch: {role}")
    return {"seconds": wall, "cpu_seconds": cpu, "peak_rss_bytes": peak,
            "report_sha256": hashlib.sha256(json.dumps(report, sort_keys=True).encode()).hexdigest(),
            "output_sha256": digest(args.worker_output) if call_args.output else None,
            "input_sha256": case["input_sha256"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--baseline-only", action="store_true")
    parser.add_argument("--implementation", type=Path, default=SUPPORT / "species_tree_constraint.py")
    parser.add_argument("--profiles", nargs="+", choices=("small", "large"), default=["small", "large"])
    parser.add_argument("--transports", nargs="+", choices=("plain", "gzip"), default=["plain", "gzip"])
    parser.add_argument("--selections", nargs="+", choices=("all", "subset"), default=["all", "subset"])
    parser.add_argument("--scopes", nargs="+", choices=("validate", "output"), default=["validate"])
    parser.add_argument("--taxa", type=int, default=128)
    parser.add_argument("--columns", type=int, default=1000000)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--worker", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--worker-output", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    args.implementation = args.implementation.resolve()
    if args.worker:
        print(json.dumps(worker(args)))
        return
    if args.taxa < 8 or args.columns < 1 or not 1 <= args.repeats <= 20:
        parser.error("--taxa must be >=8, --columns positive and --repeats 1..20")
    if not args.baseline or args.output is None or args.output.exists():
        parser.error("--baseline and a new --output are required")
    variants = {"baseline": args.baseline.resolve()}
    if not args.baseline_only:
        variants["implementation"] = args.implementation
    payload = {"python": sys.version, "platform": platform.platform(), "warmups": 1, "repeats": args.repeats,
               "scripts": {name: {"path": str(path), "sha256": digest(path)} for name, path in variants.items()},
               "scope": "Synthetic concatenated alignment validation; no IQ-TREE search; real input hashing included; "
                        "fixture creation/imports/result hashing excluded from time; VmHWM before result hashing; "
                        "fresh process per sample, warm OS cache.",
               "fixtures": {}, "samples": {}, "warmup_samples": {}, "medians": {}, "contract_probes": {}}

    def invoke(source, metadata, output):
        command = [sys.executable, str(Path(__file__).resolve()), "--implementation", str(source),
                   "--worker", str(metadata), "--worker-output", str(output)]
        return json.loads(subprocess.check_output(command, text=True))

    with tempfile.TemporaryDirectory(prefix="gg-constraint-memory-") as temporary:
        root = Path(temporary)
        for profile in dict.fromkeys(args.profiles):
            taxa, columns = (8, 1000) if profile == "small" else (args.taxa, args.columns)
            directory = root / profile
            paths, selection, constraints, inputs = prepare_profile(directory, taxa, columns, args.transports)
            for transport in dict.fromkeys(args.transports):
                for mode in dict.fromkeys(args.selections):
                    for scope in dict.fromkeys(args.scopes):
                        name = f"{profile}-{transport}-{mode}-{scope}"
                        metadata = write_case(directory / f"{name}.json", paths[transport], constraints[mode],
                                              selection if mode == "subset" else None, scope == "output", inputs)
                        payload["fixtures"][name] = {"taxa": taxa, "columns": columns, "input_sha256": inputs}
                        samples, warmups, expected = {label: [] for label in variants}, {}, None
                        for trial in range(args.repeats + 1):
                            labels = list(variants)
                            if trial % 2:
                                labels.reverse()
                            for label in labels:
                                output = directory / f"{name}-{label}-{trial}.fa"
                                sample = invoke(variants[label], metadata, output)
                                fingerprint = (sample["report_sha256"], sample["output_sha256"])
                                if expected is not None and fingerprint != expected:
                                    raise RuntimeError(f"Constraint reports/FASTA differ: {name}")
                                expected = fingerprint
                                if trial:
                                    samples[label].append(sample)
                                else:
                                    warmups[label] = sample
                                output.unlink(missing_ok=True)
                                print(name, label, "warmup" if not trial else trial,
                                      sample["seconds"], sample["peak_rss_bytes"], flush=True)
                        payload["samples"][name], payload["warmup_samples"][name] = samples, warmups
                        payload["medians"][name] = {
                            label: {metric: statistics.median(row[metric] for row in rows)
                                    for metric in ("seconds", "cpu_seconds", "peak_rss_bytes")}
                            for label, rows in samples.items()}
        for name, metadata in prepare_probes(root / "probes").items():
            results = {label: invoke(source, metadata, root / "probes" / f"{name}-{label}.fa")
                       for label, source in variants.items()}
            if len({json.dumps(result, sort_keys=True) for result in results.values()}) != 1:
                raise RuntimeError(f"Exception/output-protection contract differs: {name}")
            payload["contract_probes"][name] = results
    payload["comparison_performed"] = len(variants) > 1
    payload["equivalent"] = True if len(variants) > 1 else None
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload["medians"], indent=2))


if __name__ == "__main__":
    main()
