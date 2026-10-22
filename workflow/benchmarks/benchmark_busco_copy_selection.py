#!/usr/bin/env python3
"""Measure BUSCO extraction time/RSS and compare every output byte with a saved baseline."""

import argparse
import contextlib
import cProfile
import gzip
import hashlib
import importlib.metadata
import io
import json
import platform
import pstats
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CASES = {
    "protein": dict(species=192, markers=512, copies=4, mode="protein", duplicate_every=4),
    "cds": dict(species=64, markers=512, copies=3, mode="cds", duplicate_every=4),
    "many-copies": dict(species=8, markers=128, copies=64, mode="protein", duplicate_every=1),
}


def fixture(root, spec):
    full = root / "full"
    full.mkdir(parents=True)
    species = [f"Species_{i:04d}" for i in range(spec["species"])]
    cells = {marker: [] for marker in range(spec["markers"])}
    paths = []
    records = 0
    for label in species:
        path = root / f"{label}.fa.gz"
        paths.append(path)
        with gzip.GzipFile(filename=path, mode="wb", mtime=0) as fasta, (full / f"{label}.busco.full.tsv").open("w") as table:
            for marker in cells:
                duplicate = marker % spec["duplicate_every"] == 0
                count = spec["copies"] if duplicate else 1
                genes = []
                for copy in range(count):
                    gene = f"{label}_g{marker:05d}_{copy:03d}"
                    genes.append(gene)
                    # Different complete sequences and scores; no ambiguous/invalid cases.
                    peptide = "M" + "ACDEFGHIKLMNPQRSTVWY"[copy % 19] + "A" * 298
                    sequence = "ATG" + ("GCT" if copy % 2 else "GCC") * 299 if spec["mode"] == "cds" else peptide
                    fasta.write(f">{gene} description\n{sequence}\n".encode())
                    table.write(f"b{marker:05d}\t{'Duplicated' if duplicate else 'Complete'}\t{gene}\t{100 + copy}\t300\n")
                    records += 1
                cells[marker].append(",".join(genes))
    (root / "summary.tsv").write_text("busco_id\torthodb_url\tdescription\t" + "\t".join(species) + "\n" +
        "".join(f"b{marker:05d}\t-\t-\t" + "\t".join(values) + "\n" for marker, values in cells.items()))
    (root / "fastas.txt").write_text("".join(f"{path}\n" for path in paths))
    (root / "species.tsv").write_text("species\tfasta\n" + "".join(
        f"{label}\t{path}\n" for label, path in zip(species, paths, strict=True)))
    return records


def fingerprint(root):
    proof = hashlib.sha256()
    for path in sorted((root / "output").iterdir()) + [root / "counts.tsv", root / "selection.tsv"]:
        proof.update(path.name.encode() + b"\0")
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                proof.update(chunk)
    return proof.hexdigest()


def worker(args):
    sys.path.insert(0, str(args.support))
    import batch_extract_busco_fasta as batch
    import busco_copy_selection as selection

    spec = CASES[args.worker]
    root = args.fixture
    options = batch.parser().parse_args([
        "--summary", str(root / "summary.tsv"), "--fasta-list", str(root / "fastas.txt"),
        "--output-dir", str(root / "output"), "--suffix", ".fa", "--report", str(root / "counts.tsv"),
        "--mode", "single-copy", "--duplicate-policy", "best_score", "--require-complete",
        "--busco-full-dir", str(root / "full"), "--species-fasta-map", str(root / "species.tsv"),
        "--selection-report", str(root / "selection.tsv"), "--sequence-mode", spec["mode"],
    ])
    phases = {}
    for name in ("__init__", "extract"):
        original = getattr(selection.BestScoreSelection, name)

        def timed(self, *a, _original=original, _name=name, **kw):
            started = time.perf_counter()
            result = _original(self, *a, **kw)
            phases[_name] = time.perf_counter() - started
            return result

        setattr(selection.BestScoreSelection, name, timed)
    profiler = cProfile.Profile() if args.profile else None
    started = time.perf_counter()
    if profiler:
        profiler.enable()
    with contextlib.redirect_stdout(io.StringIO()):
        status = batch.run(options)
    seconds = time.perf_counter() - started
    if profiler:
        profiler.disable()
    if status:
        raise RuntimeError(f"Extraction returned {status}")
    result = dict(seconds=seconds, phases=phases,
                  peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss *
                  (1 if sys.platform == "darwin" else 1024), sha256=fingerprint(root))
    if profiler:
        stream = io.StringIO()
        pstats.Stats(profiler, stream=stream).sort_stats("cumulative").print_stats(25)
        result["profile"] = stream.getvalue()
    print(json.dumps(result))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-support", type=Path, required=True)
    parser.add_argument("--current-support", type=Path, default=ROOT / "workflow/support")
    parser.add_argument("--baseline-only", action="store_true")
    parser.add_argument("--cases", choices=CASES, nargs="+", default=list(CASES))
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--profile", action="store_true", help="Profile a separate untimed trial for each variant")
    parser.add_argument("--worker", choices=CASES, help=argparse.SUPPRESS)
    parser.add_argument("--support", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--fixture", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    if args.repeats < 1:
        parser.error("repeats must be positive")
    variants = {"baseline": args.baseline_support.resolve()}
    if not args.baseline_only:
        variants["current"] = args.current_support.resolve()
    result = dict(python=sys.version, platform=platform.platform(), command=sys.argv,
                  packages={name: importlib.metadata.version(name) for name in ("biopython", "pandas")},
                  warmups=1, repeats=args.repeats,
                  sources={name: {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                                  for path in [support / filename for filename in (
                                      "batch_extract_busco_fasta.py", "busco_copy_selection.py", "busco_table_utils.py",
                                      "collect_common_BUSCO_genes.py")] if path.exists()}
                           for name, support in variants.items()}, cases={})
    with tempfile.TemporaryDirectory(prefix="gg-busco-benchmark-") as temporary:
        for case in args.cases:
            root = Path(temporary) / case
            records = fixture(root, CASES[case])
            trials = {name: [] for name in variants}

            def trial(name, profile=False, *, case=case, root=root):
                command = [sys.executable, str(Path(__file__).resolve()), "--worker", case, "--support", str(variants[name]),
                           "--fixture", str(root), "--baseline-support", str(args.baseline_support), "--output", str(args.output)]
                if profile:
                    command.append("--profile")
                completed = subprocess.run(command, capture_output=True, text=True)
                if completed.returncode:
                    raise RuntimeError(f"Benchmark worker failed: {completed.stdout}\n{completed.stderr}")
                return json.loads(completed.stdout)

            reference = None
            for iteration in range(args.repeats + 1):
                order = list(variants) if iteration % 2 == 0 else list(reversed(variants))
                for name in order:
                    value = trial(name)
                    if reference is None:
                        reference = value["sha256"]
                    if value["sha256"] != reference:
                        raise RuntimeError(f"Output mismatch: {case}: {name}")
                    if iteration:
                        trials[name].append(value)
            result["cases"][case] = dict(workload=CASES[case], records=records, variants={})
            for name, values in trials.items():
                result["cases"][case]["variants"][name] = dict(
                    median_seconds=statistics.median(value["seconds"] for value in values),
                    median_peak_rss_bytes=statistics.median(value["peak_rss_bytes"] for value in values), trials=values)
                if args.profile:
                    profiled = trial(name, profile=True)
                    if profiled["sha256"] != reference:
                        raise RuntimeError(f"Profile output mismatch: {case}: {name}")
                    result["cases"][case]["variants"][name]["profile"] = profiled["profile"]
            print(case, {name: round(value["median_seconds"], 3)
                         for name, value in result["cases"][case]["variants"].items()}, flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
