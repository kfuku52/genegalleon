#!/usr/bin/env python3
"""Compare gene-count formatting with a saved baseline; exact TSV bytes must match.

Run through workflow/tests/run_in_runtime.sh with --baseline /path/to/old.py
and --output /path/to/results.json. Timings exclude fixture preparation. CLI
fixtures are prepared in the parent; CLI peak RSS covers worker imports, the
formatter run and bounded output hashing. Kernel profiles retain their in-worker
fixtures. Use --baseline-only before implementation and --cli-rows to scale only
the hog2og/sp2og CLI fixtures.
"""

import argparse
import contextlib
import hashlib
import importlib.util
import io
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

PROFILE_SHAPES = {
    "dense": (50000, 64),
    "sparse": (50000, 64),
    "wide": (10000, 256),
    "long": (10000, 32),
    "hog2og": (20000, 64),
    "sp2og": (20000, 64),
}


def load(path):
    spec = importlib.util.spec_from_file_location("formatter", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fixture(profile, rows=None):
    import numpy
    import pandas

    default_rows, species = PROFILE_SHAPES[profile]
    rows = default_rows if rows is None else rows
    values = ["gene001", "gene001, gene002", "gene001, gene002, gene003", ""]
    if profile == "sparse":
        values = ["", "", "", "gene001"]
    if profile == "long":
        values = [", ".join("gene%05d" % i for i in range(6000)), "gene1", "", "gene2, gene3"]
        # Sparse giant cells exercise the former fixed-width-string blow-up.
        offsets = (numpy.arange(rows) % 997).tolist()
        data = {
            "sp%d" % col: [values[0] if pos == col else values[1 + pos % 3] for pos in offsets]
            for col in range(species)
        }
    else:
        offsets = (numpy.arange(rows) % len(values)).tolist()
        data = {"sp%d" % col: [values[(pos + col) % len(values)] for pos in offsets] for col in range(species)}
    df = pandas.DataFrame({"Orthogroup": ["HOG%07d" % row for row in range(rows)], **data})
    return df, list(data)


def write_cli_fixture(profile, path, rows):
    df, _ = fixture(profile, rows=rows)
    if profile == "hog2og":
        df = df.rename(columns={"Orthogroup": "HOG"})
        df["OG"] = "OG0000001"
        df["Gene Tree Parent Clade"] = "N0"
    else:
        df = df.rename(columns={"Orthogroup": "group_id"})
        df["group_id"] = range(len(df))
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def worker(args):
    module = load(args.script)
    if args.profile in {"hog2og", "sp2og"}:
        with tempfile.TemporaryDirectory(prefix="gg-formatter-benchmark-") as temporary:
            root = Path(temporary)
            started = time.perf_counter()
            started_cpu = time.process_time()
            with contextlib.redirect_stdout(io.StringIO()):
                module.run(
                    SimpleNamespace(file_orthogroup_table=str(args.source), mode=args.profile, dir_out=str(root / "output"))
                )
            seconds = time.perf_counter() - started
            cpu_seconds = time.process_time() - started_cpu
            proof = hashlib.sha256()
            for name in ("Orthogroups.tsv", "Orthogroups.GeneCount.tsv", "README.txt"):
                with (root / "output" / name).open("rb") as handle:
                    for block in iter(lambda: handle.read(1024 * 1024), b""):
                        proof.update(block)
    else:
        df, species = fixture(args.profile)
        started = time.perf_counter()
        started_cpu = time.process_time()
        output = module.build_gene_count_table(df, species)
        seconds = time.perf_counter() - started
        cpu_seconds = time.process_time() - started_cpu
        proof = hashlib.sha256(output.to_csv(sep="\t", index=False).encode())
    print(
        json.dumps(
            {
                "seconds": seconds,
                "cpu_seconds": cpu_seconds,
                "sha256": proof.hexdigest(),
                "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                * (1 if sys.platform == "darwin" else 1024),
            }
        )
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--baseline-only", action="store_true")
    parser.add_argument(
        "--implementation",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "support" / "orthogroup_table_formatter.py",
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--profiles", nargs="+", choices=PROFILE_SHAPES, default=list(PROFILE_SHAPES))
    parser.add_argument("--cli-rows", type=int, help="Override row count for hog2og/sp2og CLI fixtures only.")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--script", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--profile", choices=PROFILE_SHAPES, help=argparse.SUPPRESS)
    parser.add_argument("--source", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        if args.profile in {"hog2og", "sp2og"} and args.source is None:
            parser.error("CLI workers require --source")
        worker(args)
        return
    if args.baseline is None or args.output is None:
        parser.error("--baseline and --output are required")
    if not 1 <= args.repeats <= 20:
        parser.error("--repeats must be 1..20")
    if args.cli_rows is not None and args.cli_rows < 1:
        parser.error("--cli-rows must be positive")
    import numpy
    import pandas

    script_labels = ("baseline",) if args.baseline_only else ("baseline", "implementation")
    results = {
        "python": sys.version,
        "platform": platform.platform(),
        "pandas": pandas.__version__,
        "numpy": numpy.__version__,
        "warmup": 1,
        "repeats": args.repeats,
        "cache": "warm OS caches; fresh process per sample",
        "peak_rss_scope": "CLI: worker imports, formatter and bounded hashing; kernels: includes fixture creation",
        "workloads": {
            profile: {
                "rows": args.cli_rows if profile in {"hog2og", "sp2og"} and args.cli_rows is not None
                else PROFILE_SHAPES[profile][0],
                "species": PROFILE_SHAPES[profile][1],
            }
            for profile in args.profiles
        },
        "scripts": {
            label: {
                "path": str(getattr(args, label)),
                "sha256": hashlib.sha256(getattr(args, label).read_bytes()).hexdigest(),
            }
            for label in script_labels
        },
        "samples": {},
        "inputs": {},
    }
    with tempfile.TemporaryDirectory(prefix="gg-formatter-fixtures-") as temporary:
        for profile in args.profiles:
            source = None
            if profile in {"hog2og", "sp2og"}:
                source = Path(temporary) / profile / "N0.tsv"
                write_cli_fixture(profile, source, args.cli_rows)
                with source.open("rb") as handle:
                    results["inputs"][profile] = {
                        "bytes": source.stat().st_size,
                        "sha256": hashlib.file_digest(handle, "sha256").hexdigest(),
                    }
            samples = {label: [] for label in script_labels}
            for trial in range(args.repeats + 1):
                labels = script_labels if trial % 2 == 0 else tuple(reversed(script_labels))
                for label in labels:
                    command = [
                        sys.executable,
                        str(Path(__file__).resolve()),
                        "--worker",
                        "--script",
                        str(getattr(args, label).resolve()),
                        "--profile",
                        profile,
                    ]
                    if source is not None:
                        command.extend(["--source", str(source)])
                    sample = json.loads(subprocess.check_output(command, text=True))
                    print(profile, label, "warmup" if trial == 0 else trial, sample, flush=True)
                    if trial:
                        samples[label].append(sample)
            assert len({sample["sha256"] for group in samples.values() for sample in group}) == 1, profile
            results["samples"][profile] = samples
    results["medians"] = {
        profile: {
            label: {
                metric: statistics.median(sample[metric] for sample in samples)
                for metric in ("seconds", "cpu_seconds", "peak_rss_bytes")
            }
            for label, samples in groups.items()
        }
        for profile, groups in results["samples"].items()
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, indent=2) + "\n")
    print(json.dumps(results["medians"], indent=2))


if __name__ == "__main__":
    main()
