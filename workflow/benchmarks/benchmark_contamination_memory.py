#!/usr/bin/env python3
"""Compare contamination FASTA retention through the real main without network."""

import argparse
import contextlib
import gzip
import hashlib
import importlib.util
import io
import json
import os
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

SUPPORT = Path(__file__).resolve().parents[1] / "support"
OUTPUTS = ("clean_sequences.fa", "lineage_compatibility.tsv")
IO_COUNTERS = {"rchar": "io_rchar_bytes", "wchar": "io_wchar_bytes", "read_bytes": "io_read_bytes",
               "write_bytes": "io_write_bytes", "cancelled_write_bytes": "io_cancelled_write_bytes"}


class FakeNcbi:
    """Fixed two-class taxonomy; no local taxonomy database or downloads."""

    ranks = {1: "no rank", 2: "domain", 10: "class", 11: "species", 20: "class", 21: "species"}
    names = {1: "root", 2: "Eukaryota", 10: "FocalClass", 11: "Genus species",
             20: "OtherClass", 21: "Other species"}

    def get_rank(self, taxids):
        return {taxid: self.ranks[taxid] for taxid in taxids}

    def get_name_translator(self, names):
        return {name: [11] for name in names if name == "Genus species"}

    def get_lineage(self, taxid):
        return [1, 2, 10, 11]

    def get_taxid_translator(self, taxids):
        return {taxid: self.names[taxid] for taxid in taxids}


def digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def process_io():
    """Read Linux process counters outside the measured main interval."""
    if not sys.platform.startswith("linux"):
        return None
    try:
        counters = dict(line.split(":", 1) for line in Path("/proc/self/io").read_text().splitlines())
        return {name: int(counters[name]) for name in IO_COUNTERS}
    except (OSError, KeyError, ValueError):
        return None


def create_fixture(root, records, length, scenario, codec):
    fasta = root / ("input.fa.gz" if codec == "gzip" else "input.fa")
    retained = retained_bases = total_bases = 0
    with contextlib.ExitStack() as stack:
        if codec == "gzip":
            raw = stack.enter_context(fasta.open("wb"))
            compressed = stack.enter_context(gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0))
            sequences = stack.enter_context(io.TextIOWrapper(compressed, encoding="ascii"))
        else:
            sequences = stack.enter_context(fasta.open("w"))
        taxonomy = stack.enter_context((root / "taxonomy.tsv").open("w"))
        fx = stack.enter_context((root / "fx.tsv").open("w"))
        fx.write("#id\tlength\tGC\n")
        for number in range(records):
            identifier = f"seq{number:07d}"
            size = max(1, length - number % 97)
            sequence = ("ACGTN" * (size // 5 + 1))[:size]
            sequences.write(f">{identifier} description {number}\n{sequence}\n")
            fx.write(f"{identifier}\t{size}\t40\n")
            compatible = number % 1000 == 0 if scenario == "reject-most" else number % 20 != 19
            unclassified = compatible and number % 20 == 18
            if compatible:
                retained += 1
                retained_bases += size
            total_bases += size
            taxid, name, lineage = ((0, "unclassified", "") if unclassified else
                (11, "Genus species", "1;2;10;11") if compatible else
                (21, "Other species", "1;2;20;21"))
            taxonomy.write(f"{identifier}\t{taxid}\t-\t{name}\t-\t-\t-\t-\t{lineage}\n")
    paths = [fasta, root / "taxonomy.tsv", root / "fx.tsv"]
    return {"records": records, "nominal_length": length, "sequence_bases": total_bases,
            "retained_records": retained, "retained_bases": retained_bases,
            "inputs": {path.name: {"bytes": path.stat().st_size, "sha256": digest(path)} for path in paths}}


def worker(args):
    sys.path.insert(0, str(SUPPORT))
    specification = importlib.util.spec_from_file_location("measured_contamination", args.source)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    module.ncbi_taxa = lambda dbfile=None: FakeNcbi()
    os.environ.pop("GG_TAXONOMY_DBFILE", None)
    root = args.worker_root
    fasta = root / ("input.fa.gz" if args.codec == "gzip" else "input.fa")
    sys.argv = ["remove_contaminated_sequences.py", "--fasta_file", str(fasta),
                "--mmseqs2taxonomy_tsv", str(root / "taxonomy.tsv"),
                "--fx2tab_tsv", str(root / "fx.tsv"), "--species_name", "Genus_species",
                "--rank", "class", "--rename_seq", "yes", "--verbose", "no"]
    with tempfile.TemporaryDirectory(prefix="gg-contamination-output-") as temporary:
        output = Path(temporary)
        os.chdir(output)
        with (output / "stdout.log").open("w") as log, contextlib.redirect_stdout(log):
            io_before = process_io()
            before = resource.getrusage(resource.RUSAGE_SELF)
            started = time.perf_counter()
            module.main()
            seconds = time.perf_counter() - started
            after = resource.getrusage(resource.RUSAGE_SELF)
            io_after = process_io()
        io_available = io_before is not None and io_after is not None
        result = {"seconds": seconds,
                  "user_seconds": after.ru_utime - before.ru_utime,
                  "system_seconds": after.ru_stime - before.ru_stime,
                  "cpu_seconds": after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime,
                  "peak_rss_mib": after.ru_maxrss / (1024 * 1024 if sys.platform == "darwin" else 1024),
                  "io_source": "/proc/self/io" if io_available else "unavailable",
                  **{metric: io_after[name] - io_before[name] if io_available else None
                     for name, metric in IO_COUNTERS.items()},
                  "outputs": {name: {"sha256": digest(output / name), "bytes": (output / name).stat().st_size}
                              for name in OUTPUTS}}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SUPPORT / "remove_contaminated_sequences.py")
    parser.add_argument("--baseline-source", type=Path)
    parser.add_argument("--baseline-only", action="store_true")
    parser.add_argument("--scenario", choices=("retain-most", "reject-most", "small"), default="retain-most")
    parser.add_argument("--codec", choices=("plain", "gzip"), default="plain")
    parser.add_argument("--records", type=int, default=10000)
    parser.add_argument("--length", type=int, default=40960)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--worker-root", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if min(args.records, args.length, args.repeats) < 1:
        parser.error("records, length and repeats must be positive")
    args.source = args.source.resolve()
    if args.worker_root:
        print(json.dumps(worker(args)))
        return
    if args.output is None or args.output.exists():
        parser.error("--output must name a new file")
    variants = {"current": args.source}
    if args.baseline_source:
        variants = {"baseline": args.baseline_source.resolve(), **variants}
    if args.baseline_only:
        variants = {"baseline": (args.baseline_source or args.source).resolve()}
    records = min(args.records, 32) if args.scenario == "small" else args.records
    length = min(args.length, 2048) if args.scenario == "small" else args.length
    trials = {name: [] for name in variants}
    warmups = {}
    with tempfile.TemporaryDirectory(prefix="gg-contamination-fixture-") as temporary:
        root = Path(temporary)
        fixture = create_fixture(root, records, length, args.scenario, args.codec)
        for trial in range(args.repeats + 1):
            names = list(variants)
            if trial % 2:
                names.reverse()
            for name in names:
                command = [sys.executable, str(Path(__file__).resolve()), "--source", str(variants[name]),
                           "--worker-root", str(root), "--codec", args.codec]
                result = json.loads(subprocess.check_output(command, text=True))
                if trial:
                    trials[name].append(result)
                else:
                    warmups[name] = result
                print(f"{'warmup' if trial == 0 else trial} {name}: "
                      f"{result['seconds']:.3f}s, {result['peak_rss_mib']:.2f} MiB", flush=True)
    fingerprints = {json.dumps(item["outputs"], sort_keys=True)
                    for samples in trials.values() for item in samples}
    fingerprints.update(json.dumps(item["outputs"], sort_keys=True) for item in warmups.values())
    if len(fingerprints) != 1:
        raise RuntimeError("Contamination FASTA/lineage outputs differ")
    metrics = ("seconds", "cpu_seconds", "user_seconds", "system_seconds", "peak_rss_mib",
               *IO_COUNTERS.values())
    medians = {}
    for name, samples in trials.items():
        medians[name] = {}
        for metric in metrics:
            measured = [item[metric] for item in samples if item[metric] is not None]
            medians[name][metric] = statistics.median(measured) if measured else None
    payload = {"scenario": args.scenario, "codec": args.codec, "fixture": fixture,
               "warmup_count": 1, "repeats": args.repeats, "warmups": warmups, "trials": trials,
               "equivalent_outputs": True, "python": sys.version, "platform": platform.platform(),
               "source_sha256": {name: digest(path) for name, path in variants.items()},
               "medians": medians,
               "scope": "Synthetic contamination inputs and fixed two-class taxonomy; real main runs in fresh "
                        "processes. Fixture generation, imports and output hashing excluded from wall/CPU time. "
                        "Peak RSS includes imports and main; no concurrent benchmarks/builds should run. "
                        "Linux /proc/self/io deltas include a bounded counter-read overhead; rchar/wchar "
                        "count logical bytes, read_bytes/write_bytes reflect storage accounting rather than "
                        "physical device traffic, and cancelled_write_bytes identifies cancelled dirty writes. "
                        "I/O counters are null with io_source=unavailable on unsupported systems. "
                        "Comparison includes FASTA parsing, temporary storage and both output writes; "
                        "does not establish whole-workflow or production filesystem performance."}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload["medians"], sort_keys=True))


if __name__ == "__main__":
    main()
