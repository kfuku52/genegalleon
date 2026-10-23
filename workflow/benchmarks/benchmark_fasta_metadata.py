#!/usr/bin/env python3
"""Compare real FASTA metadata consumers without fixture-generation RSS.

Fixtures are written by the controller before any fresh measured worker starts.
The chosen real annotation package is loaded and its source identity asserted
before API/main timing. Imports, fixture creation and result hashing are excluded
except that the unchanged full-record control hashes each yielded sequence while
consuming it. Plan-main runs the real local planning CLI's main, including fresh
input hashing and JSON publication, with imports outside the timer. No network,
curated data, persistent digest reuse or real analysis outputs are involved.
"""

import argparse
import contextlib
import gzip
import hashlib
import importlib
import importlib.util
import io
import json
import platform
import random
import statistics
import subprocess
import sys
import tarfile
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SUPPORT = ROOT / "workflow/support"
MIB = 1024 * 1024
PROFILES = {
    "small-raw": ("small", "raw"), "small-gz": ("small", "gz"),
    "chr16-raw": ("chr16", "raw"), "chr16-gz": ("chr16", "gz"),
    "chr64-raw": ("chr64", "raw"), "chr64-gz": ("chr64", "gz"),
    "multi-raw": ("multi", "raw"), "multi-gz": ("multi", "gz"),
    "cds20k-raw": ("cds20k", "raw"), "cds20k-gz": ("cds20k", "gz"),
    "chr16-tar": ("chr16", "tar"),
}


def sha256_path(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def peak_mib():
    for line in Path("/proc/self/status").read_text().splitlines():
        if line.startswith("VmHWM:"):
            return int(line.split()[1]) / 1024
    raise RuntimeError("Linux VmHWM is required")


def create_fasta(path, kind):
    if kind == "small":
        sizes, prefix = [300] * 32, "chr"
    elif kind == "chr16":
        sizes, prefix = [16 * MIB, 256], "chr"
    elif kind == "chr64":
        sizes, prefix = [64 * MIB, 256], "chr"
    elif kind == "multi":
        sizes, prefix = [16 * MIB] * 4, "chr"
    else:
        sizes, prefix = [1500] * 20000, "gene"
    rng = random.Random(20261023)
    # More than a deflate window of seeded DNA avoids a homopolymer compression
    # toy while allowing bounded-memory streamed input creation.
    lines = ["".join(rng.choices("ACGTN", weights=(24, 24, 24, 24, 4), k=80)).encode() + b"\n"
             for _ in range(2048)]
    with path.open("wb") as handle:
        for index, size in enumerate(sizes):
            handle.write(f">{prefix}{index} literal description\n".encode())
            left, line_index = size, index % len(lines)
            while left:
                line = lines[line_index % len(lines)]
                width = min(80, left)
                handle.write(line[:width] + b"\n")
                left -= width
                line_index += 1
    return len(sizes)


def prepare_fixture(root, profile):
    kind, transport = PROFILES[profile]
    species = root / profile / "Species_one"
    species.mkdir(parents=True)
    raw = species / "source.fa"
    records = create_fasta(raw, kind)
    if transport == "gz":
        final = raw.with_suffix(".fa.gz")
        with raw.open("rb") as source, final.open("wb") as compressed:
            with gzip.GzipFile(filename="", fileobj=compressed, mode="wb", compresslevel=6, mtime=0) as target:
                while block := source.read(MIB):
                    target.write(block)
        raw.unlink()
    elif transport == "tar":
        final = raw.with_suffix(".fa.tar.gz")
        with final.open("wb") as compressed:
            with gzip.GzipFile(filename="", fileobj=compressed, mode="wb", compresslevel=6, mtime=0) as target:
                with tarfile.open(fileobj=target, mode="w") as archive, raw.open("rb") as source:
                    member = tarfile.TarInfo("member.fa")
                    member.size = raw.stat().st_size
                    archive.addfile(member, source)
        raw.unlink()
    else:
        final = raw
    matching = "chr3" if kind == "multi" else "chr0"
    gff = species / "annotation.gff"
    gff.write_text(f"{matching}\tsrc\tgene\t1\t30\t.\t+\t.\tID=g1\n")
    return {"input_dir": str(species.parent), "fasta": str(final), "gff": str(gff),
            "records": records, "fasta_bytes": final.stat().st_size,
            "fasta_sha256": sha256_path(final), "gff_sha256": sha256_path(gff)}


def load_annotation_package(common_path):
    sys.path.insert(0, str(SUPPORT))
    package = common_path.parent
    spec = importlib.util.spec_from_file_location("format_species_annotation", package / "__init__.py",
                                                 submodule_search_locations=[str(package)])
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    common = importlib.import_module("format_species_annotation.common")
    actual = Path(common.count_fasta_records.__code__.co_filename).resolve()
    assert actual == common_path.resolve(), (actual, common_path)
    return common


def serializable(value):
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {key: serializable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [serializable(item) for item in value]
    return value


def worker(args):
    fixture = json.loads(args.fixture.read_text())
    common_path = args.implementation.resolve()
    common = load_annotation_package(common_path)
    filename = Path(fixture["fasta"])
    if args.scope == "discovery":
        tasks = importlib.import_module("format_species_annotation.tasks")
        assert Path(tasks.fasta_looks_like_gff_genome.__code__.co_filename).resolve() == common_path

        def call():
            return tasks.discover_generic_species_dir_tasks("direct", Path(fixture["input_dir"]))
    elif args.scope == "count":
        def call():
            return common.count_fasta_records(filename)
    elif args.scope == "records":
        def call():
            digest, count = hashlib.sha256(), 0
            for header, sequence in common.iter_fasta_records(filename):
                digest.update(header.encode())
                digest.update(b"\0")
                digest.update(sequence.encode())
                digest.update(b"\0")
                count += 1
            return count, digest.hexdigest()
    else:
        planner = importlib.import_module("plan_input_generation_tasks")
        stdout, stderr = io.StringIO(), io.StringIO()
        temporary = tempfile.TemporaryDirectory(prefix="gg-fasta-plan-")
        output = Path(temporary.name) / "plan.json"
        sys.argv = ["plan_input_generation_tasks.py", "--provider", "direct", "--input-dir", fixture["input_dir"],
                    "--outfile", str(output), "--strict"]

        def call():
            with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                return planner.main()
    start_wall, start_cpu = time.perf_counter(), time.process_time()
    result = call()
    wall, cpu = time.perf_counter() - start_wall, time.process_time() - start_cpu
    rss = peak_mib()
    if args.scope == "plan-main":
        result = {"exit": result, "plan": output.read_text(),
                  "stdout": stdout.getvalue().replace(str(output), "<plan-output>"), "stderr": stderr.getvalue()}
        temporary.cleanup()
    signature = hashlib.sha256(json.dumps(serializable(result), sort_keys=True, ensure_ascii=True).encode()).hexdigest()
    print(json.dumps({"wall_s": wall, "cpu_s": cpu, "peak_mib": rss, "output_sha256": signature,
                      "source": str(common_path), "source_sha256": sha256_path(common_path),
                      "profile": args.profile, "scope": args.scope, "input": fixture}))


def controller(args):
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fixtures_root = args.fixtures or args.output.parent / (args.output.stem + "-fixtures")
    if args.fixtures is None:
        fixtures_root.mkdir()
    fixtures = {}
    for profile in args.profiles.split(","):
        manifest = fixtures_root / f"{profile}.json"
        if args.fixtures is None:
            data = prepare_fixture(fixtures_root, profile)
            manifest.write_text(json.dumps(data))
        else:
            data = json.loads(manifest.read_text())
            assert sha256_path(Path(data["fasta"])) == data["fasta_sha256"]
            assert sha256_path(Path(data["gff"])) == data["gff_sha256"]
        fixtures[profile] = manifest
    records = []
    sources = {"before": args.baseline}
    if not args.baseline_only:
        sources["after"] = args.implementation
    for scope in args.scopes.split(","):
        for profile, fixture in fixtures.items():
            for pair in range(-1, args.pairs):
                order = list(sources) if pair % 2 == 0 else list(reversed(sources))
                current = {}
                for variant in order:
                    command = [sys.executable, str(Path(__file__).resolve()), "--worker", "--implementation",
                               str(sources[variant]), "--fixture", str(fixture), "--profile", profile, "--scope", scope]
                    completed = subprocess.run(command, text=True, capture_output=True)
                    if completed.returncode:
                        print(completed.stderr, end="", file=sys.stderr)
                        completed.check_returncode()
                    sample = json.loads(completed.stdout)
                    sample.update(pair=pair, variant=variant)
                    current[variant] = sample
                    records.append(sample)
                    print(json.dumps(sample), flush=True)
                if not args.baseline_only:
                    assert current["before"]["output_sha256"] == current["after"]["output_sha256"], (scope, profile, pair)
                args.output.write_text(json.dumps({"records": records}, indent=2) + "\n")
    summary = []
    for scope in args.scopes.split(","):
        for profile in fixtures:
            samples = [item for item in records if item["scope"] == scope and item["profile"] == profile and item["pair"] >= 0]
            row = {"scope": scope, "profile": profile}
            for variant in sources:
                row[variant] = {metric: statistics.median(item[metric] for item in samples if item["variant"] == variant)
                                for metric in ("wall_s", "cpu_s", "peak_mib")}
            if not args.baseline_only:
                ratios = []
                for pair in range(args.pairs):
                    selected = {item["variant"]: item for item in samples if item["pair"] == pair}
                    ratios.append(selected["before"]["wall_s"] / selected["after"]["wall_s"])
                row["paired_speedup_range"] = [min(ratios), max(ratios)]
            summary.append(row)
    report = {"environment": {"python": sys.version, "platform": platform.platform()},
              "configuration": {"pairs": args.pairs, "scopes": args.scopes, "profiles": args.profiles},
              "records": records, "summary": summary}
    for manifest in fixtures.values():
        data = json.loads(manifest.read_text())
        assert sha256_path(Path(data["fasta"])) == data["fasta_sha256"]
        assert sha256_path(Path(data["gff"])) == data["gff_sha256"]
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--baseline-only", action="store_true")
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--implementation", type=Path, default=SUPPORT / "format_species_annotation/common.py")
    parser.add_argument("--fixture", type=Path)
    parser.add_argument("--fixtures", type=Path, help="Reuse and freshly verify controller-created fixture manifests.")
    parser.add_argument("--profile", choices=PROFILES)
    parser.add_argument("--scope", choices=("discovery", "count", "records", "plan-main"))
    parser.add_argument("--profiles", default=",".join(PROFILES))
    parser.add_argument("--scopes", default="discovery,count")
    parser.add_argument("--pairs", type=int, default=5)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    if not arguments.worker and (arguments.baseline is None or arguments.output is None or arguments.pairs < 1):
        parser.error("--baseline, --output and positive --pairs are required")
    worker(arguments) if arguments.worker else controller(arguments)
