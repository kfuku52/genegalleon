#!/usr/bin/env python3
"""Compare complete HGT taxonomy migration CLIs with exact output bytes.

Run through workflow/tests/run_in_runtime.sh. Each sample uses a fresh worker,
imports the requested script and runs its real argument parser/main function.
Fixture creation and copies are outside timing; imports, all table reads,
enrichment and atomic replacements are included. Peak RSS covers the worker,
including the harness imports. The taxonomy DB is disabled: fixtures include
existing taxonomy values, as in migrations of previously enriched tables.
"""

import argparse
import csv
import hashlib
import json
import platform
import resource
import runpy
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

PROFILES = {"small": (32, 10), "screen": (2000, 10), "large": (4000, 10)}
TABLES = ("gene.tsv", "branch.tsv", "orthogroup.tsv")
SUPPORT = Path(__file__).resolve().parents[1] / "support"


def write_table(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fixture(root, profile):
    groups, per_group = PROFILES[profile]
    fields = ["orthogroup", "gene_id", "gene_taxon", "besthit_organism", "besthit_taxid",
              "recipient_phylum", "donor_phylum", "recipient_taxonomy", "donor_taxonomy", "old_score"]

    def genes():
        # Interleave groups: grouping must retain first-seen gene order.
        for slot in range(per_group):
            for group in range(groups):
                yield dict(orthogroup=f"OG{group:07d}", gene_id=f"gene{group}_{slot}",
                           gene_taxon="Example recipient", besthit_organism=f"Example donor {slot % 3}",
                           besthit_taxid="", recipient_phylum="Recipient phylum",
                           donor_phylum=f"Donor phylum {slot % 3}",
                           recipient_taxonomy="phylum:Recipient phylum",
                           donor_taxonomy=f"phylum:Donor phylum {slot % 3}", old_score=str(slot))

    write_table(root / TABLES[0], fields, genes())
    write_table(root / TABLES[1], ["orthogroup", "branch_id", "candidate_genes", "old_value"],
                (dict(orthogroup=f"OG{group:07d}", branch_id=str(group),
                      candidate_genes="; ".join(f"gene{group}_{slot}" for slot in reversed(range(per_group))),
                      old_value="preserve") for group in reversed(range(groups))))
    write_table(root / TABLES[2], ["orthogroup", "old_value"],
                (dict(orthogroup=f"OG{group:07d}", old_value="preserve") for group in reversed(range(groups))))


def hashes(root):
    return {name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in TABLES}


def worker(args):
    sys.path.insert(0, str(SUPPORT))
    sys.argv = [str(args.script), "--branch_tsv", str(args.fixture_root / "branch.tsv"),
                "--gene_tsv", str(args.fixture_root / "gene.tsv"),
                "--orthogroup_tsv", str(args.fixture_root / "orthogroup.tsv"), "--taxonomy_dbfile", ""]
    started = time.perf_counter()
    started_cpu = time.process_time()
    runpy.run_path(str(args.script), run_name="__main__")
    wall_seconds = time.perf_counter() - started
    cpu_seconds = time.process_time() - started_cpu
    # Linux VmHWM is the kernel-maintained peak, not a sampled current RSS.
    status = Path("/proc/self/status")
    if status.is_file():
        peak_rss_bytes = int(next(line.split()[1] for line in status.read_text().splitlines()
                                  if line.startswith("VmHWM:"))) * 1024
    else:
        peak_rss_bytes = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (
            1 if sys.platform == "darwin" else 1024)
    result = dict(wall_seconds=wall_seconds, cpu_seconds=cpu_seconds, peak_rss_bytes=peak_rss_bytes,
                  output_sha256=hashes(args.fixture_root))
    args.worker_output.write_text(json.dumps(result) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--implementation", type=Path, default=SUPPORT / "add_hgt_taxonomy_columns.py")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--profiles", nargs="+", choices=PROFILES, default=list(PROFILES))
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--script", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--fixture-root", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--worker-output", type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    if args.baseline is None or args.output is None:
        parser.error("--baseline and --output are required")
    if not 1 <= args.repeats <= 20:
        parser.error("--repeats must be 1..20")
    implementations = {"baseline": args.baseline.resolve(), "implementation": args.implementation.resolve()}
    report = dict(python=sys.version, platform=platform.platform(), command=sys.argv,
                  warmups=1, repeats=args.repeats, taxonomy_db="disabled; existing rank/lineage fields populated",
                  scripts={label: dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())
                           for label, path in implementations.items()},
                  shared_dependencies={name: hashlib.sha256((SUPPORT / name).read_bytes()).hexdigest()
                                       for name in ("score_hgt_candidates.py", "scaffold_taxonomy.py")}, profiles={})
    with tempfile.TemporaryDirectory(prefix="gg-taxonomy-migration-benchmark-") as temporary:
        temporary = Path(temporary)
        for profile in args.profiles:
            source = temporary / profile
            source.mkdir()
            fixture(source, profile)
            reference = None
            samples = {label: [] for label in implementations}
            for trial in range(args.repeats + 1):
                order = list(implementations) if trial % 2 == 0 else list(reversed(implementations))
                for label in order:
                    work = temporary / f"{profile}-{label}-{trial}"
                    shutil.copytree(source, work)
                    metrics = work / "metrics.json"
                    completed = subprocess.run(
                        [sys.executable, str(Path(__file__).resolve()), "--worker", "--script",
                         str(implementations[label]), "--fixture-root", str(work),
                         "--worker-output", str(metrics)], check=True, capture_output=True)
                    row = json.loads(metrics.read_text())
                    payload = {name: (work / name).read_bytes() for name in TABLES}
                    payload.update(stdout=completed.stdout, stderr=completed.stderr)
                    row["stdout_sha256"] = hashlib.sha256(completed.stdout).hexdigest()
                    row["stderr_sha256"] = hashlib.sha256(completed.stderr).hexdigest()
                    if reference is None:
                        reference = payload
                    else:
                        for name, expected in reference.items():
                            assert payload[name] == expected, (profile, label, trial, name)
                    print(profile, label, "warmup" if trial == 0 else trial,
                          row["wall_seconds"], row["cpu_seconds"], row["peak_rss_bytes"], flush=True)
                    if trial:
                        samples[label].append(row)
                    shutil.rmtree(work)
            report["profiles"][profile] = dict(orthogroups=PROFILES[profile][0],
                                                genes=PROFILES[profile][0] * PROFILES[profile][1],
                                                input_sha256=hashes(source), samples=samples,
                                                medians={label: {metric: statistics.median(row[metric] for row in rows)
                                                                 for metric in ("wall_seconds", "cpu_seconds", "peak_rss_bytes")}
                                                         for label, rows in samples.items()})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({profile: data["medians"] for profile, data in report["profiles"].items()}, indent=2))


if __name__ == "__main__":
    main()
