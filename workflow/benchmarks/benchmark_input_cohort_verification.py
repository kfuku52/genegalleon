"""Compare native view preparation/verification, including current hash reads."""

import argparse
import hashlib
import importlib.util
import json
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SUPPORT = REPO / "workflow/support"
CANDIDATE = SUPPORT / "input_generation_cohort.py"
sys.path.insert(0, str(SUPPORT))


def module(path):
    s = importlib.util.spec_from_file_location("bench_cohort_" + hashlib.md5(str(path).encode()).hexdigest(), path)
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


if len(sys.argv) > 1 and sys.argv[1] == "--worker":
    sys.argv.pop(1)
    import resource
    import time
    import tracemalloc

    import input_generation_array_state as state

    m = module(Path(sys.argv[1]))
    workspace = Path(sys.argv[2])
    plan = Path(sys.argv[3])
    phase = sys.argv[4]
    calls = []
    original = state.digest

    def counted(path):
        calls.append(Path(path).stat().st_size)
        return original(path)

    state.digest = counted
    m.digest = counted
    tracemalloc.start()
    started = time.perf_counter()
    if phase == 'prepare':
        m.prepare(workspace, plan)
        doc = json.loads((workspace / 'cohort_view.json').read_text())
    else:
        doc = m.verify_view(workspace, plan)
    seconds = time.perf_counter() - started
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    print(
        json.dumps(
            dict(
                seconds=seconds,
                peak_tracemalloc_bytes=peak,
                peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                hash_reads=len(calls),
                hash_bytes=sum(calls),
                output_sha256=hashlib.sha256(json.dumps(doc, sort_keys=True).encode()).hexdigest(),
            )
        )
    )
    raise SystemExit(0)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--baseline", type=Path, required=True, help="Saved baseline input_generation_cohort.py")
parser.add_argument('--phase', choices=['verify', 'prepare'], default='verify')
args = parser.parse_args()
BASE = args.baseline.resolve(strict=True)
m = module(BASE)
state = __import__("input_generation_array_state")
atomic_json, digest = state.atomic_json, state.digest

with tempfile.TemporaryDirectory(prefix="cohort-benchmark-") as tmp:
    tmp = Path(tmp)
    rows = []
    for i in range(2):
        sp = "Bench_species" + str(i)
        workspace = tmp / ("donor" + str(i))
        root = workspace / "output/input_generation"
        (root / "species_cds").mkdir(parents=True)
        cds = root / "species_cds" / (sp + ".fa")
        cds.write_text(">" + sp + "_g1\nATGAAATAA\n")
        genome = root / "species_genome" / (sp + ".fa")
        genome.parent.mkdir()
        with genome.open("wb") as f:
            f.write(b">chr1\n")
            for _ in range(64):
                f.write(b"ACGT" * (1024 * 1024 // 4))
        row = dict(
            species_prefix=sp,
            provider="local",
            species_key=sp,
            cds_output_path=str(cds),
            gff_output_path="",
            genome_output_path=str(genome),
        )
        plan = root / "tmp/task_plan.json"
        atomic_json(plan, dict(version=2, task_count=1, tasks=[dict(species_prefix=sp)]))
        shard = root / "tmp/species_summary_shards/1.tsv"
        m.write_table(shard, [row])
        atomic_json(
            Path(str(plan) + ".completed/1.json"),
            dict(
                task_index=1,
                species_prefix=sp,
                plan_sha256=digest(plan),
                files={str(p): digest(p) for p in [cds, genome, shard]},
            ),
        )
        rows.append(
            dict(
                species=sp,
                source_workspace=str(workspace),
                plan_path=str(plan),
                plan_sha256=digest(plan),
                cds_path=str(cds),
                gff_path="",
                genome_path=str(genome),
            )
        )
    selection = tmp / "selection.tsv"
    m.write_table(selection, rows)
    plan = tmp / "plan.json"
    atomic_json(plan, m.freeze(selection, 2))
    view = tmp / "view"
    if args.phase == 'verify':
        m.prepare(view, plan)
    results = {name: [] for name in ["baseline", "candidate"]}
    for repeat in range(4):
        for name, path in [("baseline", BASE), ("candidate", CANDIDATE)]:
            target = view if args.phase == 'verify' else tmp / ('view-' + name + '-' + str(repeat))
            output = subprocess.check_output(
                [sys.executable, str(Path(__file__).resolve()), "--worker", str(path), str(target), str(plan), args.phase], text=True
            )
            value = json.loads(output)
            if repeat:
                results[name].append(value)
    all_rows = [v for group in results.values() for v in group]
    assert len({v["output_sha256"] for v in all_rows}) == 1
    report = dict(
        workload="Two 64 MiB FASTA payloads; three hardlinked names per genome; " + args.phase + ", no biological QC timing",
        runtime=sys.version,
        platform=sys.platform,
        repetitions=3,
        outputs_identical=True,
        results=results,
        medians={
            k: {
                metric: statistics.median(v[metric] for v in values)
                for metric in ["seconds", "peak_tracemalloc_bytes", "peak_rss_kib", "hash_reads", "hash_bytes"]
            }
            for k, values in results.items()
        },
    )
    print(json.dumps(report, indent=2))
