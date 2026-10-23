"""Compare PGLS correction and complete live/archive summary database builds.

Save the baseline helper, then pass its path with --before. Runs use alternating
fresh workers; fixture generation and output validation are outside each timer.
Invoke this script through run_in_runtime.sh, without simultaneous workloads.
"""

import argparse
import ast
import contextlib
import hashlib
import importlib.util
import json
import os
import platform
import sqlite3
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
SUPPORT = REPO_ROOT / "workflow/support"
sys.path.insert(0, str(SUPPORT))


def source_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def make_frame(families, wide=True, failed=False, many_pairs=False):
    """48 represented fits/family, using actual producer diagnostic columns."""
    positions = np.arange(families * 48)
    family = positions // 48
    frame = pd.DataFrame({
        "tree_id": [f"OG{i:07d}" for i in family],
        "analysis_method": np.where(positions % 48 < 24, "rsc", "species_nwkit"),
        "aggregation": np.asarray(["sum", "mean", "max"])[positions // 8 % 3],
        "analysis_id": "a", "model_id": "m",
        "response": [f"response_{i % 4}" + (f"_{f}" if many_pairs else "")
                     for i, f in zip(positions // 2, family, strict=True)],
        "source_term": np.where(positions % 2, "temperature", "precipitation"),
        "term": np.where(positions % 2, "temperature", "precipitation"),
        "term_test": "coefficient",
        "p_value": (positions % 997 + 1) / 998,
        "inference_status": np.where((positions % 19 == 0) | failed, "failed", "ok"),
    })
    frame.loc[positions % 101 == 0, "term"] = "(Intercept)"
    for divisor, value in ((211, np.nan), (223, np.inf), (227, -0.1), (229, 1.1)):
        frame.loc[positions % divisor == 0, "p_value"] = value
    if wide:
        tree = ast.parse((SUPPORT / "reconciled_speciation_contrast.py").read_text())
        columns = next(ast.literal_eval(node.value) for node in tree.body
                       if isinstance(node, ast.Assign)
                       and any(isinstance(target, ast.Name) and target.id == "STAT_TREE_RESULT_COLUMNS"
                               for target in node.targets))
        extras = {}
        for number, column in enumerate(columns):
            if column in frame:
                continue
            if column.endswith("parameter_status"):
                extras[column] = "estimated"
            elif "optimizer_converged" in column:
                extras[column] = "yes"
            elif column.endswith(("_model", "_message", "_warning", "_policy", "_mode", "_method")):
                extras[column] = "diagnostic"
            else:
                extras[column] = (positions + number) / 7
        frame = pd.concat([frame, pd.DataFrame(extras, index=frame.index)], axis=1)
    return frame.sample(frac=1, random_state=17)


def load_helper(source):
    spec = importlib.util.spec_from_file_location("pgls_multiplicity", source)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def peak_rss_kib():
    return int(next(line.split()[1] for line in Path("/proc/self/status").read_text().splitlines()
                    if line.startswith("VmHWM:")))


def database_fingerprint(path):
    digest = hashlib.sha256()
    frames = {}
    with sqlite3.connect(path) as conn:
        schema = conn.execute("SELECT type, name, sql FROM sqlite_master ORDER BY type, name").fetchall()
        digest.update(json.dumps(schema, ensure_ascii=False).encode())
        for kind, name, _ in schema:
            if kind != "table":
                continue
            quoted = '"' + name.replace('"', '""') + '"'
            columns = conn.execute(f"PRAGMA table_info({quoted})").fetchall()
            digest.update(json.dumps(columns, ensure_ascii=False).encode())
            order = ",".join('"' + column[1].replace('"', '""') + '"' for column in columns)
            # SQLite relations have no implicit row order. The existing threaded
            # tree/branch importer may insert completed files in either order.
            frame = pd.read_sql_query(f"SELECT * FROM {quoted} ORDER BY {order}", conn)
            digest.update(frame.to_csv(sep="\t", index=False).encode())
            frames[name] = frame
    return digest.hexdigest(), frames


def worker(args):
    helper = load_helper(args.source)
    root = args.worker
    os.chdir(root)
    if args.case.startswith("db-"):
        import generate_orthogroup_database as database

        database.write_association_table = helper.write_association_table
        target = root / "result.sqlite3"
        sys.argv = [str(SUPPORT / "generate_orthogroup_database.py"), "--overwrite", "1",
                    "--dbpath", str(target), "--dir_gene_family", str(args.fixture),
                    "--dir_stat_tree", str(args.fixture / "stat_tree"),
                    "--dir_stat_branch", str(args.fixture / "stat_branch"), "--ncpu", "1"]
        started, cpu_started = time.perf_counter(), time.process_time()
        database.main()
        elapsed, cpu = time.perf_counter() - started, time.process_time() - cpu_started
        peak = peak_rss_kib()
        output_hash, frames = database_fingerprint(target)
        pd.to_pickle(frames, root / "output.pkl")
    else:
        profile = {"small": (32, True, False, False), "normal": (1000, True, False, False),
                   "large": (4000, True, False, False), "narrow": (1000, False, False, False),
                   "failed": (1000, True, True, False), "many-pairs": (1000, True, False, True)}[args.case]
        frame = make_frame(*profile)
        original = frame.copy(deep=True)
        started, cpu_started = time.perf_counter(), time.process_time()
        output = helper.adjust_associations(frame)
        elapsed, cpu = time.perf_counter() - started, time.process_time() - cpu_started
        peak = peak_rss_kib()
        pd.testing.assert_frame_equal(frame, original, check_exact=True)
        assert output.index.equals(pd.RangeIndex(len(output)))
        output.to_pickle(root / "output.pkl")
        output_hash = hashlib.sha256(output.to_csv(sep="\t", index=False).encode()).hexdigest()
    (root / "measurement.json").write_text(json.dumps({
        "wall_seconds": elapsed, "cpu_seconds": cpu, "peak_rss_kib": peak,
        "output_sha256": output_hash,
    }) + "\n")


def prepare_database(root, families, archived):
    from gene_family_output_store import archive_completed_outputs, family_context

    frame = make_frame(families)
    for family, group in frame.groupby("tree_id", sort=False):
        tables = [("pgls_comparison", "_comparison.tsv", group),
                  ("stat_tree", "_stat.tree.tsv", pd.DataFrame([dict(num_branch=3, num_spe=1, num_dup=0, num_sp=2)])),
                  ("stat_branch", "_stat.branch.tsv", pd.DataFrame([dict(branch_id=0, node_name="n0", num_sp=2, so_event="S")]))]
        for directory, suffix, table in tables:
            path = root / directory / (family + suffix)
            path.parent.mkdir(parents=True, exist_ok=True)
            table.to_csv(path, sep="\t", index=False)
        plot = root / "tree_plot" / (family + "_tree_plot.pdf")
        plot.parent.mkdir(exist_ok=True)
        plot.write_bytes(b"%PDF-test\n")
    if archived:
        counts = root / "counts.tsv"
        ids = sorted(frame.tree_id.unique())
        pd.DataFrame({"Orthogroup": ids, "Total": 1}).to_csv(counts, sep="\t", index=False)
        families, extractor = family_context("orthogroup", genecount=counts)
        with open(os.devnull, "w") as quiet, contextlib.redirect_stdout(quiet):
            archive_completed_outputs(root, "orthogroup", families, extractor)
        assert not list((root / "pgls_comparison").glob("*.tsv"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before", type=Path)
    parser.add_argument("--after", type=Path, default=SUPPORT / "pgls_multiplicity.py")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--cases", nargs="+", default=["small", "normal", "large", "narrow", "failed", "many-pairs", "db-live", "db-archive"])
    parser.add_argument("--pairs", type=int, default=5)
    parser.add_argument("--worker", type=Path)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--case")
    parser.add_argument("--fixture", type=Path)
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    if not args.before or not args.output:
        parser.error("--before and --output are required")
    args.output = args.output.resolve()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    sources = {"before": args.before.resolve(), "after": args.after.resolve()}
    result = {"source_sha256": {name: source_hash(path) for name, path in sources.items()},
              "benchmark_sha256": source_hash(Path(__file__)),
              "environment": {"python": sys.version, "platform": platform.platform(),
                              "pandas": pd.__version__, "numpy": np.__version__},
              "diagnostics": "actual STAT_TREE_RESULT_COLUMNS; 48 fits/family, 2 methods, 3 aggregations, 4 responses, 2 predictors",
              "pairs": args.pairs, "warmup_pairs": 1, "cases": {}}
    with tempfile.TemporaryDirectory(prefix="gg-pgls-benchmark-") as temporary:
        temporary = Path(temporary)
        for case in args.cases:
            fixture = temporary / f"fixture-{case}"
            if case.startswith("db-"):
                prepare_database(fixture, 256, case == "db-archive")
            measured = []
            for pair in range(args.pairs + 1):
                measurements = {}
                for name in (("before", "after") if pair % 2 == 0 else ("after", "before")):
                    root = temporary / f"{case}-{pair}-{name}"
                    root.mkdir()
                    log = args.output.with_name(f"{args.output.stem}-{case}-{pair}-{name}.log")
                    with log.open("w") as handle:
                        subprocess.run([sys.executable, str(Path(__file__).resolve()), "--worker", str(root),
                                        "--source", str(sources[name]), "--case", case, "--fixture", str(fixture)],
                                       check=True, stdout=handle, stderr=subprocess.STDOUT)
                    measurements[name] = json.loads((root / "measurement.json").read_text())
                before = pd.read_pickle(temporary / f"{case}-{pair}-before/output.pkl")
                after = pd.read_pickle(temporary / f"{case}-{pair}-after/output.pkl")
                if isinstance(before, dict):
                    assert before.keys() == after.keys()
                    for key in before:
                        pd.testing.assert_frame_equal(before[key], after[key], check_exact=True)
                else:
                    pd.testing.assert_frame_equal(before, after, check_exact=True)
                assert measurements["before"]["output_sha256"] == measurements["after"]["output_sha256"]
                if pair:
                    measured.append(measurements)
                result["cases"][case] = measured
                args.output.write_text(json.dumps(result, indent=2) + "\n")
                print(json.dumps({"case": case, "pair": pair, "warmup": not pair, **measurements}), flush=True)


if __name__ == "__main__":
    main()
