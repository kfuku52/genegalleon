#!/usr/bin/env python3
"""Per-family native CAFE rate LRT with a parametric bootstrap.

H0: one family-specific lambda on the whole tree.
H1: background lambda plus a separate lambda on one nominated branch.
All likelihoods, reconstructions and simulations are computed by cafe5. The
observed change sign comes from the input CAFE reconstruction, not from lambda.
"""
from __future__ import annotations

import argparse
import csv
from contextlib import contextmanager
import fcntl
import hashlib
import io
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from typing import Sequence

from Bio import Phylo

CONTRACT = "native_cafe_base_family_lrt_bootstrap_v1"
SCORE_TOLERANCE = 2e-5


class NativeCafeError(RuntimeError):
    pass


class OutputBusyError(NativeCafeError):
    pass


@contextmanager
def output_lock(directory: Path):
    # Keep one writer for the complete result set, including its cache and audit.
    with (directory / ".writer.lock").open("a") as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise OutputBusyError(f"Another native CAFE analysis owns {directory}") from exc
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def sha256(path: Path) -> str:
    with path.open("rb") as handle:
        h = hashlib.sha256()
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def atomic_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")
    temp.replace(path)


def write_table(path: Path, rows: list[dict], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    with temp.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({k: "NA" if r.get(k) is None else r.get(k) for k in fields} for r in rows)
    temp.replace(path)


def read_table(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if not fields or len(fields) != len(set(fields)):
            raise ValueError(f"Missing or duplicate header in {path}")
        rows = list(reader)
    if any(None in row or any(v is None for v in row.values()) for row in rows):
        raise ValueError(f"Ragged table: {path}")
    return fields, rows


def integer(value: str, label: str, signed: bool = False) -> int:
    if not re.fullmatch(r"-?\d+" if signed else r"\d+", value):
        raise ValueError(f"Expected an integer for {label}: {value!r}")
    return int(value)


@dataclass
class CafeTree:
    newick: str
    lambda_newick: str
    species: list[str]
    labels: list[str]
    root_label: str
    target_label: str
    parent_label: str
    tip_labels: dict[str, str]
    root_sides: list[list[str]]
    branch_map: list[dict]


def read_cafe_tree(path: Path, target: str, dated_tree: Path | None = None) -> CafeTree:
    # CAFE writes one fully labelled tree per family. Stream only the first;
    # a generic Nexus loader would load every reconstructed family into memory.
    pieces = []
    with path.open() as handle:
        for line in handle:
            if not pieces:
                match = re.match(r"\s*TREE\s+\S+\s*=\s*(.*)", line, re.I)
                if match:
                    pieces.append(match[1])
            else:
                pieces.append(line)
            if pieces and ";" in pieces[-1]:
                break
    if not pieces or not "".join(pieces).strip().endswith(";"):
        raise ValueError("Expected a complete CAFE TREE statement in the ASR Nexus file.")
    tree = Phylo.read(io.StringIO("".join(pieces)), "newick")
    labels = {}
    species = {}
    parents = {}
    for node in tree.find_clades(order="preorder"):
        name = node.name or ""
        match = re.fullmatch(r"(.*<\d+>)\*?_(\d+)", name)
        if not match:
            raise ValueError(f"Unsupported CAFE ASR node label: {name!r}")
        label = match[1]
        if label in labels.values():
            raise ValueError(f"Duplicate CAFE branch label: {label}")
        labels[node] = label
        if node.is_terminal():
            sp = re.sub(r"<\d+>$", "", label)
            # CAFE's native Newick reader cannot safely round-trip quoted names.
            if not re.fullmatch(r"[A-Za-z0-9_.-]+", sp):
                raise ValueError(f"Unsupported native CAFE species name: {sp!r}")
            if sp in species.values():
                raise ValueError(f"Duplicate species in CAFE ASR tree: {sp}")
            species[node] = sp
        elif len(node.clades) != 2:
            raise ValueError("The native CAFE comparison requires a bifurcating tree.")
        for child in node.clades:
            parents[child] = node
        if node is not tree.root and (node.branch_length is None or
                not math.isfinite(node.branch_length) or node.branch_length <= 0):
            raise ValueError(f"CAFE branch {label} requires a positive finite duration.")
    if dated_tree:
        original = Phylo.read(dated_tree, "newick")
        original_nodes = {}
        for node in original.find_clades():
            key = frozenset(t.name for t in node.get_terminals())
            if key in original_nodes or (not node.is_terminal() and len(node.clades) != 2):
                raise ValueError("The dated tree must be bifurcating with unique species/clades.")
            original_nodes[key] = node
        asr_keys = {frozenset(species[t] for t in n.get_terminals()) for n in labels}
        if asr_keys != set(original_nodes):
            raise ValueError("Dated tree topology/species do not match the native CAFE ASR tree.")
        for node in labels:
            if node is tree.root:
                continue
            length = original_nodes[frozenset(species[t] for t in node.get_terminals())].branch_length
            if length is None or not math.isfinite(length) or length <= 0 or not math.isclose(
                    length, node.branch_length, rel_tol=1e-5, abs_tol=1e-12):
                raise ValueError("Dated tree durations disagree with the native CAFE ASR tree.")
            # Use original precision, retaining the CAFE ordering/branch identity.
            node.branch_length = length
    targets = [n for n, label in labels.items() if label == target]
    if len(targets) != 1 or targets[0] is tree.root:
        raise ValueError("Target must exactly name a non-root CAFE branch.")
    nominated = targets[0]

    def render(node, lambda_tree=False):
        body = species[node] if node.is_terminal() else "(" + ",".join(render(c, lambda_tree) for c in node.clades) + ")"
        if lambda_tree:
            return body + (":2" if node is nominated else ":1")
        return body + ("" if node is tree.root else f":{node.branch_length:.17g}")

    mapping = [{"branch": label, "parent": labels[parents[n]] if n in parents else "",
                "length": n.branch_length if n is not tree.root else None,
                "lambda_group": 2 if n is nominated else 1,
                "descendant_species": ",".join(species[t] for t in n.get_terminals())}
               for n, label in labels.items()]
    return CafeTree(render(tree.root) + ";\n", render(tree.root, True) + ";\n",
                    list(species.values()), list(labels.values()), labels[tree.root],
                    target, labels[parents[nominated]],
                    {sp: labels[n] for n, sp in species.items()},
                    [[species[t] for t in c.get_terminals()] for c in tree.root.clades], mapping)


@dataclass
class Fit:
    nll: float
    lambdas: list[float]
    iterations: int
    run_directory: str


def parse_fit(directory: Path, expected_lambdas: int, max_iterations: int) -> Fit:
    log = (directory / "cafe.log").read_text()
    scores = re.findall(r"^Final -lnL:\s*(\S+)\s*$", log, re.M)
    iterations = re.findall(r"^Completed\s+(\d+)\s+iterations\s*$", log, re.M)
    lines = (directory / "results" / "Base_lambda_per_family.txt").read_text().strip().splitlines()
    if len(scores) != 1 or len(iterations) != 1 or len(lines) != 1:
        raise NativeCafeError(f"Incomplete native single-family optimization output: {directory}")
    parts = lines[0].split("\t")
    if len(parts) != 2 or parts[0] != "family":
        raise NativeCafeError(f"Unexpected native family ID in {directory}")
    lambdas = [float(x.strip()) for x in parts[1].split(",")]
    nll, niter = float(scores[0]), int(iterations[0])
    if len(lambdas) != expected_lambdas or not math.isfinite(nll) or any(
            not math.isfinite(x) or x < 0 for x in lambdas):
        raise NativeCafeError(f"Invalid native likelihood/lambda in {directory}")
    if niter >= max_iterations:
        raise NativeCafeError(f"Native optimization hit the iteration limit in {directory}")
    return Fit(nll, lambdas, niter, str(directory))


def lr_statistic(null: Fit, alternative: Fit) -> float:
    tolerance = SCORE_TOLERANCE * max(1., abs(null.nll), abs(alternative.nll))
    improvement = null.nll - alternative.nll
    if improvement < -tolerance:
        raise NativeCafeError("The fitted alternative is worse than its nested null; optimization is unresolved.")
    return max(0., 2 * improvement)


def bootstrap_pvalue(observed: float, replicates: Sequence[float]) -> float:
    if not math.isfinite(observed) or observed < 0 or not replicates or any(
            not math.isfinite(x) or x < 0 for x in replicates):
        raise ValueError("Bootstrap requires a finite statistic and all requested finite replicate statistics.")
    tolerance = SCORE_TOLERANCE * max(1., observed)
    return (1 + sum(x >= observed - tolerance for x in replicates)) / (len(replicates) + 1)


class NativeCafe:
    def __init__(self, tree: CafeTree, workdir: Path, executable: str, restarts: int,
                 max_iterations: int, cores: int, timeout: int, error_model: Path | None = None):
        resolved = shutil.which(executable)
        if not resolved:
            raise NativeCafeError(f"Native CAFE executable not found: {executable}")
        self.executable = Path(resolved).resolve()
        self.tree = tree
        self.workdir = workdir.resolve()
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.restarts, self.max_iterations, self.cores, self.timeout = restarts, max_iterations, cores, timeout
        self.identity = {"contract": CONTRACT, "cafe_sha256": sha256(self.executable),
                         "adapter_sha256": sha256(Path(__file__).resolve()),
                         "tree": tree.newick, "lambda_tree": tree.lambda_newick,
                         "restarts": restarts, "max_iterations": max_iterations, "cores": cores,
                         "error_model_sha256": sha256(error_model) if error_model else None}
        # Each identity has its own cache. No old-method or changed-input reuse.
        key = hashlib.sha256(json.dumps(self.identity, sort_keys=True).encode()).hexdigest()
        self.cache = self.workdir / "runs" / key
        self.cache.mkdir(parents=True, exist_ok=True)
        self.treefile = self.cache / "tree.nwk"
        self.lambdafile = self.cache / "lambda_tree.nwk"
        self.treefile.write_text(tree.newick)
        self.lambdafile.write_text(tree.lambda_newick)
        self.errorfile = None
        if error_model:
            self.errorfile = self.cache / "error_model.txt"
            shutil.copyfile(error_model, self.errorfile)
        atomic_json(self.cache / "identity.json", self.identity)

    def command(self, directory: Path, args: list[str]) -> None:
        command = [str(self.executable), *args, "-t", str(self.treefile), "-c", str(self.cores),
                   "-o", str(directory / "results")]
        if self.errorfile:
            command.append("-e" + str(self.errorfile))
        atomic_json(directory / "command.json", {"argv": command, **self.identity})
        with (directory / "cafe.log").open("w") as log:
            try:
                completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                                           cwd=directory, timeout=self.timeout, check=False)
            except subprocess.TimeoutExpired as exc:
                raise NativeCafeError(f"Native CAFE timed out; see {directory / 'cafe.log'}") from exc
        if completed.returncode:
            raise NativeCafeError(f"Native CAFE failed with status {completed.returncode}; see {directory / 'cafe.log'}")

    def input_file(self, directory: Path, counts: Sequence[int]) -> Path:
        path = directory / "family.tsv"
        write_table(path, [dict(zip(["Desc", "FamilyID", *self.tree.species], ["(null)", "family", *counts]))],
                    ["Desc", "FamilyID", *self.tree.species])
        return path

    def cached_run(self, key: dict, action, required: Sequence[str]) -> Path:
        required = [*required, "request.json", "family.tsv"]
        digest = hashlib.sha256(json.dumps(key, sort_keys=True).encode()).hexdigest()
        parent = self.cache / digest
        parent.mkdir(exist_ok=True)
        manifest = parent / "complete.json"
        if manifest.exists():
            try:
                saved = json.loads(manifest.read_text())
                name = saved["directory"]
                hashes = saved["sha256"]
                if (not isinstance(name, str) or Path(name).name != name or
                        not name.startswith("attempt_") or not isinstance(hashes, dict)):
                    raise ValueError("Invalid cache completion record")
                directory = parent / name
                if all((directory / file).is_file() and sha256(directory / file) == hashes.get(file)
                       for file in required):
                    return directory
            except (ValueError, KeyError, TypeError) as exc:
                raise NativeCafeError(f"Malformed native cache completion record: {manifest}") from exc
            raise NativeCafeError(f"Cached native outputs changed: {parent}; use a new output directory.")
        directory = Path(tempfile.mkdtemp(prefix="attempt_", dir=parent))
        atomic_json(directory / "request.json", key)
        action(directory)
        if any(not (directory / file).is_file() for file in required):
            raise NativeCafeError(f"Native CAFE did not write required outputs: {directory}")
        atomic_json(manifest, {"directory": directory.name,
                              "sha256": {file: sha256(directory / file) for file in required}})
        return directory

    def fit(self, counts: Sequence[int], alternative: bool) -> Fit:
        fits = []
        for restart in range(self.restarts):
            def run(directory):
                path = self.input_file(directory, counts)
                args = ["-b", "-i", str(path), f"-I{self.max_iterations}"]
                if alternative:
                    args += ["-y", str(self.lambdafile)]
                self.command(directory, args)
                parse_fit(directory, 2 if alternative else 1, self.max_iterations)
            directory = self.cached_run({"action": "fit", "counts": list(counts),
                                         "alternative": alternative, "restart": restart}, run,
                                        ["cafe.log", "results/Base_lambda_per_family.txt", "command.json"])
            fits.append(parse_fit(directory, 2 if alternative else 1, self.max_iterations))
        fits.sort(key=lambda x: x.nll)
        tolerance = SCORE_TOLERANCE * max(1., abs(fits[0].nll))
        if len(fits) < 2 or fits[1].nll - fits[0].nll > tolerance:
            raise NativeCafeError("Independent native optimizations disagree; inspect native CAFE convergence and fit settings before rerunning.")
        return fits[0]

    def null_root(self, counts: Sequence[int], fit: Fit) -> tuple[int, str]:
        if fit.lambdas[0] <= 0:
            raise NativeCafeError("Native CAFE cannot simulate a nonpositive fitted lambda for a nonzero LRT.")
        def run(directory):
            path = self.input_file(directory, counts)
            self.command(directory, ["-i", str(path), "-l", format(fit.lambdas[0], ".17g")])
        directory = self.cached_run({"action": "null_reconstruction", "counts": list(counts),
                                     "lambda": fit.lambdas[0]}, run,
                                    ["cafe.log", "results/Base_count.tab", "results/Base_results.txt", "command.json"])
        _, rows = read_table(directory / "results/Base_count.tab")
        if len(rows) != 1 or rows[0].get("FamilyID") != "family" or self.tree.root_label not in rows[0]:
            raise NativeCafeError(f"Unexpected native null reconstruction: {directory}")
        # This reconstruction profiles ancestral counts at the fitted null rate;
        # bootstrap generation is conditional on this plug-in ancestral root.
        root = integer(rows[0][self.tree.root_label], "null reconstructed root")
        if root <= 0:
            raise NativeCafeError("A positive reconstructed root is required for native CAFE simulation.")
        scores = re.findall(r"^Score \(-lnL\):\s*(\S+)", (directory / "cafe.log").read_text(), re.M)
        if not scores or not math.isfinite(float(scores[-1])) or abs(float(scores[-1]) - fit.nll) > SCORE_TOLERANCE * max(1., abs(fit.nll)):
            raise NativeCafeError("Native fixed-null likelihood does not reproduce the optimized null.")
        return root, str(directory)

    def simulate(self, counts: Sequence[int], null: Fit, root: int, replicates: int) -> tuple[list[list[int]], str]:
        def run(directory):
            path = self.input_file(directory, counts)
            rootfile = directory / "root_distribution.tsv"
            rootfile.write_text(f"{root}\t{replicates}\n")
            self.command(directory, [f"-s{replicates}", "-i", str(path), "-l", format(null.lambdas[0], ".17g"),
                                     "-f", str(rootfile)])
        directory = self.cached_run({"action": "simulate", "counts": list(counts), "lambda": null.lambdas[0],
                                     "root": root, "replicates": replicates}, run,
                                    ["cafe.log", "results/simulation.txt", "results/simulation_truth.txt", "command.json",
                                     "root_distribution.tsv"])
        if "Failed to create a family" in (directory / "cafe.log").read_text():
            raise NativeCafeError("Native CAFE reported a failed simulation; no replicates may be discarded.")
        fields, rows = read_table(directory / "results/simulation.txt")
        if set(fields[2:]) != set(self.tree.species) or len(rows) != replicates:
            raise NativeCafeError("Native simulation species or replicate count does not match the request.")
        result = [[integer(row[sp], f"simulated {sp}") for sp in self.tree.species] for row in rows]
        for values in result:
            by_species = dict(zip(self.tree.species, values))
            if not all(any(by_species[sp] > 0 for sp in side) for side in self.tree.root_sides):
                raise NativeCafeError("Native simulation failed its root-presence condition; no selective replacement is allowed.")
        return result, str(directory)


def compare_family(engine: NativeCafe, counts: Sequence[int], replicates: int, audit: Path) -> dict:
    null, alternative = engine.fit(counts, False), engine.fit(counts, True)
    statistic = lr_statistic(null, alternative)
    record = {"null_nll": null.nll, "alternative_nll": alternative.nll,
              "lambda_null": null.lambdas[0], "lambda_background": alternative.lambdas[0],
              "lambda_target": alternative.lambdas[1], "lrt": statistic,
              "null_fit": null.run_directory, "alternative_fit": alternative.run_directory}
    atomic_json(audit / "observed.json", record)
    # Since every bootstrap statistic is >=0 by definition, T=0 has P=1
    # exactly. No native simulation or artificial pseudo-count is needed.
    if statistic == 0:
        return {**record, "p_value": 1., "bootstrap_completed": 0,
                "bootstrap_root": None, "p_method": "zero_statistic", "status": "tested"}
    root, root_path = engine.null_root(counts, null)
    samples, simulation_path = engine.simulate(counts, null, root, replicates)
    atomic_json(audit / "generation.json", {"root": root, "null_reconstruction": root_path,
                                            "simulation_directory": simulation_path})
    bootstrap_rows = []
    for i, sample in enumerate(samples):
        sim_null, sim_alt = engine.fit(sample, False), engine.fit(sample, True)
        value = lr_statistic(sim_null, sim_alt)
        bootstrap_rows.append({"replicate": i + 1, "null_nll": sim_null.nll, "alternative_nll": sim_alt.nll,
                               "lrt": value, "null_fit": sim_null.run_directory, "alternative_fit": sim_alt.run_directory})
        write_table(audit / "bootstrap.tsv", bootstrap_rows,
                    ["replicate", "null_nll", "alternative_nll", "lrt", "null_fit", "alternative_fit"])
        print(f"  bootstrap {i + 1}/{replicates}", flush=True)
    return {**record, "p_value": bootstrap_pvalue(statistic, [r["lrt"] for r in bootstrap_rows]),
            "bootstrap_completed": len(bootstrap_rows), "bootstrap_root": root,
            "p_method": "native_parametric_bootstrap", "status": "tested"}


RESULT_FIELDS = ["FamilyID", "target_change", "direction", "lambda_null", "lambda_background", "lambda_target",
                 "null_nll", "alternative_nll", "lrt", "p_value", "bootstrap_completed", "bootstrap_root",
                 "p_method", "status", "error", "audit_directory"]


def run(args) -> int:
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    with output_lock(output):
        try:
            return run_locked(args, output)
        except (ValueError, NativeCafeError, OSError) as exc:
            atomic_json(output / "metadata.json", {"contract": CONTRACT, "status": "failed", "error": str(exc)})
            print(f"ERROR: {exc}", file=sys.stderr)
            return 1


def run_locked(args, output: Path) -> int:
    write_table(output / "family_lrt.tsv", [], RESULT_FIELDS)
    tree = read_cafe_tree(args.asr_tree, args.target_branch, args.tree)
    headers, rows = read_table(args.counts)
    if set(headers) != {"FamilyID", *tree.labels}:
        raise ValueError("CAFE count-table columns do not match the ASR tree labels.")
    change_headers, changes = read_table(args.changes)
    if set(change_headers) != {"FamilyID", *tree.labels}:
        raise ValueError("CAFE change-table columns do not match the ASR tree labels.")
    def index(records):
        result = {}
        for row in records:
            key = row["FamilyID"]
            if not key or key in result:
                raise ValueError("CAFE FamilyID must be nonempty and unique.")
            result[key] = row
        return result
    by_id, change_by_id = index(rows), index(changes)
    family_headers, family_rows = read_table(args.family_ids)
    if family_headers != ["FamilyID"]:
        raise ValueError("The family selection table requires exactly one FamilyID column.")
    families = list(index(family_rows))
    if not families or not set(families) <= by_id.keys() or not set(families) <= change_by_id.keys():
        raise ValueError("Every requested family must occur in both CAFE tables.")
    engine = NativeCafe(tree, output, args.cafe, args.fit_restarts, args.max_iterations, args.cores,
                        args.timeout, args.error_model)
    write_table(output / "branch_map.tsv", tree.branch_map,
                ["branch", "parent", "length", "lambda_group", "descendant_species"])
    metadata = {**engine.identity, "inputs": {str(path.resolve()): sha256(path) for path in
                ([args.counts, args.changes, args.asr_tree, args.family_ids] + ([args.tree] if args.tree else []))},
                "target_branch": args.target_branch, "bootstrap_replicates": args.bootstrap_replicates,
                "minimum_bootstrap_p": 1 / (args.bootstrap_replicates + 1),
                "null_model": "one_lambda_per_family", "alternative_model": "background_and_target_lambda_per_family",
                "root_generation": "native_null_reconstruction_plugin",
                "direction_source": "input_CAFE_change_table", "cafe_rng": "native_unseeded_draws_and_fits_preserved",
                "family_order": families, "status": "running"}
    atomic_json(output / "metadata.json", metadata)
    results = []
    for i, family in enumerate(families):
        print(f"CAFE specificity: family {i + 1}/{len(families)} ({family})", flush=True)
        # Never interpolate externally supplied FamilyIDs into paths or shell code.
        audit = output / "families" / hashlib.sha256(family.encode()).hexdigest()
        audit.mkdir(parents=True, exist_ok=True)
        for name in ("observed.json", "generation.json", "bootstrap.tsv"):
            (audit / name).unlink(missing_ok=True)
        record = {"FamilyID": family, "audit_directory": str(audit), "status": "failed"}
        try:
            row = by_id[family]
            node_counts = {label: integer(row[label], f"{family}:{label}") for label in tree.labels}
            delta = integer(change_by_id[family][tree.target_label], f"{family} target change", signed=True)
            if delta != node_counts[tree.target_label] - node_counts[tree.parent_label]:
                raise ValueError("CAFE target change disagrees with its reconstructed parent/child counts.")
            record.update(target_change=delta, direction="increase" if delta > 0 else "decrease" if delta < 0 else "unchanged")
            counts = [node_counts[tree.tip_labels[sp]] for sp in tree.species]
            count_map = dict(zip(tree.species, counts))
            if not all(any(count_map[sp] > 0 for sp in side) for side in tree.root_sides):
                raise ValueError("Requested family does not satisfy native CAFE root-presence ascertainment.")
            record.update(compare_family(engine, counts, args.bootstrap_replicates, audit))
        except (ValueError, NativeCafeError, OSError) as exc:
            record.update(status="failed", p_value=None, error=str(exc))
        results.append(record)
        write_table(output / "family_lrt.tsv", results, RESULT_FIELDS)
        if record["status"] != "tested":
            # Fail closed. Do not redefine the GO background after a numerical
            # or simulation failure. Completed native jobs can be resumed.
            metadata.update(status="failed", failed_family=family, error=record["error"])
            atomic_json(output / "metadata.json", metadata)
            print(record["error"], file=sys.stderr)
            return 1
    metadata["status"] = "complete"
    atomic_json(output / "metadata.json", metadata)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ("counts", "changes", "asr-tree", "family-ids", "output-dir"):
        parser.add_argument("--" + option, type=Path, required=True)
    parser.add_argument("--tree", type=Path, help="Original dated tree; checked against CAFE ASR and used at full precision.")
    parser.add_argument("--target-branch", required=True)
    parser.add_argument("--cafe", default="cafe5")
    parser.add_argument("--bootstrap-replicates", type=int, default=999)
    parser.add_argument("--fit-restarts", type=int, default=5)
    parser.add_argument("--max-iterations", type=int, default=1000)
    parser.add_argument("--cores", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=3600, help="Seconds allowed per native CAFE invocation.")
    parser.add_argument("--error-model", type=Path)
    args = parser.parse_args()
    if args.bootstrap_replicates < 1 or args.fit_restarts < 2 or args.max_iterations < 1 or args.cores < 1 or args.timeout < 1:
        parser.error("Require bootstrap >= 1, fit restarts >= 2, and positive iterations, cores and timeout.")
    try:
        return run(args)
    except (ValueError, NativeCafeError, OSError) as exc:
        # A competing writer owns its metadata; never overwrite it on lock failure.
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
