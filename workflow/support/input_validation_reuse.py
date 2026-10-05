"""Reuse native species QC only after fresh content verification of its proofs."""

import csv
import hashlib
import importlib.metadata
import importlib.util
import json
import sys
from functools import lru_cache
from pathlib import Path

from input_generation_array_state import digest_paths, load_plan, receipt_path
from input_generation_stage_resume import checkpoint_path, context, parameters


@lru_cache(maxsize=1)
def implementation_identity():
    """Bind QC to support code and the installed sequence/GFF implementations."""
    sources = []
    support = Path(__file__).resolve().parent
    sources.extend(("support/" + str(path.relative_to(support)), path) for path in support.rglob("*.py"))
    versions = {"python": list(sys.version_info[:3])}
    for package in ("kffractbias", "Bio", "pysam"):
        spec = importlib.util.find_spec(package)
        if spec is None:
            versions[package] = "unavailable"
            continue
        for directory in spec.submodule_search_locations or ():
            root = Path(directory)
            for path in root.rglob("*"):
                if path.is_file() and path.suffix in {".py", ".so"}:
                    sources.append((package + "/" + str(path.relative_to(root)), path))
    for package in ("biopython", "pysam", "kffractbias", "pandas", "polars"):
        try:
            versions[package + "_version"] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package + "_version"] = "unavailable"
    # Source location/mount names do not change the identity. Compiled sequence
    # code and installed GFF parser content do; different runtimes revalidate.
    from input_generation_array_state import digest
    files = {name: digest(path) for name, path in sorted(sources)}
    return hashlib.sha256(json.dumps({"files": files, "versions": versions}, sort_keys=True).encode()).hexdigest()


def add_arguments(parser):
    parser.add_argument("--reuse-validation-root", type=Path,
                        help="Native input-generation root containing verified worker QC.")
    parser.add_argument("--reuse-task-plan", type=Path,
                        help="Immutable native plan for --reuse-validation-root.")
    parser.add_argument("--format-contract-version", help="Expected formatter contract for QC reuse.")


def canonical(path):
    return str(Path(path).resolve(strict=True))


def clean_row(row):
    return {key: str(value or "") for key, value in row.items() if not key.startswith("_")}


def reuse_one(task, plan_path, root, contract, stage, index, plan, missing_limit):
    """Verify all sources/outputs together, without reusing old content digests."""
    documents = {}

    def document(path):
        path = canonical(path)
        raw = Path(path).read_bytes()
        documents[path] = hashlib.sha256(raw).hexdigest()
        return json.loads(raw)

    receipt = document(receipt_path(plan_path, index))
    if (receipt.get("task_index") != index or receipt.get("species_prefix") != task["species_prefix"]
            or not isinstance(receipt.get("files"), dict) or not receipt["files"]):
        return None
    expected = {}

    def expect(path, sha):
        path = canonical(path)
        if not isinstance(sha, str) or (path in expected and expected[path] != sha):
            raise ValueError("Conflicting QC source identity")
        expected[path] = sha

    for path, sha in receipt["files"].items():
        expect(path, sha)
    format_task, settings, meta, format_paths = context(plan_path, index, root, "format")
    _, _, _, validation_paths = context(plan_path, index, root, "validate")
    if format_task != plan["tasks"][index - 1] or settings.get("run_validate_inputs") != "1":
        return None
    if stage == "mapping" and str(int(task.get("strict", False))) != str(settings.get("strict", "0")):
        return None
    # The plan/settings and worker metadata are not accepted merely because
    # their filenames or the output files still exist.
    for path in (root / "tmp/task_meta_shards" / f"{index}.json",
                 Path(str(plan_path) + ".settings.json")):
        if canonical(path) not in expected:
            return None
    for label, paths in (("format", format_paths), ("validate", validation_paths)):
        saved = document(checkpoint_path(root, task["species_prefix"], label))
        if (saved.get("schema_version") != 1 or saved.get("stage") != label
                or saved.get("species") != task["species_prefix"]
                or saved.get("parameters") != parameters(settings, label, contract)
                or set(saved.get("files", {})) != set(paths)):
            return None
        for role, path in paths.items():
            item = saved["files"][role]
            if canonical(item["path"]) != canonical(path) or canonical(path) not in expected:
                return None
            expect(path, item["sha256"])
    if canonical(task["cds_file"]) != canonical(meta["cds_output_path"]):
        return None
    if stage == "mapping":
        if canonical(task["gff_file"]) != canonical(meta["gff_output_path"]):
            return None
        genome = meta.get("genome_output_path")
        if bool(genome) != bool(task.get("genome_file")):
            return None
        if genome and canonical(task["genome_file"]) != canonical(genome):
            return None
    with Path(validation_paths["summary"]).open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1 or clean_row(rows[0]) != clean_row(task["summary_row"]):
        return None
    qc_path = validation_paths["mapping_qc" if stage == "mapping" else "ownership_qc"]
    qc = document(qc_path)
    if qc.get("species_checked") != 1 or qc.get("species_passed") != 1:
        return None
    if qc.get("validation_options") != {"missing_limit": missing_limit}:
        return None
    if qc.get("validation_implementation") != implementation_identity():
        return None
    result = qc.get("species_results", {}).get(task["species_prefix"])
    if (not isinstance(result, dict) or result.get("ok") is not True
            or result.get("species_prefix") != task["species_prefix"]
            or result.get("stats_ready") is not True or not isinstance(result.get("stats"), dict)
            or not isinstance(result.get("message"), str)):
        return None
    plan_sha = receipt.get("plan_sha256")
    expect(plan_path, plan_sha)
    for path, sha in documents.items():
        expect(path, sha)
    # Each large file is read only once in this verification boundary. The
    # batch fences every file and symlink target against changes during reads.
    observed = digest_paths(expected)
    if observed != expected:
        return None
    return {**result, "index": task["index"]}


def partition(tasks, args, parser, stage):
    supplied = (args.reuse_validation_root, args.reuse_task_plan, args.format_contract_version)
    if not any(supplied):
        return [], tasks
    if not all(supplied) or not args.species_summary:
        parser.error("QC reuse requires --reuse-validation-root, --reuse-task-plan, "
                     "--format-contract-version and --species-summary together")
    try:
        plan_path = args.reuse_task_plan.resolve(strict=True)
        root = args.reuse_validation_root.resolve(strict=True)
        plan = load_plan(plan_path)
        indices = {task["species_prefix"]: index for index, task in enumerate(plan["tasks"], 1)}
    except (OSError, ValueError, KeyError, TypeError):
        return [], tasks
    reused, pending = [], []
    for task in tasks:
        # Mapping tasks assembled from directory scans have no source summary
        # and therefore cannot inherit native source-ownership proofs.
        candidate = {**task, "species_prefix": task.get("species_prefix")
                     or task.get("summary_row", {}).get("species_prefix")}
        result = None
        try:
            if candidate["species_prefix"] in indices and "summary_row" in candidate:
                result = reuse_one(candidate, plan_path, root, args.format_contract_version,
                                   stage, indices[candidate["species_prefix"]], plan, args.missing_limit)
        except (OSError, ValueError, KeyError, TypeError, IndexError, AttributeError):
            pass
        if result is None:
            pending.append(task)
        else:
            print("Reused verified {} QC: {}".format(stage, candidate["species_prefix"]))
            reused.append(result)
    return reused, pending
