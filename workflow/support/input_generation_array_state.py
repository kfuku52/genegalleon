"""Immutable array plans and verified, atomic per-worker completion receipts."""
import argparse
import csv
import hashlib
import io
import json
import os
import stat
import tempfile
from pathlib import Path


def digest(path):
    result = hashlib.sha256()
    with os.fdopen(os.open(path, os.O_RDONLY | os.O_NONBLOCK), "rb") as handle:
        before = os.fstat(handle.fileno())
        if not stat.S_ISREG(before.st_mode):
            raise ValueError("Expected a regular file for hashing: " + str(path))
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(chunk)
        after = os.fstat(handle.fileno())
        if (before.st_size, before.st_mtime_ns, before.st_ctime_ns) != (after.st_size, after.st_mtime_ns, after.st_ctime_ns):
            raise OSError("File changed while hashing: " + str(path))
    return result.hexdigest()


def atomic_json(path, value, immutable=False):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    data = json.dumps(value, sort_keys=True, indent=2) + "\n"
    with tempfile.NamedTemporaryFile(mode="w", dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
        handle.write(data)
    try:
        if immutable:
            try:
                os.link(temporary, path)
            except FileExistsError:
                if path.read_text() != data:
                    raise ValueError("Plan already exists with different inputs; use a new output workspace: " + str(path)) from None
        else:
            os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def safe_component(value):
    return (isinstance(value, str) and bool(value) and not value.startswith(".")
            and not any(char in value for char in ("/", "\\"))
            and not any(ord(char) < 32 or ord(char) == 127 for char in value))


def load_plan(path):
    plan = json.loads(Path(path).read_text())
    tasks = plan["tasks"]
    species = [task["species_prefix"] for task in tasks]
    if any(not safe_component(name) for name in species):
        raise ValueError("Task species prefixes must be safe, non-hidden filename components")
    if not tasks or len(set(species)) != len(tasks) or plan["task_count"] != len(tasks):
        raise ValueError("Invalid task plan: empty tasks, duplicate species, or inconsistent count")
    return plan


def receipt_path(plan, index):
    return Path(str(plan) + ".completed") / (str(index) + ".json")


def verify_receipt(plan, index, plan_sha256=None):
    try:
        receipt = json.loads(receipt_path(plan, index).read_text())
        if receipt["plan_sha256"] != (plan_sha256 or digest(plan)) or receipt["task_index"] != index:
            return False
        if not isinstance(receipt["files"], dict) or not receipt["files"]:
            return False
        return all(Path(path).is_file() and digest(path) == value for path, value in receipt["files"].items())
    except (OSError, ValueError, KeyError, TypeError, AttributeError):
        return False


def frozen_input_hashes(plan_path, plan, index):
    task = plan["tasks"][index - 1]
    expected = dict(task.get("input_sha256", {}))
    if "manifest_row" in task:
        cached = json.loads((Path(str(plan_path) + ".tasks") / (str(index) + ".json")).read_text())
        if cached.get("plan_sha256") != digest(plan_path) or cached.get("task_index") != index:
            raise ValueError("Resolved download cache belongs to another plan/task")
        expected.update(cached["task"]["input_sha256"])
    return expected


def export_manifest(plan, outfile):
    rows = [dict(task.get("manifest_row") or {"provider": task["provider"],
            "id": task["species_key"], "species_key": task["species_key"]}) for task in plan["tasks"]]
    fields = ["provider", "id", "species_key"] + sorted({key for row in rows for key in row} - {"provider", "id", "species_key"})
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
    path = Path(outfile)
    data = buffer.getvalue().encode("utf-8")
    if path.is_file() and path.read_bytes() == data:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as handle:
        temp = Path(handle.name)
        handle.write(data)
    try:
        os.replace(temp, path)
    finally:
        temp.unlink(missing_ok=True)


def prepared(plan):
    try:
        marker = json.loads(Path(str(plan) + ".prepared.json").read_text())
        if marker["plan_sha256"] != digest(plan) or marker["settings_sha256"] != digest(str(plan) + ".settings.json"):
            return False
        return all(digest(path) == value for path, value in marker.get("files", {}).items())
    except (OSError, ValueError, KeyError, TypeError, AttributeError):
        return False


def claim_workspace(plan, workspace, create=False, output_dirs=()):
    """A workspace's shard namespace belongs to exactly one immutable plan."""
    markers = [Path(workspace) / ".array-plan.json"]
    markers.extend(Path(str(Path(path).resolve()) + ".gg-input-generation-owner.json") for path in output_dirs)
    identity = {"task_plan": str(Path(plan).resolve()), "plan_sha256": digest(plan),
                "workspace": str(Path(workspace).resolve())}
    # All output locks are held by the core before checking or creating claims.
    for marker in markers:
        if marker.exists() and json.loads(marker.read_text()) != identity:
            raise ValueError("Output location belongs to another array plan: " + str(marker))
        if not create and not marker.is_file():
            raise ValueError("This output location is not prepared for this task plan: " + str(marker))
    if create:
        for marker in markers:
            atomic_json(marker, identity, immutable=True)


def output_lock_paths(paths):
    for path in sorted({str(Path(path).resolve()) for path in paths}):
        if "\n" in path:
            raise ValueError("Newlines are unsupported in array output directory paths")
        lock = Path(path + ".gg-input-generation.lock")
        lock.parent.mkdir(parents=True, exist_ok=True)
        yield str(lock)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("complete", "pending", "verify", "invalidate", "configure", "prepared", "check-prepared", "claim-workspace", "index", "export-manifest", "output-lock-paths"))
    parser.add_argument("--task-plan", required=True)
    parser.add_argument("--task-index", type=int)
    parser.add_argument("--file", action="append", default=[])
    parser.add_argument("--setting", action="append", default=[])
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--workspace")
    parser.add_argument("--outfile")
    args = parser.parse_args()
    if args.action == "output-lock-paths":
        print("\n".join(output_lock_paths(args.file)))
        return
    if args.action == "configure":
        path = Path(args.task_plan + ".settings.json")
        settings = dict(item.split("=", 1) for item in args.setting)
        settings["input_sha256"] = {str(Path(p).resolve()): digest(p) if Path(p).is_file() else None for p in args.file}
        if not args.prepare and not path.exists():
            parser.error("Array settings are missing; run array_prepare first")
        atomic_json(path, settings, immutable=True)
        return
    plan = load_plan(args.task_plan)
    if args.action == "export-manifest":
        if not args.outfile:
            parser.error("--outfile is required")
        export_manifest(plan, args.outfile)
        return
    if args.action == "claim-workspace":
        if not args.workspace:
            parser.error("--workspace is required")
        claim_workspace(args.task_plan, args.workspace, args.prepare, args.file)
        return
    if args.action == "prepared":
        atomic_json(args.task_plan + ".prepared.json", {"plan_sha256": digest(args.task_plan),
                    "settings_sha256": digest(args.task_plan + ".settings.json"),
                    "files": {str(Path(p).resolve()): digest(p) for p in args.file}})
        return
    if args.action == "check-prepared":
        raise SystemExit(0 if prepared(args.task_plan) else 1)
    if args.action == "pending":
        plan_sha256 = digest(args.task_plan)
        print(",".join(str(i) for i in range(1, plan["task_count"] + 1) if not verify_receipt(args.task_plan, i, plan_sha256)))
        return
    index = args.task_index
    if index is None or not 1 <= index <= plan["task_count"]:
        parser.error("A valid task index is required")
    if args.action == "index":
        print(index)
        return
    if args.action == "verify":
        raise SystemExit(0 if verify_receipt(args.task_plan, index) else 1)
    path = receipt_path(args.task_plan, index)
    if args.action == "invalidate":
        path.unlink(missing_ok=True)
        return
    if not args.file or any(not Path(p).is_file() or Path(p).stat().st_size == 0 for p in args.file):
        parser.error("Completion requires nonempty output files")
    expected_inputs = frozen_input_hashes(args.task_plan, plan, index)
    if any(digest(p) != expected for p, expected in expected_inputs.items()):
        parser.error("Raw inputs changed while the task was running")
    atomic_json(path, {"plan_sha256": digest(args.task_plan), "task_index": index,
                       "species_prefix": plan["tasks"][index-1]["species_prefix"],
                       "files": {**{str(Path(p).resolve()): digest(p) for p in args.file}, **expected_inputs}})


if __name__ == "__main__":
    main()
