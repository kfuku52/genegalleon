#!/usr/bin/env python3
"""Measure lossless rescue storage on frozen, complete producer publications.

Run through the normal runtime wrapper. Inputs are read-only; only the supplied
canonical benchmark directory and its separately owned scratch are written.
Execution FASTA/GFF files are inventoried, never copied. First-process and warm
reads do not imply a cold operating-system cache.
"""
import argparse
import hashlib
import json
import os
import platform
import re
import resource
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from collections import Counter
from contextlib import contextmanager
from pathlib import Path, PurePosixPath

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
from input_generation_array_state import atomic_json  # noqa: E402
from rescue_prediction_cache import stream_json_array  # noqa: E402

MEMBERS = {"models": "models.json", "partial": "partial_models.json",
           "revision": "revision_candidates.json"}
SHA = re.compile(r"[0-9a-f]{64}\Z")
QC_FIELDS = ("problems", "cds", "sequence", "frameshift", "coverage", "identity",
             "assembly_ambiguous_bases", "terminal_completion", "partial_evidence")


def identity(path):
    info = Path(path).stat()
    return [info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns]


def digest(path):
    path = Path(path)
    before = identity(path)
    result = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(block)
    if identity(path) != before:
        raise OSError("File changed while hashing: " + str(path))
    return result.hexdigest()


def member(root, name):
    if (not isinstance(name, str) or not name or "\\" in name
            or any(ord(c) < 32 for c in name)):
        raise ValueError("Unsafe producer member")
    relative = PurePosixPath(name)
    if relative.is_absolute() or str(relative) != name or any(p in {".", ".."} for p in relative.parts):
        raise ValueError("Unsafe producer member")
    path = Path(root)
    for part in relative.parts:
        path = path / part
        if path.is_symlink():
            raise ValueError("Symlink producer member")
    if not path.is_file():
        raise ValueError("Missing producer member: " + name)
    return path


def inventory(root):
    """Metadata only, no symlink traversal; charge each physical file once."""
    pending, seen = [Path(root)], set()
    result = {"bytes": 0, "allocated_bytes": 0, "regular_files": 0,
              "directories": 0, "symlinks": 0, "hardlink_aliases": 0}
    while pending:
        directory = pending.pop()
        with os.scandir(directory) as entries:
            for entry in entries:
                if entry.is_symlink():
                    result["symlinks"] += 1
                elif entry.is_dir(follow_symlinks=False):
                    result["directories"] += 1
                    pending.append(Path(entry.path))
                elif entry.is_file(follow_symlinks=False):
                    info = entry.stat(follow_symlinks=False)
                    result["regular_files"] += 1
                    token = info.st_dev, info.st_ino
                    if token in seen:
                        result["hardlink_aliases"] += 1
                    else:
                        seen.add(token)
                        result["bytes"] += info.st_size
                        result["allocated_bytes"] += info.st_blocks * 512
    return result


def freeze_inputs(root, species, *, expected_plan, expected_receipt):
    root = Path(root).resolve(strict=True)
    if not re.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.-]*", species) or species in {".", ".."}:
        raise ValueError("Unsafe species")
    if not SHA.fullmatch(str(expected_plan)) or not SHA.fullmatch(str(expected_receipt)):
        raise ValueError("Independent expected plan and producer receipt digests are required")
    directory = root / "rescued" / species
    plan, receipt_path = member(root, "plan.json"), member(directory, "receipt.json")
    if digest(plan) != expected_plan or digest(receipt_path) != expected_receipt:
        raise ValueError("Frozen producer plan/receipt differs")
    receipt = json.loads(receipt_path.read_text())
    if receipt.get("key", {}).get("plan") != expected_plan or receipt["key"].get("species") != species:
        raise ValueError("Producer receipt belongs to another plan/species")
    files = {}
    for kind, name in MEMBERS.items():
        if name not in receipt["files"]:
            if kind == "models" or (directory / name).exists() or (directory / name).is_symlink():
                raise ValueError("Unbound producer model member: " + name)
            files[kind] = None
            continue
        path = member(directory, name)
        expected = receipt["files"][name]
        if not SHA.fullmatch(str(expected)):
            raise ValueError("Invalid producer member checksum")
        files[kind] = {"path": str(path), "resolved_path": str(path.resolve(strict=True)),
                       "sha256": expected, "identity": identity(path)}
    frozen = {"root": str(root), "species": species, "directory": str(directory), "files": files,
              "plan": {"path": str(plan), "sha256": expected_plan, "identity": identity(plan)},
              "receipt": {"path": str(receipt_path), "sha256": expected_receipt, "identity": identity(receipt_path)},
              "baseline_capacity": inventory(directory),
              "legacy_model_bytes": sum(spec["identity"][2] for spec in files.values() if spec is not None),
              "baseline_capacity_scope": "metadata-only whole worker; execution inputs not copied"}
    check_inputs(frozen)
    return frozen


def check_inputs(frozen):
    for spec in [frozen["plan"], frozen["receipt"], *frozen["files"].values()]:
        if spec is None:
            continue
        path = Path(spec["path"])
        if identity(path) != spec["identity"] or ("resolved_path" in spec and str(path.resolve(strict=True)) != spec["resolved_path"]):
            raise OSError("Frozen input changed: " + str(path))
    for key in ("plan", "receipt"):
        if digest(frozen[key]["path"]) != frozen[key]["sha256"]:
            raise ValueError("Frozen input metadata checksum differs")


def encoded(value, *, sort=True):
    return json.dumps(value, sort_keys=sort, ensure_ascii=False, allow_nan=False,
                      separators=(",", ":")).encode("utf-8") + b"\n"


def configuration_sha(config):
    return hashlib.sha256(encoded({k: v for k, v in config.items() if k != "frozen_sha256"})).hexdigest()


def process_io():
    path = Path("/proc/self/io")
    if not path.exists():
        return None
    return {key: int(value) for key, value in (line.split(":", 1) for line in path.read_text().splitlines())}


def io_difference(before, after):
    return None if before is None or after is None else {key: after[key] - before[key] for key in before}


def owned_input_copies(config, destination):
    """Copy only model arrays to newly owned inodes; verify every source byte."""
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=False)
    start, cpu, before_io = time.perf_counter(), time.process_time(), process_io()
    check_inputs(config["inputs"])
    copied = json.loads(json.dumps(config["inputs"]))
    proof = {}
    for kind, spec in config["inputs"]["files"].items():
        if spec is None:
            continue
        original, target = Path(spec["path"]), destination / MEMBERS[kind]
        original_identity = identity(original)
        value = hashlib.sha256()
        with original.open("rb") as source, target.open("xb") as output:
            for block in iter(lambda: source.read(1024 * 1024), b""):
                value.update(block)
                output.write(block)
            output.flush()
            os.fsync(output.fileno())
        if identity(original) != original_identity or value.hexdigest() != spec["sha256"] or digest(target) != spec["sha256"]:
            raise ValueError("Owned benchmark input copy differs")
        target.chmod(0o444)
        token = identity(target)
        if token[:2] == original_identity[:2]:
            raise ValueError("Owned benchmark copy aliases a producer inode")
        copied["files"][kind] = {**spec, "path": str(target), "resolved_path": str(target.resolve(strict=True)), "identity": token}
        proof[kind] = {"source": str(original), "source_identity": original_identity,
                       "destination": str(target.resolve(strict=True)), "identity": token,
                       "sha256": spec["sha256"], "bytes": token[2]}
    check_inputs(config["inputs"])
    receipt = {"schema": 1, "owner_frozen_sha256": config["frozen_sha256"], "files": proof,
               "wall_seconds": time.perf_counter() - start, "cpu_seconds": time.process_time() - cpu,
               "io_delta": io_difference(before_io, process_io()), "capacity": inventory(destination)}
    atomic_json(destination / ".benchmark-copy-owner.json", receipt, immutable=True)
    return copied, receipt


def advise_owned_copies(destination, expected_owner):
    """DONTNEED is advisory, inode-specific, and forbidden on producer files."""
    destination = Path(destination).resolve(strict=True)
    receipt = json.loads((destination / ".benchmark-copy-owner.json").read_text())
    if receipt.get("owner_frozen_sha256") != expected_owner:
        raise ValueError("Benchmark-owned input copy owner differs")
    if not hasattr(os, "posix_fadvise") or not hasattr(os, "POSIX_FADV_DONTNEED"):
        raise RuntimeError("Owned-inode POSIX advice is unavailable")
    advised = []
    for spec in receipt["files"].values():
        path = Path(spec["destination"])
        if (path.is_symlink() or path.parent != destination or identity(path) != spec["identity"]
                or spec["identity"][:2] == spec["source_identity"][:2] or digest(path) != spec["sha256"]):
            raise ValueError("Benchmark copy ownership/content changed before advice")
        descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
        try:
            info = os.fstat(descriptor)
            if [info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns] != spec["identity"]:
                raise ValueError("Benchmark copy inode changed before advice")
            os.posix_fadvise(descriptor, 0, 0, os.POSIX_FADV_DONTNEED)
        finally:
            os.close(descriptor)
        advised.append({"path": str(path), "device_inode": spec["identity"][:2], "bytes": spec["bytes"]})
    return {"advice": "POSIX_FADV_DONTNEED", "owned_inodes_only": advised,
            "cache_regime": "owned copied-input DONTNEED advice attempted; operating-system cold cache is not guaranteed"}


class Summary:
    """Constant-space fingerprints retain values, record order and field order."""
    def __init__(self):
        self.count = 0
        self.canonical = hashlib.sha256()
        self.ordered = hashlib.sha256()
        self.qc = hashlib.sha256()
        self.accepted = hashlib.sha256()
        self.accepted_count = 0
        self.statuses = Counter()

    def add(self, record):
        if not isinstance(record, dict):
            raise ValueError("Model record is not an object")
        canonical = encoded(record)
        self.canonical.update(canonical)
        self.ordered.update(encoded(record, sort=False))
        self.qc.update(encoded({key: record[key] for key in QC_FIELDS if key in record}))
        status = record.get("status")
        self.statuses[str(status)] += 1
        if status == "accepted":
            self.accepted.update(canonical)
            self.accepted_count += 1
        self.count += 1

    def result(self):
        return {"records": self.count, "canonical_sha256": self.canonical.hexdigest(),
                "ordered_fields_sha256": self.ordered.hexdigest(), "qc_sequence_cds_sha256": self.qc.hexdigest(),
                "accepted_records": self.accepted_count, "accepted_full_fields_sha256": self.accepted.hexdigest(),
                "status_counts": dict(self.statuses)}


class SourceRows:
    def __init__(self, spec):
        self.spec, self.summary, self.complete = spec, Summary(), False
        self.raw_sha = hashlib.sha256()

    def __iter__(self):
        if self.spec is None:
            self.complete = True
            return
        path = Path(self.spec["path"])
        if identity(path) != self.spec["identity"]:
            raise OSError("Frozen model input changed")
        reader = stream_json_array(path, hasher=self.raw_sha)
        try:
            for record in reader:
                self.summary.add(record)
                yield record
            if self.raw_sha.hexdigest() != self.spec["sha256"]:
                raise ValueError("Frozen model input checksum differs")
            if identity(path) != self.spec["identity"]:
                raise OSError("Frozen model input changed while reading")
            self.complete = True
        finally:
            reader.close()


def scan_rows(rows):
    summary = Summary()
    for record in rows:
        summary.add(record)
    return summary.result()


class DiskSampler:
    def __init__(self, root, period=1.0):
        self.roots = [Path(p) for p in root] if isinstance(root, (list, tuple)) else [Path(root)]
        self.period = period
        self.maximum = {"bytes": 0, "allocated_bytes": 0, "regular_files": 0}
        self.stop = threading.Event()
        self.error = None

    def sample(self):
        try:
            value = {key: 0 for key in self.maximum}
            for root in self.roots:
                if root.exists():
                    current = inventory(root)
                    for key in value:
                        value[key] += current[key]
            for key in self.maximum:
                self.maximum[key] = max(self.maximum[key], value[key])
        except FileNotFoundError:
            pass  # Own atomic rename/temporary cleanup; next sample includes target.
        except OSError as error:
            self.error = str(error)

    def loop(self):
        while not self.stop.wait(self.period):
            self.sample()

    def __enter__(self):
        self.sample()
        self.thread = threading.Thread(target=self.loop, daemon=True)
        self.thread.start()
        return self

    def __exit__(self, *args):
        self.stop.set()
        self.thread.join()
        self.sample()


class PlannedInterruption(RuntimeError):
    pass


@contextmanager
def owned_temporary(directory):
    owned = Path(directory) / ".benchmark-tmp"
    owned.mkdir(exist_ok=True)
    previous_env, previous_dir = os.environ.get("TMPDIR"), tempfile.tempdir
    os.environ["TMPDIR"], tempfile.tempdir = str(owned), str(owned)
    try:
        yield
    finally:
        if previous_env is None:
            os.environ.pop("TMPDIR", None)
        else:
            os.environ["TMPDIR"] = previous_env
        tempfile.tempdir = previous_dir


def interrupt_after(rows, count):
    for index, record in enumerate(rows):
        if index == count:
            raise PlannedInterruption("Owned benchmark interruption")
        yield record
    raise PlannedInterruption("Owned benchmark interruption at input end")


def stamp_store(directory, config):
    store = Path(directory) / "model_store"
    files = {p.relative_to(directory).as_posix(): digest(p) for p in sorted(store.rglob("*")) if p.is_file()}
    atomic_json(Path(directory) / "receipt.json", {"key": {"plan": config["inputs"]["plan"]["sha256"],
                "species": config["inputs"]["species"], "benchmark_frozen_sha256": config["frozen_sha256"]},
                "files": files}, immutable=True)


def worker(config, mode, directory):
    import rescue_model_store as store
    if config["frozen_sha256"] != configuration_sha(config):
        raise ValueError("Frozen benchmark configuration changed")
    if (digest(Path(store.__file__)) != config["module_sha256"]
            or digest(Path(__file__)) != config["driver_sha256"]):
        raise ValueError("Benchmark implementation changed after freeze")
    check_inputs(config["inputs"])
    if "original_inputs" in config:
        check_inputs(config["original_inputs"])
    start, cpu = time.perf_counter(), time.process_time()
    before_io = process_io()
    before_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result = {"mode": mode, "cache_regime": "independent process; operating-system cold cache not guaranteed"}
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    # Partial readers may spool refs/envelopes in SQLite. Charge that temporary
    # work to this child, rather than the running dataset's shared TMPDIR.
    with owned_temporary(directory), DiskSampler(directory) as disk:
        if mode == "legacy_read":
            sources = {kind: SourceRows(spec) for kind, spec in config["inputs"]["files"].items()}
            for rows in sources.values():
                for _ in rows:
                    pass
            result["summaries"] = {kind: rows.summary.result() for kind, rows in sources.items()}
        elif mode in {"write", "restart"}:
            if mode == "restart":
                rows = SourceRows(config["inputs"]["files"]["models"])
                try:
                    store.write_model_store(directory, interrupt_after(rows, config["restart_after"]),
                                            codec=config["codec"], shard_bytes=config["shard_bytes"])
                except PlannedInterruption:
                    if (directory / "model_store").exists() or list(directory.glob(".model-store-*")):
                        raise RuntimeError("Interrupted model-store publication survived") from None
                else:
                    raise RuntimeError("Expected owned interruption did not occur")
                result["interrupted_prefix_records"] = rows.summary.count
                result["interrupted_publication_absent"] = True
            sources = {kind: SourceRows(spec) for kind, spec in config["inputs"]["files"].items()}
            manifest = store.write_model_store(directory, iter(sources["models"]),
                partial_models=iter(sources["partial"]), revisions=iter(sources["revision"]),
                codec=config["codec"], shard_bytes=config["shard_bytes"])
            if not all(rows.complete for rows in sources.values()):
                raise RuntimeError("Model-store writer did not exhaust all frozen source records")
            stamp_store(directory, config)
            result["summaries"] = {kind: rows.summary.result() for kind, rows in sources.items()}
            result["manifest_counts"] = manifest["counts"]
            result["codec"] = manifest["codec"]
        elif mode == "store_read":
            readers = {"models": store.iter_models, "partial": store.iter_partial_models,
                       "revision": store.iter_revision_models, "accepted": store.iter_accepted_models}
            result["summaries"], result["frozen_keys"] = {}, {}
            for kind, reader in readers.items():
                key = store.frozen_model_store_key(directory, kind=kind)
                result["summaries"][kind] = scan_rows(reader(directory, frozen_key=key))
                result["frozen_keys"][kind] = key
        else:
            raise ValueError("Unknown benchmark worker mode")
    check_inputs(config["inputs"])
    if "original_inputs" in config:
        check_inputs(config["original_inputs"])
    result.update(wall_seconds=time.perf_counter() - start, cpu_seconds=time.process_time() - cpu,
                  io_delta=io_difference(before_io, process_io()),
                  peak_rss_kib_linux=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  initial_peak_rss_kib_linux=before_rss, maximum_owned_disk_sampled=disk.maximum,
                  disk_sample_period_seconds=disk.period, disk_sampling_error=disk.error,
                  final_capacity=inventory(directory))
    return result


def assert_parity(baseline, current):
    for kind in MEMBERS:
        if current["summaries"][kind] != baseline["summaries"][kind]:
            raise ValueError("Full ordered model/QC/sequence/provenance parity differs: " + kind)
    if "accepted" in current["summaries"]:
        accepted = current["summaries"]["accepted"]
        models = baseline["summaries"]["models"]
        if (accepted["records"] != models["accepted_records"]
                or accepted["canonical_sha256"] != models["accepted_full_fields_sha256"]):
            raise ValueError("Selective accepted path/isoform/provenance parity differs")


def publish_store(source, destination):
    """Never overwrite canonical output; copy only this owned compact worker."""
    source, destination = Path(source), Path(destination)
    if destination.exists() or destination.is_symlink():
        raise ValueError("Canonical benchmark store already exists")
    before = {p.relative_to(source).as_posix(): (identity(p), digest(p))
              for p in source.rglob("*") if p.is_file()}
    if any(p.is_symlink() for p in source.rglob("*")):
        raise ValueError("Unexpected benchmark publication symlink")
    temporary = Path(tempfile.mkdtemp(prefix=".publishing-" + destination.name + "-", dir=destination.parent))
    shutil.copytree(source, temporary, symlinks=False, dirs_exist_ok=True)
    for relative, (token, expected) in before.items():
        original, copied = source / relative, temporary / relative
        # A just-written NFS destination may finish a ctime update during its
        # first read. Only that field may differ, and every byte must still equal
        # the immutable source digest. Source identities remain fully strict.
        copied_before = identity(copied)
        value = hashlib.sha256()
        with copied.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                value.update(block)
        if (identity(original) != token or identity(copied)[:-1] != copied_before[:-1]
                or value.hexdigest() != expected):
            raise ValueError("Canonical compact publication differs")
    copied_members = {p.relative_to(temporary).as_posix() for p in temporary.rglob("*") if p.is_file()}
    if copied_members != set(before) or destination.exists() or destination.is_symlink():
        raise ValueError("Canonical compact publication generation differs")
    temporary.rename(destination)
    return {"files": len(before), "bytes": sum(token[2] for token, _ in before.values()),
            "member_sha256": {name: value[1] for name, value in before.items()}}


def run_child(config_path, mode, directory, output):
    command = [sys.executable, str(Path(__file__).resolve()), "--worker", mode,
               "--config", str(config_path), "--worker-directory", str(directory)]
    log = Path(output).with_suffix(".log")
    if Path(output).exists() or log.exists():
        raise ValueError("Benchmark phase evidence already exists; use a new namespace")
    with log.open("w") as handle:
        child = subprocess.Popen(command, stdout=handle, stderr=subprocess.STDOUT)
        atomic_json(Path(output).parent / "progress.json", {"mode": mode, "pid": child.pid,
                    "start_monotonic": time.monotonic(), "native_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                    "command": command, "state": "running"})
        returncode = child.wait()
    if returncode:
        raise RuntimeError(f"Storage benchmark {mode} failed; preserved diagnostics: {log}")
    result = json.loads(log.read_text().splitlines()[-1])
    atomic_json(output, result, immutable=True)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path)
    parser.add_argument("--species")
    parser.add_argument("--frozen-consumer-plan", type=Path)
    parser.add_argument("--expected-source-plan-sha256")
    parser.add_argument("--expected-source-receipt-sha256")
    parser.add_argument("--output-directory", type=Path)
    parser.add_argument("--scratch-directory", type=Path)
    parser.add_argument("--codec", choices=("auto", "gzip", "zstd"), default="auto")
    parser.add_argument("--shard-bytes", type=int, default=8 * 1024 * 1024)
    parser.add_argument("--repeats", type=int, default=2)
    parser.add_argument("--restart-after", type=int, default=128)
    parser.add_argument("--freeze-only", action="store_true")
    parser.add_argument("--owned-copy-cold-advice", action="store_true",
                        help="Measure advice/first and warmed reads of newly created, fully verified model-array copies only")
    parser.add_argument("--config", type=Path)
    parser.add_argument("--worker", choices=("legacy_read", "write", "restart", "store_read"))
    parser.add_argument("--worker-directory", type=Path)
    args = parser.parse_args(argv)
    if args.worker:
        config = json.loads(args.config.read_text())
        print(json.dumps(worker(config, args.worker, args.worker_directory), sort_keys=True))
        return
    if not all((args.source_root, args.species, args.output_directory, args.scratch_directory)):
        parser.error("source root/species and distinct canonical output/scratch directories are required")
    if args.repeats < 1 or args.restart_after < 1 or args.shard_bytes < 1:
        parser.error("positive repeats/restart-after/shard-bytes required")
    expected_plan, expected_receipt = args.expected_source_plan_sha256, args.expected_source_receipt_sha256
    consumer = None
    if args.frozen_consumer_plan:
        consumer = {"path": str(args.frozen_consumer_plan.resolve(strict=True)), "sha256": digest(args.frozen_consumer_plan)}
        frozen_cache = json.loads(args.frozen_consumer_plan.read_text())["request"]["prediction_cache"]
        expected_plan = frozen_cache["plan_sha256"]
        expected_receipt = frozen_cache["species"][args.species]["receipt_sha256"]
        if args.source_root.resolve(strict=True) != Path(frozen_cache["root"]).resolve(strict=True):
            raise ValueError("Frozen consumer cache root differs")
    inputs = freeze_inputs(args.source_root, args.species, expected_plan=expected_plan, expected_receipt=expected_receipt)
    source, output, scratch = Path(inputs["root"]), args.output_directory.resolve(), args.scratch_directory.resolve()
    if (source == output or source in output.parents or source == scratch or source in scratch.parents
            or output == scratch or output in scratch.parents or scratch in output.parents):
        raise ValueError("Benchmark destinations overlap immutable input or each other")
    output.mkdir(parents=True, exist_ok=True)
    scratch.mkdir(parents=True, exist_ok=True)
    import rescue_model_store as store
    config = {"schema": 1, "inputs": inputs, "consumer": consumer, "codec": args.codec,
              "shard_bytes": args.shard_bytes, "restart_after": args.restart_after,
              "owned_copy_cold_advice": args.owned_copy_cold_advice,
              "module_sha256": digest(Path(store.__file__)), "driver_sha256": digest(Path(__file__)),
              "python": sys.version, "platform": platform.platform()}
    config["frozen_sha256"] = configuration_sha(config)
    config_path = output / "frozen_inputs.json"
    atomic_json(config_path, config, immutable=True)
    atomic_json(scratch / ".benchmark-owned.json", {"frozen_sha256": config["frozen_sha256"], "canonical": str(output)}, immutable=True)
    if args.freeze_only:
        print(json.dumps({"config": str(config_path), "baseline_capacity": inputs["baseline_capacity"]}, indent=2))
        return
    if (output / "result.json").exists() or (scratch / "worker" / "model_store").exists():
        raise ValueError("Completed benchmark exists; use a new measurement namespace")
    with DiskSampler([scratch, output]) as total_disk:
        copy_proof, advice, warmed = None, None, None
        if args.owned_copy_cold_advice:
            effective_inputs, copy_proof = owned_input_copies(config, scratch / "owned_input_copies")
            atomic_json(output / "owned_input_copy.json", copy_proof, immutable=True)
            advice = advise_owned_copies(scratch / "owned_input_copies", config["frozen_sha256"])
            atomic_json(output / "owned_input_advice.json", advice, immutable=True)
            config = {**config, "original_inputs": inputs, "inputs": effective_inputs,
                      "original_frozen_sha256": config["frozen_sha256"], "copy_receipt": copy_proof}
            config["frozen_sha256"] = configuration_sha(config)
            config_path = output / "advised_input_config.json"
            atomic_json(config_path, config, immutable=True)
        baseline = run_child(config_path, "legacy_read", scratch / "legacy_read", output / "legacy_read.json")
        if advice is not None:
            baseline["cache_regime"] = advice["cache_regime"]
            warmed = run_child(config_path, "legacy_read", scratch / "legacy_read", output / "legacy_read_warmed.json")
            assert_parity(baseline, warmed)
            warmed["cache_regime"] = "repeated read of same benchmark-owned copied inputs; warm candidate"
        creation = run_child(config_path, "restart", scratch / "worker", output / "creation.json")
        assert_parity(baseline, creation)
        reads = []
        for index in range(args.repeats):
            value = run_child(config_path, "store_read", scratch / "worker", output / f"store_read_{index}.json")
            assert_parity(baseline, value)
            value["cache_regime"] = "first independent process" if not index else "repeated independent process; warm candidate"
            reads.append(value)
        copied = publish_store(scratch / "worker", output / "canonical_worker")
        canonical = run_child(config_path, "store_read", output / "canonical_worker", output / "canonical_read.json")
        assert_parity(baseline, canonical)
    check_inputs(inputs)
    report = {"schema": 1, "native_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
              "frozen_sha256": config["frozen_sha256"], "legacy": baseline, "creation": creation,
              "owned_input_copy": copy_proof, "owned_input_advice": advice, "legacy_warmed": warmed,
              "reads": reads, "canonical_read": canonical, "canonical_publication": copied,
              "baseline_capacity": inputs["baseline_capacity"], "all_record_values_and_orders_equal": True,
              "legacy_model_stream_bytes": inputs["legacy_model_bytes"],
              "maximum_total_owned_disk_sampled": total_disk.maximum,
              "total_disk_sampling_error": total_disk.error,
              "body_duplication_evidence": {"models": creation["manifest_counts"]["models"],
                    "distinct_storage_bodies": creation["manifest_counts"]["bodies"],
                    "interpretation": "Exact storage-body sharing only; not a DNA-QC memo hit rate or locality measurement."},
              "scientific_proof": "Exact full raw/checked records, QC flags, CDS/sequences, paths, isoform/support and representative provenance preserved; accepted selective reader independently equal.",
              "consumer_pipeline_and_busco_reexecution": "Not performed by storage-only benchmark; root/consumer integration validates complete pipeline and identical BUSCO inputs separately.",
              "limitations": "Legacy capacity is existing producer inventory; legacy timings measure read/verification, not historical generation. Writer timing includes input audits and controlled prefix restart. Disk peak is sampled; first/warm labels do not guarantee OS cold cache."}
    atomic_json(output / "result.json", report, immutable=True)
    print(json.dumps({"result": str(output / "result.json"), "codec": creation["codec"],
                      "all_record_values_and_orders_equal": True}, indent=2))


if __name__ == "__main__":
    main()
