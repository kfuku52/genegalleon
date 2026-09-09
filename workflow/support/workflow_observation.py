#!/usr/bin/env python3
"""Optional, scheduler-neutral run observations. Never an execution authority."""
from __future__ import annotations

import argparse
import contextlib
import contextvars
import hashlib
import json
import os
import signal
import stat
import subprocess
import sys
import time
import uuid
from pathlib import Path

SCHEMA = "genegalleon-observation-v1"
_OBSERVED_PATHS = contextvars.ContextVar("genegalleon_observed_paths", default=None)


def observing_files():
    return _OBSERVED_PATHS.get() is not None


def path_signature(path):
    try:
        value = Path(path).lstat()
    except FileNotFoundError:
        return None
    return value.st_dev, value.st_ino, value.st_mode, value.st_size, value.st_mtime_ns, value.st_ctime_ns


def observe_path(path):
    paths = _OBSERVED_PATHS.get()
    if paths is None:
        return
    path = Path(os.path.abspath(path))
    current = path_signature(path)
    if paths.setdefault(path, current) != current:
        raise ValueError(f"Path changed during observation: {path}")


@contextlib.contextmanager
def observe_files():
    """Detect changes to observed paths without writing locks or rescanning data."""
    if _OBSERVED_PATHS.get() is not None:
        yield
        return
    paths = {}
    token = _OBSERVED_PATHS.set(paths)
    try:
        yield
        for path, expected in paths.items():
            if path_signature(path) != expected:
                raise ValueError(f"Path changed during observation; retry: {path}")
    finally:
        _OBSERVED_PATHS.reset(token)


def read_regular_json(path, limit=8 * 1024 * 1024):
    path = Path(path)
    observe_path(path)
    descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    with os.fdopen(descriptor, "rb") as handle:
        before = os.fstat(handle.fileno())
        if not stat.S_ISREG(before.st_mode):
            raise ValueError(f"JSON must be a regular file: {path}")
        raw = handle.read(limit + 1)
        after = os.fstat(handle.fileno())
    def signature(value):
        return value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns
    if signature(before) != signature(after) or signature(after) != signature(path.lstat()):
        raise ValueError(f"JSON changed during observation: {path}")
    if len(raw) > limit:
        raise ValueError(f"JSON exceeds {limit} bytes: {path}")
    return strict_json_loads(raw)


def strict_json_loads(raw):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f"Duplicate JSON key: {key}")
            result[key] = value
        return result

    def constant(value):
        raise ValueError(f"Non-finite JSON value: {value}")

    return json.loads(raw, object_pairs_hook=pairs, parse_constant=constant)


def json_digest(payload):
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def atomic_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    directory_fd = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    temporary = "." + path.name + "." + uuid.uuid4().hex
    try:
        descriptor = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600, dir_fd=directory_fd)
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, sort_keys=True, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path.name, src_dir_fd=directory_fd, dst_dir_fd=directory_fd)
        os.fsync(directory_fd)
    finally:
        try:
            os.unlink(temporary, dir_fd=directory_fd)
        except FileNotFoundError:
            pass
        os.close(directory_fd)


def report_error(code, *, step, affected_input=None, retryability="unknown", detail=None):
    """Called at the owning failure boundary; never classify arbitrary log text."""
    directory = os.environ.get("GG_OBSERVATION_ATTEMPT_DIR")
    if not directory:
        return
    try:
        atomic_json(Path(directory) / ("error-" + uuid.uuid4().hex + ".json"), {
            "schema": SCHEMA, "observed_at_ns": time.time_ns(), "step": step,
            "error_code": code, "affected_input": affected_input,
            "retryability": retryability, "detail": detail,
        })
    except (OSError, ValueError) as exc:
        print(f"GeneGalleon error observation unavailable: {exc}", file=sys.stderr)


def contract_arguments(args):
    # Export exact declared inputs/outputs/parameters, using the same public
    # provenance arguments. No shell evaluation is required to inspect later.
    scalar = ("manifest", "step", "family_id", "logical_root", "workspace_root")
    repeated = ("input", "input_gene_family_store", "input_gene_family_subdir", "input_gene_family_artifact",
                "input_logical_directory", "output", "output_logical_directory", "optional_output",
                "output_fasta_type", "recover_output", "parameter")
    argv = []
    for name in scalar:
        value = str(getattr(args, name))
        if name in {"manifest", "logical_root", "workspace_root"}:
            value = os.path.abspath(value)
        argv.extend(["--" + name.replace("_", "-"), value])
    for name in repeated:
        for value in getattr(args, name):
            if name not in {"parameter", "output_fasta_type"}:
                label, raw = value.split("=", 1)
                parts = raw.split("::") if name in {"input_gene_family_subdir", "input_gene_family_artifact"} else [raw]
                parts[0] = os.path.abspath(parts[0])
                value = label + "=" + "::".join(parts)
            argv.extend(["--" + name.replace("_", "-"), value])
    argv.extend(["--stale-policy", getattr(args, "stale_policy", "stop")])
    return argv


def record_contract_result(args, result):
    directory = os.environ.get("GG_OBSERVATION_ATTEMPT_DIR")
    if not directory or getattr(args, "dry_run", False):
        return
    try:
        try:
            argv = contract_arguments(args)
        except (ValueError, TypeError):
            argv = []
        manifest_digest = None
        try:
            manifest_digest = json_digest(read_regular_json(args.manifest))
        except (OSError, ValueError):
            pass  # Still replace a previous successful receipt with this result.
        key = json_digest([args.family_id, args.step])
        atomic_json(Path(directory) / ("contract-" + key + ".json"), {
            "schema": "genegalleon-contract-observation-v1", "attempt_id": Path(directory).name,
            "observed_at_ns": time.time_ns(), "family_id": args.family_id, "step": args.step,
            "operation": args.command, "exit_code": result, "argv": argv,
            "workspace_root": os.path.abspath(args.workspace_root),
            "manifest_sha256": manifest_digest,
        })
    except (OSError, ValueError) as exc:
        print(f"GeneGalleon contract observation unavailable: {exc}", file=sys.stderr)


def runtime_manifest(workflow, config_keys, workflow_root, version_override=None):
    """Capture only registered, resolved settings, never the entire environment."""
    root = Path(workflow_root)
    version = root.parent / "VERSION"
    files = [root / "gg_common_params.sh", root / "core" / (workflow + "_core.sh"),
             root / (workflow + "_entrypoint.sh")]
    fingerprints = {}
    for path in files:
        if path.is_file():
            fingerprints[str(path.relative_to(root))] = hashlib.sha256(path.read_bytes()).hexdigest()
    config = {key: os.environ[key] for key in config_keys if key in os.environ}
    return {
        "schema": "genegalleon-runtime-v1", "workflow": workflow,
        "genegalleon_version": version_override or (version.read_text().strip() if version.is_file() else None),
        "workflow_files_sha256": fingerprints,
        "effective_config_sha256": json_digest(config),
        "effective_config": config,
        "config_scope": "registered entrypoint settings forwarded to the core; core-local derived defaults excluded",
        "unset_config_keys": sorted(set(config_keys) - config.keys()),
        "runtime_id": os.environ.get("GG_RESOURCE_RUNTIME_ID"),
        "container_digest": os.environ.get("GG_OBSERVATION_CONTAINER_DIGEST"),
        "container_digest_source": "caller-supplied; not independently verified",
        "scheduler": {key: os.environ.get(key) for key in (
            "GG_JOB_ID", "GG_ARRAY_JOB_ID", "GG_ARRAY_TASK_ID", "GG_TASK_CPUS", "GG_MEM_TOTAL_GB")},
    }


def run(args):
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        raise ValueError("a command is required")
    if any(not 0 <= code <= 255 for code in args.accept_exit_code):
        raise ValueError("accepted exit codes must be between 0 and 255")
    attempt = uuid.uuid4().hex
    directory = args.directory / attempt
    record = {"schema": SCHEMA, "attempt_id": attempt, "workflow": args.workflow,
              "started_at_ns": time.time_ns(), "execution_state": "started",
              "completion_state": "unverified", "exit_code": None,
              "error_code": None, "retryability": "unknown", "liveness": "unknown",
              "accepted_exit_codes": sorted({0, *args.accept_exit_code})}
    env = dict(os.environ)
    script = sys.stdin.buffer.read() if args.stdin_script else None
    writable = False
    try:
        record["runtime"] = runtime_manifest(args.workflow, args.config_key, args.workflow_root, args.genegalleon_version)
        record["runtime"]["executed_core_sha256"] = hashlib.sha256(script).hexdigest() if script is not None else None
        atomic_json(directory / "run.json", record)
        env["GG_OBSERVATION_ATTEMPT_DIR"] = str(directory.absolute())
        writable = True
    except (OSError, ValueError) as exc:
        env.pop("GG_OBSERVATION_ATTEMPT_DIR", None)
        print(f"GeneGalleon run observation unavailable: {exc}", file=sys.stderr)
    # Preserve inherited streams, descriptors (including workflow/scratch locks),
    # and the child's exit status. Forward termination to the complete subtree.
    handlers = {}
    child = None
    pending = []

    def forward(signum, _frame):
        if child is None:
            pending.append(signum)
        else:
            try:
                os.killpg(child.pid, signum)
            except ProcessLookupError:
                pass

    try:
        for signum in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
            handlers[signum] = signal.signal(signum, forward)
        try:
            child = subprocess.Popen(command, env=env, start_new_session=True, close_fds=False,
                                     stdin=subprocess.PIPE if script is not None else None)
            for signum in pending:
                forward(signum, None)
            if script is not None:
                child.communicate(input=script)
                result = child.returncode
            else:
                result = child.wait()
            result = 128 - result if result < 0 else result
        except OSError as exc:
            result = 127 if isinstance(exc, FileNotFoundError) else 126
            record["error_code"] = "command_start_failed"
            print(f"GeneGalleon command start failed: {exc}", file=sys.stderr)
    finally:
        for signum, handler in handlers.items():
            signal.signal(signum, handler)
    accepted = result in record["accepted_exit_codes"]
    record.update(finished_at_ns=time.time_ns(), exit_code=result, execution_accepted=accepted,
                  execution_state="exited", error_code=record["error_code"] or (None if accepted else "unknown"))
    if writable:
        try:
            atomic_json(directory / "run.json", record)
        except (OSError, ValueError) as exc:
            print(f"GeneGalleon final observation unavailable: {exc}", file=sys.stderr)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--workflow", required=True)
    parser.add_argument("--workflow-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--config-key", action="append", default=[])
    parser.add_argument("--genegalleon-version", help="version read by the submitting entrypoint")
    parser.add_argument("--stdin-script", action="store_true", help="fingerprint and forward the exact core script")
    parser.add_argument("--accept-exit-code", action="append", default=[], type=int,
                        help="additional existing success status; never rewrites the process exit code")
    parser.add_argument("command", nargs=argparse.REMAINDER)
    return run(parser.parse_args())


if __name__ == "__main__":
    sys.exit(main())
