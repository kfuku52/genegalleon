#!/usr/bin/env python3
"""Shared lock metadata helpers used by Python and shell workflow code."""

import argparse
import json
import os
import signal
import socket
import subprocess
import sys
import tempfile
import threading
import time
import uuid
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

from shared_namespace_lock import NamespaceLockError, namespace_lock

SHARED_LOCK_FORMAT = "shared-lock-v3"
FIELD_SEPARATOR = "\x1f"


@dataclass
class LockOwnership:
    token: str
    heartbeat_state: tuple | None = None


@contextmanager
def lock_mutation_guard(lock_path):
    # A persistent namespace guard serializes *all* changes to a lease. Do not
    # unlink its coordination state: that would let contenders lock different
    # generations of the guard. This also works on Lustre localflock mounts.
    path = Path(lock_path)
    previous_mask = signal.pthread_sigmask(signal.SIG_BLOCK, {signal.SIGHUP, signal.SIGINT, signal.SIGTERM})
    try:
        with namespace_lock(path.with_name("." + path.name + ".guard"), exclusive=True):
            yield
    finally:
        signal.pthread_sigmask(signal.SIG_SETMASK, previous_mask)


def lock_pid_is_alive(pid):
    try:
        os.kill(int(pid), 0)
        return True
    except (OSError, TypeError, ValueError):
        return False


def lock_hostname():
    return socket.gethostname()


def lock_boot_id():
    boot_id_path = Path("/proc/sys/kernel/random/boot_id")
    try:
        if boot_id_path.is_file():
            return boot_id_path.read_text(encoding="utf-8").strip()
    except OSError:
        pass
    try:
        completed = subprocess.run(
            ["sysctl", "-n", "kern.bootsessionuuid"],
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError:
        return ""
    if completed.returncode != 0:
        return ""
    return completed.stdout.strip()


def parse_lock_metadata(raw):
    metadata = {
        "format": "",
        "pid": None,
        "hostname": "",
        "boot_id": "",
        "created_at": "",
        "token": "",
    }
    if raw == "":
        return metadata
    try:
        payload = json.loads(raw)
    except json.JSONDecodeError:
        return metadata
    if not isinstance(payload, dict):
        return metadata
    metadata["format"] = str(payload.get("format", "") or "")
    pid = payload.get("pid")
    if isinstance(pid, int):
        metadata["pid"] = pid
    else:
        try:
            metadata["pid"] = int(str(pid).strip())
        except (TypeError, ValueError):
            metadata["pid"] = None
    metadata["hostname"] = str(payload.get("hostname", "") or "")
    metadata["boot_id"] = str(payload.get("boot_id", "") or "")
    metadata["created_at"] = str(payload.get("created_at", "") or "")
    metadata["token"] = str(payload.get("token", "") or "")
    return metadata


def read_lock_snapshot(lock_path):
    """Read metadata and stat from the same open file, including during handoff."""
    try:
        with Path(lock_path).open(encoding="utf-8") as handle:
            raw = handle.read().strip()
            stat_result = os.fstat(handle.fileno())
    except FileNotFoundError:
        return parse_lock_metadata(""), None
    return parse_lock_metadata(raw), stat_result


def read_lock_metadata(lock_path):
    return read_lock_snapshot(lock_path)[0]


def lock_record_fields(lock_path):
    lock_path = Path(lock_path)
    metadata, stat_result = read_lock_snapshot(lock_path)
    if stat_result is None:
        mtime = ""
        device = ""
        inode = ""
    else:
        mtime = str(stat_result.st_mtime)
        device = str(stat_result.st_dev)
        inode = str(stat_result.st_ino)
    return (
        metadata.get("format", "") or "",
        "" if metadata.get("pid", None) is None else str(metadata.get("pid")),
        metadata.get("hostname", "") or "",
        metadata.get("boot_id", "") or "",
        metadata.get("created_at", "") or "",
        mtime,
        device,
        inode,
        metadata.get("token", ""),
    )


def format_lock_owner_summary(lock_path, metadata=None, stat_result=None):
    lock_path = Path(lock_path)
    if metadata is None:
        metadata, stat_result = read_lock_snapshot(lock_path)
    if stat_result is None:
        age_text = "unknown"
    else:
        age = int(time.time() - stat_result.st_mtime)
        if age < 0:
            age = 0
        age_text = "{}s".format(age)
    return "host={}, pid={}, created_at={}, boot_id={}, heartbeat_age={}, lock={}".format(
        metadata.get("hostname", "") or "unknown",
        metadata.get("pid", None) if metadata.get("pid", None) is not None else "unknown",
        metadata.get("created_at", "") or "unknown",
        metadata.get("boot_id", "") or "unknown",
        age_text,
        lock_path,
    )


def stale_lock_reason(lock_path, stale_seconds, metadata=None, stat_result=None):
    lock_path = Path(lock_path)
    if metadata is None:
        metadata, stat_result = read_lock_snapshot(lock_path)
    current_boot_id = lock_boot_id()
    is_same_host_boot = (
        metadata.get("hostname", "") != ""
        and metadata.get("boot_id", "") != ""
        and current_boot_id != ""
        and metadata.get("hostname", "") == lock_hostname()
        and metadata.get("boot_id", "") == current_boot_id
    )
    pid = metadata.get("pid", None)
    if is_same_host_boot and pid is not None and not lock_pid_is_alive(pid):
        return "same_host_same_boot_dead_pid"
    if stat_result is None:
        return ""
    age = int(time.time() - stat_result.st_mtime)
    if age < 0:
        age = 0
    if age < int(stale_seconds):
        return ""
    return "heartbeat_timeout"


def _remove_lock_if_unchanged(lock_path, stat_result):
    """Called only while holding lock_mutation_guard."""
    lock_path = Path(lock_path)
    _, current_stat = read_lock_snapshot(lock_path)
    if current_stat is None or stat_result is None:
        return False
    if current_stat.st_dev != stat_result.st_dev or current_stat.st_ino != stat_result.st_ino:
        return False
    try:
        lock_path.unlink()
    except FileNotFoundError:
        return False
    return True


def write_lock_file(lock_path, pid=None, hostname=None, boot_id=None, token=None):
    lock_path = Path(lock_path)
    payload = {
        "format": SHARED_LOCK_FORMAT,
        "pid": int(os.getpid() if pid is None else pid),
        "hostname": lock_hostname() if hostname is None else str(hostname),
        "boot_id": lock_boot_id() if boot_id is None else str(boot_id),
        "created_at": time.time(),
        "token": token or uuid.uuid4().hex,
    }
    with lock_mutation_guard(lock_path):
        # Publish complete metadata atomically so signal cleanup can always
        # identify ownership as soon as the lock becomes visible.
        fd, temporary = tempfile.mkstemp(prefix="." + lock_path.name + ".", dir=lock_path.parent)
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as handle:
                json.dump(payload, handle, separators=(",", ":"))
                handle.write("\n")
            os.link(temporary, lock_path)
        finally:
            Path(temporary).unlink(missing_ok=True)
    return payload["token"]


def try_create_lock(lock_path, pid=None, hostname=None, boot_id=None, token=None):
    try:
        return write_lock_file(lock_path, pid=pid, hostname=hostname, boot_id=boot_id, token=token)
    except FileExistsError:
        return None


def reclaim_if_stale(lock_path, stale_seconds):
    lock_path = Path(lock_path)
    with lock_mutation_guard(lock_path):
        metadata, stat_result = read_lock_snapshot(lock_path)
        reason = stale_lock_reason(lock_path, stale_seconds, metadata, stat_result)
        if reason == "":
            return None
        owner_summary = format_lock_owner_summary(lock_path, metadata, stat_result)
        if _remove_lock_if_unchanged(lock_path, stat_result):
            return reason, owner_summary
    return None


def update_owned_lock(lock_path, token, *, release=False):
    """Never touch or remove a lease acquired by a replacement owner."""
    if not token:
        return False
    with lock_mutation_guard(lock_path):
        metadata, stat_result = read_lock_snapshot(lock_path)
        if stat_result is None or metadata.get("token") != token:
            return False
        if release:
            return _remove_lock_if_unchanged(lock_path, stat_result)
        os.utime(lock_path, None)
        return True


def start_lock_heartbeat(lock_path, interval_seconds, token):
    lock_path = Path(lock_path)
    stop_event = threading.Event()

    def _heartbeat():
        while not stop_event.wait(float(interval_seconds)):
            try:
                if not update_owned_lock(lock_path, token):
                    return
            except OSError:
                continue

    thread = threading.Thread(
        target=_heartbeat,
        name="shared-lock-heartbeat-{}".format(lock_path.name),
        daemon=True,
    )
    thread.start()
    return stop_event, thread


def acquire_lock(
    lock_path,
    stale_seconds,
    timeout_seconds,
    poll_seconds,
    context,
    warning_callback=None,
    message_label="[shared-lock]",
    heartbeat_interval_seconds=None,
):
    lock_path = Path(lock_path)
    wait_started = time.monotonic()
    wait_logged = False
    lock_path.parent.mkdir(parents=True, exist_ok=True)

    def warn(message):
        if warning_callback is not None:
            warning_callback(message)

    while True:
        try:
            token = write_lock_file(lock_path)
        except FileExistsError:
            recovered = reclaim_if_stale(lock_path, stale_seconds)
            if recovered is not None:
                reason, owner_summary = recovered
                warn(
                    "{} recovered stale lock for {}: {} ({}; {})".format(
                        message_label, context, lock_path, reason, owner_summary
                    )
                )
                continue
            owner_summary = format_lock_owner_summary(lock_path)
            if not wait_logged:
                warn("{} waiting for shared lock for {} ({})".format(message_label, context, owner_summary))
                wait_logged = True
            if time.monotonic() - wait_started >= float(timeout_seconds):
                raise RuntimeError(
                    "{} timed out waiting for shared lock for {} ({})".format(message_label, context, owner_summary)
                ) from None
            time.sleep(float(poll_seconds))
            continue
        ownership = LockOwnership(token)
        if heartbeat_interval_seconds is not None:
            ownership.heartbeat_state = start_lock_heartbeat(lock_path, heartbeat_interval_seconds, token)
        return ownership


def release_lock(lock_path, ownership):
    lock_path = Path(lock_path)
    if ownership.heartbeat_state is not None:
        stop_event, thread = ownership.heartbeat_state
        stop_event.set()
        thread.join(timeout=1.0)
    return update_owned_lock(lock_path, ownership.token, release=True)


def _cmd_read_metadata(args):
    print(FIELD_SEPARATOR.join(lock_record_fields(args.lock_file)))
    return 0


def _cmd_owner_summary(args):
    print(format_lock_owner_summary(args.lock_file))
    return 0


def _cmd_try_create(args):
    if try_create_lock(args.lock_file, pid=args.pid, hostname=args.hostname, boot_id=args.boot_id, token=args.token):
        return 0
    return 1


def _cmd_release(args):
    update_owned_lock(args.lock_file, args.token, release=True)
    return 0


def _cmd_heartbeat(args):
    return 0 if update_owned_lock(args.lock_file, args.token) else 1


def _cmd_reclaim_if_stale(args):
    result = reclaim_if_stale(args.lock_file, args.stale_seconds)
    if result is None:
        return 1
    reason, owner_summary = result
    print("{}; {}".format(reason, owner_summary))
    return 0


def build_arg_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    read_parser = subparsers.add_parser("read-metadata")
    read_parser.add_argument("lock_file")
    read_parser.set_defaults(func=_cmd_read_metadata)

    summary_parser = subparsers.add_parser("owner-summary")
    summary_parser.add_argument("lock_file")
    summary_parser.set_defaults(func=_cmd_owner_summary)

    create_parser = subparsers.add_parser("try-create")
    create_parser.add_argument("lock_file")
    create_parser.add_argument("--pid", type=int)
    create_parser.add_argument("--hostname")
    create_parser.add_argument("--boot-id")
    create_parser.add_argument("--token", required=True)
    create_parser.set_defaults(func=_cmd_try_create)

    for operation, callback in (("release", _cmd_release), ("heartbeat", _cmd_heartbeat)):
        owner_parser = subparsers.add_parser(operation)
        owner_parser.add_argument("lock_file")
        owner_parser.add_argument("--token", required=True)
        owner_parser.set_defaults(func=callback)

    token_parser = subparsers.add_parser("new-token")
    token_parser.set_defaults(func=lambda args: print(uuid.uuid4().hex) or 0)

    reclaim_parser = subparsers.add_parser("reclaim-if-stale")
    reclaim_parser.add_argument("lock_file")
    reclaim_parser.add_argument("--stale-seconds", type=int, required=True)
    reclaim_parser.set_defaults(func=_cmd_reclaim_if_stale)

    return parser


def main(argv=None):
    def stop(signum, frame):
        raise SystemExit(128 + signum)

    # Unwind the short namespace critical section on ordinary job signals.
    # SIGKILL/node loss remains fail-closed, as with other namespace guards.
    for signum in (signal.SIGHUP, signal.SIGINT, signal.SIGTERM):
        signal.signal(signum, stop)
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except (NamespaceLockError, OSError) as exc:
        print(f"Shared lock operation failed: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
