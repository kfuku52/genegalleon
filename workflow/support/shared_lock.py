#!/usr/bin/env python3
"""Shared lock metadata helpers used by Python and shell workflow code."""

import argparse
import json
import os
import socket
import subprocess
import threading
import time
from pathlib import Path

SHARED_LOCK_FORMAT = "shared-lock-v2"
FIELD_SEPARATOR = "\x1f"


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


def read_lock_metadata(lock_path):
    lock_path = Path(lock_path)
    metadata = {
        "format": "",
        "pid": None,
        "hostname": "",
        "boot_id": "",
        "created_at": "",
    }
    try:
        raw = lock_path.read_text(encoding="utf-8").strip()
    except OSError:
        return metadata
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
    return metadata


def read_lock_stat(lock_path):
    try:
        return Path(lock_path).stat()
    except OSError:
        return None


def lock_record_fields(lock_path):
    lock_path = Path(lock_path)
    metadata = read_lock_metadata(lock_path)
    stat_result = read_lock_stat(lock_path)
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
    )


def format_lock_owner_summary(lock_path, metadata=None, stat_result=None):
    lock_path = Path(lock_path)
    if metadata is None:
        metadata = read_lock_metadata(lock_path)
    if stat_result is None:
        stat_result = read_lock_stat(lock_path)
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
    if not lock_path.exists():
        return ""
    if metadata is None:
        metadata = read_lock_metadata(lock_path)
    if stat_result is None:
        stat_result = read_lock_stat(lock_path)
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


def remove_lock_if_unchanged(lock_path, stat_result):
    lock_path = Path(lock_path)
    current_stat = read_lock_stat(lock_path)
    if current_stat is None or stat_result is None:
        return False
    if current_stat.st_dev != stat_result.st_dev or current_stat.st_ino != stat_result.st_ino:
        return False
    try:
        lock_path.unlink()
    except FileNotFoundError:
        return False
    except OSError:
        return False
    return True


def write_lock_file(lock_path, pid=None, hostname=None, boot_id=None):
    lock_path = Path(lock_path)
    payload = {
        "format": SHARED_LOCK_FORMAT,
        "pid": int(os.getpid() if pid is None else pid),
        "hostname": lock_hostname() if hostname is None else str(hostname),
        "boot_id": lock_boot_id() if boot_id is None else str(boot_id),
        "created_at": time.time(),
    }
    fd = os.open(str(lock_path), os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, separators=(",", ":"))
            handle.write("\n")
    except Exception:
        try:
            lock_path.unlink()
        except OSError:
            pass
        raise


def try_create_lock(lock_path, pid=None, hostname=None, boot_id=None):
    try:
        write_lock_file(lock_path, pid=pid, hostname=hostname, boot_id=boot_id)
    except FileExistsError:
        return False
    return True


def reclaim_if_stale(lock_path, stale_seconds):
    lock_path = Path(lock_path)
    if not lock_path.exists():
        return None
    metadata = read_lock_metadata(lock_path)
    stat_result = read_lock_stat(lock_path)
    reason = stale_lock_reason(lock_path, stale_seconds, metadata, stat_result)
    if reason == "":
        return None
    owner_summary = format_lock_owner_summary(lock_path, metadata, stat_result)
    if remove_lock_if_unchanged(lock_path, stat_result):
        return reason, owner_summary
    return None


def start_lock_heartbeat(lock_path, interval_seconds):
    lock_path = Path(lock_path)
    stop_event = threading.Event()

    def _heartbeat():
        while not stop_event.wait(float(interval_seconds)):
            try:
                os.utime(lock_path, None)
            except FileNotFoundError:
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
            write_lock_file(lock_path)
        except FileExistsError:
            metadata = read_lock_metadata(lock_path)
            stat_result = read_lock_stat(lock_path)
            reason = stale_lock_reason(lock_path, stale_seconds, metadata, stat_result)
            if reason != "":
                owner_summary = format_lock_owner_summary(lock_path, metadata, stat_result)
                if remove_lock_if_unchanged(lock_path, stat_result):
                    warn(
                        "{} recovered stale lock for {}: {} ({}; {})".format(
                            message_label, context, lock_path, reason, owner_summary
                        )
                    )
                    continue
            owner_summary = format_lock_owner_summary(lock_path, metadata, stat_result)
            if not wait_logged:
                warn("{} waiting for shared lock for {} ({})".format(message_label, context, owner_summary))
                wait_logged = True
            if time.monotonic() - wait_started >= float(timeout_seconds):
                raise RuntimeError(
                    "{} timed out waiting for shared lock for {} ({})".format(
                        message_label, context, owner_summary
                    )
                )
            time.sleep(float(poll_seconds))
            continue
        if heartbeat_interval_seconds is None:
            return None
        return start_lock_heartbeat(lock_path, heartbeat_interval_seconds)


def release_lock(lock_path, heartbeat_state=None):
    lock_path = Path(lock_path)
    if heartbeat_state is not None:
        stop_event, thread = heartbeat_state
        stop_event.set()
        thread.join(timeout=1.0)
    try:
        lock_path.unlink()
    except FileNotFoundError:
        return
    except OSError:
        return


def _cmd_read_metadata(args):
    print(FIELD_SEPARATOR.join(lock_record_fields(args.lock_file)))
    return 0


def _cmd_owner_summary(args):
    print(format_lock_owner_summary(args.lock_file))
    return 0


def _cmd_try_create(args):
    if try_create_lock(args.lock_file, pid=args.pid, hostname=args.hostname, boot_id=args.boot_id):
        return 0
    return 1


def _cmd_remove_if_unchanged(args):
    class _Stat:
        st_dev = int(args.device)
        st_ino = int(args.inode)

    if remove_lock_if_unchanged(args.lock_file, _Stat()):
        return 0
    return 1


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
    create_parser.set_defaults(func=_cmd_try_create)

    remove_parser = subparsers.add_parser("remove-if-unchanged")
    remove_parser.add_argument("lock_file")
    remove_parser.add_argument("device")
    remove_parser.add_argument("inode")
    remove_parser.set_defaults(func=_cmd_remove_if_unchanged)

    reclaim_parser = subparsers.add_parser("reclaim-if-stale")
    reclaim_parser.add_argument("lock_file")
    reclaim_parser.add_argument("--stale-seconds", type=int, required=True)
    reclaim_parser.set_defaults(func=_cmd_reclaim_if_stale)

    return parser


def main(argv=None):
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
