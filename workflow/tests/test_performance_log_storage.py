"""Original performance bodies remain readable with bounded live file counts."""

import json
import os
import sys
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"
if str(SUPPORT) not in sys.path:
    sys.path.insert(0, str(SUPPORT))
import metadata_store as metadata  # noqa: E402
import performance_log as log  # noqa: E402


def frame(pid):
    return {
        "schema_version": 1,
        "phase": "fixture",
        "pid": pid,
        "host": "fixture",
        "job_id": "12345",
        "array_task_id": "1",
        "seconds": 1.25,
        "status": "ok",
    }


def emit(root, pid):
    raw = (json.dumps(frame(pid), sort_keys=True) + "\n").encode()
    log.append(root, pid, raw)
    return raw


def test_live_stripes_collection_generations_and_restore_preserve_all_records(tmp_path):
    root = tmp_path / "versions/performance"
    for pid in range(1, 257):
        emit(root, pid)
    assert len(list(root.glob("*.jsonl"))) == 64
    before = list(log.records(root))
    source = {p.name: (p.read_bytes(), p.stat().st_mode & 0o777, p.stat().st_mtime_ns) for p in root.glob("*.jsonl")}
    first = metadata.pack(root, writers_stopped=True, delete_originals=True, operation_id="first")
    assert first["removed_files"] == 64 and not list(root.glob("*.jsonl"))
    assert list(log.records(root)) == before
    for pid in range(257, 513):
        emit(root, pid)
    assert len(list(root.glob("*.jsonl"))) == 64
    assert all(not p.name.startswith("trace-initial-") for p in root.glob("*.jsonl"))
    combined = list(log.records(root))
    assert sorted(combined, key=lambda r: r["pid"]) == [frame(i) for i in range(1, 513)]
    metadata.pack(root, writers_stopped=True, delete_originals=True, operation_id="second")
    assert list(log.records(root)) == combined
    metadata.restore(root, writers_stopped=True)
    for name, (raw, mode, mtime) in source.items():
        assert (root / name).read_bytes() == raw
        assert (root / name).stat().st_mode & 0o777 == mode
        assert (root / name).stat().st_mtime_ns == mtime
    assert list(log.records(root)) == combined


def test_legacy_pid_jsonl_and_custom_directories_keep_their_layout_and_bytes(tmp_path):
    root = tmp_path / "versions/performance"
    root.mkdir(parents=True)
    old = root / "123.jsonl"
    raw = b'{"pid":123,"phase":"old","status":"exit"}\n'
    old.write_bytes(raw)
    (root / "notes.tsv").write_text("keep\tstatistics\n")
    metadata.pack(root, writers_stopped=True, delete_originals=True)
    assert metadata.archived_bytes(old) == raw and (root / "notes.tsv").read_text() == "keep\tstatistics\n"
    assert list(log.records(root)) == [json.loads(raw)]
    custom = tmp_path / "custom"
    assert emit(custom, 123) == (custom / "123.jsonl").read_bytes()
    with pytest.raises(ValueError, match="registered"):
        metadata.pack(custom, writers_stopped=True)


def test_readers_reject_changed_or_invalid_bodies_and_stripe_redirection(tmp_path):
    root = tmp_path / "versions/performance"
    emit(root, 1)
    stripe = root / "trace-initial-01.jsonl"
    outside = tmp_path / "keep"
    outside.write_bytes(b"keep")
    stripe.unlink()
    stripe.symlink_to(outside)
    with pytest.raises(OSError):
        emit(root, 1)
    assert outside.read_bytes() == b"keep"
    stripe.unlink()
    os.link(outside, stripe)
    with pytest.raises(ValueError, match="Unsafe"):
        emit(root, 1)
    assert outside.read_bytes() == b"keep"


def test_unrelated_performance_directory_is_not_a_registered_metadata_root(tmp_path):
    root = tmp_path / "performance"
    root.mkdir()
    (root / "1.jsonl").write_text("{}\n")
    with pytest.raises(ValueError, match="registered"):
        metadata.pack(root, writers_stopped=True)


def test_same_stripe_concurrent_writers_preserve_complete_json_records(tmp_path):
    from concurrent.futures import ThreadPoolExecutor

    root = tmp_path / "versions/performance"
    with ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(lambda i: emit(root, 64 * i + 1), range(1, 81)))
    assert len(list(root.glob("*.jsonl"))) == 1
    assert sorted(record["pid"] for record in log.records(root)) == [64 * i + 1 for i in range(1, 81)]


def test_package_import_and_enabled_telemetry_preserve_the_public_import_boundary(tmp_path):
    import os
    import subprocess

    root = tmp_path / "versions" / "performance"
    code = (
        "from workflow.support.performance_metrics import write_row;"
        "from workflow.support.performance_log import records;"
        "write_row({'phase':'package-import','elapsed':1});"
        "import os,json;print(json.dumps(list(records(os.environ['GG_PERFORMANCE_DIR']))))"
    )
    environment = {key: value for key, value in os.environ.items() if key != "PYTHONPATH"}
    environment["GG_PERFORMANCE_DIR"] = str(root)
    result = subprocess.run(
        [sys.executable, "-c", code], cwd=SUPPORT.parents[1], env=environment, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr
    rows = json.loads(result.stdout)
    assert len(rows) == 1 and rows[0]["phase"] == "package-import" and rows[0]["elapsed"] == 1
