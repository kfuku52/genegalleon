from __future__ import annotations

import importlib.util
import stat
import zipfile
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "container" / "scripts" / "extract_notung_jar.py"
SPEC = importlib.util.spec_from_file_location("extract_notung_jar", SCRIPT_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_extract_notung_jar_selects_latest_and_publishes_atomically(tmp_path: Path):
    archive_path = tmp_path / "notung.zip"
    destination = tmp_path / "bin" / "Notung.jar"
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.writestr("Notung-2.9.1.jar", b"old")
        archive.writestr("release/Notung-2.9.10.jar", b"latest")

    result = MODULE.extract_notung_jar(archive_path, destination)

    assert result == destination
    assert destination.read_bytes() == b"latest"
    assert stat.S_IMODE(destination.stat().st_mode) == 0o755
    assert not list(destination.parent.glob(".Notung.jar.partial.*"))


def test_extract_notung_jar_rejects_symlink_candidate(tmp_path: Path):
    archive_path = tmp_path / "notung.zip"
    link = zipfile.ZipInfo("Notung-2.9.1.jar")
    link.create_system = 3
    link.external_attr = (stat.S_IFLNK | 0o777) << 16
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.writestr(link, b"../../payload")

    with pytest.raises(MODULE.NotungArchiveError, match="was not found"):
        MODULE.extract_notung_jar(archive_path, tmp_path / "Notung.jar")


def test_extract_notung_jar_rejects_special_file_candidate(tmp_path: Path):
    archive_path = tmp_path / "notung.zip"
    fifo = zipfile.ZipInfo("Notung-2.9.1.jar")
    fifo.create_system = 3
    fifo.external_attr = (stat.S_IFIFO | 0o600) << 16
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.writestr(fifo, b"")

    with pytest.raises(MODULE.NotungArchiveError, match="was not found"):
        MODULE.extract_notung_jar(archive_path, tmp_path / "Notung.jar")


def test_extract_notung_jar_failure_preserves_existing_destination(tmp_path: Path):
    archive_path = tmp_path / "notung.zip"
    archive_path.write_bytes(b"not a zip")
    destination = tmp_path / "Notung.jar"
    destination.write_bytes(b"existing")

    with pytest.raises(MODULE.NotungArchiveError, match="Failed to extract"):
        MODULE.extract_notung_jar(archive_path, destination)

    assert destination.read_bytes() == b"existing"
    assert not list(tmp_path.glob(".Notung.jar.partial.*"))


def test_container_builds_do_not_extract_notung_with_unzip():
    dockerfile = (REPO_ROOT / "container" / "Dockerfile").read_text(encoding="utf-8")
    apptainer = (REPO_ROOT / "container" / "apptainer_local_build.def.template").read_text(
        encoding="utf-8"
    )

    for text in (dockerfile, apptainer):
        assert "unzip -q" not in text
        assert "extract_notung_jar.py" in text
