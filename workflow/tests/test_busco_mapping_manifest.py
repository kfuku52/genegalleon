"""BUSCO's standard manifest works without an HTML directory listing."""

import hashlib
import io
import os
import re
import shlex
import subprocess
import tarfile
import urllib.error
import urllib.request
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"
BASE = "https://busco-data.ezlab.org/v5/data/"
DOMAINS = ("archaea", "bacteria", "eukaryota")


def inline(function):
    file = "03_species_helpers.sh" if function.startswith("gg_fetch") else "04_busco_runtime.sh"
    text = (SUPPORT / "gg_util" / file).read_text().split(function + "() {", 1)[1]
    return text.split("<<'PY'\n", 1)[1].split("\nPY\n", 1)[0]


def manifest_fixture(version="12"):
    rows = []
    archives = {}
    for domain in DOMAINS:
        name = f"mapping_taxids-busco_dataset_name.{domain}_odb{version}"
        member = f"{name}.2025-01-15.txt"
        body = f"123\t{domain}_odb{version}\n".encode()
        stream = io.BytesIO()
        with tarfile.open(fileobj=stream, mode="w:gz") as tar:
            info = tarfile.TarInfo(member)
            info.size = len(body)
            tar.addfile(info, io.BytesIO(body))
        archive = stream.getvalue()
        archives[BASE + "placement_files/" + member + ".tar.gz"] = archive
        rows.append(f"{name}.txt\t2025-01-15\t{hashlib.md5(archive).hexdigest()}\t{domain}\tplacement_files")
    return "\n".join(rows) + "\n", archives


def responses(monkeypatch, manifest, archives=None):
    calls = []

    def urlopen(url, timeout):
        assert timeout == 120
        calls.append(url)
        if url == BASE + "file_versions.tsv":
            return io.BytesIO(manifest.encode())
        if url == BASE + "placement_files/":
            return io.BytesIO(b"")  # S3 directory marker: HTTP 200 with zero bytes.
        return io.BytesIO((archives or {})[url])

    monkeypatch.setattr(urllib.request, "urlopen", urlopen)
    return calls


def download(monkeypatch, tmp_path, version="12"):
    monkeypatch.setenv("GG_BUSCO_MAPPING_DIR", str(tmp_path / "mappings"))
    monkeypatch.setenv("GG_BUSCO_MAPPING_STAMP", str(tmp_path / "ready.tsv"))
    monkeypatch.setenv("GG_BUSCO_MAPPING_ODB_VERSION", version)
    exec(compile(inline("_download_busco_dataset_mapping_files_locked"), "<mapping-download>", "exec"), {})


def test_version_uses_standard_manifest_when_directory_response_is_empty(monkeypatch, capsys):
    manifest, _ = manifest_fixture()
    calls = responses(monkeypatch, manifest)
    exec(compile(inline("gg_fetch_latest_busco_mapping_odb_version"), "<mapping-version>", "exec"), {})
    assert capsys.readouterr().out == "12\n"
    assert calls == [BASE + "file_versions.tsv"]


def test_version_keeps_numeric_common_integer_odb_selection(monkeypatch, capsys):
    old, _ = manifest_fixture("9")
    latest, _ = manifest_fixture("12")
    fractional, _ = manifest_fixture("12.2")
    partial, _ = manifest_fixture("13")
    responses(monkeypatch, old + latest + fractional + partial.splitlines()[0] + "\n")
    exec(inline("gg_fetch_latest_busco_mapping_odb_version"), {})
    assert capsys.readouterr().out == "12\n"


def test_missing_domain_still_fails(monkeypatch):
    manifest, _ = manifest_fixture()
    responses(monkeypatch, "\n".join(manifest.splitlines()[:2]))
    with pytest.raises(SystemExit, match="No common BUSCO ODB"):
        exec(inline("gg_fetch_latest_busco_mapping_odb_version"), {})


@pytest.mark.parametrize("bad", ["not-a-date", "2025-99-99"])
def test_bad_manifest_date_is_rejected(monkeypatch, bad):
    manifest, _ = manifest_fixture()
    responses(monkeypatch, manifest.replace("2025-01-15", bad, 1))
    with pytest.raises(SystemExit, match="Invalid BUSCO placement manifest"):
        exec(inline("gg_fetch_latest_busco_mapping_odb_version"), {})


def test_duplicate_conflicting_record_is_rejected(monkeypatch):
    manifest, _ = manifest_fixture()
    responses(monkeypatch, manifest + manifest.splitlines()[0].replace("2025-01-15", "2025-01-16") + "\n")
    with pytest.raises(SystemExit, match="Conflicting BUSCO placement manifest"):
        exec(inline("gg_fetch_latest_busco_mapping_odb_version"), {})


def test_download_uses_manifest_dates_and_preserves_stamp_and_mapping_names(monkeypatch, tmp_path):
    manifest, archives = manifest_fixture()
    calls = responses(monkeypatch, manifest, archives)
    download(monkeypatch, tmp_path)
    assert calls == [BASE + "file_versions.tsv", *archives]
    stamp = (tmp_path / "ready.tsv").read_text().splitlines()
    assert stamp[:2] == ["odb_version\t12", "domain\tarchive\tmapping_file"]
    assert len(stamp) == 5
    for domain in DOMAINS:
        name = f"mapping_taxids-busco_dataset_name.{domain}_odb12.2025-01-15.txt"
        assert (tmp_path / "mappings" / name).read_text() == f"123\t{domain}_odb12\n"


def test_archive_checksum_mismatch_never_creates_ready_stamp(monkeypatch, tmp_path):
    manifest, archives = manifest_fixture()
    fields = manifest.splitlines()[0].split("\t")
    responses(monkeypatch, manifest.replace(fields[2], "0" * 32, 1), archives)
    with pytest.raises(SystemExit, match="checksum mismatch"):
        download(monkeypatch, tmp_path)
    assert not (tmp_path / "ready.tsv").exists()


def test_download_missing_domain_never_stamps_partial_cache(monkeypatch, tmp_path):
    manifest, archives = manifest_fixture()
    responses(monkeypatch, "\n".join(manifest.splitlines()[:2]), archives)
    with pytest.raises(SystemExit, match="domain/version not found"):
        download(monkeypatch, tmp_path)
    assert not (tmp_path / "ready.tsv").exists()


def test_http_failure_is_reported_without_changing_odb(monkeypatch):
    def fail(*args, **kwargs):
        raise urllib.error.HTTPError(BASE + "file_versions.tsv", 503, "unavailable", {}, None)

    monkeypatch.setattr(urllib.request, "urlopen", fail)
    with pytest.raises(urllib.error.HTTPError):
        exec(inline("gg_fetch_latest_busco_mapping_odb_version"), {})


def test_failed_remote_and_empty_cache_preserve_original_error(tmp_path):
    script = (
        "set -euo pipefail\nsource " + shlex.quote(str(SUPPORT / "gg_util.sh")) + "\n"
        "gg_fetch_latest_busco_mapping_odb_version() { echo 'provider manifest unavailable' >&2; return 1; }\n"
        "ensure_busco_dataset_mapping_files " + shlex.quote(str(tmp_path))
    )
    result = subprocess.run(["bash", "-c", script], capture_output=True, text=True,
                            env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
    assert result.returncode != 0
    assert "provider manifest unavailable" in result.stderr
    assert "Failed to determine a BUSCO placement mapping ODB version" in result.stderr


def test_named_lineage_still_does_not_fetch_placement_manifest():
    text = (SUPPORT / "gg_util/04_busco_runtime.sh").read_text()
    function = text.split("gg_resolve_busco_lineage() {", 1)[1]
    assert function.index('"${requested_lc}" != "auto"') < function.index("ensure_busco_dataset_mapping_files")
    assert re.search(r'printf.*normalized_requested.*\n\s*return 0', function)
