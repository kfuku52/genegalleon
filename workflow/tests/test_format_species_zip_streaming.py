from __future__ import annotations

import gzip
import io
import sys
import zipfile
from pathlib import Path

import pytest

SUPPORT_DIR = Path(__file__).resolve().parents[1] / "support"
if str(SUPPORT_DIR) not in sys.path:
    sys.path.insert(0, str(SUPPORT_DIR))

from format_species_download import locking, targets  # noqa: E402
from format_species_download.cache_validation import (  # noqa: E402
    GzipValidationCache,
    gzip_validation_key_for_target,
)


def _zip_bytes(member_name: str, payload: bytes) -> bytes:
    output = io.BytesIO()
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr(member_name, payload)
    return output.getvalue()


def _reject_whole_member_reads(*args, **kwargs):
    raise AssertionError("ZIP members must be streamed with ZipFile.open, not read whole")


def _archive_response(payload):
    response = io.BytesIO(payload)
    response.status = 200
    response.headers = {"Content-Length": str(len(payload)), "Content-Type": "application/zip"}
    return response


def test_direct_archive_member_download_streams_zip_payload(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    url = "http://127.0.0.1/archive.zip"
    member_name = "release/genome.fasta"
    payload = b">chr1\n" + (b"ACGT" * 4096) + b"\n"
    destination = tmp_path / "download" / "genome.fasta"
    response_payload = _zip_bytes(member_name, payload)
    monkeypatch.setattr(locking, "urlopen", lambda *args, **kwargs: _archive_response(response_payload))
    monkeypatch.setattr(locking, "acquire_download_lock", lambda *args, **kwargs: None)
    monkeypatch.setattr(locking, "release_download_lock", lambda *args, **kwargs: None)
    monkeypatch.setattr(locking.zipfile.ZipFile, "read", _reject_whole_member_reads)

    changed = locking.download_url_to_file(
        url,
        destination,
        {},
        10,
        False,
        False,
        60,
        [],
        "streaming-test",
        archive_member=member_name,
    )

    assert changed is True
    assert destination.read_bytes() == payload


def test_ncbi_datasets_member_download_streams_zip_payload(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    member_name = "ncbi_dataset/data/GCF_TEST/cds_from_genomic.fna"
    payload = b">gene1\n" + (b"ATGC" * 4096) + b"\n"
    response_payload = _zip_bytes(member_name, payload)
    destination = tmp_path / "cds.fna.gz"
    monkeypatch.setattr(targets, "acquire_download_lock", lambda *args, **kwargs: None)
    monkeypatch.setattr(targets, "release_download_lock", lambda *args, **kwargs: None)
    monkeypatch.setattr(locking, "urlopen", lambda *args, **kwargs: _archive_response(response_payload))
    monkeypatch.setattr(targets.zipfile.ZipFile, "read", _reject_whole_member_reads)
    validation_cache = GzipValidationCache(tmp_path / "validation-cache")
    validation_key = gzip_validation_key_for_target(destination, tmp_path, "ncbi-datasets://GCF_TEST/CDS")

    changed = targets.download_ncbi_datasets_file_from_id(
        "GCF_TEST",
        "CDS",
        destination,
        {},
        10,
        False,
        False,
        60,
        [],
        "streaming-test",
        validation_cache,
        validation_key,
        str(destination.relative_to(tmp_path)),
        "ncbi-datasets://GCF_TEST/CDS",
    )

    assert changed is True
    with gzip.open(destination, "rb") as handle:
        assert handle.read() == payload
    assert validation_cache.diagnostics()["validation_cache_records"] == 1
