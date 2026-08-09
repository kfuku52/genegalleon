from __future__ import annotations

import gzip
import hashlib
import io
import sys
import zipfile
from pathlib import Path

import pytest

SUPPORT_DIR = Path(__file__).resolve().parents[1] / "support"
if str(SUPPORT_DIR) not in sys.path:
    sys.path.insert(0, str(SUPPORT_DIR))

from format_species_download import locking, targets  # noqa: E402


def _zip_bytes(member_name: str, payload: bytes) -> bytes:
    output = io.BytesIO()
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr(member_name, payload)
    return output.getvalue()


def _reject_whole_member_reads(*args, **kwargs):
    raise AssertionError("ZIP members must be streamed with ZipFile.open, not read whole")


def test_direct_archive_member_download_streams_zip_payload(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    url = "http://127.0.0.1/archive.zip"
    member_name = "release/genome.fasta"
    payload = b">chr1\n" + (b"ACGT" * 4096) + b"\n"
    destination = tmp_path / "download" / "genome.fasta"
    archive_name = "archive.zip"
    archive_hash = hashlib.sha256(url.encode("utf-8")).hexdigest()[:16]
    archive_cache = destination.parent / ".archive_cache"
    archive_cache.mkdir(parents=True)
    cached_archive = archive_cache / f"{archive_hash}__{archive_name}"
    cached_archive.write_bytes(_zip_bytes(member_name, payload))
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
    monkeypatch.setattr(targets, "urlopen", lambda *args, **kwargs: io.BytesIO(response_payload))
    monkeypatch.setattr(targets.zipfile.ZipFile, "read", _reject_whole_member_reads)

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
    )

    assert changed is True
    with gzip.open(destination, "rb") as handle:
        assert handle.read() == payload
