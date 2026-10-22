"""Real loopback transport tests for declared opaque archive parts and models."""
import gzip
import hashlib
import json
import socket
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))

from annevo_runtime_download import fetch_verified
from format_species_download.locking import download_url_to_file, prepare_partial_download
from test_download_reliability import server


def identity(body):
    return {"expected_sha256": hashlib.sha256(body).hexdigest(), "expected_size": len(body)}


def download(url, destination, **options):
    return download_url_to_file(url, destination, {}, 3, False, False, 60, [], "byte-contract fixture", **options)


@pytest.fixture(autouse=True)
def bounded_local_transport(monkeypatch):
    monkeypatch.setenv("GG_DOWNLOAD_RETRY_BASE_SECONDS", "0")
    monkeypatch.setenv("GG_DOWNLOAD_ATTEMPTS", "2")
    monkeypatch.setenv("GG_TEST_ALLOW_ONLY_LOOPBACK_HTTP", "1")
    monkeypatch.delenv("GG_INPUT_DOWNLOAD_LIMIT_DIR", raising=False)
    monkeypatch.delenv("GG_DOWNLOAD_SHARED_CACHE_DIR", raising=False)


def response(body, calls, content_type="application/octet-stream"):
    def reply(handler):
        calls.append(handler.headers.get("Range"))
        handler.send_response(200)
        handler.send_header("Content-Type", content_type)
        handler.send_header("Content-Length", str(len(body)))
        handler.send_header("ETag", '"same"')
        handler.end_headers()
        handler.wfile.write(body)
    return reply


@pytest.mark.parametrize("body", [gzip.compress(b"opaque runtime archive body", mtime=0)[:-6], b"plain bytes", b"PK\x03\x04ZIP fragment"])
def test_exact_declared_bytes_are_fetched_once_reused_and_never_normalized(tmp_path, body):
    calls = []
    target = tmp_path / "runtime.tar.gz"
    with server(response(body, calls)) as url:
        assert download(url, target, **identity(body)) is True
        assert download(url, target, **identity(body)) is False
    assert calls == [None]
    assert target.read_bytes() == body
    assert not list(tmp_path.glob("*.corrupt.*"))


@pytest.mark.parametrize("bad", ["SHA", "size"])
def test_wrong_declared_hash_or_size_never_publishes_a_destination(tmp_path, bad):
    body = gzip.compress(b"runtime archive body")[:-6]
    expected = identity(body)
    expected["expected_sha256" if bad == "SHA" else "expected_size"] = "a" * 64 if bad == "SHA" else len(body) + 1
    target = tmp_path / "part0000"
    with server(response(body, [])) as url, pytest.raises((OSError, ValueError), match="expected|integrity"):
        download(url, target, **expected)
    assert not target.exists()


@pytest.mark.parametrize("expected", [{"expected_sha256": "a" * 64}, {"expected_size": 12},
                                     {"expected_sha256": "bad", "expected_size": 12},
                                     {"expected_sha256": "a" * 64, "expected_size": False},
                                     {"expected_sha256": "a" * 64, "expected_size": 0}])
def test_both_strong_identity_fields_are_mandatory_before_any_network_io(tmp_path, expected):
    calls = []
    with server(response(b"body", calls)) as url, pytest.raises(ValueError, match="both"):
        download(url, tmp_path / "destination", **expected)
    assert calls == []
    assert not list(tmp_path.iterdir())


def test_byte_contract_retains_safe_range_resume_and_exact_final_hash(tmp_path):
    body = gzip.compress(b"runtime split archive body")[:-6]
    split = len(body) // 2
    calls = []
    def reply(handler):
        ranged = handler.headers.get("Range")
        calls.append(ranged)
        handler.send_response(206 if ranged else 200)
        handler.send_header("ETag", '"same"')
        payload = body[split:] if ranged else body
        handler.send_header("Content-Length", str(len(payload)))
        if ranged:
            assert ranged == "bytes={}-".format(split)
            assert handler.headers.get("If-Range") == '"same"'
            handler.send_header("Content-Range", "bytes {}-{}/{}".format(split, len(body) - 1, len(body)))
        handler.end_headers()
        handler.wfile.write(payload if ranged else body[:split])
        if not ranged:
            handler.wfile.flush()
            handler.connection.shutdown(socket.SHUT_RDWR)
            handler.connection.close()
    with server(reply) as url:
        download(url, tmp_path / "part0", **identity(body))
    assert calls == [None, "bytes={}-".format(split)]
    assert (tmp_path / "part0").read_bytes() == body


def test_changed_byte_contract_discards_partial_identity_from_another_declared_object(tmp_path):
    old, body = b"old object", b"new exact object"
    calls = []
    target = tmp_path / "part0"
    with server(response(body, calls)) as url:
        partial, _ = prepare_partial_download(target, url, False, **identity(old))
        partial.write_bytes(old[:4])
        Path(str(partial) + ".identity.json").write_text(json.dumps({"etag": '"same"', "total": len(old)}))
        download(url, target, **identity(body))
    assert calls == [None]
    assert target.read_bytes() == body


def test_shared_cache_raw_fragments_reuse_copy_and_do_not_alias_structural_or_other_hash_modes(tmp_path, monkeypatch):
    cache = tmp_path / "shared"
    monkeypatch.setenv("GG_DOWNLOAD_SHARED_CACHE_DIR", str(cache))
    body = gzip.compress(b"runtime archive body")[:-6]
    calls = []
    with server(response(body, calls)) as url:
        first, second = tmp_path / "planA/part0", tmp_path / "planB/part0"
        download(url, first, **identity(body))
        download(url, second, **identity(body))
        assert len(calls) == 1 and first.stat().st_ino != second.stat().st_ino
        first.write_bytes(b"downstream edit")
        download(url, first, **identity(body))
        assert len(calls) == 1 and first.read_bytes() == second.read_bytes() == body
        with pytest.raises(OSError):
            download(url, tmp_path / "structural/part0")
        assert len(calls) == 3
        with pytest.raises(OSError):
            download(url, tmp_path / "other_hash/part0", expected_sha256="a" * 64, expected_size=len(body))
        assert len(calls) == 5
    assert len(list(cache.iterdir())) == 3
    receipts = list(cache.rglob("receipt*.json"))
    assert len(receipts) == 1
    assert json.loads(receipts[0].read_text())["byte_integrity"] == {"sha256": identity(body)["expected_sha256"], "size": len(body)}


def test_default_compressed_and_html_scientific_contracts_stay_strict(tmp_path):
    body = gzip.compress(b"<html>provider error</html>")
    with server(response(body, [])) as url, pytest.raises(OSError):
        download(url, tmp_path / "science.gz")
    html = b"<html>declared opaque artifact</html>"
    calls = []
    with server(response(html, calls, "text/html")) as url:
        with pytest.raises(ValueError, match="HTML"):
            download(url, tmp_path / "science")
        download(url, tmp_path / "opaque", **identity(html))
    assert (tmp_path / "opaque").read_bytes() == html
    assert not (tmp_path / "science").exists()


def test_annevo_http_assets_use_the_exact_byte_transport_for_split_gzip(tmp_path):
    body = gzip.compress(b"ANNEVO runtime archive content")[:-8]
    calls = []
    with server(response(body, calls)) as url:
        assert fetch_verified(url, tmp_path / "bundle.part0000", identity(body)["expected_sha256"], len(body)) is True
        assert fetch_verified(url, tmp_path / "bundle.part0000", identity(body)["expected_sha256"], len(body)) is False
    assert calls == [None]
    assert (tmp_path / "bundle.part0000").read_bytes() == body
