"""Fault injection for resumable downloads; all HTTP stays on loopback."""
import gzip
import sys
import threading
from contextlib import contextmanager
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.error import HTTPError

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'support'))
from format_species_download.locking import download_url_to_file
from format_species_download.targets import pick_ncbi_datasets_member_name
from format_species_network import download_retry_delay


@contextmanager
def server(reply):
    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            reply(self)
        def log_message(self, *args):
            pass
    http = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
    thread = threading.Thread(target=http.serve_forever, daemon=True)
    thread.start()
    try:
        yield f'http://127.0.0.1:{http.server_port}/data'
    finally:
        http.shutdown()
        http.server_close()
        thread.join()


def download(url, target):
    return download_url_to_file(url, target, {}, 3, False, False, 60, [], 'test')


@pytest.fixture(autouse=True)
def no_delays(monkeypatch):
    monkeypatch.setenv('GG_DOWNLOAD_RETRY_BASE_SECONDS', '0')
    monkeypatch.setenv('GG_DOWNLOAD_ATTEMPTS', '3')
    monkeypatch.setenv('GG_TEST_ALLOW_ONLY_LOOPBACK_HTTP', '1')
    monkeypatch.delenv('GG_INPUT_DOWNLOAD_LIMIT_DIR', raising=False)


def test_chunked_ranges_produce_exact_file(tmp_path, monkeypatch):
    payload = b'>sequence\nACGTACGTACGT\n'
    requests = []
    monkeypatch.setenv('GG_DOWNLOAD_RANGE_CHUNK_BYTES', '7')
    def reply(h):
        requests.append((h.headers.get('Range'), h.headers.get('If-Range')))
        start, end = map(int, h.headers['Range'][6:].split('-'))
        end = min(end, len(payload)-1)
        h.send_response(206)
        h.send_header('ETag', '"v1"')
        h.send_header('Content-Range', f'bytes {start}-{end}/{len(payload)}')
        h.send_header('Content-Length', str(end-start+1))
        h.end_headers()
        h.wfile.write(payload[start:end+1])
    with server(reply) as url:
        download(url, tmp_path/'sequence')
    assert (tmp_path/'sequence').read_bytes() == payload
    assert requests[1] == ('bytes=7-13', '"v1"')
    assert not list(tmp_path.glob('*.part*'))


def test_changed_etag_restarts_without_mixing(tmp_path, monkeypatch):
    monkeypatch.setenv('GG_DOWNLOAD_RANGE_CHUNK_BYTES', '4')
    seen = []
    def reply(h):
        seen.append(h.headers.get('If-Range'))
        if len(seen) == 1:
            h.send_response(206)
            h.send_header('ETag', '"old"')
            h.send_header('Content-Range', 'bytes 0-3/8')
            body = b'AAAA'
        else:
            h.send_response(200) # If-Range mismatch
            h.send_header('ETag', '"new"')
            body = b'BBBBBBBB'
        h.send_header('Content-Length', str(len(body)))
        h.end_headers()
        h.wfile.write(body)
    with server(reply) as url:
        download(url, tmp_path/'data')
    assert (tmp_path/'data').read_bytes() == b'BBBBBBBB'
    assert seen == [None, '"old"']


def test_midstream_disconnect_resumes(tmp_path):
    payload = b'A' * (1024 * 1024 + 17)
    seen = []
    def reply(h):
        seen.append(h.headers.get('Range'))
        start = int(h.headers.get('Range', 'bytes=0-')[6:].split('-')[0])
        h.send_response(206 if start else 200)
        h.send_header('ETag', '"same"')
        if start:
            h.send_header('Content-Range', f'bytes {start}-{len(payload)-1}/{len(payload)}')
        h.send_header('Content-Length', str(len(payload)-start))
        h.end_headers()
        h.wfile.write(payload[start:] if start else payload[:1024*1024])
        h.close_connection = True
    with server(reply) as url:
        download(url, tmp_path/'data')
    assert (tmp_path/'data').read_bytes() == payload
    assert seen == [None, 'bytes=1048576-']


@pytest.mark.parametrize('kind', ['html_header', 'html_body', 'bad_range', 'corrupt_gzip'])
def test_invalid_response_never_published(tmp_path, monkeypatch, kind):
    if kind == 'bad_range':
        monkeypatch.setenv('GG_DOWNLOAD_RANGE_CHUNK_BYTES', '4')
    def reply(h):
        h.send_response(206 if kind == 'bad_range' else 200)
        if kind == 'bad_range':
            h.send_header('Content-Range', 'bytes 9-12/13')
        if kind == 'html_header':
            h.send_header('Content-Type', 'text/html')
        body = gzip.compress(b'data')[:-5] if kind == 'corrupt_gzip' else b'<html>failure</html>'
        h.send_header('Content-Length', str(len(body)))
        h.end_headers()
        h.wfile.write(body)
    with server(reply) as url, pytest.raises((ValueError, OSError)):
        download(url, tmp_path/'12345')
    assert not (tmp_path/'12345').exists()


def test_retry_budget_is_not_multiplied(tmp_path):
    calls = []
    def reply(h):
        calls.append(1)
        h.send_response(503)
        h.send_header('Retry-After', '0')
        h.end_headers()
    with server(reply) as url, pytest.raises(HTTPError):
        download(url, tmp_path/'data')
    assert len(calls) == 3


def test_retry_after_and_exponential_cap(monkeypatch):
    monkeypatch.setattr('format_species_network.random.uniform', lambda low, high: high)
    assert download_retry_delay(1, 5) == 5
    assert download_retry_delay(4, 5) == 40
    assert download_retry_delay(100, 5) == 300
    error = HTTPError('http://localhost', 429, '', {'Retry-After': '120'}, None)
    assert download_retry_delay(1, 5, error) == 120


def test_datasets_member_is_accession_scoped():
    names = ['ncbi_dataset/data/GCA_1.1/genomic.fna', 'ncbi_dataset/data/GCA_2.1/genomic.fna']
    assert pick_ncbi_datasets_member_name(names, 'GENOME', 'GCA_2.1') == names[1]
    assert pick_ncbi_datasets_member_name(names, 'GENOME', 'GCA_3.1') == ''
    with pytest.raises(ValueError, match='ambiguous'):
        pick_ncbi_datasets_member_name(names, 'GENOME')


def test_shared_cache_reuses_and_repairs_corruption(tmp_path, monkeypatch):
    monkeypatch.setenv('GG_DOWNLOAD_SHARED_CACHE_DIR', str(tmp_path/'cache'))
    monkeypatch.setenv('GG_DOWNLOAD_EVENT_DIR', str(tmp_path/'events'))
    calls = []
    payload = gzip.compress(b'>seq\nACGT\n')
    def reply(h):
        calls.append(1)
        h.send_response(200)
        h.send_header('Content-Length', str(len(payload)))
        h.end_headers()
        h.wfile.write(payload)
    with server(reply) as url:
        download(url, tmp_path/'plan1'/'data.gz')
        download(url, tmp_path/'plan2'/'data.gz')
        assert len(calls) == 1
        cached = next((tmp_path/'cache').rglob('payload.gz'))
        cached.write_bytes(b'broken')
        download(url, tmp_path/'plan3'/'data.gz')
    assert len(calls) == 2
    assert (tmp_path/'plan1'/'data.gz').read_bytes() == payload
    assert (tmp_path/'plan3'/'data.gz').read_bytes() == payload
    assert len(list((tmp_path/'events').glob('*.json'))) == 3


def test_missing_validator_falls_back_to_full_response(tmp_path, monkeypatch):
    monkeypatch.setenv('GG_DOWNLOAD_RANGE_CHUNK_BYTES', '4')
    seen = []
    def reply(h):
        seen.append(h.headers.get('Range'))
        ranged = h.headers.get('Range') is not None
        h.send_response(206 if ranged else 200)
        if ranged:
            h.send_header('Content-Range', 'bytes 0-3/8')
        body = b'AAAA' if ranged else b'BBBBBBBB'
        h.send_header('Content-Length', str(len(body)))
        h.end_headers()
        h.wfile.write(body)
    with server(reply) as url:
        download(url, tmp_path/'data')
    assert seen == ['bytes=0-3', None]
    assert (tmp_path/'data').read_bytes() == b'BBBBBBBB'


def test_metadata_retry_after_without_limiter(monkeypatch):
    from format_species_network import guarded_urlopen
    delays = []
    monkeypatch.setattr('format_species_network.time.sleep', delays.append)
    calls = []
    def reply(h):
        calls.append(1)
        h.send_response(429 if len(calls) == 1 else 200)
        h.send_header('Retry-After', '2')
        h.end_headers()
        h.wfile.write(b'ok')
    with server(reply) as url, guarded_urlopen(url) as response:
        assert response.read() == b'ok'
    assert delays == [2]


def test_http_403_preserves_identified_partial(tmp_path):
    import hashlib
    import json
    target = tmp_path/'data'
    partial = tmp_path/'data.part'
    calls = []
    def reply(h):
        calls.append(1)
        h.send_response(403)
        h.end_headers()
    with server(reply) as url:
        partial.write_bytes(b'AAAA')
        (tmp_path/'data.part.urlsha256').write_text(hashlib.sha256(url.encode()).hexdigest())
        (tmp_path/'data.part.identity.json').write_text(json.dumps({'etag': '"v1"', 'total': 8}))
        with pytest.raises(HTTPError):
            download(url, target)
    assert partial.read_bytes() == b'AAAA'
    assert len(calls) == 1


def test_corrupt_zip_is_not_cached_as_success(tmp_path, monkeypatch):
    monkeypatch.setenv('GG_DOWNLOAD_SHARED_CACHE_DIR', str(tmp_path/'cache'))
    def reply(h):
        h.send_response(200)
        body = b'PK\x03\x04broken'
        h.send_header('Content-Length', str(len(body)))
        h.end_headers()
        h.wfile.write(body)
    with server(reply) as url, pytest.raises(OSError):
        download(url, tmp_path/'archive.zip')
    assert not (tmp_path/'archive.zip').exists()
    assert not list((tmp_path/'cache').rglob('receipt*.json'))


def test_insufficient_disk_space_does_not_publish(tmp_path, monkeypatch):
    from collections import namedtuple
    usage = namedtuple('usage', 'total used free')
    monkeypatch.setattr('format_species_download.locking.shutil.disk_usage', lambda path: usage(1, 1, 0))
    def reply(h):
        h.send_response(200)
        h.send_header('Content-Length', '4')
        h.end_headers()
        h.wfile.write(b'ACGT')
    with server(reply) as url, pytest.raises(OSError) as error:
        download(url, tmp_path/'data')
    assert error.value.errno == 28
    assert not (tmp_path/'data').exists()


def test_overwrite_corrupt_download_preserves_previous_file(tmp_path):
    target = tmp_path/'data.gz'
    original = gzip.compress(b'previous valid file')
    target.write_bytes(original)
    def reply(h):
        h.send_response(200)
        h.send_header('Content-Length', '7')
        h.end_headers()
        h.wfile.write(b'invalid')
    with server(reply) as url, pytest.raises(OSError):
        download_url_to_file(url, target, {}, 3, False, True, 60, [], 'test')
    assert target.read_bytes() == original


def test_zero_byte_gzip_rejected(tmp_path):
    def reply(h):
        h.send_response(200)
        h.send_header('Content-Length', '0')
        h.end_headers()
    with server(reply) as url, pytest.raises(OSError):
        download(url, tmp_path/'data.gz')
    assert not (tmp_path/'data.gz').exists()


def test_etag_disappearing_does_not_mix_generations(tmp_path, monkeypatch):
    monkeypatch.setenv('GG_DOWNLOAD_RANGE_CHUNK_BYTES', '4')
    calls = []
    def reply(h):
        calls.append(1)
        ranged = h.headers.get('Range')
        start = int(ranged[6:].split('-')[0]) if ranged else 0
        body = b'AAAA' if len(calls) == 1 else (b'BBBB' if ranged else b'BBBBBBBB')
        h.send_response(206 if ranged else 200)
        if ranged:
            h.send_header('Content-Range', f'bytes {start}-{start+3}/8')
        if len(calls) == 1:
            h.send_header('ETag', '"old"')
        h.send_header('Last-Modified', 'Wed, 16 Sep 2026 00:00:00 GMT')
        h.send_header('Content-Length', str(len(body)))
        h.end_headers()
        h.wfile.write(body)
    with server(reply) as url:
        download(url, tmp_path/'data')
    assert (tmp_path/'data').read_bytes() == b'BBBBBBBB'


def test_enospc_preserves_valid_partial(tmp_path, monkeypatch):
    import json

    from format_species_download import locking
    def interrupted(url, partial, *args):
        partial.write_bytes(b'valid partial')
        Path(str(partial)+'.identity.json').write_text(json.dumps({'etag': '"v1"', 'total': 99}))
        raise OSError(28, 'disk full')
    monkeypatch.setattr(locking, 'download_url_to_partial', interrupted)
    with pytest.raises(OSError):
        download('http://localhost/data', tmp_path/'data')
    assert (tmp_path/'data.part').read_bytes() == b'valid partial'


def test_shared_concurrent_plans_download_once(tmp_path, monkeypatch):
    import json
    from concurrent.futures import ThreadPoolExecutor
    monkeypatch.setenv('GG_DOWNLOAD_SHARED_CACHE_DIR', str(tmp_path/'cache'))
    monkeypatch.setenv('GG_DOWNLOAD_EVENT_DIR', str(tmp_path/'events'))
    monkeypatch.setenv('GG_DOWNLOAD_LOCK_POLL_SECONDS', '0.01')
    calls = []
    def reply(h):
        calls.append(1)
        h.send_response(200)
        h.send_header('Content-Length', '4')
        h.end_headers()
        h.wfile.write(b'ACGT')
    with server(reply) as url, ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(download, url, tmp_path/f'plan{i}'/'data') for i in range(2)]
        assert [f.result(timeout=10) for f in futures] == [True, True]
    assert len(calls) == 1
    events = [json.loads(p.read_text()) for p in (tmp_path/'events').glob('*.json')]
    assert sorted(event['status'] for event in events) == ['downloaded', 'materialized']


def test_shared_receipt_unknown_schema_requires_refresh(tmp_path, monkeypatch):
    import json
    monkeypatch.setenv('GG_DOWNLOAD_SHARED_CACHE_DIR', str(tmp_path/'cache'))
    calls = []
    def reply(h):
        calls.append(1)
        h.send_response(200)
        h.send_header('Content-Length', '4')
        h.end_headers()
        h.wfile.write(b'ACGT')
    with server(reply) as url:
        download(url, tmp_path/'one')
        receipt = next((tmp_path/'cache').rglob('receipt*.json'))
        value = json.loads(receipt.read_text())
        value['schema_version'] = 999
        receipt.write_text(json.dumps(value))
        download(url, tmp_path/'two')
    assert len(calls) == 2


def test_archive_member_corrupt_cache_recovers(tmp_path):
    import io
    import zipfile
    payload = io.BytesIO()
    with zipfile.ZipFile(payload, 'w') as archive:
        archive.writestr('data.fa', '>seq\nACGT\n')
    body = payload.getvalue()
    calls = []
    def reply(h):
        calls.append(1)
        h.send_response(200)
        h.send_header('Content-Length', str(len(body)))
        h.end_headers()
        h.wfile.write(body)
    with server(reply) as url:
        def member(name):
            return download_url_to_file(url, tmp_path/name, {}, 3, False, False, 60, [], 'test', archive_member='data.fa')
        member('one.fa')
        cached = [p for p in (tmp_path/'.archive_cache').iterdir() if p.is_file() and '__' in p.name and not p.name.endswith(('.lock', '.guard'))]
        assert len(cached) == 1
        cached[0].write_bytes(b'PK\x03\x04broken')
        member('two.fa')
    assert len(calls) == 2
    assert (tmp_path/'two.fa').read_text() == '>seq\nACGT\n'


def test_archive_member_validation_keeps_existing_destination(tmp_path):
    import io
    import zipfile
    payload = io.BytesIO()
    with zipfile.ZipFile(payload, 'w') as archive:
        archive.writestr('data.gz', b'bad gzip')
    body = payload.getvalue()
    target = tmp_path/'data.gz'
    previous = gzip.compress(b'old')
    target.write_bytes(previous)
    def reply(h):
        h.send_response(200)
        h.send_header('Content-Length', str(len(body)))
        h.end_headers()
        h.wfile.write(body)
    with server(reply) as url, pytest.raises(OSError):
        download_url_to_file(url, target, {}, 3, False, True, 60, [], 'test', archive_member='data.gz')
    assert target.read_bytes() == previous


@pytest.mark.parametrize('body', [b'\xef\xbb\xbf<!--error--><html>failure</html>', gzip.compress(b'<html>error</html>')])
def test_disguised_html_is_not_data(tmp_path, body):
    def reply(h):
        h.send_response(200)
        h.send_header('Content-Length', str(len(body)))
        h.end_headers()
        h.wfile.write(body)
    with server(reply) as url, pytest.raises(OSError):
        download(url, tmp_path/'1234')
    assert not (tmp_path/'1234').exists()


@pytest.mark.parametrize('headers', [{'Content-Length': '-1'}, {'Content-Length': 'oops'}, {'Content-Range': 'bytes 0-3/8'}])
def test_invalid_full_response_framing(tmp_path, headers):
    def reply(h):
        h.send_response(200)
        for key, value in headers.items():
            h.send_header(key, value)
        h.end_headers()
        h.wfile.write(b'ACGT')
    with server(reply) as url, pytest.raises(ValueError):
        download(url, tmp_path/'data')
    assert not (tmp_path/'data').exists()


def test_shared_cache_rejects_recursive_destination(tmp_path, monkeypatch):
    monkeypatch.setenv('GG_DOWNLOAD_SHARED_CACHE_DIR', str(tmp_path/'cache'))
    with pytest.raises(ValueError, match='outside the shared cache'):
        download('http://localhost/data', tmp_path/'cache'/'output')


@pytest.mark.parametrize('network_failure', [False, True])
def test_ncbi_fallback_distinguishes_body_disconnect_from_disk_error(tmp_path, monkeypatch, network_failure):
    from http.client import IncompleteRead

    from format_species_download import manifest
    calls = []
    def fail(*args, **kwargs):
        raise IncompleteRead(b'', 4) if network_failure else OSError(28, 'disk full')
    def fallback(**kwargs):
        calls.append(kwargs['source_id'])
        kwargs['destination'].write_text('>seq\nACGT\n')
        return True
    monkeypatch.setattr(manifest, 'download_url_to_file', fail)
    monkeypatch.setattr(manifest, 'download_ncbi_datasets_file_from_id', fallback)
    result = manifest.execute_download_target_job({
        'provider': 'ncbi', 'source_id': 'GCA_123.1', 'species_key': 'Test_species',
        'label': 'GENOME', 'url': 'https://ftp.ncbi.nlm.nih.gov/test.fa', 'target': tmp_path/'data.fa',
    }, {}, 3, False, 60)
    assert bool(calls) == network_failure
    assert bool(result['failed']) != network_failure
