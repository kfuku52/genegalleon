import multiprocessing
import os
import sys
import time
from pathlib import Path
from queue import Empty

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT))
from format_species_network import guarded_urlopen, retry_after_seconds  # noqa: E402
from input_download_limiter import Admission, database_key  # noqa: E402


def _hold(root, queue):
    os.environ["GG_INPUT_DOWNLOAD_LIMIT_DIR"] = root
    os.environ["GG_INPUT_MAX_CONCURRENT_DOWNLOADS_NCBI"] = "2"
    os.environ["GG_INPUT_REQUEST_INTERVAL_NCBI"] = "0.12"
    permit = Admission("https://ftp.ncbi.nlm.nih.gov/test").acquire()
    queue.put((os.getpid(), time.monotonic()))
    time.sleep(60)
    permit.close()


def test_shared_slots_rate_spacing_and_process_death(tmp_path):
    context = multiprocessing.get_context("spawn")
    queue = context.Queue()
    workers = [context.Process(target=_hold, args=(str(tmp_path), queue)) for _ in range(3)]
    try:
        for worker in workers:
            worker.start()
        first, second = queue.get(timeout=10), queue.get(timeout=10)
        assert abs(second[1] - first[1]) >= 0.10
        with pytest.raises(Empty):
            queue.get(timeout=0.25)
        holder = next(worker for worker in workers if worker.pid == first[0])
        holder.kill()
        holder.join(5)
        third = queue.get(timeout=10)
        assert third[0] not in (first[0], second[0])
    finally:
        for worker in workers:
            if worker.is_alive():
                worker.kill()
            worker.join(5)


def test_shared_cooldown_and_mismatched_policy(tmp_path, monkeypatch):
    monkeypatch.setenv("GG_INPUT_DOWNLOAD_LIMIT_DIR", str(tmp_path))
    monkeypatch.setenv("GG_INPUT_REQUEST_INTERVAL_NCBI", "0")
    first = Admission("https://api.ncbi.nlm.nih.gov/").acquire()
    first.cooldown(0.2)
    first.close()
    start = time.monotonic()
    Admission("https://ftp.ncbi.nlm.nih.gov/").acquire().close()
    assert time.monotonic() - start >= 0.18
    monkeypatch.setenv("GG_INPUT_MAX_CONCURRENT_DOWNLOADS_NCBI", "7")
    with pytest.raises(ValueError, match="policy mismatch"):
        Admission("https://www.ncbi.nlm.nih.gov/").acquire()


def test_logical_database_aliases_and_direct_hosts():
    assert database_key("https://ftp.ncbi.nlm.nih.gov/a") == database_key("https://api.ncbi.nlm.nih.gov/b")
    assert database_key("https://plants.ensembl.org/a") == database_key("https://ftp.ensemblgenomes.org/b")
    assert database_key("https://example.org/") != database_key("https://other.org/")
    assert database_key("https://notncbi.nlm.nih.gov/") != "ncbi"
    assert retry_after_seconds("2", 0) == 2
    assert retry_after_seconds("invalid", 2) == 4


def test_response_holds_permit_until_close(tmp_path, monkeypatch):
    import format_species_network as network
    monkeypatch.setenv("GG_INPUT_DOWNLOAD_LIMIT_DIR", str(tmp_path))
    monkeypatch.setenv("GG_INPUT_MAX_CONCURRENT_DOWNLOADS_DIRECT", "1")
    monkeypatch.setenv("GG_INPUT_REQUEST_INTERVAL_DIRECT", "0")
    monkeypatch.setenv("GG_INPUT_DOWNLOAD_LIMIT_WAIT", "0.1")
    class Response:
        def read(self):
            return b"stream"
        def close(self):
            pass
    monkeypatch.setattr(network, "limited_urlopen", lambda request, *args, **kwargs: network.open_with_permit(lambda req: Response(), request))
    with guarded_urlopen("http://localhost/test") as response:
        assert response.read() == b"stream"
        with pytest.raises(TimeoutError):
            Admission("http://localhost/test").acquire()
    Admission("http://localhost/test").acquire().close()


def test_http_redirect_and_retry_after_share_destination_bucket(tmp_path, monkeypatch):
    import threading
    from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

    from format_species_network import guarded_urlopen
    monkeypatch.setenv("GG_INPUT_DOWNLOAD_LIMIT_DIR", str(tmp_path))
    monkeypatch.setenv("GG_INPUT_REQUEST_INTERVAL_DIRECT", "0")
    # Use separate loopback hosts to test that redirects acquire a new bucket.
    requests = []
    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            requests.append((self.path, time.monotonic()))
            if self.path == "/redirect":
                self.send_response(302)
                self.send_header("Location", "http://127.0.0.1:" + str(self.server.server_port) + "/target")
                self.end_headers()
            elif len([path for path, _ in requests if path == "/target"]) == 1:
                self.send_response(429)
                self.send_header("Retry-After", "0.2")
                self.end_headers()
            else:
                self.send_response(200)
                self.end_headers()
                self.wfile.write(b"complete")
        def log_message(self, *args):
            pass
    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        with guarded_urlopen("http://localhost:" + str(server.server_port) + "/redirect", timeout=3) as response:
            assert response.read() == b"complete"
        target_times = [started for path, started in requests if path == "/target"]
        assert len(target_times) == 2
        assert target_times[1] - target_times[0] >= 0.18
        import json
        states = [json.loads(path.read_text()) for path in tmp_path.glob("*/state.json")]
        assert {state["policy"]["database"] for state in states} == {"host:localhost", "host:127.0.0.1"}
        assert all(not list(path.parent.glob("*.tmp")) for path in tmp_path.glob("*/state.json"))
    finally:
        server.shutdown()
        server.server_close()
        thread.join(3)


def test_admission_metadata_lock_obeys_wait_timeout(tmp_path, monkeypatch):
    import fcntl
    monkeypatch.setenv("GG_INPUT_DOWNLOAD_LIMIT_DIR", str(tmp_path))
    monkeypatch.setenv("GG_INPUT_DOWNLOAD_LIMIT_WAIT", "0.1")
    permit = Admission("https://ftp.ncbi.nlm.nih.gov/example")
    with open(permit.directory / "admission.lock", "a+b") as blocker:
        fcntl.flock(blocker, fcntl.LOCK_EX)
        started = time.monotonic()
        with pytest.raises(TimeoutError, match="metadata lock"):
            permit.acquire()
        assert time.monotonic() - started < 1
    permit.acquire().close()


@pytest.mark.parametrize("value", ["inf", "nan", "-inf"])
def test_nonfinite_retry_after_uses_bounded_backoff(value):
    assert retry_after_seconds(value, 2) == 4


def test_provider_cdn_redirect_keeps_logical_database_and_cooldown(tmp_path, monkeypatch):
    import json
    import threading
    from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

    from format_species_network import request_provider
    monkeypatch.setenv("GG_INPUT_DOWNLOAD_LIMIT_DIR", str(tmp_path))
    monkeypatch.setenv("GG_INPUT_REQUEST_INTERVAL_NCBI", "0")
    class Handler(BaseHTTPRequestHandler):
        count = 0
        def do_GET(self):
            if self.path == "/origin":
                self.send_response(302)
                self.send_header("Location", f"http://127.0.0.1:{self.server.server_port}/cdn")
            else:
                Handler.count += 1
                self.send_response(503 if Handler.count == 1 else 200)
                if Handler.count == 1:
                    self.send_header("Retry-After", "0.15")
            self.end_headers()
            self.wfile.write(b"data")
        def log_message(self, *args):
            pass
    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        with request_provider("refseq"):
            with guarded_urlopen(f"http://localhost:{server.server_port}/origin", timeout=2) as response:
                assert response.read() == b"data"
        states = [json.loads(path.read_text()) for path in tmp_path.glob("*/state.json")]
        assert len(states) == 1 and states[0]["policy"]["database"] == "ncbi"
        assert states[0]["cooldown"] > 0
        # Scope must not leak into unrelated direct requests on a reused thread.
        with guarded_urlopen(f"http://localhost:{server.server_port}/cdn", timeout=2) as response:
            assert response.read() == b"data"
        assert len(list(tmp_path.glob("*/state.json"))) == 2
    finally:
        server.shutdown()
        server.server_close()
        thread.join(3)


def test_cross_origin_redirect_does_not_forward_credentials(tmp_path, monkeypatch):
    import threading
    from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
    from urllib.request import Request
    monkeypatch.setenv("GG_INPUT_DOWNLOAD_LIMIT_DIR", str(tmp_path))
    monkeypatch.setenv("GG_INPUT_REQUEST_INTERVAL_DIRECT", "0")
    observed = []
    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            observed.append((self.path, self.headers.get("Authorization"), self.headers.get("Cookie")))
            if self.path == "/origin":
                self.send_response(302)
                self.send_header("Location", f"http://localhost:{self.server.server_port}/same")
            elif self.path == "/same":
                self.send_response(302)
                self.send_header("Location", f"http://127.0.0.1:{self.server.server_port}/other")
            else:
                self.send_response(200)
            self.end_headers()
        def log_message(self, *args):
            pass
    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        request = Request(f"http://localhost:{server.server_port}/origin",
                          headers={"Authorization": "Bearer test-only", "Cookie": "test-only-cookie"})
        with guarded_urlopen(request, timeout=2) as response:
            response.read()
        assert observed == [("/origin", "Bearer test-only", "test-only-cookie"),
                            ("/same", "Bearer test-only", "test-only-cookie"), ("/other", None, None)]
    finally:
        server.shutdown()
        server.server_close()
        thread.join(3)
