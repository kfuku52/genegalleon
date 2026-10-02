import io
import sys
from pathlib import Path
from urllib.error import HTTPError
from urllib.request import Request

import pytest

SUPPORT_DIR = Path(__file__).resolve().parents[1] / "support"
if str(SUPPORT_DIR) not in sys.path:
    sys.path.insert(0, str(SUPPORT_DIR))

import format_species_network as network  # noqa: E402
from format_species_network import assert_url_allowed_by_test_guard  # noqa: E402


def test_network_guard_allows_local_and_file_urls(monkeypatch):
    monkeypatch.setenv("GG_TEST_ALLOW_ONLY_LOOPBACK_HTTP", "1")

    assert_url_allowed_by_test_guard("http://127.0.0.1:1234/example")
    assert_url_allowed_by_test_guard("http://localhost:1234/example")
    assert_url_allowed_by_test_guard("file:///tmp/example.fa")


def test_network_guard_blocks_external_urls(monkeypatch):
    monkeypatch.setenv("GG_TEST_ALLOW_ONLY_LOOPBACK_HTTP", "1")

    with pytest.raises(RuntimeError, match="External network access is disabled"):
        assert_url_allowed_by_test_guard(Request("https://example.com/data"))


@pytest.mark.parametrize("code,attempts", [(404, 1), (429, 2), (503, 2)])
def test_terminal_http_error_preserves_body_and_retry_closes_consumed_errors(monkeypatch, code, attempts):
    errors = []
    def fail(url, *args, **kwargs):
        error = HTTPError(url, code, "failed", {"Retry-After": "0"}, io.BytesIO(b'{"message":"provider detail"}'))
        errors.append(error)
        raise error
    monkeypatch.setattr(network, "limited_urlopen", fail)
    with pytest.raises(HTTPError) as raised:
        network.guarded_urlopen("http://localhost/test", retry_attempts=attempts)
    assert len(errors) == attempts
    assert all(error.closed for error in errors[:-1])
    with raised.value as terminal:
        assert terminal.read() == b'{"message":"provider detail"}'
    assert raised.value.closed
