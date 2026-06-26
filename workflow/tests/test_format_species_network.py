import sys
from pathlib import Path
from urllib.request import Request

import pytest

SUPPORT_DIR = Path(__file__).resolve().parents[1] / "support"
if str(SUPPORT_DIR) not in sys.path:
    sys.path.insert(0, str(SUPPORT_DIR))

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
