import ipaddress
import os
from urllib.parse import urlparse
from urllib.request import urlopen as _urlopen

NETWORK_GUARD_ENV = "GG_TEST_ALLOW_ONLY_LOOPBACK_HTTP"


def test_network_guard_enabled():
    return os.environ.get(NETWORK_GUARD_ENV, "").strip() not in ("", "0", "false", "False")


def request_url(request_or_url):
    if hasattr(request_or_url, "get_full_url"):
        return request_or_url.get_full_url()
    if hasattr(request_or_url, "full_url"):
        return request_or_url.full_url
    return str(request_or_url)


def is_loopback_hostname(hostname):
    host = str(hostname or "").strip()
    if host == "":
        return False
    if host.lower() == "localhost":
        return True
    try:
        return ipaddress.ip_address(host).is_loopback
    except ValueError:
        return False


def assert_url_allowed_by_test_guard(request_or_url):
    if not test_network_guard_enabled():
        return
    url = request_url(request_or_url)
    parsed = urlparse(url)
    if parsed.scheme in ("", "file", "data"):
        return
    if is_loopback_hostname(parsed.hostname):
        return
    raise RuntimeError(
        "External network access is disabled by {}: {}".format(
            NETWORK_GUARD_ENV,
            url,
        )
    )


def guarded_urlopen(request_or_url, *args, **kwargs):
    assert_url_allowed_by_test_guard(request_or_url)
    return _urlopen(request_or_url, *args, **kwargs)
