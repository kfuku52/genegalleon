import ipaddress
import math
import os
import time
from contextlib import contextmanager
from contextvars import ContextVar
from email.utils import parsedate_to_datetime
from functools import wraps
from urllib.error import HTTPError
from urllib.parse import urlparse
from urllib.request import FTPHandler, HTTPHandler, HTTPRedirectHandler, HTTPSHandler, build_opener
from urllib.request import urlopen as _urlopen

from input_download_limiter import Admission, database_key, limit_directory

NETWORK_GUARD_ENV = "GG_TEST_ALLOW_ONLY_LOOPBACK_HTTP"


_REQUEST_PROVIDER = ContextVar("input_request_provider", default=None)


def normalize_request_provider(provider):
    key = str(provider or "").lower()
    if key in ("refseq", "genbank"):
        key = "ncbi"
    elif key.startswith("ensembl"):
        key = "ensembl"
    return key if key not in ("", "all", "direct", "local") else None


def request_database(url, provider=None):
    destination = database_key(url)
    return destination if not destination.startswith("host:") else normalize_request_provider(provider) or destination


def set_request_provider(provider):
    _REQUEST_PROVIDER.set(normalize_request_provider(provider))


@contextmanager
def request_provider(provider):
    token = _REQUEST_PROVIDER.set(None)
    set_request_provider(provider)
    try:
        yield
    finally:
        _REQUEST_PROVIDER.reset(token)


def isolated_request_provider(function):
    @wraps(function)
    def wrapped(*args, **kwargs):
        with request_provider(None):
            return function(*args, **kwargs)
    return wrapped


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
    if not limit_directory():
        return _urlopen(request_or_url, *args, **kwargs)
    # Transport handlers acquire permits for each hop. urllib closes a redirect
    # response before opening its destination, releasing the original permit.
    for attempt in range(4):
        try:
            return limited_urlopen(request_or_url, *args, **kwargs)
        except HTTPError as exc:
            if exc.code not in (429, 503):
                exc.close()
                raise
            delay = retry_after_seconds(exc.headers.get("Retry-After"), attempt)
            try:
                permit = getattr(exc.fp, "permit", None) or Admission(exc.geturl())
                permit.cooldown(delay)
            finally:
                exc.close()
            if attempt == 3:
                raise
            time.sleep(delay)


def limited_urlopen(request_or_url, *args, **kwargs):
    logical_database = _REQUEST_PROVIDER.get()

    def open_hop(open_request, request):
        nonlocal logical_database
        destination = database_key(request_url(request))
        if not destination.startswith("host:"):
            logical_database = destination
        effective = logical_database or destination
        return open_with_permit(open_request, request, database=effective)

    class SafeRedirect(HTTPRedirectHandler):
        def redirect_request(self, request, response, code, message, headers, newurl):
            redirected = super().redirect_request(request, response, code, message, headers, newurl)
            old = urlparse(request_url(request))
            new = urlparse(newurl)
            old_origin = (old.scheme, old.hostname, old.port or (443 if old.scheme == "https" else 80))
            new_origin = (new.scheme, new.hostname, new.port or (443 if new.scheme == "https" else 80))
            if redirected is not None and old_origin != new_origin:
                for header in ("Authorization", "Cookie", "Proxy-authorization"):
                    redirected.remove_header(header)
            return redirected

    class LimitedHTTP(HTTPHandler):
        def http_open(self, request):
            return open_hop(super().http_open, request)

    class LimitedHTTPS(HTTPSHandler):
        def https_open(self, request):
            return open_hop(super().https_open, request)

    class LimitedFTP(FTPHandler):
        def ftp_open(self, request):
            return open_hop(super().ftp_open, request)

    context = kwargs.pop("context", None)
    opener = build_opener(SafeRedirect(), LimitedHTTP(), LimitedHTTPS(context=context), LimitedFTP())
    return opener.open(request_or_url, *args, **kwargs)


def open_with_permit(open_request, request, database=None):
    assert_url_allowed_by_test_guard(request)
    permit = Admission(request_url(request), database=database).acquire()
    try:
        return LimitedResponse(open_request(request), permit)
    except BaseException:
        permit.close()
        raise


def retry_after_seconds(value, attempt):
    try:
        seconds = float(value)
        return max(0.0, seconds) if math.isfinite(seconds) else float(2 ** attempt)
    except (TypeError, ValueError):
        try:
            return max(0.0, parsedate_to_datetime(value).timestamp() - time.time())
        except (TypeError, ValueError, OverflowError):
            return float(2 ** attempt)


class LimitedResponse:
    def __init__(self, response, permit):
        self.response = response
        self.permit = permit

    def __getattr__(self, name):
        return getattr(self.response, name)

    def __iter__(self):
        return iter(self.response)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()

    def close(self):
        try:
            self.response.close()
        finally:
            self.permit.close()
