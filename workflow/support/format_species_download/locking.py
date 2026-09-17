"""Download runtime implementation: locking."""

import errno
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import tarfile
import tempfile
import time
import uuid
import zipfile
from http.client import IncompleteRead
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import unquote, urlparse
from urllib.request import Request

from format_species_common import (
    is_transient_network_error,
)
from format_species_constants import (
    DEFAULT_DOWNLOAD_ATTEMPTS,
    DEFAULT_DOWNLOAD_LOCK_ACQUIRE_TIMEOUT_SECONDS,
    DEFAULT_DOWNLOAD_LOCK_HEARTBEAT_SECONDS,
    DEFAULT_DOWNLOAD_LOCK_POLL_SECONDS,
    DEFAULT_DOWNLOAD_LOCK_STALE_SECONDS,
    DEFAULT_DOWNLOAD_RETRY_BASE_SECONDS,
    DOWNLOAD_DIAGNOSTIC_KEYS,
)
from format_species_network import download_retry_delay
from format_species_network import guarded_urlopen as urlopen
from shared_lock import acquire_lock, release_lock

from .local import (
    quarantine_corrupt_gzip,
    quarantine_existing_file,
    validate_gzip_with_cache,
)

RAR4_SIGNATURE = b"Rar!\x1a\x07\x00"
RAR5_SIGNATURE = b"Rar!\x1a\x07\x01\x00"
CONTENT_RANGE_PATTERN = re.compile(r"^bytes\s+([0-9]+)-([0-9]+)/([0-9]+|[*])$", re.IGNORECASE)
UNSATISFIED_CONTENT_RANGE_PATTERN = re.compile(r"^bytes\s+[*]/([0-9]+)$", re.IGNORECASE)


def _fsync_directory(path):
    try:
        descriptor = os.open(path, os.O_RDONLY)
    except OSError:
        return
    try:
        os.fsync(descriptor)
    except OSError:
        pass
    finally:
        os.close(descriptor)


def is_rar_archive(path):
    with open(path, "rb") as handle:
        signature = handle.read(len(RAR5_SIGNATURE))
    return signature.startswith(RAR4_SIGNATURE) or signature.startswith(RAR5_SIGNATURE)


def extract_rar_archive_member(archive_path, archive_member, destination):
    bsdtar = shutil.which("bsdtar")
    if bsdtar is None:
        raise RuntimeError(
            "RAR archive extraction requires 'bsdtar'; use the GeneGalleon container runtime"
        )
    with open(destination, "wb") as out:
        completed = subprocess.run(
            [bsdtar, "-xOf", str(archive_path), "--", str(archive_member)],
            stdin=subprocess.DEVNULL,
            stdout=out,
            stderr=subprocess.PIPE,
            check=False,
        )
    if completed.returncode == 0:
        return
    try:
        destination.unlink()
    except FileNotFoundError:
        pass
    except OSError:
        pass
    detail = completed.stderr.decode("utf-8", errors="replace").strip()
    if detail == "":
        detail = "bsdtar exited with status {}".format(completed.returncode)
    raise OSError(
        "failed to extract RAR archive member '{}': {}".format(archive_member, detail)
    )


def resolve_download_lock_stale_seconds():
    raw = os.environ.get("GG_DOWNLOAD_LOCK_STALE_SECONDS", "").strip()
    if raw == "":
        return DEFAULT_DOWNLOAD_LOCK_STALE_SECONDS
    try:
        value = int(raw)
    except ValueError:
        return DEFAULT_DOWNLOAD_LOCK_STALE_SECONDS
    if value < 1:
        return 1
    return value


def resolve_download_lock_heartbeat_seconds():
    raw = os.environ.get("GG_DOWNLOAD_LOCK_HEARTBEAT_SECONDS", "").strip()
    if raw == "":
        return DEFAULT_DOWNLOAD_LOCK_HEARTBEAT_SECONDS
    try:
        value = int(raw)
    except ValueError:
        return DEFAULT_DOWNLOAD_LOCK_HEARTBEAT_SECONDS
    if value < 1:
        return 1
    return value


def resolve_download_attempts():
    raw = os.environ.get("GG_DOWNLOAD_ATTEMPTS", "").strip()
    if raw == "":
        return DEFAULT_DOWNLOAD_ATTEMPTS
    try:
        value = int(raw)
    except ValueError:
        return DEFAULT_DOWNLOAD_ATTEMPTS
    if value < 1:
        return 1
    return value


def resolve_download_retry_base_seconds():
    raw = os.environ.get("GG_DOWNLOAD_RETRY_BASE_SECONDS", "").strip()
    if raw == "":
        return DEFAULT_DOWNLOAD_RETRY_BASE_SECONDS
    try:
        value = float(raw)
    except ValueError:
        return DEFAULT_DOWNLOAD_RETRY_BASE_SECONDS
    if value < 0:
        return 0.0
    return value


def sleep_before_download_retry(attempt, base_seconds, error=None):
    delay = download_retry_delay(attempt, base_seconds, error)
    if delay > 0:
        time.sleep(delay)


def partial_download_paths(target):
    partial_path = Path(str(target) + ".part")
    url_hash_path = Path(str(partial_path) + ".urlsha256")
    return partial_path, url_hash_path


def discard_partial_download(partial_path, url_hash_path):
    for path in (partial_path, url_hash_path, Path(str(partial_path) + ".identity.json"),
                 Path(str(partial_path) + ".identity.json.tmp")):
        try:
            path.unlink()
        except FileNotFoundError:
            pass
        except OSError:
            pass


def prepare_partial_download(target, url, overwrite):
    partial_path, url_hash_path = partial_download_paths(target)
    expected_url_hash = hashlib.sha256(str(url).encode("utf-8")).hexdigest()
    if overwrite:
        discard_partial_download(partial_path, url_hash_path)
    elif partial_path.exists():
        try:
            stored_url_hash = url_hash_path.read_text(encoding="utf-8").strip()
        except (FileNotFoundError, OSError):
            stored_url_hash = ""
        if stored_url_hash != expected_url_hash:
            discard_partial_download(partial_path, url_hash_path)
    elif url_hash_path.exists():
        discard_partial_download(partial_path, url_hash_path)
    url_hash_path.write_text(expected_url_hash + "\n", encoding="utf-8")
    return partial_path, url_hash_path


def parse_content_range(value):
    match = CONTENT_RANGE_PATTERN.match(str(value or "").strip())
    if match is None:
        return None, None
    start = int(match.group(1))
    total_text = match.group(3)
    total = None if total_text == "*" else int(total_text)
    return start, total


def parse_unsatisfied_content_range_total(value):
    match = UNSATISFIED_CONTENT_RANGE_PATTERN.match(str(value or "").strip())
    if match is None:
        return None
    return int(match.group(1))


def response_status_code(response):
    status = getattr(response, "status", None)
    if status is None:
        status = response.getcode()
    try:
        return int(status)
    except (TypeError, ValueError):
        return 0


def response_content_length(response):
    raw = response.headers.get("Content-Length")
    try:
        value = int(raw)
    except (TypeError, ValueError):
        if raw is not None:
            raise ValueError("invalid Content-Length") from None
        return None
    if value < 0:
        raise ValueError("negative Content-Length")
    return value


def download_url_to_partial(url, partial_path, headers, timeout, warnings, lock_context):
    """Download sequential ranges, preserving only bytes from one representation."""
    identity_path = Path(str(partial_path) + ".identity.json")
    try:
        identity = json.loads(identity_path.read_text())
        if (not isinstance(identity, dict)
                or not isinstance(identity.get("etag", ""), str)
                or not isinstance(identity.get("last_modified", ""), str)
                or (identity.get("total") is not None and (type(identity["total"]) is not int or identity["total"] < 0))):
            identity = {}
    except (OSError, ValueError):
        identity = {}
    if partial_path.exists() and partial_path.stat().st_size and not (identity.get("etag") or identity.get("last_modified")):
        warnings.append("{} restarting partial download without a representation validator".format(lock_context))
        partial_path.unlink()
    request_headers = {k: v for k, v in (headers or {}).items()
                       if k.lower() not in ("range", "if-range", "accept-encoding")}
    request_headers["Accept-Encoding"] = "identity"
    # Bound CNGB connection lifetimes; other providers retain streaming by default.
    host = urlparse(url).hostname or ""
    default_chunk = 128 * 1024 * 1024 if host == "cngb.org" or host.endswith(".cngb.org") else 0
    try:
        chunk_size = max(0, int(os.environ.get("GG_DOWNLOAD_RANGE_CHUNK_BYTES", default_chunk)))
    except ValueError:
        chunk_size = default_chunk
    while True:
        resume_from = partial_path.stat().st_size if partial_path.exists() else 0
        current_headers = dict(request_headers)
        if resume_from or chunk_size:
            end = str(resume_from + chunk_size - 1) if chunk_size else ""
            current_headers["Range"] = "bytes={}-{}".format(resume_from, end)
        validator = identity.get("etag") or identity.get("last_modified")
        if resume_from and validator:
            current_headers["If-Range"] = validator
        try:
            # The outer file loop owns retries, including errors while reading bodies.
            response_context = urlopen(Request(url, headers=current_headers), timeout=timeout, retry_attempts=1)
        except HTTPError as exc:
            # A size alone cannot establish that a partial file is still current.
            if exc.code == 416 and resume_from:
                exc.close()
                partial_path.unlink(missing_ok=True)
                identity_path.unlink(missing_ok=True)
                identity = {}
                raise IncompleteRead(b"", 1) from exc
            raise
        with response_context as response:
            status = response_status_code(response)
            content_type = response.headers.get("Content-Type", "").split(";", 1)[0].lower()
            if content_type in ("text/html", "application/xhtml+xml"):
                raise ValueError("download returned an HTML page instead of a data file")
            if response.headers.get("Content-Encoding", "identity").lower() not in ("", "identity"):
                raise ValueError("download returned unsupported Content-Encoding")
            expected = response_content_length(response)
            total = expected
            append = False
            if status == 206:
                match = CONTENT_RANGE_PATTERN.fullmatch(response.headers.get("Content-Range", "").strip())
                if match is None or match.group(3) == "*":
                    raise ValueError("invalid or unknown-total Content-Range")
                start, end, total = map(int, match.groups())
                if start != resume_from or end < start or end >= total or (expected is not None and expected != end - start + 1):
                    raise ValueError("inconsistent Content-Range")
                expected = end - start + 1
                if chunk_size and expected > chunk_size:
                    raise ValueError("response exceeded requested Range")
                new_validator = response.headers.get("ETag") if identity.get("etag") else response.headers.get("Last-Modified")
                if resume_from and ((validator and validator != new_validator)
                                    or (identity.get("total") is not None and identity["total"] != total)):
                    partial_path.unlink(missing_ok=True)
                    identity_path.unlink(missing_ok=True)
                    raise IncompleteRead(b"", total)
                append = resume_from > 0
            elif response.headers.get("Content-Range"):
                raise ValueError("Content-Range on a non-partial response")
            elif status not in (0, 200):
                raise ValueError("unexpected download response status {}".format(status))
            response_etag = response.headers.get("ETag", "")
            response_validator = (response_etag and not response_etag.startswith("W/")) or response.headers.get("Last-Modified")
            if status == 206 and "Range" not in current_headers:
                raise ValueError("unsolicited partial response")
            if status == 206 and not response_validator:
                # Multiple requests cannot safely share bytes without a validator.
                partial_path.unlink(missing_ok=True)
                identity_path.unlink(missing_ok=True)
                chunk_size = 0
                identity = {}
                continue
            if not append:
                if resume_from:
                    warnings.append("{} server did not honor Range; restarting the download".format(lock_context))
                resume_from = 0
            else:
                warnings.append("{} resuming partial download at byte {}".format(lock_context, resume_from))
            # Check writable storage and space before consuming a potentially large body.
            required = expected or 0
            if required > shutil.disk_usage(partial_path.parent).free:
                raise OSError(28, "insufficient space for download response")
            etag = response.headers.get("ETag", "")
            identity = {"etag": etag if not etag.startswith("W/") else "",
                        "last_modified": response.headers.get("Last-Modified", ""), "total": total}
            # Truncate before publishing a replacement identity.
            with open(partial_path, "ab" if append else "wb") as out:
                identity_tmp = Path(str(identity_path) + ".tmp")
                identity_tmp.write_text(json.dumps(identity))
                identity_tmp.replace(identity_path)
                response_bytes = 0
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    if expected is not None and response_bytes + len(chunk) > expected:
                        raise ValueError("response exceeded declared length")
                    out.write(chunk)
                    response_bytes += len(chunk)
                out.flush()
                os.fsync(out.fileno())
            if expected is not None and response_bytes != expected:
                raise IncompleteRead(b"", max(0, expected - response_bytes))
            if status == 206 and partial_path.stat().st_size < total:
                continue
            return


def normalize_uncompressed_gzip(partial_path, warnings, lock_context):
    """Compress a complete plain-text response when the target is ``.gz``.

    Public download endpoints sometimes return an uncompressed FASTA/GFF body
    even though the manifest has no usable filename and the normalized target
    is assigned a ``.gz`` suffix.  Keep the target contract stable by
    normalizing only after the complete response has been downloaded; a gzip
    response is left byte-for-byte unchanged.
    """
    partial_path = Path(partial_path)
    if partial_path.stat().st_size == 0:
        raise OSError("downloaded file is empty")
    with open(partial_path, "rb") as source:
        if source.read(2) == b"\x1f\x8b":
            return False
    temporary = Path("{}.gziptmp.{}".format(partial_path, os.getpid()))
    try:
        with open(partial_path, "rb") as source, gzip.open(temporary, "wb") as compressed:
            shutil.copyfileobj(source, compressed, length=1024 * 1024)
            compressed.flush()
        with open(temporary, "rb") as handle:
            os.fsync(handle.fileno())
        os.replace(temporary, partial_path)
        _fsync_directory(partial_path.parent)
    except Exception:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass
        except OSError:
            pass
        raise
    warnings.append("{} normalized uncompressed response to gzip: {}".format(lock_context, partial_path))
    return True


def finalize_partial_download(partial_path, url_hash_path, destination):
    partial_path.replace(destination)
    discard_partial_download(partial_path, url_hash_path)


def empty_download_diagnostics():
    return {key: 0 for key in DOWNLOAD_DIAGNOSTIC_KEYS}


def scan_download_cache_diagnostics(download_root):
    counts = {
        "partial_tmp": 0,
        "corrupt": 0,
        "locks": 0,
    }
    root = Path(download_root)
    if not root.exists():
        return counts
    try:
        iterator = root.rglob("*")
        for path in iterator:
            if ".gg-gzip-validation" in path.parts:
                continue
            try:
                if not path.is_file():
                    continue
            except OSError:
                continue
            name = path.name
            if ".tmp." in name or name.endswith(".part"):
                counts["partial_tmp"] += 1
            if ".corrupt." in name:
                counts["corrupt"] += 1
            if name.endswith(".lock"):
                counts["locks"] += 1
    except OSError:
        return counts
    return counts


def summarize_download_diagnostics(
    preexisting_cache,
    final_cache,
    warnings,
    download_jobs_count,
    failed_downloads_count,
    validation_diagnostics=None,
):
    diagnostics = empty_download_diagnostics()
    diagnostics["cache_preexisting_partial_tmp"] = int(preexisting_cache.get("partial_tmp", 0))
    diagnostics["cache_preexisting_corrupt"] = int(preexisting_cache.get("corrupt", 0))
    diagnostics["cache_preexisting_locks"] = int(preexisting_cache.get("locks", 0))
    diagnostics["cache_final_partial_tmp"] = int(final_cache.get("partial_tmp", 0))
    diagnostics["cache_final_corrupt"] = int(final_cache.get("corrupt", 0))
    diagnostics["cache_final_locks"] = int(final_cache.get("locks", 0))
    diagnostics["download_jobs"] = int(download_jobs_count)
    diagnostics["failed_downloads"] = int(failed_downloads_count)
    diagnostics.update({key: int(value) for key, value in (validation_diagnostics or {}).items()
                        if key in diagnostics})
    for warning in warnings:
        text = str(warning or "").lower()
        if "failed transiently; retrying" in text:
            diagnostics["transient_retries"] += 1
        if "resuming partial download at byte" in text:
            diagnostics["range_resumes"] += 1
        if "downloaded corrupt gzip" in text and "retrying" in text:
            diagnostics["corrupt_download_retries"] += 1
        if "found corrupt gzip cache" in text:
            diagnostics["corrupt_cache_recoveries"] += 1
        if "[download-lock] recovered stale lock" in text:
            diagnostics["stale_locks_recovered"] += 1
        if "[download-lock] waiting for shared lock" in text:
            diagnostics["lock_waits"] += 1
    return diagnostics


def format_download_diagnostics_line(diagnostics):
    values = empty_download_diagnostics()
    values.update(diagnostics or {})
    return (
        "Download diagnostics: "
        "cache_preexisting partial_tmp={cache_preexisting_partial_tmp},corrupt={cache_preexisting_corrupt},locks={cache_preexisting_locks}; "
        "cache_final partial_tmp={cache_final_partial_tmp},corrupt={cache_final_corrupt},locks={cache_final_locks}; "
        "download_jobs={download_jobs}; "
        "failed_downloads={failed_downloads}; "
        "retries transient={transient_retries},corrupt_gzip={corrupt_download_retries},range_resumes={range_resumes}; "
        "corrupt_cache_recoveries={corrupt_cache_recoveries}; "
        "validation_cache hits={validation_cache_hits},misses={validation_cache_misses},records={validation_cache_records}; "
        "stale_locks recovered={stale_locks_recovered},waits={lock_waits}"
    ).format(**values)


def resolve_download_lock_acquire_timeout_seconds():
    raw = os.environ.get("GG_DOWNLOAD_LOCK_ACQUIRE_TIMEOUT_SECONDS", "").strip()
    if raw == "":
        return DEFAULT_DOWNLOAD_LOCK_ACQUIRE_TIMEOUT_SECONDS
    try:
        value = int(raw)
    except ValueError:
        return DEFAULT_DOWNLOAD_LOCK_ACQUIRE_TIMEOUT_SECONDS
    if value < 1:
        return 1
    return value


def resolve_download_lock_poll_seconds():
    raw = os.environ.get("GG_DOWNLOAD_LOCK_POLL_SECONDS", "").strip()
    if raw == "":
        return DEFAULT_DOWNLOAD_LOCK_POLL_SECONDS
    try:
        value = float(raw)
    except ValueError:
        return DEFAULT_DOWNLOAD_LOCK_POLL_SECONDS
    if value <= 0:
        return 0.1
    return value


def acquire_download_lock(lock_path, stale_seconds, warnings, lock_context):
    return acquire_lock(
        lock_path,
        stale_seconds=stale_seconds,
        timeout_seconds=resolve_download_lock_acquire_timeout_seconds(),
        poll_seconds=resolve_download_lock_poll_seconds(),
        context=lock_context,
        warning_callback=warnings.append,
        message_label="[download-lock]",
        heartbeat_interval_seconds=resolve_download_lock_heartbeat_seconds(),
    )


def release_download_lock(lock_path, ownership):
    release_lock(lock_path, ownership)


def _download_url_to_file(
    url,
    destination,
    headers,
    timeout,
    dry_run,
    overwrite,
    lock_stale_seconds,
    warnings,
    lock_context,
    archive_member="",
    validation_cache=None,
    validation_key="",
    validation_relative_target="",
):
    if dry_run:
        return False
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    # A real create catches read-only mounts and ACLs before requesting a body.
    with tempfile.TemporaryFile(dir=destination.parent):
        pass
    request_headers = dict(headers or {})
    if not any(key.lower() == "user-agent" for key in request_headers):
        host = urlparse(url).hostname or ""
        request_headers["User-Agent"] = ("Mozilla/5.0 (compatible; genegalleon-input-generation)"
                                         if host == "figshare.com" or host.endswith(".figshare.com")
                                         else "genegalleon-input-generation")
    archive_member_text = str(archive_member or "").strip()
    archive_cache_path = None
    if archive_member_text == "":
        lock_path = Path(str(destination) + ".lock")
    else:
        archive_name = Path(unquote(urlparse(url).path)).name
        if archive_name == "":
            archive_name = "archive.bin"
        archive_hash = hashlib.sha256(json.dumps([str(url), sorted((k.lower(), v) for k, v in request_headers.items())]).encode()).hexdigest()
        archive_cache_dir = destination.parent / ".archive_cache"
        archive_cache_dir.mkdir(parents=True, exist_ok=True)
        archive_cache_path = archive_cache_dir / "{}__{}".format(archive_hash, archive_name)
        lock_path = Path(str(destination) + ".lock")
    ownership = acquire_download_lock(lock_path, lock_stale_seconds, warnings, lock_context)
    tmp = Path(str(destination) + ".tmp.{}".format(os.getpid()))
    try:
        if destination.exists() and destination.stat().st_size > 0 and not overwrite:
            if quarantine_corrupt_gzip(
                destination,
                warnings,
                lock_context,
                validation_cache=validation_cache,
                validation_key=validation_key,
                source_url=url,
                archive_member=archive_member_text,
                relative_target=validation_relative_target,
            ):
                pass
            else:
                return False
        if destination.exists() and destination.stat().st_size > 0 and not overwrite:
            return False
        attempts = resolve_download_attempts()
        retry_base_seconds = resolve_download_retry_base_seconds()
        if archive_member_text == "":
            partial_path, partial_url_hash_path = prepare_partial_download(destination, url, overwrite)
            last_validation_error = None
            for attempt in range(1, attempts + 1):
                try:
                    download_url_to_partial(
                        url,
                        partial_path,
                        request_headers,
                        timeout,
                        warnings,
                        lock_context,
                    )
                    if destination.name.lower().endswith(".gz"):
                        normalize_uncompressed_gzip(partial_path, warnings, lock_context)
                    validation_error = validate_gzip_with_cache(
                        partial_path,
                        expected_path=destination,
                        validation_cache=validation_cache,
                        validation_key=validation_key,
                        source_url=url,
                        archive_member=archive_member_text,
                        relative_target=validation_relative_target,
                    )
                    if validation_error is None:
                        last_validation_error = None
                        finalize_partial_download(partial_path, partial_url_hash_path, destination)
                        _fsync_directory(destination.parent)
                        break
                    last_validation_error = validation_error
                    quarantined = quarantine_existing_file(partial_path, warnings, lock_context, validation_error)
                    if attempt < attempts:
                        warnings.append(
                            "{} downloaded corrupt gzip to {}; retrying attempt {}/{} ({})".format(
                                lock_context,
                                quarantined,
                                attempt + 1,
                                attempts,
                                validation_error,
                            )
                        )
                        sleep_before_download_retry(attempt, retry_base_seconds)
                        continue
                    raise OSError("downloaded gzip failed integrity check: {}".format(validation_error))
                except Exception as exc:
                    transient_error = is_transient_network_error(exc)
                    if attempt < attempts and transient_error:
                        warnings.append(
                            "{} download attempt {}/{} failed transiently; retrying ({})".format(
                                lock_context,
                                attempt,
                                attempts,
                                exc,
                            )
                        )
                        sleep_before_download_retry(attempt, retry_base_seconds, exc)
                        continue
                    if not transient_error and not isinstance(exc, HTTPError) and getattr(exc, "errno", None) not in (errno.ENOSPC, errno.EDQUOT):
                        discard_partial_download(partial_path, partial_url_hash_path)
                    raise
            if last_validation_error is not None:
                raise OSError("downloaded gzip failed integrity check: {}".format(last_validation_error))
        else:
            # Reuse the same validation/retry path for archives themselves.
            _download_url_to_file(
                url, archive_cache_path, request_headers, timeout, False, overwrite,
                lock_stale_seconds, warnings, lock_context + " archive",
            )
            archive_lock = Path(str(archive_cache_path) + ".lock")
            archive_ownership = acquire_download_lock(archive_lock, lock_stale_seconds, warnings, lock_context)
            try:
                if zipfile.is_zipfile(archive_cache_path):
                    with zipfile.ZipFile(archive_cache_path) as archive:
                        with archive.open(archive_member_text, "r") as source, open(
                            tmp, "wb"
                        ) as out:
                            shutil.copyfileobj(source, out, length=1024 * 1024)
                            out.flush()
                            os.fsync(out.fileno())
                elif is_rar_archive(archive_cache_path):
                    extract_rar_archive_member(archive_cache_path, archive_member_text, tmp)
                    with open(tmp, "rb") as handle:
                        os.fsync(handle.fileno())
                else:
                    with tarfile.open(archive_cache_path, "r:*") as archive:
                        extracted = archive.extractfile(archive_member_text)
                        if extracted is None:
                            raise KeyError(archive_member_text)
                        with extracted, open(tmp, "wb") as out:
                            shutil.copyfileobj(extracted, out, length=1024 * 1024)
                            out.flush()
                            os.fsync(out.fileno())
            except (zipfile.BadZipFile, tarfile.ReadError, EOFError) as exc:
                quarantine_existing_file(archive_cache_path, warnings, lock_context, exc)
                raise
            finally:
                release_download_lock(archive_lock, archive_ownership)
            validation_error = validate_gzip_with_cache(
                tmp,
                expected_path=destination,
                validation_cache=validation_cache,
                validation_key=validation_key,
                source_url=url,
                archive_member=archive_member_text,
                relative_target=validation_relative_target,
            )
            if validation_error is not None:
                quarantine_existing_file(tmp, warnings, lock_context, validation_error)
                raise OSError("downloaded archive member gzip failed integrity check: {}".format(validation_error))
            os.replace(tmp, destination)
            _fsync_directory(destination.parent)
    except Exception:
        try:
            tmp.unlink()
        except FileNotFoundError:
            pass
        except OSError:
            pass
        raise
    finally:
        release_download_lock(lock_path, ownership)
    return True


def download_url_to_file(
    url, destination, headers, timeout, dry_run, overwrite, lock_stale_seconds,
    warnings, lock_context, archive_member="", validation_cache=None,
    validation_key="", validation_relative_target="",
):
    """Download with optional cross-plan cache and private structured outcomes."""
    from format_species_network import request_database

    from .shared_cache import materialize_shared_download
    arguments = dict(
        url=url, destination=destination, headers=headers, timeout=timeout,
        dry_run=dry_run, overwrite=overwrite, lock_stale_seconds=lock_stale_seconds,
        warnings=warnings, lock_context=lock_context, archive_member=archive_member,
        validation_cache=validation_cache, validation_key=validation_key,
        validation_relative_target=validation_relative_target,
    )
    started = time.monotonic()
    warning_start = len(warnings)
    error = None
    changed = False
    fetched = True
    try:
        if not dry_run and os.environ.get("GG_DOWNLOAD_SHARED_CACHE_DIR", "").strip():
            Path(destination).parent.mkdir(parents=True, exist_ok=True)
            destination_lock = Path(str(destination) + ".lock")
            ownership = acquire_download_lock(destination_lock, lock_stale_seconds, warnings, lock_context)
            try:
                changed, fetched = materialize_shared_download(_download_url_to_file, arguments)
            finally:
                release_download_lock(destination_lock, ownership)
        else:
            changed = _download_url_to_file(**arguments)
        return changed
    except Exception as exc:
        error = exc
        raise
    finally:
        event_directory = os.environ.get("GG_DOWNLOAD_EVENT_DIR", "").strip()
        if event_directory and not dry_run:
            event = {
                "schema_version": 1,
                "database": request_database(url),
                "url_sha256": hashlib.sha256(str(url).encode()).hexdigest(),
                "status": "failed" if error else ("downloaded" if changed and fetched else ("materialized" if changed else "reused")),
                "http_status": error.code if isinstance(error, HTTPError) else None,
                "error_class": type(error).__name__ if error else None,
                "elapsed_seconds": round(time.monotonic() - started, 3),
                "retries": sum("retrying" in str(w) for w in warnings[warning_start:]),
            }
            try:
                event["bytes"] = Path(destination).stat().st_size if not error else 0
                root = Path(event_directory)
                root.mkdir(parents=True, exist_ok=True)
                # Independent files avoid append races between scheduler processes.
                name = root / (uuid.uuid4().hex + ".json")
                temporary = name.with_suffix(".tmp")
                with open(temporary, "x", encoding="utf-8") as out:
                    json.dump(event, out)
                    out.write("\n")
                    out.flush()
                    os.fsync(out.fileno())
                temporary.replace(name)
            except OSError as log_error:
                warnings.append("{} could not write download event ({})".format(lock_context, type(log_error).__name__))
