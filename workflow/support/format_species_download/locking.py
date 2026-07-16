"""Download runtime implementation: locking."""

import hashlib
import os
import tarfile
import time
import zipfile
from pathlib import Path
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
from format_species_network import guarded_urlopen as urlopen
from shared_lock import acquire_lock, release_lock

from .local import (
    gzip_integrity_error,
    quarantine_corrupt_gzip,
    quarantine_existing_file,
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


def sleep_before_download_retry(attempt, base_seconds):
    delay = float(base_seconds) * float(attempt)
    if delay > 0:
        time.sleep(delay)


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
            try:
                if not path.is_file():
                    continue
            except OSError:
                continue
            name = path.name
            if ".tmp." in name:
                counts["partial_tmp"] += 1
            if ".corrupt." in name:
                counts["corrupt"] += 1
            if name.endswith(".lock"):
                counts["locks"] += 1
    except OSError:
        return counts
    return counts


def summarize_download_diagnostics(
    preexisting_cache, final_cache, warnings, download_jobs_count, failed_downloads_count
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
    for warning in warnings:
        text = str(warning or "").lower()
        if "failed transiently; retrying" in text:
            diagnostics["transient_retries"] += 1
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
        "retries transient={transient_retries},corrupt_gzip={corrupt_download_retries}; "
        "corrupt_cache_recoveries={corrupt_cache_recoveries}; "
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


def release_download_lock(lock_path, heartbeat_state=None):
    release_lock(lock_path, heartbeat_state)


def download_url_to_file(
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
):
    if dry_run:
        return False
    request_headers = dict(headers or {})
    if "User-Agent" not in request_headers:
        request_headers["User-Agent"] = "genegalleon-input-generation"
    archive_member_text = str(archive_member or "").strip()
    archive_cache_path = None
    archive_tmp = None
    if archive_member_text == "":
        lock_path = Path(str(destination) + ".lock")
    else:
        archive_name = Path(unquote(urlparse(url).path)).name
        if archive_name == "":
            archive_name = "archive.bin"
        archive_hash = hashlib.sha256(str(url).encode("utf-8")).hexdigest()[:16]
        archive_cache_dir = destination.parent / ".archive_cache"
        archive_cache_dir.mkdir(parents=True, exist_ok=True)
        archive_cache_path = archive_cache_dir / "{}__{}".format(archive_hash, archive_name)
        lock_path = Path(str(archive_cache_path) + ".lock")
    heartbeat_state = acquire_download_lock(lock_path, lock_stale_seconds, warnings, lock_context)
    tmp = Path(str(destination) + ".tmp.{}".format(os.getpid()))
    try:
        if destination.exists() and destination.stat().st_size > 0 and not overwrite:
            if quarantine_corrupt_gzip(destination, warnings, lock_context):
                pass
            else:
                return False
        if destination.exists() and destination.stat().st_size > 0 and not overwrite:
            return False
        attempts = resolve_download_attempts()
        retry_base_seconds = resolve_download_retry_base_seconds()
        if archive_member_text == "":
            last_validation_error = None
            for attempt in range(1, attempts + 1):
                try:
                    request = Request(url, headers=request_headers)
                    with urlopen(request, timeout=timeout) as response, open(tmp, "wb") as out:
                        while True:
                            chunk = response.read(1024 * 1024)
                            if not chunk:
                                break
                            out.write(chunk)
                    tmp.replace(destination)
                    validation_error = gzip_integrity_error(destination)
                    if validation_error is None:
                        last_validation_error = None
                        break
                    last_validation_error = validation_error
                    quarantined = quarantine_existing_file(destination, warnings, lock_context, validation_error)
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
                    try:
                        tmp.unlink()
                    except FileNotFoundError:
                        pass
                    except OSError:
                        pass
                    if attempt < attempts and is_transient_network_error(exc):
                        warnings.append(
                            "{} download attempt {}/{} failed transiently; retrying ({})".format(
                                lock_context,
                                attempt,
                                attempts,
                                exc,
                            )
                        )
                        sleep_before_download_retry(attempt, retry_base_seconds)
                        continue
                    raise
            if last_validation_error is not None:
                raise OSError("downloaded gzip failed integrity check: {}".format(last_validation_error))
        else:
            archive_tmp = Path(str(archive_cache_path) + ".tmp.{}".format(os.getpid()))
            if overwrite or not archive_cache_path.exists() or archive_cache_path.stat().st_size == 0:
                for attempt in range(1, attempts + 1):
                    try:
                        request = Request(url, headers=request_headers)
                        with urlopen(request, timeout=timeout) as response, open(archive_tmp, "wb") as out:
                            while True:
                                chunk = response.read(1024 * 1024)
                                if not chunk:
                                    break
                                out.write(chunk)
                        archive_tmp.replace(archive_cache_path)
                        break
                    except Exception as exc:
                        try:
                            archive_tmp.unlink()
                        except FileNotFoundError:
                            pass
                        except OSError:
                            pass
                        if attempt < attempts and is_transient_network_error(exc):
                            warnings.append(
                                "{} archive download attempt {}/{} failed transiently; retrying ({})".format(
                                    lock_context,
                                    attempt,
                                    attempts,
                                    exc,
                                )
                            )
                            sleep_before_download_retry(attempt, retry_base_seconds)
                            continue
                        raise
            if zipfile.is_zipfile(archive_cache_path):
                with zipfile.ZipFile(archive_cache_path) as archive:
                    payload = archive.read(archive_member_text)
            else:
                with tarfile.open(archive_cache_path, "r:*") as archive:
                    extracted = archive.extractfile(archive_member_text)
                    if extracted is None:
                        raise KeyError(archive_member_text)
                    payload = extracted.read()
            with open(tmp, "wb") as out:
                out.write(payload)
            tmp.replace(destination)
            validation_error = gzip_integrity_error(destination)
            if validation_error is not None:
                quarantine_existing_file(destination, warnings, lock_context, validation_error)
                raise OSError("downloaded archive member gzip failed integrity check: {}".format(validation_error))
    except Exception:
        try:
            tmp.unlink()
        except FileNotFoundError:
            pass
        except OSError:
            pass
        if archive_tmp is not None:
            try:
                archive_tmp.unlink()
            except FileNotFoundError:
                pass
            except OSError:
                pass
        raise
    finally:
        release_download_lock(lock_path, heartbeat_state)
    return True
