"""Download runtime implementation: locking."""

import hashlib
import os
import re
import shutil
import subprocess
import tarfile
import time
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
from format_species_network import guarded_urlopen as urlopen
from shared_lock import acquire_lock, release_lock

from .local import (
    gzip_integrity_error,
    quarantine_corrupt_gzip,
    quarantine_existing_file,
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


def sleep_before_download_retry(attempt, base_seconds):
    delay = float(base_seconds) * float(attempt)
    if delay > 0:
        time.sleep(delay)


def partial_download_paths(target):
    partial_path = Path(str(target) + ".part")
    url_hash_path = Path(str(partial_path) + ".urlsha256")
    return partial_path, url_hash_path


def discard_partial_download(partial_path, url_hash_path):
    for path in (partial_path, url_hash_path):
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
        return None
    if value < 0:
        return None
    return value


def download_url_to_partial(url, partial_path, headers, timeout, warnings, lock_context):
    resume_from = partial_path.stat().st_size if partial_path.exists() else 0
    request_headers = dict(headers or {})
    if resume_from > 0:
        request_headers["Range"] = "bytes={}-".format(resume_from)
    request = Request(url, headers=request_headers)
    try:
        response_context = urlopen(request, timeout=timeout)
    except HTTPError as exc:
        if int(getattr(exc, "code", 0)) == 416 and resume_from > 0:
            total = parse_unsatisfied_content_range_total(exc.headers.get("Content-Range"))
            if total is not None and total == resume_from:
                return
        raise

    with response_context as response:
        status = response_status_code(response)
        content_range_start, content_range_total = parse_content_range(response.headers.get("Content-Range"))
        append_response = resume_from > 0 and status == 206 and content_range_start == resume_from
        if append_response:
            warnings.append("{} resuming partial download at byte {}".format(lock_context, resume_from))
            file_mode = "ab"
        else:
            if resume_from > 0:
                warnings.append(
                    "{} server did not honor Range bytes={}-; restarting the download".format(
                        lock_context,
                        resume_from,
                    )
                )
            resume_from = 0
            content_range_total = None
            file_mode = "wb"

        expected_response_bytes = response_content_length(response)
        response_bytes = 0
        with open(partial_path, file_mode) as out:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                out.write(chunk)
                response_bytes += len(chunk)

        if expected_response_bytes is not None and response_bytes != expected_response_bytes:
            raise IncompleteRead(b"", max(0, expected_response_bytes - response_bytes))
        final_size = partial_path.stat().st_size
        if content_range_total is not None and final_size != content_range_total:
            raise IncompleteRead(b"", max(0, content_range_total - final_size))


def finalize_partial_download(partial_path, url_hash_path, destination):
    partial_path.replace(destination)
    try:
        url_hash_path.unlink()
    except FileNotFoundError:
        pass


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
                    partial_path.replace(destination)
                    validation_error = gzip_integrity_error(destination)
                    if validation_error is None:
                        last_validation_error = None
                        discard_partial_download(partial_path, partial_url_hash_path)
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
                        sleep_before_download_retry(attempt, retry_base_seconds)
                        continue
                    if not transient_error:
                        discard_partial_download(partial_path, partial_url_hash_path)
                    raise
            if last_validation_error is not None:
                raise OSError("downloaded gzip failed integrity check: {}".format(last_validation_error))
        else:
            if overwrite or not archive_cache_path.exists() or archive_cache_path.stat().st_size == 0:
                archive_partial_path, archive_url_hash_path = prepare_partial_download(
                    archive_cache_path,
                    url,
                    overwrite,
                )
                for attempt in range(1, attempts + 1):
                    try:
                        download_url_to_partial(
                            url,
                            archive_partial_path,
                            request_headers,
                            timeout,
                            warnings,
                            lock_context,
                        )
                        finalize_partial_download(
                            archive_partial_path,
                            archive_url_hash_path,
                            archive_cache_path,
                        )
                        break
                    except Exception as exc:
                        transient_error = is_transient_network_error(exc)
                        if attempt < attempts and transient_error:
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
                        if not transient_error:
                            discard_partial_download(archive_partial_path, archive_url_hash_path)
                        raise
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
            os.replace(tmp, destination)
            _fsync_directory(Path(destination).parent)
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
        raise
    finally:
        release_download_lock(lock_path, heartbeat_state)
    return True
