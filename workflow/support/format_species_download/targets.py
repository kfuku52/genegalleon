"""Download runtime implementation: targets."""

import gzip
import io
import os
import shutil
import time
import zipfile
from pathlib import Path
from urllib.parse import quote, urlparse
from urllib.request import Request

from format_species_common import (
    parse_positive_int,
)
from format_species_constants import (
    DEFAULT_GLOBAL_DOWNLOAD_WORKERS,
    NCBI_DATASETS_INCLUDE_BY_LABEL,
    PROVIDER_DEFAULT_MAX_CONCURRENT_DOWNLOADS,
)
from format_species_network import guarded_urlopen as urlopen
from format_species_provider_config import (
    ENSEMBL_LIKE_PROVIDERS,
    PROVIDERS,
)
from format_species_provider_urls import (
    append_cookie_header,
    fetch_jgi_session_cookie,
    resolve_ncbi_datasets_base_url,
)

from .locking import (
    acquire_download_lock,
    release_download_lock,
)


def default_download_filename(provider, species_key, label, url, archive_member=""):
    archive_member_text = str(archive_member or "").strip()
    if archive_member_text != "":
        base = Path(archive_member_text).name
    else:
        parsed = urlparse(url)
        base = Path(parsed.path).name
    if base == "":
        if label in ("cds", "genome"):
            ext = ".fa.gz"
        elif label == "gbff":
            ext = ".gbff.gz"
        else:
            ext = ".gff3.gz"
        base = "{}.{}{}".format(species_key, label, ext)
    if provider in ENSEMBL_LIKE_PROVIDERS:
        if not base.startswith(species_key + "."):
            base = "{}.{}".format(species_key, base)
    return base


def parse_http_headers(
    http_header_values,
    auth_bearer_token_env,
    auth_cookie_env="",
    jgi_login_env="",
    jgi_password_env="",
    timeout=120.0,
):
    headers = {}
    for raw in http_header_values:
        text = raw.strip()
        if text == "":
            continue
        if ":" not in text:
            raise ValueError("Invalid --http-header value (missing ':'): {}".format(raw))
        key, value = text.split(":", 1)
        key = key.strip()
        value = value.strip()
        if key == "":
            raise ValueError("Invalid --http-header value (empty key): {}".format(raw))
        headers[key] = value
    if auth_bearer_token_env != "":
        token = os.environ.get(auth_bearer_token_env, "").strip()
        if token == "":
            raise ValueError(
                "Environment variable for bearer token is empty or undefined: {}".format(auth_bearer_token_env)
            )
        headers["Authorization"] = "Bearer {}".format(token)
    if auth_cookie_env != "":
        cookie_value = os.environ.get(auth_cookie_env, "").strip()
        if cookie_value == "":
            raise ValueError("Environment variable for cookie header is empty or undefined: {}".format(auth_cookie_env))
        headers["Cookie"] = append_cookie_header(headers.get("Cookie", ""), cookie_value)
    if jgi_login_env != "" or jgi_password_env != "":
        if jgi_login_env == "" or jgi_password_env == "":
            raise ValueError("--jgi-login-env and --jgi-password-env must be provided together")
        login_value = os.environ.get(jgi_login_env, "").strip()
        password_value = os.environ.get(jgi_password_env, "").strip()
        if login_value == "":
            raise ValueError("Environment variable for JGI login is empty or undefined: {}".format(jgi_login_env))
        if password_value == "":
            raise ValueError("Environment variable for JGI password is empty or undefined: {}".format(jgi_password_env))
        jgi_cookie = fetch_jgi_session_cookie(
            login=login_value,
            password=password_value,
            timeout=float(timeout),
            headers=headers,
        )
        headers["Cookie"] = append_cookie_header(headers.get("Cookie", ""), jgi_cookie)
    return headers


def resolve_parallel_jobs(requested_jobs):
    jobs = parse_positive_int(requested_jobs, 0)
    if jobs > 0:
        return jobs
    env_slots = os.environ.get("GG_TASK_CPUS", "").strip()
    if env_slots == "":
        env_slots = os.environ.get("NSLOTS", "").strip()
    if env_slots != "":
        return parse_positive_int(env_slots, DEFAULT_GLOBAL_DOWNLOAD_WORKERS)
    return DEFAULT_GLOBAL_DOWNLOAD_WORKERS


def resolve_provider_download_limits(global_jobs):
    limits = {}
    for provider in PROVIDERS:
        default_limit = PROVIDER_DEFAULT_MAX_CONCURRENT_DOWNLOADS.get(provider, 1)
        env_name = "GG_INPUT_MAX_CONCURRENT_DOWNLOADS_{}".format(provider.upper())
        env_value = os.environ.get(env_name, "").strip()
        limit = parse_positive_int(env_value, default_limit)
        limits[provider] = max(1, min(int(global_jobs), int(limit)))
    return limits


def resolve_download_urls_from_templates(provider, source_id, species_key, row):
    cds_template = str(row.get("cds_url_template", "") or "").strip()
    gff_template = str(row.get("gff_url_template", "") or "").strip()
    genome_template = str(row.get("genome_url_template", "") or "").strip()
    if cds_template == "" and gff_template == "" and genome_template == "":
        return None
    mapping = {"id": source_id, "species_key": species_key, "provider": provider}
    resolved = {}
    try:
        if cds_template != "":
            resolved["cds_url"] = cds_template.format(**mapping)
        if gff_template != "":
            resolved["gff_url"] = gff_template.format(**mapping)
        if genome_template != "":
            resolved["genome_url"] = genome_template.format(**mapping)
    except KeyError as exc:
        raise ValueError("template placeholder is not supported: {}".format(exc))
    return resolved


def pick_ncbi_datasets_member_name(member_names, label):
    if label == "CDS":
        suffixes = ("/cds_from_genomic.fna", "_cds_from_genomic.fna")
    elif label == "GFF":
        suffixes = ("/genomic.gff", "_genomic.gff", "/genomic.gff3", "_genomic.gff3")
    elif label == "GBFF":
        suffixes = ("/genomic.gbff", "_genomic.gbff", "/genomic.gbff3", "_genomic.gbff3")
    elif label == "GENOME":
        suffixes = ("/genomic.fna", "_genomic.fna", "/genomic.fa", "_genomic.fa")
    else:
        raise ValueError("Unsupported NCBI datasets label: {}".format(label))

    matches = [name for name in member_names if any(name.endswith(suffix) for suffix in suffixes)]
    if len(matches) == 0:
        return ""
    return sorted(matches)[0]


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


def write_download_stream(destination, source):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(str(destination) + ".tmp.{}".format(os.getpid()))
    try:
        if destination.name.lower().endswith(".gz"):
            with gzip.open(tmp, "wb") as out:
                shutil.copyfileobj(source, out, length=1024 * 1024)
        else:
            with open(tmp, "wb") as out:
                shutil.copyfileobj(source, out, length=1024 * 1024)
        with open(tmp, "rb") as handle:
            os.fsync(handle.fileno())
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


def write_download_payload(destination, payload_bytes):
    with io.BytesIO(payload_bytes) as source:
        write_download_stream(destination, source)


def download_ncbi_datasets_file_from_id(
    source_id,
    label,
    destination,
    headers,
    timeout,
    dry_run,
    overwrite,
    lock_stale_seconds,
    warnings,
    lock_context,
):
    include_annotation_type = NCBI_DATASETS_INCLUDE_BY_LABEL.get(label, "")
    if include_annotation_type == "":
        raise ValueError("No datasets include_annotation_type is defined for label: {}".format(label))
    if dry_run:
        return False

    lock_path = Path(str(destination) + ".lock")
    heartbeat_state = acquire_download_lock(lock_path, lock_stale_seconds, warnings, lock_context)
    tmp_zip = Path(str(destination) + ".datasets.tmp.{}".format(os.getpid()))
    try:
        if destination.exists() and destination.stat().st_size > 0 and not overwrite:
            return False

        base = resolve_ncbi_datasets_base_url()
        datasets_url = "{}/genome/accession/{}/download?include_annotation_type={}".format(
            base,
            quote(source_id, safe=""),
            quote(include_annotation_type, safe=""),
        )

        req_headers = dict(headers)
        if "User-Agent" not in req_headers:
            req_headers["User-Agent"] = "genegalleon-input-generation"
        if "Accept" not in req_headers:
            req_headers["Accept"] = "application/zip"
        if "Accept-Encoding" not in req_headers:
            req_headers["Accept-Encoding"] = "identity"

        last_error = None
        for attempt in range(1, 4):
            try:
                request = Request(datasets_url, headers=req_headers)
                with urlopen(request, timeout=timeout) as response, open(tmp_zip, "wb") as out:
                    while True:
                        chunk = response.read(1024 * 1024)
                        if not chunk:
                            break
                        out.write(chunk)

                with zipfile.ZipFile(tmp_zip) as archive:
                    member_name = pick_ncbi_datasets_member_name(archive.namelist(), label)
                    if member_name == "":
                        raise ValueError(
                            "datasets archive did not contain expected {} member for id {}".format(label, source_id)
                        )
                    with archive.open(member_name, "r") as source:
                        write_download_stream(destination, source)
                last_error = None
                break
            except Exception as exc:
                last_error = exc
                try:
                    tmp_zip.unlink()
                except FileNotFoundError:
                    pass
                except OSError:
                    pass
                if attempt < 3:
                    time.sleep(float(attempt))
                    continue
        if last_error is not None:
            raise last_error
    except Exception:
        try:
            tmp_zip.unlink()
        except FileNotFoundError:
            pass
        except OSError:
            pass
        raise
    finally:
        try:
            tmp_zip.unlink()
        except FileNotFoundError:
            pass
        except OSError:
            pass
        release_download_lock(lock_path, heartbeat_state)
    return True
