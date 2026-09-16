# Reliable input downloads

Input generation retains one bounded retry budget per downloaded file, including
body interruptions. `GG_DOWNLOAD_ATTEMPTS` defaults to four. Transient failures
use exponential backoff with jitter, capped at 300 seconds before applying a
server's `Retry-After` minimum. Permanent HTTP errors are not repeatedly retried.
The shared database limiter still applies cooldowns and per-database concurrency.
Metadata requests also honor `Retry-After` when the limiter is disabled.

Before reading a body, the downloader checks destination writability, reported
response size against available filesystem space, content type, encoding and
Range consistency. This is a per-response check, not a reservation or a quota
query; concurrent writers and user quotas can still exhaust storage. Metadata is
read from the GET response, avoiding HEAD-only failures or a duplicate request.
Existing manifest resolution checks local source paths before download execution.

A `.part.identity.json` sidecar records the strong ETag (or Last-Modified) and
expected total size. Resume uses `If-Range`; a full response replaces the partial
file, and inconsistent ranges cannot be appended. Legacy partials without a
validator restart safely. Servers without validators use a single full response.
HTTP authorization/unavailability and disk/quota exhaustion preserve identified
partial files for a later run. A missing or changed resume validator causes a
fresh transfer, even when another validator type is present. A 416 response triggers a fresh request rather than declaring success
from the byte count alone.

CNGB transfers use sequential 128 MiB ranges when the server supports validated
ranges. Set `GG_DOWNLOAD_RANGE_CHUNK_BYTES=0` to stream, or a positive byte count
to select a range size for all providers. This bounds connection lifetimes;
it does not imply measured improvements in throughput or failure rate. A server
ignoring Range falls back to its full response. Existing
`GG_INPUT_MAX_CONCURRENT_DOWNLOADS_<PROVIDER>` controls remain available; CNGB's
default is one. No extra parallel chunk workers are created.

Figshare uses an identifiable browser-compatible default User-Agent, while an
explicit header (case-insensitive) takes precedence. This is a compatibility
workaround for observed 403 responses, not a bypass for private data. Remove it
if controlled provider tests establish that the ordinary agent is consistently
accepted. HTML responses are rejected, and gzip is detected by its signature as
well as its extension, including numeric Figshare filenames. Empty files and
compressed or BOM/comment-prefixed HTML are rejected. Version-2 validation
receipts invalidate earlier checks so cached data receives these stronger checks.

NCBI Datasets uses the same resumable transport and retains its downloaded ZIP.
Extraction requires the requested accession/version and an unambiguous member;
ZIP member CRC and gzip integrity are checked before final publication. Failed
overwrite attempts retain the previous complete destination. Generic archive
members share the same validated, resumable archive download path. Local disk failures do not trigger
another network fallback. Existing provider resolvers supply official alternate
locations; the downloader never substitutes a different assembly version to make
a failed bundle succeed.

## Reuse across retry workspaces

Opt in by exporting these variables before starting the input job. The standard
input entrypoint forwards `GG_DOWNLOAD_*` settings into the container:

```bash
export GG_DOWNLOAD_SHARED_CACHE_DIR=/shared/project/input-download-cache
export GG_DOWNLOAD_EVENT_DIR=/shared/project/input-download-events
```

The cache keys URL, archive member and request headers without storing their
plaintext. A completed payload has a SHA-256/size receipt checked on every reuse;
an invalid or older receipt requires a fresh download. The materialized copy
is checked against the receipt before atomic publication. This verifies local cache integrity, not
publisher authenticity. Partial files survive job retries. Publication uses locks
and atomic replacement; destinations are independent copies so downstream edits
cannot alter cached bytes. Allow storage for both cache and materialized outputs.
Only use a private project directory visible and writable inside the container.
Output destinations must be outside the shared-cache directory.
These settings are opt-in and do not change existing running jobs or frozen plans.

The event directory contains one schema-version-1 JSON file per transport outcome:
`database`, `url_sha256`, `status` (`downloaded`, `materialized`, `reused`, `failed`), `http_status`,
`error_class`, `elapsed_seconds`, `retries`, and final `bytes`. URLs and headers are
omitted. `materialized` means a shared-cache copy with no new network download;
`reused` means the destination already matched. Event files appear atomically.
These are transport outcomes, not species completion or workflow success;
manifest cache skips before transport are not included. Compare counts within the
same collection scope, rather than treating these events as all historical jobs.

## Validation

`workflow/tests/test_download_reliability.py` injects loopback HTTP failures for
midstream disconnects, changed validators, invalid ranges, HTML, corrupt
extensionless gzip, retry-budget exhaustion, shared-cache corruption and
accession-scoped extraction. It also checks failed-overwrite preservation,
quota/disk partial preservation, concurrent cache publication, stale receipts,
corrupted archive recovery and credential stripping on redirects with and
without the shared limiter. Run these with the existing download, limiter and
cache-validation tests in the GeneGalleon container, as described in
[agent runtime validation](agent-runtime-validation.md).
