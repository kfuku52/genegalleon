"""Fair download dispatch with database admission and bounded concurrency."""

import os
import time
from collections import defaultdict, deque
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from urllib.parse import urlparse

from format_species_network import request_database
from input_download_limiter import Admission

from .targets import resolve_provider_download_limits


def dispatch_download_jobs(download_jobs, max_workers, headers, timeout, overwrite, lock_stale_seconds, execute):
    """Fairly dispatch database queues without parking the pool behind one DB.

    Readiness is advisory: the transport still acquires the shared permit at
    every actual request/redirect. Concurrent dispatchers may race for a slot.
    """
    queues = defaultdict(deque)
    for job in download_jobs:
        key = request_database(job["url"], job["provider"]) if urlparse(job["url"]).scheme in ("http", "https", "ftp") else "local"
        queues[key].append(job)
    keys = deque(queues)
    active = defaultdict(int)
    futures = {}
    provider_limits = resolve_provider_download_limits(max_workers)
    permits = {key: Admission(queue[0]["url"], database=key) for key, queue in queues.items() if key != "local"}
    limits = {key: (max_workers if key == "local" else min(
        provider_limits.get("direct" if key.startswith("host:") else key, 1),
        permits[key].limit if permits[key].directory is not None else max_workers)) for key in keys}
    wait_timeout = float(os.environ.get("GG_INPUT_DOWNLOAD_LIMIT_WAIT", "3600"))
    last_progress = time.monotonic()
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        while keys or futures:
            for _ in range(len(keys)):
                key = keys.popleft()
                queue = queues[key]
                if len(futures) < max_workers and active[key] < limits[key] and (key == "local" or permits[key].ready()):
                    job = queue.popleft()
                    future = pool.submit(execute, job, headers, timeout, overwrite, lock_stale_seconds)
                    futures[future] = key
                    active[key] += 1
                    last_progress = time.monotonic()
                if queue:
                    keys.append(key)
            if futures:
                done, _ = wait(futures, timeout=0.05, return_when=FIRST_COMPLETED)
                for future in done:
                    key = futures.pop(future)
                    active[key] -= 1
                    last_progress = time.monotonic()
                    try:
                        yield future.result()
                    except Exception as exc:
                        yield {"errors": ["Unhandled download worker error: {}".format(exc)]}
            elif keys:
                if time.monotonic() - last_progress >= wait_timeout:
                    raise TimeoutError("Timed out waiting for database download queues: " + ", ".join(keys))
                time.sleep(0.05)
