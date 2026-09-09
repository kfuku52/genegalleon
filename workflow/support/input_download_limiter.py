"""Cross-process request admission on a shared filesystem with working flock locks.

Never delete lock files while clients may be running. Kernel-held slot locks are
released on process exit; timestamps limit request starts, not lock ownership.
"""
import fcntl
import hashlib
import json
import math
import os
import time
from contextlib import contextmanager
from pathlib import Path
from urllib.parse import urlparse

DATABASE_DOMAINS = {
    "ncbi": ("ncbi.nlm.nih.gov",),
    "ensembl": ("ensembl.org", "ensemblgenomes.org"),
    "coge": ("genomevolution.org",),
    "cngb": ("cngb.org",),
    "gwh": ("big.ac.cn", "ngdc.cncb.ac.cn"),
    "ddbj": ("ddbj.nig.ac.jp",),
}


def database_key(url):
    host = (urlparse(url).hostname or "").lower().rstrip(".")
    for key, domains in DATABASE_DOMAINS.items():
        if any(host == domain or host.endswith("." + domain) for domain in domains):
            return key
    return "host:" + host


@contextmanager
def locked(path, deadline):
    with open(path, "a+b") as handle:
        while True:
            try:
                fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
                break
            except BlockingIOError:
                if time.monotonic() >= deadline:
                    raise TimeoutError("Timed out waiting for download admission metadata lock") from None
                time.sleep(0.05)
        yield


def limit_directory():
    return os.environ.get("GG_INPUT_DOWNLOAD_LIMIT_DIR", "").strip() or os.environ.get("download_limit_dir", "").strip()


class Admission:
    def __init__(self, url, database=None):
        self.slot = None
        root = limit_directory()
        self.directory = None
        if not root or urlparse(url).scheme not in ("http", "https", "ftp"):
            return
        self.key = database or database_key(url)
        suffix = "DIRECT" if self.key.startswith("host:") else self.key.upper()
        self.limit = int(os.environ.get("GG_INPUT_MAX_CONCURRENT_DOWNLOADS_" + suffix, "2"))
        self.interval = float(os.environ.get("GG_INPUT_REQUEST_INTERVAL_" + suffix, "0.4"))
        self.wait_timeout = float(os.environ.get("GG_INPUT_DOWNLOAD_LIMIT_WAIT", "3600"))
        if self.limit < 1 or not math.isfinite(self.interval) or self.interval < 0 or not math.isfinite(self.wait_timeout) or self.wait_timeout <= 0:
            raise ValueError("Invalid shared download limit for " + self.key)
        self.directory = Path(root).expanduser().resolve() / hashlib.sha256(self.key.encode()).hexdigest()
        self.directory.mkdir(parents=True, exist_ok=True)
        self.state_path = self.directory / "state.json"
        self.policy = {"database": self.key, "limit": self.limit, "interval": self.interval}

    def state(self):
        if not self.state_path.exists():
            return {"policy": self.policy, "next_start": 0, "cooldown": 0}
        state = json.loads(self.state_path.read_text())
        if state["policy"] != self.policy:
            raise ValueError("Shared download policy mismatch for " + self.key)
        return state

    def save(self, state):
        temp = self.directory / "state.tmp"
        temp.write_text(json.dumps(state), encoding="utf-8")
        os.replace(temp, self.state_path)

    def acquire(self):
        if self.directory is None:
            return self
        deadline = time.monotonic() + self.wait_timeout
        while time.monotonic() < deadline:
            with locked(self.directory / "admission.lock", deadline):
                state = self.state()
                if not self.state_path.exists():
                    self.save(state)
                now = time.time()
                if now >= max(state["next_start"], state["cooldown"]):
                    for index in range(self.limit):
                        slot = open(self.directory / (str(index) + ".slot"), "a+b")
                        try:
                            fcntl.flock(slot, fcntl.LOCK_EX | fcntl.LOCK_NB)
                        except BlockingIOError:
                            slot.close()
                            continue
                        self.slot = slot
                        state["next_start"] = now + self.interval
                        try:
                            self.save(state)
                        except BaseException:
                            self.close()
                            raise
                        return self
            time.sleep(0.05)
        raise TimeoutError("Timed out waiting for shared download admission: " + self.key)

    def cooldown(self, seconds):
        if self.directory is not None:
            with locked(self.directory / "admission.lock", time.monotonic() + self.wait_timeout):
                state = self.state()
                state["cooldown"] = max(state["cooldown"], time.time() + seconds)
                self.save(state)

    def close(self):
        if self.slot is not None:
            self.slot.close()
            self.slot = None
