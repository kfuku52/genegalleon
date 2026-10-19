"""Finite recovery for interrupted bodies that make validated byte progress."""
import math
import os
import time


def positive_setting(name, default):
    value = float(os.environ.get(name, default))
    if not math.isfinite(value) or value <= 0:
        raise ValueError(name + " must be finite and positive")
    return value


class DownloadRecovery:
    def __init__(self, attempts):
        self.attempts = attempts
        self.failures = 0
        self.started = self.last_progress = time.monotonic()
        self.total_seconds = positive_setting("GG_DOWNLOAD_MAX_TOTAL_SECONDS", 21600)
        self.stall_seconds = positive_setting("GG_DOWNLOAD_STALL_SECONDS", 900)
        self.min_bytes = positive_setting("GG_DOWNLOAD_MIN_PROGRESS_BYTES", 1024 * 1024)
        self.identity = None
        self.anchor = 0

    def check(self):
        now = time.monotonic()
        if now - self.started >= self.total_seconds:
            raise TimeoutError("download recovery exceeded total elapsed-time limit")
        if now - self.last_progress >= self.stall_seconds:
            raise TimeoutError("download recovery exceeded meaningful-progress stall limit")

    def remaining(self):
        self.check()
        now = time.monotonic()
        return min(self.total_seconds - (now - self.started),
                   self.stall_seconds - (now - self.last_progress))

    def response(self, identity, offset):
        validator = identity.get("etag") or identity.get("last_modified")
        signature = (validator, identity.get("total")) if validator else None
        if signature is None or signature != self.identity:
            # Changing representations or repeatedly restarting at zero must
            # not itself replenish a retry or elapsed-time budget.
            self.identity = signature
            self.anchor = offset
        elif offset < self.anchor:
            self.anchor = offset

    def progress(self, offset):
        self.check()
        if offset - self.anchor >= self.min_bytes:
            self.anchor = offset
            self.last_progress = time.monotonic()
            if self.identity is not None:
                self.failures = 0

    def retry(self):
        self.check()
        self.failures += 1
        return self.failures < self.attempts
