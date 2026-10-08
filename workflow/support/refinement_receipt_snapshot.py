"""Command-local, full-byte-first proofs of immutable stage dependencies.

The ordinary receipt verifier establishes each completed generation. Later uses
fence its file and pathname identities; they never repair or absorb a change.
No proof survives a command, and no global checksum/verifier is overridden.
"""
import hashlib
import json
import os
import stat
from pathlib import Path

try:
    from input_generation_array_state import digest
except ImportError:
    from .input_generation_array_state import digest


class ReceiptChanged(ValueError):
    """An observed dependency generation changed, including during verification."""


def _identity(path, *, member=False):
    info = path.lstat()
    common = info.st_dev, info.st_ino, info.st_mode
    if member or stat.S_ISLNK(info.st_mode):
        common += info.st_size, info.st_mtime_ns, info.st_ctime_ns
    if stat.S_ISLNK(info.st_mode):
        common += (os.readlink(path),)
    return common


class _ReceiptProof:
    def __init__(self, directory):
        self.directory = directory
        self.components = {}
        self.targets = {}

    def observe(self, paths):
        """Capture all identities before the first content verification."""
        self.check()
        try:
            roots = self.directory, self.directory.resolve(strict=True)
            for path in paths:
                resolved = path.resolve(strict=True)
                self.targets[path] = resolved
                for component in (*path.parents, path, *resolved.parents, resolved):
                    # Ancestors outside this publication can gain unrelated
                    # siblings. Its own directories cannot change membership.
                    member = component in (path, resolved) or any(component.is_relative_to(root) for root in roots)
                    current = _identity(component, member=member)
                    previous = self.components.get(component)
                    if previous is not None and previous != (current, member):
                        raise ReceiptChanged('Stage dependency pathname changed: ' + str(component))
                    self.components[component] = current, member
        except OSError as exc:
            raise ReceiptChanged('Stage dependency unavailable: ' + str(self.directory)) from exc
        self.check()

    def check(self):
        try:
            for path, (expected, member) in self.components.items():
                if _identity(path, member=member) != expected:
                    raise ReceiptChanged('Stage dependency content changed: ' + str(path))
            for path, target in self.targets.items():
                if path.resolve(strict=True) != target:
                    raise ReceiptChanged('Stage dependency resolved target changed: ' + str(path))
        except OSError as exc:
            raise ReceiptChanged('Stage dependency unavailable') from exc

    def receipt(self, path, expected):
        self.observe([path])
        try:
            actual = digest(path)
        except (OSError, ValueError) as exc:
            raise ReceiptChanged('Stage dependency receipt verification failed: ' + str(path)) from exc
        self.check()
        try:
            raw = path.read_bytes()
        except OSError as exc:
            raise ReceiptChanged('Stage dependency receipt unavailable: ' + str(path)) from exc
        self.check()
        # Fence before parsing, so a mid-verification malformed replacement is
        # an observed change and cannot be mistaken for an initial repair case.
        if hashlib.sha256(raw).hexdigest() != actual or actual != expected:
            raise ReceiptChanged('Stage dependency receipt changed: ' + str(path))
        value = json.loads(raw)
        if (not isinstance(value, dict) or 'key' not in value
                or not isinstance(value.get('files'), dict) or not value['files']
                or not all(isinstance(name, str) and isinstance(sha, str)
                           for name, sha in value['files'].items())):
            raise ValueError('Invalid stage dependency receipt: ' + str(path))
        for name in value['files']:
            relative = Path(name)
            if (not name or relative.is_absolute() or '..' in relative.parts
                    or relative in (Path('.'), Path('receipt.json'))):
                raise ValueError('Unsafe stage dependency receipt member: ' + name)
        self.observe([self.directory / name for name in value['files']])
        return value


class ReceiptSnapshot:
    """Explicit command scope; each dependency is checked independently.

    Failed first verification retains no successful proof. The caller's ordinary
    stage owns restart/repair of its own output. Dependencies still must already
    be complete, as in the original dependency guard.
    """

    def __init__(self, verifier):
        self.verifier = verifier
        self.entries = {}
        self.closed = False
        self.failure = None

    def _open(self):
        if self.closed:
            raise ValueError('Receipt proof is closed')
        if self.failure is not None:
            raise ReceiptChanged('Stage dependency proof previously changed') from self.failure

    def __enter__(self):
        self._open()
        return self

    def __exit__(self, exc_type, _exc, _traceback):
        try:
            if exc_type is None:
                self.check_all()
        finally:
            self.close()

    def close(self):
        self.entries.clear()
        self.verifier = None
        self.closed = True

    def verify(self, receipts):
        self._open()
        try:
            self._verify(receipts)
        except ReceiptChanged as exc:
            self.failure = exc
            raise

    def _verify(self, receipts):
        for filename, expected in receipts.items():
            path = Path(filename).absolute()
            if path in self.entries:
                self.check({path: expected})
                continue
            proof = _ReceiptProof(path.parent)
            receipt = proof.receipt(path, expected)
            valid = self.verifier(path.parent, receipt['key'])
            proof.check()
            if not valid:
                raise ValueError('Stage dependency content changed: ' + str(path.parent))
            self.entries[path] = expected, proof, dict(receipt['files'])

    def read_json_member(self, filename, expected, member):
        """Read a receipted small projection, with full first proof and fences.

        Return presence separately from JSON null. An unreceipted sidecar is
        never used as evidence; absence permits the caller's legacy path.
        """
        self.verify({filename: expected})
        path = Path(filename).absolute()
        _previous, proof, files = self.entries[path]
        if member not in files:
            return False, None
        try:
            raw = (path.parent / member).read_bytes()
            proof.check()
            if hashlib.sha256(raw).hexdigest() != files[member]:
                raise ReceiptChanged('Stage dependency member changed: ' + member)
        except OSError as exc:
            self.failure = exc
            raise ReceiptChanged('Stage dependency member unavailable: ' + member) from exc
        except ReceiptChanged as exc:
            self.failure = exc
            raise
        return True, json.loads(raw)

    def check(self, receipts):
        self._open()
        try:
            self._check(receipts)
        except ReceiptChanged as exc:
            self.failure = exc
            raise

    def _check(self, receipts):
        for filename, expected in receipts.items():
            path = Path(filename).absolute()
            if path not in self.entries:
                raise ValueError('Stage dependency has no completed proof: ' + str(path))
            previous, proof, _files = self.entries[path]
            if expected != previous:
                raise ReceiptChanged('Stage dependency receipt generation changed: ' + str(path))
            proof.check()

    def check_all(self):
        self._open()
        for _expected, proof, _files in self.entries.values():
            proof.check()
