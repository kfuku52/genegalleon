"""Full-byte-first input proofs owned by one explicit refinement command.

The preparation primitive already fences file bytes, resolved targets, pathname
ancestors, permissions and symlinks. This wrapper binds the whole plan and keeps
failures sticky; no checksum or stat-authorized proof survives the command.
"""
import copy
import os
from pathlib import Path

try:
    from rescue_prepared_snapshot import _Proof
except ImportError:
    from .rescue_prepared_snapshot import _Proof


class InputChanged(ValueError):
    """An input or frozen plan generation changed during the command."""


class RefinementInputSnapshot:
    def __init__(self):
        self.root = None
        self.plan = None
        self.frozen = None
        self.plan_proof = _Proof()
        self.groups = {}
        self.path_groups = {}
        self.expected = {}
        self.hashes = {}
        self.hash_owners = {}
        self.physical_owners = {}
        self.target_owners = {}
        self.closed = False
        self.failure = None

    def _open(self):
        if self.closed:
            raise ValueError('Refinement input proof is closed')
        if self.failure is not None:
            raise InputChanged('Refinement input proof previously changed') from self.failure

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
        self.groups.clear()
        self.path_groups.clear()
        self.expected.clear()
        self.hashes.clear()
        self.hash_owners.clear()
        self.physical_owners.clear()
        self.target_owners.clear()
        self.root = self.plan = self.frozen = self.plan_proof = None
        self.closed = True

    def _failed(self, exc):
        self.failure = exc
        if isinstance(exc, InputChanged):
            raise exc
        raise InputChanged('Frozen refinement input or plan changed') from exc

    def reject_plan_read(self, exc):
        self._open()
        self._failed(exc)

    def verify_plan(self, root, value):
        """Read small plan bytes afresh and reject semantic or path replacement."""
        self._open()
        try:
            root = Path(root).absolute()
            if self.root is not None and root != self.root:
                raise InputChanged('Refinement input proof belongs to a different plan')
            parsed, _sha = self.plan_proof.json(root / 'plan.json')
            if self.root is None:
                if parsed != value:
                    raise InputChanged('Frozen refinement plan changed during execution')
                self.root, self.plan, self.frozen = root, value, copy.deepcopy(value)
                self.expected = dict(value['request']['files'])
                # Assign shared logical paths once. The caller still passes its
                # exact names/global selection; groups do not broaden it.
                for name, source in value['request']['sources'].items():
                    for field in ('fasta', 'gff', 'genome'):
                        self.path_groups.setdefault(source[field], name)
            if parsed != self.frozen or value != self.frozen or self.plan != self.frozen:
                raise InputChanged('Frozen refinement plan changed during execution')
        except (OSError, ValueError) as exc:
            self._failed(exc)

    def verify_files(self, files):
        """Verify only the mapping selected by ordinary load(), bytes first."""
        self._open()
        try:
            if self.root is None:
                raise ValueError('Refinement input proof has no frozen plan')
            selected = {}
            for path, expected in files.items():
                if path not in self.expected or expected != self.expected[path]:
                    raise InputChanged('Frozen refinement input changed')
                selected.setdefault(self.path_groups.get(path), {})[path] = expected
            for group, entries in selected.items():
                if group not in self.groups:
                    proof = _Proof()
                    # Preserve FreshDigestBatch physical-identity dedup across
                    # scopes, but independently fence each new alias pathname.
                    proof.batch.hashes = self.hashes
                    self.groups[group] = proof
                proof = self.groups[group]
                owners = set()
                for path in entries:
                    info = os.stat(path)
                    identity = info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns
                    owners.update((self.hash_owners.get(identity),
                                   self.physical_owners.get(identity[:2]),
                                   self.target_owners.get(str(Path(path).resolve(strict=True)))))
                for owner in owners - {None}:
                    owner.check()
                actual = proof.read(entries)
                if any(actual[str(Path(path).absolute())] != expected for path, expected in entries.items()):
                    raise InputChanged('Frozen refinement input changed')
                for identity, target in proof.batch.paths.values():
                    self.hash_owners.setdefault(identity, proof)
                    self.physical_owners.setdefault(identity[:2], proof)
                    self.target_owners.setdefault(target, proof)
        except (OSError, ValueError) as exc:
            self._failed(exc)

    def check_all(self):
        self._open()
        if self.root is None:
            return
        try:
            parsed, _sha = self.plan_proof.json(self.root / 'plan.json')
            if parsed != self.frozen or self.plan != self.frozen:
                raise InputChanged('Frozen refinement plan changed during execution')
            for proof in self.groups.values():
                proof.check()
        except (OSError, ValueError) as exc:
            self._failed(exc)
