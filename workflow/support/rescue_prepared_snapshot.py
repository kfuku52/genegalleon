"""Fresh, command-local proofs of immutable rescue preparation publications.

Only completed prepared inputs are retained. Mutable synteny/cache outputs keep
their ordinary full verification. A changed generation is an error, never a
reason to refresh this proof. No stat-authorized checksum survives the command.
"""
import copy
import hashlib
import json
import os
import stat
from pathlib import Path

try:
    from input_generation_array_state import FreshDigestBatch
except ImportError:
    from .input_generation_array_state import FreshDigestBatch


class _Proof:
    """Full bytes first, with pathname, ancestor and permission fences."""

    def __init__(self):
        self.batch = FreshDigestBatch()
        self.components = {}

    @staticmethod
    def identity(path):
        info = path.lstat()
        common = info.st_dev, info.st_ino, info.st_mode
        if stat.S_ISLNK(info.st_mode):
            return common + (info.st_size, info.st_mtime_ns, info.st_ctime_ns, os.readlink(path))
        return common

    def check(self):
        try:
            self.batch.check()
            for path, identity in self.components.items():
                if self.identity(path) != identity:
                    raise OSError('Frozen preparation pathname changed: ' + str(path))
        except OSError as exc:
            raise PreparationChanged(str(exc)) from exc

    def read(self, paths):
        self.check()
        paths = [Path(path).absolute() for path in paths]
        for path in paths:
            resolved = path.resolve(strict=True)
            for component in (*path.parents, path, *resolved.parents, resolved):
                if component not in self.components:
                    self.components[component] = self.identity(component)
        try:
            result = self.batch.read(paths)
        except OSError as exc:
            raise PreparationChanged(str(exc)) from exc
        self.check()
        return result

    def json(self, path):
        path = Path(path).absolute()
        expected = self.read([path])[str(path)]
        try:
            raw = path.read_bytes()
        except OSError as exc:
            raise PreparationChanged(str(exc)) from exc
        self.check()
        if hashlib.sha256(raw).hexdigest() != expected:
            raise PreparationChanged('Frozen preparation JSON changed: ' + str(path))
        return json.loads(raw), expected


class PreparationChanged(OSError):
    """A change after observing an identity must not become an initial repair."""


class PreparedSnapshot:
    """One CLI's independently verified plan and per-species preparation proofs.

    Pair checks touch only their two species, rather than the growing collection
    of previously used species. check_all is the publication/end fence. Missing
    or initially unverified stages use ordinary recovery before their output
    proof is established. Changes to an established generation fail closed.
    """

    def __init__(self, root, plan, producer):
        self.root = Path(root).absolute()
        self.plan = plan
        self.producer = producer
        self.frozen = copy.deepcopy(plan)
        self.plan_proof = _Proof()
        self.species = {}
        self.sources = {}
        self.closed = False
        parsed, self.plan_hash = self.plan_proof.json(self.root / 'plan.json')
        if parsed != self.frozen:
            raise ValueError('Frozen rescue plan changed during execution')
        source_paths = {source[key] for source in self.frozen['request']['sources'].values()
                        for key in ('fasta', 'gff', 'genome')}
        other_files = {path: expected for path, expected in self.frozen['request']['files'].items()
                       if path not in source_paths}
        self._verify_files(self.plan_proof, other_files)

    @staticmethod
    def _verify_files(proof, expected):
        hashes = proof.read(expected)
        if any(hashes[str(Path(path).absolute())] != value for path, value in expected.items()):
            raise ValueError('Frozen rescue input changed')
        return hashes

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
        self.species.clear()
        self.sources.clear()
        self.plan_proof = None
        self.frozen = None
        self.closed = True

    def _open(self):
        if self.closed:
            raise ValueError('Preparation proof is closed')

    def assert_context(self, root, plan=None):
        self._open()
        if Path(root).absolute() != self.root or plan is not None and plan is not self.plan:
            raise ValueError('Preparation proof belongs to a different plan')

    def check(self, names=()):
        self._open()
        self.plan_proof.check()
        # Only global comparison semantics and the requested species are checked
        # here. The complete in-memory plan is checked at publication/end.
        for key in ('parameters', 'tools'):
            if self.plan['request'].get(key) != self.frozen['request'].get(key):
                raise ValueError('Frozen rescue plan changed during execution')
        for name in set(names):
            source = self.frozen['request']['sources'][name]
            if (self.plan['request']['sources'][name] != source
                    or any(self.plan['request']['files'].get(source[key])
                           != self.frozen['request']['files'][source[key]]
                           for key in ('fasta', 'gff', 'genome'))):
                raise ValueError('Frozen rescue plan changed during execution')
            if name in self.species:
                entry = self.species[name]
                entry['proof'].check()
                info = entry['directory'].stat()
                if (info.st_dev, info.st_ino, info.st_mtime_ns, info.st_ctime_ns) != entry['directory_identity']:
                    raise OSError('Frozen prepared directory changed: ' + name)
            if name in self.sources:
                self.sources[name].check()

    def check_all(self):
        self.check(self.species.keys() | self.sources.keys())
        if self.plan != self.frozen:
            raise ValueError('Frozen rescue plan changed during execution')

    def plan_digest(self):
        self.check()
        return self.plan_hash

    def check_job(self, job):
        self.check((job['a'], job['b']))
        index = job.get('index')
        if (type(index) is not int or not 1 <= index <= len(self.frozen['synteny_jobs'])
                or job != self.frozen['synteny_jobs'][index - 1]):
            raise ValueError('Frozen rescue comparison changed during execution')

    def verify_source(self, name):
        self.check([name])
        if name not in self.sources:
            source = self.frozen['request']['sources'][name]
            proof = _Proof()
            self._verify_files(proof, {source[key]: self.frozen['request']['files'][source[key]]
                                       for key in ('fasta', 'gff', 'genome')})
            self.sources[name] = proof
        self.check([name])

    def prepared(self, name):
        self.check([name])
        if name in self.species:
            return self.species[name]['directory']
        self.verify_source(name)
        directory = self.root / 'prepared' / name
        journal = self.root / '.locks' / ('prepared__' + name + '.publish.json')
        entry = None
        if not journal.exists():
            try:
                entry = self._completed(directory, name)
            except PreparationChanged:
                raise
            except (OSError, ValueError):
                # No output proof has been retained yet. Preserve the original
                # first-use restart/recovery policy and its producer lock.
                pass
        if entry is None:
            self.producer.prepared(self.root, self.plan, name, source_snapshot=self)
            self.check([name])
            entry = self._completed(directory, name)
        self.species[name] = entry
        self.check([name])
        return directory

    def _completed(self, directory, name):
        proof = _Proof()
        info = directory.stat()
        identity = info.st_dev, info.st_ino, info.st_mtime_ns, info.st_ctime_ns
        receipt, receipt_hash = proof.json(directory / 'receipt.json')
        key = {'plan': self.plan_hash, 'species': name}
        if (not isinstance(receipt, dict) or receipt.get('key') != key
                or not isinstance(receipt.get('files'), dict) or not receipt['files']
                or not all(isinstance(path, str) and isinstance(value, str)
                           for path, value in receipt['files'].items())):
            raise ValueError('Prepared annotation incomplete or corrupted: ' + name)
        for relative in receipt['files']:
            path = Path(relative)
            if path.is_absolute() or '..' in path.parts or path == Path('receipt.json'):
                raise ValueError('Unsafe prepared receipt path: ' + relative)
        hashes = self._verify_files(proof, {directory / path: value for path, value in receipt['files'].items()})
        return {'directory': directory, 'directory_identity': identity,
                'proof': proof, 'receipt_hash': receipt_hash, 'hashes': hashes}

    def receipt_digest(self, name):
        self.prepared(name)
        return self.species[name]['receipt_hash']

    def file_digest(self, name, filename):
        directory = self.prepared(name)
        path = str(directory / filename)
        if path not in self.species[name]['hashes']:
            raise ValueError('Prepared file absent from receipt: ' + filename)
        return self.species[name]['hashes'][path]
