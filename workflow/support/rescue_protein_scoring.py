"""Ordered, restartable protein scoring for a frozen genome-wide rescue search."""
import hashlib
import json
import mmap
import sqlite3
import sys
import zlib
from array import array
from concurrent.futures import ProcessPoolExecutor
from contextlib import ExitStack, closing, contextmanager
from functools import partial
from itertools import islice
from multiprocessing import get_context
from struct import Struct

from Bio.Align import _pairwisealigner
from cdskit import protein_similarity

if __package__:
    from . import gene_model_selection as selection
else:
    import gene_model_selection as selection


CHECKPOINT_PAIRS = 16_384
_PAIR_IDS = Struct('!II')
_SCORE = Struct('!5dI???')
_MAX_ALIGNMENT_CELLS = 25_000_000
_CONTEXT = None


def _layout(path, proteins):
    """Verify the search FASTA and index it once for all read-only workers."""
    offsets, sizes = array('Q'), array('I')
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for token, protein in proteins.items():
            header = stream.readline()
            offset = stream.tell()
            sequence = stream.readline()
            if header != f'>{token}\n'.encode() or sequence != (protein + '\n').encode():
                raise ValueError(f'Genomewide score FASTA differs at {token} in {path}')
            offsets.append(offset)
            sizes.append(len(sequence) - 1)
            digest.update(header)
            digest.update(sequence)
        if stream.read(1):
            raise ValueError(f'Unexpected additional records in genomewide score FASTA {path}')
    return offsets, sizes, digest.hexdigest()


@contextmanager
def _scorer(paths, layouts):
    with ExitStack() as resources:
        references = []
        for path, (offsets, sizes, _) in zip(paths, layouts, strict=True):
            with path.open('rb') as stream:
                data = resources.enter_context(mmap.mmap(stream.fileno(), 0, access=mmap.ACCESS_READ))
            references.append((data, offsets, sizes))
        aligner = selection._aligner()
        yield references, aligner


def _initialise(paths, layouts):
    global _CONTEXT
    # Keep the mapping alive until the worker exits.
    manager = _scorer(paths, layouts)
    _CONTEXT = manager, manager.__enter__()


def _protein(reference, index):
    data, offsets, sizes = reference
    start = offsets[index]
    return data[start:start + sizes[index]].decode().upper().removesuffix('*')


def _score_batch(ids, context=None):
    if context is None:
        context = _CONTEXT[1]
    (queries, subjects), aligner = context
    score_pair = partial(selection._score_signatures, max_alignment_cells=_MAX_ALIGNMENT_CELLS,
                         aligner_factory=lambda: aligner)
    result = bytearray(len(ids) // _PAIR_IDS.size * _SCORE.size)
    last_query = None
    for row, (q, s) in enumerate(_PAIR_IDS.iter_unpack(ids)):
        # Most adjacent pairs use the same native protein.
        if last_query != q:
            query = (_protein(queries, q), (), True)
            last_query = q
        donor = (_protein(subjects, s), (), True)
        score = score_pair(query, donor)
        _SCORE.pack_into(result, row * _SCORE.size, score.score, score.identity, score.coverage_a,
                         score.coverage_b, score.junction_similarity, score.aligned_residues,
                         score.internal_loss_a, score.internal_loss_b, score.bounded)
    return result


def _tasks(ids, query_sizes, donor_sizes):
    """Keep both the number of pairs and estimated alignment work per task bounded."""
    start = cells = 0
    for i, (q, s) in enumerate(_PAIR_IDS.iter_unpack(ids)):
        cost = query_sizes[q] * donor_sizes[s]
        if cost > _MAX_ALIGNMENT_CELLS:
            # Allow for the optional terminal stop in the FASTA lengths.
            cost = 10_000 if (query_sizes[q] - 1) * (donor_sizes[s] - 1) > _MAX_ALIGNMENT_CELLS else cost
        end = i * _PAIR_IDS.size
        if end > start and (end - start >= 256 * _PAIR_IDS.size or cells + cost > 20_000_000):
            yield ids[start:end]
            start, cells = end, 0
        cells += cost
    yield ids[start:]


@contextmanager
def _open_checkpoint(path, fasta_digests, key):
    path.parent.mkdir(parents=True, exist_ok=True)
    # The scientific implementation and the on-disk format are both part of
    # the checkpoint identity. CPUs and task boundaries within a group are not.
    identity = json.dumps({
        'format': 1, 'key': key, 'fasta': fasta_digests, 'python': sys.version,
        'kernel': [hashlib.file_digest(opened, 'sha256').hexdigest()
                   for opened in _kernel_files()],
    }, sort_keys=True)
    with closing(sqlite3.connect(path, timeout=60)) as db:
        db.execute('PRAGMA journal_mode=DELETE').fetchone()
        with db:
            db.execute('CREATE TABLE IF NOT EXISTS scope (id INTEGER PRIMARY KEY CHECK (id = 1), identity TEXT)')
            db.execute('INSERT OR IGNORE INTO scope VALUES (1, ?)', (identity,))
            if db.execute('SELECT identity FROM scope WHERE id=1').fetchone() != (identity,):
                raise ValueError(f'Genomewide score checkpoint inputs changed: {path}')
            db.execute('CREATE TABLE IF NOT EXISTS scores (position INTEGER PRIMARY KEY, '
                       'pairs_sha BLOB, scores_sha BLOB, payload BLOB)')
        yield db


def _kernel_files():
    for path in (__file__, selection.__file__, protein_similarity.__file__, _pairwisealigner.__file__):
        with open(path, 'rb') as stream:
            yield stream


def _cached(db, position, pairs_sha, pairs_count):
    if db is None:
        return None
    row = db.execute('SELECT pairs_sha, scores_sha, payload FROM scores WHERE position=?', (position,)).fetchone()
    if row is None:
        return None
    old_pairs, scores_sha, payload = row
    if old_pairs != pairs_sha:
        raise ValueError(f'Genomewide score checkpoint pair IDs changed at group {position}')
    if hashlib.sha256(payload).digest() != scores_sha:
        raise ValueError(f'Damaged genomewide score checkpoint at group {position}')
    size = pairs_count * _SCORE.size
    decoder = zlib.decompressobj()
    scores = decoder.decompress(payload, size + 1)
    if len(scores) != size or not decoder.eof or decoder.unused_data:
        raise ValueError(f'Invalid genomewide score checkpoint size at group {position}')
    return scores


def _save(db, position, pairs_sha, scores):
    payload = zlib.compress(scores, level=1)
    with db:
        db.execute('INSERT OR IGNORE INTO scores VALUES (?, ?, ?, ?)',
                   (position, pairs_sha, hashlib.sha256(payload).digest(), payload))


def scored_pairs(pairs, queries, subjects, directory, cpus, *, checkpoint=None, checkpoint_key=None):
    """Read the existing FASTAs, returning exact scores in their original pair order.

    Only compact pair IDs cross worker boundaries. Completed groups commit into
    one compressed checkpoint, and retries can use a different worker count.
    """
    pairs = iter(pairs)
    group = list(islice(pairs, CHECKPOINT_PAIRS))
    if not group:
        return
    paths = directory / 'native.pep', directory / 'donors.pep'
    layouts = tuple(_layout(path, proteins) for path, proteins in zip(paths, (queries, subjects), strict=True))
    qids = {token: i for i, token in enumerate(queries)}
    sids = {token: i for i, token in enumerate(subjects)}
    pool = serial = None
    position = 0
    with ExitStack() as resources:
        db = (resources.enter_context(_open_checkpoint(checkpoint, [layout[2] for layout in layouts], checkpoint_key))
              if checkpoint is not None else None)

        def start(group, position):
            nonlocal pool, serial
            ids = bytearray(len(group) * _PAIR_IDS.size)
            for i, (q, s) in enumerate(group):
                _PAIR_IDS.pack_into(ids, i * _PAIR_IDS.size, qids[q], sids[s])
            pairs_sha = hashlib.sha256(ids).digest() if db is not None else None
            saved = _cached(db, position, pairs_sha, len(group))
            if saved is not None:
                return saved, None, pairs_sha
            if cpus > 1 and len(group) >= 512 and pool is None:
                pool = resources.enter_context(ProcessPoolExecutor(max_workers=cpus, mp_context=get_context('spawn'),
                                                                   initializer=_initialise, initargs=(paths, layouts)))
            if pool is not None:
                tasks = _tasks(ids, layouts[0][1], layouts[1][1])
                return None, pool.map(_score_batch, tasks), pairs_sha
            if serial is None:
                serial = resources.enter_context(_scorer(paths, layouts))
            return None, (_score_batch(ids, serial) for ids in (ids,)), pairs_sha

        work = start(group, position)
        while group:
            scores, pending, pairs_sha = work
            if pending is not None:
                scores = b''.join(pending)
                if len(scores) != len(group) * _SCORE.size:
                    raise ValueError('Genomewide scorer changed the pair count')
                if db is not None:
                    _save(db, position, pairs_sha, scores)
            next_group = list(islice(pairs, CHECKPOINT_PAIRS))
            position += len(group)
            if next_group:
                # Overlap the next block of CPU work with genomic validation.
                work = start(next_group, position)
            for (q, s), fields in zip(group, _SCORE.iter_unpack(scores), strict=True):
                yield q, s, selection.PairScore(*fields)
            group = next_group
