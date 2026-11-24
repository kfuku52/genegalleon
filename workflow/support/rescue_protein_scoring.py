"""Ordered, restartable protein scoring for a frozen genome-wide rescue search."""
import hashlib
import json
import mmap
import sqlite3
import sys
import zlib
from array import array
from collections import deque
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from contextlib import ExitStack, closing, contextmanager
from functools import partial
from itertools import islice
from multiprocessing import get_context
from struct import Struct

from Bio.Align import _pairwisealigner
from cdskit import protein_similarity

if __package__:
    from . import gene_model_selection as selection
    from . import rescue_alignment_kernel as alignment
else:
    import gene_model_selection as selection
    import rescue_alignment_kernel as alignment


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
def _scorer(paths, layouts, kernel):
    with ExitStack() as resources:
        references = []
        for path, (offsets, sizes, _) in zip(paths, layouts, strict=True):
            with path.open('rb') as stream:
                data = resources.enter_context(mmap.mmap(stream.fileno(), 0, access=mmap.ACCESS_READ))
            references.append((data, offsets, sizes))
        yield references, kernel


def _initialise(paths, layouts):
    global _CONTEXT
    # Keep the mapping alive until the worker exits.
    manager = _scorer(paths, layouts, None)
    _CONTEXT = manager, manager.__enter__()


def _protein(reference, index):
    data, offsets, sizes = reference
    start = offsets[index]
    return data[start:start + sizes[index]].decode().upper().removesuffix('*')


def _score_batch(ids, context=None):
    if context is None:
        context = _CONTEXT[1]
    (queries, subjects), kernel = context
    aligner = selection._aligner()
    score_pair = partial(selection._score_signatures, max_alignment_cells=_MAX_ALIGNMENT_CELLS,
                         aligner_factory=lambda: aligner)
    result = bytearray(len(ids) // _PAIR_IDS.size * _SCORE.size)
    last_query = None
    stats = alignment.buffer() if kernel is not None else None
    for row, (q, s) in enumerate(_PAIR_IDS.iter_unpack(ids)):
        # Most adjacent pairs use the same native protein.
        if last_query != q:
            query = _protein(queries, q)
            encoded_query = query.encode('ascii') if kernel is not None and query.isascii() else None
            last_query = q
        donor = _protein(subjects, s)
        n, m = len(query), len(donor)
        if (encoded_query is not None and donor.isascii() and query != donor
                and min(n, m) >= alignment.MIN_FAST_LENGTH and max(n, m) <= alignment.MAX_FAST_LENGTH
                and n * m <= _MAX_ALIGNMENT_CELLS):
            status = kernel(encoded_query, n, donor.encode('ascii'), m, stats)
            if status != 0:
                raise ValueError(f'Native protein alignment failed for lengths {n} and {m} (status {status})')
            score = alignment.score_counts(n, m, *stats)
        else:
            score = score_pair((query, (), True), (donor, (), True))
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
            yield start, end, cells
            start, cells = end, 0
        cells += cost
    yield start, len(ids), cells


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
    paths = [__file__, selection.__file__, protein_similarity.__file__, _pairwisealigner.__file__, alignment.__file__]
    if alignment.compiled_kernel() is not None:
        paths.append(alignment._LIBRARY)
    for path in paths:
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


def scored_pairs(pairs, queries, subjects, directory, cpus, *, checkpoint=None, checkpoint_key=None, worker_budget=None):
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
    pool = context = None
    kernel = alignment.compiled_kernel()
    if worker_budget is not None and kernel is None:
        raise ValueError('Shared native worker budget requires an installed SIMD kernel')
    with ExitStack() as resources:
        db = (resources.enter_context(_open_checkpoint(checkpoint, [layout[2] for layout in layouts], checkpoint_key))
              if checkpoint is not None else None)

        def start(group, position):
            nonlocal pool, context
            ids = bytearray(len(group) * _PAIR_IDS.size)
            for i, (q, s) in enumerate(group):
                _PAIR_IDS.pack_into(ids, i * _PAIR_IDS.size, qids[q], sids[s])
            pairs_sha = hashlib.sha256(ids).digest() if db is not None else None
            saved = _cached(db, position, pairs_sha, len(group))
            if saved is not None:
                return saved, None, pairs_sha
            if cpus > 1 and len(group) >= 512 and pool is None:
                if kernel is not None:
                    # ctypes releases the GIL in the SIMD kernel. Threads share
                    # the FASTA mappings and offsets and need no score IPC.
                    context = resources.enter_context(_scorer(paths, layouts, kernel))
                    pool = resources.enter_context(ThreadPoolExecutor(max_workers=cpus))
                else:
                    pool = resources.enter_context(ProcessPoolExecutor(
                        max_workers=cpus, mp_context=get_context('spawn'), initializer=_initialise,
                        initargs=(paths, layouts)))
            if pool is None:
                if context is None:
                    context = resources.enter_context(_scorer(paths, layouts, kernel))
                return None, partial(_score_batch, ids, context), pairs_sha
            tasks = sorted(_tasks(ids, layouts[0][1], layouts[1][1]), key=lambda task: -task[2])
            worker_context = context if kernel is not None else None
            if worker_budget is None:
                pending = [(a, b, pool.submit(_score_batch, ids[a:b], worker_context)) for a, b, _ in tasks]
            else:
                pending = [(a, b, pool.submit(worker_budget.score, _score_batch, ids[a:b], worker_context))
                           for a, b, _ in tasks]

            def collect():
                output = bytearray(len(group) * _SCORE.size)
                for a, b, future in pending:
                    start_byte = a // _PAIR_IDS.size * _SCORE.size
                    end_byte = b // _PAIR_IDS.size * _SCORE.size
                    result = future.result()
                    if len(result) != end_byte - start_byte:
                        raise ValueError('Genomewide scoring task changed the pair count')
                    output[start_byte:end_byte] = result
                return output

            return None, collect, pairs_sha

        ready = deque()
        position = 0

        def submit(group):
            nonlocal position
            ready.append((position, group, *start(group, position)))
            position += len(group)

        submit(group)
        if next_group := list(islice(pairs, CHECKPOINT_PAIRS)):
            submit(next_group)
        while ready:
            position_in_cache, group, scores, collect, pairs_sha = ready.popleft()
            if collect is not None:
                scores = collect()
                if len(scores) != len(group) * _SCORE.size:
                    raise ValueError('Genomewide scorer changed the pair count')
                if db is not None:
                    _save(db, position_in_cache, pairs_sha, scores)
            # Keep two groups in flight while the caller performs genomic QC.
            if next_group := list(islice(pairs, CHECKPOINT_PAIRS)):
                submit(next_group)
            for (q, s), fields in zip(group, _SCORE.iter_unpack(scores), strict=True):
                yield q, s, selection.PairScore(*fields)
