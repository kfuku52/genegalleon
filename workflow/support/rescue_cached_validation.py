"""Ordered, bounded CPU parallelism for pristine cached predictor validation.

Every record still uses the same genomic, terminal and boundary checks. Only
read-only indexed genome handles and local validation memos live in workers.
"""
import os
import time
from collections import deque
from concurrent.futures import ProcessPoolExecutor
from itertools import chain, islice
from multiprocessing import get_context


def _r():
    try:
        import rescue_gene_models
    except ImportError:
        from . import rescue_gene_models
    return rescue_gene_models


def checked_prediction(model, memo, genome, code, params, protein, *, engine=None):
    r = engine or _r()
    raw = {key: model[key] for key in r.RAW_FIELDS if key in model}
    raw['cds'] = [list(block) for block in model['cds']]
    raw['evidence'] = dict(model['evidence'])
    checked = memo.validate(model)
    checked = r.complete_terminals(checked, genome, code, params, protein, r.validate_model)
    checked = r.revalidate_boundaries(checked, genome, code, params, protein, r.validate_model)
    checked['problems'] = sorted(set(checked['problems']) | set(
        model.get('annevo_problems', []) + model.get('oriongeno_problems', [])))
    checked['raw_prediction'] = raw
    return r.check_interval(checked)


_CONTEXT = None


def _initialise(genome_path, code, params):
    global _CONTEXT
    import pysam
    r = _r()
    genome = pysam.FastaFile(genome_path)
    memo = r.RawValidationMemo(r.validate_model, genome, code, params,
                               candidate_homology=r.homology_problems)
    _CONTEXT = genome, code, params, memo


def _batch(rows):
    genome, code, params, memo = _CONTEXT
    result = [checked_prediction(model, memo, genome, code, params, protein) for model, protein in rows]
    return result, os.getpid(), memo.diagnostics()


def validate_cached(models, genome, code, params, proteins, cpus, serial, audit):
    """Yield in source order with at most two small pending batches per CPU.

    Small jobs retain the ordinary in-process path. Spawned workers open the
    parent's already staged/indexed genome read-only; no live Fasta/SQLite handle
    crosses the process boundary. Candidate donor proteins travel with the batch.
    """
    if type(cpus) is not int or cpus < 1:
        raise ValueError('Cached prediction validation CPUs must be positive')
    started = time.perf_counter()
    models = iter(models)
    initial = list(islice(models, 512)) if cpus > 1 else []
    audit.update(schema=1, records=0, workers=0, batch_size=64, max_pending_batches=0,
                 decisions='Unchanged genomic, homology, terminal, boundary and interval checks; original record order')
    if cpus == 1 or len(initial) < 512:
        audit['workers'] = 1
        for model in chain(initial, models):
            yield serial(model)
            audit['records'] += 1
        audit['wall_seconds'] = time.perf_counter() - started
        return
    workers = min(cpus, (len(initial) + 63) // 64)
    pending, snapshots = deque(), {}
    audit['workers'] = workers
    with ProcessPoolExecutor(max_workers=workers, mp_context=get_context('spawn'),
            initializer=_initialise, initargs=(os.fsdecode(genome.filename), code, params)) as pool:
        def consume(item):
            expected, future = item
            results, pid, snapshot = future.result()
            if len(results) != expected:
                raise ValueError('Cached prediction validation batch changed record count')
            snapshots[pid] = snapshot
            return results
        remaining = chain(initial, models)
        while batch := list(islice(remaining, 64)):
            rows = [(model, proteins[model['evidence']['donor']][model['evidence']['query']]) for model in batch]
            pending.append((len(rows), pool.submit(_batch, rows)))
            audit['max_pending_batches'] = max(audit['max_pending_batches'], len(pending))
            if len(pending) >= workers * 2:
                for result in consume(pending.popleft()):
                    yield result
                    audit['records'] += 1
        while pending:
            for result in consume(pending.popleft()):
                yield result
                audit['records'] += 1
    fields = ('calls', 'hits', 'misses', 'bypasses', 'evictions', 'oversized', 'exceptions')
    audit['worker_memo_counts'] = {field: sum(value[field] for value in snapshots.values()) for field in fields}
    audit['worker_memo_peak_accounted_bytes_total'] = sum(value['peak_accounted_bytes'] for value in snapshots.values())
    audit['wall_seconds'] = time.perf_counter() - started
