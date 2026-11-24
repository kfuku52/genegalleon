"""Ordered, bounded CPU parallelism for pristine cached predictor validation.

Every record still uses the same genomic, terminal and boundary checks. Only
read-only indexed genome handles and local validation memos live in workers.
"""
import os
import time
from collections import deque
from concurrent.futures import ProcessPoolExecutor
from contextlib import ExitStack
from itertools import batched, chain, islice
from multiprocessing import get_context
from threading import Condition


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
    results = []
    for model, protein in rows:
        result = checked_prediction(model, memo, genome, code, params, protein)
        # The caller already owns this evidence. Returning it over IPC would
        # create a separate copy for each batch, held until rescue publication.
        evidence = model['evidence']
        if result['evidence'] == evidence:
            del result['evidence']
        raw = result['raw_prediction']
        if raw['evidence'] == evidence:
            del raw['evidence']
        results.append(result)
    return results, os.getpid(), memo.diagnostics()


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
            inputs, future = item
            results, pid, snapshot = future.result()
            if len(results) != len(inputs):
                raise ValueError('Cached prediction validation batch changed record count')
            snapshots[pid] = snapshot
            for result, model in zip(results, inputs, strict=True):
                evidence = model['evidence']
                if 'evidence' not in result:
                    result['evidence'] = evidence
                if 'evidence' not in result['raw_prediction']:
                    result['raw_prediction']['evidence'] = dict(evidence)
                yield result
        remaining = chain(initial, models)
        while batch := list(islice(remaining, 64)):
            rows = [(model, proteins[model['evidence']['donor']][model['evidence']['query']]) for model in batch]
            pending.append((batch, pool.submit(_batch, rows)))
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


def native_validation_workers(cpus):
    """Bound child processes for donor-specific QC."""
    return min(2, (cpus - 1) // 3)


class NativeWorkerBudget:
    """Share all child CPU slots between exact scoring threads and repair jobs."""

    def __init__(self, slots, repair_workers):
        if type(slots) is not int or type(repair_workers) is not int or not 1 <= repair_workers <= slots:
            raise ValueError('Invalid shared native worker budget')
        self.slots = slots
        self.repair_workers = repair_workers
        self.available = slots
        self.repairs = 0
        self.waiting_repair = False
        self.condition = Condition()

    def score(self, function, *args):
        with self.condition:
            self.condition.wait_for(lambda: self.available and not (
                self.waiting_repair and self.repairs < self.repair_workers))
            self.available -= 1
        try:
            return function(*args)
        finally:
            with self.condition:
                self.available += 1
                self.condition.notify_all()

    def _finish_repair(self, _future=None):
        with self.condition:
            self.repairs -= 1
            self.available += 1
            self.condition.notify_all()

    def submit_repair(self, pool, function, rows):
        with self.condition:
            self.waiting_repair = True
            try:
                self.condition.wait_for(lambda: self.available and self.repairs < self.repair_workers)
                self.available -= 1
                self.repairs += 1
            finally:
                self.waiting_repair = False
                self.condition.notify_all()
        try:
            future = pool.submit(function, rows)
        except BaseException:
            self._finish_repair()
            raise
        future.add_done_callback(self._finish_repair)
        return future


def native_worker_allocation(cpus):
    """Use the shared budget on the thread-based SIMD path."""
    if __package__:
        from .rescue_alignment_kernel import compiled_kernel
    else:
        from rescue_alignment_kernel import compiled_kernel

    repairs = native_validation_workers(cpus)
    if not repairs:
        return cpus, None
    if compiled_kernel() is None:
        return max(1, cpus - repairs - 1), None
    return cpus - 1, NativeWorkerBudget(cpus - 1, repairs)


class _StagedValidation:
    def __init__(self, checked):
        self.checked = checked

    def validate(self, _model):
        return self.checked


def _finish_native_batch(rows):
    genome, code, params, _memo = _CONTEXT
    results = []
    for raw, checked, protein in rows:
        output = checked_prediction(raw, _StagedValidation(checked), genome, code, params, protein)
        results.append((output, output['evidence'] is raw['evidence']))
    return results


def validate_native(models, genome, code, params, proteins, cpus, memo, audit, *, budget=None, repair_workers=None):
    """Share genomic memo hits in the parent and send only donor repair to workers."""
    if __package__:
        from . import rescue_boundary_revalidation, rescue_terminal_completion
    else:
        import rescue_boundary_revalidation
        import rescue_terminal_completion

    workers = native_validation_workers(cpus) if repair_workers is None else repair_workers
    if type(workers) is not int or not 0 <= workers < cpus or (budget and workers != budget.repair_workers):
        raise ValueError('Invalid native repair worker count')
    audit.update(workers=0, offloaded=0, shared_worker_budget=budget.slots if budget else None)
    pending = deque()
    started = time.perf_counter()

    def emit(rows, indexes, future):
        if future is not None:
            for index, (output, shared_evidence) in zip(indexes, future.result(), strict=True):
                raw, _previous = rows[index]
                evidence = raw['evidence']
                if shared_evidence:
                    output['evidence'] = evidence
                output['raw_prediction']['evidence'] = dict(evidence)
                rows[index] = raw, output
        yield from rows

    with ExitStack() as stack:
        executor = None
        for batch in batched(models, 128):
            rows, indexes, repair = [], [], []
            for model in batch:
                evidence = model['evidence']
                protein = proteins[evidence['donor']][evidence['query']]
                checked = memo.validate(model)
                problems = set(checked.get('problems', ()))
                offload = (workers and len(batch) == 128 and not checked.get('frameshift') and (
                    (problems & rescue_terminal_completion._TERMINAL_PROBLEMS
                     and not problems - rescue_terminal_completion._RECHECKABLE)
                    or (problems and not problems - rescue_boundary_revalidation._RECHECKABLE)))
                if offload:
                    indexes.append(len(rows))
                    repair.append((model, checked, protein))
                    rows.append((model, None))
                else:
                    output = checked_prediction(model, _StagedValidation(checked), genome, code, params, protein)
                    rows.append((model, output))
            future = None
            if repair:
                if executor is None:
                    executor = stack.enter_context(ProcessPoolExecutor(
                        max_workers=workers, mp_context=get_context('spawn'), initializer=_initialise,
                        initargs=(os.fsdecode(genome.filename), code, params)))
                    audit['workers'] = workers
                audit['offloaded'] += len(repair)
                future = (executor.submit(_finish_native_batch, repair) if budget is None
                          else budget.submit_repair(executor, _finish_native_batch, repair))
            pending.append((rows, indexes, future))
            if len(pending) >= max(1, 2 * workers):
                yield from emit(*pending.popleft())
        while pending:
            yield from emit(*pending.popleft())
    audit['wall_seconds_including_scoring'] = time.perf_counter() - started
