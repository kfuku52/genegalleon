"""Bound queued I/O work without retaining every completed future."""

import weakref
from collections import deque
from itertools import islice
from queue import SimpleQueue


def ordered_bounded_results(executor, function, arguments, max_pending):
    """Yield input-order results while bounding queued and completed work.

    Arguments are tuples, as for bounded_results. The limit includes completed
    results waiting behind a slow earlier task. Closing cancels queued work;
    the caller still owns the executor and waits for active tasks.
    """
    if max_pending < 1:
        raise ValueError("max_pending must be positive")
    arguments = iter(arguments)
    pending = deque()
    try:
        for args in islice(arguments, max_pending):
            pending.append(executor.submit(function, *args))
        while pending:
            future = pending.popleft()
            yield future.result()
            del future
            for args in islice(arguments, 1):
                pending.append(executor.submit(function, *args))
    finally:
        for future in pending:
            future.cancel()


def bounded_results(executor, function, arguments, max_pending):
    """Yield completion-order results from at most max_pending argument tuples.

    Worker and argument-iterator errors propagate. Closing the iterator cancels
    queued work; the caller owns the executor and waits for active I/O to finish.
    """
    if max_pending < 1:
        raise ValueError("max_pending must be positive")
    arguments = iter(arguments)
    pending = set()
    completed = SimpleQueue()
    completed_ref = weakref.ref(completed)

    def record_completed(future):
        # An abandoned iterator must not be retained by active futures, or form
        # a queue -> completed-future -> callback -> queue reference cycle.
        queue = completed_ref()
        if queue is not None:
            queue.put(future)

    def submit(args):
        future = executor.submit(function, *args)
        pending.add(future)
        future.add_done_callback(record_completed)

    try:
        for args in islice(arguments, max_pending):
            submit(args)
        while pending:
            future = completed.get()
            pending.remove(future)
            yield future.result()
            for args in islice(arguments, 1):
                submit(args)
    finally:
        for future in pending:
            future.cancel()
