"""Bound queued I/O work without retaining every completed future."""

import weakref
from itertools import islice
from queue import SimpleQueue


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
