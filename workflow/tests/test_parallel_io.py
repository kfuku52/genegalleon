"""Queued metadata I/O stays bounded and fails without scheduling the remainder."""
from concurrent.futures import Future, ThreadPoolExecutor
from threading import Event

import pytest

from workflow.support.parallel_io import bounded_results


@pytest.mark.parametrize('limit', [1, 4])
def test_bounded_results_does_not_consume_all_arguments(limit):
    release = Event()
    filled = Event()
    consumed = []

    def arguments():
        for value in range(100):
            consumed.append(value)
            if len(consumed) == limit:
                filled.set()
            yield (value,)

    def worker(value):
        assert release.wait(10)
        return value * 2

    with ThreadPoolExecutor(max_workers=2) as io_pool, ThreadPoolExecutor(max_workers=1) as consumer:
        results = bounded_results(io_pool, worker, arguments(), max_pending=limit)
        first = consumer.submit(next, results)
        try:
            assert filled.wait(10)
            assert consumed == list(range(limit))
        finally:
            release.set()
        output = [first.result(timeout=10), *results]
    assert sorted(output) == list(range(0, 200, 2))
    assert consumed == list(range(100))


class ManualExecutor:
    def __init__(self):
        self.futures = []

    def submit(self, function, *arguments):
        future = Future()
        if not self.futures:
            try:
                future.set_result(function(*arguments))
            except Exception as exc:
                future.set_exception(exc)
        self.futures.append(future)
        return future


def test_close_cancels_pending_work_without_reading_remainder():
    pool = ManualExecutor()
    arguments = iter((value,) for value in range(100))
    results = bounded_results(pool, lambda value: value, arguments, 4)
    assert next(results) == 0
    results.close()
    assert len(pool.futures) == 4
    assert all(future.cancelled() for future in pool.futures[1:])
    assert next(arguments) == (4,)


def test_active_worker_can_finish_after_iterator_is_closed(caplog):
    pool = ManualExecutor()
    results = bounded_results(pool, lambda value: value, ((value,) for value in range(100)), 4)
    assert next(results) == 0
    assert pool.futures[1].set_running_or_notify_cancel()
    results.close()
    pool.futures[1].set_result(1)
    assert pool.futures[1].result() == 1
    assert not caplog.records
    assert all(future.cancelled() for future in pool.futures[2:])


def test_worker_failure_propagates_and_cancels_queued_work():
    pool = ManualExecutor()

    def fail(value):
        raise OSError('broken metadata')

    with pytest.raises(OSError, match='broken metadata'):
        list(bounded_results(pool, fail, ((value,) for value in range(100)), 4))
    assert len(pool.futures) == 4
    assert all(future.cancelled() for future in pool.futures[1:])


def test_argument_failure_cancels_submitted_work():
    pool = ManualExecutor()

    def arguments():
        yield (0,)
        yield (1,)
        raise ValueError('invalid inventory')

    with pytest.raises(ValueError, match='invalid inventory'):
        list(bounded_results(pool, lambda value: value, arguments(), 4))
    assert pool.futures[1].cancelled()


@pytest.mark.parametrize('limit', [0, -1])
def test_invalid_queue_limit_is_rejected(limit):
    with pytest.raises(ValueError, match='positive'):
        list(bounded_results(ManualExecutor(), lambda: None, [], limit))


def test_empty_inventory_needs_no_work():
    pool = ManualExecutor()
    assert list(bounded_results(pool, lambda: None, [], 4)) == []
    assert pool.futures == []
