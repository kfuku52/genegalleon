"""Shared species-validation scheduling and file discovery; validators own scientific checks."""

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

if __package__:
    from .parallel_io import bounded_results
else:
    from parallel_io import bounded_results


def has_known_suffix(name, suffixes):
    lower = name.lower()
    return any(lower.endswith(suffix) for suffix in suffixes)


def list_nonhidden_files(directory, suffixes):
    out = []
    for path in sorted(Path(directory).iterdir()):
        if not path.is_file():
            continue
        if path.name.startswith("."):
            continue
        if has_known_suffix(path.name, suffixes):
            out.append(path)
    return out


def resolve_nthreads(args):
    if args.nthreads is not None:
        return max(1, int(args.nthreads))
    if args.ncpu is not None:
        return max(1, int(args.ncpu))
    return 1


def run_validation_tasks(tasks, missing_limit, nthreads, validator, executor_error):
    if not tasks:
        return []
    if nthreads <= 1 or len(tasks) <= 1:
        return [validator(task=task, missing_limit=missing_limit) for task in tasks]
    ordered_results = {}

    def validate(task):
        try:
            return task, validator(task, missing_limit), None
        except Exception as exc:
            return task, None, exc

    with ThreadPoolExecutor(max_workers=min(nthreads, len(tasks))) as executor:
        for task, result, exc in bounded_results(
            executor, validate, ((task,) for task in tasks), max_pending=2 * nthreads,
        ):
            if exc is not None:
                result = executor_error(task, exc)
            ordered_results[result["index"]] = result
    return [ordered_results[index] for index in sorted(ordered_results)]
