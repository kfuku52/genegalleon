"""Hash every published regular output using one directory traversal."""
import os
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path

try:
    from input_generation_array_state import digest
except ImportError:
    from .input_generation_array_state import digest


def output_paths(root, excluded=("genome.fa", "genome.mpi")):
    """Yield regular files once, retaining cached readdir types and no symlinks."""
    root = Path(root)
    pending = [root]
    while pending:
        directory = pending.pop()
        with os.scandir(directory) as entries:
            for entry in entries:
                if entry.is_dir(follow_symlinks=False):
                    pending.append(Path(entry.path))
                elif entry.name not in excluded and entry.is_file(follow_symlinks=False):
                    yield Path(entry.path)


def hash_outputs(root, workers=1, excluded=("genome.fa", "genome.mpi"), hash_function=digest):
    """Retain full content/stat guards, skip symlinks, and bound pending hashes.

    Results are consumed as completed so large FASTA outputs cannot block the
    other workers from verifying small interval files.
    """
    return hash_paths(root, output_paths(root, excluded), workers=workers, hash_function=hash_function)


def hash_paths(root, paths, workers=1, hash_function=digest):
    """Verify a receipt's listed paths without rescanning its directories."""
    if type(workers) is not int or workers < 1:
        raise ValueError("Output hash workers must be a positive integer")
    root = Path(root)
    if workers == 1:
        return {str(path.relative_to(root)): hash_function(path) for path in paths}
    result, pending = {}, {}
    with ThreadPoolExecutor(max_workers=workers) as executor:
        def consume():
            completed, _ = wait(pending, return_when=FIRST_COMPLETED)
            for future in completed:
                relative = pending.pop(future)
                result[relative] = future.result()
        for path in paths:
            pending[executor.submit(hash_function, path)] = str(path.relative_to(root))
            if len(pending) >= workers * 2:
                consume()
        while pending:
            consume()
    return result
