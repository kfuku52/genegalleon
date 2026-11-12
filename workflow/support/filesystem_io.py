"""Small filesystem primitives shared by otherwise independent publication protocols."""

import os
import tempfile
from contextlib import contextmanager
from pathlib import Path


@contextmanager
def scratch_build_directory(staging: Path):
    """Build on temporary storage when it is on a different filesystem."""
    scratch = Path(tempfile.gettempdir())
    if scratch.stat().st_dev == staging.stat().st_dev:
        yield staging
        return
    with tempfile.TemporaryDirectory(prefix="gg-build-", dir=scratch) as directory:
        yield Path(directory)


def fsync_directory(path: Path) -> None:
    try:
        descriptor = os.open(path, os.O_RDONLY)
    except OSError:
        return
    try:
        os.fsync(descriptor)
    except OSError:
        pass
    finally:
        os.close(descriptor)
