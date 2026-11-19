"""Compare selection and plot outputs for wide CSUBST site-summary tables.

Run under GeneGalleon with --support pointing at the baseline or updated support
directory. Run in a fresh process per measurement; compare the digest and row
counts as well as time and peak process RSS.
"""
import argparse
import contextlib
import logging
import os
import resource
import sys
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory
from time import perf_counter

import numpy as np
import pandas as pd


def selection(wrapper, *, sparse=False):
    rng = np.random.default_rng(119)
    high_parts, low_parts = [], []
    families, per_family = (2400, 4) if sparse else (4, 2400)
    for _ in range(families):
        groups = np.array([rng.choice(400, 4, replace=False) for _ in range(per_family)])
        high_parts.append(groups)
        pairs = np.concatenate([groups[:600, :2], groups[300:900, 2:], groups[:, [0, 3]][::-1]])
        rng.shuffle(pairs)
        # Include non-subsets as well as pairs within known quartets.
        pairs[::2, 1] += 400
        low_parts.append(pairs)
    tables = []
    for parts in (high_parts, low_parts):
        frame = pd.DataFrame(np.vstack(parts), columns=[f"branch_id_{i+1}" for i in range(parts[0].shape[1])])
        frame["orthogroup"] = np.repeat(np.arange(len(parts)), len(parts[0]))
        tables.append(frame)
    digest = sha256()
    rows = []
    state = {"trait": {}}
    with open(os.devnull, "w") as sink, contextlib.redirect_stdout(sink):
        start = perf_counter()
        for frame in tables:
            selected, state = wrapper.skip_lower_order(frame, len(frame.columns) - 1, "trait", state)
            rows.append(len(selected))
            digest.update(pd.util.hash_pandas_object(selected).to_numpy().tobytes())
        seconds = perf_counter() - start
    return seconds, rows, digest.hexdigest()


def figures(wrapper):
    import matplotlib
    import matplotlib.pyplot as plt
    # Remove file metadata randomness, not any plotted content.
    os.environ["SOURCE_DATE_EPOCH"] = "0"
    matplotlib.rcParams["svg.hashsalt"] = "119"
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    rng = np.random.default_rng(119)
    frame = pd.DataFrame(rng.uniform(0, 12, size=(80000, 64)))
    digest = sha256()
    with TemporaryDirectory() as td:
        start = perf_counter()
        for i in range(8):
            output = Path(td) / str(i)
            wrapper.plot_scatter(frame, 0, 1 + i, str(output), polygon_xmin=2, polygon_ymin=5)
            for suffix in (".pdf", ".svg"):
                digest.update(output.with_suffix(suffix).read_bytes())
        seconds = perf_counter() - start
    return seconds, len(plt.get_fignums()), digest.hexdigest()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=("selection", "sparse", "figures"))
    parser.add_argument("--support", type=Path, default=Path(__file__).resolve().parents[1] / "support")
    args = parser.parse_args()
    sys.path.insert(0, str(args.support.resolve()))
    import csubst_site_wrapper as wrapper

    seconds, count, digest = (figures(wrapper) if args.case == "figures"
                              else selection(wrapper, sparse=args.case == "sparse"))
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    mib = rss / (1024 ** 2 if sys.platform == "darwin" else 1024)
    print(f"{args.case}\t{seconds:.4f}\t{mib:.1f}\t{count}\t{digest}")
