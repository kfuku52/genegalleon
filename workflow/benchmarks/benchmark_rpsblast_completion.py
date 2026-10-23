#!/usr/bin/env python3
"""Compare actual old/new shell postprocessing; require byte-identical results."""
import argparse
import hashlib
import json
import platform
import statistics
import subprocess
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CORE = "workflow/core/gg_gene_evolution_core.sh"
FIELDS = "qacc sacc pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen stitle"


def block(source):
    begin = source.index("  if ! rpsblast ")
    marker = "  # Complete query coverage" if "  # Complete query coverage" in source else "  genes=()"
    start = source.index(marker, begin)
    end = source.index('\n  cp_out "${og_id}.rpsblast.tsv"', start)
    return 'set -euo pipefail\nog_id=OG\noutfmt="' + FIELDS + '"\n' + source[start:end]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", required=True)
    parser.add_argument("--queries", nargs="+", type=int, default=[200, 1000])
    parser.add_argument("--trials", type=int, default=3)
    args = parser.parse_args()
    assert args.trials >= 1 and all(0 < n <= 10000 for n in args.queries)
    scripts = [block(subprocess.check_output(["git", "show", args.base + ":" + CORE], cwd=ROOT, text=True)), block((ROOT / CORE).read_text())]
    rows = []
    for count in args.queries:
        timings = [[], []]
        hashes = set()
        with tempfile.TemporaryDirectory(prefix="gg-rps-completion-") as temporary:
            folder = Path(temporary)
            (folder / "ungapped_translated_cds.fas").write_text("".join(f">Species_g{i:05d}\n" + "MACDE" * 60 + "\n" for i in range(count)))
            hits = "".join(f"Species_g{i:05d}\tpfam\t80\t20\t0\t0\t1\t20\t1\t20\t1e-8\t42\t300\t40\tdomain\n" for i in range(0, count, 2))
            for trial in range(args.trials + 1):
                for version in ([0, 1] if trial % 2 == 0 else [1, 0]):
                    (folder / "OG.rpsblast.tmp.tsv").write_text(hits)
                    started = time.perf_counter()
                    result = subprocess.run(["bash", "-c", scripts[version]], cwd=folder, capture_output=True, check=True)
                    elapsed = time.perf_counter() - started
                    hashes.add(hashlib.sha256((folder / "OG.rpsblast.tsv").read_bytes() + result.stdout).hexdigest())
                    if trial:
                        timings[version].append(elapsed)
        assert len(hashes) == 1, "Scientific table or query log changed"
        medians = [statistics.median(t) for t in timings]
        rows.append(dict(queries=count, trials=args.trials, seconds=timings, median_seconds=medians, speedup=medians[0] / medians[1], full_table_and_logs_identical=True, sha256=hashes.pop()))
    print(json.dumps(dict(platform=platform.platform(), seqkit=subprocess.check_output(["seqkit", "version"], text=True).strip(), baseline=args.base, results=rows), indent=2))


if __name__ == "__main__":
    main()
