#!/usr/bin/env python3
"""Compare bounded native shell pools with identical uneven worker durations.

This replays scheduling with sleep workers; it does not run BUSCO or MAFFT.
An optional private duration file contains a ``workers`` list with ``seconds``.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

SCRIPT = r'''
set -euo pipefail
source "$1"
root=$2
slots=$3
while IFS=$'\t' read -r index duration; do
  wait_until_jobn_le "${slots}"
  (
    printf '%s\tstart\t%s\n' "${EPOCHREALTIME}" "${index}" >> "${root}/events.tsv"
    sleep "${duration}"
    printf 'result-%s\n' "${index}" > "${root}/${index}.out"
    printf '%s\tend\t%s\n' "${EPOCHREALTIME}" "${index}" >> "${root}/events.tsv"
  ) &
  gg_background_register "$!"
done < "${root}/durations.tsv"
wait_for_background_jobs
[[ ${#GG_BACKGROUND_PIDS[@]} -eq 0 ]]
'''


def worker(helper: Path, durations: list[float], slots: int) -> dict:
    with tempfile.TemporaryDirectory(prefix="gg-background-pool-") as directory:
        root = Path(directory)
        (root / "durations.tsv").write_text("".join(
            f"{i:04d}\t{duration:.6f}\n" for i, duration in enumerate(durations)
        ))
        start = time.perf_counter()
        subprocess.run(["bash", "-c", SCRIPT, "bash", str(helper), str(root), str(slots)], check=True)
        seconds = time.perf_counter() - start
        proof = hashlib.sha256()
        outputs = sorted(root.glob("*.out"))
        assert len(outputs) == len(durations)
        for i, output in enumerate(outputs):
            assert output.read_text() == f"result-{i:04d}\n"
            proof.update(output.name.encode())
            proof.update(output.read_bytes())
        events = sorted(line.split("\t") for line in (root / "events.tsv").read_text().splitlines())
        active, maximum = set(), 0
        for _timestamp, kind, index in events:
            if kind == "start":
                assert index not in active
                active.add(index)
                maximum = max(maximum, len(active))
            else:
                active.remove(index)
        assert not active and maximum <= slots
        usage = resource.getrusage(resource.RUSAGE_CHILDREN)
        return {"seconds": seconds, "max_active": maximum, "outputs_sha256": proof.hexdigest(),
                "largest_child_peak_rss_bytes": usage.ru_maxrss * (1 if sys.platform == "darwin" else 1024)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-helper", type=Path, required=True)
    parser.add_argument("--current-helper", type=Path,
                        default=Path(__file__).resolve().parents[1] / "support/gg_util/06_workspace_validation.sh")
    parser.add_argument("--durations", type=Path)
    parser.add_argument("--slots", type=int, default=4)
    parser.add_argument("--trials", type=int, default=5)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.slots < 1 or not 1 <= args.trials <= 10:
        parser.error("Require positive slots and 1..10 trials")
    durations = [0.6, *([0.05] * 8)] * 4
    if args.durations:
        raw = json.loads(args.durations.read_text())
        values = [float(row["seconds"]) for row in raw["workers"]]
        if not values or min(values) <= 0:
            parser.error("Worker durations must be nonempty and positive")
        durations = [value / max(values) * 0.3 for value in values]
    if args.worker:
        print(json.dumps(worker(args.current_helper.resolve(), durations, args.slots)))
        return
    helpers = {"baseline": args.baseline_helper.resolve(), "current": args.current_helper.resolve()}
    result = {"python": sys.version, "platform": platform.platform(), "slots": args.slots,
              "workers": len(durations), "workload": "scaled recorded durations" if args.durations else "uneven synthetic durations",
              "durations_sha256": hashlib.sha256(json.dumps(durations).encode()).hexdigest(),
              "helper_sha256": {name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in helpers.items()},
              "samples": {name: [] for name in helpers}}
    for trial in range(args.trials + 1):
        order = list(helpers.items())
        if trial % 2:
            order.reverse()
        for name, helper in order:
            command = [sys.executable, str(Path(__file__).resolve()), "--worker", "--baseline-helper", str(helper),
                       "--current-helper", str(helper), "--slots", str(args.slots), "--output", str(args.output)]
            if args.durations:
                command.extend(["--durations", str(args.durations.resolve())])
            sample = json.loads(subprocess.check_output(command, text=True))
            print(name, "warmup" if trial == 0 else trial, sample, flush=True)
            if trial:
                result["samples"][name].append(sample)
    assert len({row["outputs_sha256"] for samples in result["samples"].values() for row in samples}) == 1
    result["medians"] = {name: {key: statistics.median(row[key] for row in samples)
                                for key in ("seconds", "largest_child_peak_rss_bytes")}
                         for name, samples in result["samples"].items()}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result["medians"]))


if __name__ == "__main__":
    main()
