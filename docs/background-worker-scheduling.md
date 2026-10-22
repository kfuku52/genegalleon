# Background worker scheduling

The shared shell pool used by BUSCO, MAFFT and other independent workers now
reuses a slot when any registered worker finishes. Previously it waited for
the oldest registered PID, leaving slots idle behind a slow worker even when
later workers had finished. The worker commands, launch order, CPU/memory caps,
stage order and scientific parameters are unchanged. Every registered child's
exit status still contributes to the final failure result; unrelated children
are neither waited on nor reaped. Completed children are collected before
waiting, including completion between the job snapshot and `wait`.

The completion-identified wait uses Bash 5.1 or newer, as supplied by the tested
GeneGalleon runtimes. Earlier supported Bash versions retain their existing
FIFO wait behavior. A running immutable workflow does not acquire this change
by updating another checkout.

## Comparable scheduling measurements

Linux arm64 Docker, Python 3.12.15, Bash 5.2.21, one warmup and five alternating
trials compared the helper from `cd032844` with the changed helper in the same
immutable runtime. Four slots processed 36 sleep workers: four repetitions of
one 0.6-second worker followed by eight 0.05-second workers.

| Measurement | FIFO | Completion-based wait |
|---|---:|---:|
| Median wall time | 2.701 s | 1.322 s |
| Wall-time range | 2.696–2.717 s | 1.319–1.330 s |
| Largest-child peak RSS, median | 18.59 MiB | 18.58 MiB |

Wall time fell 51.1%. Each trial checked identical complete worker output bytes,
all workers collected, and at most four simultaneous workers. RSS is the
largest individual child recorded by `RUSAGE_CHILDREN`, including the benchmark
harness; it is neither aggregate simultaneous RSS nor a scientific-tool memory
measurement. The fixture measures scheduling with uneven durations, not BUSCO,
alignment quality, IQ-TREE search or complete workflow throughput.

The same fixture also ran in the pinned amd64 SIF used by the observed analysis
(Linux 6.11, Bash 5.2.21). One warmup and three alternating trials measured
2.709 to 1.339 seconds, with identical worker outputs and the same four-slot
bound. This is a SIF helper comparison, not a complete scientific workflow run.

Reproduce from the checkout without concurrent builds or tests:

```bash
mkdir -p tmp/background-pool
git show cd032844:workflow/support/gg_util/06_workspace_validation.sh > tmp/background-pool/baseline.sh
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_background_pool.py \
  --baseline-helper tmp/background-pool/baseline.sh --trials 5 \
  --output tmp/background-pool/comparison.json
```

An optional `--durations FILE` reads a private JSON object containing a `workers`
list of positive `seconds` values. It scales the maximum to 0.3 seconds and
preserves list order. Native BUSCO-reported durations retained in completion
order describe a timing distribution, not a reconstruction of species launch
order. Filesystem modification times from another clock must not be combined
with process log timestamps to infer worker duration.
