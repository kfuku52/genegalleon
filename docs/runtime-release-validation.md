# Immutable runtime verification

The explicit-release bootstrap reads and hashes every managed runtime file on
each invocation. Sorted paths share open directory descriptors; each child open
rejects symlinks, and final path/identity checks reject directory replacement,
file changes and manifest changes during verification. Only the current directory
ancestry stays open. File reads are bounded to 1 MiB chunks. Root aliases retain
their existing logical-manifest/physical-directory interpretation.

Verification measures successive groups of eight reads. When their mean cost
reaches 10 ms, independent file reads overlap with at most four threads, further
bounded by positive `GG_TASK_CPUS`, `SLURM_CPUS_PER_TASK`, `NSLOTS` and `PBS_NP`
limits. Fast storage and one-CPU allocations stay serial. Each directory stays
open until its ordered results are consumed; final identity fences still cover
all files and directories. The read buffers therefore use at most 4 MiB.

This reduces repeated path traversal on shared filesystems without replacing
content verification with saved timestamps or weakening the managed inventory.
Scheduler startup preflight and subsequent native startup retain separate proofs.

## Bounded measurement

A Linux/arm64 Docker run with Python 3.12 used a 476-file fixture: 26 native
support files and 450 1,040-byte files across five shared directories. After one
warmup, five alternating baseline/candidate trials returned identical workflow
paths and hashed identical bytes, including contract/snapshot payloads.

| Metric | Previous verifier | Descriptor traversal |
| --- | ---: | ---: |
| Median validation time | 0.0361 s | 0.0161 s |
| Python `stat`/`lstat` calls | 7,072 | 989 |
| Median process peak RSS | 32,488 KiB | 32,488 KiB |

With 12 ms of simulated input-delivery latency per nonempty hash input, three
alternating measured trials after warmup took 6.995 / 1.877 seconds, with equal
outputs and hashed bytes. This controlled delay exercises the adaptive branch;
it is not a measurement of NFS or disk latency. Add `--input-delay-ms 12` to the
command below to reproduce it.

The filesystem-call counts cover instrumented Python `os.stat`/`os.lstat`, not
every kernel operation or descriptor `fstat`. This warm synthetic measurement
does not establish NFS startup time, SIF compatibility or scientific workflow
speed. The known image was used deliberately offline because the measured
bootstrap uses host Python standard-library operations; runtime freshness was
not claimed. Real startup qualification remains required before submission.

To reproduce with a saved previous `gg_entrypoint_bootstrap.sh` mounted into the
selected validation runtime:

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_runtime_release_validation.py \
  --baseline /data/baseline_gg_entrypoint_bootstrap.sh --files 450 --repeats 5
```

The benchmark records fresh hashed bytes and checks output equality. Its default
2,000-file case also exercises larger runtime inventories. Regression tests cover
regular content tampering, symlinked leaves/ancestors, nonregular files, directory
replacement, earlier-file changes during a later read, and manifest mutation.
Latency fixtures also verify the one-, two- and four-thread bounds and detection
of slow reads that begin after an initially warm part of the inventory.
