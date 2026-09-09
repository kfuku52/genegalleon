# IQ-TREE persistent-session validation (2026-09-09)

The persistent connection reuses the loaded alignment, fixed model, topology and
allocated likelihood buffers. It reduces repeated-call overhead, most clearly
for small alignments. These measurements do **not** show a substantial speedup
for large codon alignments, or a speed advantage over NWKIT's native engine.

## Controlled likelihood/score benchmark

Both modes used the same locally extended IQ-TREE executable, the same simulated
alignment and branch-length requests, `GY{0.5,2}+FQ+G4{1}`, and one worker thread.
The benchmark alternated mode order over three repetitions. Each repetition
started a fresh Python process, fitted the initial unclocked tree, warmed the
likelihood kernel once, and timed 20 distinct uncached branch-length updates.
The table reports medians; ranges are the minimum and maximum of three runs.
Setup is separate from the timed requests and includes the unclocked prefit.

| Tips × codons | Subprocess: 20 requests (s) | Persistent: 20 requests (s) | Subprocess setup (s) | Persistent setup (s) |
| --- | ---: | ---: | ---: | ---: |
| 4 × 150 | 0.550 (0.526–0.950) | 0.256 (0.254–0.279) | 0.095 | 0.098 |
| 16 × 150 | 2.702 (2.647–2.728) | 2.384 (2.363–2.390) | 0.865 | 0.900 |
| 16 × 1,500 | 18.563 (18.159–19.562) | 18.355 (17.861–18.982) | 8.253 | 8.207 |
| 64 × 1,500 | 86.612 (85.716–86.900) | 86.725 (85.967–88.221) | 44.247 | 43.435 |

The small case took approximately half as long for repeated evaluations. For
64 tips the evaluation ranges overlapped, and setup was comparable;
there is no demonstrated meaningful overall gain for that 20-request workload.
Data retention alone does not remove the cost of codon likelihood/score
calculation. Long-lived workers also retain their memory between requests.

Largest-child peak RSS medians were approximately 142 MiB for the first three
workloads and 227 MiB for the last, with no substantial reduction in persistent
mode. This is Unix `RUSAGE_CHILDREN.ru_maxrss`, including fork/exec peaks; it is
not the sum of simultaneously resident Python/worker memory or a measurement
of IQ-TREE's steady-state allocations alone.

All 24 trials checked log-likelihood and score equivalence. The largest absolute
log-likelihood difference was `4.78e-7`; scores agreed within `rtol=1e-5` and
`atol=2e-4`, accounting for the subprocess IQ2MC export's six-significant-digit
scores. Persistent responses use 17 significant digits. The modes use the same
conditional model and do not substitute different frequency estimators.

## Dating-stage comparison

One controlled paired integration run used the existing GeneGalleon fixture:
four sequences, 150 simulated codons (seed 91), default `GY+F3X4+G4`, fitted
nuisance parameters and rate variance, fixed species ages, and 95% conditional
profile intervals. Wall time includes the real shell stage, NWKIT startup,
model fitting, dating, intervals, PDF generation and atomic bundle publication;
fixture generation is outside the timer.

| Connection | Complete dating stage (s) |
| --- | ---: |
| Subprocess | 112.116 |
| Persistent | 53.152 |

The frozen model string and interval statuses were identical. Maximum absolute
differences were `1.35e-5` in point ages (species-root age 10), `5.99e-10` in
interval endpoints and `8.31e-8` in estimated branch rates. The integration test
requires matching clade IDs and agreement within `rtol=1e-4`, `atol=1e-4`.
This is one paired end-to-end check, not a multi-family throughput benchmark;
do not generalize its approximately twofold improvement to large families.

To reproduce this pair from the GeneGalleon checkout:

```sh
docker run --rm \
  -v "$PWD:/gg:ro" \
  -v /tmp/iqtree-session-benchmark:/results \
  -w /gg local/genegalleon:iqtree-before-speedup \
  python -m pytest -q \
  workflow/tests/test_native_tree_dating.py::test_iqtree_session_matches_subprocess_profile \
  --basetemp=/results/stage-pair
```

Each mode's temporary fixture directory contains `stage_wall_seconds.txt` and
the dated output bundle. The full dating-stage test file also covers default
Gamma, FreeRate, the native codon engines and preservation of previous outputs
on failure.

## Numerical and integration validation

The timed worker includes a cancellation check in its spectral likelihood
calculation. Ill-conditioned requests are recomputed within IQ-TREE using its
original-state pruning kernels and a stable matrix exponential of the same
fitted generator. Neither branch lengths nor the requested model are changed.
This path retains the process and loaded data, but requires additional work.
The ordinary-input timings above do not measure the cost of that path.

A codon regression with two sibling branches near `8e-11` and `6e-11` checks
likelihoods and scores against an independent SciPy matrix exponential, then
returns to ordinary branch lengths in the same worker. On the original failure
reproducer the negative log likelihood was `1863.0148201671896`, compared with
`1863.0148201671898` from the independent calculation. The pre-fix spectral
worker returned a finite but incorrect value; merely checking finiteness was
insufficient. Empirical-model frequency normalization is also checked by
forcing the original-state path through all 39 IQ-TREE model/derivative tests.

The IQ-TREE producer protocol test, all 136 targeted NWKIT tests, and all three
affected GeneGalleon IQ-TREE stage tests passed on the corrected implementation.
NWKIT's full checks had passed with 3,052 tests and 37 skips before this C++
numerical correction; the affected runtime tests were rerun afterward. Lint,
type, security, maintainability and distribution-reproducibility checks passed.
GeneGalleon's broader fast suite reported 1,618 passes and three failures in
unrelated input-format/array-workflow changes; those files were not changed to
make this feature pass. SIF execution was unavailable on this macOS host.

The measured image has been retained as `local/genegalleon:iqtree-before-speedup`.
The development alias now advances with subsequent changes; see
[the follow-up optimization measurements](iqtree-speedup-validation.md).

## Environment and reproduction

- Docker Linux ARM64 (`Linux-6.12.76-linuxkit-aarch64`, glibc 2.39), on macOS.
- Python 3.12.14; IQ-TREE reports 3.1.4 with the local session extension.
- IQ-TREE binary SHA-256: `a7bd4e334b63c0ba83eccad793d795e330d111745ca10008fc65e5ba10139131`.
- Runtime image: `local/genegalleon:iqtree-before-speedup`, image SHA-256
  `78338e630bdbb90a65be12ca007c8a16576d265588c7aa12d1a7a3a4b6aa78dd`.
- Alignment simulation: seed 918, balanced gene trees, homogeneous GY94 with
  equal codon frequencies. The evaluation model includes four Gamma categories;
  this benchmark tests transport equivalence, not model adequacy or CI coverage.

Build the local overlay as described in [Gene-tree dating](gene-tree-dating.md).
From the matching NWKIT source checkout, run:

```sh
mkdir -p /tmp/iqtree-session-benchmark
docker run --rm \
  -v "$PWD:/nwkit:ro" \
  -v /tmp/iqtree-session-benchmark:/results \
  -w /nwkit local/genegalleon:iqtree-before-speedup \
  python tools/benchmark_radte_iqtree.py --output /results/benchmark.json
```

The JSON contains every timed trial, all likelihood/gradient outputs and binary
identity. Avoid concurrent CPU-intensive work during timing. These results cover
synthetic codon inputs, a single thread and this ARM64 build; they establish
neither native SIF compatibility nor performance on x86, other thread counts,
other models or an entire biological dataset.
