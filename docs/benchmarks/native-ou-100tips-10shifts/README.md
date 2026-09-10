# 100-tip / 10-shift paired development benchmark

A [new seven-method comparison](../native-ou-100tips-10shifts-ic/RESULTS.md)
uses the corrected pBIC backend and supersedes this historical baseline.

**Post-run correction:** the installed kfl1ou backend failed the
`ou-optimum-information-v1` pBIC capability probe. See
[the audit](pbic-postrun-audit.json). These measurements describe the old,
uncorrected pBIC backend and must not be interpreted as a comparison against
corrected pBIC. The baseline requires rerunning after its capability check passes.
The preflight was omitted from this benchmark's direct R adapter; package
version 3.0.9 alone did not establish corrected behavior.

This directory contains the frozen specification, independent simulator, method
adapters, input datasets, per-run logs, resource measurements and scientific
figure for an initial comparison of kfl1ou and NWKIT native shift selection.
It is not a production-adoption study.

See the [results and figure](RESULTS.md) and [individual-run metrics](metrics.csv).

## Design

Six datasets use one balanced 100-tip tree, three random seeds and two shift
strengths. Each has ten distinct, disjoint shifted clades of 6–7 tips. Each
shift's optimum is a random sign times the effect scale times a uniform
0.8–1.2 multiplier. Optima are distinct; no regime-merging search is requested.
Within each seed, both strengths share the tree, locations and Gaussian
innovations. All trees have height one, positive branches and fixed root state.
The true OU alpha is 3 and the marginal evolutionary tip variance is 1.
There is one trait, no missingness and no observation error.

The simulator implements the Gaussian branch recurrence directly, independently
of either inference package. Both methods estimate alpha and process variance
from the same observed tips; neither receives the planted shifts or true
covariance. The maximum allowed shift count is 10, equal to the truth. This is
an upper-cap benchmark and does not test recovery when the truth is unknown
relative to a larger search cap.

kfl1ou uses pBIC and its existing default candidate search. NWKIT uses the native
heuristic (24 candidate branches, 48 final refits, 2,000 screening evaluations,
beam width 2) and full-search plug-in bootstrap selection with B=19 at level
0.05. This is the minimum B resolving that level and is a coarse development
setting, below the GeneGalleon native default of 199. Calibration is part of
NWKIT's model selection, not an optional support calculation. Both exclude
support bootstrap and plotting from the timed work. The model selection rules
and optimizer bounds remain method-specific; this is a comparison of configured
procedures, not equivalent algorithms or a pure programming-language speed test.

## Measurement and accuracy

Methods run sequentially, with order alternating across paired datasets and
BLAS/OpenMP limited to one thread. Wall time includes process startup, input,
inference and result serialization. Linux `wait4` supplies process peak RSS in
MiB (2^20 bytes); it does not measure the entire Docker VM. The kernel returns
RSS even for a killed timed-out run. Timeout time and memory are censored, not
successful completion measurements.

Shift identities are sorted descendant-tip sets, avoiding NWKIT/R branch-number
mismatches. Exact-branch precision is TP/(TP+FP), recall is TP/10 and F1 is
2TP/(2TP+FP+FN). Precision is defined as zero for an empty predicted set.
RMSE compares fitted tip means with the simulator's noise-free expected means,
not the noisy observations. These metrics do not assess historical shift time,
parameter identifiability, confidence-interval coverage or false-positive rates
under the global null.

Accuracy is computed only for completed selections, with completion denominators
reported separately. An initial uncalibrated NWKIT candidate is saved for
search diagnosis; it is never substituted for a completed native selection.
Three replicates per strength are exploratory and do not justify population
confidence intervals or broad reliability claims.

The initial resource ceiling was five minutes. Once the first strong-effect
native run reached that ceiling, the ceiling was extended to twenty minutes
for its unchanged rerun and all pending runs, to obtain completed accuracy
measurements. Search budgets, calibration B, data and seeds were unchanged.
`first-pass-results.json` and `resource-extension.json` preserve that transition.
The first pass's early, interrupted second strong-effect native run is retained
as a partial JSON but is not treated as a completed selection or a method failure.
Final resource summaries use the completed follow-up measurement, not the sum
of an initial timeout and its rerun. The initial kfl1ou invocation also exposed
a harness omission (required postorder tree ordering); it was corrected using
`ape::reorder.phylo` and rerun, with the preflight error retained separately.

## Reproduction

Run in the recorded local GeneGalleon image with this directory mounted at
`/bench`; the image identity and numerical versions are recorded separately.
The image contains kfl1ou 3.0.9 and the native NWKIT development snapshot.
No source mount overrides either package during inference. For a clean rerun,
copy only the five Python/R runner and summary files plus `specification.json`
to a new directory (the Python files are `benchmark.py`, `native.py`,
`complete_comparison.py` and `summarize.py`, plus `kfl1ou.R`). Run the commands
below from that new directory.

```sh
docker run --rm -v "$PWD:/bench" local/genegalleon:native-ou-dev \
  python /bench/benchmark.py generate

docker run --rm -v "$PWD:/bench" \
  -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 -e MKL_NUM_THREADS=1 \
  -e VECLIB_MAXIMUM_THREADS=1 local/genegalleon:native-ou-dev \
  python /bench/benchmark.py run

docker run --rm -v "$PWD:/bench" \
  -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 -e MKL_NUM_THREADS=1 \
  -e VECLIB_MAXIMUM_THREADS=1 local/genegalleon:native-ou-dev \
  python /bench/complete_comparison.py

docker run --rm -v "$PWD:/bench" local/genegalleon:native-ou-dev \
  python /bench/summarize.py
```

Measurement scripts write their own results and the follow-up script preserves
the first pass only once. `benchmark.py generate` recreates inputs in a new directory; it deliberately
refuses to overwrite existing datasets. The source generator, input hashes and
raw results are retained. No SIF execution was performed.
