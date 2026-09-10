# Native OU information-criterion comparison

See the [measured results and figure](RESULTS.md) and [individual-run metrics](metrics.csv).

A subsequent [AIC candidate-path study](../native-ou-aic-improvement/DECISION.md)
uses these six datasets for development and evaluates 30 new seeds per condition.
The present results remain the frozen original-search comparison.

This reruns the 100-tip, 10-shift experiment with seven configured procedures:
NWKIT native bootstrap, pBIC, AIC and BIC; kfl1ou pBIC, AIC and BIC. Ordinary AIC
is distinct from AICc. The repeated kfl1ou AIC in the request is counted once.

The previous benchmark used an uncorrected pBIC backend and remains archived in
[the original directory](../native-ou-100tips-10shifts/README.md). This comparison
requires the corrected pBIC capability probe and cross-backend score checks to
pass **before** timing. See [validation](validation/pbic-attestation.json) and
[fixed-alpha score agreement](validation/native-kfl-score-agreement.json).

## Experiment

The six inputs are byte-identical to the original experiment: one balanced
100-tip tree, three seeds, two effect scales (2 and 6), one trait, ten disjoint
shifted clades of 6–7 tips. True alpha times height is 3 and marginal process tip
variance is 1. The fixed-root Gaussian simulator propagates branch innovations
independently of both inference packages. All branches are positive; there is
no missingness, observation error or convergence. Locations and Gaussian noise
are paired across effect strengths within each seed.

Every procedure estimates its own alpha and process variance. All have a
10-shift cap, equal to truth; this does not test behavior with an unknown truth
below a larger cap. NWKIT uses its heuristic candidate search with pool 24,
48 refits, 2,000 screening evaluations and beam width 2. All four native methods
share those settings. The three information criteria score every refitted
candidate, including layouts that lose on likelihood within a shift count.
The likelihood-based candidate search is not a certified global IC optimizer.
kfl1ou retains its default candidate search and optimizer bounds.

Native bootstrap selection repeats the full search at B=19 and level 0.05.
This is the minimum B resolving 0.05 and a coarse development setting; the
workflow default is 199. IC methods do not run calibration draws. No method
runs optional support bootstrap in these measurements. Selection rules,
optimization bounds and candidate generation differ between packages, so
runtime ratios compare configured procedures, not language or kernel speed.

## Measurements and accuracy

The 42 runs execute sequentially in the same frozen GeneGalleon Docker image,
with BLAS/OpenMP restricted to one thread and cyclic method-order rotation
across datasets. Every run has a 1,800-second ceiling. Wall time includes
process startup, input, selection and JSON serialization. Linux `wait4`
provides CPU time and peak process RSS (MiB); this is not total Docker-VM memory.
Failures and censored timeouts remain in raw results and completion counts.
Wall-time collection polls every 0.1 seconds; differences of that order should
not be interpreted as meaningful speed differences.

Shift identities are sorted descendant-tip sets, avoiding package-specific
branch IDs. Precision is TP/(TP+FP), recall TP/10 and F1
2TP/(2TP+FP+FN). Empty predictions have precision zero. RMSE compares fitted tip
means with noise-free simulated expected means. The adjusted Rand index compares
tip partitions induced by the selected and planted shifts. Every individual
run is retained; summary accuracy is the mean over completed runs and resource
summaries are medians. With three seeds per effect, no population confidence
intervals or broad reliability claims are justified.

The design does not validate global-null false positives, historical timing,
convergence, multivariate behavior, observation errors, missing data or
production adoption. Equal fixed-layout scores do not imply equal search
procedures. Boundary information criteria also have nonregular statistical
limitations described in NWKIT's native inference documentation.

## Reproduction

[specification.json](specification.json) records the image digest and parameters;
[source-manifest.json](source-manifest.json) fingerprints the source files.
The NWKIT and kfl1ou source patches preserve the changes relative to their
recorded HEADs, including the corrected pBIC implementation.
[Validation evidence](validation/README.md) reports checks and existing failures.
[input-qa.json](input-qa.json) verifies unchanged inputs. Reuse the frozen image
or rebuild the development packages from the recorded source state in the
GeneGalleon runtime. The following commands require the benchmark directory
mounted at `/bench`, and a new output directory for another measured run:

```bash
docker run --rm -v "$PWD:/bench" local/genegalleon:native-ic-dev python /bench/validate.py
docker run --rm -v "$PWD:/bench" \
  -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 -e MKL_NUM_THREADS=1 \
  -e NUMEXPR_NUM_THREADS=1 -e VECLIB_MAXIMUM_THREADS=1 \
  local/genegalleon:native-ic-dev python /bench/run_comparison.py
```

`run_comparison.py` refuses to overwrite an existing results file. Run
`summarize.py` afterward to calculate metrics and export the scientific figure.
This is Docker validation on macOS; no SIF runtime was available or tested.
