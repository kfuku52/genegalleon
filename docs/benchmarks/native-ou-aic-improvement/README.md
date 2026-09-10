# Native OU AIC candidate-search improvement

[Decision and limitations](DECISION.md) · [Independent results](RESULTS.md) · [Shape/count scope check](extension/RESULTS.md)
· [Protocol](PROTOCOL.md) · [Validation](validation/README.md)

This experiment separates the six original development datasets from 30 new
paired seeds at each of three strengths (zero, weak, strong). Three configured
AIC procedures run sequentially: original NWKIT beam search, updated native
optimum-increment path, and kfl1ou. A further 12 datasets check five shifts below
the cap and random ultrametric tree shapes. Neither independent set is used for
tuning. See the protocol for limitations and provisional comparison margins.

## Diagnosis and development

`crossfit.json` contains both implementations' refits of the original native,
kfl1ou, and planted layouts. Native likelihoods evaluated at the R-fitted alpha
agree within 2e-5 in all 18 comparisons. Free-alpha optima also agree except for
a documented upper-bound difference in three refits across two datasets
(the planted layouts and one identical native-selected layout). Four of six native-selected configurations
already have lower AIC than the kfl1ou-selected configuration when evaluated by
the same native fitter. Consequently, simply minimizing AIC more aggressively
cannot be assumed to improve recovery.

The original 24-branch pool omits some kfl1ou-selected branches in five of six
cases and covers 0–10 of the planted branches. Development ablations are saved
as `pilot-*.json`, including ineffective or adverse quick-scale and local-refine
variants. The chosen path retains joint lasso configurations, using
unstandardized OU optimum-increment columns and updated covariance. Penalized
coefficients are never final estimates. It uses the same likelihood, alpha
bounds, ordinary AIC formula and 48-fit budget. Its candidate restriction can
select a higher-AIC model with better recovery; it is not a global optimizer.

`variants.py` and `path_variant.py` preserve development prototypes only. The
independent candidate adapter imports the frozen candidate implementation in
`candidate/nwkit`, never those prototypes. The source manifest fingerprints all
Python package files; `freeze.json` records the pre-validation freeze. After
measurement, `final-integration.patch` restores the original auto dispatch while
retaining explicit native-path. The numerical implementation is unchanged; see
[the integration decision](DECISION.md).

## Reproduction

Use the GeneGalleon Docker image recorded in `image-id.txt` (development tag
`local/genegalleon:native-ic-dev`) and mount this directory at `/bench`. Python source
snapshots are supplied as baseline-source.tar.gz and candidate-source.tar.gz
for original and candidate imports; extract them in the fresh working directory. Do not
reuse completed output directories: the runners refuse to overwrite results.
For a clean run, copy the protocol, source snapshots and scripts into a new
working directory, then generate the datasets there:

```bash
tar -xzf baseline-source.tar.gz
tar -xzf candidate-source.tar.gz
docker run --rm -v "$PWD:/bench" local/genegalleon:native-ic-dev python /bench/benchmark.py
docker run --rm -v "$PWD:/bench" local/genegalleon:native-ic-dev python /bench/generate_extension.py
docker run --rm -v "$PWD:/bench" \
  -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 -e MKL_NUM_THREADS=1 \
  -e NUMEXPR_NUM_THREADS=1 -e VECLIB_MAXIMUM_THREADS=1 \
  local/genegalleon:native-ic-dev bash -c 'python /bench/run.py && python /bench/run.py extension'
docker run --rm -v "$PWD:/bench" local/genegalleon:native-ic-dev python /bench/summarize.py
docker run --rm -v "$PWD:/bench" local/genegalleon:native-ic-dev python /bench/summarize_extension.py
```

`check_baseline.py` additionally needs the original six-case benchmark mounted
read-only at `/old`. It confirms original selected clades, likelihoods and tip
means before any new timing. Linux `wait4` records CPU and peak process RSS;
wall-time polling resolution is 0.1 seconds. Every process has a 1,800-second cap.
These are configured-procedure comparisons with different outputs, not
language/kernel-speed measurements. This is Docker validation, not SIF validation.

The frozen R backend can also be reconstructed from the [previous benchmark
source manifest](../native-ou-100tips-10shifts-ic/source-manifest.json) and its
kfl1ou-source.patch. Source snapshots are experiment provenance, not upstream
version defaults. Apply final-integration.patch to the candidate snapshot only
when reproducing the delivered default-dispatch behavior; the benchmark adapter
explicitly selects the same path either way.
