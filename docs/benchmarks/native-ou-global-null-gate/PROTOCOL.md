# Global-null AIC gate: development pilot

Prespecified method: T = max(0, AIC(null) - min searched AIC), simulated under
the fitted null, with complete native-path replay and covariance refitting.
Retain the AIC winner only when the plus-one bootstrap p-value is <= 0.05.
No conditional false-branch control is claimed.

Inputs are the independently simulated, frozen fixtures in
`../native-ou-aic-improvement/data` and `../native-ou-aic-improvement/extension/data`.
Input hashes are saved per result. All have 100 tips, true alpha*height 3 and
unit process tip variance. Search cap 10, refit budget 48, screening budget 2000,
150 lasso iterations; all other native defaults retained. No tuning from results.

- Diagnosis: all 30 null fixtures (seeds 28101–28130), fixing alpha to its true
  value while estimating process variance; compare ordinary AIC to estimated-alpha AIC.
- Paired gate pilot: those 30 null fixtures; six 10-shift fixtures (effects 2/6,
  seeds 28101–28103); all twelve 5-shift extension fixtures (balanced/random,
  effects 2/6, seeds 29101–29103).
- Inner B=19, level 0.05, bootstrap seed = outer seed + 71000. This minimum
  resolution is a development pilot, not a 5% calibration certification. The
  independent outer fixtures do not use the fitting package's simulator.
- Outcomes: probability of any selected shift, mean selected count, exact-clade
  true/false positives, precision, recall, tip-mean RMSE, search/gate wall time.
  Undefined precision for no selected branches is stored as null.
- Four worker processes, one BLAS/OpenMP thread each; timings are under this
  concurrency and must not be compared with earlier isolated-process timings.

Run from the GeneGalleon repository:

```sh
docker run --rm -v /Users/kf/repos/nwkit:/src:ro -v /Users/kf/repos/genegalleon:/gg \
  -e PYTHONPATH=/src -e OPENBLAS_NUM_THREADS=1 -e OMP_NUM_THREADS=1 -e MKL_NUM_THREADS=1 \
  -w /gg local/genegalleon:native-ic-dev python docs/benchmarks/native-ou-global-null-gate/run.py
```

The script refuses to overwrite result files. The main pilot estimates nuisance
parameters anew for every inner draw. It includes the unchanged AIC result from
the exact same observed search for paired comparison. No resampling failures may
be dropped. Docker validation does not establish SIF compatibility.

Additional default-count smoke check: run `run.py --default-check` to use B=199
on the first main fixture at each effect (0, 2, 6; seed 28101). These three
prespecified representative cases check the default draw count; they do not
estimate its false-positive rate or power. Their filenames end in `-B199.json`
and they are excluded from the B=19 pilot summaries. Some timing measurements
overlap container test-suite execution; treat wall times as indicative costs,
not controlled performance comparisons. Peak memory was not measured.
