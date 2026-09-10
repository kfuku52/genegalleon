# 1,000 tips / 100 true shifts: native AICc search timing

Compare NWKIT `native-path` with its existing beam search, both using AICc.
This measures configured search procedures, including unknown alpha/process
variance estimation. No null gate, calibration or support bootstrap is run.
The software implementation is unchanged for this measurement.

## Inputs and scope

Six independent branch-recursion simulations: one balanced 1,000-tip tree,
100 distinct, disjoint shifted clades (7–8 tips each), one trait, fixed-root OU,
true alpha*height 3, process tip variance 1, effect scales 2/6 and seeds
30101–30103. Locations and noise are paired across effect strengths within a
seed. `input-qa.json` verifies counts, disjointness, taxa, input hashes and paired
noise. The simulator does not call the inference package.

All branches are positive and the tree height is one. There are no measurement
errors, missing values or shared regimes. The search cap is 100, equal to truth;
this is not an unknown-count or null false-positive validation study. Recovery
uses exact descendant-tip clades, not package-specific branch IDs.

## Search budgets

Both methods use the documented large-search configuration:

```
criterion=AICc
max_shifts=100
candidate_pool=128
refit_budget=220
screening_budget=20000
beam_width=1
lasso_iterations=150
search_memory_mb=512
```

The small default pool/refit settings cannot accommodate 100 shifts, so the
existing explicit large-search configuration is used. Pool and beam width govern
beam search; path search considers all branches. Both share the refit/screening
limits, but path points and beam layout screens are different units of work.
The comparison reports actual candidate counts and the largest fitted model.
The 512-MiB search-memory guard is not a limit on total process RSS.

Finite-budget candidate paths can stop before convergence; their status is
retained in the raw outputs. Final coefficients are refitted without shrinkage.
The selected models and recovery can differ between methods, so timing ratios
must not be interpreted as equivalent-output algorithmic speedups.

## Measurement and reproduction

Primary measurements are sequential fresh processes with alternating method
order across the six datasets. One 100-tip run per method first warms library
and file caches; those two recorded warmups are excluded from the primary
statistics. Three independent datasets per effect provide timing repetition;
this is not a broad statistical-power study. Wall time includes startup, reading,
search and JSON serialization. Linux `wait4` records CPU and peak process RSS.
All BLAS/OpenMP thread variables are one. The assistant runs no concurrent test
suite or other benchmark; the host itself is shared. Each run has an 1,800-second
ceiling; timeouts and failures remain in raw results.

The same GeneGalleon Docker image and frozen NWKIT source serve both methods.
`image-id.txt`, `source.tar.gz` and `source-manifest.json` capture them.
`protocol.json` records the experiment settings. From this repository, mount
NWKIT at `/src` and the repository at `/gg` and run:

```sh
python docs/benchmarks/native-ou-1000tips-100shifts-aicc/benchmark.py
python docs/benchmarks/native-ou-1000tips-100shifts-aicc/run.py
python docs/benchmarks/native-ou-1000tips-100shifts-aicc/summarize.py
python docs/benchmarks/native-ou-1000tips-100shifts-aicc/verify.py
```

Use `PYTHONPATH=/src` inside `local/genegalleon:native-ic-dev`, with
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
VECLIB_MAXIMUM_THREADS=1`. Generation and measurement refuse existing outputs;
replay in a clean copy. The two cache-warmup outputs in this study were produced
with the prior AICc adapter on `effect6-seed27101`, once per method.
This is Docker validation; SIF was not tested.
