# Admission-update memory

The admission patcher previously serialized every changed model into a second
in-memory dictionary before writing any shard. A large OrionGeno replay failed
at that allocation under a 64 GiB address-space limit. The patcher now serializes
one changed record at a time and locates each shard's changed ordinals by binary
search. Scientific admission checks, record order and file formats are unchanged.

A synthetic fixture uses 20,000 changed records, 12,006 nt per coding sequence,
348 shards and gzip. Fixture construction and input decoding precede the timed
patch call. Each measurement uses a fresh Python process; peak RSS includes the
caller's complete decoded change set. Three runs use the same immutable Docker
image, Python 3.12.15, input generator and compression settings.

| Median | Before | After |
| --- | ---: | ---: |
| Peak RSS (KiB) | 678,680 | 397,128 |
| Patch wall time (s) | 2.492 | 2.040 |

Peak RSS decreased by 41.5%; elapsed time decreased by 18.1%. These are synthetic
patch-stage measurements, not a prediction of complete genome-analysis savings.
All six runs produced the same manifest SHA256, including all 352 output file
hashes: `b51d9b36481d67da0030811c7624bfef206d0fc156688d7df2f341da7de405db`.
The baseline module was extracted from commit
`47e24782c8534d34834c576f2c256f660e98e368`; its SHA256 is
`c98b3485fa6615be0bf54ac346a0dd804f50ead183d5192a4dbd4c7a5e3ab7a2`.

The retained fixture supports gzip and zstd. Run both implementations through
the same runtime wrapper, supplying the extracted baseline module as `MODULE`:

```bash
GENEGALLEON_DOCKER_EXTRA_BINDS=/absolute/baseline-directory \
GG_TEST_RUNTIME=docker bash workflow/tests/run_in_runtime.sh \
  python workflow/tests/benchmark_rescue_model_store_patch.py \
  --implementation MODULE --records 20000 --repeats 3 --codec gzip
```

Then repeat with `--implementation workflow/support/rescue_model_store_patch.py`.
The measured image ID is
`sha256:594a13182c692adfa393107910c67cf57591ff8127d7eaad7724f5ea552e83af`.
A regression test also bounds the number of simultaneously retained serialized
envelopes while checking the complete decoded output.
