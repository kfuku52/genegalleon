# Rescue dependency-key verification

Rescue finalization already keeps a command-local receipt proof, but comparison
key construction used a fresh full preparation check for each comparison pair.
It now passes the same proof into those nested preparation checks. Each prepared
publication is read and verified fully on its first use, then its original file,
permission, ancestor and resolved-path identities remain fenced until exit.

This uses the existing `ReceiptSnapshot` boundary. It creates no persistent
checksum cache, changes no scientific selection, and keeps the same complete
ordered comparison and rescue keys. Calls without a snapshot retain ordinary
fresh full-byte verification. Combining the independent preparation and receipt
contexts is rejected explicitly.

## Measurement

The synthetic workload has 23 species, 158 comparison pairs, one 8 MiB prepared
payload per species and 20,000 frozen plan metadata entries. It constructs every
species' rescue key twice in one receipt scope, as in the two finalization guards.
One warmup per side precedes three alternating paired trials; each trial runs in
a fresh Python process. The fixture contains no project sequences or models.

On Linux arm64 / Python 3.12.15 in Docker:

| Metric | Before | After |
| --- | ---: | ---: |
| Median dependency-key wall time | 6.907 s | 1.389 s |
| Hashed bytes | 12,063,780,403 | 1,460,579,891 |
| Content reads | 4,878 | 3,614 |
| Median maximum RSS | 71,618,560 bytes | 71,340,032 bytes |

The measured operation is 4.97 times faster (79.9% less wall time), with 87.9%
fewer hashed bytes. Memory use is essentially unchanged. These figures do not
measure complete rescue or refinement, whose other computation, verification and
I/O remain. They do not measure an HPC filesystem or an already running job.

All eight warmup/measured processes produced complete ordered key SHA-256
`ab42e2cbf9241048a4cd63b06046d0e6b11f9bc0b9284b0577e777a57054c1f3`.
The [full measurements](benchmarks/rescue-dependency-keys-2026-10-22.json) retain
individual trials, the exact baseline, runtime identity and command.

Reproduce with a fresh GeneGalleon runtime:

```bash
GG_TEST_RUNTIME=docker bash workflow/tests/run_in_runtime.sh \
  python -B workflow/benchmarks/benchmark_rescue_dependency_keys.py \
  --compare e8f3f7ab9b7467a7cfda9588c4ee2301933ad10e \
  --repeats 3 --warmups 1
```

Regressions check exact keys, one full preparation proof per command, fresh
proofs in a new command, unchanged ordinary library calls, initial corruption,
and refusal of observed content, receipt, permission, replacement and symlink
changes. Existing real-tool rescue and refinement checks remain required.
