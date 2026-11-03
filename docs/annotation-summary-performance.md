# Annotation import and branch-summary performance

Version 0.8.479 removes repeated work in the native annotation DB and
reconciliation-summary helpers. This does not change the scientific models,
thresholds, BH scope, schema, output file count or input verification.

GO tuples are validated once for identical complete ID/aspect/term/evidence
values within one import. The LRU holds at most 256 immutable results; keys
longer than 4,096 characters bypass it. Changed tuple values still undergo the
full validation, and each input is still hashed while it is parsed. The cache
does not survive a command. Wide annotation rows also avoid Python generator
overhead when collecting unchanged field values and batch sizes.

For summaries, product and unique-protein JSON and total extant-gene counts
are computed once per family. Event counts and BH results remain branch-local.
The TSV writer uses the established column order directly, and SQLite inserts
use batches of 1,000 rows within the same transaction. The final partial batch,
missing XML, quoted product strings and failed publication remain covered by
regression tests.

## Comparable measurements

On 2026-11-03, compare the two helpers against commit
`f20e6381ffd4419f4a083bfb7859b1a9ffe3539d` (0.8.478), using the same GeneGalleon
Docker runtime (arm64, Python 3.12.15, image
`sha256:f0632dbd422a3ee87e3c2a82c3612cf64413a06c932140e69c178a16ef73b708`).
Each version gets one warmup and three fresh-process measurements, alternating
baseline and candidate. Fixture creation and output-comparison hashes are outside
the measured interval. The interval includes native source-copy hashing,
annotation parsing and GO indexes, publication, native summary construction,
TSV/DB writes and manifest hashes. No scientific project data is used.

Both cases have 4,096 families, 24 species, 98,304 annotated genes, 589,824 GO
tuples, 32 additional annotation fields, 15 branches and 61,440 branch/family
rows. The first repeats GO tuples and has 24 literal products per family. The
control has a different GO tuple for every gene and empty products.

| Median wall time | Baseline | Candidate | Reduction |
| --- | ---: | ---: | ---: |
| Repeated GO, annotation import | 3.49 s | 3.01 s | 13.9% |
| Repeated GO, branch summary | 2.89 s | 2.12 s | 26.5% |
| Repeated GO, combined | 6.38 s | 5.13 s | 19.6% |
| Distinct GO control, annotation import | 3.01 s | 2.91 s | 3.4% |
| Distinct GO control, branch summary | 1.25 s | 0.97 s | 22.4% |
| Distinct GO control, combined | 4.26 s | 3.88 s | 9.0% |

Median peak RSS is essentially unchanged: 106.71 to 106.75 MiB for the first
case, and 105.06 to 105.16 MiB for the control. All five output files, including
both SQLite DBs, both TSVs and the manifest, have identical byte hashes across
all six measured runs of each case. These are synthetic helper measurements,
not a claimed speedup or ETA for the SHIROKANE project. Existing fixed job
sources do not receive this change automatically.

## Reproduce

Save the old two helper files under `tmp/annotation-baseline/` using `git show
f20e6381ffd4419f4a083bfb7859b1a9ffe3539d:workflow/support/<filename>`. Keep all
other support files from the checkout on the import path; these dependencies
were unchanged by this patch. Select a fresh runtime as described in
[development and tests](development-and-tests.md#choose-checks-for-a-change).

```bash
bash workflow/tests/run_in_runtime.sh python workflow/tests/benchmark_annotation_summary.py \
  --baseline-support tmp/annotation-baseline --families 4096 \
  --warmups 1 --samples 3 --output tmp/annotation-products.json

bash workflow/tests/run_in_runtime.sh python workflow/tests/benchmark_annotation_summary.py \
  --baseline-support tmp/annotation-baseline --families 4096 \
  --warmups 1 --samples 3 --distinct-go --empty-products \
  --output tmp/annotation-control.json
```

Add `--profile` only for a separate diagnosis run; profiling changes the timing
overhead. The driver rejects any complete output-byte mismatch.
