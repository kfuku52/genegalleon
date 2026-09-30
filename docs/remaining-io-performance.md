# Verification, database and plot I/O performance

GeneGalleon 0.8.40 removes repeated work without changing scientific settings,
public result schemas, SQLite contents, archive layout or plot geometry.

- `verify` shares declared-source digests across required steps and exact-attempt
  preflights. Standalone preflights also share their repeated contract reads.
  Each reused digest is independently recomputed and manifest evidence reread
  before returning; persistent cache entries cannot establish verification.
  Logical source identities, optional-output types/absence, archive generations,
  family observations and exact-attempt receipts remain enforced. Worker-thread
  final reads retain the explicit runtime mapping and read-only logical store.
- Store and subdirectory fingerprints enumerate metadata without first hashing
  raw contents. Observational fingerprints stream each member once, retaining
  source-change and archive-generation checks. Ordinary runtime fingerprints
  retain their existing persistent raw-file digest-cache behavior.
- Database input size lookup uses metadata rather than computing SHA256. Small
  raw and ZIP TSVs use the same bounded 256 KiB read path. Larger files retain
  their full-file type-discovery pass, including late text and nullable Boolean
  promotion. ZIP readers are pooled under a complete build snapshot, and the
  database is published only after the final archive fence. Create, replace and
  append modes keep their existing atomic publication and namespace guard.
- Tree plotting checks `ggimage` in the rendering R process rather than starting
  another R process. Its absence still skips plotting; other rendering failures
  fail the stage. Domain tables and trimmed/untrimmed FASTA data are parsed once
  per unique input in that family. Signature checks reject changes on reuse and
  before writing the PDF. Existing panel order, styling and parameters remain.

## Measurements, 2026-09-30

Baseline: `b60b10fed483e6306a583984b5a9968e20646cbc` (0.8.39).
Both revisions used the same `local/genegalleon:dev` Linux/arm64 Docker image
(`4174483dc4df`), Python 3.12, on macOS Docker Desktop. No additional CPU/memory
limit was applied. Each case used a fresh process, one discarded warmup and three
measured trials; the table reports median operation wall time. Fixture generation
and output comparison are outside the measured interval. All cases were serialized;
this task ran no tests or other benchmark builds concurrently.

| Case | Before (s) | After (s) | Speedup |
| --- | ---: | ---: | ---: |
| Exact-attempt verification, five steps sharing a 64 MiB input | 0.4520 | 0.0955 | 4.73x |
| Database, 256 families x 128 branches, raw | 0.6751 | 0.5481 | 1.23x |
| Same database, ZIP-held inputs | 6.0116 | 0.5647 | 10.64x |
| Observational whole-store fingerprint, 64 MiB FASTA plus statistics | 0.1655 | 0.0665 | 2.49x |
| PDF, tree panel and dependency check | 1.5157 | 1.1699 | 1.30x |
| PDF, tree/domain/alignment panels, shared trimmed/untrimmed input | 1.8543 | 1.6865 | 1.10x |

All trials produced equal verification decisions and scientific contract payloads,
equal hashes of sorted SQLite contents and table schemas, equal store digests,
and byte-identical PDF contents after removing only creation/modification dates
(including compressed drawing commands, fonts and page sizes).

Median parent-process peak RSS changed by less than 1 MiB in each Python case
(including fixture creation). The peak R-child RSS increased from about 219 MiB
to 229 MiB: the single rendering process now retains the validated `ggimage`
namespace. Parsed plot inputs are retained for that family and released when its
R process exits; larger alignments need their own memory measurement. The JSON
records `RUSAGE_SELF` and `RUSAGE_CHILDREN` separately; these are not the sum of
simultaneously live process memory, and child measurements can include fixture
conversion. These small synthetic cases do not establish production HPC speedups
or global BH-FDR memory behavior.

## Reproduce

Run from the checkout root with the same fresh GeneGalleon runtime for both
revisions. Keep the private directory untracked and do not commit its results.

```bash
gg_perf_dir="$(mktemp -d "$PWD/.gg-performance-XXXXXX")"
git archive b60b10fed483e6306a583984b5a9968e20646cbc workflow/support |
  tar -x -C "$gg_perf_dir"
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_remaining_io.py \
  --support-root "$gg_perf_dir/workflow/support" --output "$gg_perf_dir/before.json"
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_remaining_io.py \
  --output "$gg_perf_dir/after.json"
```

Use `--case` to select one of the table cases and `--repeats` for more trials.
Compare every `output_sha256` before interpreting timing differences. The helper
builds only private synthetic data; it neither reads project data nor submits jobs.
