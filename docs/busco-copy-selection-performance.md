# BUSCO copy-selection performance

This comparison measures the raw `best_score` batch extraction helper, including
candidate validation, FASTA extraction, and both TSV reports. It excludes BUSCO
discovery, workflow finalization/compression, alignment, and tree inference.

## Recorded comparison

On 2026-10-22, Docker Linux/arm64, Python 3.12.15, Biopython 1.88 and pandas 3.0.6
used image `sha256:86245ca742d3aff03d4e5fc53c9f0ca116ec55b51f046f92a29f7b4182fd03e2`.
Each variant had one warmup and three measured trials, with alternating execution
order and a fresh worker process per trial. The timings below are medians without
profiling; peak RSS is the median worker peak. No other builds or tests ran during
the comparison. Fixtures and outputs lived on the container's local temporary
filesystem, with gzip input FASTAs.

| Workload | Species × markers | Duplicated fraction / copies | Input sequences | Before → after seconds | Speedup | Before → after peak MiB |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| Protein, 300 residues | 192 × 512 | 25% / 4 | 172,032 | 5.559 → 2.068 | 2.69× | 184.5 → 131.7 |
| CDS, 900 nucleotides | 64 × 512 | 25% / 3 | 49,152 | 2.884 → 1.684 | 1.71× | 107.6 → 55.8 |
| Many-copy stress case, 300 residues | 8 × 128 | 100% / 64 | 65,536 | 2.245 → 0.642 | 3.50× | 135.8 → 85.9 |

Every FASTA and both TSV reports matched byte-for-byte between variants in all
warmup, measured and separately profiled trials. The 64-copy case exercises a
scaling limit rather than a typical biological copy-number distribution.
Synthetic sequences have low complexity; these results do not establish
production NAS throughput, SIF performance, or a whole-workflow speedup.

The [recorded JSON](benchmarks/busco-copy-selection-2026-10-22.json) includes raw
trials, source hashes, package versions, output fingerprints and profiles.
The baseline was the implementation saved immediately before this review;
the local source snapshot remains in
`workspace/output/benchmarks/busco-copy-selection-review-20261022/baseline-support`.

## Measured bottlenecks and changes

The protein profile repeatedly resolved the same output paths when reopening
files through the bounded handle pool. Validated paths are now cached per marker;
every reopen still uses `O_NOFOLLOW`. Candidate rows are indexed by species,
BUSCO and normalized source gene, removing repeated scans of unrelated species
and quadratic searches within many-copy groups. A lightweight shared identifier
module avoids importing pandas for copy selection, and string translation counts
ambiguous residues in native code with the same residue definition.

The ranking rule and scientific thresholds are unchanged. Stronger input/output
validation is included in the after measurements. Remaining CDS time is dominated
by Biopython translation, whose genetic-code behavior is retained. Output writing
and audit serialization remain significant; candidate audit metadata is retained
in memory, while FASTAs are scanned once and only current winning duplicate
sequences are buffered per species.

## Reproduce a comparison

Before changing the implementation, save its support directory. The comparison
requires that saved code; copying the already optimized code gives a self-comparison.

```bash
baseline_support=$(mktemp -d /tmp/gg-busco-baseline.XXXXXX)
cp -a workflow/support/. "${baseline_support}/"
```

After making a change, use the normal runtime wrapper and run without concurrent
builds or tests:

```bash
GENEGALLEON_DOCKER_EXTRA_BINDS="${baseline_support}" \
  bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_busco_copy_selection.py \
  --baseline-support "${baseline_support}" --repeats 3 --profile \
  --output tmp/busco-copy-comparison.json
```

The harness creates deterministic fixtures in its own temporary directory, runs
variants against the same output paths, and fails if any output fingerprint
differs. `--baseline-only` records an initial measurement, and `--cases protein`,
`--cases cds` or `--cases many-copies` selects a workload. A usable GeneGalleon
runtime and a baseline with the requested `best_score` CLI are required.
