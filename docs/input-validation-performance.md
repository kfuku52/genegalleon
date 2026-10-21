# Input validation performance

CDS reconstruction indexes FASTA reference lengths once per opened genome.
Every fragment still checks reference existence, nonnegative start, positive
span and end bounds before fetching bases. Strand, block ordering, overlap
rejection, sequence content and genomic reconstruction rules are unchanged.

Array finalization can restore worker mapping/selection QC instead of repeating
the same validators. Reuse requires the current format and validation contracts,
matching options, support-code/runtime identity and source summary, complete format/validation checkpoints,
and QC included in the worker completion receipt. One fresh hashing batch per
species in the combined native final-QC invocation covers the plan, settings, worker metadata, all receipt
sources/outputs, checkpoints and QC. File identities and symlink targets are
fenced across both validators and checked again before accepting combined QC;
persistent metadata-only digest reuse is not used. Independent public validator
CLIs retain their separate verification boundaries. Implementation identity is
computed once in the combined process, with its full existing dependency scope.
Species identity, original/formatted genomes and source-ownership proof are
required. Failed, absent, legacy or incompatible proofs run ordinary validation.
The validators retain per-species results and diagnostics in additive
`species_results`, `validation_options` and `validation_implementation` fields
in their JSON QC outputs. Implementation identity covers support Python sources,
installed sequence/GFF parser source and compiled libraries, and runtime versions.
Older JSON remains readable but cannot supply these cached diagnostics.

Array completion verification, complete species-set checks, CDS format checks,
and publication after all required stages still execute. Single-mode validation
is unchanged. Direct validator calls can opt in with all three arguments:

```bash
--reuse-validation-root /data/workspace/output/input_generation \
--reuse-task-plan /data/workspace/output/input_generation/tmp/task_plan.json \
--format-contract-version FORMAT_CONTRACT_FROM_CURRENT_CORE
```

## Bounded measurements

Measurements on 2026-10-05 used Docker Linux/arm64, Python 3.12.14,
pysam 0.23.3 and immutable runtime
`sha256:8a6d69ec10defff91117d59591eeaabdbd76b41218525a7c038de815af54f87e`.
Each comparison alternated the implementations in fresh processes, with one
warmup and three measured trials. Baseline source was committed `0dbb271`
(version 0.8.135), including the separately developed staging/resume improvements.

| Workload | Before median | After median | Equivalence |
|---|---:|---:|---|
| 89,579 references; 1,000 CDS reconstructions with four blocks and both strands | 1.437 s | 0.222 s | Identical reconstructed-CDS SHA256 |
| Final mapping/selection QC; two native worker fixtures, each with 64 MiB genome | 2.789 s | 0.399 s | Identical aggregate QC and source ownership |

Reference indexing increased process peak RSS from about 52.1 to 57.2 MiB
(roughly 5 MiB). Final-QC trials both peaked at about 71.4 MiB. These are
bounded synthetic comparisons, not full-finalize or production NAS speedups.
The reference benchmark includes first FASTA index/open and fetches but excludes
fixture generation and GFF parsing. Final-QC fixtures use the native workers
and real validators/GFF reader, with fake BUSCO/seqkit and local taxonomy during
fixture setup only; setup and imports are excluded from measured validation.
Implementation identity is recomputed for each validator, as in the core's
separate validator processes. Native performance instrumentation is enabled
in the final-QC benchmark.

Final-QC reuse read about 315 MB for fresh SHA256 checks in these raw-FASTA
fixtures and avoided reference scans. It can increase hashing traffic versus
ordinary QC that does not otherwise consume every raw source. NAS throughput,
compression, input size and reconstruction needs can change that tradeoff;
measure real workloads before predicting total elapsed time. Checksum byte
counts exclude ordinary parser reads and are not total I/O estimates.

## Reproduction

Keep a baseline support tree and normaliser file visible to the same runtime.
Run without concurrent builds, benchmarks or tests. Output paths must be new.

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_cds_reference_lookup.py \
  --baseline-source /data/baseline/workflow/support/cds_model_normalisation.py \
  --output tmp/reference-comparison.json

bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_final_validation.py \
  --baseline-support /data/baseline/workflow/support \
  --output tmp/final-validation-comparison
```

Use the runtime wrapper's extra-bind mechanism for the baseline. Results retain
all measured timings, Linux process peak RSS, workload, environment and logical
output hashes. The final-QC benchmark also retains the native fixture workspace
and proofs. Contract mismatch, missing diagnostic/genome proof and changed
summary are covered by the input-generation end-to-end tests; ordinary validator
and source-ownership tests still exercise fresh validation.

## Sharing current worker-proof verification

The 2026-10-05 comparison against v0.8.139 (`1140bbd`) in the same immutable
Docker runtime used two worker fixtures with 64 MiB genomes each, one warmup
and three alternating fresh processes. Each variant generated its own native
proofs so its implementation identity matched. Both restored all four species
QC results; complete aggregate and `species_results` JSON matched after
excluding only the deliberately changed implementation identity.

The median final-QC time fell from 0.393 to 0.185 seconds, and SHA-256 bytes
from 314,523,308 to 157,299,949. This measures the two validators and proof reads,
excluding fixture creation and module imports; it is not a complete finalize,
NAS or real BUSCO measurement. Merge receipt verification remains independent.
Use `benchmark_native_efficiency.py --mode final` as documented in
[input-generation performance](input-generation-performance.md).

Tests retain ordinary fallback for absent/old QC, changed implementation/options,
summary and genome evidence, plus rejection of changes between the combined
validators. Source summary rows are parsed from exact bytes bound to the proof.

## Completion receipt aliases

Worker completion records can contain both a raw source's resolved path and its
original alias. Completion verification now hashes each physical file once per
receipt, including symlink, hard-link and explicitly mapped `/workspace` aliases.
Each declared path must still match its own expected SHA-256. The receipt and all
member identities/targets are checked again after reading the last output, so a
change to an earlier source cannot be accepted. No checksum is retained between
calls. Existing receipt schemas, retry selection and scientific validators are
unchanged; a supplied plan hash retains its existing caller-owned binding.

A 2026-10-16 comparison against v0.8.276 (`7f2fa685`) used Docker Linux/arm64,
Python 3.12.15 and runtime
`sha256:a28e1760f1cf2cca40798943bc8e48d594b160147c00baf0441820ada07e41a7`.
Eight synthetic receipts each covered a 64 MiB source and a 1 MiB output.
One warmup preceded five alternating fresh-process trials in a warm filesystem.

| Receipt workload | Before median | After median | SHA-256 bytes before / after |
|---|---:|---:|---:|
| Source plus one alias | 0.519 s | 0.274 s | 1,082,133,880 / 545,267,336 |
| Source without aliases | 0.273 s | 0.273 s | 545,262,968 / 545,266,304 |

The alias workload reduced verification time by 47% and hash bytes by 50%.
Median process peak RSS was 24.3 / 24.6 MiB, including fixture setup. The control
has effectively unchanged elapsed time and small additional receipt reads for
mutation detection. Both implementations accepted the same valid receipts,
rejected changed output content and left receipt bytes unchanged. These figures
measure receipt verification, not ORF prediction, whole-job time or NAS throughput.
Regression tests additionally reject changes to earlier files, receipt metadata,
plan content and alias targets during verification.

Keep the baseline support tree visible to the runtime and use a new output path:

```bash
GENEGALLEON_DOCKER_EXTRA_BINDS=/data/baseline \
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_native_efficiency.py \
  --mode receipt --baseline-repo /data/baseline \
  --species 8 --genome-mib 64 --trials 5 --output tmp/receipt-comparison.json
```

Repeat with `--receipt-aliases 0` and a different output path for the control.

## Genome sequence identity

OrionGeno, ANNEVO and native candidate admission compute a canonical identity
from contig names and DNA. Wrapped FASTA previously normalized, validated,
encoded and updated SHA-256 separately for every sequence line. The reader now
batches about 64 KiB of stripped sequence fragments before those operations.
Each contig still receives its own exact length and sequence hash; sorting,
header identity, case/whitespace semantics and the canonical JSON stay unchanged.
The parser still reads the complete stream, checks gzip CRC/truncation, rejects
invalid DNA and duplicate/empty contigs, and checks the open-file/path identities
before returning. No saved digest or metadata-only reuse is introduced.

On 2026-10-21, baseline v0.8.337 (`f7e1e559`) and the updated reader ran in the
same Docker Linux/arm64 runtime, Python 3.12.15:
`sha256:594a13182c692adfa393107910c67cf57591ff8127d7eaad7724f5ea552e83af`.
One warmup preceded three alternating fresh-process trials per implementation.
All canonical output hashes matched, including the plain/gzip identity.

| Synthetic workload, 60 columns | Wall seconds before / after | CPU seconds before / after | Peak RSS MiB before / after |
|---|---:|---:|---:|
| 128 MiB DNA, eight contigs | 1.083 / 0.527 | 1.082 / 0.527 | 19.92 / 19.92 |
| Same DNA, gzip | 1.165 / 0.643 | 1.165 / 0.643 | 20.13 / 20.13 |
| 8 MiB DNA, 4,096 contigs | 0.073 / 0.039 | 0.073 / 0.039 | 21.26 / 21.26 |

These warm-cache measurements cover parsing and identity computation, excluding
fixture generation, imports and inference. The repeated-pattern gzip fixture
compresses strongly; it does not model production disk/NAS throughput, native
GPU prediction time or whole-job speedup. A single unusually long FASTA line
can still determine peak memory; batching avoids accumulating entire contigs
across their wrapped lines.

Keep the baseline source visible to the runtime and choose a new output path:

```bash
GENEGALLEON_DOCKER_EXTRA_BINDS=/data/baseline \
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_genome_sequence_identity.py \
  --baseline-source /data/baseline/genome_sequence_identity.py \
  --mib 128 --trials 3 --output tmp/genome-identity-comparison.json
```

The report retains input/source hashes, workload, environment, every measured
wall/CPU/RSS sample and canonical output hashes. Regression tests cover Unicode
header and legacy case behavior, wrapped/gzip equivalence, invalid symbols
across batch boundaries, complete gzip trailers and late input replacement or
same-size edits with a restored mtime.
