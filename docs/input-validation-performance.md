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
