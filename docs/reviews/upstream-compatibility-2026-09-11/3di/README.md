# 3Di implementation validation — 2026-09-11

The implementation connects GeneGalleon's search, scan, normal sites and scan
candidate sites to one full-CDS family bundle. The existing statistical options
are preserved. This initial validation used local CSUBST changes. Subsequent
published revisions and clean-image checks are recorded in the
[publication validation](../publication/README.md).

## Runtime

Validation used a Linux ARM64 GeneGalleon Docker runtime based on
`local/genegalleon:nwkit-shift-covariance-dev`
(`sha256:c8c69075b5d591499bf3195ec912ef5b5e1aa1cd163f03864639ed6a669e9429`).
The current GeneGalleon workflow was mounted, the current treevis source was
installed, and modified CSUBST was built and reinstalled as a compiled wheel.
The other owned Python tools were subsequently rebuilt from their current
upstream source heads for the final integration run:

| Program | Source used |
| --- | --- |
| GeneGalleon | `10fc03bd9c40c51cc635c5254e6b099c5ad75323` plus local changes |
| CSUBST | `cde05f34e5ea51c3f0f1bae721ccd45243e46f38` plus local changes |
| CDSKIT | `f7ee05ea43451929da1db4d32e4f6cb5c02719ec`, version 0.31.1 |
| NWKIT | `1996588e561d10aebc47c4daf90be836bbebd315`, version 0.43.16 |
| AMALGKIT | `21ce1417410ec05aea409370d496ce01896671d0`, version 0.16.88 |

The runtime has Python 3.12.14, NumPy 1.26.4, SciPy 1.17.1, pandas 3.0.5,
PyTorch 2.13.0, Transformers 5.17.0 and PEFT 0.20.0. The 3Di dependency/API check
and `pip check` passed. Backend settings came from CSUBST: ESM3Di-35m, direct
ASR and the uniform GTR structural model.

This was a source-wheel update of an existing GeneGalleon Docker image, not a
clean complete multi-platform image build. Native Apptainer definition/staging
was tested, but Apptainer/Singularity is unavailable on this macOS host, so SIF
execution and AMD64 execution were not tested in this initial run. The default
remote-source build required publication of the CSUBST changes before this
feature could pass its runtime check.
No dependency-side fallback was added to GeneGalleon.

## Owning CSUBST fixes

- Native 3Di joint scan now uses the independent structural generator,
  stationary frequencies and branch lengths. Codon synonymous counts retain
  their own codon fit and rate mixture. Unsupported native analytical,
  parametric and bridge paths are still rejected.
- Joint search acquires both expected reducers before early state release, so
  computing S does not fingerprint already-freed 3Di states.
- 3Di state-cache matching retains content hashes while allowing identical
  files to move. IQ-TREE is resolved through PATH before its executable path
  is recorded, and genetic code is part of the cache context.

The scan test compares native observed tensors with the separate joint
endpoint implementation, exposure with a direct matrix exponential, and
synonymous counts with independent codon inference. Early release produces
the same expected-count table as retained states. Cache tests cover relocation,
changed content, changed genetic code and PATH resolution.

## Checks

| Check | Result |
| --- | --- |
| CSUBST affected unit/integration tests with compiled extensions | 185 passed |
| CSUBST fast unit/CLI lane | 1,821 passed; 2 optional Gemmi reference tests skipped; 7 deselected by the lane |
| CSUBST lint, repository hygiene, documentation and configured mypy targets | Passed |
| GeneGalleon input/command/report/cache/build-staging regression tests | 232 passed |
| Real 3Di integration with current owned tools | 2 cases passed with workflow defaults (genetic codes 1 and 2) |
| GeneGalleon static lane | 397 passed |
| GeneGalleon Ruff and configured ShellCheck on changed shell implementation | Passed |
| Current-source treevis R tests | Passed |
| 3Di imports/API contract and installed dependency consistency | Passed |

The separate real-predictor suite uses eight retained tips, 100 full codons and
90 trimmed codons. It builds the actual core ASR fragment and full bundle,
archives and relocates that bundle, then executes the actual core search/scan
commands. Direct CSUBST calls produce matching search/scan tables. Both genetic
codes 1 and 2 are exercised. The code-2 input also includes a codon gap, an
ambiguous codon and an excluded extra tip. Invalid frames and duplicate IDs
are rejected in the small input tests.

The final suite uses the workflow's existing codon defaults, `ECMK07+F+R4`
for code 1 and `GY+F+R4` for code 2. The mitochondrial fixture converts
AGA/AGG to sense codons and includes both TGA and TGG tryptophan codons.
During fixture development, IQ-TREE 3.1.4 emitted all-NaN ancestral
probabilities with `GY+F+R4` when one tryptophan codon was unobserved,
including after internal stops were removed. It reported numerical underflow
and exited successfully; CSUBST correctly rejected the invalid probabilities.
The [minimal reproduction](iqtree-nan-reproducer/README.md) preserves this
IQ-TREE numerical limitation, which remains unresolved. No pseudocount,
statistical-model substitution or relaxed probability validation was added.

Normal sites generates its combined report through `process_index`; scan
candidate sites generates its focused report through `analyze_candidate`.
Synthetic sequence IDs have no public structure accession, so the test disables
PDB lookup while retaining the real model inference and R/PDF generation.
The reports were rendered and their full-CDS and structural site panels
inspected. Site-label wrapping was tightened for adjacent single-site panels.

For the code-1 fixture, raw scan site 75 is public full-CDS codon 76, corresponding
to trimmed codon 66. Branches 0 (`a`) and 7 (`e`) both have structural state `D`
at that site, while their amino acids are `E` and `H`. This confirms that the
structural panel uses the predicted states and full positions rather than an
amino-acid recoding or the same numeric position in the trimmed alignment.
These synthetic fixtures validate integration, not biological significance.

## Resource reuse and CPU measurement

Two simultaneous first downloads into a new shared resource directory completed.
One process acquired the CSUBST download lock; the other waited and reused the
validated files. `csubst download --no_download yes` then succeeded. The
integration suite runs subsequent analyses with Hugging Face/Transformers
network access disabled and a required state cache.

Three fresh-directory CPU runs used the eight-tip, 100-codon fixture, two
threads, predownloaded model weights, `--sa_cache no`, `--sa_state_cache no`
and `--sa_no_download yes`. Each `inspect` run includes fresh codon/structural
fits and prediction. Wall times were 13.71, 14.91 and 13.65 seconds; peak child
RSS was 841,276–841,796 KiB (about 822 MiB). These are short-input measurements,
with no before/after speed claim or extrapolation to long proteins or GPUs.

The exact command and measurements are in `cpu-measurements.json`; `measure.py`
reproduces each run with repository mounted at `/review`, writable outputs at
`/audit`, resources at `/audit/downloads-parallel`, and `REPLICATE=1`, `2` or `3`.
