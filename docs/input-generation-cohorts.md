# Aggregate verified formatted input cohorts

`workflow/support/input_generation_cohort.py` combines selected successful native
workers from multiple input-generation workspaces. It preserves compressed CDS,
GFF and genome bytes, sequence IDs, transcript ORF audits, experimental projection
proofs and cloned-model history. It never reruns the formatter or collapses
experimental sequence groups. GFF and genome remain optional for species whose
reviewed source has no such output.

The selection TSV has `species`, `source_workspace`, `plan_path`, `plan_sha256`,
`cds_path`, `gff_path`, `genome_path` and optional `cloned_model_audit_path`.
Each selected species must have exactly one matching native plan task and a
completion receipt binding its summary shard and all selected outputs/evidence.
A published multispecies summary is not required for unrelated donor species.

```bash
python workflow/support/input_generation_cohort.py plan \
  --selection selected_sources.tsv --expected-species 383 --plan cohort_plan.json
python workflow/support/input_generation_cohort.py prepare \
  --plan cohort_plan.json --workspace /data/project/cohort/workspace
```

Prepare takes native donor locks and stages an unpublished regular-file view
using hardlinks where supported and bounded streaming copies across filesystems.
After staging, it freshly verifies every receipt-bound source, linked/copied
output, summary and cloned audit snapshot before publishing `cohort_view.json`.
Each physical file is hashed once within that boundary. Failed copies, unexpected
species files and plan/source changes cannot publish a successful view.
Donor metadata and evidence remain required: retain those workspaces and treat
all linked inputs as read-only. No stat-only or persisted checksum establishes
success. A partially prepared identical view can resume; conflicting files or
changed selection/proofs require investigation.

The view provides standard `input/species_cds`, `input/species_gff` and
`input/species_genome` directories, native BUSCO/fx2tab artifacts under
`output/input_generation`, the selected original summary rows, byte-identical
cloned audit snapshots and a `cohort_view.json` provenance index. Original summary
paths remain in their donor namespaces so source-gene and ORF checks retain their
meaning. BUSCO artifact contracts remain usable against byte-identical files.

Run `qc --task-index N --plan cohort_plan.json --workspace ...` once per species,
using the standard longest-CDS/source-ownership and strict CDS/GFF/reference
validators. When saved rows contain `/workspace` paths, mount that exact donor
workspace at `/workspace` inside the same qualified container. Other project
paths must remain visible at their original absolute locations. An absent or
wrong donor mount is an error; ORF audits are never rewritten to waive checks.

For a relocated raw source, `plan --source-aliases aliases.json` accepts an
explicit JSON map from original absolute input paths to physical replacements.
Only receipt-bound raw CDS/GFF/genome inputs qualify, and their SHA-256 must
match the original receipt. Original summaries and receipt identities remain
in the plan alongside the relocation history; QC uses the recorded current
paths. Missing outputs, evidence, changed bytes and unused aliases are errors.

Cohort QC independently rederives historical transcript ORFs while retaining
their original audit bytes and IDs. Source/output hashes, CDSKIT implementation,
genetic code and selection policy must agree. The original header grouping is
reproduced for audits predating declared Trinity gene grouping. Every rederived
ID and sequence must match; the result records both implementation identities
and the grouping policy, without certifying biological gene uniqueness.
Ordinary formatter reuse and validation remain strict against the current
contract; this historical path only establishes independent rederivation.

Run `finalize --plan ... --workspace ...` only after every QC worker succeeds.
It verifies the exact cohort/view and each species result, writes `cohort_qc.json`
and a BUSCO `marker_occupancy.tsv` (S/D/F/M, `?` for ambiguous Complete rows,
`.` for absent marker evidence). Marker inventories/lineages must agree; taxa
with no single-copy markers are reported without automatic exclusion. Scheduler
acceptance and a prepared view do not certify final QC completion.

## Verification cost

Fresh input and view verification share one invocation-local `FreshDigestBatch`.
Each physical file is hashed once even when several view paths hardlink to it.
The stat/target fences remain live through read-only native QC and are rechecked
before publication; another invocation always hashes the content again. Native
source changes during QC prevent a success receipt.

The bounded verification benchmark uses two 64 MiB FASTA payloads, three
hardlinked names per payload, one warmup and three repetitions. Docker
Linux/arm64 (`0e46be37b7dd`), comparing the initial cohort implementation with
this change, produced identical returned documents. Median verification time
was 0.176 s → 0.076 s; hashed bytes were 384 MiB → 128 MiB. Peak Python
allocations stayed about 2.03 MiB. This measures content verification only,
not biological QC or total workflow speed.

To reproduce with a saved baseline module visible inside the container:

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_input_cohort_verification.py \
  --baseline /data/benchmark/baseline_input_generation_cohort.py
```

Preparation uses the same final content boundary after link creation. Compared
with v0.8.394, the same warmup/three-repetition Docker fixture produced identical
view manifests with median preparation time 0.133 s → 0.081 s and hashed bytes
256 MiB → 128 MiB. Peak Python allocations stayed about 2.04 MiB. Copied outputs
are now explicitly hashed before publication as well; distinct physical copies
still require their own content read. This measures assembly only, excluding
biological QC. Reproduce with the same command plus `--phase prepare` and a saved
v0.8.394 baseline module.

For a code-fix retry, `qc` and `finalize` accept `--qc-directory` beneath the
same cohort workspace. Use the same directory for every member and its finalizer;
previous QC records stay intact. This does not skip validation or rebuild inputs.
The finalizer still checks every selected species and freshly verifies all inputs.

Coding-GFF validation streams the filtered rows directly to the table reader,
retaining alignment-ancestor closure and historical row order. It no longer
writes an expanded GFF to temporary storage. In a Linux/arm64 Docker fixture
with 250,000 CDS rows, three alternating runs after a warmup produced identical
rows, dtypes and checksums: expanded temporary writes fell from 14,314,862 bytes
to zero, median peak RSS from 205,636 to 182,692 KiB. Median parsing times were
0.566 and 0.537 s; this small fixture does not establish whole-workflow speedup.
Reproduce with `workflow/benchmarks/benchmark_coding_gff_stream.py --baseline`
pointing to the previous `gff2genestat.py`, inside the GeneGalleon runtime.
