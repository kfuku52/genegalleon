# Genome prediction with ANNEVO

`gg_genome_annotation_entrypoint.sh` defaults to [OrionGeno](oriongeno-genome-annotation.md).
Set `genome_annotation_backend="none"` for existing CDS-based behavior, or
`genome_annotation_backend="annevo"` to predict from the FASTA files
in `workspace/input/species_genome`. No curated GFF or CDS is required for that
route. Each array task uses the sorted genome file list; set the array size to
that list's size. Scheduler resource directives are not changed automatically.

The selected backend prepares a dedicated cached ANNEVO Python environment,
source on first use. Only the selected lineage model is downloaded and verified;
all species share the environment and immutable model cache. It runs within the main GG
container, with native ANNEVO and its dependencies isolated from GG's base
Python environment. With `genome_annotation_backend="none"`, no ANNEVO cache
lookup, download or GPU request occurs.

The qualified conda-pack environment has 255-byte binary prefix fields. Its final
installation prefix below `workspace/downloads/annotation/annevo/runtimes/` must
fit that limit, including the cache key, generation and `environment` directory.
Choose a shorter workspace path before first installation when needed; installed
generations cannot be renamed after relocation. Failed installation records and
generations remain available for diagnosis. Remove this path constraint only
after qualifying a payload whose binary relocation fields accept longer prefixes.

## Configure the stage

Edit the existing entrypoint's editable block, or use its registered overrides:

```bash
GG_GENOME_ANNOTATION_GENOME_ANNOTATION_BACKEND=annevo \
  bash workflow/gg_genome_annotation_entrypoint.sh
```

Run this inside an allocation with a CUDA GPU. The scheduler must supply
`CUDA_VISIBLE_DEVICES`, or an explicit GPU selection must be supplied before
running locally. The ANNEVO path requires actual CUDA execution; it does not
fall back to CPU prediction. Apptainer/Singularity receives `--nv`. The Docker
adapter translates this to `--gpus device=<selection>`, accepting only GPU
indices or UUIDs. A different explicit Docker selection or forwarded CUDA
visibility is rejected. `GG_CONTAINER_DOCKER_GPU_DEVICES` can provide an explicit
selection when no scheduler selection exists; it cannot override a different
`CUDA_VISIBLE_DEVICES` allocation. The adapter never requests all GPUs
implicitly. Docker receives the physical selection and maps CUDA visibility
inside the container to `0..N-1`; for example, host `2,4` becomes container
`0,1`. SIF execution retains the scheduler's original CUDA selection.

| Entrypoint option | Default | Purpose |
| --- | --- | --- |
| `genome_annotation_backend` | `oriongeno` | `none` retains CDS analysis; `annevo` and `oriongeno` enable native genome prediction. |
| `annevo_model_profile` | `Magnoliopsida` | Explicit native profile, or `auto` for per-species NCBI ancestor selection. |
| `annevo_taxid_map` | blank | Explicit TSV with `species` and positive `taxid` columns. |
| `annevo_genome_metadata` | blank | Trusted TSV with `species`, `taxid`, `genome_sha256`, and optional `scientific_name`. |
| `annevo_taxonomy_db` | blank | Read-only ETE NCBI SQLite snapshot; blank uses the existing workspace taxonomy cache. |
| `annevo_runtime_manifest` | blank | Public runtime registry on first use, or an explicit qualified manifest locator. |
| `annevo_no_download` | `0` | `1` requires an existing verified runtime and model cache. |
| `annevo_batch_size` | `32` | Native GPU prediction batch size. |
| `annevo_sequence_threshold` | `25` | Native sequence-group threshold; whole contigs remain intact. |

The other overrides follow the same `GG_GENOME_ANNOTATION_` prefix, for example
`GG_GENOME_ANNOTATION_ANNEVO_NO_DOWNLOAD=1`. A verified `representative_inputs`
view cannot be combined with new ANNEVO prediction. Native windows, flank sizes,
minimum intron/protein criteria and decoding remain those of the qualified
runtime. The qualified ANNEVO annotation/candidate pipeline supports genetic
code 1; other `genetic_code` settings fail before taxonomy or runtime acquisition
and GPU prediction. The existing CDS-based backend retains its genetic-code
support. BUSCO, annotation and contamination settings retain their meaning.

## Select a model automatically

Set `annevo_model_profile="auto"` and supply `annevo_taxid_map` or
`annevo_genome_metadata`. The `species` column must exactly match the selected
workflow species key. One map can contain all species in a mixed array:

```tsv
species	taxid
Homo_sapiens	9606
Gallus_gallus	9031
Danio_rerio	7955
Drosophila_melanogaster	7227
Saccharomyces_cerevisiae	4932
Oryza_sativa	4530
```

Genome metadata must come from a reviewed assembly/provider record and bind its
TaxID to the SHA-256 of the actual input FASTA bytes (compressed bytes for `.gz`).
A checksum binds the record to an input; it does not establish the biological
accuracy of a provider's TaxID. GeneGalleon does not promote the older formatter's
name-derived TaxID summary to authoritative genome metadata. Species labels,
qualified names, assembly descriptions and nearest-clade guesses are never used
to fill missing TaxIDs. Explicit maps are user-curated authoritative evidence.
If both sources provide a TaxID, they must agree after NCBI merged-ID resolution.
An optional scientific name must have exactly one NCBI identity in the genome's
ancestor chain. Duplicate species rows, bad hashes, conflicting TaxIDs/names,
clade-only TaxIDs, incomplete ancestor chains and unresolved/deleted TaxIDs stop
prediction with a reason in `model_selection.json.failed.json`.

The existing shared ETE taxonomy cache is acquired only when `auto` needs it and
downloads are enabled. A cached database or explicit snapshot is reused without
refreshing it. `annevo_no_download=1` requires the taxonomy database, runtime and
selected model to be present. To adopt newer taxonomy, replace the snapshot
explicitly and use the normal stale-artifact policy. Active SQLite transactions
are rejected before and after selection and at the prediction boundaries; use a
quiescent snapshot. Selection reads the hashed SQLite main file in one read
transaction, excluding WAL changes, and rejects changes before publishing its
receipt. CLI outputs cannot replace the genome, TaxID evidence, or taxonomy DB.

| NCBI ancestor | TaxID | Native model and `-l` | Overlapping prediction |
| --- | --- | --- | --- |
| Mammalia | 40674 | Mammalia | enabled |
| Aves | 8782 | Aves | native default |
| Actinopteri | 186623 | Actinopteri | enabled |
| Insecta | 50557 | Insecta | native default |
| Fungi | 4751 | Fungi | native default |
| Magnoliopsida | 3398 | Magnoliopsida | native default |

Magnoliopsida means NCBI flowering plants, including monocots and dicots.
Gymnosperms, ferns, non-insect arthropods and other unsupported lineages receive
an explicit error; no alternative model is substituted. ANNEVO itself supplies
window/flank sizes for the selected `-l`; the execution receipt records the
actual native configuration and prediction source hash. Mammalia and Actinopteri
use `--overlap_pred` as recommended by the [official ANNEVO source](https://github.com/xjtu-omics/ANNEVO).
All six names remain available as explicit selections. The default stays
`Magnoliopsida`, preserving existing explicit workflows and receipts.

## Outputs and restart

Generated data stays under
`workspace/output/genome_prediction/annevo/<species>/`:

- the untouched native `<species>.gff3`, native execution receipt and prediction
  provenance;
- `formatted/<content-key>/species_cds`, `species_gff`, `species_genome`, plus
  standard formatter summaries and source provenance;
- a separate formatter provenance record;
- for `auto`, `model_selection.json` and a species-local `runtime.json`. The
  selection record includes input/evidence/database hashes, original and merged
  TaxIDs, the full ancestor chain, matched model anchor, reason and overlap mode.
  The native receipt embeds the selection and records model/source/runtime hashes,
  argv, GPU proof and inference conditions. Failed resolution retains the previous
  complete selection and prediction.
- `candidate_inputs.tsv` and its provenance for optional
  [rescue/refinement additions and corrections](annevo-model-improvement.md).

The standard input formatter extracts CDS from the native GFF and genome,
retaining its existing gene grouping, phase, identifier and longest-CDS rules.
The stage passes those generated files to the existing CDS resolution and enabled
annotation/QC steps. It does not copy predictions into curated `workspace/input`
files or declare BUSCO complete from prediction exit alone.

A verified unchanged prediction is reused, including after interruption between
prediction and formatting. The standard `artifact_stale_policy=stop` rejects a
changed genome, runtime/model receipt or native parameters. Automatic
selection also binds taxonomy/evidence/selector hashes and the resolved profile.
An unchanged offline rerun preserves selection bytes and reuses prediction;
changes to taxonomy, evidence, inputs or models invalidate it. To adopt newer
model/runtime metadata, use an explicit new manifest locator or the existing
`annevo_runtime.py prepare --refresh` operation before the next stage run;
ordinary reruns retain the first-use snapshot for offline reproducibility. Species locks
serialize the same species; auto selection and runtime receipts are per species,
while the existing environment and model locks serialize shared first use. An
explicitly selected `rebuild` regenerates the affected stage; previous complete formatted
content generations remain available. Failed native prediction preserves its
previous complete GFF and native receipt. Cache/model availability alone does
not satisfy prediction or downstream completion checks.

ANNEVO's component uses a non-commercial license; retain its license and notices
with the cached source/model. See [runtime bundles](annevo-runtime-bundles.md)
for source, model, dependency and distribution provenance. GG's base environment
and its own license remain independently recorded.

For the opt-in six-clade GPU test and its bounded input manifest, see
[development validation](development-and-tests.md#optional-annevo-clade-gpu-qualification).

## Pipeline several genomes in one GPU allocation

For rescue/refinement over many genomes, run `workflow/support/annevo_pipeline.py`
inside the GPU-enabled GG container. Prepare the optional runtime/model once with
`annevo_runtime.py prepare`, then create a TSV with `species`, `genome`, `receipt`
columns. `receipt` names the installed runtime receipt, not a prediction receipt.
An optional fourth `selection_receipt` column retains an automatic model selection.
Paths can be absolute or relative to the TSV. The program invokes native
`prediction.py` on one GPU and native `decoding.py` in two independent CPU slots:

```bash
python workflow/support/annevo_pipeline.py --inputs genomes.tsv --output predictions \
  --predict-threads 4 --decoder-workers 2 --decode-threads 4 --max-tmp-gib 128
```

Reserve at least 14 CPUs for this configuration: four predictor threads, two
data-loader workers and two four-thread decoders. One GPU allocation is sufficient;
independent Slurm jobs requesting an exclusive GPU cannot share that allocation.
The ordinary annotation entrypoint and scientific defaults are unchanged. Completed
species publish native GFF/receipt pairs; after all species succeed, the combined
`candidate_inputs.tsv` feeds the existing rescue/refinement `annevo_inputs` option.
Unchanged pairs are verified and reused. Changed inputs/parameters or corrupt pairs
stop execution; select a new output directory for a new prediction generation.
Failures retain logs and never publish an incomplete campaign table.

Complete contigs and native thresholds are preserved. Temporary space is reserved
before prediction using 60 bytes/base for native probabilities plus FASTA/metadata
headroom. The program bounds pending work and cleans each intermediate after its
decoder finishes. A genome exceeding the budget stops with an explicit error;
the sequence-group threshold does not cap the complete HDF5 file size.

`workflow/benchmarks/benchmark_annevo_pipeline.py` compares two alternating pairs
of serial/native annotation and pipelined execution on exactly the same TSV,
runtime, model and GPU. It records wall time, GPU samples, native receipts and GFF
feature/byte equality. It runs workloads only inside an existing allocation;
resource planning and scheduler submission remain separate. Rescue admission and
representative selection require a separate comparison using the frozen donors
and actual production gates.

The whole-genome admission benchmark accepts `--compare-annevo-gff` and
`--compare-receipt` alongside its original native prediction. It nominates and
validates both actual predictions against the same miniprot placements, donor
catalogs and masks, then checks individual accepted cases, representative choices
and selected GFF/CDS bytes. Independent execution-receipt and raw GFF artifact
identities are retained in the report. Native parallel decoding can reorder
sequence-number comments; these identities are excluded from decision equality
only after every native annotation field and model/genome input agrees.
Indexing and donor placement finish before full annotation catalogs are retained,
so their independent memory peaks do not overlap in bounded validation runtimes.
Its placement graph remains explicitly
approximate; a fixed alternative annotation is a regression control, not
independent biological truth.

The 2026-10-14 whole-genome comparison used Aethionema arabicum (203 Mb),
Brachypodium stacei (234 Mb) and Geranium dissectum (568 Mb), one RTX 6000 Ada,
16 allocated CPUs and 64 GiB RAM. Both routes used the same installed runtime,
Magnoliopsida model, batch size 32 and four decoding threads per species.

| Matched trial | Serial | Pipeline |
| --- | --- | --- |
| 1, serial first | 32.53 min | 27.50 min |
| 2, pipeline first | 29.86 min | 26.15 min |
| Median | 31.19 min | 26.83 min |

The measured reduction was 14.0% (1.16x speedup), with a sampled peak GPU memory
of 5.33 GiB and peak intermediate reservation of 50.8 GiB. All nine GFF annotation
fields matched across all four runs. Whole-gene and internal-exon masking plus
unmasked controls retained identical individual admissions, representative choices
and selected GFF/CDS. These admission assays used an approximate placement graph;
Geranium's fixed Helixer annotation was a regression control. In this input order,
GPU prediction overlapped CPU decoding, but only one of the two decoder slots was
active at a time. The large final genome and probability-file I/O limited the
gain. This measurement does not establish steady-state throughput for 500 species.

Additional pipeline options preserve the original configuration when omitted:

| Option | Default | Purpose |
| --- | --- | --- |
| --order | input | largest-first schedules longer assemblies earlier; candidate TSV rows retain the original input order. |
| --compression | none | lzf passes native --comp for lossless probability storage; storage mode is recorded in prediction parameters. |
| --jit-cache-root | omitted | A private shared cache for native Numba decoding; omitted retains per-species temporary caches. |

For a comparison, use a new output directory rather than reusing predictions
made with different storage parameters. Scheduling verifies input hashes before
reusing assembly lengths. Shared JIT caches bind the installed runtime, source
generation, CPU identity and Numba CPU overrides. Native Numba additionally
checks compiler/CPU compatibility. Cache roots must be private, owned by the
calling user and separate from inputs and the immutable installed environment.
Select a local filesystem that preserves private POSIX permissions; output
shares that create permissive directories despite a `0700` request are rejected.
Changed identities select a new cache; unsafe members or modified identity
records stop execution. No native program, model, prediction window, probability
precision or admission threshold is changed by these options.

For isolated effects and resource measurements, run
workflow/benchmarks/benchmark_annevo_pipeline_optimizations.py in an existing GPU
allocation with --inputs, a fresh --output, and --probability-proof-genome naming
a qualified small matched-profile FASTA. It first verifies actual native plain/LZF
probability bits, shapes, dtypes and metadata on that fixture outside timing.
Two untimed native GPU/decode probes then verify cold/warm JIT reuse, GFF
equivalence and resource sampling before full-genome trials. The harness keeps
compiled caches under `--tmp-root`, separate from the GFF output share. Preflight
and timed caches are separate, so the first optimized trial starts cold; timed
trials share their cache until the harness exits and removes it.
Four cumulative configurations run twice in reverse order: original settings,
ordering only, ordering plus LZF, and ordering plus LZF plus shared JIT.
Reports retain per-species phases, cache load/save events, GPU samples and
one-second process-tree RSS/PSS. Sampled peaks are lower bounds and do not include
whole-job page-cache/cgroup accounting. Large-genome output equivalence and actual
rescue/refinement admission still require their recorded comparisons.
