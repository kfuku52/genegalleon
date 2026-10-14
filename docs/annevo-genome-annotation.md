# Genome prediction with ANNEVO

`gg_genome_annotation_entrypoint.sh` keeps its existing CDS-based behavior by
default. Set `genome_annotation_backend="annevo"` to predict from the FASTA files
in `workspace/input/species_genome`. No curated GFF or CDS is required for that
route. Each array task uses the sorted genome file list; set the array size to
that list's size. Scheduler resource directives are not changed automatically.

The selected backend prepares a dedicated cached ANNEVO Python environment,
source and Magnoliopsida model on first use. It runs within the main GG
container, with native ANNEVO and its dependencies isolated from GG's base
Python environment. With `genome_annotation_backend="none"`, no ANNEVO cache
lookup, download or GPU request occurs.

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
| `genome_annotation_backend` | `none` | `none` retains CDS analysis; `annevo` enables genome prediction. |
| `annevo_model_profile` | `Magnoliopsida` | Qualified native flowering-plant profile. |
| `annevo_runtime_manifest` | blank | Public runtime registry on first use, or an explicit qualified manifest locator. |
| `annevo_no_download` | `0` | `1` requires an existing verified runtime and model cache. |
| `annevo_batch_size` | `32` | Native GPU prediction batch size. |
| `annevo_sequence_threshold` | `25` | Native sequence-group threshold; whole contigs remain intact. |

The other overrides follow the same `GG_GENOME_ANNOTATION_` prefix, for example
`GG_GENOME_ANNOTATION_ANNEVO_NO_DOWNLOAD=1`. A verified `representative_inputs`
view cannot be combined with new ANNEVO prediction. Native windows, flank sizes,
minimum intron/protein criteria and decoding remain those of the qualified
runtime. Existing genetic-code, BUSCO, annotation and contamination settings
retain their meaning.

## Outputs and restart

Generated data stays under
`workspace/output/genome_prediction/annevo/<species>/`:

- the untouched native `<species>.gff3`, native execution receipt and prediction
  provenance;
- `formatted/<content-key>/species_cds`, `species_gff`, `species_genome`, plus
  standard formatter summaries and source provenance;
- a separate formatter provenance record.
- `candidate_inputs.tsv` and its provenance for optional
  [rescue/refinement additions and corrections](annevo-model-improvement.md).

The standard input formatter extracts CDS from the native GFF and genome,
retaining its existing gene grouping, phase, identifier and longest-CDS rules.
The stage passes those generated files to the existing CDS resolution and enabled
annotation/QC steps. It does not copy predictions into curated `workspace/input`
files or declare BUSCO complete from prediction exit alone.

A verified unchanged prediction is reused, including after interruption between
prediction and formatting. The standard `artifact_stale_policy=stop` rejects a
changed genome, runtime/model receipt or native parameters. An explicitly
selected `rebuild` regenerates the affected stage; previous complete formatted
content generations remain available. Failed native prediction preserves its
previous complete GFF and native receipt. Cache/model availability alone does
not satisfy prediction or downstream completion checks.

ANNEVO's component uses a non-commercial license; retain its license and notices
with the cached source/model. See [runtime bundles](annevo-runtime-bundles.md)
for source, model, dependency and distribution provenance. GG's base environment
and its own license remain independently recorded.
