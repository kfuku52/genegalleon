# OrionGeno prediction and representative UTR completion

Genome annotation defaults to `genome_annotation_backend="oriongeno"`.
Set it to `annevo` for the existing ANNEVO route, or `none` for CDS-based
annotation, including annotation of an existing verified representative bundle.
Both predictors read `workspace/input/species_genome` and retain curated inputs.
Each array task handles one sorted genome; provision the array for all species.

## Native prediction

Run the genome annotation entrypoint in an explicitly allocated CUDA GPU job.
The existing SIF `--nv` / Docker GPU adapter forwards that allocation; the
predictor never requests every GPU implicitly. The first environment preparation
imports native CUDA/Triton extensions and also needs this GPU allocation and
driver mount. OrionGeno bootstraps its own
Linux x86-64 environment from the moving source branch declared in
`container/source_branches.env`, using the upstream environment and native-wheel
requirements. Source and model licenses remain with those assets. The base GG
Python environment is unchanged. See the [upstream installation guidance](https://github.com/BGIResearch/OrionGeno)
for its GPU/driver requirements and non-commercial license.

| Entrypoint variable | Default | Meaning |
| --- | --- | --- |
| `genome_annotation_backend` | `oriongeno` | `none`, `annevo`, or `oriongeno`. |
| `oriongeno_model_profile` | `plants` | Explicit checkpoint: `plants`, `mammals`, `birds`, `fish`, `other_vertebrates`, `arthropods`, `other_invertebrates`, or `fungi`. |
| `oriongeno_source_ref` | declared branch | Moving branch read from `container/source_branches.env` on the host and forwarded into the container. |
| `oriongeno_batch_size` | `1` | Native prediction batch size. |
| `oriongeno_no_download` | `0` | `1` requires verified local runtime/model assets. |
| `oriongeno_python` | blank | Optional external native Python installation. |
| `oriongeno_checkpoint` | blank | Optional local checkpoint directory; otherwise fetch only the selected checkpoint. |
| `oriongeno_source_commit` | blank | Source provenance required for an external installation. This is a one-off override, never a repository source default. |

Scoped overrides use `GG_GENOME_ANNOTATION_`, for example
`GG_GENOME_ANNOTATION_GENOME_ANNOTATION_BACKEND=annevo`. The explicit default
profile is for plants; choose the correct profile for other lineages. Different
lineages can be submitted as separate tasks with their own profile overrides.
The integrated native/candidate path requires nuclear genetic code 1.

Source, checkpoint and installed-runtime identities are frozen on first use.
Unchanged reruns verify and reuse them offline; use a new cache/receipt location
to adopt another source or checkpoint snapshot. Changed or incomplete artifacts
stop the stage under the usual stale-artifact policy. Failed native execution
retains logs and the previous complete GFF/receipt pair.
External Python paths retain their virtual-environment prefix. Preparation and
prediction publication use the shared-filesystem namespace lock; all concurrent
jobs sharing those paths must use the same lock protocol.

For separate species workspaces, set `GG_COMMON_DOWNLOADS_ROOT` to one existing,
owned physical absolute directory and bind that directory at the same path in
the container. OrionGeno environments/checkpoints and downloaded BUSCO lineages
then share the existing asset locks and provenance checks. The default remains
each workspace's `downloads` directory. Keep one source/model snapshot per cache;
use a new cache for an intentional upgrade.

## Separate acquisition, GPU and CPU allocations

`workflow/support/genome_annotation_cohort.py` runs one frozen member with
`--manifest COHORT.json --manifest-sha256 SHA --task N --phase acquire|predict|busco`.
It runs on the host and invokes the installed native entrypoints. It does not
submit jobs. A scheduler can connect arrays with Slurm `aftercorr` so one slow
download does not delay other species, and run BUSCO without reserving a GPU.

The version-1 JSON contains `project_root`, `runtime_root`, `image_path`,
`work_root`, `downloads_root`, `busco_lineage` (`embryophyta_odb12`) and `members`.
All roots are existing physical absolute paths under the project, which is bound
into the container at the same path. Each member has a consecutive 1-based `task`,
unique `species`, `filename` (`SPECIES.genome.fa` or `.fa.gz`) and `source`:
either `local_path` plus `sha256`, or a public `url` with optional `md5` and
`archive_member`. MD5 applies to a directly downloaded file, not an archive member.

Acquisition records the actual genome hash and hard-links its input. Prediction
uses a separate species workspace and the shared native runtime/model cache.
The CPU phase audits current prediction/formatting provenance before linking the
formatted triplet into its own BUSCO workspace. A missing, changed or incorrectly
labelled BUSCO artifact fails the phase even if its entrypoint returned zero.
Original annotations remain separate. Resource limits, duplicate exclusion and
accepted-job verification belong to the scheduler submitting these phases.

Outputs are below `output/genome_prediction/oriongeno/SPECIES/`: native GFF3,
execution receipt, prediction/formatter provenance, formatted GFF/CDS/genome,
and `candidate_inputs.tsv`. These outputs do not establish annotation accuracy
or BUSCO completeness by themselves.
Native GFFs can encode UTRs through exon sequence outside CDSs without separate
UTR feature rows. Representative completion derives those intervals and emits
explicit UTR rows in its final GFF views.

## Connect the whole cohort to representative selection

Input generation defaults `gene_model_oriongeno_inputs` to `auto`: completed
`output/genome_prediction/oriongeno/SPECIES/candidate_inputs.tsv` tables for the
selected cohort are discovered automatically when rescue/refinement plans are
prepared. A plan freezes the tables available then; newly generated predictions
require a new plan/output directory. If none exist, existing-annotation UTRs are
still processed and the audit records missing native predictions. Set the option
to blank to disable native import, or to an explicit table/prediction directory.

For an explicit portable cohort table, run inside the GG container:

```bash
python workflow/support/oriongeno_runtime.py collect \
  --prediction-root workspace/output/genome_prediction/oriongeno \
  --output workspace/output/oriongeno_candidates.tsv
```

Supply that table to input generation:

```bash
GG_INPUT_GENE_MODEL_ORIONGENO_INPUTS=/absolute/path/oriongeno_candidates.tsv \
GG_INPUT_RUN_GENE_MODEL_REFINEMENT=1 \
  bash workflow/gg_input_generation_entrypoint.sh
```

The independent rescue/refinement planners accept `--oriongeno-inputs` too.
Refinement inherits a verified OrionGeno table from a rescue publication.
Each row binds `species`, native `gff`, and its native execution `receipt`.
The target assembly, GFF and receipt are checked at plan/execution/publication
boundaries. Old `gene_model_annevo_inputs` / `--annevo-inputs` remain supported.
Native coding candidates pass the same independent donor, genomic, ownership
and representative-adoption gates as existing candidates. A predictor cannot
provide its own independent structural confirmation.

## Longest projected UTR per side

After the final representative CDS is selected, every species and every
exported representative passes UTR completion. UTR length does not affect CDS
selection. Candidate transcript ends come from existing annotation and matched
OrionGeno models, including models with different CDS starts/stops/internal
exons. Native models first match a unique selected gene with an identical CDS
exon path on the same contig and strand, then use unique same-strand coding
overlap when no exact path exists. Merged or copy-ambiguous matches stay unassigned.

For each donor transcript and each side, keep its outer transcript end fixed
and retain its exonic sequence outside the selected CDS boundary. This can
extend a UTR into the donor's former CDS, or consume part/all of its former UTR.
Introns never count toward UTR length. Existing annotations without explicit
UTR rows can supply UTRs through their own exon/CDS difference. A donor without
a UTR on a side does not manufacture one there.

Compare projected spliced lengths separately for 5-prime and 3-prime UTRs.
Choose the longest candidate per side; equal lengths prefer existing annotation,
then the selected model, then OrionGeno, with stable transcript-ID tie breaking.
Different sides may use different donors. Each side retains one complete donor
exon path rather than a union of competing UTRs. CDS coordinates, phase, sequence,
protein, gene IDs and representative choice remain unchanged.
When no UTR is available and the exon path is unchanged, completion preserves
the GFF bytes while still exporting the transcript sequence and audit.

When a boundary falls in a donor intron, retain the projected UTR and mark the
new terminal connection, its genomic motif and whether it is canonical in the
audit. Motifs do not provide RNA support or add a UTR-loss criterion. A CDS
reaching the donor's outer transcript end removes the UTR on that side. Missing
source UTRs, unmapped predictions and unsupported source models are also audited.
An overlap with another coding locus is not used to silently trim a UTR.

The resulting model is a derived prediction and may combine different source
transcripts. Longer UTRs are a coverage policy, not a claim of greater accuracy,
experimental support or confirmed complete transcript ends.

The representative bundle adds:

- updated selected/full/analysis GFF views with consistent exon and UTR parents;
- `species_transcript/SPECIES.fa`, spliced sequences keyed by selected transcript
  IDs, separate from unchanged CDS/protein FASTAs;
- `utr_audit/SPECIES.json`, per-side candidates, source identities, adjusted
  blocks/lengths, boundary shifts and newly inferred junctions;
- `utr_summary.tsv`, cohort denominators and 5-prime/3-prime/both-side coverage;
- transcript paths/hashes in `effective/inputs.tsv` and all new files in the
  publication receipt. Exact original annotation remains in `source_annotation`.

FASTA headers identify the projection policy; presence of a transcript sequence
does not certify full-length UTRs. Missing/unexportable models stay represented
in the audit and coverage denominator. With no OrionGeno table, completion still
compares existing source isoforms and records unavailable native matches.
