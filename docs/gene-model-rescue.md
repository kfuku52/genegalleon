# Synteny-guided gene-model rescue

`workflow/support/rescue_gene_models.py` is an independent, restartable CLI.
Input generation optionally runs it after formatting, initial BUSCO and species
taxonomy, before augmented CDS/GFF inputs are used in OrthoFinder. Enable it with
`GG_INPUT_RUN_GENE_MODEL_RESCUE=1`; it is off by default. Evidence-based repairs
of existing CDS models already run during ordinary input formatting, before
longest-isoform selection, independently of this flag. That step exports matching
CDS/GFF corrections and an adjacent `*.cds-normalisation.json` audit while
preserving acquired sources; see [input conventions](input-conventions.md).
The rescue flag controls synteny searches for additional models.

```mermaid
flowchart LR
  A[Formatted CDS, GFF, genomes] --> B[Initial BUSCO and initial tree]
  B --> C[Freeze five common references and three nearest donors]
  C --> J[Audit and admit existing protein anchors]
  J --> D[Deduplicated pair comparisons and unquota self synteny]
  D --> E[Two-anchor candidate intervals]
  E --> F[miniprot and optional GeMoMa]
  F --> G[ORF and conflict validation]
  G --> H[Augmented CDS and GFF]
  H --> I[BUSCO for changed species and downstream orthogroups]
```

## Inputs and reference selection

Provide exactly one formatted CDS FASTA, GFF3 and genome FASTA per species.
Species prefixes and IDs follow [input conventions](input-conventions.md).
CDS-to-GFF mapping is complete and unambiguous; representative isoforms are
selected with the existing synteny mapper. Existing models undergo automatic
anchor admission before comparisons. Original files are never rewritten.
Genome contig names must be unique, and every original coding-feature span and
selected anchor must lie within a matching genome contig, including species with
no candidates. FASTA indexing warnings fail the attempt and remain in its logs.
These checks detect coordinate incompatibility; they do not establish exact
sequence agreement between every original CDS and its genomic exons.

### Existing-model anchor admission

Ordinary stop-free CDS use their original translation under the species' genetic
code. For an invalid translation, two narrowly supported normalisations are
available as a fallback for legacy or externally formatted inputs, applied only
to the temporary anchor protein:

* UTR contamination requires agreement with the same model's annotated spliced
  exons, and a stop-free, in-frame genomic CDS contained within that transcript
  and the supplied sequence. At most two terminal formatter-added Ns and the
  legacy conversion of ambiguity codes to N are allowed in identity comparisons.
  If the supplied sequence matches an annotated CDS, that model takes precedence;
  a shorter ORF in another isoform cannot reclassify its internal stop as UTR.
* Partial CDS require agreement with the same model's genomic CDS and a coherent
  chain of annotated phases. Standard GFF3 skip-count phases and complementary
  frame fields are distinguished from their exon-length recurrence. Ambiguous
  single-block models require unanimous informative evidence elsewhere in the
  source GFF. Missing, inconsistent or mixed evidence cannot establish a frame.
  Only complete codons after the annotated offset are translated; omitted
  terminal bases are recorded. Newly rescued models never establish a source
  file's phase convention.

Translation exceptions, annotated pseudogenes, conflicting reconstructions and
unexplained internal stops are withheld from both anchors and donor queries.
Stops are neither deleted nor replaced with X, and a stop-free alternative frame
alone is insufficient evidence. These records remain in exported CDS/GFF and in
the original-feature overlap checks. Withholding a model is not a gene-loss call.
If no usable anchors remain, preparation fails with its admission audit intact.
New predictions still require the complete ORF and other checks below.

Admission also recognises GeneGalleon's formatted NCBI `GeneID123` identifiers
against `GeneID:123` in `Dbxref`, and GWH gene `Accession` identifiers. These
exact aliases select the public annotation mapper's feature/attribute pair;
full, unambiguous CDS-to-GFF mapping remains required. The general pairwise
synteny command retains its existing strict translation and mapping behaviour.

Admission JSON records the genomic transcript, phase decision, sequence hashes
and reasons for every normalised or excluded representative; a TSV covers all
representatives. The helper's identity and source genome hashes are frozen in
the plan, and admission files participate in prepared/finalized receipts.

All first-pass BUSCO short summaries must contain `C:`, dataset, BUSCO version,
mode and marker count (`n:`). Comparisons require the same dataset/version/mode,
marker count and dataset creation date when present. Matching lineage names
alone cannot distinguish database revisions. Complete BUSCO is
S+D, so WGD duplication does not lower the score. NWKIT `sample --method max-pd`
selects five common references among species with completeness at least 90%.
BUSCO rank breaks equal phylogenetic-diversity gains; it does not define a
weighted quality-diversity objective. Fewer than five eligible species is an
error; change the documented threshold explicitly if appropriate.

Each target also receives its three nearest species by tree distance, with
BUSCO completeness and species name as tie-breakers. Remove self and deduplicate
shared donors and unordered pairs. Every species receives an additional raw,
unquota self comparison. At 500 species this gives at most 4,000 pair assignments
before deduplication and 500 self jobs, compared with 124,750 all-pairs jobs.

Use an initial external tree with `GG_INPUT_GENE_MODEL_RESCUE_TREE`, or the
generated `output/species_taxonomy/taxonomy_tree.nwk` with the default `auto`.
The CLI requires an explicit `--tree`; it does not discover downstream inferred
trees. A tree without positive branch lengths uses unit edges and records that
choice. An NCBI taxonomy topology supplies no divergence-time information.
If any non-root edge has a positive length, every non-root edge must have an
explicit length. Partially specified lengths are rejected rather than treating
missing edges as zero distance. Root stem lengths do not affect this decision.

The plan freezes all source hashes, genetic codes, thresholds, reference lists,
tool identities (including executable and alignment/annotation source hashes)
and comparison indices. Post-rescue BUSCO never changes this
selection. Changed inputs/settings/tools require a new output directory.

## Candidate and model evidence

Pair comparisons use JCVI/DIAMOND; self comparisons use JCVI/LAST alignment and
raw `kffractbias.selfscan`, with no quota screening. Self alignment removes
same-ID hits, preserves high-identity paralogs and disables tandem collapsing.
Pair alignment also disables tandem collapsing. The selfscan dependency uses
four anchors and chaining distance 20; `--min-anchors` and `--distance` tune
pair comparisons. A successfully scanned pair without anchors is complete
empty evidence, not a failed job or a gene-loss call.
For each chromosome pair within a block, adjacent target anchors define
intervals. Multiple donor anchors at either flank are matched in donor rank
order, maximising matched paths and minimising rank distances; all equally good
matches are retained. This keeps parallel copies when a lifted block combines
WGD segments on the same chromosome. It is a local pairing rule, not an
orthology assignment or a global copy quota. Both target and donor flanking IDs
are recorded. Donor genes between the matched flanks nominate searches.
Both directions of a self block are examined. A donor gene can nominate multiple WGD intervals;
another retained homolog elsewhere does not remove a candidate.

Each distinct interval is searched independently, batching its donor proteins.
This prevents stronger paralogs in other intervals suppressing the local hit.
The default maximum interval is 200 kb and maximum intron is 20 kb; both are
configurable. A large/rearranged/poorly anchored region can remain unexamined.
No candidate or no alignment therefore does not imply gene loss.

The default engine is [miniprot](https://github.com/lh3/miniprot), including its
embedded PAF evidence. Automatic acceptance requires at least 95% donor protein
coverage and 50% identity, an intact start and terminal stop, consistent CDS
phases, canonical GT–AG/GC–AG/AT–AC splice junctions, and no assembly ambiguities,
frameshifts or internal stops. If the predictor omits the terminal stop from CDS,
the export extends only when the actual next genomic codon is a stop. Minus-strand
models are reconstructed in transcription order. Coverage is the fraction of
donor residues aligned to genomic residues, including split codons, excluding
query insertions; an alignment spanning both ends is insufficient. Both the
aligned fraction and query-span fraction are retained in model evidence. These are conservative
defaults, not calibrated species-independent sensitivity estimates.

Unresolved queries optionally receive a whole-genome search. Every engine,
including searches with padding, requires the assembled CDS to stay inside the
expected two-anchor interval. Outside hits remain unresolved evidence, including
possible relocations and other WGD copies. No disruption is repaired or masked. Models
overlapping any original gene/transcript/CDS, or conflicting with another new
model, are withheld. Such conflicting predictions do not mark a query resolved
and suppress optional refinement. Queries recovered by the whole-genome fallback
are rechecked before GeMoMa dispatch. Identical coordinates from multiple donors
are consolidated with their supporting records. Accepted IDs are stable hashes of genomic exon
coordinates; neither orthogroup IDs nor run order define them.

Optional `--gemoma-jar` / `GG_INPUT_GENE_MODEL_RESCUE_GEMOMA_JAR` runs the
[GeMoMa pipeline](https://www.jstacs.de/index.php/GeMoMa-Docs) for selected donor
transcripts and target regions, using Java and tblastn. It is an explicit local
jar dependency, not downloaded at execution. Select its compatible Java with
`--gemoma-java` / `GG_INPUT_GENE_MODEL_RESCUE_GEMOMA_JAVA`. The official 1.9 jar
uses a JavaScript engine; its absence was reproduced with Java 24, and a
separately selected Java 8 was verified. This requirement should be rechecked when the jar is
updated, rather than changing the container's default JVM. Its CDS models undergo the same
ORF/conflict checks and a donor-protein alignment for coverage/identity. Use
standard genetic code 1 for this optional refinement. Ambiguous reference
transcript mappings and references requiring complementary phase normalisation
are recorded and skipped by GeMoMa; miniprot evidence remains.

## Commands and outputs

```bash
python workflow/support/rescue_gene_models.py plan \
  --cds-dir workspace/output/input_generation/species_cds \
  --gff-dir workspace/output/input_generation/species_gff \
  --genome-dir workspace/output/input_generation/species_genome \
  --busco-dir workspace/output/input_generation/species_cds_busco_short \
  --tree initial_species_tree.nwk --output workspace/output/input_generation/gene_model_rescue

python workflow/support/rescue_gene_models.py run \
  --output workspace/output/input_generation/gene_model_rescue --cpus 4
```

`run` executes all comparison jobs, all species rescue jobs and finalization.
Alternatively use `synteny --task-index N`, `rescue --task-index N` and
`finalize` independently. Indices are one-based and frozen. `status` reports
pending jobs. `--cpus` is the total threads for one sequential CLI worker;
parallelism across jobs is supplied by the scheduler.

| Output | Meaning |
|---|---|
| `plan.json`, `references.tsv`, `selection.nwk` | Frozen sources, selection and sparse job graph |
| `prepared/SPECIES/` | Admitted proteins, BED, ID map and original-source metadata |
| `prepared/SPECIES/genes.anchor_admission.json`, `.tsv` | Existing-model decisions, evidence and counts |
| `synteny/comparison_NNNNNN/` | Raw anchors, all blocks, commands and logs |
| `rescued/SPECIES/candidates.json` | Donor gene, flanks, interval, orientation and comparison |
| `rescued/SPECIES/models.json`, `audit.tsv` | Accepted, duplicate-support and unresolved models with reasons |
| `effective/SPECIES/` | Validated per-species augmented inputs, ready for BUSCO workers |
| `augmented/species_cds`, `augmented/species_gff` | Original records plus accepted new models |
| `augmented/inputs.tsv`, `summary.json` | Explicit effective input paths and model counts |
| `augmented/anchor_admission/SPECIES.json`, `.tsv` | Admission audits after full export validation |
| `qc_report/before_after.tsv` | Initial/post-rescue BUSCO and changes, after the input-generation QC stage |

Exported FASTA keeps all original IDs, headers and sequences. GFF keeps all
original features, including any embedded FASTA. New transcript/gene IDs follow
the source mapping granularity. Finalization checks full CDS-to-GFF mapping and
frozen genome hashes. Point downstream input directories to the **augmented**
CDS/GFF directories together; using first-pass CDS silently omits rescued models.
The usual input-generation outputs themselves are not automatically installed
into `workspace/input`, and the rescue stage preserves that convention.

Input generation reruns BUSCO only for changed species, in a separate `qc/`
namespace, reusing first-pass tables for unchanged species. The independent CLI
can validate externally generated post-rescue summaries with
`qc --busco-dir DIR --output DIR`. That command checks comparability and writes
the before/after report; it does not invoke BUSCO itself.
QC report inputs and BUSCO summaries are hashed before reading and checked again
before publication. Worker completion likewise freezes model/export/QC files
before checking comparability and rechecks files and dependencies before writing
its receipt. Concurrent replacement therefore cannot publish a completed receipt
for evidence different from that checked.

Per-job locks and hashed receipts support retries. Results are published only
after successful computation; failed attempts retain their logs in `.failed`
directories and leave previous completed output intact. Genomic indexes are
scratch within one species attempt. Successful interval/model evidence persists.
A publication journal records replacement of an existing result and the expected
new job key. After SIGKILL,
the next attempt under the same job lock restores the previous directory or
retains a verified new directory with the matching key and removes its backup.
An invalid new result, including one with another job's receipt, is quarantined.
Unsafe or malformed journal paths fail without recovery writes.
This covers process interruption; power-loss durability is not established.
Comparison receipts bind prepared annotations; model receipts bind comparisons
and prepared annotations. Dependencies are checked again under the stage lock
before reuse and after computation, so replacements during a job invalidate
publication. Finalization also checks frozen inputs after export.
Finalization verifies comparison and prepared-annotation outputs as well as
model receipts. Malformed receipts and receipts belonging to other jobs are
pending work. Effective-input validation errors fail before QC dispatch; final
aggregation reuses verified worker QC even if `overwrite=1`. Array retries
check the full dependency chain, including comparisons repaired successfully
before the retry, instead of trusting a worker receipt alone.

## Array input generation

```bash
python workflow/gg_input_generation_array.py \
  --task-plan /path/workspace/output/input_generation/tmp/task_plan.json \
  --rescue --cpus 4 --memory 32G --max-running 8 --submit
```

This adds `rescue_synteny` comparison arrays, `rescue_models` species arrays and
`rescue_finalize` after the initial prepare/worker/finalize chain. Initial
finalization is waited on so reference selection and task counts exist before
submission. Use `--rescue-output` for a custom output location. All settings use
the normal `GG_INPUT_*` forwarding registry. External trees, jars and Java
installations must be reachable at their configured paths inside the runtime;
bind an external Java installation as a directory, including its libraries.
`--retry` checks receipts and refuses
submission while rescue jobs are active or scheduler status is unavailable.
Each species worker also performs its post-rescue BUSCO and publishes a receipt
only after both model export and QC succeed; shared finalization verifies/reuses
these outputs. A failed BUSCO therefore cannot make a model worker look complete.
Worker receipt file paths are relative, so the host helper can verify them
even when the container uses different mount paths. The rescue output location
has its own phase lock, including when multiple workspaces share a custom path.
The helper previews submissions without `--submit`. UGE/PBS users can submit
these same core modes manually with their usual one-based task IDs.

This stage recovers supported models and records unresolved evidence. Orthology,
copy-specific loss and ancestral-copy reconciliation belong to downstream gene
trees and synteny analyses. Self synteny alone does not date a WGD or identify a
missing ancestral copy definitively.

## Validation and scale

The runtime tests mask intact/disrupted models, exercise actual WGD pair/self
comparisons on shared and separate chromosomes, and reconstruct split codons in both strands and all three intron
phases. They also test transcript-level CDS/GFF export, corrupted dependencies,
failed publication, dependency replacement during prediction, padded interval
boundaries, conflict-aware refinement, worker QC/retry and a balanced 500-tip reference-selection
plan. The 500-tip test uses tiny synthetic sources; it does not benchmark 500
real plant genomes or establish biological sensitivity/false-positive rates.
Optional GeMoMa 1.9 was additionally exercised with Java 8 on a hidden intact
single-exon model and a minus-strand, two-exon model with a split codon.
Linux arm64 Docker has been tested; SIF execution has not.
Additional tests cover actual SIGKILL at three publication boundaries, corrupted
new outputs and unsafe recovery journals, four concurrent processes publishing
one job, compressed inputs and literal numeric/missing-like contig names, partial
tree lengths, incompatible genome coordinates, duplicate contigs and QC input
replacement. Local donor-flank matching is also checked against exhaustive
ordered assignments on small randomized cases.
Mixed genetic codes are exercised with a code-4 target and code-1 donors,
including an internal TGA that codes for tryptophan in the target, in both local
and whole-genome miniprot searches.

The interval-overlap microbenchmark below measures consolidation only, excluding
synteny and alignment. On `local/genegalleon:rescue-dev-20261003` (image
`sha256:52af7d7cd2756a53fcad19b1b5534e81ef1e6203c7d5f3c1e5a1d0b00cb4773a`),
three repetitions before/after the interval-index change gave median wall times
3.709/0.0122 s (about 300-fold for this operation). Maximum process RSS was
62.5/66.6 MiB, including imports, inputs and copies. The serialized outputs had
the identical SHA256
`eed4fa36ca712a6788d8ac70e35fca007a789c30fab817179f5459fe30c9970c`.
The indexed method trades additional interval storage for fewer comparisons;
these numbers do not predict whole-workflow speed.

```bash
GG_CONTAINER_DOCKER_IMAGE=local/genegalleon:rescue-dev-20261003 \
bash workflow/tests/run_in_runtime.sh python - <<'PY'
import copy, hashlib, json, resource, time
from workflow.support.rescue_gene_models import consolidate
existing = [{'seqid': 'chr1', 'start': i * 100, 'end': i * 100 + 60}
            for i in range(30000)]
models = [{'seqid': 'chr1', 'strand': '+',
           'cds': [[4000000 + i * 100, 4000000 + i * 100 + 60, 0]],
           'problems': [], 'query': str(i), 'evidence': {'i': i}}
          for i in range(2000)]
for _ in range(3):
    work = copy.deepcopy(models)
    start = time.perf_counter()
    output = consolidate(work, existing, 'Plant_species')
    print({'seconds': time.perf_counter() - start,
           'peak_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
           'sha256': hashlib.sha256(json.dumps(output, sort_keys=True).encode()).hexdigest()})
PY
```
