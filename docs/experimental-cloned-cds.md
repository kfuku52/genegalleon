# Prefer an experimental CDS with a reviewed genomic projection

`gg_input_generation` supports `cloned_cds_mode=experimental` alongside the
unchanged default `genomic` and transcriptome-only `sequence` modes. It adopts
the supplied CDS verbatim, replaces its reviewed source gene models and records
the differences from the unchanged reference genome. This is an explicit
experimental transcript projection, not a genome-identical annotation.

Use the same `cloned_cds_fasta`, SHA256 and species-specific ID selector as
[ordinary clone replacement](input-conventions.md#cloned-cds-model-replacement).
Each selected CDS must still contain an intact, unambiguous coding frame and
terminal stop. A biologically partial sequence can satisfy that syntactic
check; declare its known missing end explicitly.

## Reviewed targets

The manifest column `cloned_cds_projection_targets` contains a JSON object
with exactly one entry per selected clone. Obtain candidate coordinates from
the native genomic mode's `audit.json:candidates` and inspect the source GFF,
or independently verify the coordinates. This mode does not choose a relaxed
best hit automatically. Unknown or competing loci still need coordinate review.

```json
{
  "clone_A": {
    "seqid": "chr1",
    "strand": "+",
    "cds": [[100, 400, 0], [484, 910, 0]],
    "source_gene_ids": ["N_fragment", "C_fragment"],
    "partial_5prime": false,
    "partial_3prime": false
  }
}
```

Coordinates are **zero-based, half-open**, phases use GFF3 skip counts, and
blocks are checked in transcript order. All named source genes must overlap
the projected coding path on that contig and strand. An omitted overlapping
coding gene rejects publication. Coding sequence outside the reviewed path
requires an explicit partial-end declaration, so a shorter experimental
record cannot silently delete an annotated N terminus.

An unannotated locus uses `"source_gene_ids": [], "new_locus": true`. It must
have no overlapping source coding model. A new model is added while other
source genes remain. The declaration is a reviewed projection, not evidence
that the reference locus is intact or expressed.

## One model per gene

`cloned_cds_priority_ids` is a JSON array in descending preference order, for
example `["clone_A", "reference_accession"]`. Clones sharing a source gene
or overlapping coding coordinates form one locus group. Different sequences
in a group require an explicit priority; otherwise publication stops. Exactly
one model is emitted, with all other supplied clone IDs retained as aliases
and `same_gene_lower_priority` decisions. Identical sequences can be deduplicated
without a preference. The existing longer-only policy is evaluated first and
still retains an equal or shorter guarded record.

Do not supply genomic-search batch/unresolved options or sequence-only BLAST
thresholds in this mode. Targets, priority, clone hashes and the projection
implementation are frozen in the normal task plan. Changes need a new prepare
plan/workspace. Standard single/array formatting, receipts, QC and resume apply.

## Outputs and verification

GFF genes, transcripts and projected CDS carry
`gg_sequence_origin=experimental_clone` and `gg_projection_*` attributes. These
bind the experimental CDS, reconstructed reference CDS, reference coordinate
blocks, known partial ends and an exact nucleotide difference CIGAR. `M` means
an exact nucleotide match, `X` a substitution, `I` experimental-only bases and
`D` reference-only bases, along the spliced reference in transcript direction.
The CDS intervals and phases remain **reference projections**; reconstructing
them alone need not recover the experimental sequence.

The native audit records substitutions/indels with both sequence offsets,
old-to-new gene identities, priority decisions and reference/experimental
lengths. The bundle retains every supplied clone, removed source CDS and GFF,
and the projected reference CDS. FASTA archive roles with no applicable records
are omitted from `outputs`; `archive_record_counts` records their zero count.
No dummy sequence or weaker completion check is used. Those archived full source models remain
available for N-terminal/localization comparisons; a partial experimental CDS
must not be reported as a verified full-length localization input.

Formatting and exported CDS/GFF QC independently reconstruct the actual
reference intervals and validate both sequence hashes, coordinate hashes and
every CIGAR operation. Native summary QC also checks the expected projection
count against the curation audit. Missing, changed or partially removed proofs
cannot turn a difference into an ordinary exact match.

The genomic isoform catalog preserves experimental projection provenance and
the supplied CDS separately. A reference reconstruction with a sequence
disagreement remains ineligible as an ordinary genomic anchor; experimental
protein validity never repairs a disrupted reference genome implicitly.

Original inputs and the reference genome are preserved. All source ownership,
shared-feature, identifier-collision and atomic-publication guards remain.
