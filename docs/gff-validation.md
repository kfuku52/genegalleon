# GFF coordinates and coding structure

`gff2genestat.py` selects a longest CDS model before calculating structure.
For GFF3 CDS features attached directly to a gene, distinct CDS IDs are distinct
models; repeated rows with the same CDS ID are parts of one model. Models under
an explicit transcript remain grouped by transcript. Equal-length alternatives
retain the existing deterministic selection and warning.
When a root feature's structural `ID` resolves to an input gene, it takes
precedence over display aliases such as `Name`. Aliases are consulted only when
that structural ID cannot be resolved, preventing obsolete gene names from
merging distinct annotated loci.

`exception=ribosomal slippage` is accepted only when all parts declare it and
share one explicit CDS ID and Parent. Overlapping parts remain intact.
`feature_size` is the sum of distinct annotated block lengths, including reused
bases, not their genomic union or a claim of translation correctness. The
`start`/`end` fields are the genomic envelope. `splice_mode=ribosomal-slippage`
distinguishes these records; `num_intron` and `cds_first_phase` are unavailable.
Small gaps between translation parts are not inferred to be spliceosomal introns.
Tree plots omit their conventional exon/intron structure and intron-site
correspondence; scaffold taxonomy can still use the validated chromosome.
Free-text notes alone do not authorize overlapping coordinates.

Source-audited overlapping CDS parts may carry the local
`gg_source_overlap=confirmed` marker. This preserves the source coordinates
and repeated bases without misclassifying the record as ribosomal slippage or
pseudogene. The marker is accepted only when every part of one explicit CDS
model carries it; mixed or unmarked overlaps still fail.

Explicit `pseudo=true` CDS parts similarly retain their genomic location with
`splice_mode=pseudogene` and `phase_status=pseudogene`. They must share one
CDS ID; partial/mixed pseudogene declarations fail. Their intron count and global
reading frame are unavailable, and tree structure/site panels omit them. This
does not reclassify them as functional coding genes or remove their existing
identifiers from scaffold composition.

`--phase-policy strict` (default) rejects conflicting phases. With
`--phase-policy report`, valid coordinates are retained, `phase_status` is
`conflicting`, and `cds_first_phase` is missing. No phase is repaired or inferred
from a length match. Intron-site correspondence excludes these records. Other
values of `phase_status` are `consistent`, `missing`, `ribosomal-slippage`, `pseudogene`, `source-overlap`, and
`not_evaluated` (non-CDS features). Invalid coordinates, mixed coordinate systems,
and invalid phase values still fail. `consistent` describes internal phase
consistency, not verification against a genome or protein sequence.
The output also reports `cds_partial` as `none`, `5prime`, `3prime`,
`5prime+3prime`, or `unknown`. With `--validate-cds-length`, fuzzy termini are
checked against the corresponding phase and the 0--2 bases needed to complete
a terminal codon; non-partial records retain exact length validation.

Genome annotation uses report mode so scaffold evidence does not depend on a
usable global reading frame. It also uses `--require-matches`, which rejects a
zero-match result before publication. Partial coverage is not zero support;
compare mapped IDs with the input CDS set when reviewing evidence. Input
generation retains strict validation. `--validate-cds-length` remains an
independent check and does not imply nucleotide sequence identity.

Preserve existing FASTA identifiers when repairing a project with downstream
results. Repair an exact, audited GFF identifier correspondence rather than
enabling heuristic suffix matching. Retain original inputs and a change log.

Cross-contig CDS parts with one explicit transcript Parent and complete unique
`number=1..N` order are represented as `splice_mode=ordered-fragments`. Each
contig's coordinates must still agree with its strand and the declared order.
The source blocks and junction positions are retained, without assuming
trans-splicing, a common scaffold, or a genomic intron count. Conventional
structure and intron-site plots omit these models. Missing/contradictory part
order and conflicting phases still follow the existing validation policy.

## Automatic CDS resolution

Genome annotation now validates CDS candidates before sequence analyses. It
retains the original input FASTA/GFF and writes the selected CDS, coordinate
traits and per-gene decisions under `workspace/output/species_cds_resolved`.
Genome/gene evolution use a content-bound view of these resolved inputs when
available. Run genome annotation before those analyses to resolve new inputs;
a stale resolution is an error and must be regenerated, never silently reused.

A supplied CDS is preferred if its translation passes. When a reference genome
is available, exact GFF CDS blocks provide a second candidate. Explicit fragment
order and strand are retained; overlapping UTRs never subtract coding bases.
GFF phase constrains the reference candidate. Otherwise CDSKit selects padding,
and GeneGalleon checks the resulting reading frame independently. Padding is
recorded, never used to mask internal stops. A shifted frame with multiple
stop-free alternatives is unresolved. Internal ambiguous bases and
context-dependent stop codons are also unresolved; terminal partial codons may
retain padding. The configured genetic code is used throughout.

If both candidates pass, biological length before padding breaks the tie only
for an exact in-frame extension; other disagreements retain the supplied CDS and are reported.
If neither passes, that gene is excluded, while other genes continue.
This is an analysis admission policy, not a claim that an excluded gene is a
pseudogene. The report retains source hashes, phase, padding, internal-stop
counts, selection reason and coordinate usability. A CDS without exact
coordinate correspondence remains usable for sequence analyses, but its
intron/exon structure is unavailable. Conflicting UTRs alone disable UTR
information, not an independently validated CDS.

The standalone `cds_resolution.py --cds FASTA --gff GFF --genome FASTA
--output-dir DIR --genetic-code 1` command uses the same policy. GFF and genome
arguments are optional; they are required for coordinate-derived rescue.
Existing artifact-provenance checks remain active: adopting a new resolution
may require the explicit workflow `artifact_stale_policy=rebuild` setting.
