# CDS and intron structure in tree plots

The default gene-family tree plot includes a structure column between the intron
count and protein domains when GFF-derived block coordinates are available.
Boxes show CDS in the tip-label color and explicitly annotated UTRs in a
paler tint of that color. UTR boxes are also thinner. Lines show intervening
introns. Unannotated UTRs are not inferred from transcript spans or other isoforms. Each row runs from the 5′ to 3′ coding end, including minus-strand
genes. Boxes, lines, and abbreviated labels beside pairwise heatmaps use the
same tip-label colors as the main labels.

The GFF summary now retains feature_blocks (semicolon-separated, 1-based
inclusive start-end pairs in transcript order), feature_type, gff_transcript_id,
utr_blocks, and cds_first_phase. UTRs must reference the selected transcript
exactly; CDS/UTR overlap is rejected. Fully unknown CDS phase stays missing; consistent later phases can establish a missing initial phase.
Coordinates come from the same selected transcript used for intron counts:
the existing transcript-ID matching and longest-CDS selection policy. Isoforms are not merged. Existing GFF summaries
need regeneration before the structure column can appear.

Use --panelN=gene_structure,compressed,23 with stat_branch2tree_plot.r
to select the default 23 mm column. CDS and UTR widths retain their respective sequence lengths;
each intron of length L bp occupies 100 * ln(1 + L/100) display units.
The same transformation and horizontal scale apply to every row. The caption
identifies the compression; mixed display units are not presented as genomic bp.
Use gene_structure,linear,23 for a common uncompressed genomic scale instead; the axis uses kb for spans of at least 1,000 bp, otherwise bp.
Width is physical and does not depend on the figure width. Intron counts do
not expand the requested width; positional numbers are retained only in the
underlying correspondence data.

GFF-missing tips remain blank; an observed intron-free CDS appears as one box.
Protein-domain intron marks are hidden when the structure column is present,
regardless of panel ordering. Without a structure column the existing marks
remain available. The separate intron-count column is retained.

## Intron correspondence within the structure column

The default CDS-mode plot annotates introns directly inside the exon/intron
structure column. No separate correspondence matrix is added. Use
`--panelN=gene_structure,compressed,23,untrimmed_cds_alignment.fa` to request
this explicitly; gzip FASTA is supported. Omit the alignment argument to draw
structure alone. Linear mode supports the same fourth argument.

The data-panel width defaults to 23 mm (about one third of the earlier AHA
figure's 69 mm column). Intron numbers and their leader lines are not drawn,
and the width no longer grows with the number of introns. CDS/UTR boxes retain
their tip-label colors and stay centered on each row. Positional site IDs
(I001, I002, etc.) remain in the underlying correspondence data.
Background bands anchor at actual intron positions:

- Darker bands (opacity 0.16) connect exact alignment-position and known-phase
  matches.
- Lighter bands (opacity 0.05) indicate position-only candidates within **3 aligned
  nucleotides** (one codon), including unknown or different phase.

For each intron, the next row in tip order containing a nearby candidate is
considered. Only reciprocal unique nearest matches are connected; ties are not
resolved arbitrarily and an ambiguous row is not bypassed. An unknown-phase
occurrence interrupts a dark band, replacing adjoining links with light bands.
A gap-spanning or ambiguous-base boundary has no resolved position and is not
connected. Band intensity is a display category, not a probability. Near-position
links do not merge correspondence numbers or establish common ancestry. A band
crossing another row does not assign an intron to that intervening gene; only
its endpoints are matched observations.

An observed CDS intron offset is mapped to its two flanking nucleotide positions
in the **untrimmed** CDS alignment. Both must be unambiguous A/C/G/T bases and
adjacent alignment columns for a resolved position. Coding phase is
`(CDS offset - first CDS phase) mod 3`. Different known phases are separate
groups even at the same alignment position. Unknown phase is treated as a
position-only candidate, never as a resolved phase match.

A gap-spanning intron is not snapped to a nearby site. Missing GFF/alignment
data are excluded from correspondence assignment. A supplied alignment with a CDS-length mismatch stops structure plotting.
The observed intron-count column continues to count CDS introns, as before.

These are positional correspondence candidates, not a reconstruction of intron
ancestry or proof of homology. The approach follows the position/phase comparison
principle described by [GenePainter](https://genepainter.motorprotein.de/help),
with explicit uncertainty instead of merging across alignment gaps.

## Validation details

Only positive genomic gaps between CDS blocks count as CDS introns; adjacent
CDS fragments remain separate annotated boxes but introduce no intron or site ID.
Every known CDS block phase is checked against the cumulative CDS length in
transcript order. Conflicting phases are rejected. A missing initial phase can
be inferred from later known phases only when they imply one consistent frame;
fully unknown phase remains unknown. Empty intron families produce empty site
and connection tables without failing.

The workflow validates selected GFF CDS lengths against the nucleotide input
(`gff2genestat.py --validate-cds-length`); protein-mode inputs are not subject to
this nucleotide length check. A structure plot supplied with an alignment also
stops on a CDS-length mismatch, so an excluded match cannot silently leave an
incorrect count or structure in a PDF. Missing GFF remains missing rather than
zero. Length equality alone cannot establish sequence identity.

For the legacy CoGe nested-mRNA export, `repair_coge_transcripts.py` provides an
explicit, sequence-verified input repair. It requires a single linear chain,
shared CoGe feature ID, nonoverlapping CDS/exon pairs, matching total CDS length,
and agreement of every non-N CDS base with reconstructed genomic sequence. It
writes a new GFF and a JSON audit, preserving source inputs. It does not merge
ordinary alternative transcripts or infer missing phase.
