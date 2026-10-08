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

Source-audited CDS models with `splice_mode=source-overlap` retain their
coordinates for sequence and scaffold analyses but are omitted from the
ordinary exon/intron drawing, because repeated genomic bases do not define a
single linear genomic geometry. They are reported as excluded diagnostics
when intron correspondence is requested.

## Gene cluster membership

The membership column has a graphical legend: colored circles and connecting
lines mark clusters containing at least two tips from this gene family; pale
gray circles mark singletons. Clusters are formed separately within each species
and chromosome/scaffold, splitting when the gap from the rightmost covered
coordinate to the next gene exceeds the displayed maximum distance in bp.
Overlapping and nested intervals remain together; a chain can span more than
the maximum gap. Different clusters within a species use
different shades of its tip-label color, within bounded lightness ranges so
large cluster counts do not collapse to repeated white/black symbols. Very pale
tip colors start from a darker shade to distinguish clusters from singletons.
A thin pale gray background line joins
the species' tips across cluster boundaries. Missing taxon, coordinates, or
scaffold assignments remain blank. These symbols describe family-gene proximity,
not HGT direction, conserved gene order, or host-scaffold background support.

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

## Explicit trans-splicing

CDS rows that all declare `exception=trans-splicing` may use `part=1`, `part=2`,
etc., or the `part=X/Y` form to specify their complete transcript order.
[NCBI documents this ordering attribute](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/file-formats/annotation-files/about-ncbi-gff3/#unofficial-attributes).
Missing, conflicting, or incomplete part orders are rejected. GeneGalleon does
not infer trans-splicing from inconsistent coordinates alone.

The summary preserves numeric `feature_blocks` in that order and adds
`splice_mode`, `feature_block_sequences` (semicolon-separated, URL-escaped
sequence IDs), `feature_block_strands`, and `transcript_junction_positions`
(cumulative CDS lengths before each join). These fields also accompany ordinary
CDS rows. Genome-to-CDS extraction uses each block's own strand before joining.
A mixed-strand transcript has scalar `strand=?`; a multi-contig transcript has
no single chromosome/start/end envelope. Per-block coordinates remain complete.

Trans-spliced CDS remain in sequence analyses and structure plots. Their joins
are drawn as dashed vertical boundaries on cumulative CDS coordinates, without
an invented genomic gap. Because the transcript-level exception alone does not
classify each join as a cis intron, `num_intron` is missing and
`intron_positions` is empty for these transcripts. Cis-intron correspondence
reports `trans_splicing` explicitly rather than treating them as intron-free.
CDS length checks still apply. Trans-spliced transcripts with explicit UTRs are
rejected until their UTR order can also be represented unambiguously.

## Shared Pfam filter for focused HGT

The focused HGT PDF's second page compares each focal and displayed neighbor
gene's saved MMseqs2 LCA classification/host match with its separate Swiss-Prot
best-hit product and taxonomy. Missing and unresolved classifications remain
explicit. Overlapping and intron-hosting genes are included in the bounded
neighborhood. See [context inputs and provenance](host-scaffold-taxonomy.md).

Focused category-1 HGT tables and plots now require at least one donor/recipient
gene pair sharing an exact query Pfam accession, at least 50% shared-domain
query coverage on each protein, and a shorter/longer protein length ratio of
at least 0.5 by default. Both genes must be retained
descendants linked to the same modeled event and individually satisfy the
existing class-level scaffold-background thresholds. One qualifying pair
retains the event; other eligible descendants remain in its context tables.
Events, genes and independent orthogroups are counted separately.

`hgt_summary_focus_min_shared_pfam_coverage=0.5` sets this inclusive fraction;
`0` restores any shared query Pfam. For each exact pair, take the union of saved
query amino-acid intervals for its shared accessions and divide by each query's
own protein length. Overlaps count once. Coverage is not pairwise alignment
coverage; separate pairs cannot supply separate sides of a passing decision.
Event/pair audits retain both lengths, covered amino acids, coverage and a
traceable best pair. Review flags identify repeat/generic-binding-only matches,
differing domain sets and proteins shorter than 100 aa. These flags are not
exclusions or proof of a partial gene model. The enabled protein-length-ratio
criterion independently excludes pairs below its threshold. Domain-set
differences alone do not establish incompatible architecture.
For a passing event, the recorded best pair is selected from its passing pairs;
a bilateral no-hit exception retains both gene IDs with unmeasured coverage.

The filter reads saved `rpsblast/<OG>_rpsblast.tsv` (also the existing
`<OG>.rpsblast.tsv` convention), including ZIP-backed families. It does not
use best-hit `pfam_ids`, borrow neighbor domains, infer domain architecture
equivalence, or run additional sequence searches. Exact Pfam sharing is an
additional candidate criterion, not proof of HGT or complete sequence quality.

In `gg_gene_summary`, `hgt_summary_focus_require_shared_pfam=1` enables the
Pfam requirement by default. Independently,
`hgt_summary_focus_require_length_ratio=1` requires positive measured lengths
and `hgt_summary_focus_min_length_ratio=0.5` sets the inclusive shorter/longer
ratio. Both enabled rules must pass for the same exact pair. Query-protein
lengths are measured in amino acids; family statistics, best-hit lengths and
alignment lengths are not substitutes. Missing lengths do not pass an enabled
length rule. Disable only the length rule to reproduce the preceding Pfam-only
cohort, or disable both requirements for the earlier scaffold-only cohort.
`hgt_summary_focus_allow_both_no_pfam=0` excludes pairs without a shared hit.
Setting this option to `1` additionally allows pairs where **both genes have
explicit searched-no-hit rows**. A no-hit on only one side, disjoint detected
domains, and missing search/query records still fail. These correspond to
`focus_hgt_traits.py --require_shared_pfam 1 --allow_both_no_pfam 0`.
The bilateral no-hit opt-in is an explicit exception to domain coverage;
coverage stays unmeasured, and the enabled length rule still applies. Saved
hits receive no additional E-value cutoff.

The exact-pair step evaluates the shared input cohort once; the combined taxonomy/trait step selects each
trait's category-1 recipients from that same passing set. Existing Pfam-named audit
files retain their paths when only the length rule is enabled. The bundle root
contains `pfam_events.tsv` and `pfam_event_audit.tsv`, `pfam_pair_audit.tsv` and
`pfam_gene_audit.tsv`, recording all input events before trait selection, exact
gene pairs, rejection reasons, detected accessions and source hashes.
Per-trait audit paths are retained as trait-specific views of those shared decisions.
Selected event tables add `pfam_*` decision/count columns; per-tip tables, native tree
PDFs and the three aggregate figures use the same filtered event set. This
filter also runs when plotting is disabled. Missing evidence remains unknown
and does not qualify through the bilateral no-hit option.

Python callers may supply booleans or explicit `0`/`1` (`false`/`true`) flags;
ambiguous flag values fail before output generation. Direct native-tree and
summary exporters also validate event IDs, exact event-gene identities and
retained lineage eligibility. Summary product labels use only selected events,
and an explicit native no-hit record cannot be overridden by a stale link hit.
