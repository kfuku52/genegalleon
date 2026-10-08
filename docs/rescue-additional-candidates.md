# Additional rescue exploration and verified prediction reuse

Two flanking synteny anchors remain the first candidate source. For the selected
nearest and common donor references, the existing pairwise protein comparisons
also nominate prepared donor genes without a full target match, with partial
matches, or with ambiguous full matches. A full match requires the configured
minimum protein identity and donor coordinate coverage. Multiple disjoint target
proteins are never combined to claim a full match. Matches within the configured
`cscore` fraction of the best full-match bit score remain copy ambiguous.
Distinct donor genes uniquely pointing to the same full target match also remain
copy ambiguous. This reciprocal competition does not establish orthology.

Existing interval candidates, including separate WGD-copy windows, are preserved.
Additional candidates bypass interval extraction and join the deduplicated genome
protein search. The default bound is 20,000 extra queries per target species;
deferred query identities and screen counts are recorded rather than silently
discarded. A deterministic round robin allocates two selections to nearest
references for each common-reference selection, balancing donors within each
group. Strong partial matches have two turns for each absent and ambiguous
representation turn; low-identity matches are considered after these queues
empty. Match identity ranks queries within a bucket. Disabling genome fallback
disables this exploration.
Comparisons and prepared donor proteins must be receipt verified by
the caller before nomination. Every prediction passes genome, genetic-code,
ORF, splice, frameshift, coverage and identity checks again.

Outside-block and genome-only predictions remain proposals unless an intact,
mutually compatible genomic coding path is supported by at least two independent
donor species. Two genes from one donor species supply one species of support;
self evidence does not supply an independent donor species. A donor query with
intact predictions at multiple genomic loci cannot supply unique placement
support. Incompatible paths and ambiguous copies remain visible as proposals.
An admitted model is annotated as **unanchored**. Neither orthology nor recovery
of an expected WGD copy is inferred from this annotation.

`frozen_prediction_cache_key` freezes the original plan and completed worker
receipt hashes, with their models, candidate definitions and genome-query mapping
hashes. `verify_prediction_cache` freshly verifies these contents, original and
current source genome/CDS/GFF hashes, genetic codes, prepared query proteins,
miniprot version and binary hash, and the search parameters. Changes to QC or
acceptance rules are permitted because decisions are recomputed; changes to
prediction inputs fail explicitly. Receipt members cannot escape their producer
directory.

The cache yields only genomic predictor fields and the current unchanged
candidate evidence, stripping old sequence, QC problems, support, status,
selection and model IDs. Old local-search coverage and old genome-search coverage
are separate: an interval formerly accepted but newly unresolved may still need
a genome search. An additional current query with an identical full prepared
protein may reuse a proven old genome search, including a proven empty result.
Interval alignments never authorise genome-search reuse. Query names and evidence
are rebound to the current candidate without borrowing the old synteny placement;
combined genome-query coverage is preserved for subsequent verified reuse.
New producers retain pristine predictor coordinates before terminal completion so
later policy changes cannot inherit completed coordinates as raw predictions.
JSON results are decoded and hashed in bounded chunks. This
avoids a second full JSON string allocation while preserving a fresh content
verification and a final file-identity fence before publication.
