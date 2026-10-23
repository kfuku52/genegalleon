# Reviewed experimental CDS without genome coordinates

Use `cloned_cds_mode=sequence` for species whose input contains only CDS or
transcriptome sequences. The default search and ambiguity handling remain
unchanged. An explicit `cloned_cds_sequence_groups` JSON manifest column can
resolve reviewed clone competition or admit an experimentally validated CDS
with no supported hit in that sequence catalog.

```json
{
  "candidate_TPS": {
    "clone_ids": ["clone1", "clone2"],
    "source_cds_ids": ["old_fragment"],
    "unplaced": false
  }
}
```

The first `clone_ids` member is the reviewed canonical representative. The
remaining sequences are retained verbatim in `alternative_cloned_cds.fa` and
linked to that representative in `audit.json`. This ordering is an explicit
input choice, not evidence of different function or a confirmed genomic locus.
The candidate group remains genomically unresolved. All members must match
exactly the declared source CDS set under the existing search thresholds.
Those source sequences are archived and replaced once. Alternative clones do
not become extra genes in the canonical set; downstream analyses may include
the archived alternatives with their candidate-group labels.

To include every experimental sequence in phylogeny and localization inputs,
set `"analysis_policy": "all_clones"` in that group. The default
`"representative_only"` preserves the behavior above. With `all_clones`, each
clone receives a distinct analysis ID and remains in the formatted species CDS;
the reviewed source fragments are still archived and removed only once.
Non-first members have status `retain_unresolved_clone` and also remain in the
alternative archive. The first member is only the curation representative,
without priority over other members in downstream analysis.

All reviewed groups record `copy_number_status=unresolved` and a null
`confirmed_gene_copy_count`. The receipt-bound `sequence_group_members.tsv`
links clone IDs, formatted analysis IDs, candidate groups and inclusion status.
Carry these labels into trees and localization results. Two retained sequences
do not establish two gene copies: do not interpret sequence/FASTA counts as
confirmed genomic copy numbers or infer a duplication event from this group
without independent locus evidence. This option does not distinguish alleles
from recent paralogs. Changing it requires a new frozen task plan.

For a validated sequence without a supported catalog match:

```json
{
  "unplaced_TPS": {
    "clone_ids": ["experimental_TPS"],
    "source_cds_ids": [],
    "unplaced": true
  }
}
```

This adds one canonical experimental CDS and preserves all existing source
sequences. It requires a completed uncapped search with no supported match.
The audit records `no_supported_match_in_sequence_catalog`, null genomic coordinates and
`not_assessed_no_genome`. Missing assembly regions, expression differences and
divergent reference individuals are possible explanations, not established
origins. This sequence-only route rejects paired GFF/GBFF/genome inputs and
does not fabricate GFF coordinates.

Group names use letters, digits, `_`, `.` and `-`, starting with a letter.
Every selected clone must belong to exactly one group; source CDS ownership
cannot overlap groups. Groups require `clone_ids`, `source_cds_ids` and boolean
`unplaced`; the only optional key is `analysis_policy`.
Explicit `cloned_cds_ids` are required.
Unplaced groups cannot use the longer-only rule because no comparator exists.
Missing source IDs, changed matches, unsafe ambiguous alignments, capped
searches and failed longer-only comparisons stop before publication.

The derived input archives all selected experimental CDS, the reviewed group
specification, removed source sequences and any alternatives. Every artifact
enters the ordinary native receipt. Changing group membership or representative
order requires a new frozen task plan; checksum validation and stage reuse
remain unchanged. For a reviewed projection onto an available reference genome,
use [experimental genomic projections](experimental-cloned-cds.md).
