# NOTUNG replacement with NWKIT

GeneGalleon no longer invokes or installs NOTUNG. Only its rooting and
reconciliation paths are replaced. GeneRax, GRAMPA, SeqKit and other analysis
programs retain their roles.

## Configuration migration

| Previous setting | Current setting |
| --- | --- |
| `tree_rooting_method=notung` | `tree_rooting_method=reconciliation` |
| `run_notung_reconcil=1` | `run_reconciliation=1` |
| `run_busco_dupaware_notung_root_dna=1` | `run_busco_dupaware_reconciliation_root_dna=1` |
| `run_busco_dupaware_notung_root_pep=1` | `run_busco_dupaware_reconciliation_root_pep=1` |
| `notung_jar`, `NOTUNG_DOWNLOAD_*`, `NOTUNG_ZIP_SHA256` | Remove; no replacement is needed |

Gene-tree rooting still defaults to MAD. Reconciliation rooting uses
`reconciliation_duplication_cost=1.5` and `reconciliation_loss_cost=1` by
default. These are the historical NOTUNG weights. The previously available
NWKIT rooting option used NWKIT's defaults of 1 and 1: explicitly set
`reconciliation_duplication_cost=1` to retain that weighting.

The two weights must be finite and nonnegative, with at least one positive.
Species labels use the configured parser, regular expression or mapping TSV.
Both the rooted species tree and gene-tree internal topology must be binary;
unmatched species and unsupported multifurcations fail rather than being
silently pruned or arbitrarily resolved. Transfer inference is not introduced.

## Root candidates and selection

NWKIT's `root --method reconciliation --candidates-out PATH` exports every
optimal physical edge as a Newick tree collection. Candidates are ordered
canonically by topology and tip names. There is no 1,000-candidate truncation.
The selected root and collection retain the input's pairwise tip distances.
Different ordering from NOTUNG can change the selected tree when roots tie.

For BUSCO trees, the candidate-generation stages write
`busco_reconciliation_dna/<busco>.busco.roots.nwk` and the corresponding
`busco_reconciliation_pep/` files. The existing `run_busco_dupaware_root_dna`
and `run_busco_dupaware_root_pep` stages then select:

1. MAD, if its root edge is an optimal reconciliation edge.
2. Midpoint, if its root edge is an optimal reconciliation edge.
3. The first canonical optimal reconciliation root.

Missing, empty or invalid candidate collections are errors. The selected
`.root.nwk`, `.root.txt`, `.root.tsv` and `.root.pdf` bundle is retained; the
comparison plot still compares MAD and midpoint. Candidate/root generation
uses the configured species mapping and weights.

For gene families, reconciliation rooting writes
`root_candidates/<family>_roots.nwk`. Query-family orthogroup extraction
considers the optimal roots, including when the rooting stage is reused from
cache. Missing or empty collections stop extraction and preserve existing
outputs; they never reduce the search silently to a single selected root.
This replaces the ephemeral NOTUNG directory.

## Reconciliation and downstream consumers

`run_reconciliation` writes `reconciliation/<family>_reconciliation.tsv` using
`nwkit reconcile --event-source lca --unmatched error`. It contains stable clade
identifiers, event assignments and `implied_losses` for each event. Losses are
counted on topological species-tree edges below the gene root, not from the
species root down to the gene root. Leaf losses are zero. NWKIT leaves loss
counts undefined for non-LCA event sources or unmapped events.

Without GeneRax, RADTE consumes this TSV and the same active rooted gene tree
directly. In query2family mode, reconciliation and dating use the extracted
rooted tree and alignment, rather than mixing an original full-family tree
with an extracted alignment.
GeneRax NHX dating remains available. Summary statistics now use
`reconciliation_num_dup`, `reconciliation_num_speciation` and
`reconciliation_num_loss`, replacing NOTUNG-specific `ntg_*` fields. They do
not report transfer or polytomy counts as if those models had been fitted.

Legacy NOTUNG ZIPs are not translated into native TSVs and remain untouched.
Update saved configurations and rerun enabled candidate/reconciliation stages.
New paths and artifact contracts prevent treating legacy ZIPs as NWKIT results.
The root candidate collection is included in orthogroup extraction provenance.

## Dependency deployment and validation

The sibling NWKIT checkout adds optimal-root collection export and LCA loss
columns. A standard container continues to install NWKIT from the moving branch
in `container/source_branches.env`; local dependency edits must be published
there before a remote-source standard build can include them. No source SHA is
pinned into GeneGalleon defaults.

The container build runs `check_nwkit_reconciliation.py` and rejects NWKIT
installations without the required exports or with incorrect fixture counts.

Runtime validation exercises native root collections, BUSCO selection,
reconciliation summaries, RADTE's native TSV route, and failure preservation.
Docker validation does not establish native Apptainer/SIF compatibility.

Validated on 2026-09-11 with the local ARM64 Docker image
`local/genegalleon:nwkit-reconciliation-dev`, derived from
`local/genegalleon:cdskit-stats-dev` with a wheel built from the local NWKIT
checkout, updated validation scripts, and `Notung.jar` removed:

- All 56 required runtime checks passed, including the new reconciliation
  export probe and the existing IQ-TREE worker/CLI agreement check.
- All 12 recorded NOTUNG reference fixtures matched for the complete set of
  optimal root splits and for duplication/loss counts. The native regression
  tests use the recorded data without requiring NOTUNG.
- 30 focused integration tests passed, including query extraction over all
  candidate roots, reuse from cache, extracted-tree reconciliation, RADTE's
  native/GeneRax routes, and failure preservation.
- 351 broader workflow/support tests, 397 declared static-suite tests, and 39
  final container/configuration tests passed. These suites overlap.
- 294 NWKIT root, reconciliation, CLI contract, RADTE and contrast tests passed.
- ShellCheck, changed-Python-file Ruff checks, configuration-schema validation,
  and whitespace checks passed.
- The new export probe rejected the older NWKIT image as expected, preventing
  a build from silently accepting an incompatible dependency.

This validates a locally derived runtime, not a fresh full multi-platform image
build. Native amd64 and Apptainer/SIF execution were not performed. Neither
repository was pushed as part of this change.


### Follow-up review

The commit review additionally covered empty output-path validation, missing
and empty cached root collections, and failures while writing both candidate
and selected-tree outputs. Query-family extraction now rejects absent or empty
candidate collections instead of silently using only the selected root.
NWKIT avoids an unnecessary second copy of each candidate tree, and its export
documentation clarifies that an already matching root keeps its input position.

A Docker test snapshot containing only the reconciliation changes passed 354
NWKIT root, reconciliation, CLI, RADTE and contrast tests. The final export
suite passed all 15 tests, including injected failures during candidate and
selected-tree writes and independently enumerated LCA scores under four cost
settings. Unrelated in-progress shift changes were excluded from that snapshot
and from the reconciliation commit.

The follow-up GeneGalleon run passed 375 focused workflow, configuration,
resource and development-tooling tests in the Docker runtime, including the
native/GeneRax RADTE routes. The declared static suite passed 397 tests, the
required runtime validator passed, and changed-file Ruff, the repository's
configured ShellCheck and whitespace checks passed. These suites overlap.
