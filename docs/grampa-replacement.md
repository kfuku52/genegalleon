# GRAMPA Replacement with NWKIT

GeneGalleon's BUSCO DNA, BUSCO protein and orthogroup polyploidy stages now
invoke `nwkit mul-reconcile`, not `grampa.py`. The native engine searches
GRAMPA-style auto-/allopolyploid MUL-tree hypotheses and minimizes exact
duplication-plus-loss parsimony, including species-root-to-gene-root losses.
It is separate from the [WGD count/synteny/Ks analysis](wgd-ssd.md).

## Configuration and Outputs

Existing `run_busco_dupaware_grampa_dna`, `run_busco_dupaware_grampa_pep`,
`run_orthogroup_grampa` and `grampa_h1` settings retain their names and values.
An empty `grampa_h1` still disables these stages. H1 species sets and numeric
postorder internal-node selectors retain their meaning, including multiple
space-separated H1 selectors. Gene/species trees
must be rooted and binary. Input preparation now uses a tree parser and the
configured species parser/regex/TSV; malformed or unmatched tips fail instead
of being skipped. The old underscore-to-hyphen naming in prepared inputs and
summary columns is retained, with normalization collisions rejected.

Output directories and `grampa_summary.tsv` columns remain unchanged, together
with `grampa_out.txt`, `grampa_det.txt`, `grampa_checknums.txt` and
`best_mul_tree.nwk`. New `nwkit_mul_reconcile.json` records the native method,
all tied hypotheses, and limitations. The detailed `maps` field retains
annotated Newick; additional `node.maps` contains auditable JSON node mappings.
Check-table group counts now mean
ambiguous tips rather than GRAMPA's collapsed groups. Gene/tree IDs are
run-local, and the first minimum-score candidate supplies the summary.
The prepared inputs and literal filename inventory are published with the
results. Each run uses an isolated temporary directory; all nine related
outputs use the existing recoverable bundle publisher. Failed computation or
publication preserves previous results. Empty inputs with existing results
fail rather than record those results as a new analysis.

No candidate-dependent gene filtering or map truncation is performed. Resource
limits fail the whole analysis rather than compare hypotheses on different
gene subsets. All three artifact contracts record the engine, NWKIT identity,
input/summary adapters, species settings and the complete output bundle, so
old GRAMPA caches are stale;
use the existing stale-artifact policy to authorize rebuilding them. Existing
results and curated inputs are not automatically converted or deleted.

## Interpretation and Deployment

This is D+L parsimony, not a calibrated test or donor posterior. Fixed trees,
missing genes, ILS, HGT and gene-tree uncertainty can affect the ranking.
The exact solver can improve upon the original grouping heuristics: a verified
whole-root auto fixture scores 13 versus the original program's 15, with 13
confirmed by complete enumeration using the original LCA scoring function.
The author's 25-tree manual dataset matches all 10 tested candidate scores.

Standard images obtain the command from the moving NWKIT branch declared in
`container/source_branches.env`; an installed NWKIT without the required
exports cannot run these stages. No upstream default version or SHA
has been pinned. Both architecture-specific Conda manifests and required-tool
lists no longer request GRAMPA. Runtime validation requires a real native MUL
reconciliation, score reconstruction from node mappings, and model metadata;
an installed NWKIT lacking this command fails the build.

See the [NWKIT command guide](https://github.com/kfuku52/nwkit/blob/master/docs/guides/MUL_RECONCILE.md)
and its [validation report](https://github.com/kfuku52/nwkit/blob/master/docs/validation/MUL_RECONCILE_VALIDATION.md).
The workflow's native MUL integration tests exercise the installed dependency
and recoverable publication through the normal runtime checker. See
[development validation](development-and-tests.md) for Docker and SIF commands;
a passing Docker run alone does not establish SIF compatibility.
