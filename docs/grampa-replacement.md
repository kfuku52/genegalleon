# GRAMPA Replacement with NWKIT

GeneGalleon's BUSCO DNA, BUSCO protein and orthogroup polyploidy stages now
invoke `nwkit mul-reconcile`, not `grampa.py`. The native engine searches
GRAMPA-style auto-/allopolyploid MUL-tree hypotheses and minimizes exact
duplication-plus-loss parsimony, including species-root-to-gene-root losses.
It is separate from the [WGD count/synteny/Ks analysis](wgd-ssd.md).
That analysis can optionally [join per-family MUL node diagnostics](wgd-ssd.md#optional-mul-node-diagnostics)
without relabeling WGD/SSD evidence or reusing this genome-stage result bundle.

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

### Optional DL + ILS Comparison

Setting `grampa_locus_model` to an explicit NWKIT locus-model JSON enables
an additional **experimental** `--score-model locus-mc` analysis for each
enabled BUSCO/orthogroup stage. Empty (the default) keeps D+L only. The JSON
supplies every rate, Ne, detection probability, selection limit and budget;
GeneGalleon does not invent these settings. NWKIT supports explicit
`detection-rb` and `hybrid-rb` conditional integration as well as the older
histogram modes; see its [locus guide](https://github.com/kfuku52/nwkit/blob/master/docs/guides/MUL_LOCUS_MC.md).

Also supply `grampa_locus_species_tree`: a rooted binary ultrametric tree
with branch lengths explicitly in **generations**, the same species and
rooted topology as the workflow species tree. Ordinary workflow dating
units are not assumed to be generations and are never automatically
converted. The adapter matches child order to preserve numeric H1/H2 IDs.
Model detection keys use the original species names; the adapter normalizes
them together with tree labels, rejects collisions, and preserves both input
files. Keep model/tree inputs outside the result directory.

`grampa_h1` must identify one supported non-root clade for this mode.
`grampa_locus_h2` optionally restricts its donors (empty searches all supported
parents). `grampa_locus_bootstrap=0` disables null calibration; a positive
value is replicates per null point. `grampa_locus_null_calibration` accepts
`plug-in` (default) or `grid-supremum` (requires positive replicates). These
parameters also support `GG_GENOME_EVOLUTION_` scoped environment overrides.

The separate `locus/` result directory contains `scores.tsv`, `families.tsv`,
`checks.tsv`, `results.json`, normalized `input_model.json` and
`species_tree.nwk`, plus `null_search.tsv` when requested. These are **not**
GRAMPA-format scores or inferred biological point trees. The unchanged D+L
summary and experimental directory publish in one recoverable transaction;
failure in either analysis/publication preserves the entire previous bundle.
All three stage contracts hash the original model/tree, scope, calibration
settings and result directory. Disabling the option retains previous locus
results with a warning; they are not recorded as part of the D+L-only run.

No families are silently removed to fit the model's observation limit. In
particular, the existing orthogroup size defaults (5-50) are incompatible
with a model declaring at most four observed tips; choose a scientifically
justified input selection explicitly. Additional workflow curation and
IQ-TREE/rooting error are not reproduced by the topology-only null simulator.
Its P-values do not calibrate those real-data processes. Wide MC bounds,
missing support and resource failures remain explicit; GeneGalleon never
converts these results to an automatic WGD call. This integration does not
change the default scientific model or its thresholds.

### D+L Deployment

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
