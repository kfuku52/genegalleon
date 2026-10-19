# Optional HGT direction filters

Direction filtering is **off by default**. It affects the trait-focused view,
not the saved input cohort, reconciliation, scientific inference, or global
candidate tables. Enable it explicitly in `gg_gene_summary_entrypoint.sh`:

```bash
hgt_summary_focus_direction_filter="custom"
hgt_summary_focus_direction_rules="/absolute/path/direction-rules.json"
```

A versioned JSON file defines the two endpoints independently. For example,
select fungal donors and let the existing category-1 trait select recipients:

```json
{
  "schema_version": 1,
  "donor": {"kind": "taxonomy", "rank": "kingdom", "values": ["Fungi"], "match": "within"},
  "recipient": {"kind": "any"}
}
```

| Selector | Fields | Meaning |
| --- | --- | --- |
| `taxonomy` | `rank`, `values`, optional `match` | Match the exact saved host-species rank column; use `species_rank` for the species rank. `hgt_summary_focus_species_taxonomy` supplies the table. No network lookup or best-hit substitution. |
| `tips` | `values`, optional `match` | Classify every full-tree descendant by membership in the explicit species-tip list. All-positive internal branches qualify automatically. |
| `branches` | `values`, optional `match` | Select exact named branches of the analysis tree, including terminal or internal incoming branches. Does not expand the list to ancestors or descendants. |
| `any` | no other fields | Leave this endpoint unconstrained; mapping must still be resolved. |

`match` defaults to `within`; `outside` requires all descendants to be outside
that endpoint's taxonomy/tip set. Mixed, unknown, empty and unmapped branches
are withheld, including for exclusion rules. A missing rank is unknown, never
an automatic nonmember. Multi-valued rank cells use the saved JSON rank array.
Names are exact; duplicate/unknown fields and ambiguous or absent tip/branch
identities stop the run. Explicit tip/branch rules require no taxonomy table.
Supplying rules with an off filter is an error, not an implicit opt-in.

Both endpoints must be specified. For a plant-only recipient set, replace the
`any` selector with `{"kind":"tips","values":["Species_one","Species_two"]}`
using the exact tree labels. A shared ancestor whose complete descendant set is
in that list qualifies. This is branch selection, not ancestral trait-state
reconstruction. Category-1 trait selection remains a separate subsequent stage.
The same ancestral event counts once by event ID in the aggregate; its per-tip
context reports must not be summed as independent acquisitions.

The existing Pfam, protein-length and scaffold evidence rules are unchanged.
They remain separate from direction and trait membership. Mixed/missing evidence
is not converted to a passing score, and no new support threshold is imposed.

To render only the focused reports, explicitly set
`hgt_summary_plot_scope="focused"` with `run_hgt_summary_plots=1` and
`run_hgt_trait_focus=1`. This produces the focused tables, three summary figures
and qualifying native family/context PDFs, while omitting the all-candidate
summary and family PDFs. The default `all` retains the original rendering scope.
Global candidate tables remain available in either scope. Raw/ZIP gene-family
results, alignments and reconciliations are reused; no new inference is run.

## Identity and audit

The full rooted descendant-clade sets of the reconciliation and analysis trees
must form a unique bijection. Internal labels may differ; missing species,
topology differences and ambiguous unary clades are not repaired by pruning or
name guessing. Original `generax_transfer`, GeneRax endpoint labels and event
IDs are retained. Separate `analysis_donor_node` / `analysis_recipient_node`
columns identify the verified analysis-tree endpoints. Filters, trait branches
and display edges use those IDs; gene-tree event matching retains original IDs.

`direction_rules.json`, `direction_event_audit.tsv`,
`direction_species_branches.tsv` and the focused manifest record the selectors,
full descendant labels, decisions, source hashes and surviving cohort. Changing
the rule file invalidates reuse, including when plots are disabled. Mutation
during generation fails before result publication. The filtering-flow figure
labels this as `Configured donor → recipient`, separate from the pair and trait
stages. Preserve original input tables and unfiltered reports.

`any` remains a compatible spelling of `off`. The explicit
`non_arthropoda_to_insecta` and `non_arthropoda_to_arthropoda` values remain
compatibility presets with their original rules and audit columns; general
projects should use `custom`. No project species or fungal assumption is a
GeneGalleon default. The standalone `focus_hgt_traits.py` CLI exposes the same
`--direction_filter` and `--direction_rules` options.
