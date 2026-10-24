# Species-tree topology constraints

Concatenated protein and DNA IQ-TREE searches optionally accept an unrooted,
multifurcating Newick constraint. Existing runs remain unconstrained by default.
Set these entrypoint options (or their `GG_GENOME_EVOLUTION_` uppercase overrides):

```bash
species_tree_topology_constraint="/project/constraints/reviewed.nwk"
species_tree_topology_constraint_scope="complete"
undated_species_tree="iqtree_pep"
run_astral_pep=0
run_astral_dna=0
```

Tip labels must exactly match alignment sequence IDs. Empty names, duplicates,
unknown tips, omitted tips in `complete` mode, and constraints without a
nontrivial unrooted split fail before IQ-TREE runs. Set scope `partial` explicitly
to constrain only a subset; other taxa can attach freely. Polytomies remain free.
The final tree is checked for the selected taxa and all required splits before
publication. A stale search restarts IQ-TREE with `-redo`, so an old checkpoint
cannot silently resume a different constraint. The constraint file contents and
scope enter the existing artifact provenance contract. Changing them invalidates
the concatenated tree and its consumers, leaving BUSCO, extraction, alignments
and individual gene trees unchanged.

`species_tree_iqtree_taxa` optionally names a text file of exact alignment IDs,
one per line. It selects rows from the existing concatenation without realigning
or filtering columns. The constraint is validated against that selected set.
This allows a reviewed taxon exclusion while preserving earlier BUSCO/MSA work.
It does not change the species input cohort or the individual BUSCO gene trees.
ASTRAL must be disabled and an IQ-TREE summary selected when using either option;
stale-artifact handling also cannot automatically re-enable ASTRAL in this mode.

Prepare a constraint separately with `nwkit constrain`, or supply a reviewed
Newick. For minimal higher-group assumptions, use a manually reviewed backbone
with `--backbone ncbi_user`; NCBI supplies species-to-group membership.
`--backbone ncbi` imposes NCBI's broader hierarchy, which may differ from the
intended hypothesis. Inspect missing/ambiguous assignments in the NWKIT log;
GeneGalleon's default complete-taxon validation rejects dropped taxa. A full
reviewed species-level backbone can instead use `--backbone user` and exact-ID
parsing (`--species-parser legacy --species-regex '^(.+)$'`). No NCBI lookup is
performed implicitly during inference.

Constraints specify unrooted topology, not the biological root. Configure
`species_tree_rooting` separately, for example `outgroup,Species_A,Species_B`.
The per-search `topology_constraint.json` and `topology_constraint.result.json`
record input hashes, selected taxa and assumed splits. Bootstrap support on an
imposed split is conditional on that assumption, not independent evidence for it.
