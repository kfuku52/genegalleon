# Host-scaffold taxonomy context for HGT

This measures whether the **annotated CDS background of a scaffold** matches
the host's taxonomy. It is independent of shared-neighbor synteny and is neither
an HGT probability nor proof that an assembly is correct. No threshold or
automatic filtering is applied.

## Generation

`gg_genome_annotation_entrypoint.sh` has `run_scaffold_taxonomy=1` by default.
The stage runs only if the existing `species_gff_info` and raw
`species_cds_mmseqs2taxonomy` outputs are present. It does not enable MMseqs2
searches or GFF extraction. Enable `run_collect_gff_info` and
`run_cds_mmseqs2taxonomy` separately when their inputs need generating.
It executes before contamination removal and never uses the cleaned FASTA.

The existing ETE/NCBI taxonomy database resolves the species name uniquely.
Unresolvable/ambiguous host names are errors, not a guessed classification.
For a manual invocation, `scaffold_taxonomy.py --host-taxid` supplies an explicit
host anchor. The pipeline intentionally does not borrow the contamination-removal
target, which can be a broad taxon rather than the actual host.

Outputs in `output/species_scaffold_taxonomy/`:

- `<species>_gene_taxonomy.tsv`: one CDS identifier × rank, with `species`,
  `gene_id`, `scaffold`, `locus_id`, `count_unit`, `rank`, `host_taxid`, `label`.
- `<species>_scaffold_taxonomy.tsv`: one scaffold × rank with `species`,
  `scaffold`, `rank`, `host_taxid`, and the metrics below.

Ranks: domain (including NCBI superkingdom), phylum, class, order, family,
genus, species. Labels are `compatible`, `incompatible`, or `unresolved`.
A coarser LCA cannot establish compatibility at a finer rank. For example,
Eukaryota alone is unresolved at phylum. Missing host ranks, missing assignments,
and unknown taxids are also unresolved. A bacterial LCA without a phylum is
incompatible at domain but unresolved at phylum; consult both ranks.

The GFF's explicit transcript-to-gene parent chain groups isoforms into loci;
conflicting isoform classifications make the locus unresolved. Without that
relationship the counting unit is the input CDS ID (`count_unit=cds_id`), not
an inferred locus. Such inputs may overcount isoforms; use one representative
CDS per locus or provide the complete GFF. Missing coordinates and trans-spliced
genes are excluded from this single-scaffold measurement, not assigned to a
guessed scaffold. The denominator is **coordinate-mapped input CDS loci**, not
all DNA, all GFF features, or unannotated genes.

## Metrics

| Suffix | Meaning |
| --- | --- |
| `total_count` | All counted loci, including unresolved |
| `compatible_count` | Loci matching the host taxid at this rank |
| `incompatible_count` | Loci assigned a different taxid at this rank |
| `unresolved_count` | Loci without a usable rank comparison |
| `cds_id_count` | Subset of total lacking explicit GFF locus mapping; counted by CDS ID |
| `classified_fraction` | (compatible + incompatible) / total |
| `compatible_fraction` | compatible / (compatible + incompatible) |
| `compatible_all_fraction` | compatible / total |

Zero denominators produce blank fractions, not zero support. A high compatible
fraction with few classified background loci is weak information; inspect counts
and coverage alongside fractions. Short scaffolds and database representation
bias can limit interpretation. Close-relative transfers may be taxonomically
indistinguishable at broad ranks. Contamination or chimeric assembly cannot be
ruled out from this metric alone.

## HGT integration

`gg_hgt_core.sh` automatically loads these outputs during HGT evaluation.
Gene candidates retain their existing row unit and get `host_scaffold_*` columns.
`host_scaffold_<rank>_*` measures the whole scaffold;
`host_scaffold_background_<rank>_*` excludes the union of **all HGT candidate
loci across orthogroups** on that scaffold, not just the focal gene. All isoforms
of an excluded locus are removed. This conservatively accommodates multi-gene
transfers, but the exclusion set can also include donor-side candidate descendants;
it describes candidate-free background, not a reconstructed transferred block.

Branch candidates use only candidate genes whose species are descended from the
GeneRax `Y@donor@recipient` recipient label in the species tree. Counts are pooled
over distinct `(species, scaffold)` pairs, then fractions are recalculated;
multiple genes on one scaffold do not multiply its weight. Large scaffolds have
more weight than small ones. This is context of sampled present-day descendants,
not reconstructed ancestral scaffold structure. Unresolved recipient labels do
not fall back to mixing donor and recipient descendants.

`host_scaffold_status=measured` means data were joined, **not** that they passed
any quality threshold. Empty backgrounds have zero counts and blank fractions;
missing tables or missing recipient mapping retain explicit status and blank
measurements. Orthogroup summary is unchanged; it does not sum repeated scaffold
context across independent event rows. `output/hgt/README.md` documents every
new branch/gene column automatically. Artifact provenance includes the taxonomy
database, GFF, raw taxonomy, helper code, and species tree at their consuming stages.
