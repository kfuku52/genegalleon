# Shared NWKIT and CDSKIT functionality

GeneGalleon keeps workflow decisions, source evidence, stable identifiers and
artifact publication. Generic tree and sequence operations live in NWKIT and
CDSKIT. The compatibility modules in `workflow/support/` preserve existing
commands and import names; they require the corresponding upstream APIs.

| Audit group | Library implementation | GeneGalleon consumer |
| --- | --- | --- |
| N01 topology export | NWKIT `tool_trees.topology_newick` | `write_newick_topology.py` |
| N02 support/root transfer | NWKIT `branch_features.map_internal_support_by_split`, `transfer.transfer_properties`, `root.transfer_root_with_taxon_mode` | `orthogroup_statistics.py` |
| N03 clade/CI/CAFE correspondence | NWKIT `numbered_tree.numbered_clade_keys`, `tree_features.projected_clade_candidates` | `tree_topology_helpers.R`, `cafe_branch_map.py` |
| N04 parent/LCA/distance index | NWKIT `parent_index.ParentIndex`, `tree_features.patristic_distances` | `query_gene_orthologs.py`, `busco_reference_quality.py` |
| N05 extended species names | NWKIT `SpeciesParser("extended-taxonomic")`, `taxonomic_labels` | `species_labeling.py`, common species-parser configuration |
| N06 ASTRAL annotations | NWKIT `astral`, `tree_formats`, `convert` | `extract_astral_support_label.py` |
| N07 plain tool trees | NWKIT `tool_trees`, `convert --comments drop` | `prepare_generax_species_tree.py`, `prepare_generax_ufboot_target.py` |
| N08 recPhyloXML | NWKIT `recphyloxml_records.read_event_records` | `reconciliation_database.py` |
| N09 branch features | NWKIT `branch_features` | `orthogroup_statistics.py` |
| N10 eligible roots | NWKIT `root_candidates.select_eligible_root` | `species_tree_guided_gene_tree_rooting.py` |
| N11 calibration units | NWKIT `calibration_scale` | `mcmctree_time_scale.py` |
| N12 tree rendering | NWKIT `dated_presentation`, optional `r_render` backend | `dated_tree_presentation.py`, `nwk2pdf.r` |
| N13 trait preparation | NWKIT public helpers in `contrast` and `ordinary_regression` | `species_tree_pgls.py` |
| N14 constraint splits | NWKIT `tree_features.unrooted_splits`/`verify_constraint`, CDSKIT exact selection | `species_tree_constraint.py` |
| N15 foreground clades | NWKIT `tree_features.foreground_components` | `resolve_binary_trait_foregrounds.py` |
| N16 chain diagnostics | NWKIT `posterior_diagnostics` with R `posterior` | MCMCtree diagnostic R wrappers |
| C01 codon meaning | CDSKIT `codonutil`, `translate` | CDS normalization/catalog/resolution and synteny dS inputs |
| C02 strict FASTA | CDSKIT `fasta_contract`, `validate --profile` | `fasta_sequence_contract.py` |
| C03 GFF/CDS geometry | CDSKIT `coding_gff_stream`, `gff_attributes`, `gff_feature_structure`, `cds_features`, `transcript_coordinates` | annotation, CDS catalog and transcript UTR adapters |
| C04 streaming selection | CDSKIT `sequence_selection` | `fasta_sequence_store.py`, `species_tree_constraint.py` |
| C05 shared labels | CDSKIT `label_mapping`, NWKIT `rename` | `shorten_fasta_newick_names.py` and its shell compatibility entrypoint |
| C06 model resolution | CDSKIT `latest_model`, `localize_models` | `resolve_cdskit_localize_model.py` |
| C07 coordinate mapping | CDSKIT `coordinates`, R coordinate resources | treevis alignment/intron data |
| C08 protein/junction similarity | CDSKIT `protein_similarity`, `coordinates.pairwise_identity` | gene-model pair scoring and treevis identity data |
| C09 alignment metrics | CDSKIT `alignment_metrics` | orthogroup lengths and IQ-TREE/mapnh F3X4 frequencies |

Extended species parsing is opt-in in NWKIT; its `legacy` and `taxonomic`
presets retain their behavior. GeneGalleon's unset parser default is
`extended-taxonomic`; an explicitly configured parser or mapping still wins.
Cultivar groups, subgroups and interspecific hybrid names remain distinct IDs.
Host-side cohort acquisition handles frozen species IDs and exact filename
bindings without interpreting taxonomy. Metadata-only provenance operations
do not import the sequence-analysis libraries; scientific parsing and explicit
FASTA validation use the corresponding runtime APIs.

Support bounds, partial-label rejection and root-edge agreement remain strict.
Raw XML evidence preserves event order and incomplete transfer observations;
it is distinct from NWKIT's resolved reconciliation document, which requires
transfer pairing. Database IDs and schemas do not adopt NWKIT traversal IDs.

Codon translation does not infer complete-CDS termination. Context-dependent
codes 27/28/31 require explicit terminal context. Conserved isoform admission,
source-specific corrections, statistical testing families and thresholds stay
in GeneGalleon. The strict FASTA profile also checks ungapped codon lengths.

Coordinate mapping reports unique/ambiguous/unmatched correspondence. The
legacy alignment display explicitly uses the leftmost complete map; callers
can request unresolved output or rejection instead. Intron position, phase,
missingness and ancestry evidence remain separate.

R coordinate helpers are packaged with CDSKIT; the R `ggtree` renderer and
MCMCtree `posterior` diagnostics are optional NWKIT backends. GeneGalleon's
container supplies these R dependencies. Docker validation does not establish
SIF compatibility. Use the existing [development check entrypoint](development-and-tests.md)
for affected callers and serialization, and benchmark changed algorithms with
the same inputs and output checks.
