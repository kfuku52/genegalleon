# Host-scaffold taxonomy context for HGT

See [GFF validation](gff-validation.md) for coordinate validity, phase status,
translation exceptions, and repairing pre-existing project inputs.

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

The host anchor first uses a unique NCBI taxid explicitly recorded by the input
GFF's `##species` directive or source `region` `Dbxref=taxon:` attributes.
Conflicting source taxids are errors. When no source taxid is recorded, the
existing ETE/NCBI taxonomy database must resolve the species name uniquely.
An explicit `Genus_sp_unknown` label may resolve to one verified genus taxid;
species-level composition remains unresolved. Other unresolvable/ambiguous
host names are errors, not a guessed classification.
For a manual invocation, `scaffold_taxonomy.py --host-taxid` supplies an explicit
host anchor. For a scheduled single-species annotation task, set
`GG_GENOME_ANNOTATION_SCAFFOLD_HOST_TAXID` to that NCBI TaxID; the value is
recorded in the stage provenance and overrides an ambiguous name. The pipeline
intentionally does not borrow the contamination-removal
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
any quality threshold. Candidates whose species cannot be mapped to a tree tip
are counted in `host_scaffold_unresolved_taxon_gene_count`; a branch with such
candidates is `partial` even when its known recipient genes all have coordinates.
Empty backgrounds have zero counts and blank fractions;
missing tables or missing recipient mapping retain explicit status and blank
measurements. Orthogroup summary is unchanged; it does not sum repeated scaffold
context across independent event rows. `output/hgt/README.md` documents every
new branch/gene column automatically. Artifact provenance includes the taxonomy
database, GFF, raw taxonomy, helper code, and species tree at their consuming stages.

Input gene-taxonomy tables must have every rank exactly once per gene, a
consistent scaffold/locus/counting unit across ranks, consistent host taxids,
and reconciled isoform labels. Corrupt or incomplete tables raise an error;
missing rank rows are not silently omitted from the denominator. Explicit
GFF loci and fallback CDS IDs use separate counting namespaces.
Numeric internal species-tree names are retained verbatim (including leading
zeros and long identifiers). Duplicate names or space/underscore alias
collisions are rejected rather than assigned to an arbitrary branch.

## Event-resolved donor and recipient context

With `run_hgt_candidate_summary=1`, `gg_gene_summary` also writes
`hgt_transfer_events.tsv` and `hgt_transfer_event_genes.tsv` in the selected
HGT output directory. Existing branch/gene/orthogroup tables retain their
schemas and meaning. No threshold or automatic filtering is introduced.

The event table has one row per family × gene-tree branch × transfer-token
position, including unresolved events. It matches `Y@donor@recipient` and the
branch's exact descendant-gene set to one GeneRax XML `branchingOut` / child
`transferBack` pair. Missing or ambiguous matches remain unresolved with
blank measurements. Raw and ZIP-backed reconciliations use the same logical
gene-family store reader. XML species-tree branch labels and descendant sets
are authoritative; a supplied species tree must agree with their named clades.

Each side follows its continuation lineage and excludes later `transferBack`
descendants and leaves outside the original species branch from aggregation.
The event-gene table retains these excluded links with their reasons. It joins
existing gene-level scaffold measurements by family and gene, verifying the
XML leaf species against the gene's own host. Best-hit-derived `donor_*`
taxonomy columns in older tables are never used to assign transfer roles.

`donor_host_scaffold_*` and `recipient_host_scaffold_*` pool counts over unique
`(species, scaffold)` pairs and recalculate fractions; repeated copies on the
same scaffold do not multiply its background. The existing candidate-free
background exclusion set is retained. Per-gene measurements, locus identifiers,
count units, and available auxiliary evidence are preserved in the link table.
Missing measurements stay blank. `measured` describes availability, not a pass.

Terminal species branches use `extant_terminal_genome` evidence; internal
branches use `extant_descendant_proxy`. Neither reconstructs ancestral scaffold
structure or proves physical integration. Event counts, unique genes, and
unique families must be reported separately. The event's UFBoot is read only
from the same `stat_branch` family/branch/descendant-set match; terminal branches
and missing values remain unavailable, with no generic-support substitution.

The native species-tree summary plot consumes the event table when available.
For a project-filtered figure, pass its selected rows to
`plot_hgt_summary.py --transfer_event_tsv PATH` alongside the existing
branch/gene overview inputs. Each row counts as one event; unique families are
counted separately. The drawing and direction/count conventions are unchanged.

To highlight an observed binary tip trait, also pass
`--species_trait PATH --transfer_tree_highlight_trait COLUMN`. Positive tip
labels and incoming species branches are orange; an internal branch is colored
only when all descendant tips have an observed `1`. Missing descendants prevent
highlighting. This is a display rule for homogeneous clades, not inferred
ancestral states. Links entering a highlighted recipient branch use the same
orange and render above other links. Each direction has one constant-width arrow;
reciprocal directions use separate curves and recipient-specific colors.
Arrow opacity defaults to `0.55` and can be changed with
`--transfer_arrow_alpha` (range `0` to `1`). Event counts, selection and distance
measurements are unchanged.
Ordinary arrows use one blue color; highlighted arrows use the trait color.
Endpoint-node distance does not control arrow color and has no colorbar.

## Category-1 focused results

When HGT summaries are enabled, `run_hgt_trait_focus=1` (the default) also writes
`hgt/trait_focus/index.tsv` and result bundles for each eligible trait, each
observed category-1 tip, and each all-positive internal recipient branch.
Binary and declared categorical traits are eligible; schema-free columns must
contain only observed 0/1 values. Declared continuous numeric/text traits and
observation/quality columns are excluded. Categorical values other than 1 are
other categories, not evidence of biological absence; missing stays unknown.
Hash-bound trait schemas and observation contracts are validated.

The exporter uses native event/context tables by default. To return a previously
filtered project cohort, set `hgt_summary_focus_event_tsv` and optionally
`hgt_summary_focus_event_gene_tsv` to its event and enriched link TSVs. Existing
support, direction, scaffold, product-name and quality columns are preserved;
no additional scientific threshold or organism-specific filter is introduced.
`run_hgt_trait_focus=0` disables these outputs; disabling all HGT summaries does
not start an HGT analysis merely because the focus default is on.

Each bundle contains events, branch pairs, event-gene links, donor/recipient gene
tables and counts. Internal donor/recipient IDs have their full species-tree tip
lists alongside them. Per-tip `direct_events.tsv` and
`ancestral_recipient_events.tsv` separate direct transfers from qualifying
ancestor-branch context. An ancestral event can occur in multiple tip reports;
the aggregate counts it once by event ID. Do not sum per-tip totals as independent
acquisitions. Mixed/missing recipient clades remain outside the focused cohort
with explicit reasons. Empty targets are retained as zero-result reports.

With `run_hgt_summary_plots=1`, each trait aggregate's `plots/` contains three
single-page PDFs: `filtering_flow.pdf`, `orthogroup_species_distribution.pdf`
and `donor_recipient_counts.pdf`. The distribution includes an existing protein
product/annotation label, its source gene and evidence basis in the accompanying
TSV. Counts distinguish event IDs, gene-tree tip copies and orthogroups. All
bundles retain their directed edge TSVs; focused species-tree and per-recipient
PDFs are not produced.

Each trait aggregate also has a `tree_plot/` folder with one native gene-tree PDF
per qualifying orthogroup. Orange diamonds mark the exact transfer nodes;
labels give the matched `support_generax_ufboot`, including internal nodes.
This presentation subset requires UFBoot >=90 and at least one retained
event-linked gene on **each** side with candidate-free class scaffold background:
>=10 classified units, >=50% classification coverage and >=90% host compatibility.
Those recipient genes have orange tip labels. The parent focused tables are
unchanged. `tree_plot/event_node_audit.tsv` retains every requested event with
selection/withholding reasons; terminal nodes, missing support and unresolved
branch/token matches do not qualify. Background context is not conserved gene
order or proof of physical integration. Standalone `focus_hgt_traits.py` calls
enable this folder by supplying `--gene_family_root` with existing raw/ZIP-backed
family outputs; the workflow supplies this path automatically when plots run.

Each individual PDF has two pages. Page 1 replays the panels and arguments
recorded by `gg_gene_evolution` (including domain, structure, alignment,
localization and available similarity panels). Older results use their saved
tree-plot parameter/input provenance; unavailable optional inputs remain missing.
`renderer_settings.json` records the settings source and input availability.
Page 2 shows the exact HGT node/UFBoot and existing GFF neighborhoods in two
columns: donor descendants in blue, recipient descendants in orange, and nearby annotations
in gray. It draws at most three distinct genes per side on one page, with
shown/total/omitted counts. Eligible event-linked genes with failing or unavailable
scaffold evidence can appear with explicit status and pale hatched focal blocks;
this does not change event selection or the individually supported gene counts.
Display priority is passing scaffold support, available GFF, background coverage,
host compatibility, then gene ID. Repeated side/gene links draw one track while
retaining each exact event reference in the audit. All genomic tracks share
one linear kb axis centered on the focal feature midpoint. The shared window
includes each focal feature plus at least 20 kb flanks; asterisks mark neighboring
features extending beyond the window. CDS and available UTR
blocks retain genomic intron lengths; missing structure and trans-splicing are
explicitly unavailable. Neighbor annotations do not establish their taxonomy or
conserved gene order. `context_gene_audit.tsv` records displayed and omitted
genes, selection rank/reason, measured evidence and count units, complete gene
IDs, structures, exact event branches and shared axis limits. The workflow supplies
existing `species_gff_info`; standalone calls use `--gff_info_root`.

Under each genomic track, page 2 lists the focal gene and every displayed
neighbor with its own protein product, Swiss-Prot best-hit organism/accession,
and kingdom, phylum, class, order, family and genus. The table uses the same
left-to-right neighbor numbers as the track; focal rows are colored by side.
GFF products and best-hit protein-name predictions have explicit labels.
Missing hits, products and ranks remain unavailable; a neighbor never inherits
the focal gene's annotation. These hit ranks describe the annotation hit and
do not identify the modeled GeneRax donor or establish host scaffold support.

Standalone calls accept `--context_annotations_tsv`, an existing normalized
per-gene TSV with unique full `gene_id`, `orthogroup`, `protein_product_name`
(own GFF product), `protein_product_status`, `besthit_accession`,
`besthit_organism`, and `besthit_{kingdom,phylum,class,order,family,genus}` columns.
Optional fields include `swissprot_best_hit_protein_name`, hit taxid, GFF/hit/taxonomy
source paths and SHA-256 hashes, product feature IDs and selected-transcript/locus
scope. Join each neighbor to its **own** family and exact leaf hit, then resolve
that same hit's taxid using the existing taxonomy database. This input adds
annotation only and does not request a sequence search or taxonomy update.
Without it, exact focal-family leaf hit annotations are used when available;
unavailable neighbor annotations and taxonomic ranks stay missing.
`context_annotation_audit.tsv` records every displayed gene, its annotation
inputs and exact context/event references; the input TSV hash is in the manifest.
The one-page canvas expands to fit complete annotations while retaining the
three-gene-per-side display limit and the same linear genomic scale per track.

`hgt_summary_focus_filter_audit_tsv` (standalone `--filter_audit_tsv`) accepts an
optional existing project event-level direction/UFBoot audit, including gzip TSV.
When supplied, the filtering-flow PDF also shows the upstream modeled,
non-Insecta-to-Insecta and verified UFBoot >=90 counts. Without it, upstream
counts are not inferred from an already filtered input table. No sequence,
phylogenetic or annotation-search analysis is run for these figures.

`manifest.json` binds inputs and every output to SHA-256 hashes. Focus bundles
are replaced as a managed unit; invalid inputs do not publish a partial bundle.
