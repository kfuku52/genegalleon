# Pairwise genome synteny

`gg_genome_evolution` can infer gene-based synteny between two species with
JCVI/MCscan and DIAMOND, then draw a chromosome ribbon plot and a dotplot.
The ribbon plot serves the chromosome-comparison use case of a GENESPACE
riparian plot; it does not construct GENESPACE synteny-constrained orthogroups
or pan-gene sets. No species tree, BUSCO database or OrthoFinder run is required.

## Inputs and first run

Provide one matching annotation release per species in `workspace/input/species_gff`
and protein FASTA files in `species_protein`, or CDS FASTA files in `species_cds`.
Use the [species naming conventions](input-conventions.md). `auto` prefers a
protein source for each species and translates CDS where no protein source exists.
CDS translation uses `genetic_code` and the existing per-species
`species_genetic_code/species_genetic_code.tsv` overrides. Non-triplet CDS,
internal stops, duplicate FASTA IDs, ambiguous ID aliases and conflicting
annotation coordinates are errors.

Create `workspace/input/synteny_pairs.tsv` with these required columns:

```tsv
analysis_id	target_species	query_species
triphyophyllum_ancistrocladus	Triphyophyllum_peltatum	Ancistrocladus_abbreviatus
```

The same row is supplied in `workspace/input/synteny_pairs.example.tsv`.
`analysis_id` must be unique. Both species must differ. Each row is one comparison;
rows run sequentially inside this single-task workflow using its allocated CPUs.

```bash
GG_GENOME_EVOLUTION_GENOME_EVOLUTION_MODE=synteny \
  bash workflow/gg_genome_evolution_entrypoint.sh
```

`genome_evolution_mode=synteny` enables this stage and exits before species-tree
and orthogroup setup. Existing outputs from those stages are preserved. In
ordinary `all` mode, set `run_pairwise_synteny=1` to add the stage to the pipeline;
its default is `0`. Scheduler resource directives are unchanged.

Optional pair-table columns are `target_fasta`, `target_gff`, `query_fasta`,
`query_gff`, `target_feature`, `target_attribute`, `query_feature`,
`query_attribute`, `target_seqids`, and `query_seqids`. Paths are absolute or
workspace-relative. Multiple same-species source files require explicit paths;
explicit FASTA paths also require `synteny_sequence_mode=protein` or `cds`.
GFF feature/attribute overrides must be supplied together. Seqids are exact,
comma-separated chromosome names in the desired ribbon display order.
Unlisted chromosomes remain in the analysis and full-genome dotplot.

FASTA IDs match GFF attributes exactly or after removal of the exact species
prefix. Ambiguous matches fail. The annotation mapper bundled in kfFractBias
resolves known gene loci; the longest available protein/CDS per locus is selected,
with identifier ordering breaking length ties. Unknown locus relationships are
reported in the summary instead of being guessed. All original FASTA IDs,
selected representatives, normalized JCVI IDs, excluded isoforms and unmapped
records remain in the ID mapping tables. The default required mapping fraction
is `1`; inspect source-release and identifier mismatches before changing it.

## Analysis and outputs

Default JCVI filters are `synteny_cscore=0.7`, `synteny_min_anchors=4`, and
`synteny_search_distance=20`. There is no quota filter: many-to-many blocks
remain available, including duplications. The result alone does not establish
whole-genome duplication or ancestral chromosome structure. Tandem-duplicate
filtering follows JCVI's algorithm. Both seed and lifted anchors are retained.

Outputs are under `workspace/output/genome_evolution/synteny/`:

| Directory | Contents |
| --- | --- |
| `analysis/<analysis_id>/` | BED and protein inputs, ID maps, seed/lifted anchors, `anchors.tsv`, `blocks.tsv`, `summary.json`, commands and logs |
| `plots/<analysis_id>/` | `karyotype` and `dotplot` in PDF/SVG/PNG, `display.tsv`, reproducible layout/seqids and copied plotting inputs |

Both plots use **gene rank**, rather than physical chromosome length. BED/block
coordinates use 0-based, half-open intervals. Ribbons are colored by the target
chromosome; dotplot colors distinguish block orientation. Chromosome labels and
their original orientation are preserved; natural ordering puts Chr2 before
Chr10. Display exclusions are recorded, and the dotplot retains every anchor
without downsampling. Short/unplaced scaffolds are included by default, so
fragmented assemblies may require an explicit ribbon seqid list for readability.

The summary records source hashes, selected/unmapped gene counts, syntenic gene
fractions, block/anchor counts, genetic codes, algorithm parameters, tool versions
and annotation-mapper identity. A run with no blocks fails before publication.
Analysis and plots have separate standard artifact-provenance contracts and
recoverable bundle publication. A stage lock prevents concurrent runs from
mixing analysis and plots. Completed analysis is retained if a plot fails.

## Redraw and resume

An unchanged rerun reuses both analysis and plots. A changed display order or
output format invalidates only the plots. Use this after changing `*_seqids`:

```bash
GG_GENOME_EVOLUTION_GENOME_EVOLUTION_MODE=synteny \
GG_GENOME_EVOLUTION_SYNTENY_PLOT_ONLY=1 \
GG_GENOME_EVOLUTION_SYNTENY_PLOT_FORMATS=pdf,svg,png \
artifact_stale_policy=rebuild \
  bash workflow/gg_genome_evolution_entrypoint.sh
```

Plot-only requires a complete, current analysis. Changed FASTA/GFF inputs or
analysis settings cannot silently use the previous analysis. The normal
`artifact_stale_policy=stop` stops on drift; an explicitly selected `rebuild`
regenerates the affected stage. Plot-only does not rebuild analysis. Review
`summary.json`, both figures and `display.tsv` before interpreting real data.

For plot-only, stale analysis is rejected even with `artifact_stale_policy=reuse`.
Input hashes are rechecked before publication to detect changes during a run.
