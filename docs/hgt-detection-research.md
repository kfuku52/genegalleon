# HGT Detection Research Notes for GeneGalleon

Date: 2026-03-17

## Executive Summary

The strongest overall strategy for horizontal gene transfer (HGT) detection in GeneGalleon is:

1. Use species-tree-aware gene-tree reconciliation as the primary decision engine.
2. Use fast sequence-similarity filters only for candidate generation and prioritization.
3. Add aggressive eukaryote-focused validation layers to suppress false positives.

For GeneGalleon specifically, the best immediate upgrade path is to build around the existing GeneRax integration and switch the reconciliation model from duplication-loss only to duplication-transfer-loss for HGT-focused runs.

## What GeneGalleon Already Has

GeneGalleon already contains most of the structural pieces needed for a serious HGT workflow:

- per-family alignments and IQ-TREE inference
- rooted species trees and orthogroup mode
- GeneRax integration
- branch/table summarization that already imports `generax_event`
- GFF-derived intron and positional annotations
- expression support
- synteny neighborhood summaries
- upstream contamination filtering in the annotation/transcriptome stages

Current bottlenecks in the shipped defaults:

- `run_generax=0` by default
- `generax_rec_model="UndatedDL"` by default, so transfers are not modeled in default runs
- tree visualization defaults to `species_overlap`, not GeneRax events
- NOTUNG calls use `--infertransfers false`

This means GeneGalleon already has the plumbing for transfer-aware analysis, but HGT is not yet the primary operating mode.

## Method Comparison

| Method family | Accuracy / practical performance | Main strengths | Main failure modes | Fit for GeneGalleon |
| --- | --- | --- | --- | --- |
| Gene-tree/species-tree reconciliation with DTL | Best overall when alignments and species tree are reasonably good; strongest direct HGT signal | Explicitly models duplication, transfer, and loss | Sensitive to gene-tree error, bad rooting, sparse taxon sampling, and biological confounders such as ILS/hybridization | Excellent |
| Cost-based reconciliation with uncertainty handling | Usually weaker as a final decision engine than full ML, but very useful for robustness checks | Can account for rooting/topology uncertainty and alternative optimal reconciliations | Still inherits gene-tree uncertainty; event costs can influence calls | Good as optional confirmation |
| Best-hit / taxonomic-distribution screening | Good speed and often good recall for first-pass screening | Cheap, scalable, robust candidate generation | Weak as standalone evidence; sensitive to sampling gaps and donor/recipient database bias | Excellent as prefilter |
| Hybrid best-hit + phylogeny pipelines | Better than either component alone in large-scale screening contexts | Practical compromise between scale and precision | Still depends on downstream phylogenetic validation quality | Very good |
| Composition-only methods | Useful mainly for recent prokaryotic transfers | Fast and reference-light | Weak for ameliorated transfers and especially weak in eukaryotes | Poor as primary method |
| Community/metagenome HGT tools | Useful for MAG/metagenome settings | Works without isolate-quality references in some cases | Assembly fragmentation limits recent-HGT recovery; not designed for clean phylogenomic family analysis | Low for current GeneGalleon use cases |

## What the Primary Literature Says

### 1. DTL reconciliation is the highest-value core method

The clearest evidence comes from the GeneRax paper. GeneRax is a species-tree-aware maximum-likelihood method that jointly uses sequence evolution and duplication-transfer-loss reconciliation. In the authors' simulations it recovered trees closest to the true tree in most runs, and it was also fast at scale when starting from alignments rather than precomputed bootstrap trees.

Implication for GeneGalleon:

- HGT detection should be centered on GeneRax-style DTL reconciliation, not on BLAST score heuristics alone.
- Orthogroup mode is the natural place to run it.

Relevant source:

- GeneRax paper: <https://pubmed.ncbi.nlm.nih.gov/32502238/>
- GeneRax project: <https://github.com/BenoitMorel/GeneRax>

### 2. Uncertainty-aware reconciliation is important for support, not just speed

RANGER-DTL 2.0 is explicitly designed to handle uncertain rooting, topological uncertainty, multiple optimal reconciliations, and alternative event costs. This is valuable because raw HGT calls are often unstable when the gene tree is weakly resolved.

Implication for GeneGalleon:

- GeneRax should be the main engine.
- A second-pass uncertainty-aware reconciliation layer is attractive for confidence scoring on high-value candidates.

Relevant source:

- RANGER-DTL 2.0: <https://pubmed.ncbi.nlm.nih.gov/29688310/>

### 3. Best-hit methods are useful, but mainly as candidate generators

HGTector was built specifically to improve first-pass HGT screening over naive BLAST heuristics. Its key advantage is that it uses hit-distribution statistics across predefined phylogenetic groups, which makes it less fragile to gene loss, rate variation, and database noise than simple top-hit logic.

Implication for GeneGalleon:

- A HGTector-like score is worth adding as a fast prefilter.
- It should not be treated as final evidence.

Relevant source:

- HGTector: <https://pubmed.ncbi.nlm.nih.gov/25159222/>

### 4. Hybrid pipelines outperform single-signal pipelines in practice

MetaCHIP combines best-match and phylogenetic approaches. Although it targets community-level microbial data rather than clean eukaryotic orthogroups, the architectural lesson is important: fast similarity screens plus phylogenetic confirmation scale well and recover transfers across different divergence levels better than a single signal alone.

Implication for GeneGalleon:

- Use a two-stage workflow: cheap genome-wide screening, then reconciliation on prioritized families.

Relevant source:

- MetaCHIP: <https://pubmed.ncbi.nlm.nih.gov/30832740/>

### 5. Composition-only approaches are a poor default for eukaryotic datasets

The fungal data-quality study is particularly important for GeneGalleon because many intended use cases are eukaryotic. That study found that both composition-based and phylogeny-based automated detection had limited statistical power in complex fungal datasets, with extreme gene-tree variation and heavy dependence on data quality.

Implication for GeneGalleon:

- Composition-based HGT calls should never be used as final evidence in eukaryotes.
- Even phylogenetic calls need stronger validation layers than in prokaryotes.

Relevant source:

- Genomic Data Quality Impacts Automated Detection of Lateral Gene Transfer in Fungi: <https://pubmed.ncbi.nlm.nih.gov/28235827/>

### 6. False positives from contamination and assembly problems are a major eukaryotic risk

The tardigrade re-analysis is the canonical warning: an apparent extreme HGT signal largely disappeared after contamination-aware reassessment.

Implication for GeneGalleon:

- contamination control is not a side issue; it is part of HGT detection itself
- HGT candidates from short scaffolds, low-support loci, or suspicious neighborhood context need extra skepticism

Relevant source:

- Tardigrade reassessment: <https://pubmed.ncbi.nlm.nih.gov/27173901/>

### 7. Large-scale eukaryotic HGT studies rely on phylogenomics plus contextual filtering

The ochrophyte phylogenomic study shows what successful large-scale eukaryotic HGT work looks like: dense phylogenomic screening, functional interpretation, and a focus on the evolutionary pattern across clades rather than isolated top-hit anomalies.

Implication for GeneGalleon:

- the right target is a family-aware phylogenomic HGT workflow, not a per-gene anomaly detector
- orthogroup-level summaries and branch-level aggregation are therefore the correct product shape

Relevant source:

- Ochrophyte phylogenomics: <https://pubmed.ncbi.nlm.nih.gov/33419955/>

## Recommended Detection Architecture for GeneGalleon

### Tier 1: Fast candidate generation

Add a cheap screen before full reconciliation:

- donor-vs-ingroup best-hit contrast using existing BLAST/MMseqs outputs
- taxonomic hit-distribution score in the style of HGTector
- optional compositional outlier score, but only as a weak auxiliary feature

This stage should optimize recall and runtime, not final precision.

### Tier 2: Primary confirmation with GeneRax DTL

For HGT-focused runs:

- set `run_generax=1`
- set `generax_rec_model="UndatedDTL"`
- use orthogroup mode by default
- expose GeneRax H events directly in candidate tables
- make `treevis_event_method="generax"` the default for HGT report mode

This is the single most important upgrade.

### Tier 3: Eukaryote-focused false-positive suppression

A candidate should be down-ranked or rejected if several of the following are suspicious:

- no expression support
- no introns in an otherwise intron-rich lineage
- weak or broken local synteny
- locus sits on a contamination-prone or low-context scaffold
- donor support disappears after broader taxon sampling
- gene-tree support is weak or unstable

GeneGalleon already has most of these signals except explicit contig-level QC for final HGT reports.

### Tier 4: Confidence support from uncertainty-aware reconciliation

For top candidates, optionally run an uncertainty-aware confirmation layer:

- RANGER-DTL 2.0 across alternative rootings / bootstrap trees
- consensus support for whether transfer remains necessary

This is slower, so it should be reserved for a narrowed candidate set.

## Recommended Scoring Model

Instead of a binary call, produce an evidence tier:

- `Tier A`: strong DTL support plus clean genomic context plus expression/synteny support
- `Tier B`: DTL support but one contextual weakness
- `Tier C`: best-hit and taxonomic-distribution anomaly only
- `Rejected`: contamination-like or low-support candidates

Suggested score components:

- GeneRax transfer event present
- branch support on the relevant gene-tree split
- donor-vs-recipient hit contrast
- taxonomic patchiness
- expression detected
- intron pattern plausible for host integration
- synteny retained
- contamination warning flags absent

## Practical Recommendations by Priority

### Highest priority

1. Add an explicit HGT mode that flips GeneRax to `UndatedDTL`.
2. Emit transfer-aware summary tables from existing `generax_event` fields.
3. Add a fast HGTector-like prefilter using current search outputs.

### Medium priority

1. Add composite confidence scoring.
2. Add stricter contamination and scaffold-context reporting for candidate loci.
3. Add candidate ranking reports at the orthogroup and gene levels.

### Later priority

1. Add optional RANGER-DTL-style uncertainty confirmation.
2. Evaluate whether dated/time-consistent reconciliation should be added using the existing dated species tree outputs.

Inference:

- Because GeneGalleon already builds dated species trees, a future time-consistent transfer filter is likely to improve precision relative to undated DTL alone, especially for deep eukaryotic datasets. This is an inference from the current pipeline structure plus the known limitations of undated reconciliation, not a direct benchmark result from the sources listed above.

## Methods I Would Not Make Primary

- composition-only HGT detection
- naive top-hit or Alien Index only
- NOTUNG as the main HGT caller in the current GeneGalleon configuration

These can still be useful as auxiliary signals or for niche datasets, but they should not define the final decision rule.

## Bottom-Line Recommendation

If the goal is the best practical HGT performance in GeneGalleon, the right design is:

- HGTector-like candidate generation
- GeneRax `UndatedDTL` as the main confirmation engine
- eukaryote-specific QC using contamination, intron, expression, and synteny signals
- optional uncertainty-aware confirmation for high-confidence calls

That combination is the best match to both the published literature and GeneGalleon's existing architecture.
