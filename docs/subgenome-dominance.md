# Subgenome retention and expression contrasts

An optional `gg_genome_evolution` stage consumes completed synteny/phasing and
expression results. Set `genome_evolution_mode="subgenome"` for an independent
run, or `run_subgenome_dominance=1` in `all` mode. It does not assign subgenomes
from gene number or expression, and does not rerun quantification or synteny.

`subgenome_manifest` defaults to `workspace/input/subgenome_analyses.tsv`.
Its TSV columns are `analysis_id`, `species`, `mapping_file`, and optional
`retention_file`, `homoeolog_file`, `expression_file`, `samples_file`,
`mapping_validated` and `retention_validated` (both 0 by default). Input paths are absolute or relative to the
manifest. Analysis IDs must be unique safe path components.
`expression_unit` defaults to `TPM`; `FPKM` is also accepted for within-sample
copy ratios. Other units, including counts and log transforms, are rejected.

Input tables:

| File | Required TSV columns | Contract |
| --- | --- | --- |
| mapping | gene_id, group_id, subgenome, assignment_scope, assignment_basis, evidence | One row per gene. Scope is local or global. Basis is synteny, phylogeny, kmer or curated, supported independently of retention/expression. At least two labels per group. |
| retention | group_id, block_id, locus_id, subgenome, callable, retained | One row per ancestral/reference opportunity per expected subgenome, including explicit uncallable rows. Boolean values are 0/1. All sides of a locus share one nonoverlapping block. |
| homoeolog | pair_id, block_id, gene_id | One independently inferred homoeolog gene per subgenome per set. Genes cannot occur in multiple sets. Tandem copies and ambiguous many-to-many assignments must be resolved upstream. |
| expression | gene_id, sample columns | Raw, finite, nonnegative TPM or FPKM for the mapped representative transcript. Do not use OG sums, log-transformed values, counts without effective-length correction, or inconsistent annotations. |
| samples | column, tissue, biological_id | Columns must match the expression table. Technical replicates sharing a tissue/biological_id are averaged. Tissue-specific biological IDs are required. |

Native kfFractBias `*.genes.tsv` is an intermediate, not a callable-opportunity
table: zero outside the expected collinear interval is not evidence of loss.
Define matched ancestral intervals independently; exclude assembly gaps and
unresolved annotation/alignability. Compare the same reference loci on both
sides, using outgroup and window-size sensitivity analyses. Self-synteny
retention alone is conditional on surviving genes and cannot establish loss.
Set `retention_validated=1` only after those callability and annotation checks.
Otherwise the result is explicitly exploratory syntelog detection, and a zero
must not be described as proven ancestral gene loss.

Each local group is analysed separately. Labels A/B in distinct local groups
must not be pooled or interpreted as parental subgenomes. `global` requires
independent phasing across groups; it adds `ALL_GROUPS` contrasts using pooled
loci with group-specific resampling blocks. Ancestral locus IDs must be unique
across groups for global retention. Local assignments leave genome-wide
dominance untested. More than two labels are
supported by pairwise contrasts. Display ordering of a synteny plot is not a
validated assignment.

For retention the effect is the paired difference in retained fractions.
For expression it is the mean log2(TPM_A/TPM_B), after averaging technical
replicates and then biological-sample log ratios. Both copies must be positive
in every biological sample of a tissue; exclusions are recorded in
`expression_coverage.tsv`. No pseudocount is added. The separate
`expression_detection_difference` metric compares the paired fraction of
biological samples with positive expression. Missing gene IDs remain unassessable; genuine
numeric zeros contribute to detection, not a log ratio. Analyse
mapping ambiguity separately, using unique diagnostic reads or simulation
before setting `mapping_validated=1`.

Statistics resample **nonoverlapping independently defined blocks**, not
overlapping sliding windows. At least three blocks are required for a 95%
cluster bootstrap interval and a two-sided block sign-flip test. With up to 16
blocks the sign-flip test is exact; otherwise it uses Monte Carlo draws. BH
adjustment is within each metric/analysis across groups, contrasts and tissues.
These intervals describe loci conditional on the sampled tissues, annotation
and mapping; they do not measure between-population biological uncertainty.
Block independence/exchangeability is a scientific assumption to check upstream.
Non-significance is not evidence that dominance is absent.

Outputs under `workspace/output/genome_evolution/subgenome_dominance/<analysis_id>`:
`statistics.tsv`, `expression_pairs.tsv`, `expression_coverage.tsv`,
`summary.json`, and PNG/SVG contrasts when estimable. The stage retains input
SHA256s, seed, NumPy version, input/implementation provenance, and uses the
existing transaction lock, stale-policy and atomic publication mechanisms.
Missing retention or expression data remain explicit `not_estimable` results.

`contrasts_absolute.png` and `.svg` also show A/B-invariant magnitudes: the
absolute retention/detection difference or absolute **mean** log2 expression
ratio, rather than a mean of absolute per-gene ratios. Confidence intervals are
images of the signed 95% intervals under the absolute-value transform; an
interval spanning zero has a lower bound of zero. These plots discard direction
and do not show whether retention and expression favour the same copy.
The signed statistics and their P/q values are retained unchanged.

Reference: [Saul et al. (2023), Nature Plants](https://www.nature.com/articles/s41477-023-01562-2).
Its retention and expression evidence should be assessed separately; a gene
copy with higher expression is not automatically the dominant subgenome.
