# BUSCO quality diagnostics and trait correlations

Enable `run_copy_number_quality_diagnostics=1` in the genome-evolution
entrypoint, or use the environment override:

```bash
GG_GENOME_EVOLUTION_RUN_COPY_NUMBER_QUALITY_DIAGNOSTICS=1 \
bash workflow/gg_genome_evolution_entrypoint.sh
```

This optional stage prepares the shared copy-number matrix even when CAFE and
the ordinary trait-PGLS stage are disabled. It uses the dated species tree,
the existing shared species BUSCO **short summaries**, and the numeric columns
of `species_trait.tsv`. No BUSCO job is started by this setting alone. Existing
sequence/tree producer flags retain their usual behavior.

The output directory is
`workspace/output/genome_evolution/orthogroup_copy_number/quality_diagnostics/`.
`manifest.tsv` lists the files belonging to the current run; older or unrelated
files in the directory are not implicitly current. Failed fits have explicit
statuses, and a failed publication preserves the previously published bundle.

## Inputs and scope

- `copy_number_quality_busco_table="auto"` reads shared short summaries. An
  explicit path is workspace-relative (or absolute) and must contain a first
  species column and `busco_complete_pct`, in **percent units from 0 to 100**.
  Optional `lineage`, `mode`, `busco_version`, and `source` columns document
  comparability. The stage records missing metadata rather than claiming it has
  verified comparability. Known mixed lineages, modes, or versions among the
  analyzed species are rejected. Use comparable runs before interpreting their
  quality gradient; a common lineage still needs to be appropriate for the taxa.
- Completeness is C = S + D. Duplicated BUSCOs are complete; C is not the
  single-copy percentage. Missing BUSCO measurements remain missing, never zero.
  These percentages do not estimate family-specific detection probabilities.
- Species join by complete labels with spaces normalized to underscores;
  duplicates and ambiguous short-summary filenames fail. Trait and quality
  rows may be in any order. Tree species define the analyzed population.
- `orthogroup_copy_number_trait` selects the biological trait columns; BUSCO
  completeness is always added to the plot and tested separately. The same
  hash-bound `.schema.json` and `.metadata.json` rules as ordinary PGLS apply:
  declared categorical/text columns are excluded and reported in
  `excluded_traits.tsv`; undeclared invalid numeric values fail validation.
  Numeric 0/1 traits may be plotted. GBIF observation columns require explicit
  selection and valid acquisition metadata; ineligible observations are masked.
  Preserve the generated input audit when interpreting or sharing the plots.
  A conflicting `busco_complete_pct` already in the trait table fails rather
  than silently replacing either source.
- `orthogroup_copy_number_trait_family_ids`,
  `orthogroup_copy_number_trait_family_file`, and
  `orthogroup_copy_number_trait_max_families`
  set the family scope. Set it before inspecting biological trait results.
  Unknown/non-finite, negative, or fractional copy counts are rejected.
- The inherited minimum species setting must be at least four. This is a
  computational minimum, not a guarantee of calibrated small-sample inference.

Example external table:

```tsv
species	busco_complete_pct	lineage	mode	busco_version
Arabidopsis_thaliana	98.7	embryophyta_odb12	proteins	6.0.0
Amborella_trichopoda	96.2	embryophyta_odb12	proteins	6.0.0
```

The shared protein/CDS gene-set assessment is used for gene-catalog diagnostics.
Genome-mode completeness need not measure completeness of the annotations
actually counted; transcriptome completeness also reflects expression and
sampling. Quality flags cannot distinguish these mechanisms.

## Outputs and interpretation

| File | Meaning |
| --- | --- |
| `busco_quality.tsv` | Tree-matched percentages, sources, metadata and missing/comparability statuses |
| `traits_with_busco.tsv` | Numeric trait table with BUSCO completeness joined by species |
| `trait_correlations.pdf`, `.svg` | Descriptive correlation matrix, including pairwise species counts |
| `trait_correlations.tsv` | Every ordered trait pair, correlation, pairwise n, method and status |
| `family_busco_associations.tsv` | Family-wise Brownian PGLS: `BUSCO completeness ~ log1p(copy number)`; raw P, BH q and quality flag |
| `trait_quality_sensitivity.tsv` | Biological trait PGLS comparisons, with the family quality flag joined |
| `sensitivity_cohorts.tsv` | Species membership for each biological trait's sensitivity cohorts |
| `cafe_family_quality.tsv` | Original CAFE family results with extant-count quality flags attached, when CAFE runs alongside diagnostics |
| `pgls_quality.tsv` | Original trait-PGLS rows and P/q values with family quality flags attached, when trait PGLS runs alongside diagnostics |
| `selected_species_traits.tsv` and sidecars, `species_trait_input.json`, `trait_selection.tsv`, `excluded_traits.tsv` | Selected inputs, type/provenance audit and exclusions when a trait table is supplied |
| `parameters.tsv`, `manifest.tsv` | Settings and authoritative current output members |

The matrix defaults to Spearman correlation. Use
`copy_number_quality_correlation_method="pearson"` for Pearson correlation.
Correlations use pairwise complete data, require at least three species and
variation in both variables, and are **not phylogenetically corrected tests**.
They have no significance stars. Missing/unavailable cells are gray and labeled
NA, and constant traits are not assigned an artificial diagonal correlation of 1.

The family diagnostic uses NWKIT Brownian REML GLS with an intercept and tests
the copy-number coefficient. Completeness is a bounded response; this Gaussian
diagnostic is exploratory and its calibration is not guaranteed, especially
with small samples or scores near 100%. BH covers every prespecified family,
including the planned multiplicity of unavailable fits (whose own q remains NA).
`orthogroup_copy_number_trait_alpha`, default 0.05, controls the diagnostic flag:

- `quality_associated`: BUSCO association passes the diagnostic BH cutoff.
- `no_quality_association_detected`: a usable test did not pass; not evidence
  that the family is free of annotation/detection error.
- `not_assessable`: missing quality data, insufficient species/variation, or
  unavailable inference. CAFE families outside the diagnostic scope instead
  receive `not_assessed`.

No family, PGLS association, CAFE branch, or GO input is removed based on these
flags. BUSCO loci overlap the counts being analyzed and completeness can vary
biologically: an association is a **quality warning, not a false-positive verdict**.
The attached CAFE diagnostic concerns extant family counts. It neither retests
branches nor corrects CAFE P-values, and CAFE's original columns remain intact.
This stage does not estimate a CAFE error model or rerun CAFE on a pruned tree.

## PGLS sensitivity comparisons

With `copy_number_quality_sensitivity=1` (default within the enabled stage),
the same family/trait set receives four separately labeled phylogenetic fits:

1. `baseline_all`: `trait ~ log1p(copy number)` on all available trait species.
2. `baseline_busco_observed`: the same model on species with BUSCO measurements.
3. `busco_adjusted`: `trait ~ log1p(copy number) + BUSCO completeness` on exactly
   the same cohort as (2).
4. `high_completeness`: the baseline model restricted to C at or above
   `copy_number_quality_high_completeness` (default 95%). Choose this threshold
   before inspecting associations; 95% is not a universal quality guarantee.

Compare (2) against (3) to separate adjustment from missing-species effects.
Each variant has its own BH column `quality_qvalue`, covering all planned
family × biological trait tests for that variant. This does not control error
rates for picking whichever variant is most significant. Coefficients and
standard errors remain available for comparison; a P-value crossing 0.05 is
not an automatic interpretation rule. The original root/shared Brownian history
is retained under subsetting. Collinear adjustment or exhausted residual degrees
of freedom are reported as not estimable.

Sensitivity fits inherit `orthogroup_copy_number_trait_response_families`
(default Gaussian). Explicit `trait=binomial`, `trait=poisson`, and
`trait=negative-binomial` mappings use the same NWKIT response families and
links as ordinary trait PGLS. A numeric 0/1 column alone does not select logistic
regression. Standalone callers use `--response_families=trait=binomial` (or the
other mappings). The shared adapter validates raw counts and applies `log1p`
exactly once. The ordinary PGLS outputs and their global BH family
remain unchanged. No measurement-error model is inferred from BUSCO, and the
covariate adjustment is not a causal or observation-error correction.

## Standalone plotting or diagnostics

For plots alone, run inside a GeneGalleon container with no copy-number argument:

```bash
Rscript workflow/support/copy_number_quality_diagnostics.r \
  --file_sptree=/data/dated_species_tree.nwk \
  --file_trait=/data/species_trait.tsv \
  --file_busco=/data/busco_quality.tsv \
  --outdir=/data/trait_correlations
```

Add `--file_copy_number=/data/orthogroup_copy_number.tsv` for the family flags
and sensitivity fits. Alternatively supply
`--busco_short_dir=/data/species_busco_short` instead of `--file_busco`.
An explicit `--file_cafe_results=/data/Gamma_family_results.txt` joins an existing
CAFE family table; the caller must ensure it belongs to the same count analysis.
Likewise, `--file_pgls_results=/data/orthogroup_copy_number_trait_pgls.tsv` annotates
an existing PGLS table without changing its tests or adjusted P-values.
Use `--sensitivity=0` for family diagnostics and plotting only.

Validation uses `workflow/tests/test_busco_quality_metadata.py` and
`workflow/tests/test_copy_number_quality_diagnostics.R` in a GeneGalleon runtime.
Tests cover parsing/joining, missing/invariant data, protocol conflicts, independent
matrix-GLS coefficient/SE calculations, identical adjustment cohorts, preserved
subtree covariance, CAFE annotations without filtering, and publication failure.
They do not establish biological specificity or empirical FDR calibration.

With F families and T numeric biological traits, sensitivity analysis can request
up to F × (1 + 4T) NWKIT fits; unavailable/invariant cohorts are skipped. Use a
prespecified trait/family scope for a pilot, or set sensitivity to 0 when only
the quality flags and matrix are needed. Wall time depends on species count and
the fitting runtime; no performance benchmark is claimed here.

### Validation record (2026-09-10)

Validation used a clean snapshot of main plus this change, excluding unrelated
working changes. All diagnostic tests and the full declared R validation passed.

- Broad Python run: **2,506 passed**. Four calibration tests failed because the
  initial runtime lacked the R package `posterior`; these were environment
  failures outside the diagnostic implementation.
- The aligned Docker runtime, `local/genegalleon:busco-quality-review-20260910`,
  combines the NWKIT installation from `native-ou-dev` with the R environment
  from `gbif-main-review-dev`. **74 targeted Python tests passed**, including all
  four calibration failures, BUSCO parsing and real core/CLI execution, mixed
  response families, trait schemas, GBIF provenance/masking, HGT integration,
  and owned runtime contracts.
- All **16 commands** in the declared R runner passed in the aligned runtime,
  including the treevis package check, original copy-number PGLS, independent
  adjusted-GLS coefficient/SE checks, and every treevis test. The runtime update
  also replaced an older installed treevis that rejected trans-spliced inputs.
- Full Python Ruff and tracked-shell ShellCheck checks passed on the host;
  affected shell syntax and whitespace checks passed. Runtime behavior was
  checked in Docker. The correlation PDF layout was visually checked during
  implementation; PDF fonts are embedded.

SIF execution and large empirical calibration simulations were not performed.
Docker validation does not establish SIF compatibility, biological specificity,
or empirical false-discovery calibration.
