# CAFE GO enrichment

`go_enrichment_method="event"` remains the default. Its event counts, numerical
results for the same inputs, GO candidate set and output paths are unchanged.
The optional `cafe_branch_flags` method consumes standard, unmodified CAFE
outputs. It does not fit new models, simulate families, run an additional CAFE
process, or require a custom CAFE executable.

## Exploratory target-restricted branch flags

```bash
run_go_enrichment=1
target_branch_go="A<1>"              # Exact non-root CAFE column label
change_direction_go="both"          # increase | decrease | both
go_enrichment_method="cafe_branch_flags"
go_family_alpha=0.05
```

A family is selected when all of the following hold:

1. Its **existing CAFE family-wide P value**, adjusted by BH across all
   GO-annotated families in the requested GO categories, is below
   `go_family_alpha`.
2. The target branch has a nonzero reconstructed change and native branch
   probability below 0.05, matching the legacy change-flag convention.
3. No other non-root branch has that flag with a nonzero change, of **either
   sign**. For example, a target gain with a flagged loss elsewhere is excluded.
4. The target change has the requested sign. `both` reports gains and losses
   separately, using the same family BH correction.

This is an operational screen for **flags observed only on the target branch**.
It is not a target-versus-background rate test. A significant target and a
non-significant other branch do not establish a significant difference between
those branches. Unflagged branches may still have changed; their detection power
can differ with duration, family size and reconstruction uncertainty. The method
cannot claim that a selected family evolved exclusively on the target branch.

CAFE branch probabilities are used only as native flags. They are not relabelled
as rate-contrast P values and are not combined into an invented specificity P
value. The reported family P values remain the native family-wide test values.
CAFE's finite simulation resolution, fitted-model uncertainty and possible
reported zero P values remain limitations; a reported zero is not an exact
probability of zero.

## Native inputs and missing reports

The existing change, branch-probability, orthogroup and annotation files are
used, plus matching native `*_family_results.txt` and `*_asr.tre` files. Both
`Gamma_*` and `Base_*` prefixes are supported. The ASR tree identifies the root
exactly so it is not counted as another branch. Family IDs, columns, integer
changes and probability ranges are checked for consistency.

CAFE normally writes branch-probability rows only for families passing its
family-wide reporting threshold. A missing row is recorded as **not reported**,
not as zero flagged branches. Families failing the family BH screen may lack
that row and remain in the GO background. If a BH-selected family has no branch
report, or a reported row has missing non-root probabilities, the analysis stops.
It does not silently remove that family or treat missing evidence as absence.
Provide complete standard CAFE reports with an appropriate native reporting
threshold before attempting that analysis.

## GO hypotheses and background

The GO candidate IDs for each direction are exactly those the legacy event
analysis would test for the same inputs and categories. This set is formed
before the new family screen. A retained GO term with no selected families gets
P=1. No extra GO terms are introduced by this option.

Each annotated family counts once. Selected families are compared with the
remaining annotated families using one-sided Fisher enrichment, followed by BH
within each direction/category run. Duplicate gene/GO annotations do not multiply
family counts. `both` uses separate GO corrections for gains and losses; it does
not provide joint two-direction FDR control.

Because the requested GO candidate set is selected from the observed data and
family/GO hypotheses are dependent, the resulting tables are **exploratory**.
This procedure does not establish end-to-end GO FDR control or resolve the
legacy GO-ascertainment concern. Confirmatory interpretation would require a
separately justified analysis and calibration.

## Outputs and direct use

Core-workflow outputs are separated under `go_enrichment/cafe_branch_flags/`:

- `family_branch_flags.tsv`: native family P and BH-adjusted P, target change and
  branch probability, report availability, target flag, other flagged branches,
  direction, selection and reason for exclusion.
- `branch_flags_metadata.tsv`: the screen, threshold, family/GO scope and explicit
  exploratory interpretation.
- `enrichment_significant_*_all_go.tsv` and `*_significant_go.tsv`: unique-family
  GO tables, with `n_selected_*` and `n_background_*` counts. A no-discovery result
  has a header-only significant table.
- Existing `orthogroup_*significant*` tables keep their original CAFE event
  meaning; use `family_branch_flags.tsv` for the new selected-family list.

The core's artifact provenance includes the native reports, annotations, adapter,
settings and primary results. Published optional summaries are invalidated before
reading new inputs, so a failed rerun cannot leave an old summary appearing to be
current. Input CAFE files are never modified. Reprocessing the same inputs is
cheap and deterministic; no optimizer or simulation cache is needed.

```bash
Rscript workflow/support/cafe_go_enrichment.r \
  Gamma_change.tab Gamma_branch_probabilities.tab gene_ids.tsv \
  reference.annotation.tsv outdir 'A<1>' both BP,MF,CC \
  cafe_branch_flags 0.05
```

The original eight-argument invocation retains `event` mode. The prior local
`cafe_lrt` implementation and its bootstrap/restart settings have been removed;
that name is rejected rather than silently assigned a different meaning. Old
LRT output directories are not reused. Neither a CAFE source patch nor upstream
publication of one is required for this replacement.

## Validation

Run in a normal GeneGalleon container containing unmodified CAFE:

```bash
Rscript workflow/tests/test_cafe_go_enrichment.R
python -m pytest -q --gg-strict-runtime workflow/tests/test_cafe_branch_flags_runtime.py
```

The R regressions check family BH, both signs, exclusion of flags on other
branches (including the opposite sign), normal missing reports, malformed input,
unique-family GO counts, empty selections and unchanged legacy numerical/CLI
results. The runtime test generates standard CAFE output using its ordinary
fixed-lambda option, consumes those files, and checks gain/loss screening,
reprocessing and source-file integrity. This tests execution, not statistical
power or end-to-end calibration. Docker validation does not establish SIF
compatibility.

The native output/reporting semantics are described in the
[CAFE5 documentation](https://github.com/hahnlab/CAFE5) and implemented in its
[standard report generation](https://github.com/hahnlab/CAFE5/blob/master/src/execute.cpp)
and [branch-probability calculation](https://github.com/hahnlab/CAFE5/blob/master/src/gene_family_reconstructor.cpp).
