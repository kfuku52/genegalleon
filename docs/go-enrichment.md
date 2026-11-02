# CAFE GO enrichment

The native workflow accepts `go_annotation_file` as an absolute path to an
existing annotation table when annotations were produced in another workspace.
Set the ordinary `annotation_species` to the exact orthogroup species column;
the workflow forwards it explicitly. Direct R calls accept an optional eleventh
`annotation_species` argument. Without that argument, both native
`Species_annotation.tsv` and legacy `Species.annotation.tsv` names are accepted.
An unmatched reference species stops instead of silently losing annotations.
For a table outside the workflow workspace, set
`GG_CONTAINER_PROJECT_ROOT_BIND=/path/to/project:/path/to/project` to the
project directory containing both the workspace and the annotation table.
SHIROKANE, NIG and audrey1 honor this explicit same-path project mount;
unrelated host directories are not added automatically.

Native UniProt GO cross-references retain one aspect, term and evidence entry
per GO ID. The GO reader accepts aligned C/F/P or CC/MF/BP aspects and native or
legacy semicolon spacing. Older metadata that independently deduplicated aspects
has lost those associations; rebuild it from the matching cached DAT, rather than
guessing the missing associations. Completed best hits can be reused:

```bash
python workflow/support/build_uniprot_metadata_from_dat.py \
  --uniprot_dat frozen/uniprot_sprot.dat.gz --outfile fixed.meta.tsv.gz
python workflow/support/reformat_uniprot_diamond.py \
  --existing_annotation_tsv completed/Arabidopsis_thaliana_annotation.tsv \
  --uniprot_meta_tsv fixed.meta.tsv.gz --outdir refreshed --audit_json refresh.json
```

Repeat `--existing_annotation_tsv` for a batch. Only the four GO columns are
refreshed; gene order, accepted accessions, search statistics, orthogroups and
all other annotations are preserved. No homology search is repeated.

When species sampling changes, map a previous CAFE target by descendant taxa:

```bash
python workflow/support/cafe_branch_map.py \
  --old-tree old/Gamma_asr.tre --new-tree new/Gamma_asr.tre \
  --old-target '<9>' --new-branch-table new/Gamma_change.tab \
  --output branch_mapping.json
```

This reads the first native tree, preserves full species keys, and reports added
or missing taxa and all matching non-root clades. A unique match may include
explicitly reported added species. A split edge, missing target species, or
non-monophyletic target gives exit 2 and no selected branch. Resolve biological
intent before analysis in those cases; old numeric branch IDs alone are not
stable identifiers between runs. The mapping does not infer orthology or
ancestral gene content.

`go_enrichment_method="event"` remains the default. Its event counts, numerical
results for the same inputs, GO candidate set and output paths are unchanged.
The optional `cafe_branch_flags` method consumes standard, unmodified CAFE
outputs. It does not fit new models, simulate families, run an additional CAFE
process, or require a custom CAFE executable.

For new analyses, use `go_enrichment_method="cafe_branch_flags_all_go"`.
It uses the same native family screen, but fixes the GO hypotheses from **all
annotated background families before examining target events**. Zero-hit terms
remain in the correction, each family counts once, and `both` uses one BH
correction across gains and losses together. This removes target-observed GO
filtering and repeated-family counting; native screening and annotation
dependencies still make the results exploratory rather than calibrated causal
or rate-contrast evidence. Legacy `event` now emits an explicit warning, while
retaining its numerical results and paths.

The fixed-universe outputs use `go_enrichment/cafe_branch_flags_all_go/` with
the same table schemas as `cafe_branch_flags`. `branch_flags_metadata.tsv`
records `all_annotated_background_go` and `BH_all_requested_directions`.
Switching methods selects a separate output directory and changes the artifact
contract; use the normal stale-artifact review/rebuild procedure rather than
renaming old results into the new directory. The sections below describe the
original target-restricted `cafe_branch_flags` mode unless stated otherwise.

## Fixed event GO universe and Tarone sensitivity

Use `event_all_go` to compare multiplicity corrections while retaining the
legacy **event** population and one-sided Fisher tests. This is separate from
the target-exclusive family screen in `cafe_branch_flags_all_go`. Freeze the GO
catalog once, before examining a target branch or direction:

```bash
Rscript workflow/support/cafe_go_enrichment.r --freeze-event-go-universe \
  saved/Gamma_change.tab saved/Orthogroups.selected.tsv \
  annotations/Arabidopsis_thaliana_annotation.tsv frozen_event_go.tsv \
  BP,MF,CC Arabidopsis_thaliana
```

The catalog contains every GO annotated on a family present in the CAFE change
table, restricted to the requested aspects and fixed reference annotation.
Only family identities are used from the change table; branch changes,
probabilities, target and direction do not determine the catalog. This is not
the complete ontology, and no new GO ancestor propagation is added. Existing
catalogs cannot be overwritten. Preserve the catalog, source-input hashes and
its own hash as part of the analysis contract.

```bash
Rscript workflow/support/cafe_go_enrichment.r \
  saved/Gamma_change.tab saved/Gamma_branch_probabilities.tab \
  saved/Orthogroups.selected.tsv annotations/Arabidopsis_thaliana_annotation.tsv \
  results/branch_37 '<37>' decrease BP,MF,CC event_all_go .05 \
  Arabidopsis_thaliana /absolute/path/frozen_event_go.tsv
```

For the workflow entrypoint, select `go_enrichment_method="event_all_go"` and
set `go_fixed_universe_file` to the absolute catalog path. Use the same catalog
for all target branches. Each invocation verifies the complete catalog against
its eligible families and annotation; mismatches fail rather than silently
redefining the background. The method accepts one direction per invocation.

All frozen GO rows remain in `_all_go.tsv`, including zero-target and zero-event
terms with one-sided P=1. `p_value_adjusted` is BH across the **full** catalog.
Existing event-only raw P values and event counts are preserved, but BH values
can increase because the hypothesis family has changed. `_significant_go.tsv`
contains full-catalog BH discoveries. `_tarone_significant_go.tsv` contains the
separate Tarone discoveries at alpha=0.05, and `_all_go.tsv` additionally records
minimum attainable P, its log, testability, Tarone decisions and full-catalog
Bonferroni results.

For each GO, the smallest attainable one-sided Fisher P is the hypergeometric
mass at the maximum feasible overlap, conditional on its two-by-two margins.
The implemented Tarone (1990) variant chooses the smallest integer K such that
`M(K) = count(p_min < alpha/K) <= K`, then rejects only tests with
`P < alpha/K`. Strict inequalities match the existing P<alpha convention.
`event_all_go_metadata.tsv` records the total GO count, K, the distinct testable
count M(K), threshold and `M(K)*alpha/K` union bound. K and M(K) are not
interchangeable. No Tarone q value is fabricated: the decision is specific to
the recorded alpha and variant.

Tarone targets **FWER**, whereas BH targets **FDR**; a reduced testable count
does not guarantee more discoveries than BH. Bounds require valid conditional
tests. Freezing GO does not repair repeated-family dependence in the legacy
event model or turn its screen into a branch-rate comparison. Treat these
outputs as a sensitivity analysis, preserve the original BH results, and compare
original BH, fixed-catalog BH and Tarone separately. No CAFE refit is needed.

Reference: R. E. Tarone (1990), *A modified Bonferroni method for discrete data*,
[Biometrics 46:515–522](https://pubmed.ncbi.nlm.nih.gov/2364136/).
The integer-K formulation is also described in
[Gilbert (2005), section 3.2](https://faculty.washington.edu/peterg/Vaccine2006/articles/GilbertPaperAppliedStatistics.2004.pdf).

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
changes and probability ranges are checked for consistency. The branch-probability
path must match the change-file prefix; mixing reports from separate runs or
models is rejected even when their family IDs and columns happen to agree.

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
results. The runtime tests generate standard Base and Gamma CAFE outputs using ordinary
fixed-parameter options, consume those files, and check gain/loss screening,
rejection of other-branch flags, preservation of unselected GO candidates,
reprocessing, source-file integrity and rejection of mismatched report paths.
Full CLI regressions also cover empty branch reports and internal targets. This tests execution, not statistical
power or end-to-end calibration. Docker validation does not establish SIF
compatibility.

The native output/reporting semantics are described in the
[CAFE5 documentation](https://github.com/hahnlab/CAFE5) and implemented in its
[standard report generation](https://github.com/hahnlab/CAFE5/blob/master/src/execute.cpp)
and [branch-probability calculation](https://github.com/hahnlab/CAFE5/blob/master/src/gene_family_reconstructor.cpp).
