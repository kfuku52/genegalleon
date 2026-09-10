# CAFE GO enrichment

`go_enrichment_method="event"` remains the default. It compares CAFE-significant
FamilyID × branch events in the target versus other branches. Its target-observed
GO selection, numerical results and output paths are unchanged. The optional
native model comparison below uses CAFE outputs; no manual branch-length table
is required.

## Native CAFE family comparison

This option requires a CAFE build with corrected Nelder–Mead contraction logic.
The stock CAFE 5.1 executable failed the native integration check: independent
fits of a bootstrap family did not agree. Its inside-contraction expression and
second-worst/reflection branching were incorrect; the same defects were present
in the upstream source inspected during implementation. The correction belongs
in CAFE, and is not embedded or worked around in GeneGalleon. Until that change
is available in your runtime, the optional mode may stop on an unresolved fit.
The default event mode does not use these new per-family optimizations.

```bash
run_go_enrichment=1
target_branch_go="A<1>"         # Exact non-root column in Gamma_change.tab
change_direction_go="both"     # increase | decrease | both (both requires cafe_lrt)
go_enrichment_method="cafe_lrt"
go_family_alpha=0.05
go_cafe_bootstrap_replicates=999
go_cafe_fit_restarts=5
go_cafe_max_iterations=1000
```

For each GO-annotated family in the CAFE change table, GeneGalleon asks **CAFE**
to fit two nested birth-death models to the extant family sizes:

| Model | Free change-rate parameters |
|---|---|
| Null | One family-specific lambda for the whole tree |
| Alternative | A background lambda and a separate lambda on the nominated branch |

The null lets a family change rapidly everywhere. The alternative must explain
additional target-branch behavior; a family is not selected merely because it
changes on many branches. Counts are observations at extant species, not
independent reconstructed FamilyID × branch events. Native CAFE performs the
ancestral-state likelihood calculation on the full tree.

The LRT statistic is `max(0, 2 * (NLL_null - NLL_alternative))`, where CAFE
reports **negative** log likelihoods. A materially worse alternative is an
optimization failure, not a value to clip away. Both rate decreases and rate
increases remain in the family audit table. GO candidates require a significantly
**faster** target turnover rate (`lambda_target > lambda_background`) and the
requested sign of the reconstructed copy-number change.

Lambda is a birth/death **turnover** parameter, not a signed copy-number growth
rate. An accelerated branch can have a gain or a loss. `Gamma_change.tab`
supplies that direction and magnitude; it is not inferred from lambda's sign.
With `both`, gains and losses are reported separately, using the same family
P values and a single family BH correction before splitting by direction.
A statistically different but slower target rate is recorded as `rate_shift=slower`
and is not called an accelerated gain or loss.

### Native inputs and fitting

- `Gamma_asr.tre` supplies native branch IDs, topology and durations. Target
  identity is exact; internal nodes are not inferred from substrings. The core
  also passes the original dated tree: its clades and durations must match the
  ASR tree, and its full-precision lengths are used. Direct use of ASR outputs
  alone is supported, at the precision CAFE wrote into that file.
- `Gamma_count.tab` supplies observed tip counts. Parent/child counts are
  checked against the target's `Gamma_change.tab` entry.
- `Gamma_branch_probabilities.tab` retains its existing role in defining the
  legacy GO candidate set. These probabilities are **not** substituted for the
  target-versus-background LRT.
- A matching `Gamma_error_model.txt`, if present, is passed unchanged to all
  native fits, reconstructions and simulations. No error model is silently fitted.
- CAFE's native `-b` single-family optimization mode supplies each fitted
  lambda and likelihood, including the `-y` two-rate topology for the alternative.
  Each native invocation contains exactly one family. This uses CAFE's documented
  experimental per-family fitting facility; it does not port the model to Python.
  Two of the independent restart scores must agree near the best score. Missing
  results, an iteration cap, nonfinite values, or disagreement stop the analysis.

These comparisons use a **Base model fitted separately per family**, allowing
between-family rate variation through separate nuisance lambdas. They do not
reuse the global Gamma likelihood or attempt to estimate a gamma mixture from
one family. The original Gamma reconstruction remains the declared source of
observed direction; model comparison results do not override its signs.

### Bootstrap calibration

For every nonzero observed LRT:

1. Run native CAFE at the fitted null lambda to reconstruct the null root size.
2. Generate `go_cafe_bootstrap_replicates` families using native CAFE `-s`, that
   lambda, the same dated tree/error model, and a root distribution concentrated
   at the reconstructed null root.
3. Refit **both** models for **every** simulated family, with the same native
   optimization/restart procedure, and recompute the LRT.
4. Compute `(1 + number of simulated LRTs >= observed LRT) / (B + 1)`; ties
   include a numerical tolerance. An exactly zero observed statistic has P=1
   without simulation because all possible LRTs are nonnegative.

This is a plug-in parametric bootstrap conditional on the estimated null root,
tree and observation model. It is not an exact unconditional test: root/date/model
uncertainty and model misspecification still require calibration and sensitivity
analysis. Native CAFE's root-presence ascertainment is retained in fits and
simulations. A missing/failed simulated family is not dropped or replaced by a
favorable draw; the analysis stops with an audit record.

The minimum nontrivial P value is **1/(B+1)**, independent of the number of tree
branches. With B=999 it is 0.001. Family BH can demand finer resolution when many
families are tested; 999 is not guaranteed to be sufficient. Fix the replicate
budget before a confirmatory run. Repeatedly increasing it only for attractive
results introduces a further selection issue. The code does not use an
asymptotic chi-square fallback or relax thresholds to produce discoveries.

### GO universe and interpretation

Family BH includes every requested GO-annotated family, without preliminary
screening by target CAFE significance or reconstructed direction. No failed
family is silently removed from the testing universe or GO background.

For GO enrichment, each tested family counts **once**. Selected accelerated
gain/loss families are compared to the remaining tested annotated families.
As requested, the candidate GO IDs for each direction are those the legacy event
analysis would test for the same inputs/categories. That set is determined before
LRT selection. No new GO terms are added; a retained term with no selected families
has enrichment P=1. GO BH is performed within the declared target/direction/category
run. With `both`, the output adds a `direction` column and retains separate GO BH
corrections for gains and losses; it does not claim a joint two-direction GO FDR.

Preserving the legacy GO ascertainment rule means scientific-review item 01 is
not generally resolved by this option. Family selection power may also depend on
size/annotation, and GO hypotheses overlap. A calibrated family test followed by
GO Fisher/BH does not, by itself, establish end-to-end GO FDR control. Validate
that full selection pipeline under relevant evolutionary and annotation scenarios
before using results as confirmatory biological evidence.

### Runtime, evidence and resume

This method is substantially more expensive than reading existing branch P values.
For F distinct count patterns, B replicates and R restarts, budget roughly
`2 * F * (B+1) * R` native optimizations plus one null reconstruction/simulation
per nonzero observed statistic. Identical count patterns and completed native jobs
are cached. Estimate runtime from a small representative pilot using its native
logs; do not assume a runtime from tree branch count alone. Logs can also be large.

The CAFE executable in the validated runtime does not expose a random-seed option.
GeneGalleon does not invent one or substitute another simulator. It preserves raw
simulated tip/truth tables, every fit log, arguments, executable/input hashes and
replicate statistics. Re-running the same request reuses verified native outputs;
changing the tree, executable, error model or fitting settings uses a new cache.
Changing the adapter source also uses a new cache. Completion records verify
both native inputs and outputs; malformed or changed records fail explicitly.
This supports exact replay of the saved draws/results, not seed-based reproduction
of a fresh run. Failed attempts are preserved and do not receive a completion mark.
The native adapter holds an exclusive output-directory lock while running, so
a competing adapter cannot overwrite its cache, audit, or status. A new family
attempt clears its prior summary audit while preserving the raw native runs.

Results are written separately under `go_enrichment/cafe_lrt/`:

- `family_specificity.tsv`: native LRT, raw/adjusted P, null/background/target
  lambdas, turnover-rate shift, original copy-number change and direction,
  selection, bootstrap count/root and audit path.
- `native_cafe/family_lrt.tsv`, `metadata.json`, `branch_map.tsv`: model results,
  declared assumptions, source identity, run status and branch-to-lambda mapping.
- `native_cafe/families/*/`: observed-fit references, bootstrap generation evidence
  and all replicate LRTs. `native_cafe/runs/*/`: native inputs, commands, logs,
  reconstructions, simulated data and checksummed completion records.
- `specificity_metadata.tsv`: family/GO scope, thresholds and P-value resolution.
- `enrichment_significant_*_all_go.tsv` and `*_significant_go.tsv`: family-count
  GO tables. `n_specific_*` and `n_other_*` refer to unique tested families.
  `orthogroup_in_target` lists selected families; degenerate odds ratios are NA.
  No discoveries produce a header-only significant table, not a failed run.
- Existing `orthogroup_*significant*` tables keep their original CAFE significance
  meaning; they are not the LRT-selected family list.

Changing mode does not overwrite the event-mode outputs. Method, adapters, native
executable, CAFE inputs, dated tree, settings and primary results participate in
workflow provenance. Failed recomputation removes the prior published optional
summaries before reading new inputs, so even an input-validation error cannot
leave them looking like a current successful analysis.

## Direct use and validation

The R interface retains the original eight arguments. Optional arguments follow:

```bash
Rscript workflow/support/cafe_go_enrichment.r \
  Gamma_change.tab Gamma_branch_probabilities.tab gene_ids.tsv \
  reference.annotation.tsv outdir 'A<1>' both BP,MF,CC \
  cafe_lrt 0.05 999 5 1000 4 dated_tree.nwk
```

The additional numbers are family alpha, bootstrap replicates, fit restarts,
maximum native iterations and cores. The dated tree argument is optional for
direct runs from native output files alone. Existing eight-argument calls retain
the default event analysis.

Run in a GeneGalleon container:

```bash
python -m pytest -q workflow/tests/test_cafe_branch_specificity.py
Rscript workflow/tests/test_cafe_go_enrichment.R
python -m pytest -q --gg-strict-runtime workflow/tests/test_cafe_branch_specificity_runtime.py
```

The native runtime smoke uses one bootstrap draw to exercise the full path and
replay; it cannot demonstrate statistical power or scientific calibration.
Regression coverage checks exact branch matching, model construction, likelihood
sign, refitting both models, ties/resolution, failure handling, direction selection,
unique family counts, unchanged legacy GO behavior and artifact replay.

### Implementation validation (2026-09-10)

Docker/Linux arm64 validation passed the Python adapter tests, R family-selection
and legacy numerical/CLI regressions, and the real native CAFE gain/loss bootstrap
and replay integration. The integration used one simulated family per direction
and eight independent fits per model; all 64 native optimizations completed.
This is execution coverage, not evidence of power or calibrated GO FDR. SIF
compatibility and a production-scale bootstrap were not tested.

The pre-commit review also passed 367 Python adapter/shell/config checks and
50 CI/validation-runner checks, the R regressions and parsing of 29 R scripts,
a fresh native gain/loss/replay run,
and a native ancestral-target comparison using the default five restarts.
Failure-path regressions cover invalid inputs, corrupt completion records,
changed cached inputs, competing native writers, and stale family audits.

The CAFE correction was made in a separate owning-dependency checkout, with both
regressions failing before the fix and all 207 CAFE tests passing afterward.
Eight independent fits of the previously failing `[2, 3, 3]` tip-count dataset
then agreed within 4e-9 in negative log likelihood. No convergence threshold
was relaxed. The CAFE changes are local and have not been published upstream.
They are committed as `9784f3c` in the owning CAFE checkout. Its 207 tests now
include 591 assertions and pass in normal and randomized order.

A local development runtime containing that corrected executable is available as
`local/genegalleon:cafe-lrt-dev-20260910`. To use it with the entrypoint:

```bash
GG_CONTAINER_RUNTIME=docker \
GG_CONTAINER_DOCKER_IMAGE=local/genegalleon:cafe-lrt-dev-20260910 \
bash workflow/gg_genome_evolution_entrypoint.sh
```

Set the GO options in the entrypoint/config first. This image name identifies
the local validation build; it is not a new repository default or a published
image. Other runtimes need the owning CAFE correction before the optional method
can be considered validated there. In particular, the required native integration
check can fail in the stock runtime; it is deliberately not skipped or weakened
to accommodate the uncorrected dependency. Upstream publication and incorporation
into standard container builds remain separate release work.

Primary references:

- [CAFE5 native model, simulation, multi-lambda and per-family options](https://github.com/hahnlab/CAFE5).
- [CAFE5 native estimator implementation](https://github.com/hahnlab/CAFE5/blob/master/src/execute.cpp).
- [R Fisher exact-test documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/fisher.test.html).
- [R multiple-testing documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html).
