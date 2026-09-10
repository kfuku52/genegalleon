# Native AIC improvement experiment

Authorized after the six paired 100-tip/10-shift comparisons. This work changes
NWKIT only if independent validation supports adoption. Existing results remain
immutable. No new information-criterion penalty or alpha restriction is chosen
using the planted truth.

## Diagnostic gate

Refit the NWKIT, kfl1ou and planted layouts in both implementations. Check native
likelihood at kfl1ou's fitted alpha independently, distinguish alpha bounds from
location search, and audit original-pool coverage. The planted layout is a
reference only; never pass it to candidate generation.

## Development ablations

Use the original six inputs only. Separate covariance-scale profiling in quick
ranking, IC-prioritized iterative local refinement, and updated-covariance branch
screening. Keep the initial pool/refit/screen/beam budgets 24/48/2000/2. Then
consider larger budgets if the low-cost variants are promising. Record all
variants, including negative results. These are development comparisons, not
unbiased performance estimates.

## Independent comparison

Freeze the selected implementation before inspecting new results. Use seeds
28101–28130, one balanced 100-tip tree and effect strengths 0, 2 and 6, otherwise
the original independent branch-recursion simulation. Strength-zero datasets
have zero true shifts (the random draw order remains paired). Maximum shifts
remains 10 in all three procedures. This tests an unknown count below the cap
for null data, in addition to the original nonnull design.

Compare frozen original NWKIT AIC, candidate NWKIT AIC and kfl1ou AIC on exactly
the same inputs in sequential fresh processes, one BLAS/OpenMP thread, cyclic
method order. Retain errors/timeouts and wall time, CPU time, peak process RSS,
selected clades, fitted tip means and score. No concurrent numerical work during
these timings. Runtime cap: 1800 seconds per process.

Primary criteria are paired F1 and fitted-mean RMSE at strengths 2 and 6;
secondary are shift count, global-null false-positive frequency, time and peak
RSS. Report means and paired bootstrap uncertainty over seeds, pairing strengths
when aggregating. No inferential F1 metric for the global null; report any-shift
frequency and count instead. Runtime is a configured-procedure comparison, not
an equivalent-output kernel-speed claim.

Adoption requires a meaningful accuracy gain over original NWKIT without clear
regression in the other nonnull condition or null behavior. To describe the
candidate as comparable to kfl1ou, use provisional engineering margins: F1
difference >= -0.02 and RMSE ratio <= 1.05; report uncertainty and do not turn
point-estimate margins into a claim of established statistical equivalence.
Prefer runtime in the kfl1ou range, measuring the tradeoff rather than assuming
it. If these criteria fail, preserve the baseline default and document what
failed. Any tuning after opening independent results requires new validation
seeds. Additional tree shapes/counts are required before broad generalization.

## Frozen candidate decision (before independent runs)

The joint optimum-increment path was selected from development data. Its six-case
weak/strong mean F1 is 0.367/0.433 and mean RMSE 0.653/0.595. This is a development
result, not the independent estimate. The 150-iteration and 300-iteration pilot
settings yielded the same selected layouts; retain the existing 150 default.
The production implementation reproduces the prototype's selected layouts.

The quick-scale-only and refinement-order variants did not consistently help.
Covariance-updated branch pools improved one weak case but left substantial mean
error. All pilot outputs are retained. No larger refit budget was adopted.

The final path uses two covariance rounds, 80 penalty points per round before
early termination, all branches with unstandardized OU-effect columns, up to
150 proximal iterations per point, cap 10 and total refit budget 48 including
the covariance-only seed. Screening budget 2000 caps lambda points; this work
unit differs from the baseline's 2000 quick layout scores. The time/memory
measurements compare the actual configured procedures.

## Prespecified scope check

After the primary experiment, run a small extension with five true shifts below
the unchanged cap of ten. Use seeds 29101–29103, effect scales 2 and 6, and both
the original balanced topology and independently generated random ultrametric
100-tip topologies. Shifted clades are disjoint and contain 3–10 tips. These 12
additional datasets (36 runs) are exploratory scope checks, not a replacement
for the 30-seed primary comparison. Do not tune from them.
