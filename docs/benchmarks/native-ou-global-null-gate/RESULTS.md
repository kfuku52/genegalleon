# Global-null AIC gate: reduced null detections, substantial power cost

The optional gate reduces false detections in this development pilot, but is
not suitable for automatic adoption. It tests whether any shifts exist and then
retains the entire AIC winner; it does not correct the selected shift count or
test individual branches. Code and outputs are implemented; statistical adoption
is pending.

## Diagnosis

On all 30 no-shift, 100-tip fixtures, ordinary native-path AIC with estimated
alpha selects ten shifts. Fixing alpha at the simulation truth (alpha*height=3)
still selects shifts in 28/30 datasets, with mean K=8.73. Alpha estimation alone
therefore cannot explain the false detections. This comparison does not decompose
the effects of process-variance estimation and searching many branch layouts.

## Paired B=19 pilot

All 48 paired cases completed, replaying the full search for all 19 bootstrap
draws in every case. At nominal level 0.05, null false detections fall from 30/30
to 3/30. The observed 10% rate has a 95% Wilson interval of 3.46–25.62%; **5%
control has not been demonstrated**. Those three rejections each have p=0.05,
the minimum p-value resolvable with 19 draws. Mean null K falls from 10 to 1,
and mean tip-mean RMSE falls from 0.717 to 0.180.

The following are exact-clade recovery metrics. Precision and recall pool
true/false positives across datasets; no-selection outcomes contribute zero
true positives to recall. K=5 pools balanced and random trees (three seeds each);
K=10 uses the first three main-study seeds. These subsets are small and were
chosen before running the gate, not to estimate population power accurately.

| Truth / effect | n | Any shifts, AIC → gate | Precision, AIC → gate | Recall, AIC → gate | Mean false branches, AIC → gate | Mean RMSE, AIC → gate |
|---|---:|---:|---:|---:|---:|---:|
| No shifts | 30 | 30 → 3 | — | — | 10.00 → 1.00 | 0.717 → 0.180 |
| 5 shifts / 2 | 6 | 6 → 1 | 0.237 → 0.000 | 0.467 → 0.000 | 7.50 → 1.67 | 0.789 → 0.646 |
| 5 shifts / 6 | 6 | 6 → 5 | 0.450 → 0.440 | 0.900 → 0.733 | 5.50 → 4.67 | 0.715 → 0.866 |
| 10 shifts / 2 | 3 | 3 → 1 | 0.367 → 0.500 | 0.367 → 0.167 | 6.33 → 1.67 | 0.690 → 1.013 |
| 10 shifts / 6 | 3 | 3 → 2 | 0.625 → 0.600 | 0.500 → 0.400 | 3.00 → 2.67 | 1.413 → 1.625 |

For weak five-shift data, the sole retained model has ten entirely false branch
locations. For strong five-shift data, mean false branches still exceed four
after gating. Fewer selected datasets must not be interpreted as accurate branch
localization. The gate often sacrifices true-shift detection, and prediction
error worsens in three of the four nonnull strata.

Mean observed search time was 1.55 seconds and mean additional gate time was
28.71 seconds for B=19. These timings include four-worker execution and overlap
with tests/default-count checks. They indicate added computational cost only;
this is not a controlled speed comparison and peak memory was not measured.

## Default B=199 representative checks

The first main fixture at each effect (seed 28101) also completed with all 199
inner draws and newly estimated covariance in every search:

| Truth / effect | Bootstrap p | AIC K | Gated K | Additional gate wall seconds |
|---|---:|---:|---:|---:|
| No shifts | 0.200 | 10 | 0 | 327.5 |
| 10 shifts / 2 | 0.010 | 10 | 10 | 287.1 |
| 10 shifts / 6 | 0.005 | 10 | 10 | 277.5 |

These confirm default-count execution and the expected accept/reject outputs
on three representative cases. They cannot estimate default-count FPR or power.
All bootstrap p-values were independently recomputed from saved statistics.
Together with the pilot, 1,509 inner full-search replicates completed without
failures (48×19 + 3×199). Timing caveats above apply.

## Interpretation and next evidence

Keep `--global-null-gate` opt-in. Defaults and the AIC formula are unchanged.
The calibrated decision, null-selected predictions, resume configuration, and
nested support replay are included in the implementation. The test remains a
plug-in parametric bootstrap with no uniform composite-null guarantee.

The next adoption study needs many independent outer simulations with an
adequate inner draw count, varied true alpha/variance/tree shapes and
measurement patterns, and explicit lower bounds on useful-shift recovery.
Conditional over-selection requires a separate complexity or branch-selection
procedure; the present gate cannot supply it. No penalty or threshold was tuned
from this pilot.

See [protocol](PROTOCOL.md), [raw summary](summary.json), and
[independent metric audit](independent-metric-qa.json). The latter reconstructs
all clade counts and RMSEs from previously archived AIC predictions and fresh
no-shift fits. This study reuses development fixtures; it is not held-out evidence.
