# Chart contract

Question: compare the runtime, process memory and shift-recovery accuracy of seven configured methods on six identical simulated datasets.
Takeaway: report the measured speed/accuracy tradeoffs after all runs complete; no prespecified winner.
Surface: standalone scientific PNG and SVG, rendered with Matplotlib, with complete metrics CSV and reproducible source.
Family: horizontal categorical dot plots, six panels (time, RSS, precision, recall, F1, expected-mean RMSE). Seven named method rows and two effect-strength series; three independently simulated seeds per series.
Data sufficiency: 42 planned runs. Show every completed observation; no boxplot or inferred confidence interval with n=3. Identify failures/timeouts and completion denominators separately.
Scales: time logarithmic and labelled; RSS/RMSE start at zero; precision/recall/F1 use the same 0..1 limits. Every panel repeats method labels and order.
Palette: hard two-root cap, blue #0072B2 and orange #D55E00, representing the two effect scales. Circular versus triangular markers and vertical offsets preserve distinction in grayscale. Neutral grids and statistic ticks.
Summary: median time/RSS and mean accuracy, labelled in the caption. CPU time, selected counts and partition ARI remain available in tables/CSV.
Footprint: 15 by 12 inches; 180-dpi PNG and editable SVG. Inspect actual exported PNG for overlap, clipping and legibility before delivery.
Limits: one tree, one trait, 10-shift cap equal to truth, no convergence/missing/error. Bootstrap B=19 is coarse; scores and search/optimizer defaults differ. No global-null or production reliability claim.
