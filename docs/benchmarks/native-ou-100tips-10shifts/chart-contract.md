# Benchmark chart contract

Question: at 100 tips and ten planted shifts, what computation and completed
selection accuracy do the two current methods deliver?

Surface: reproducible scientific figure (Matplotlib PNG and SVG) displayed in the
conversation, with machine-readable paired results alongside it. Static small
multiples show the requested metrics simultaneously. No performance conclusion
is prespecified: annotate only observed results after all runs finish.

Grain: six paired datasets, three seeds at each of two effect scales. Preserve
individual points. Use method colors plus marker shapes; distinguish censored
runs with upward triangles. Display wall seconds (log scale if needed), peak
RSS MiB, exact-branch precision/recall/F1, and true tip-mean RMSE. Label accuracy
as conditional on completion; explicitly report success denominators. Do not
assign zero accuracy to a timeout or substitute the uncalibrated native model.

Palette: accessible blue/orange method colors, neutral grid and labels. Export
at readable laptop width; inspect PNG for label collisions, censoring markers,
units and completion denominators. No confidence interval or general accuracy
claim from three replicates per condition. Save all point-level values as CSV.
