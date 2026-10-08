# HGT PDF rendering performance

The focused HGT renderer replays the saved gene-evolution plot settings for page
one and adds donor/recipient genomic context on page two. Rendering existing
results does not require rerunning tree inference or HGT detection.

Profiling large saved-input plots identified four expensive operations:

- Constructing a data frame for every gene-by-intron-site cell.
- Scanning every later gene row for each intron connection.
- Constructing a data frame for every alignment run and domain interval, and
  repeatedly ordering identical domain sets along a sequence.
- Scanning all compressed noncoding gaps for every plotted genomic coordinate.

The renderer now assembles cells and intervals in batches, indexes nearby intron
occurrences, reuses ordered domain sets, and uses prefix offsets plus binary
search for compressed coordinates. It preserves cell states, phase ambiguity,
reciprocal connection rules, interval ordering, coordinates and missing-data
diagnostics. Context source verification remains enabled at the original points.

## Verification and measurement

Compare frozen and current helpers in the same dependency runtime, with the same
saved branch tables, alignments, domains, genomic models, annotations, event
selection and plot settings. Run timings serially, exclude warmups, and retain
individual wall times and peak RSS. Use profiling to locate expensive functions;
distinguish profiled measurements from ordinary throughput measurements.

Require exact equality of intermediate objects and complete audit tables. Compare
both PDF pages after removing only document creation/modification timestamps;
also compare the merged PDF bytes. Preserve input digests and renderer parameters
with the private benchmark outputs rather than committing research files.

Regression tests cover intron phase/gap states, sparse cell grids, first-row and
reciprocal-nearest ties, domain stacking and gap-separated alignment runs. Context
tests check exact floating-point results at gap boundaries, large coordinates,
numeric-type behavior, and PDF/audit equality against the scalar transform. Batch
tests retain failed-output isolation and source-change rejection.

Use the existing check entrypoint and runtime freshness checks described in
[Development and Tests](development-and-tests.md#choose-checks-for-a-change).
An isolated package library in a pinned comparison image is useful for measuring
the change; fresh production-image validation is a separate requirement.

The 2026-10-08 saved-input benchmark and its raw trials are kept outside version
control. These measurements describe rendering kernels and representative PDFs;
they do not establish total runtime for a new cohort, remote acquisition speed,
or SIF compatibility.
