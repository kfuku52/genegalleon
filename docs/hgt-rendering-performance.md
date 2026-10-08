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
diagnostics.

Alignment domain keys use an interval sweep over the existing monotone
trimmed-to-untrimmed mapping. Active counts retain repeated overlapping labels;
inclusive endpoints, missing mappings and the final stacking order are unchanged.
Intron event rows and cell columns are assembled in batches, and a nongap mask is
reused without removing nucleotide validation. Domain sweeps use indexed query
groups and parallel event vectors while retaining stable boundary ordering.
Context gap types use a sorted span index and prefix maximum ends to test exact
containment and overlap, preserving each locus rather than merging loci.

Focused exports now send up to eight OGs at a time to the existing
`tree_plot_batch.r` worker. Each job evaluates gene evolution's shared
`stat_branch2tree_plot.r` source with its own saved settings and inputs; there is
no separate tree renderer. The worker restores plot options, theme, species
parser and random state between jobs, releases each job's data, and preserves
input/source mutation checks and failed-output isolation. Context and annotation
caches remain shared across the export. Event-gene links are indexed by family.

Each context page verifies the exact annotation, classification, taxonomy and
neighbor-family inputs it uses, including cached dependencies. The final export
still hashes every accumulated input before completing. This removes repeated
reads of unrelated earlier pages without replacing content hashes with file
timestamps or treating a batch receipt as proof of workflow completion.

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
reciprocal-nearest ties, domain stacking and gap-separated alignment runs. An
independent point-membership oracle checks inclusive domain endpoints, repeated
overlaps, missing mappings and complete plot layers. Context tests check exact
floating-point results at gap boundaries, large coordinates, numeric-type
behavior, gap classification against individual loci, and PDF/audit equality
against the scalar transform. Batch tests retain failed-output isolation and
source-change rejection.

`workflow/benchmarks/benchmark_focused_hgt_render.py` compares complete saved-input
exports, including the native tree PDFs, merged two-page PDFs, ordered audit
tables, replayed settings and consumed-input hashes. It runs warmups and alternating
trials in one dependency runtime. Supply frozen support directories and a private
plan containing each family's `events.tsv` and `links.tsv`:

```bash
bash workflow/tests/run_in_runtime.sh python workflow/benchmarks/benchmark_focused_hgt_render.py \
  --inputs /path/to/saved-inputs --event-plan /path/to/family-requests \
  --families OG0001 OG0002 --baseline-root /path/to/baseline/support \
  --candidate-root /path/to/candidate/support --trials 3 --output /path/to/comparison
```

Inputs must be available inside the selected runtime. Peak Python RSS and the
largest R child's RSS are reported separately; they are not summed. Regression
tests also check family/event isolation across batch boundaries and refusal of
missing, mismatched or failed worker receipts.

Use the existing check entrypoint and runtime freshness checks described in
[Development and Tests](development-and-tests.md#choose-checks-for-a-change).
An isolated package library in a pinned comparison image is useful for measuring
the change; fresh production-image validation is a separate requirement.

The 2026-10-08 saved-input benchmark and its raw trials are kept outside version
control. These measurements describe rendering kernels and representative PDFs;
they do not establish total runtime for a new cohort, remote acquisition speed,
or SIF compatibility.
