# Tree Annotation Helpers

This directory contains the local R helper code used by GeneGalleon's tree
annotation and plotting scripts.

Current consumers include:

- `workflow/support/stat_branch2tree_plot.r`
- `workflow/support/annotation_summary.r`

Notes:

- the helpers are sourced directly from `workflow/support/tree_annotation/R/`
  during workflow runs; they are not treated as a separately installed R
  package in normal GeneGalleon execution,
- the package metadata files (`DESCRIPTION`, `NAMESPACE`) remain here for local
  development and editor tooling,
- when updating plotting behavior for branch/stat annotation, start with
  `R/main.R` in this directory and the caller script that sources it.
