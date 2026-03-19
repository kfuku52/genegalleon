# Tree Visualization Helpers

This directory contains the local R helper code used by GeneGalleon's tree
visualization and plotting scripts. The package name is `genegalleon.treevis`.

Current consumers include:

- `workflow/support/stat_branch2tree_plot.r`
- `workflow/support/annotation_summary.r`

Notes:

- `R/main.R` is a loader that sources the split module files in this
  directory,
- the helpers are sourced directly from `workflow/support/treevis/R/` during
  workflow runs; the legacy `workflow/support/tree_annotation/` path is kept as
  a compatibility symlink,
- the package metadata files (`DESCRIPTION`, `NAMESPACE`) remain here for local
  development and editor tooling,
- when updating plotting behavior for branch/stat annotation, start with the
  module file that owns the function and the caller script that sources it.
