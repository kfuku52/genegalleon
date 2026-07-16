# Tree Visualization Helpers

This directory contains the local R helper code used by GeneGalleon's tree
visualization and plotting scripts. The package name is `genegalleon.treevis`.

Current consumers include:

- `workflow/support/stat_branch2tree_plot.r`
- `workflow/support/annotation_summary.r`

Notes:

- the package is installed into the GeneGalleon container image during the
  image build,
- workflow consumers load it with `library(genegalleon.treevis)`; running those
  consumers outside the GeneGalleon runtime therefore requires installing this
  package and its declared dependencies first,
- `R/` contains the implementation modules and `NAMESPACE` is the public API
  boundary,
- validate package changes with `workflow/tests/check_treevis_package.sh` and
  `workflow/tests/test_treevis_main.R` inside a GeneGalleon runtime.
