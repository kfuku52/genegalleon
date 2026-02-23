# cafe_plot_branch_id.r
# Usage:
#   Rscript cafe_plot_branch_id.r Gamma_asr.tre outdir
#
# Outputs:
#   outdir/branch_id.(pdf|svg)
cat(as.character(Sys.time()), 'Starting cafe_plot_branch_id.r\n')

suppressPackageStartupMessages({
  library(ggtree)
  library(ape)
  library(ggplot2)
})

# Input
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript cafe_plot_branch_id.r Gamma_asr.tre outdir")
}
tree_file <- args[1]
outdir <- args[2]

trees <- ape::read.nexus(tree_file)
tree <- trees[[1]]

# Trim _~ suffix from tip labels and node labels
tree$tip.label <- sub(">.*$", ">", tree$tip.label)
tree$node.label <- sub(">.*$", ">", tree$node.label)

# Visualize the CAFE tree with branch IDs
ref_fort <- ggtree::fortify(tree)
xmax <- suppressWarnings(max(ref_fort$x[is.finite(ref_fort$x)], na.rm = TRUE)); if (!is.finite(xmax)) xmax <- 1
tip_y <- sort(unique(ref_fort$y[ref_fort$isTip & is.finite(ref_fort$y)]))
y_step <- if (length(tip_y) >= 2) median(diff(tip_y)) else 1
if (!is.finite(y_step) || y_step <= 0) y_step <- 1

font_size <- 8
font_size_factor <- 0.352777778 # https://stackoverflow.com/questions/17311917/ggplot2-the-unit-of-size
p <- ggtree(tree) +
  coord_cartesian(clip = "off") +
  xlim(0, xmax * 1.35) +
  geom_text2(aes(subset = !isTip, label = label), hjust = -0.1, size = font_size * font_size_factor, color = "purple") +
  geom_tiplab(size = font_size * font_size_factor, color = "purple")

# Save
width <- 7.2
height <- 3 + Ntip(tree) * font_size / 72
dpi <- 300
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
ggsave(file.path(outdir, "branch_id.pdf"), plot = p, width = width, height = height, dpi = dpi)
ggsave(file.path(outdir, "branch_id.svg"), plot = p, width = width, height = height, dpi = dpi)
cat(as.character(Sys.time()), 'cafe_plot_branch_id.r completed successfully. Exiting\n')
