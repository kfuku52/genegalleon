# cafe_plot_each_family.r
# Usage:
#   Rscript cafe_plot_each_family.r Gamma_asr.tre Gamma_count.tab Gamma_change.tab outdir core
cat(as.character(Sys.time()), 'Starting cafe_plot_each_family.r\n')

suppressPackageStartupMessages({
  library(ggtree)
  library(treeio)
  library(ape)
  library(ggplot2)
})

# Input
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript cafe_plot_each_family.r Gamma_asr.tre Gamma_count.tab Gamma_change.tab outdir core")
}
tree_file <- args[1]
count_file <- args[2]
change_file <- args[3]
outdir <- args[4]
core <- as.integer(args[5])

# Font size
font_size <- 8
font_size_factor <- 0.352777778 # https://stackoverflow.com/questions/17311917/ggplot2-the-unit-of-size

read_tsv_base <- function(path, na = character()) {
  read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = na
  )
}

# Normalize labels:
# - drop "_{copynumber}"
# - strip "*"
norm_key <- function(lbl) {
  x <- sub("_[^_]+$", "", lbl)
  gsub("\\*", "", x)
}

# Extract species name for tip labels: remove "*" and "_{copynumber}", then drop "<...>"
get_species <- function(lbl) {
  x <- gsub("\\*", "", lbl)
  x <- sub("_[^_]+$", "", x)
  x <- sub("<.*$", "", x)
  gsub("_", " ", x)
}

# Per-OG annotation
create_annotation_df <- function(og, tree, counts_tbl, changes_tbl) {
  fort <- ggtree::fortify(tree)
  keys <- norm_key(fort$label)

  counts_og <- counts_tbl[counts_tbl$FamilyID == og, , drop = FALSE]
  changes_og <- changes_tbl[changes_tbl$FamilyID == og, , drop = FALSE]

  if (nrow(counts_og) == 0 || nrow(changes_og) == 0) {
    return(data.frame(
      isTip = fort$isTip,
      x = fort$x,
      y = fort$y,
      raw_lbl = fort$label,
      key = keys,
      species = ifelse(fort$isTip, get_species(fort$label), NA_character_),
      has_star = grepl("\\*", fort$label),
      count = rep(NA_real_, length(keys)),
      change = rep(NA_real_, length(keys)),
      count_s = rep(NA_character_, length(keys)),
      change_s = rep(NA_character_, length(keys)),
      stringsAsFactors = FALSE
    ))
  }

  cnt <- as.numeric(counts_og[1, keys, drop = TRUE])
  dlt <- as.numeric(changes_og[1, keys, drop = TRUE])

  data.frame(
    isTip = fort$isTip,
    x = fort$x,
    y = fort$y,
    raw_lbl = fort$label,
    key = keys,
    species = ifelse(fort$isTip, get_species(fort$label), NA_character_),
    has_star = grepl("\\*", fort$label),
    count = cnt,
    change = dlt,
    count_s = as.character(cnt),
    change_s = sprintf("%+d", as.integer(dlt)),
    stringsAsFactors = FALSE
  )
}

# Plot one OG
plot_cafe5_each_family <- function(og, tree, counts_tbl, changes_tbl, outdir, width, height, dpi) {
  ann_df <- create_annotation_df(og, tree, counts_tbl, changes_tbl)

  # Base plot
  p <- ggtree::ggtree(tree) +
       ggtree::theme_tree2() +
       theme(plot.margin = margin(10, 70, 10, 10), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) +
       coord_cartesian(clip = "off")

  # Geometry helpers
  xmax <- suppressWarnings(max(ann_df$x[is.finite(ann_df$x)], na.rm = TRUE)); if (!is.finite(xmax)) xmax <- 1
  ymax <- suppressWarnings(max(ann_df$y[is.finite(ann_df$y)], na.rm = TRUE)); if (!is.finite(ymax)) ymax <- 1

  # Compute vertical step between adjacent tips
  tip_y <- sort(unique(ann_df$y[ann_df$isTip & is.finite(ann_df$y)]))
  y_step <- if (length(tip_y) >= 2) median(diff(tip_y)) else 1
  if (!is.finite(y_step) || y_step <= 0) y_step <- 1

  # Tip species names
  tips <- ann_df[ann_df$isTip, , drop = FALSE]
  tips$species_x <- tips$x + 0.015 * xmax
  p <- p + geom_text(data = tips, aes(x = species_x, y = y, label = species), hjust = 0, size = font_size * font_size_factor, fontface = "italic", inherit.aes = FALSE)

  # Branch annotations
  edge_labels <- ann_df
  edge_labels$y_off <- y_step * font_size_factor * 0.8
  edge_labels$x_off <- ifelse(edge_labels$isTip, 0, -0.45)
  edge_labels$x_lab <- edge_labels$x + edge_labels$x_off
  edge_labels$y_lab <- edge_labels$y + edge_labels$y_off
  edge_labels$color_group <- ifelse(
    is.na(edge_labels$change),
    "other",
    ifelse(edge_labels$change > 0, "inc", ifelse(edge_labels$change < 0, "dec", "other"))
  )
  edge_labels$label_num <- paste0(ifelse(edge_labels$has_star, "*", ""), edge_labels$change_s, " (", edge_labels$count_s, ")")

  p <- p + geom_text(data = edge_labels, aes(x = x_lab, y = y_lab, label = label_num, color = color_group), hjust = 1, vjust = 0.5, size = font_size * font_size_factor * 0.8, inherit.aes = FALSE, show.legend = FALSE) +
    scale_color_manual(values = c(inc = "red", dec = "blue", other = "black"))

  # Title
  p <- p + ggtitle(og) +
    theme(plot.title = element_text(hjust = 0.5, size = font_size * 1.5, face = "bold"))

  # Branch scale
  p <- p + ggtree::geom_treescale(
    x = xmax * 0.02, y = - 0.5 * y_step,
    fontsize = font_size * font_size_factor
  )
  
  # Legend
  leg_lines <- c("Increase", "Decrease", "No change", "* Significant change")
  max_tip_y <- if (length(tip_y)) max(tip_y) else ymax
  line_gap <- 0.9 * y_step
  leg_base_y <- max_tip_y + 0.6 * y_step
  leg_y <- leg_base_y + rev(seq_along(leg_lines) - 1) * line_gap
  leg_x <- 0

  p <- p +
    annotate("text", x = leg_x, y = leg_y[1], label = leg_lines[1], color = "red", hjust = 0, vjust = 1, size = font_size * font_size_factor) +
    annotate("text", x = leg_x, y = leg_y[2], label = leg_lines[2], color = "blue", hjust = 0, vjust = 1, size = font_size * font_size_factor) +
    annotate("text", x = leg_x, y = leg_y[3], label = leg_lines[3], color = "black", hjust = 0, vjust = 1, size = font_size * font_size_factor) +
    annotate("text", x = leg_x, y = leg_y[4], label = leg_lines[4], color = "black", hjust = 0, vjust = 1, size = font_size * font_size_factor)

  p <- p + xlim(0, xmax * 1.35)

  # Save
  outfile_pdf <- file.path(outdir, paste0(og, ".pdf"))
  outfile_svg <- file.path(outdir, paste0(og, ".svg"))
  ggsave(outfile_pdf, plot = p, width = width, height = height, dpi = dpi)
  ggsave(outfile_svg, plot = p, width = width, height = height, dpi = dpi)
}

# -- Main --
counts_tbl <- read_tsv_base(count_file)
colnames(counts_tbl)[1] <- "FamilyID"
changes_tbl <- read_tsv_base(change_file)
colnames(changes_tbl)[1] <- "FamilyID"
trees <- ape::read.nexus(tree_file)
ogs <- intersect(counts_tbl$FamilyID, names(trees))
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

cl <- parallel::makeCluster(core, type = "PSOCK")
parallel::clusterExport(
  cl,
  varlist = c("plot_cafe5_each_family", "create_annotation_df", "norm_key", "get_species", "font_size",
              "font_size_factor", "trees", "counts_tbl", "changes_tbl", "outdir"),
  envir = environment()
)
  parallel::clusterEvalQ(
  cl,
  suppressPackageStartupMessages({
    library(ggplot2)
    library(ggtree)
    library(treeio)
    library(ape)
  })
)
invisible(
  parallel::parLapply(cl, ogs, function(og) {
    tree <- trees[[og]]
    width <- 7.2
    height <- 3 + ape::Ntip(tree) * font_size / 72
    dpi <- 300
    plot_cafe5_each_family(og, tree, counts_tbl, changes_tbl, outdir, width, height, dpi)
    NULL
  })
)
parallel::stopCluster(cl)

cat(as.character(Sys.time()), 'cafe_plot_each_family.r completed successfully. Exiting\n')
