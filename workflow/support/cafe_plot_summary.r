# cafe_plot_summary.r
# Usage:
#   Rscript cafe_plot_summary.r Gamma_asr.tre Gamma_change.tab outdir
#
# Outputs:
#   outdir/summary_all.(pdf|svg)
#   outdir/summary_significant.(pdf|svg)
cat(as.character(Sys.time()), 'Starting cafe_plot_summary.r\n')

suppressPackageStartupMessages({
  library(ggtree)
  library(treeio)
  library(ape)
  library(ggplot2)
})

script_args <- commandArgs(trailingOnly = FALSE)
script_path_arg <- grep("^--file=", script_args, value = TRUE)
script_dir <- if (length(script_path_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_path_arg[[1]]), winslash = "/", mustWork = TRUE))
} else {
  getwd()
}
source(file.path(script_dir, "tree_label_formatting.R"))

# Input
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript cafe_plot_summary.r Gamma_asr.tre Gamma_change.tab outdir")
}
tree_file <- args[1]
change_file <- args[2]
outdir <- args[3]

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

source(file.path(script_dir, "cafe_plot_labels.R"))

# Per-OG annotation
create_annotation_df <- function(og, tree, changes_tbl) {
  fort <- ggtree::fortify(tree)
  keys <- norm_key(fort$label)

  changes_og <- changes_tbl[changes_tbl$FamilyID == og, , drop = FALSE]
  if (nrow(changes_og) == 0) {
    return(data.frame(
      isTip = fort$isTip,
      x = fort$x,
      y = fort$y,
      raw_lbl = fort$label,
      key = keys,
      species = ifelse(fort$isTip, get_species(fort$label), NA_character_),
      has_star = grepl("\\*", fort$label),
      change = rep(NA_real_, length(keys)),
      stringsAsFactors = FALSE
    ))
  }

  dlt <- as.numeric(changes_og[1, keys, drop = TRUE])

  data.frame(
    isTip = fort$isTip,
    x = fort$x,
    y = fort$y,
    raw_lbl = fort$label,
    key = keys,
    species = ifelse(fort$isTip, get_species(fort$label), NA_character_),
    has_star = grepl("\\*", fort$label),
    change = dlt,
    stringsAsFactors = FALSE
  )
}

# Accumulate per-branch counts across all families
accumulate_branch_counts <- function(trees, changes_tbl, ogs, ref_tree) {
  ref_fort <- ggtree::fortify(ref_tree)
  base <- data.frame(
    key = norm_key(ref_fort$label),
    isTip = ref_fort$isTip,
    x = ref_fort$x,
    y = ref_fort$y,
    stringsAsFactors = FALSE
  )
  agg <- base
  agg$inc <- 0L
  agg$dec <- 0L
  agg$inc_sig <- 0L
  agg$dec_sig <- 0L

  if (anyNA(agg$key) || anyDuplicated(agg$key) || anyDuplicated(changes_tbl$FamilyID) ||
      anyDuplicated(names(trees)) || !all(agg$key %in% names(changes_tbl))) {
    stop("CAFE summary requires unique family and branch identities and matching change columns.")
  }
  selected <- ogs[ogs %in% names(trees)]
  if (!length(selected)) return(agg)
  values <- as.matrix(changes_tbl[match(selected, changes_tbl$FamilyID), agg$key, drop = FALSE])
  if (!is.numeric(values)) stop("CAFE branch changes must be numeric.")
  agg$inc <- as.integer(colSums(values > 0, na.rm = TRUE))
  agg$dec <- as.integer(colSums(values < 0, na.rm = TRUE))

  # Geometry is needed only for the reference tree. Native CAFE trees share
  # branch identities; their labels carry the family-specific significance flags.
  for (i in seq_along(selected)) {
    tree <- trees[[selected[i]]]
    labels <- c(tree$tip.label, tree$node.label)
    keys <- norm_key(labels)
    if (anyNA(keys) || anyDuplicated(keys) || !setequal(keys, agg$key)) {
      stop("CAFE family tree branch identities differ from the reference: ", selected[i])
    }
    positions <- match(keys[grepl("*", labels, fixed = TRUE)], agg$key)
    change <- values[i, positions]
    agg$inc_sig[positions] <- agg$inc_sig[positions] + as.integer(!is.na(change) & change > 0)
    agg$dec_sig[positions] <- agg$dec_sig[positions] + as.integer(!is.na(change) & change < 0)
  }
  agg
}

# Plot helper
build_base_plot <- function(ref_tree, ref_fort, xmax, y_step) {
  tip_df <- data.frame(
    isTip = ref_fort$isTip,
    x = ref_fort$x,
    y = ref_fort$y,
    label = ref_fort$label,
    stringsAsFactors = FALSE
  )
  tip_df <- tip_df[tip_df$isTip, , drop = FALSE]
  tip_df$species <- get_species(tip_df$label)
  tip_df$species_x <- tip_df$x + 0.015 * xmax
  tip_df <- tip_df[!is.na(tip_df$species) & tip_df$species != "", , drop = FALSE]

  ggtree::ggtree(ref_tree) +
    ggtree::theme_tree2() +
    theme(
      plot.margin = margin(10, 70, 10, 10),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank()
    ) +
    coord_cartesian(clip = "off") +
    geom_text(
      data = tip_df,
      aes(x = species_x, y = y, label = species),
      hjust = 0, size = font_size * font_size_factor, fontface = "italic", inherit.aes = FALSE,
      na.rm = TRUE
    ) +
    ggtree::geom_treescale(
      x = xmax * 0.02, y = - 0.5 * y_step,
      fontsize = font_size * font_size_factor
    )
}

# -- Main --
changes_tbl <- read_tsv_base(change_file)
colnames(changes_tbl)[1] <- "FamilyID"
trees <- ape::read.nexus(tree_file)

ogs <- intersect(changes_tbl$FamilyID, names(trees))
if (length(ogs) == 0) stop("No overlapping OG IDs between tree file and tables")

# Reference tree
ref_tree <- trees[[ogs[1]]]
ref_fort <- ggtree::fortify(ref_tree)

# Geometry helpers
xmax <- suppressWarnings(max(ref_fort$x[is.finite(ref_fort$x)], na.rm = TRUE)); if (!is.finite(xmax)) xmax <- 1
tip_y <- sort(unique(ref_fort$y[ref_fort$isTip & is.finite(ref_fort$y)]))
y_step <- if (length(tip_y) >= 2) median(diff(tip_y)) else 1
if (!is.finite(y_step) || y_step <= 0) y_step <- 1

# Output format
width <- 7.2
height <- 3 + Ntip(ref_tree) * font_size / 72
dpi <- 300

agg <- accumulate_branch_counts(trees, changes_tbl, ogs, ref_tree)

# Plot
base <- build_base_plot(ref_tree, ref_fort, xmax, y_step) + xlim(0, xmax * 1.35)
y_off <- y_step * font_size_factor * 0.8

df_all <- agg
df_all$x_lab <- df_all$x + ifelse(df_all$isTip, 0, -0.45)
df_all$y_lab <- df_all$y + y_off
df_all$lbl <- sprintf("+%d/-%d", df_all$inc, df_all$dec)
p_all <- base +
  ggtitle("Summary of Gene Family Expansion/Contraction") +
  theme(plot.title = element_text(hjust = 0.5, size = font_size * 1.5, face = "bold")) +
  geom_text(data = df_all, aes(x = x_lab, y = y_lab, label = lbl), hjust = 1, vjust = 0.5, color = "purple", size = font_size * font_size_factor * 0.8, inherit.aes = FALSE, na.rm = TRUE)

df_sig <- agg
df_sig$x_lab <- df_sig$x + ifelse(df_sig$isTip, 0, -0.45)
df_sig$y_lab <- df_sig$y + y_off
df_sig$lbl <- sprintf("+%d/-%d", df_sig$inc_sig, df_sig$dec_sig)
p_sig <- base +
  ggtitle("Summary of Significant Gene Family Expansion/Contraction") +
  theme(plot.title = element_text(hjust = 0.5, size = font_size * 1.5, face = "bold")) +
  geom_text(data = df_sig, aes(x = x_lab, y = y_lab, label = lbl), hjust = 1, vjust = 0.5, color = "purple", size = font_size * font_size_factor * 0.8, inherit.aes = FALSE, na.rm = TRUE)

# Save
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
ggsave(file.path(outdir, "summary_all.pdf"), plot = p_all, width = width, height = height, dpi = dpi)
ggsave(file.path(outdir, "summary_all.svg"), plot = p_all, width = width, height = height, dpi = dpi)
ggsave(file.path(outdir, "summary_significant.pdf"), plot = p_sig, width = width, height = height, dpi = dpi)
ggsave(file.path(outdir, "summary_significant.svg"), plot = p_sig, width = width, height = height, dpi = dpi)

cat(as.character(Sys.time()), 'cafe_plot_summary.r completed successfully. Exiting\n')
