args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- NULL
for (arg in args) {
  if (startsWith(arg, file_arg)) {
    script_path <- substring(arg, nchar(file_arg) + 1)
    break
  }
}

if (is.null(script_path) || nchar(script_path) == 0) {
  stop("Could not determine test script path from commandArgs().")
}

resolve_script_path <- function(path_in) {
  candidates <- unique(c(
    path_in,
    gsub("~\\+~", " ", path_in),
    gsub("%20", " ", path_in)
  ))
  for (candidate in candidates) {
    if (file.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }
  return(NA_character_)
}

resolved_script_path <- resolve_script_path(script_path)
repo_root <- NA_character_

if (!is.na(resolved_script_path)) {
  script_dir <- dirname(resolved_script_path)
  repo_root_candidate <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = FALSE)
  if (dir.exists(file.path(repo_root_candidate, "workflow", "script"))) {
    repo_root <- normalizePath(repo_root_candidate, winslash = "/", mustWork = TRUE)
  }
}

if (is.na(repo_root)) {
  cwd_candidate <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (dir.exists(file.path(cwd_candidate, "workflow", "script"))) {
    repo_root <- cwd_candidate
  }
}

if (is.na(repo_root)) {
  stop("Could not resolve repository root.")
}

main_r <- file.path(repo_root, "workflow", "script", "tree_annotation", "R", "main.R")
if (!file.exists(main_r)) {
  stop(sprintf("main.R not found: %s", main_r))
}

suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
source(main_r, local = .GlobalEnv)

# 1) tidy_df_tip: group order and tip order are stable and numeric conversion works.
df_tip <- data.frame(
  label = c("g2", "g1"),
  y = c(2, 1),
  B = c("2", "4"),
  A = c("1", "3"),
  stringsAsFactors = FALSE
)
df_trait <- data.frame(B = c(2, 4), A = c(1, 3), row.names = c("g2", "g1"))
tidy <- tidy_df_tip(df_tip, df_trait)
if (nrow(tidy) != 4) stop("tidy_df_tip returned unexpected row count.")
if (!all(levels(tidy$group) == c("A", "B"))) stop("tidy_df_tip group order is incorrect.")

# 2) get_rel_widths: --rel_widths override string is parsed and applied.
dummy_tip_plot <- ggplot(data.frame(label = c("a", "bb"), x = c(0, 0), y = c(1, 2))) + geom_blank(aes(x = x, y = y))
dummy_tip_plot$data$label <- factor(dummy_tip_plot$data$label, levels = dummy_tip_plot$data$label)
g <- list(tree = ggplot(), tiplabel = dummy_tip_plot)
w <- get_rel_widths(g, "tree,2,tiplabel,1")
if (!("tree" %in% names(w) && "tiplabel" %in% names(w))) stop("get_rel_widths missing expected keys.")
if (abs(unname(w["tree"]) - 2) > 1e-9) stop("get_rel_widths did not apply tree override.")
if (abs(unname(w["tiplabel"]) - 1) > 1e-9) stop("get_rel_widths did not apply tiplabel override.")

# 3) get_df_trait: relative scaling should not produce Inf/NaN for all-zero rows.
b <- data.frame(
  so_event = c("L", "L"),
  node_name = c("n1", "n2"),
  expression_t1 = c(0, 0),
  expression_t2 = c(0, 5),
  stringsAsFactors = FALSE
)
trait <- get_df_trait(b, transform = "no", scale = "rel", trait_prefix = "expression_", negative2zero = TRUE)
if (any(is.infinite(as.matrix(trait)))) stop("get_df_trait generated Inf values.")
if (any(is.nan(as.matrix(trait)))) stop("get_df_trait generated NaN values.")

# 4) get_df_N: N-slice parser returns expected rows.
df_tip_n <- data.frame(
  promoter_N = c("1,3:5", "", NA),
  y = c(10, 20, 30),
  stringsAsFactors = FALSE
)
df_n <- get_df_N(df_tip_n)
if (nrow(df_n) != 2) stop("get_df_N returned unexpected number of rows.")
if (!all(colnames(df_n) == c("x_start", "x_end", "y_start", "y_end"))) stop("get_df_N column names are incorrect.")
df_tip_n_bad <- data.frame(
  promoter_N = c("bad,7:3,8:x,9"),
  y = c(5),
  stringsAsFactors = FALSE
)
df_n_bad <- get_df_N(df_tip_n_bad)
if (nrow(df_n_bad) != 2) stop("get_df_N should ignore malformed tokens and keep valid ones.")
if (!any(abs(df_n_bad$x_start - 2.5) < 1e-9 & abs(df_n_bad$x_end - 7.5) < 1e-9)) {
  stop("get_df_N should normalize reversed ranges (e.g., 7:3).")
}

# 5) add_gene_cluster_membership: splits by intergenic distance and chromosome.
df_tip_cluster <- data.frame(
  so_event = c("L", "L", "L", "L"),
  taxon = c("Sp one", "Sp one", "Sp one", "Sp one"),
  chromosome = c("chr1", "chr1", "chr1", "chr2"),
  label = c("g1", "g2", "g3", "g4"),
  start = c(100, 140, 1000, 50),
  end = c(120, 160, 1020, 70),
  stringsAsFactors = FALSE
)
clustered <- add_gene_cluster_membership(df_tip_cluster, max_bp_membership = 100)
cid <- setNames(as.character(clustered$cluster_membership), clustered$label)
if (!(cid[["g1"]] == cid[["g2"]])) stop("Nearby genes should share the same cluster_membership.")
if (cid[["g2"]] == cid[["g3"]]) stop("Distant genes should not share the same cluster_membership.")
if (cid[["g3"]] == cid[["g4"]]) stop("Different chromosomes should not share the same cluster_membership.")

# 6) add_complete_overlap_groups: fully overlapping motifs are merged.
df_fimo_overlap <- data.frame(
  label = c("geneA", "geneA", "geneB", "geneB"),
  start = c(10, 12, 20, 22),
  end = c(30, 28, 40, 38),
  motif_altid = c("M1", "M2", "M1", "M2"),
  stringsAsFactors = FALSE
)
merged <- add_complete_overlap_groups(df_fimo_overlap, "motif_altid")
merged_ids <- unique(as.character(merged$motif_altid))
if (length(merged_ids) != 1) stop("Completely overlapping motifs should be merged into one ID.")

# 7) add_pointplot_column: multi-char replicate separator is handled correctly.
g_tree_min <- list(
  data = data.frame(
    isTip = c(TRUE, TRUE),
    label = c("g1", "g2"),
    y = c(1, 2),
    tiplab_color = c("black", "black"),
    stringsAsFactors = FALSE
  )
)
g_in <- list(tree = g_tree_min)
args_pp <- list(
  font_size = 6,
  trait_axis_label = "Expression",
  trait_point_size = 1,
  trait_colors = c("#000000"),
  margins = c(0.01, 0, 0.01, 0)
)
df_trait_pp <- data.frame(
  "A__rep1" = c(1, 2),
  "A__rep2" = c(3, 4),
  row.names = c("g1", "g2")
)
g_pp <- add_pointplot_column(g_in, args_pp, df_trait_pp, replicate_sep = "__")
pp_groups <- unique(as.character(g_pp$pointplot$data$group))
if (!(length(pp_groups) == 1 && pp_groups[1] == "A")) stop("add_pointplot_column failed to strip multi-char replicate suffix.")

# 8) add_fimo_column: NA fimo_start should be skipped safely.
g_fimo_in <- list(
  tree = list(
    data = data.frame(
      isTip = c(TRUE),
      label = c("g1"),
      y = c(1),
      fimo_start = c(NA),
      fimo_end = c(NA),
      fimo_strand = c(NA),
      fimo_alt_id = c(NA),
      fimo_qvalue = c(NA),
      promoter_available = c("Y"),
      promoter_N = c(""),
      stringsAsFactors = FALSE
    )
  )
)
g_fimo_out <- add_fimo_column(g_fimo_in, list(font_size = 6, margins = c(0, 0, 0, 0)), qname = "fimo", xmax = 100, qvalue = 0.01)
if ("fimo" %in% names(g_fimo_out)) stop("add_fimo_column should not add panel when all FIMO entries are NA.")

# 8b) add_fimo_column: entries failing qvalue cutoff should not add panel.
g_fimo_cutoff_in <- list(
  tree = list(
    data = data.frame(
      isTip = c(TRUE),
      label = c("g1"),
      y = c(1),
      fimo_start = c("10"),
      fimo_end = c("14"),
      fimo_strand = c("+"),
      fimo_alt_id = c("M1"),
      fimo_qvalue = c("0.5"),
      promoter_available = c("Y"),
      promoter_N = c(""),
      stringsAsFactors = FALSE
    )
  )
)
g_fimo_cutoff_out <- add_fimo_column(g_fimo_cutoff_in, list(font_size = 6, margins = c(0, 0, 0, 0)), qname = "fimo", xmax = 100, qvalue = 0.01)
if ("fimo" %in% names(g_fimo_cutoff_out)) stop("add_fimo_column should not add panel when all FIMO entries fail qvalue cutoff.")

# 8c) add_fimo_column: when polygon data are empty, panel should still be created safely.
g_fimo_poly_empty_in <- list(
  tree = list(
    data = data.frame(
      isTip = c(TRUE),
      label = c("g1"),
      y = c(1),
      fimo_start = c("10"),
      fimo_end = c("14"),
      fimo_strand = c("+"),
      fimo_alt_id = c("M1"),
      fimo_qvalue = c("0.0001"),
      promoter_available = c("Y"),
      promoter_N = c(""),
      stringsAsFactors = FALSE
    )
  )
)
g_fimo_poly_empty_out <- add_fimo_column(g_fimo_poly_empty_in, list(font_size = 6, margins = c(0, 0, 0, 0)), qname = "fimo", xmax = 100, qvalue = 0.01)
if (!("fimo" %in% names(g_fimo_poly_empty_out))) stop("add_fimo_column should create panel even when polygon data are empty.")

# 8d) add_fimo_column: ncpu path should be equivalent to serial output.
g_fimo_parallel_in <- list(
  tree = list(
    data = data.frame(
      isTip = c(TRUE, TRUE),
      label = c("g1", "g2"),
      y = c(1, 2),
      fimo_start = c("10;20", "12"),
      fimo_end = c("14;24", "16"),
      fimo_strand = c("+;+", "-"),
      fimo_alt_id = c("M1;M2", "M1"),
      fimo_qvalue = c("0.0001;0.0002", "0.0003"),
      promoter_available = c("Y", "Y"),
      promoter_N = c("", ""),
      stringsAsFactors = FALSE
    )
  )
)
fimo_args <- list(font_size = 6, margins = c(0, 0, 0, 0))
g_fimo_serial <- add_fimo_column(g_fimo_parallel_in, fimo_args, qname = "fimo", xmax = 100, qvalue = 0.01, ncpu = 1)
g_fimo_parallel <- add_fimo_column(g_fimo_parallel_in, fimo_args, qname = "fimo", xmax = 100, qvalue = 0.01, ncpu = 2)
if (!("fimo" %in% names(g_fimo_serial))) stop("add_fimo_column serial path did not create a panel.")
if (!("fimo" %in% names(g_fimo_parallel))) stop("add_fimo_column parallel path did not create a panel.")
serial_df <- g_fimo_serial$fimo$data[, c("label", "start", "end", "strand", "motif_altid", "count")]
parallel_df <- g_fimo_parallel$fimo$data[, c("label", "start", "end", "strand", "motif_altid", "count")]
serial_df$motif_altid <- as.character(serial_df$motif_altid)
parallel_df$motif_altid <- as.character(parallel_df$motif_altid)
serial_df <- serial_df[order(serial_df$label, serial_df$start, serial_df$end, serial_df$motif_altid), ]
parallel_df <- parallel_df[order(parallel_df$label, parallel_df$start, parallel_df$end, parallel_df$motif_altid), ]
rownames(serial_df) <- NULL
rownames(parallel_df) <- NULL
if (!isTRUE(all.equal(serial_df, parallel_df, check.attributes = FALSE))) {
  stop("add_fimo_column parallel path output differs from serial output.")
}

# 8e) add_fimo_column: semicolon-delimited tokens with spaces should be trimmed safely.
g_fimo_spaces_in <- list(
  tree = list(
    data = data.frame(
      isTip = c(TRUE),
      label = c("g1"),
      y = c(1),
      fimo_start = c("10; 20"),
      fimo_end = c("14; 24"),
      fimo_strand = c("+; -"),
      fimo_alt_id = c("M1; M2"),
      fimo_qvalue = c("0.0001; 0.2"),
      promoter_available = c("Y"),
      promoter_N = c(""),
      stringsAsFactors = FALSE
    )
  )
)
g_fimo_spaces_out <- add_fimo_column(g_fimo_spaces_in, fimo_args, qname = "fimo", xmax = 100, qvalue = 0.01, ncpu = 1)
if (!("fimo" %in% names(g_fimo_spaces_out))) stop("add_fimo_column should parse spaced semicolon-delimited tokens.")
motif_vals_spaces <- as.character(g_fimo_spaces_out$fimo$data$motif_altid)
if (!(length(motif_vals_spaces) == 1 && !grepl("\\s", motif_vals_spaces[1]))) {
  stop("add_fimo_column should trim spaces and apply qvalue cutoff correctly.")
}

# 9) add_pointplot_column: all-NA values should return without adding panel.
df_trait_na <- data.frame(
  "A__rep1" = c(NA_real_, NA_real_),
  "A__rep2" = c(NA_real_, NA_real_),
  row.names = c("g1", "g2")
)
g_pp_na <- add_pointplot_column(g_in, args_pp, df_trait_na, replicate_sep = "__")
if ("pointplot" %in% names(g_pp_na)) stop("add_pointplot_column should skip panel when all values are NA.")

# 10) add_heatmap_column: all-NA values should return without adding panel.
args_hm <- list(
  font_size = 6,
  margins = c(0.01, 0, 0.01, 0),
  panel1 = "tree,bl_rooted,support_unrooted,species,L"
)
g_hm_na <- add_heatmap_column(g_in, args_hm, df_trait_na)
if ("heatmap" %in% names(g_hm_na)) stop("add_heatmap_column should skip panel when all values are NA.")

# 11) enhance_branch_table: one subroot support should be masked.
b_min <- data.frame(
  branch_id = c(3, 1, 2),
  parent = c(-999, 3, 3),
  child1 = c(1, NA, NA),
  child2 = c(2, NA, NA),
  so_event_parent = c("S", "D", "S"),
  support_unrooted = c(100, 95, 90),
  node_name = c("n3", "n1", "n2"),
  stringsAsFactors = FALSE
)
enhanced <- enhance_branch_table(
  b_min,
  args = list(max_delta_intron_present = -0.5),
  event_method = "species_overlap"
)
subroot_support <- enhanced[enhanced$parent == 3, "support_unrooted"]
if (sum(is.na(subroot_support)) != 1) stop("enhance_branch_table should mask support for one subroot branch.")

# 12) get_df_protein_backbone: invalid coordinates are ignored safely.
df_rps_backbone <- data.frame(
  qacc = c("q1", "q1", "q2"),
  qlen = c("10", "10", "bad"),
  qstart = c(2, -5, 1),
  qend = c(4, 100, 3),
  stringsAsFactors = FALSE
)
pb_seg <- get_df_protein_backbone(df_rps_backbone, mode = "segmented")
if (!all(c("x", "xend", "label") %in% colnames(pb_seg))) stop("get_df_protein_backbone(segmented) returned unexpected columns.")
pb_one <- get_df_protein_backbone(df_rps_backbone, mode = "one_line")
if (!all(pb_one$label %in% c("q1"))) stop("get_df_protein_backbone(one_line) should skip invalid qlen rows.")
df_rps_backbone_tail <- data.frame(
  qacc = c("q1", "q1"),
  qlen = c("10", "10"),
  qstart = c(2, 8),
  qend = c(4, 9),
  stringsAsFactors = FALSE
)
pb_tail <- get_df_protein_backbone(df_rps_backbone_tail, mode = "segmented")
if (!(nrow(pb_tail) == 3 && any(pb_tail$x == 9 & pb_tail$xend == 10))) {
  stop("get_df_protein_backbone(segmented) should keep trailing non-domain segment.")
}

# 13) get_df_domain: malformed rows should not crash and should return valid schema.
df_rps_domain <- data.frame(
  qacc = c("q1", "q1", "q2"),
  qlen = c("10", "10", "0"),
  qstart = c(2, "bad", 1),
  qend = c(4, 6, 3),
  ymin = c(1, 1, 2),
  ymax = c(2, 2, 3),
  sacc = c("D1", "D2", "D3"),
  stringsAsFactors = FALSE
)
dom <- get_df_domain(df_rps_domain)
if (!all(c("xmin", "xmax", "ymin", "ymax", "sacc", "label") %in% colnames(dom))) stop("get_df_domain returned unexpected columns.")

# 14) get_df_clade + split_nested_clade: grouping and nested split are stable.
df_nearest <- data.frame(
  query = c("A_1", "A_2", "A_3", "A_4"),
  nearests = c("SPX_g1", "SPX_g1", "SPX_g2", "SPX_g2"),
  y = c(1, 2, 3, 4),
  stringsAsFactors = FALSE
)
clade <- get_df_clade(df_nearest, ortholog_prefix = "SPX_")
if (nrow(clade) != 2) stop("get_df_clade should create two grouped clades.")
df_nested <- data.frame(
  nearests = c("x", "y"),
  ymin = c(1, 2),
  ymid = c(5, 3),
  ymax = c(9, 4),
  is_self = c(FALSE, TRUE),
  num_line = c(8, 2),
  num_label = c(1, 1),
  stringsAsFactors = FALSE
)
split_out <- split_nested_clade(df_nested, margin = 0.3)
if (nrow(split_out) < 2) stop("split_nested_clade should preserve/split nested rows, not drop all.")

# 15) merge_overlap_group: contiguous overlaps are grouped, separated intervals split.
df_merge <- data.frame(
  label = c("g1", "g1", "g1"),
  motif_altid = c("M1", "M1", "M1"),
  start = c(1, 3, 20),
  end = c(4, 5, 25),
  strand = c("+", "+", "+"),
  y = c(1, 1, 1),
  stringsAsFactors = FALSE
)
mg <- merge_overlap_group(df_merge, merge_level = "TF")
if (nrow(mg) != 2) stop("merge_overlap_group should split non-overlapping motif blocks.")

# 16) get_df_polygon: should return polygons for simple two-gene motif chain.
df_poly_in <- data.frame(
  motif_altid = c("M1", "M1"),
  label = c("g1", "g2"),
  y = c(1, 2),
  ymin = c(1, 2),
  ymax = c(1.3, 2.3),
  start = c(10, 12),
  end = c(14, 16),
  stringsAsFactors = FALSE
)
poly <- get_df_polygon(df_poly_in, min_count = 1, max_count = Inf, multiple_connection = "align_from_tss")
if (nrow(poly) == 0) stop("get_df_polygon should emit polygon rows for connected motifs.")
if (!all(c("x", "y", "polygon_id", "motif_altid", "alpha_value") %in% colnames(poly))) {
  stop("get_df_polygon returned unexpected schema.")
}

# 17) get_df_polygon: empty after filtering should return empty schema safely.
poly_empty <- get_df_polygon(df_poly_in, min_count = 999, max_count = Inf, multiple_connection = "align_from_tss")
if (nrow(poly_empty) != 0) stop("get_df_polygon should return empty output when nothing passes min_count.")
if (!all(c("x", "y", "polygon_id", "motif_altid", "alpha_value") %in% colnames(poly_empty))) {
  stop("get_df_polygon empty output schema is incorrect.")
}

# 18) merge_tip_trait: preserves row order and joins by label.
df_tip_mt <- data.frame(label = c("g2", "g1"), y = c(2, 1), stringsAsFactors = FALSE)
df_trait_mt <- data.frame(v = c(10, 20), row.names = c("g1", "g2"))
mt <- merge_tip_trait(df_tip_mt, df_trait_mt)
if (!all(mt$label == c("g2", "g1"))) stop("merge_tip_trait should preserve df_tip row order.")
if (!all(mt$v == c(20, 10))) stop("merge_tip_trait joined values are incorrect.")

# 19) get_line_coordinate: returns valid chains and pairwise lines.
g_line <- list(data = data.frame(
  branch_id = c(1, 2, 3),
  branch = c(0.1, 0.2, 0.3),
  y = c(3, 2, 1),
  stringsAsFactors = FALSE
))
lc_chain <- get_line_coordinate(g_line, c(1, 2, 3), jitter = FALSE, pairwise = FALSE)
if (nrow(lc_chain) != 2) stop("get_line_coordinate(pairwise=FALSE) should return n-1 segments.")
lc_pair <- get_line_coordinate(g_line, c(1, 2, 3), pairwise = TRUE)
if (nrow(lc_pair) != 3) stop("get_line_coordinate(pairwise=TRUE) should return nC2 segments.")

# 20) propagate_tiplab_colors_to_internal_branches: paints internal branch when descendant tips share one color.
gtree_col <- list(data = data.frame(
  node = c(1, 2, 3),
  parent = c(3, 3, NA),
  isTip = c(TRUE, TRUE, FALSE),
  branch_color = c("red", "red", "black"),
  x = c(1, 1, 0),
  stringsAsFactors = FALSE
))
gtree_col2 <- propagate_tiplab_colors_to_internal_branches(gtree_col, tree = NULL)
internal_col <- gtree_col2$data$branch_color[gtree_col2$data$isTip == FALSE]
if (!(length(internal_col) == 1 && internal_col == "red")) stop("Internal branch color propagation failed.")

# 21) get_df_polygon: NA y rows should be ignored without crashing.
df_poly_na <- data.frame(
  motif_altid = c("M1", "M1", "M1"),
  label = c("g1", "g2", "g3"),
  y = c(1, NA, 3),
  ymin = c(1, NA, 3),
  ymax = c(1.3, NA, 3.3),
  start = c(10, 20, 30),
  end = c(14, 24, 34),
  stringsAsFactors = FALSE
)
poly_na <- get_df_polygon(df_poly_na, min_count = 1, max_count = Inf, multiple_connection = "align_from_tss")
if (!all(c("x", "y", "polygon_id", "motif_altid", "alpha_value") %in% colnames(poly_na))) {
  stop("get_df_polygon NA-y handling returned unexpected schema.")
}

# 22) merge_overlap_group: start/end reversal should be normalized safely.
df_merge_rev <- data.frame(
  label = c("g1", "g1"),
  motif_altid = c("M1", "M1"),
  start = c(10, 6),
  end = c(5, 12),
  strand = c("+", "-"),
  y = c(1, 1),
  stringsAsFactors = FALSE
)
mg_rev <- merge_overlap_group(df_merge_rev, merge_level = "TF")
if (nrow(mg_rev) != 1) stop("merge_overlap_group should merge overlapping reversed-coordinate intervals.")
if (!(mg_rev$start[1] == 5 && mg_rev$end[1] == 12)) {
  stop("merge_overlap_group did not normalize reversed start/end correctly.")
}

# 23) add_motif_colors: NA-only motif IDs should be handled without size mismatch errors.
df_color_na <- data.frame(motif_altid = NA_character_, stringsAsFactors = FALSE)
colored_na <- add_motif_colors(df_color_na, min_count = 2)
if (!all(c("count", "color_id") %in% colnames(colored_na))) stop("add_motif_colors should append count/color_id columns.")
if (!is.na(colored_na$count[1])) stop("add_motif_colors should keep NA count for NA motif IDs.")
df_overlap_naid <- data.frame(
  label = c("g1", "g2"),
  start = c(10, 20),
  end = c(15, 25),
  motif_altid = c(NA_character_, "M1"),
  stringsAsFactors = FALSE
)
overlap_naid <- add_complete_overlap_groups(df_overlap_naid, "motif_altid")
if (!is.na(overlap_naid$motif_altid[1])) stop("add_complete_overlap_groups should keep NA motif IDs as NA.")

# 24) merge_overlap_group: transitive overlap chain should stay in one group.
df_merge_chain <- data.frame(
  label = c("g1", "g1", "g1"),
  motif_altid = c("M1", "M1", "M1"),
  start = c(1, 2, 9),
  end = c(10, 3, 11),
  strand = c("+", "+", "+"),
  y = c(1, 1, 1),
  stringsAsFactors = FALSE
)
mg_chain <- merge_overlap_group(df_merge_chain, merge_level = "TF")
if (nrow(mg_chain) != 1) stop("merge_overlap_group should keep transitive-overlap intervals in one group.")

# 25) merge_overlap_group: reversed wide interval should bridge blocks after normalization.
df_merge_bridge <- data.frame(
  label = c("g1", "g1", "g1"),
  motif_altid = c("M1", "M1", "M1"),
  start = c(10, 12, 100),
  end = c(11, 13, 1),
  strand = c("+", "+", "+"),
  y = c(1, 1, 1),
  stringsAsFactors = FALSE
)
mg_bridge <- merge_overlap_group(df_merge_bridge, merge_level = "TF")
if (nrow(mg_bridge) != 1) stop("merge_overlap_group should merge intervals connected by a reversed wide interval.")

# 26) add_complete_overlap_groups: ncpu path should match serial output.
df_overlap_parallel <- data.frame(
  label = rep(c("geneA", "geneB", "geneC"), each = 4),
  start = c(10, 11, 10, 12, 20, 21, 20, 22, 30, 31, 30, 32),
  end = c(14, 15, 14, 16, 24, 25, 24, 26, 34, 35, 34, 36),
  motif_altid = rep(c("X1", "X2", "Y1", "Y2"), times = 3),
  stringsAsFactors = FALSE
)
overlap_serial <- add_complete_overlap_groups(df_overlap_parallel, "motif_altid", ncpu = 1)
overlap_parallel <- add_complete_overlap_groups(df_overlap_parallel, "motif_altid", ncpu = 2)
serial_cmp <- overlap_serial[, c("label", "start", "end", "motif_altid", "complete_overlap_group")]
parallel_cmp <- overlap_parallel[, c("label", "start", "end", "motif_altid", "complete_overlap_group")]
serial_cmp$motif_altid <- as.character(serial_cmp$motif_altid)
parallel_cmp$motif_altid <- as.character(parallel_cmp$motif_altid)
serial_cmp <- serial_cmp[order(serial_cmp$label, serial_cmp$start, serial_cmp$end, serial_cmp$motif_altid), ]
parallel_cmp <- parallel_cmp[order(parallel_cmp$label, parallel_cmp$start, parallel_cmp$end, parallel_cmp$motif_altid), ]
rownames(serial_cmp) <- NULL
rownames(parallel_cmp) <- NULL
if (!isTRUE(all.equal(serial_cmp, parallel_cmp, check.attributes = FALSE))) {
  stop("add_complete_overlap_groups parallel path output differs from serial output.")
}

# 27) get_df_polygon: NA motif_altid rows should be ignored safely.
df_poly_na_motif <- data.frame(
  motif_altid = c(NA_character_, NA_character_, "M1", "M1"),
  label = c("g1", "g2", "g1", "g2"),
  y = c(1, 2, 1, 2),
  ymin = c(1, 2, 1, 2),
  ymax = c(1.3, 2.3, 1.3, 2.3),
  start = c(10, 12, 20, 22),
  end = c(14, 16, 24, 26),
  stringsAsFactors = FALSE
)
poly_na_motif <- get_df_polygon(df_poly_na_motif, min_count = 1, max_count = Inf, multiple_connection = "align_from_tss")
if (any(is.na(as.character(poly_na_motif$motif_altid)))) {
  stop("get_df_polygon should skip NA motif_altid rows.")
}

# 28) get_df_polygon: ncpu path should match serial output (ignoring polygon_id labels).
df_poly_parallel <- data.frame(
  motif_altid = c("M1", "M1", "M2", "M2"),
  label = c("g1", "g2", "g1", "g2"),
  y = c(1, 2, 1, 2),
  ymin = c(1, 2, 1, 2),
  ymax = c(1.3, 2.3, 1.3, 2.3),
  start = c(10, 12, 30, 32),
  end = c(14, 16, 34, 36),
  stringsAsFactors = FALSE
)
poly_serial <- get_df_polygon(df_poly_parallel, min_count = 1, max_count = Inf, multiple_connection = "align_from_tss", ncpu = 1)
poly_parallel <- get_df_polygon(df_poly_parallel, min_count = 1, max_count = Inf, multiple_connection = "align_from_tss", ncpu = 2)
serial_cmp <- poly_serial[, c("x", "y", "motif_altid", "alpha_value")]
parallel_cmp <- poly_parallel[, c("x", "y", "motif_altid", "alpha_value")]
serial_cmp$motif_altid <- as.character(serial_cmp$motif_altid)
parallel_cmp$motif_altid <- as.character(parallel_cmp$motif_altid)
serial_cmp <- serial_cmp[order(serial_cmp$motif_altid, serial_cmp$y, serial_cmp$x, serial_cmp$alpha_value), ]
parallel_cmp <- parallel_cmp[order(parallel_cmp$motif_altid, parallel_cmp$y, parallel_cmp$x, parallel_cmp$alpha_value), ]
rownames(serial_cmp) <- NULL
rownames(parallel_cmp) <- NULL
if (!isTRUE(all.equal(serial_cmp, parallel_cmp, check.attributes = FALSE))) {
  stop("get_df_polygon parallel path output differs from serial output.")
}

cat("tree_annotation main.R tests passed.\n")
