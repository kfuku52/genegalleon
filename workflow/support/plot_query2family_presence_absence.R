#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ape, quietly = TRUE))
suppressPackageStartupMessages(library(cowplot, quietly = TRUE))
suppressPackageStartupMessages(library(ggplot2, quietly = TRUE))
suppressPackageStartupMessages(library(ggtree, quietly = TRUE))

options(stringsAsFactors = FALSE)

font_size_pt <- 8
font_size_mm <- font_size_pt * 0.352777778
support_color <- "#8c2d04"
hpd_color <- "#6b7280"
ortholog_colors <- c(undetected = "#f5f5f5", detected = "#2166ac")
busco_colors <- c(Single = "#000000", Duplicated = "#b22222", Fragmented = "#666666", Missing = "#cccccc")
heatmap_cell_pitch <- 0.68

parse_args <- function(argv) {
  out <- list()
  for (arg in argv) {
    if (!grepl("^--", arg)) {
      next
    }
    kv <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- kv[[1]]
    value <- if (length(kv) >= 2) paste(kv[-1], collapse = "=") else "1"
    out[[key]] <- value
  }
  out
}

get_arg <- function(args, key, default = "") {
  if (key %in% names(args) && nzchar(as.character(args[[key]]))) {
    return(as.character(args[[key]]))
  }
  default
}

has_nonempty_file <- function(file_path) {
  if (is.null(file_path) || length(file_path) == 0 || is.na(file_path) || !nzchar(file_path)) {
    return(FALSE)
  }
  if (!file.exists(file_path)) {
    return(FALSE)
  }
  file_info <- file.info(file_path)
  !is.na(file_info$size[[1]]) && file_info$size[[1]] > 0
}

descendant_tips <- function(tree, node) {
  n_tip <- length(tree[["tip.label"]])
  children <- tree[["edge"]][tree[["edge"]][, 1] == node, 2]
  tips <- children[children <= n_tip]
  internal_children <- children[children > n_tip]
  for (child in internal_children) {
    tips <- c(tips, descendant_tips(tree, child))
  }
  tips
}

descendant_tip_labels <- function(tree, node, tip_set = NULL) {
  labels <- tree[["tip.label"]][descendant_tips(tree, node)]
  if (!is.null(tip_set)) {
    labels <- labels[labels %in% tip_set]
  }
  labels
}

node_clade_key <- function(tree, node, tip_set = NULL) {
  tips <- sort(descendant_tip_labels(tree, node, tip_set = tip_set))
  if (length(tips) < 2) {
    return(NA_character_)
  }
  paste(tips, collapse = "\t")
}

node_split_key <- function(tree, node, tip_set = NULL) {
  if (is.null(tip_set)) {
    tip_set <- tree[["tip.label"]]
  }
  tip_set <- sort(unique(tip_set))
  tips <- sort(descendant_tip_labels(tree, node, tip_set = tip_set))
  if (length(tips) == 0 || length(tips) == length(tip_set)) {
    return(NA_character_)
  }
  other <- sort(setdiff(tip_set, tips))
  if (length(other) < length(tips)) {
    tips <- other
  } else if (length(other) == length(tips)) {
    tip_key <- paste(tips, collapse = "\t")
    other_key <- paste(other, collapse = "\t")
    if (other_key < tip_key) {
      tips <- other
    }
  }
  paste(tips, collapse = "\t")
}

transfer_internal_labels <- function(target_tree, source_tree) {
  target_tips <- target_tree[["tip.label"]]
  if (length(setdiff(target_tips, source_tree[["tip.label"]])) > 0) {
    warning("Cannot transfer support values: support tree does not contain all target tips.")
    return(target_tree)
  }

  source_nodes <- seq.int(
    length(source_tree[["tip.label"]]) + 1,
    length(source_tree[["tip.label"]]) + source_tree[["Nnode"]]
  )
  source_labels <- source_tree[["node.label"]]
  if (is.null(source_labels) || length(source_labels) == 0) {
    return(target_tree)
  }
  label_by_split <- list()
  for (idx in seq_along(source_nodes)) {
    if (idx > length(source_labels)) {
      next
    }
    label <- source_labels[[idx]]
    if (is.na(label) || !nzchar(label)) {
      next
    }
    split_key <- node_split_key(source_tree, source_nodes[[idx]], tip_set = target_tips)
    if (!is.na(split_key) && nzchar(split_key)) {
      label_by_split[[split_key]] <- label
    }
  }

  n_tip <- length(target_tree[["tip.label"]])
  target_nodes <- seq.int(n_tip + 1, n_tip + target_tree[["Nnode"]])
  target_labels <- rep("", target_tree[["Nnode"]])
  for (idx in seq_along(target_nodes)) {
    split_key <- node_split_key(target_tree, target_nodes[[idx]], tip_set = target_tips)
    if (!is.na(split_key) && nzchar(split_key) && !is.null(label_by_split[[split_key]])) {
      target_labels[[idx]] <- label_by_split[[split_key]]
    }
  }
  target_tree[["node.label"]] <- target_labels
  target_tree
}

format_support_labels <- function(labels) {
  values <- suppressWarnings(as.numeric(labels))
  values[!is.finite(values)] <- NA_real_
  if (all(is.na(values))) {
    return(rep("", length(labels)))
  }
  max_value <- max(values, na.rm = TRUE)
  if (is.finite(max_value) && max_value <= 1.0001) {
    return(ifelse(is.na(values), "", sprintf("%.2f", values)))
  }
  ifelse(is.na(values), "", sprintf("%.0f", round(values, digits = 0)))
}

extract_numeric_ci_from_label <- function(label_text) {
  if (is.na(label_text) || nchar(trimws(label_text)) == 0) {
    return(NULL)
  }
  nums <- regmatches(label_text, gregexpr("[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?", label_text, perl = TRUE))[[1]]
  nums <- suppressWarnings(as.numeric(nums))
  nums <- nums[is.finite(nums)]
  if (length(nums) < 2) {
    return(NULL)
  }
  c(min(nums[1], nums[2]), max(nums[1], nums[2]))
}

extract_ci_table <- function(ci_tree_text) {
  ci_pattern <- "\\[&95%(?:HPD)?=\\{[[:space:]]*[-+0-9.eE]+[[:space:]]*,[[:space:]]*[-+0-9.eE]+[[:space:]]*\\}\\]"
  ci_hits <- gregexpr(ci_pattern, ci_tree_text, perl = TRUE)
  hit_texts <- regmatches(ci_tree_text, ci_hits)[[1]]
  ci_tree_tagged <- ci_tree_text

  ci_table <- data.frame(
    tag = character(0),
    lower = numeric(0),
    upper = numeric(0),
    stringsAsFactors = FALSE
  )

  if (length(hit_texts) == 0) {
    cleaned <- gsub("\\[&[^]]*\\]", "", ci_tree_text, perl = TRUE)
    tree_ci <- ape::read.tree(text = cleaned)
  } else {
    for (i in seq_along(hit_texts)) {
      h <- hit_texts[i]
      mm <- regexec("\\{[[:space:]]*([-+0-9.eE]+)[[:space:]]*,[[:space:]]*([-+0-9.eE]+)[[:space:]]*\\}", h, perl = TRUE)
      cap <- regmatches(h, mm)[[1]]
      if (length(cap) < 3) {
        next
      }
      nums <- as.numeric(cap[2:3])
      if (any(is.na(nums))) {
        next
      }
      tag <- paste0("CI_TAG_", i)
      ci_tree_tagged <- sub(ci_pattern, tag, ci_tree_tagged, perl = TRUE)
      ci_table <- rbind(
        ci_table,
        data.frame(
          tag = tag,
          lower = min(nums[1], nums[2]),
          upper = max(nums[1], nums[2]),
          stringsAsFactors = FALSE
        )
      )
    }
    ci_tree_tagged <- gsub("\\[&[^]]*\\]", "", ci_tree_tagged, perl = TRUE)
    tree_ci <- ape::read.tree(text = ci_tree_tagged)
  }

  node_labels <- tree_ci[["node.label"]]
  if (!is.null(node_labels) && length(node_labels) > 0) {
    n_tip <- length(tree_ci[["tip.label"]])
    for (i in seq_along(node_labels)) {
      ci_vals <- extract_numeric_ci_from_label(node_labels[i])
      if (is.null(ci_vals)) {
        next
      }
      ci_table <- rbind(
        ci_table,
        data.frame(
          tag = paste0("NODE_LABEL_", i),
          lower = ci_vals[1],
          upper = ci_vals[2],
          node_ci = n_tip + i,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  if (nrow(ci_table) == 0) {
    return(list(tree_ci = tree_ci, ci_table = data.frame()))
  }

  if (!"node_ci" %in% colnames(ci_table)) {
    ci_table$node_ci <- NA_integer_
  }
  idx_missing_node <- is.na(ci_table$node_ci)
  if (any(idx_missing_node)) {
    n_tip <- length(tree_ci[["tip.label"]])
    node_labels <- tree_ci[["node.label"]]
    ci_table$node_ci[idx_missing_node] <- n_tip + match(ci_table$tag[idx_missing_node], node_labels)
  }
  ci_table <- ci_table[!is.na(ci_table$node_ci), c("tag", "lower", "upper", "node_ci"), drop = FALSE]
  list(tree_ci = tree_ci, ci_table = ci_table)
}

map_ci_nodes_to_mean_tree <- function(tree_ci, tree_mean, ci_table) {
  if (nrow(ci_table) == 0) {
    return(ci_table)
  }
  target_tips <- tree_mean[["tip.label"]]
  if (length(setdiff(target_tips, tree_ci[["tip.label"]])) > 0) {
    warning("Cannot transfer time CI values: CI tree does not contain all target tips.")
    return(data.frame())
  }

  n_tip_mean <- length(tree_mean[["tip.label"]])
  mean_nodes <- seq.int(n_tip_mean + 1, n_tip_mean + tree_mean[["Nnode"]])
  mean_signatures <- vapply(mean_nodes, function(x) node_clade_key(tree_mean, x, tip_set = target_tips), character(1))
  sig_to_node_mean <- stats::setNames(mean_nodes[!is.na(mean_signatures)], mean_signatures[!is.na(mean_signatures)])

  ci_table$signature <- vapply(ci_table$node_ci, function(x) node_clade_key(tree_ci, x, tip_set = target_tips), character(1))
  ci_table <- ci_table[!is.na(ci_table$signature) & ci_table$signature %in% names(sig_to_node_mean), , drop = FALSE]
  if (nrow(ci_table) == 0) {
    return(ci_table)
  }
  ci_table$node_mean <- unname(sig_to_node_mean[ci_table$signature])
  ci_table
}

read_time_ci_table <- function(tree_ci_path, tree_mean) {
  if (!has_nonempty_file(tree_ci_path)) {
    return(data.frame())
  }
  ci_tree_text <- paste(readLines(tree_ci_path, warn = FALSE), collapse = "")
  ci_parsed <- extract_ci_table(ci_tree_text)
  ci_table <- map_ci_nodes_to_mean_tree(ci_parsed$tree_ci, tree_mean, ci_parsed$ci_table)
  if (nrow(ci_table) == 0) {
    warning("No mappable 95% HPD annotations were found in: ", tree_ci_path)
  }
  ci_table
}

nice_scale_width <- function(span) {
  if (!is.finite(span) || span <= 0) {
    return(1)
  }
  raw_width <- span / 5
  base <- 10^floor(log10(raw_width))
  candidates <- c(1, 2, 5, 10) * base
  valid <- candidates[candidates <= raw_width]
  if (length(valid) == 0) {
    return(min(candidates))
  }
  max(valid)
}

format_ruler_labels <- function(x) {
  x <- as.numeric(x)
  max_value <- max(abs(x), na.rm = TRUE)
  digits <- if (is.finite(max_value) && max_value < 10) 1 else 0
  formatC(x, format = "f", digits = digits)
}

normalize_busco_species_id <- function(x) {
  x <- basename(x)
  x <- sub("\\.busco\\.full\\.tsv$", "", x)
  x <- sub("_busco\\.full\\.tsv$", "", x)
  x <- sub("\\.busco\\.short\\.txt$", "", x)
  x <- sub("_busco\\.short\\.txt$", "", x)
  x <- sub("\\.busco.*$", "", x)
  x <- sub("_busco.*$", "", x)
  x <- sub("([._-]?longestCDS)$", "", x)
  x <- sub("_sp_[^._-]+$", "_sp", x)
  x
}

match_busco_species <- function(species, tip_y_by_species) {
  species <- normalize_busco_species_id(species)
  candidates <- unique(c(
    species,
    gsub(" ", "_", species, fixed = TRUE),
    gsub("_", " ", species, fixed = TRUE)
  ))
  matched <- candidates[candidates %in% names(tip_y_by_species)]
  if (length(matched) > 0) {
    return(matched[[1]])
  }
  species
}

empty_busco_plot_table <- function() {
  data.frame(
    species = character(0),
    status = character(0),
    count = numeric(0),
    total = numeric(0),
    percent = numeric(0),
    y = numeric(0),
    stringsAsFactors = FALSE
  )
}

normalize_busco_status <- function(status) {
  status <- as.character(status)
  status[status == "Complete"] <- "Single"
  status
}

busco_counts_to_plot_rows <- function(counts, tip_y_by_species) {
  empty <- empty_busco_plot_table()
  if (nrow(counts) == 0) {
    return(empty)
  }
  counts$species <- vapply(as.character(counts$species), match_busco_species, character(1), tip_y_by_species = tip_y_by_species)
  counts$status <- normalize_busco_status(counts$status)
  counts$count <- suppressWarnings(as.numeric(counts$count))
  counts$count[!is.finite(counts$count)] <- 0
  counts <- counts[
    counts$species %in% names(tip_y_by_species) &
      counts$status %in% names(busco_colors),
    ,
    drop = FALSE
  ]
  if (nrow(counts) == 0) {
    return(empty)
  }

  species_levels <- names(tip_y_by_species)[names(tip_y_by_species) %in% unique(counts$species)]
  out <- list()
  out_idx <- 0L
  for (species in species_levels) {
    species_counts <- stats::setNames(rep(0, length(busco_colors)), names(busco_colors))
    tmp <- counts[counts$species == species, , drop = FALSE]
    if (nrow(tmp) == 0) {
      next
    }
    aggregate_counts <- stats::aggregate(count ~ status, data = tmp, FUN = sum)
    species_counts[aggregate_counts$status] <- aggregate_counts$count
    total <- sum(species_counts)
    if (!is.finite(total) || total <= 0) {
      next
    }
    x_start <- 0
    for (state in names(busco_colors)) {
      count <- species_counts[[state]]
      percent <- count / total * 100
      out_idx <- out_idx + 1L
      out[[out_idx]] <- data.frame(
        species = species,
        status = state,
        count = count,
        total = total,
        percent = percent,
        xmin_percent = x_start,
        xmax_percent = x_start + percent,
        y = unname(tip_y_by_species[[species]]),
        stringsAsFactors = FALSE
      )
      x_start <- x_start + percent
    }
  }
  if (length(out) == 0) {
    return(empty)
  }
  do.call(rbind, out)
}

read_busco_full_table_dir <- function(busco_dir_path, tip_y_by_species) {
  empty <- empty_busco_plot_table()
  files <- list.files(busco_dir_path, full.names = TRUE)
  if (length(files) == 0) {
    return(empty)
  }
  file_is_dir <- file.info(files)$isdir
  files <- files[!is.na(file_is_dir) & !file_is_dir]
  records <- list()
  record_idx <- 0L
  for (file_path in files) {
    species <- match_busco_species(basename(file_path), tip_y_by_species)
    if (!species %in% names(tip_y_by_species)) {
      next
    }
    all_lines <- readLines(file_path, warn = FALSE)
    if (!any(grepl("^# Busco id", all_lines))) {
      next
    }
    data_lines <- grep("^[^#]", all_lines, value = TRUE)
    if (length(data_lines) == 0) {
      next
    }
    split_lines <- strsplit(data_lines, "\t", fixed = TRUE)
    tmp <- data.frame(
      busco_id = vapply(split_lines, function(x) if (length(x) >= 1) x[[1]] else NA_character_, character(1)),
      status = vapply(split_lines, function(x) if (length(x) >= 2) x[[2]] else NA_character_, character(1)),
      species = species,
      stringsAsFactors = FALSE
    )
    tmp <- tmp[
      !is.na(tmp$busco_id) & nzchar(tmp$busco_id) &
        !is.na(tmp$status) & nzchar(tmp$status),
      ,
      drop = FALSE
    ]
    tmp$status <- normalize_busco_status(tmp$status)
    tmp <- tmp[tmp$status %in% names(busco_colors), , drop = FALSE]
    if (nrow(tmp) == 0) {
      next
    }
    record_idx <- record_idx + 1L
    records[[record_idx]] <- tmp[, c("species", "status", "busco_id"), drop = FALSE]
  }
  if (length(records) == 0) {
    return(empty)
  }
  df <- unique(do.call(rbind, records))
  counts <- stats::aggregate(busco_id ~ species + status, data = df, FUN = length)
  colnames(counts) <- c("species", "status", "count")
  busco_counts_to_plot_rows(counts, tip_y_by_species)
}

read_busco_count_summary_table <- function(busco, tip_y_by_species) {
  empty <- empty_busco_plot_table()
  if (nrow(busco) == 0 || ncol(busco) < 2) {
    return(empty)
  }
  normalized_cols <- tolower(gsub("[^a-z0-9]+", "_", colnames(busco)))
  species_col <- match("species", normalized_cols)
  if (is.na(species_col)) {
    return(empty)
  }
  status_cols <- vapply(names(busco_colors), function(status) {
    matches <- grep(paste0("(^|_)", tolower(status), "($|_)"), normalized_cols)
    matches <- setdiff(matches, species_col)
    if (length(matches) == 0) {
      return(NA_integer_)
    }
    matches[[1]]
  }, integer(1))
  if (all(is.na(status_cols))) {
    return(empty)
  }

  rows <- list()
  row_idx <- 0L
  for (i in seq_len(nrow(busco))) {
    species <- match_busco_species(as.character(busco[[species_col]][[i]]), tip_y_by_species)
    if (!species %in% names(tip_y_by_species)) {
      next
    }
    for (status in names(busco_colors)) {
      col_idx <- status_cols[[status]]
      count <- 0
      if (!is.na(col_idx)) {
        count <- suppressWarnings(as.numeric(busco[[col_idx]][[i]]))
      }
      if (!is.finite(count)) {
        count <- 0
      }
      row_idx <- row_idx + 1L
      rows[[row_idx]] <- data.frame(
        species = species,
        status = status,
        count = count,
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) {
    return(empty)
  }
  busco_counts_to_plot_rows(do.call(rbind, rows), tip_y_by_species)
}

read_busco_wide_gene_table <- function(busco, busco_table_path, tip_y_by_species) {
  empty <- empty_busco_plot_table()
  if (ncol(busco) < 4 || nrow(busco) == 0) {
    warning("BUSCO summary table does not contain species columns: ", busco_table_path)
    return(empty)
  }
  species_cols <- colnames(busco)[4:ncol(busco)]
  rows <- list()
  row_idx <- 0L
  for (col in species_cols) {
    species <- match_busco_species(col, tip_y_by_species)
    if (!species %in% names(tip_y_by_species)) {
      next
    }
    values <- as.character(busco[[col]])
    values[is.na(values)] <- "-"
    status <- ifelse(values == "-" | values == "", "Missing", ifelse(grepl(",", values, fixed = TRUE), "Duplicated", "Single"))
    tab <- table(factor(status, levels = names(busco_colors)))
    for (state in names(busco_colors)) {
      row_idx <- row_idx + 1L
      rows[[row_idx]] <- data.frame(
        species = species,
        status = state,
        count = as.numeric(tab[[state]]),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) {
    return(empty)
  }
  cat("BUSCO wide gene table does not encode Fragmented; plotting Fragmented as 0:", busco_table_path, "\n")
  busco_counts_to_plot_rows(do.call(rbind, rows), tip_y_by_species)
}

read_busco_summary_table <- function(busco_table_path, tip_y_by_species) {
  if (is.null(busco_table_path) || length(busco_table_path) == 0 || is.na(busco_table_path) || !nzchar(busco_table_path)) {
    return(empty_busco_plot_table())
  }
  if (dir.exists(busco_table_path)) {
    return(read_busco_full_table_dir(busco_table_path, tip_y_by_species))
  }
  if (!has_nonempty_file(busco_table_path)) {
    return(empty_busco_plot_table())
  }

  busco <- read.table(busco_table_path, sep = "\t", header = TRUE, quote = "", comment.char = "", check.names = FALSE)
  count_summary <- read_busco_count_summary_table(busco, tip_y_by_species)
  if (nrow(count_summary) > 0) {
    return(count_summary)
  }
  read_busco_wide_gene_table(busco, busco_table_path, tip_y_by_species)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
tree_path <- get_arg(args, "species_tree")
tree_ci_path <- get_arg(args, "species_tree_ci")
support_tree_path <- get_arg(args, "support_tree")
busco_table_path <- get_arg(args, "busco_table")
long_path <- get_arg(args, "long_table")
out_pdf <- get_arg(args, "out_pdf")
out_svg <- get_arg(args, "out_svg")
value_mode <- get_arg(args, "value", "presence")
plot_width <- as.numeric(get_arg(args, "width", "7.2"))
plot_height_arg <- get_arg(args, "height", "auto")

if (!nzchar(tree_path) || !file.exists(tree_path)) {
  stop("Missing --species_tree: ", tree_path)
}
if (!nzchar(long_path) || !file.exists(long_path)) {
  stop("Missing --long_table: ", long_path)
}
if (!value_mode %in% c("presence", "copy_number")) {
  stop("--value must be presence or copy_number: ", value_mode)
}
if (!is.finite(plot_width) || plot_width <= 0) {
  stop("Invalid --width: ", plot_width)
}
plot_width <- min(plot_width, 7.2)

tree <- ape::read.tree(tree_path)
tree_for_plot <- tree
if (has_nonempty_file(support_tree_path)) {
  tree_for_plot <- transfer_internal_labels(tree_for_plot, ape::read.tree(support_tree_path))
}

df <- read.table(long_path, sep = "\t", header = TRUE, quote = "", comment.char = "", check.names = FALSE)
if (nrow(df) == 0) {
  stop("Long table is empty: ", long_path)
}

tree_plot <- ggtree::ggtree(tree_for_plot, size = 0.3)
plot_data <- tree_plot$data
plot_data$y <- (plot_data$y - 1) * heatmap_cell_pitch + 1
tip_df <- subset(plot_data, isTip)
tip_levels <- as.character(tip_df$label[order(tip_df$y, decreasing = FALSE)])
tip_levels <- tip_levels[tip_levels %in% unique(as.character(df$species))]
if (length(tip_levels) == 0) {
  stop("No species overlap between tree tips and long table.")
}
tip_y_by_species <- stats::setNames(tip_df$y, as.character(tip_df$label))
y_values <- unname(tip_y_by_species[tip_levels])
y_limits <- c(min(y_values, na.rm = TRUE) - 0.55, max(y_values, na.rm = TRUE) + 0.55)

query_levels <- unique(as.character(df$query[order(df$query_order, df$query)]))
df <- df[df$species %in% tip_levels & df$query %in% query_levels, , drop = FALSE]
df$species <- factor(as.character(df$species), levels = tip_levels)
df$query <- factor(as.character(df$query), levels = query_levels)
df$y <- unname(tip_y_by_species[as.character(df$species)])
query_index_by_query <- stats::setNames(seq_along(query_levels), query_levels)
df$query_index <- unname(query_index_by_query[as.character(df$query)])

tree_height <- max(plot_data$x[plot_data$isTip], na.rm = TRUE)
tree_span <- max(plot_data$x, na.rm = TRUE) - min(plot_data$x, na.rm = TRUE)
if (!is.finite(tree_span) || tree_span <= 0) {
  tree_span <- tree_height
}
if (!is.finite(tree_height) || tree_height <= 0) {
  tree_height <- max(plot_data$x, na.rm = TRUE)
}
if (!is.finite(tree_height) || tree_height <= 0) {
  tree_height <- 1
}

ci_table <- read_time_ci_table(tree_ci_path, tree)
ci_plot_df <- data.frame()
if (nrow(ci_table) > 0) {
  node_xy <- plot_data[, c("node", "x", "y"), drop = FALSE]
  ci_plot_df <- merge(ci_table, node_xy, by.x = "node_mean", by.y = "node", all.x = FALSE, all.y = FALSE)
  if (nrow(ci_plot_df) > 0) {
    ci_plot_df$x_old <- tree_height - ci_plot_df$upper
    ci_plot_df$x_young <- tree_height - ci_plot_df$lower
    ci_plot_df$xmin <- pmin(ci_plot_df$x_old, ci_plot_df$x_young)
    ci_plot_df$xmax <- pmax(ci_plot_df$x_old, ci_plot_df$x_young)
  }
}

support_labels <- format_support_labels(plot_data$label)
root_node <- plot_data$node[plot_data$parent == plot_data$node]
is_root <- plot_data$node %in% root_node
is_subroot <- rep(FALSE, nrow(plot_data))
if (length(root_node) > 0) {
  is_subroot <- plot_data$parent %in% root_node & !plot_data$isTip
}
support_df <- plot_data[!plot_data$isTip & !is_root & !is_subroot & nzchar(support_labels), , drop = FALSE]
support_df$support_label <- support_labels[!plot_data$isTip & !is_root & !is_subroot & nzchar(support_labels)]

x_left <- 0
if (nrow(ci_plot_df) > 0) {
  x_left <- min(0, min(ci_plot_df$xmin, na.rm = TRUE))
}
x_right <- max(plot_data$x, na.rm = TRUE)
x_pad <- max((x_right - x_left) * 0.025, tree_height * 0.01, 0.001)
x_left <- x_left - x_pad
x_right <- x_right + x_pad

max_ruler_age <- tree_height
if (nrow(ci_plot_df) > 0 && any(is.finite(ci_plot_df$upper))) {
  max_ruler_age <- max(max_ruler_age, max(ci_plot_df$upper, na.rm = TRUE))
}
ruler_ages <- pretty(c(0, max_ruler_age), n = 5)
ruler_ages <- ruler_ages[ruler_ages >= 0]
ruler_breaks <- tree_height - ruler_ages
ruler_keep <- is.finite(ruler_breaks) & ruler_breaks >= x_left & ruler_breaks <= x_right
ruler_breaks <- ruler_breaks[ruler_keep]
ruler_labels <- format_ruler_labels(ruler_ages[ruler_keep])
ruler_order <- order(ruler_breaks)
ruler_breaks <- ruler_breaks[ruler_order]
ruler_labels <- ruler_labels[ruler_order]

label_df <- unique(df[, c("species", "species_display"), drop = FALSE])
label_df$species <- factor(as.character(label_df$species), levels = tip_levels)
label_df$y <- unname(tip_y_by_species[as.character(label_df$species)])
label_df <- label_df[order(label_df$y), , drop = FALSE]

if (value_mode == "copy_number") {
  df$plot_value <- suppressWarnings(as.numeric(df$copy_number))
  finite_values <- df$plot_value[is.finite(df$plot_value)]
  max_value <- if (length(finite_values) > 0) max(finite_values) else 1
  if (!is.finite(max_value) || max_value <= 0) {
    max_value <- 1
  }
  df$tile_fill <- ifelse(is.na(df$plot_value), "#d9d9d9", grDevices::colorRampPalette(c("#f5f5f5", "#2166ac"))(101)[pmax(1, pmin(101, round(df$plot_value / max_value * 100) + 1))])
} else {
  df$plot_value <- ifelse(is.na(df$presence), NA, ifelse(as.numeric(df$presence) > 0, "detected", "undetected"))
  df$plot_value <- factor(df$plot_value, levels = c("undetected", "detected"))
  df$tile_fill <- ifelse(is.na(df$plot_value), "#d9d9d9", ortholog_colors[as.character(df$plot_value)])
}

max_label_width <- max(nchar(as.character(label_df$species_display), type = "width"))
max_query_width <- max(nchar(as.character(query_levels), type = "width"))
busco_df <- read_busco_summary_table(busco_table_path, tip_y_by_species)

tree_width <- 10.5
tree_left <- 0
tree_map <- function(x) tree_left + (x - x_left) / (x_right - x_left) * tree_width
tree_right <- tree_left + tree_width
label_right <- tree_right + max(4.0, max_label_width * 0.33)
heatmap_left <- label_right + 0.02
heatmap_right <- heatmap_left + length(query_levels) * heatmap_cell_pitch
busco_left <- heatmap_right + ifelse(nrow(busco_df) > 0, 0.95, 0.15)
busco_width <- ifelse(nrow(busco_df) > 0, 3.6, 0)
busco_right <- busco_left + busco_width
x_max <- busco_right + ifelse(nrow(busco_df) > 0, 0.55, 0.15)
x_min <- tree_left - 0.25

node_xy <- plot_data[, c("node", "parent", "x", "y"), drop = FALSE]
parent_xy <- node_xy[, c("node", "x", "y"), drop = FALSE]
colnames(parent_xy) <- c("parent", "parent_x", "parent_y")
edge_df <- merge(node_xy, parent_xy, by = "parent", all.x = FALSE, all.y = FALSE)
edge_df <- edge_df[edge_df$node != edge_df$parent, , drop = FALSE]
edge_df$x <- tree_map(edge_df$x)
edge_df$parent_x <- tree_map(edge_df$parent_x)

ci_plot_df$xmin_plot <- tree_map(ci_plot_df$xmin)
ci_plot_df$xmax_plot <- tree_map(ci_plot_df$xmax)
support_df$x_plot <- tree_map(support_df$branch)
ruler_df <- data.frame(
  x = tree_map(ruler_breaks),
  label = ruler_labels,
  stringsAsFactors = FALSE
)

heatmap_df <- df
heatmap_df$x <- heatmap_left + (heatmap_df$query_index - 0.5) * heatmap_cell_pitch
heatmap_cell_half <- heatmap_cell_pitch * 0.46
heatmap_df$xmin <- heatmap_df$x - heatmap_cell_half
heatmap_df$xmax <- heatmap_df$x + heatmap_cell_half
heatmap_df$ymin <- heatmap_df$y - heatmap_cell_half
heatmap_df$ymax <- heatmap_df$y + heatmap_cell_half
query_label_df <- data.frame(
  query = query_levels,
  x = heatmap_left + (seq_along(query_levels) - 0.5) * heatmap_cell_pitch,
  stringsAsFactors = FALSE
)

busco_rect_df <- data.frame()
if (nrow(busco_df) > 0) {
  busco_rect_df <- busco_df
  busco_rect_df$xmin <- busco_left + busco_rect_df$xmin_percent / 100 * busco_width
  busco_rect_df$xmax <- busco_left + busco_rect_df$xmax_percent / 100 * busco_width
  busco_rect_df$ymin <- busco_rect_df$y - heatmap_cell_pitch * 0.36
  busco_rect_df$ymax <- busco_rect_df$y + heatmap_cell_pitch * 0.36
  busco_rect_df$fill <- busco_colors[busco_rect_df$status]
}

row_y_min <- min(y_values, na.rm = TRUE)
row_y_max <- max(y_values, na.rm = TRUE)
y_ruler <- row_y_min - 0.85
y_axis_label <- row_y_min - 2.35
y_query_label <- row_y_min - heatmap_cell_half - 0.06
y_legend <- row_y_min - 4.15
y_busco_axis <- row_y_min - 0.85
y_busco_label <- row_y_min - 2.35
y_min <- row_y_min - 6.70
y_max <- row_y_max + 0.75

legend_df <- data.frame()
if (value_mode == "presence") {
  legend_df <- data.frame(
    label = c("undetected", "detected"),
    fill = unname(ortholog_colors[c("undetected", "detected")]),
    x = heatmap_left + c(2.95, 7.35),
    stringsAsFactors = FALSE
  )
}
busco_legend_df <- data.frame()
if (nrow(busco_df) > 0) {
  busco_legend_df <- data.frame(
    label = names(busco_colors),
    fill = unname(busco_colors),
    x = busco_left + 0.35,
    y = y_legend - seq(0, by = 0.75, length.out = length(busco_colors)),
    stringsAsFactors = FALSE
  )
}

combined <- ggplot()
if (nrow(ci_plot_df) > 0) {
  combined <- combined +
    geom_segment(data = ci_plot_df, aes(x = xmin_plot, xend = xmax_plot, y = y, yend = y), linewidth = 0.8, color = hpd_color, alpha = 0.75, lineend = "butt")
}
combined <- combined +
  geom_segment(data = edge_df, aes(x = parent_x, xend = x, y = y, yend = y), linewidth = 0.28, color = "black", lineend = "square") +
  geom_segment(data = edge_df, aes(x = parent_x, xend = parent_x, y = parent_y, yend = y), linewidth = 0.28, color = "black", lineend = "square")
if (nrow(support_df) > 0) {
  combined <- combined +
    geom_text(data = support_df, aes(x = x_plot, y = y + 0.28, label = support_label), size = font_size_mm, color = support_color)
}
combined <- combined +
  geom_text(data = label_df, aes(x = label_right, y = y, label = species_display), hjust = 1, size = font_size_mm, color = "black", fontface = "italic") +
  geom_rect(data = heatmap_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = tile_fill), color = "white", linewidth = 0.18) +
  geom_text(data = query_label_df, aes(x = x, y = y_query_label, label = query), angle = 90, hjust = 1, vjust = 0.5, size = font_size_mm, color = "black") +
  geom_segment(aes(x = tree_map(x_left + x_pad), xend = tree_map(x_right - x_pad), y = y_ruler, yend = y_ruler), linewidth = 0.28, color = "black") +
  geom_segment(data = ruler_df, aes(x = x, xend = x, y = y_ruler, yend = y_ruler - 0.15), linewidth = 0.28, color = "black") +
  geom_text(data = ruler_df, aes(x = x, y = y_ruler - 0.30, label = label), size = font_size_mm, vjust = 1, color = "black") +
  annotate("text", x = mean(c(tree_left, tree_right)), y = y_axis_label, label = "Million years ago", size = font_size_mm, color = "black") +
  scale_fill_identity() +
  coord_fixed(ratio = 1, xlim = c(x_min, x_max), ylim = c(y_min, y_max), expand = FALSE, clip = "off") +
  theme_void(base_size = font_size_pt) +
  theme(plot.margin = margin(3, 3, 3, 3))

if (nrow(busco_rect_df) > 0) {
  busco_tick_df <- data.frame(
    pct = c(0, 50, 100),
    x = busco_left + c(0, 50, 100) / 100 * busco_width,
    stringsAsFactors = FALSE
  )
  combined <- combined +
    geom_rect(data = busco_rect_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill), color = "white", linewidth = 0.14) +
    geom_segment(aes(x = busco_left, xend = busco_right, y = y_busco_axis, yend = y_busco_axis), linewidth = 0.28, color = "black") +
    geom_segment(data = busco_tick_df, aes(x = x, xend = x, y = y_busco_axis, yend = y_busco_axis - 0.15), linewidth = 0.28, color = "black") +
    geom_text(data = busco_tick_df, aes(x = x, y = y_busco_axis - 0.30, label = pct), size = font_size_mm, vjust = 1, color = "black") +
    annotate("text", x = mean(c(busco_left, busco_right)), y = y_busco_label, label = "BUSCO (%)", size = font_size_mm, color = "black")
}

if (value_mode == "presence") {
  combined <- combined +
    annotate("text", x = heatmap_left + 0.05, y = y_legend, label = "Ortholog", hjust = 0, size = font_size_mm, color = "black") +
    geom_rect(data = legend_df, aes(xmin = x, xmax = x + 0.36, ymin = y_legend - 0.18, ymax = y_legend + 0.18, fill = fill), color = "black", linewidth = 0.14) +
    geom_text(data = legend_df, aes(x = x + 0.62, y = y_legend, label = label), hjust = 0, size = font_size_mm, color = "black")
}
if (nrow(support_df) > 0) {
  combined <- combined +
    annotate("text", x = tree_left + 1.8, y = y_legend, label = "00", hjust = 0, size = font_size_mm, color = support_color) +
    annotate("text", x = tree_left + 2.40, y = y_legend, label = ":", hjust = 0, size = font_size_mm, color = "black") +
    annotate("text", x = tree_left + 2.85, y = y_legend, label = "Ultrafast bootstrap support (%)", hjust = 0, size = font_size_mm, color = "black")
}
if (nrow(busco_legend_df) > 0) {
  combined <- combined +
    annotate("text", x = busco_left, y = y_legend + 0.75, label = "BUSCO", hjust = 0, size = font_size_mm, color = "black") +
    geom_rect(data = busco_legend_df, aes(xmin = x, xmax = x + 0.32, ymin = y - 0.16, ymax = y + 0.16, fill = fill), color = "white", linewidth = 0.14) +
    geom_text(data = busco_legend_df, aes(x = x + 0.42, y = y, label = label), hjust = 0, size = font_size_mm, color = "black")
}

plot_height <- if (identical(plot_height_arg, "auto")) {
  max(2.8, plot_width * (y_max - y_min) / (x_max - x_min))
} else {
  as.numeric(plot_height_arg)
}
if (!is.finite(plot_height) || plot_height <= 0) {
  stop("Invalid --height: ", plot_height_arg)
}

if (nzchar(out_pdf)) {
  dir.create(dirname(out_pdf), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_pdf, plot = combined, width = plot_width, height = plot_height, units = "in", dpi = 300, limitsize = FALSE, bg = "transparent")
  cat("Wrote PDF:", out_pdf, "\n")
}
if (nzchar(out_svg)) {
  dir.create(dirname(out_svg), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_svg, plot = combined, width = plot_width, height = plot_height, units = "in", dpi = 300, limitsize = FALSE, bg = "transparent")
  cat("Wrote SVG:", out_svg, "\n")
}
