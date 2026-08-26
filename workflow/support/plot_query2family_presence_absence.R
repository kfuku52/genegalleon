#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ape, quietly = TRUE))
suppressPackageStartupMessages(library(cowplot, quietly = TRUE))
suppressPackageStartupMessages(library(ggplot2, quietly = TRUE))
suppressPackageStartupMessages(library(ggtree, quietly = TRUE))

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE))
} else {
  getwd()
}
source(file.path(script_dir, "species_label_utils.r"), local = TRUE)

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

transfer_species_mapping_labels <- function(target_tree, mapping_tree) {
  target_tips <- as.character(target_tree[["tip.label"]])
  mapping_tips <- as.character(mapping_tree[["tip.label"]])
  if (!setequal(target_tips, mapping_tips)) {
    stop(
      "Species mapping tree tips do not match the plotted species tree: ",
      "only_in_plot=",
      paste(setdiff(target_tips, mapping_tips), collapse = ","),
      "; only_in_mapping=",
      paste(setdiff(mapping_tips, target_tips), collapse = ",")
    )
  }

  mapping_labels <- mapping_tree[["node.label"]]
  if (is.null(mapping_labels)) {
    mapping_labels <- rep("", mapping_tree[["Nnode"]])
  }
  mapping_nodes <- seq.int(
    length(mapping_tips) + 1,
    length(mapping_tips) + mapping_tree[["Nnode"]]
  )
  label_by_clade <- list()
  for (index in seq_along(mapping_nodes)) {
    label <- if (index <= length(mapping_labels)) {
      as.character(mapping_labels[[index]])
    } else {
      ""
    }
    if (is.na(label) || !nzchar(label)) {
      next
    }
    clade_key <- node_clade_key(
      mapping_tree,
      mapping_nodes[[index]],
      tip_set = target_tips
    )
    if (!is.na(clade_key) && nzchar(clade_key)) {
      label_by_clade[[clade_key]] <- label
    }
  }

  target_labels <- target_tree[["node.label"]]
  if (is.null(target_labels)) {
    target_labels <- rep("", target_tree[["Nnode"]])
  }
  target_nodes <- seq.int(
    length(target_tips) + 1,
    length(target_tips) + target_tree[["Nnode"]]
  )
  for (index in seq_along(target_nodes)) {
    clade_key <- node_clade_key(
      target_tree,
      target_nodes[[index]],
      tip_set = target_tips
    )
    if (!is.na(clade_key) && !is.null(label_by_clade[[clade_key]])) {
      target_labels[[index]] <- label_by_clade[[clade_key]]
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

nice_count_breaks <- function(max_count, target_intervals = 3) {
  max_count <- suppressWarnings(as.numeric(max_count))
  max_count <- max_count[is.finite(max_count)]
  if (length(max_count) == 0 || max(max_count) <= 0) {
    return(c(0, 1))
  }
  max_count <- ceiling(max(max_count))
  breaks <- pretty(c(0, max_count), n = target_intervals)
  breaks <- breaks[is.finite(breaks) & breaks >= 0]
  if (any(abs(breaks - round(breaks)) > 1e-8)) {
    breaks <- seq.int(0, max_count)
  }
  breaks <- sort(unique(as.numeric(round(breaks))))
  if (length(breaks) == 0 || max(breaks) < max_count) {
    step <- if (length(breaks) >= 2) breaks[[2]] - breaks[[1]] else 1
    step <- max(1, step)
    upper <- ceiling(max_count / step) * step
    breaks <- seq(0, upper, by = step)
  }
  breaks
}

nice_anchor_breaks <- function(max_anchor_count, max_labels = 6) {
  max_anchor_count <- suppressWarnings(as.numeric(max_anchor_count))
  max_anchor_count <- max_anchor_count[is.finite(max_anchor_count)]
  if (length(max_anchor_count) == 0 || max(max_anchor_count) < 1) {
    return(1)
  }
  max_anchor_count <- as.integer(ceiling(max(max_anchor_count)))
  if (max_anchor_count <= max_labels) {
    return(seq_len(max_anchor_count))
  }
  breaks <- pretty(c(1, max_anchor_count), n = max_labels)
  breaks <- as.integer(round(breaks[
    is.finite(breaks) & breaks >= 1 & breaks <= max_anchor_count
  ]))
  sort(unique(c(1L, breaks, max_anchor_count)))
}

viridis_unit_colors <- function(fraction) {
  fraction <- suppressWarnings(as.numeric(fraction))
  palette <- grDevices::hcl.colors(256, palette = "viridis")
  palette_index <- 1L + as.integer(round(pmax(0, pmin(1, fraction)) * 255))
  palette[palette_index]
}

viridis_anchor_colors <- function(anchor_count, max_anchor_count) {
  anchor_count <- suppressWarnings(as.numeric(anchor_count))
  max_anchor_count <- suppressWarnings(as.numeric(max_anchor_count))
  max_anchor_count <- max(1, max_anchor_count[is.finite(max_anchor_count)], na.rm = TRUE)
  fraction <- if (max_anchor_count <= 1) {
    rep(0.5, length(anchor_count))
  } else {
    (anchor_count - 1) / (max_anchor_count - 1)
  }
  viridis_unit_colors(fraction)
}

viridis_anchor_rail_colors <- function(anchor_count, max_anchor_count) {
  anchor_count <- suppressWarnings(as.numeric(anchor_count))
  max_anchor_count <- suppressWarnings(as.numeric(max_anchor_count))
  max_anchor_count <- max(1, max_anchor_count[is.finite(max_anchor_count)], na.rm = TRUE)
  viridis_unit_colors(anchor_count / max_anchor_count)
}

viridis_ufboot_colors <- function(ufboot) {
  ufboot <- suppressWarnings(as.numeric(ufboot))
  fraction <- pmax(0, pmin(100, ufboot)) / 100
  viridis_unit_colors(fraction)
}

inferno_ufboot_colors <- function(ufboot) {
  ufboot <- suppressWarnings(as.numeric(ufboot))
  fraction <- pmax(0, pmin(100, ufboot)) / 100
  palette <- grDevices::hcl.colors(256, palette = "Inferno")
  palette_index <- 1L + as.integer(round(fraction * 255))
  palette[palette_index]
}

format_ruler_labels <- function(x) {
  x <- as.numeric(x)
  max_value <- max(abs(x), na.rm = TRUE)
  digits <- if (is.finite(max_value) && max_value < 10) 1 else 0
  formatC(x, format = "f", digits = digits)
}

format_scale_bar_label <- function(x) {
  x <- as.numeric(x)
  if (!is.finite(x) || x == 0) {
    return("0")
  }
  digits <- max(0, -floor(log10(abs(x))))
  formatC(x, format = "f", digits = min(digits, 8))
}

normalize_busco_species_id <- function(x) {
  gg_species_label_from_filename(x)
}

match_busco_species <- function(species, tip_y_by_species) {
  species <- normalize_busco_species_id(species)
  display_species <- gg_species_display_from_key(species)
  candidates <- unique(c(
    species,
    display_species,
    gsub(" ", "_", species, fixed = TRUE),
    gsub("_", " ", species, fixed = TRUE),
    gsub(" ", "_", display_species, fixed = TRUE)
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
species_mapping_tree_path <- get_arg(args, "species_mapping_tree")
busco_table_path <- get_arg(args, "busco_table")
long_path <- get_arg(args, "long_table")
ortholog_column_path <- get_arg(args, "ortholog_column_table")
ortholog_glyph_path <- get_arg(args, "ortholog_glyph_table")
ortholog_tree_path <- get_arg(args, "ortholog_tree_table")
ortholog_synteny_path <- get_arg(args, "ortholog_synteny_table")
ortholog_ufboot_path <- get_arg(args, "ortholog_ufboot_table")
reference_species <- get_arg(args, "reference_species")
out_pdf <- get_arg(args, "out_pdf")
out_svg <- get_arg(args, "out_svg")
value_mode <- get_arg(args, "value", "presence")
tree_scale_mode <- get_arg(args, "tree_scale", "bar")
plot_width <- as.numeric(get_arg(args, "width", "7.2"))
plot_height_arg <- get_arg(args, "height", "auto")
evidence_layout <- tolower(get_arg(args, "evidence_layout", "band"))
glyph_mode <- has_nonempty_file(ortholog_column_path) && has_nonempty_file(ortholog_glyph_path)
reference_species_display <- ""
if (glyph_mode) {
  if (!nzchar(reference_species)) {
    stop("--reference_species is required with ortholog tables")
  }
  reference_species_display <- gg_species_display_from_key(reference_species)
}

if (!nzchar(tree_path) || !file.exists(tree_path)) {
  stop("Missing --species_tree: ", tree_path)
}
if (!nzchar(long_path) || !file.exists(long_path)) {
  stop("Missing --long_table: ", long_path)
}
if (!value_mode %in% c("presence", "copy_number")) {
  stop("--value must be presence or copy_number: ", value_mode)
}
if (!tree_scale_mode %in% c("bar", "ruler")) {
  stop("--tree_scale must be bar or ruler: ", tree_scale_mode)
}
if (!evidence_layout %in% c("band", "rail", "glyph", "off")) {
  stop("--evidence_layout must be band, rail, glyph, or off: ", evidence_layout)
}
evidence_strip_mode <- evidence_layout %in% c("band", "rail")
evidence_layout_label <- if (identical(evidence_layout, "band")) {
  "Evidence band"
} else if (identical(evidence_layout, "rail")) {
  "Evidence rail"
} else {
  ""
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
tree_for_species_mapping <- tree
if (glyph_mode && nzchar(species_mapping_tree_path)) {
  if (!has_nonempty_file(species_mapping_tree_path)) {
    stop("Missing --species_mapping_tree: ", species_mapping_tree_path)
  }
  tree_for_species_mapping <- transfer_species_mapping_labels(
    tree,
    ape::read.tree(species_mapping_tree_path)
  )
}

df <- read.table(long_path, sep = "\t", header = TRUE, quote = "", comment.char = "", check.names = FALSE)
if (nrow(df) == 0) {
  stop("Long table is empty: ", long_path)
}
source_df <- df
ortholog_columns <- data.frame()
ortholog_glyphs <- data.frame()
ortholog_tree_nodes <- data.frame()
ortholog_synteny <- data.frame()
ortholog_ufboot <- data.frame()
ortholog_family_count <- 0
duplication_family_mode <- FALSE
family_colors <- character(0)
if (glyph_mode) {
  ortholog_columns <- read.table(ortholog_column_path, sep = "\t", header = TRUE, quote = "", comment.char = "", check.names = FALSE)
  ortholog_glyphs <- read.table(ortholog_glyph_path, sep = "\t", header = TRUE, quote = "", comment.char = "", check.names = FALSE)
  if (has_nonempty_file(ortholog_tree_path)) {
    ortholog_tree_nodes <- read.table(ortholog_tree_path, sep = "\t", header = TRUE, quote = "", comment.char = "", check.names = FALSE)
  }
  if (has_nonempty_file(ortholog_synteny_path)) {
    ortholog_synteny <- read.table(ortholog_synteny_path, sep = "\t", header = TRUE, quote = "", comment.char = "", check.names = FALSE)
  }
  if (has_nonempty_file(ortholog_ufboot_path)) {
    ortholog_ufboot <- read.table(ortholog_ufboot_path, sep = "\t", header = TRUE, quote = "", comment.char = "", check.names = FALSE)
  }
  required_column_fields <- c(
    "column_order", "family_id", "family_order", "reference_species",
    "cds_fasta_id", "gene_id", "plot_label"
  )
  required_glyph_fields <- c(
    "species", "family_id", "family_order", "relation", "copy_number",
    "start_order", "end_order", "lane_index", "lane_count"
  )
  if (!all(required_column_fields %in% colnames(ortholog_columns))) {
    stop("Ortholog column table is missing required columns: ", ortholog_column_path)
  }
  if (!all(required_glyph_fields %in% colnames(ortholog_glyphs))) {
    stop("Ortholog glyph table is missing required columns: ", ortholog_glyph_path)
  }
  required_tree_fields <- c(
    "family_id", "family_order", "node_id", "parent_node_id", "is_tip", "event",
    "node_height", "plot_order", "mapped_species_node", "duplication_index",
    "in_reference_tree"
  )
  if (nrow(ortholog_tree_nodes) > 0 && !all(required_tree_fields %in% colnames(ortholog_tree_nodes))) {
    stop("Ortholog tree table is missing required columns: ", ortholog_tree_path)
  }
  if (nrow(ortholog_columns) == 0) {
    stop("Ortholog column table is empty: ", ortholog_column_path)
  }
  if (nrow(ortholog_glyphs) == 0) {
    stop("Ortholog glyph table is empty: ", ortholog_glyph_path)
  }
  column_order_values <- suppressWarnings(as.numeric(ortholog_columns$column_order))
  if (
    any(!is.finite(column_order_values)) ||
      any(abs(column_order_values - round(column_order_values)) > 1e-8) ||
      !identical(sort(as.integer(column_order_values)), seq_len(nrow(ortholog_columns)))
  ) {
    stop("Ortholog column_order must contain each integer from 1 through the number of columns exactly once")
  }
  ortholog_columns$column_order <- as.integer(column_order_values)
  ortholog_columns$family_id <- as.character(ortholog_columns$family_id)
  ortholog_columns$family_order <- suppressWarnings(as.numeric(ortholog_columns$family_order))
  ortholog_columns$reference_species <- as.character(ortholog_columns$reference_species)
  ortholog_columns$cds_fasta_id <- as.character(ortholog_columns$cds_fasta_id)
  ortholog_columns$gene_id <- as.character(ortholog_columns$gene_id)
  ortholog_columns$plot_label <- as.character(ortholog_columns$plot_label)
  if (any(is.na(ortholog_columns$family_id) | !nzchar(ortholog_columns$family_id))) {
    stop("Ortholog column family_id values must be non-empty")
  }
  if (any(is.na(ortholog_columns$plot_label) | !nzchar(ortholog_columns$plot_label))) {
    stop("Ortholog column plot_label values must be non-empty")
  }
  if (
    any(is.na(ortholog_columns$reference_species) | !nzchar(ortholog_columns$reference_species)) ||
      any(ortholog_columns$reference_species != reference_species) ||
      any(is.na(ortholog_columns$cds_fasta_id) | !nzchar(ortholog_columns$cds_fasta_id)) ||
      any(is.na(ortholog_columns$gene_id) | !nzchar(ortholog_columns$gene_id))
  ) {
    stop("Ortholog column reference species and gene identifiers must be non-empty and reference-consistent")
  }
  if (
    any(!is.finite(ortholog_columns$family_order)) ||
      any(ortholog_columns$family_order < 1) ||
      any(abs(ortholog_columns$family_order - round(ortholog_columns$family_order)) > 1e-8)
  ) {
    stop("Ortholog column family_order values must be positive integers")
  }
  ortholog_columns$family_order <- as.integer(ortholog_columns$family_order)
  ortholog_columns <- ortholog_columns[
    order(ortholog_columns$column_order),
    ,
    drop = FALSE
  ]
  family_pairs <- unique(
    ortholog_columns[, c("family_id", "family_order"), drop = FALSE]
  )
  if (anyDuplicated(family_pairs$family_id) || anyDuplicated(family_pairs$family_order)) {
    stop("Ortholog columns must map each family_id to exactly one unique family_order")
  }
  family_runs <- rle(ortholog_columns$family_id)$values
  if (anyDuplicated(family_runs)) {
    stop("Ortholog columns for each family must form one contiguous block")
  }
  family_order_by_id <- stats::setNames(family_pairs$family_order, family_pairs$family_id)
  family_id_by_column <- stats::setNames(
    ortholog_columns$family_id,
    as.character(ortholog_columns$column_order)
  )
  reference_species_by_column <- stats::setNames(
    ortholog_columns$reference_species,
    as.character(ortholog_columns$column_order)
  )
  reference_cds_id_by_column <- stats::setNames(
    ortholog_columns$cds_fasta_id,
    as.character(ortholog_columns$column_order)
  )
  reference_gene_id_by_column <- stats::setNames(
    ortholog_columns$gene_id,
    as.character(ortholog_columns$column_order)
  )

  ortholog_glyphs$family_id <- as.character(ortholog_glyphs$family_id)
  ortholog_glyphs$family_order <- suppressWarnings(as.numeric(ortholog_glyphs$family_order))
  ortholog_glyphs$relation <- as.character(ortholog_glyphs$relation)
  glyph_numeric_fields <- c(
    "copy_number", "start_order", "end_order", "lane_index", "lane_count"
  )
  for (field in glyph_numeric_fields) {
    ortholog_glyphs[[field]] <- suppressWarnings(as.numeric(ortholog_glyphs[[field]]))
  }
  if (
    any(is.na(ortholog_glyphs$family_id) | !nzchar(ortholog_glyphs$family_id)) ||
      any(!ortholog_glyphs$family_id %in% names(family_order_by_id))
  ) {
    stop("Ortholog glyph family_id values must match the ortholog column table")
  }
  expected_glyph_family_order <- unname(family_order_by_id[ortholog_glyphs$family_id])
  if (
    any(!is.finite(ortholog_glyphs$family_order)) ||
      any(ortholog_glyphs$family_order != expected_glyph_family_order)
  ) {
    stop("Ortholog glyph family_order values disagree with the ortholog column table")
  }
  allowed_glyph_relations <- c("specific", "shared_ancestral", "ambiguous")
  if (any(is.na(ortholog_glyphs$relation) | !ortholog_glyphs$relation %in% allowed_glyph_relations)) {
    stop("Ortholog glyph relation must be specific, shared_ancestral, or ambiguous")
  }
  if (
    any(!is.finite(ortholog_glyphs$copy_number)) ||
      any(ortholog_glyphs$copy_number < 1) ||
      any(abs(ortholog_glyphs$copy_number - round(ortholog_glyphs$copy_number)) > 1e-8)
  ) {
    stop("Ortholog glyph copy_number values must be positive integers")
  }
  for (field in c("start_order", "end_order", "lane_index", "lane_count")) {
    values <- ortholog_glyphs[[field]]
    if (any(!is.finite(values)) || any(abs(values - round(values)) > 1e-8)) {
      stop("Ortholog glyph ", field, " values must be integers")
    }
    ortholog_glyphs[[field]] <- as.integer(values)
  }
  if (
    any(ortholog_glyphs$start_order < 1) ||
      any(ortholog_glyphs$end_order > nrow(ortholog_columns)) ||
      any(ortholog_glyphs$start_order > ortholog_glyphs$end_order)
  ) {
    stop("Ortholog glyph column spans are outside the ortholog column table or reversed")
  }
  if (
    any(ortholog_glyphs$lane_index < 1) ||
      any(ortholog_glyphs$lane_count < ortholog_glyphs$lane_index)
  ) {
    stop("Ortholog glyph lanes must satisfy 1 <= lane_index <= lane_count")
  }
  glyph_start_family <- unname(
    family_id_by_column[as.character(ortholog_glyphs$start_order)]
  )
  glyph_end_family <- unname(
    family_id_by_column[as.character(ortholog_glyphs$end_order)]
  )
  if (
    any(ortholog_glyphs$family_id != glyph_start_family) ||
      any(ortholog_glyphs$family_id != glyph_end_family)
  ) {
    stop("Ortholog glyph column spans must remain within their declared family")
  }

  if (nrow(ortholog_synteny) > 0) {
    required_synteny_fields <- c(
      "family_id", "family_order", "species", "reference_species", "relation",
      "reference_cds_fasta_id", "reference_gene_id", "column_order",
      "candidate_cds_fasta_id", "glyph_copy_number",
      "glyph_start_order", "glyph_end_order", "glyph_lane_index", "glyph_lane_count",
      "synteny_status", "support_min_anchor_count", "synteny_window_radius",
      "shared_anchor_count"
    )
    if (!all(required_synteny_fields %in% colnames(ortholog_synteny))) {
      stop("Ortholog synteny table is missing required columns: ", ortholog_synteny_path)
    }
    character_synteny_fields <- c(
      "family_id", "species", "reference_species", "relation",
      "reference_cds_fasta_id", "reference_gene_id", "candidate_cds_fasta_id",
      "synteny_status"
    )
    for (field in character_synteny_fields) {
      ortholog_synteny[[field]] <- as.character(ortholog_synteny[[field]])
      ortholog_synteny[[field]][is.na(ortholog_synteny[[field]])] <- ""
    }
    numeric_synteny_fields <- c(
      "family_order", "column_order", "glyph_copy_number", "glyph_start_order",
      "glyph_end_order", "glyph_lane_index", "glyph_lane_count",
      "support_min_anchor_count", "synteny_window_radius", "shared_anchor_count"
    )
    for (field in numeric_synteny_fields) {
      ortholog_synteny[[field]] <- suppressWarnings(as.numeric(ortholog_synteny[[field]]))
    }
    if (
      any(!nzchar(ortholog_synteny$family_id)) ||
        any(!ortholog_synteny$family_id %in% names(family_order_by_id))
    ) {
      stop("Ortholog synteny family_id values must match the ortholog column table")
    }
    if (
      any(!nzchar(ortholog_synteny$species)) ||
        any(!nzchar(ortholog_synteny$reference_species)) ||
        any(ortholog_synteny$reference_species != reference_species) ||
        any(!nzchar(ortholog_synteny$reference_cds_fasta_id)) ||
        any(!nzchar(ortholog_synteny$reference_gene_id)) ||
        any(!nzchar(ortholog_synteny$candidate_cds_fasta_id))
    ) {
      stop("Ortholog synteny species and gene identifiers must be non-empty and reference-consistent")
    }
    expected_synteny_family_order <- unname(
      family_order_by_id[ortholog_synteny$family_id]
    )
    finite_integer_fields <- c(
      "family_order", "column_order", "glyph_copy_number", "glyph_start_order",
      "glyph_end_order", "glyph_lane_index", "glyph_lane_count",
      "support_min_anchor_count"
    )
    for (field in finite_integer_fields) {
      values <- ortholog_synteny[[field]]
      if (any(!is.finite(values)) || any(abs(values - round(values)) > 1e-8)) {
        stop("Ortholog synteny ", field, " values must be finite integers")
      }
      ortholog_synteny[[field]] <- as.integer(values)
    }
    if (any(ortholog_synteny$family_order != expected_synteny_family_order)) {
      stop("Ortholog synteny family_order values disagree with the ortholog column table")
    }
    expected_synteny_family_id <- unname(
      family_id_by_column[as.character(ortholog_synteny$column_order)]
    )
    if (
      any(is.na(expected_synteny_family_id)) ||
        any(ortholog_synteny$family_id != expected_synteny_family_id)
    ) {
      stop("Ortholog synteny column_order values cross family boundaries")
    }
    expected_synteny_reference_species <- unname(
      reference_species_by_column[as.character(ortholog_synteny$column_order)]
    )
    expected_synteny_reference_cds_id <- unname(
      reference_cds_id_by_column[as.character(ortholog_synteny$column_order)]
    )
    expected_synteny_reference_gene_id <- unname(
      reference_gene_id_by_column[as.character(ortholog_synteny$column_order)]
    )
    if (
      any(ortholog_synteny$reference_species != expected_synteny_reference_species) ||
        any(ortholog_synteny$reference_cds_fasta_id != expected_synteny_reference_cds_id) ||
        any(ortholog_synteny$reference_gene_id != expected_synteny_reference_gene_id)
    ) {
      stop("Ortholog synteny reference identifiers disagree with the ortholog column table")
    }
    if (
      any(ortholog_synteny$column_order < ortholog_synteny$glyph_start_order) ||
        any(ortholog_synteny$column_order > ortholog_synteny$glyph_end_order) ||
        any(ortholog_synteny$glyph_copy_number < 1) ||
        any(ortholog_synteny$glyph_lane_index < 1) ||
        any(ortholog_synteny$glyph_lane_count < ortholog_synteny$glyph_lane_index)
    ) {
      stop("Ortholog synteny rows contain invalid glyph spans, copy numbers, or lanes")
    }
    allowed_synteny_status <- c(
      "supported", "single_anchor", "no_support", "not_evaluable", "reference_self"
    )
    if (
      any(is.na(ortholog_synteny$synteny_status)) ||
        any(!ortholog_synteny$synteny_status %in% allowed_synteny_status)
    ) {
      stop("Ortholog synteny table contains an invalid synteny_status")
    }
    is_synteny_reference_self <- ortholog_synteny$candidate_cds_fasta_id ==
      ortholog_synteny$reference_cds_fasta_id
    if (
      any(is_synteny_reference_self & ortholog_synteny$synteny_status != "reference_self") ||
        any(!is_synteny_reference_self & ortholog_synteny$synteny_status == "reference_self")
    ) {
      stop("Ortholog synteny reference_self status disagrees with candidate/reference IDs")
    }
    if (any(ortholog_synteny$support_min_anchor_count != 2)) {
      stop("Ortholog synteny support_min_anchor_count must be 2")
    }
    evaluated_synteny_rows <- ortholog_synteny$synteny_status %in% c(
      "supported", "single_anchor", "no_support"
    )
    finite_window_rows <- is.finite(ortholog_synteny$synteny_window_radius)
    if (
      any(
        finite_window_rows &
          (
            ortholog_synteny$synteny_window_radius < 1 |
              abs(
                ortholog_synteny$synteny_window_radius -
                  round(ortholog_synteny$synteny_window_radius)
              ) > 1e-8
          )
      ) ||
        any(evaluated_synteny_rows & !finite_window_rows)
    ) {
      stop(
        "Evaluated ortholog synteny rows require a positive-integer ",
        "synteny_window_radius"
      )
    }
    ortholog_synteny$synteny_window_radius[finite_window_rows] <- as.integer(
      ortholog_synteny$synteny_window_radius[finite_window_rows]
    )
    if (
      any(!is.finite(ortholog_synteny$shared_anchor_count[evaluated_synteny_rows])) ||
        any(ortholog_synteny$shared_anchor_count[evaluated_synteny_rows] < 0) ||
        any(abs(
          ortholog_synteny$shared_anchor_count[evaluated_synteny_rows] -
            round(ortholog_synteny$shared_anchor_count[evaluated_synteny_rows])
        ) > 1e-8)
    ) {
      stop("Evaluated ortholog synteny shared_anchor_count values must be non-negative integers")
    }
    if (
      any(ortholog_synteny$shared_anchor_count[ortholog_synteny$synteny_status == "supported"] < 2) ||
        any(ortholog_synteny$shared_anchor_count[ortholog_synteny$synteny_status == "single_anchor"] != 1) ||
        any(ortholog_synteny$shared_anchor_count[ortholog_synteny$synteny_status == "no_support"] != 0)
    ) {
      stop("Ortholog synteny status values disagree with shared_anchor_count")
    }
    synteny_pair_keys <- paste(
      ortholog_synteny$family_id,
      ortholog_synteny$reference_cds_fasta_id,
      ortholog_synteny$candidate_cds_fasta_id,
      sep = "::"
    )
    if (anyDuplicated(synteny_pair_keys)) {
      stop("Ortholog synteny table contains duplicate candidate/reference pairs within a family")
    }
    glyph_keys <- paste(
      ortholog_glyphs$family_id,
      ortholog_glyphs$species,
      ortholog_glyphs$start_order,
      ortholog_glyphs$end_order,
      ortholog_glyphs$lane_index,
      ortholog_glyphs$lane_count,
      ortholog_glyphs$copy_number,
      sep = "::"
    )
    synteny_glyph_keys <- paste(
      ortholog_synteny$family_id,
      ortholog_synteny$species,
      ortholog_synteny$glyph_start_order,
      ortholog_synteny$glyph_end_order,
      ortholog_synteny$glyph_lane_index,
      ortholog_synteny$glyph_lane_count,
      ortholog_synteny$glyph_copy_number,
      sep = "::"
    )
    if (any(!synteny_glyph_keys %in% glyph_keys)) {
      stop("Ortholog synteny rows do not map to an ortholog glyph")
    }
  }

  if (nrow(ortholog_ufboot) > 0) {
    required_ufboot_fields <- c(
      "family_id", "family_order", "species", "reference_species", "relation",
      "reference_cds_fasta_id", "reference_gene_id", "column_order",
      "candidate_cds_fasta_id",
      "glyph_copy_number", "glyph_start_order", "glyph_end_order",
      "glyph_lane_index", "glyph_lane_count", "orthology_mrca_branch_id",
      "orthology_mrca_event", "ufboot_support_source", "decisive_branch_ufboot",
      "orthology_ufboot_status", "orthology_ufboot_unavailable_reason"
    )
    if (!all(required_ufboot_fields %in% colnames(ortholog_ufboot))) {
      stop("Ortholog UFBoot table is missing required columns: ", ortholog_ufboot_path)
    }
    character_ufboot_fields <- c(
      "family_id", "species", "reference_species", "relation",
      "reference_cds_fasta_id", "reference_gene_id", "candidate_cds_fasta_id",
      "orthology_mrca_branch_id",
      "orthology_mrca_event", "ufboot_support_source", "orthology_ufboot_status",
      "orthology_ufboot_unavailable_reason"
    )
    for (field in character_ufboot_fields) {
      ortholog_ufboot[[field]] <- as.character(ortholog_ufboot[[field]])
      ortholog_ufboot[[field]][is.na(ortholog_ufboot[[field]])] <- ""
    }
    integer_ufboot_fields <- c(
      "family_order", "column_order", "glyph_copy_number",
      "glyph_start_order", "glyph_end_order", "glyph_lane_index",
      "glyph_lane_count"
    )
    for (field in integer_ufboot_fields) {
      values <- suppressWarnings(as.numeric(ortholog_ufboot[[field]]))
      if (any(!is.finite(values)) || any(abs(values - round(values)) > 1e-8)) {
        stop("Ortholog UFBoot ", field, " values must be finite integers")
      }
      ortholog_ufboot[[field]] <- as.integer(values)
    }
    ortholog_ufboot$decisive_branch_ufboot <- suppressWarnings(
      as.numeric(ortholog_ufboot$decisive_branch_ufboot)
    )
    if (
      any(!nzchar(ortholog_ufboot$family_id)) ||
        any(!ortholog_ufboot$family_id %in% names(family_order_by_id))
    ) {
      stop("Ortholog UFBoot family_id values must match the ortholog column table")
    }
    if (
      any(!nzchar(ortholog_ufboot$species)) ||
        any(!nzchar(ortholog_ufboot$reference_species)) ||
        any(ortholog_ufboot$reference_species != reference_species) ||
        any(!nzchar(ortholog_ufboot$reference_cds_fasta_id)) ||
        any(!nzchar(ortholog_ufboot$reference_gene_id)) ||
        any(!nzchar(ortholog_ufboot$candidate_cds_fasta_id))
    ) {
      stop("Ortholog UFBoot species and gene identifiers must be non-empty and reference-consistent")
    }
    expected_ufboot_family_order <- unname(
      family_order_by_id[ortholog_ufboot$family_id]
    )
    if (any(ortholog_ufboot$family_order != expected_ufboot_family_order)) {
      stop("Ortholog UFBoot family_order values disagree with the ortholog column table")
    }
    expected_ufboot_family_id <- unname(
      family_id_by_column[as.character(ortholog_ufboot$column_order)]
    )
    if (
      any(is.na(expected_ufboot_family_id)) ||
        any(ortholog_ufboot$family_id != expected_ufboot_family_id)
    ) {
      stop("Ortholog UFBoot column_order values cross family boundaries")
    }
    expected_ufboot_reference_species <- unname(
      reference_species_by_column[as.character(ortholog_ufboot$column_order)]
    )
    expected_ufboot_reference_cds_id <- unname(
      reference_cds_id_by_column[as.character(ortholog_ufboot$column_order)]
    )
    expected_ufboot_reference_gene_id <- unname(
      reference_gene_id_by_column[as.character(ortholog_ufboot$column_order)]
    )
    if (
      any(ortholog_ufboot$reference_species != expected_ufboot_reference_species) ||
        any(ortholog_ufboot$reference_cds_fasta_id != expected_ufboot_reference_cds_id) ||
        any(ortholog_ufboot$reference_gene_id != expected_ufboot_reference_gene_id)
    ) {
      stop("Ortholog UFBoot reference identifiers disagree with the ortholog column table")
    }
    if (
      any(ortholog_ufboot$column_order < ortholog_ufboot$glyph_start_order) ||
        any(ortholog_ufboot$column_order > ortholog_ufboot$glyph_end_order) ||
        any(ortholog_ufboot$glyph_copy_number < 1) ||
        any(ortholog_ufboot$glyph_lane_index < 1) ||
        any(ortholog_ufboot$glyph_lane_count < ortholog_ufboot$glyph_lane_index)
    ) {
      stop("Ortholog UFBoot rows contain invalid glyph spans, copy numbers, or lanes")
    }
    allowed_ufboot_status <- c("evaluated", "not_evaluable", "reference_self")
    if (any(!ortholog_ufboot$orthology_ufboot_status %in% allowed_ufboot_status)) {
      stop("Ortholog UFBoot table contains an invalid orthology_ufboot_status")
    }
    allowed_ufboot_sources <- c("support_generax_ufboot", "support_unrooted")
    if (any(!ortholog_ufboot$ufboot_support_source %in% allowed_ufboot_sources)) {
      stop("Ortholog UFBoot table contains an invalid ufboot_support_source")
    }
    is_reference_self <- ortholog_ufboot$candidate_cds_fasta_id ==
      ortholog_ufboot$reference_cds_fasta_id
    if (
      any(is_reference_self & ortholog_ufboot$orthology_ufboot_status != "reference_self") ||
        any(!is_reference_self & ortholog_ufboot$orthology_ufboot_status == "reference_self")
    ) {
      stop("Ortholog UFBoot reference_self status disagrees with candidate/reference IDs")
    }
    evaluated_ufboot_rows <- ortholog_ufboot$orthology_ufboot_status == "evaluated"
    nonself_ufboot_rows <- !is_reference_self
    nonself_mrca_branch_ids <- suppressWarnings(as.numeric(
      ortholog_ufboot$orthology_mrca_branch_id[nonself_ufboot_rows]
    ))
    if (
      any(nonself_ufboot_rows & !nzchar(ortholog_ufboot$orthology_mrca_branch_id)) ||
        any(!is.finite(nonself_mrca_branch_ids)) ||
        any(abs(nonself_mrca_branch_ids - round(nonself_mrca_branch_ids)) > 1e-8) ||
        any(is_reference_self & nzchar(ortholog_ufboot$orthology_mrca_branch_id)) ||
        any(nonself_ufboot_rows & ortholog_ufboot$orthology_mrca_event != "S") ||
        any(is_reference_self & nzchar(ortholog_ufboot$orthology_mrca_event))
    ) {
      stop("Ortholog UFBoot MRCA fields are inconsistent with the orthology assignment")
    }
    if (
      any(!is.finite(ortholog_ufboot$decisive_branch_ufboot[evaluated_ufboot_rows])) ||
        any(ortholog_ufboot$decisive_branch_ufboot[evaluated_ufboot_rows] < 0) ||
        any(ortholog_ufboot$decisive_branch_ufboot[evaluated_ufboot_rows] > 100) ||
        any(is.finite(ortholog_ufboot$decisive_branch_ufboot[!evaluated_ufboot_rows]))
    ) {
      stop("Ortholog UFBoot values must be 0-100 only for evaluated rows")
    }
    if (
      any(evaluated_ufboot_rows & nzchar(ortholog_ufboot$orthology_ufboot_unavailable_reason)) ||
        any(!evaluated_ufboot_rows & !nzchar(ortholog_ufboot$orthology_ufboot_unavailable_reason))
    ) {
      stop("Ortholog UFBoot availability reasons disagree with evaluation status")
    }
    if (
      any(
        ortholog_ufboot$orthology_ufboot_status == "not_evaluable" &
          !ortholog_ufboot$orthology_ufboot_unavailable_reason %in% c(
            "missing_support", "mrca_is_root"
          )
      ) ||
        any(
          ortholog_ufboot$orthology_ufboot_status == "reference_self" &
            ortholog_ufboot$orthology_ufboot_unavailable_reason != "reference_self"
        )
    ) {
      stop("Ortholog UFBoot table contains an invalid unavailability reason")
    }
    ufboot_pair_keys <- paste(
      ortholog_ufboot$family_id,
      ortholog_ufboot$reference_cds_fasta_id,
      ortholog_ufboot$candidate_cds_fasta_id,
      sep = "::"
    )
    if (anyDuplicated(ufboot_pair_keys)) {
      stop("Ortholog UFBoot table contains duplicate candidate/reference pairs within a family")
    }
    glyph_keys <- paste(
      ortholog_glyphs$family_id,
      ortholog_glyphs$species,
      ortholog_glyphs$start_order,
      ortholog_glyphs$end_order,
      ortholog_glyphs$lane_index,
      ortholog_glyphs$lane_count,
      ortholog_glyphs$copy_number,
      sep = "::"
    )
    ufboot_glyph_keys <- paste(
      ortholog_ufboot$family_id,
      ortholog_ufboot$species,
      ortholog_ufboot$glyph_start_order,
      ortholog_ufboot$glyph_end_order,
      ortholog_ufboot$glyph_lane_index,
      ortholog_ufboot$glyph_lane_count,
      ortholog_ufboot$glyph_copy_number,
      sep = "::"
    )
    if (any(!ufboot_glyph_keys %in% glyph_keys)) {
      stop("Ortholog UFBoot rows do not map to an ortholog glyph")
    }
    ufboot_glyph_row_groups <- split(
      seq_len(nrow(ortholog_ufboot)),
      factor(ufboot_glyph_keys, levels = unique(ufboot_glyph_keys))
    )
    for (row_indices in ufboot_glyph_row_groups) {
      glyph_rows <- ortholog_ufboot[row_indices, , drop = FALSE]
      glyph_is_reference_self <-
        glyph_rows$orthology_ufboot_status == "reference_self"
      if (any(glyph_is_reference_self)) {
        if (!all(glyph_is_reference_self)) {
          stop(
            "Ortholog UFBoot rows for one glyph cannot mix reference-self ",
            "and non-self evidence"
          )
        }
        next
      }
      glyph_mrca_branch_ids <- suppressWarnings(as.numeric(
        glyph_rows$orthology_mrca_branch_id
      ))
      if (length(unique(glyph_mrca_branch_ids)) != 1) {
        stop(
          "Ortholog UFBoot rows for one glyph must share one ",
          "orthology-defining speciation branch"
        )
      }
      if (
        length(unique(glyph_rows$relation)) != 1 ||
          length(unique(glyph_rows$ufboot_support_source)) != 1 ||
          length(unique(glyph_rows$orthology_ufboot_status)) != 1 ||
          length(unique(glyph_rows$orthology_ufboot_unavailable_reason)) != 1
      ) {
        stop(
          "Ortholog UFBoot rows for one glyph must share one availability state"
        )
      }
      if (
        glyph_rows$orthology_ufboot_status[[1]] == "evaluated" &&
          length(unique(glyph_rows$decisive_branch_ufboot)) != 1
      ) {
        stop(
          "Ortholog UFBoot rows for one glyph must share one ",
          "orthology-defining branch UFBoot value"
        )
      }
    }
  }

  if (nrow(ortholog_tree_nodes) > 0) {
    ortholog_tree_nodes$family_id <- as.character(ortholog_tree_nodes$family_id)
    ortholog_tree_nodes$family_order <- suppressWarnings(
      as.numeric(ortholog_tree_nodes$family_order)
    )
    ortholog_tree_nodes$node_id <- as.character(ortholog_tree_nodes$node_id)
    ortholog_tree_nodes$event <- as.character(ortholog_tree_nodes$event)
    ortholog_tree_nodes$is_tip <- suppressWarnings(
      as.numeric(ortholog_tree_nodes$is_tip)
    )
    ortholog_tree_nodes$in_reference_tree <- suppressWarnings(
      as.numeric(ortholog_tree_nodes$in_reference_tree)
    )
    if (
      any(is.na(ortholog_tree_nodes$family_id) | !nzchar(ortholog_tree_nodes$family_id)) ||
        any(!ortholog_tree_nodes$family_id %in% names(family_order_by_id))
    ) {
      stop("Ortholog tree family_id values must match the ortholog column table")
    }
    expected_tree_family_order <- unname(
      family_order_by_id[ortholog_tree_nodes$family_id]
    )
    if (
      any(!is.finite(ortholog_tree_nodes$family_order)) ||
        any(ortholog_tree_nodes$family_order != expected_tree_family_order)
    ) {
      stop("Ortholog tree family_order values disagree with the ortholog column table")
    }
    tree_node_keys <- paste(
      ortholog_tree_nodes$family_id,
      as.character(ortholog_tree_nodes$node_id),
      sep = "::"
    )
    if (anyDuplicated(tree_node_keys)) {
      stop("Ortholog tree contains duplicate node_id values within a family")
    }
    if (any(is.na(ortholog_tree_nodes$node_id) | !nzchar(ortholog_tree_nodes$node_id))) {
      stop("Ortholog tree node_id values must be non-empty")
    }
    if (
      any(!is.finite(ortholog_tree_nodes$is_tip)) ||
        any(!ortholog_tree_nodes$is_tip %in% c(0, 1)) ||
        any(!is.finite(ortholog_tree_nodes$in_reference_tree)) ||
        any(!ortholog_tree_nodes$in_reference_tree %in% c(0, 1))
    ) {
      stop("Ortholog tree is_tip and in_reference_tree values must be 0 or 1")
    }
    compact_tree_rows <- ortholog_tree_nodes$in_reference_tree == 1
    compact_parent_ids <- as.character(
      ortholog_tree_nodes$parent_node_id[compact_tree_rows]
    )
    compact_parent_ids[is.na(compact_parent_ids)] <- ""
    compact_family_ids <- ortholog_tree_nodes$family_id[compact_tree_rows]
    compact_node_keys <- tree_node_keys[compact_tree_rows]
    compact_parent_keys <- paste(compact_family_ids, compact_parent_ids, sep = "::")
    missing_compact_parents <- nzchar(compact_parent_ids) &
      !compact_parent_keys %in% compact_node_keys
    if (any(missing_compact_parents)) {
      stop("Compact ortholog tree contains a parent_node_id absent from its family")
    }
    compact_root_count <- table(compact_family_ids[!nzchar(compact_parent_ids)])
    if (
      !setequal(names(compact_root_count), family_pairs$family_id) ||
        any(compact_root_count != 1)
    ) {
      stop("Compact ortholog tree must contain exactly one root for every family")
    }
    compact_node_height <- suppressWarnings(
      as.numeric(ortholog_tree_nodes$node_height[compact_tree_rows])
    )
    compact_plot_order <- suppressWarnings(
      as.numeric(ortholog_tree_nodes$plot_order[compact_tree_rows])
    )
    if (
      any(!is.finite(compact_node_height)) ||
        any(!is.finite(compact_plot_order))
    ) {
      stop("Compact ortholog tree nodes require finite node_height and plot_order values")
    }
  }
  ortholog_family_count <- length(unique(as.character(ortholog_columns$family_id)))
  mapped_species_nodes <- as.character(ortholog_tree_nodes$mapped_species_node)
  has_mapped_duplications <- nrow(ortholog_tree_nodes) > 0 && any(
    as.character(ortholog_tree_nodes$event) == "D" &
      !is.na(mapped_species_nodes) & nzchar(mapped_species_nodes)
  )
  duplication_family_mode <- has_mapped_duplications
  family_colors <- grDevices::hcl.colors(
    max(1, ortholog_family_count),
    palette = "Dark 3"
  )
}

tree_plot <- ggtree::ggtree(tree_for_plot, size = 0.3)
plot_data <- tree_plot$data
plot_data$y <- (plot_data$y - 1) * heatmap_cell_pitch + 1
original_node_labels <- tree_for_species_mapping$node.label
if (is.null(original_node_labels)) {
  original_node_labels <- rep("", tree$Nnode)
}
species_mapping_label_by_node <- stats::setNames(
  c(as.character(tree$tip.label), as.character(original_node_labels)),
  seq_len(length(tree$tip.label) + tree$Nnode)
)
plot_data$species_mapping_label <- unname(
  species_mapping_label_by_node[as.character(plot_data$node)]
)
plot_data$species_mapping_label[is.na(plot_data$species_mapping_label)] <- ""
if (glyph_mode && nrow(ortholog_tree_nodes) > 0) {
  duplication_event_values <- as.character(ortholog_tree_nodes$event)
  duplication_mapping_values <- as.character(ortholog_tree_nodes$mapped_species_node)
  mapped_duplication_labels <- unique(
    duplication_mapping_values[
      !is.na(duplication_event_values) & duplication_event_values == "D" &
        !is.na(duplication_mapping_values) & nzchar(duplication_mapping_values)
    ]
  )
  nonempty_species_mapping_labels <- plot_data$species_mapping_label[
    nzchar(plot_data$species_mapping_label)
  ]
  species_mapping_counts <- table(nonempty_species_mapping_labels)
  ambiguous_mapped_labels <- intersect(
    names(species_mapping_counts[species_mapping_counts > 1]),
    mapped_duplication_labels
  )
  if (length(ambiguous_mapped_labels) > 0) {
    stop(
      "Mapped duplication label occurs more than once in the species tree: ",
      paste(ambiguous_mapped_labels, collapse = ", ")
    )
  }
}
tip_df <- subset(plot_data, isTip)
tip_levels <- as.character(tip_df$label[order(tip_df$y, decreasing = FALSE)])
tip_levels <- tip_levels[tip_levels %in% unique(as.character(source_df$species))]
if (length(tip_levels) == 0) {
  stop("No species overlap between tree tips and long table.")
}
tip_y_by_species <- stats::setNames(tip_df$y, as.character(tip_df$label))
y_values <- unname(tip_y_by_species[tip_levels])
y_limits <- c(min(y_values, na.rm = TRUE) - 0.55, max(y_values, na.rm = TRUE) + 0.55)

if (glyph_mode) {
  query_levels <- paste0("query_column_", as.numeric(ortholog_columns$column_order))
  query_plot_labels <- as.character(ortholog_columns$plot_label)
  df <- expand.grid(species = tip_levels, query = query_levels, stringsAsFactors = FALSE)
  display_lookup <- unique(source_df[, c("species", "species_display"), drop = FALSE])
  display_lookup <- display_lookup[!duplicated(as.character(display_lookup$species)), , drop = FALSE]
  display_by_species <- stats::setNames(as.character(display_lookup$species_display), as.character(display_lookup$species))
  df$species_display <- unname(display_by_species[df$species])
  missing_display <- is.na(df$species_display) | !nzchar(df$species_display)
  df$species_display[missing_display] <- gg_species_display_from_key(df$species[missing_display])
} else {
  query_levels <- unique(as.character(df$query[order(df$query_order, df$query)]))
  query_plot_labels <- query_levels
  df <- df[df$species %in% tip_levels & df$query %in% query_levels, , drop = FALSE]
}
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

root_stem_needed <- FALSE
root_x_raw <- NA_real_
root_y <- NA_real_
root_stem_start_raw <- NA_real_
root_mapped_family_count <- 1
if (length(root_node) == 1) {
  root_row <- plot_data[plot_data$node == root_node[[1]], , drop = FALSE]
  root_mapping_label <- as.character(root_row$species_mapping_label[[1]])
  mapped_duplication_nodes <- character(0)
  if (glyph_mode && nrow(ortholog_tree_nodes) > 0) {
    mapped_duplication_nodes <- as.character(
      ortholog_tree_nodes$mapped_species_node[
        as.character(ortholog_tree_nodes$event) == "D"
      ]
    )
    mapped_duplication_nodes <- mapped_duplication_nodes[
      !is.na(mapped_duplication_nodes) & nzchar(mapped_duplication_nodes)
    ]
  }
  root_stem_needed <- nzchar(root_mapping_label) &&
    root_mapping_label %in% mapped_duplication_nodes
  if (root_stem_needed) {
    root_event_values <- as.character(ortholog_tree_nodes$event)
    root_mapping_values <- as.character(ortholog_tree_nodes$mapped_species_node)
    root_duplication_rows <- !is.na(root_event_values) & root_event_values == "D" &
      !is.na(root_mapping_values) & root_mapping_values == root_mapping_label
    root_mapped_family_count <- max(
      1,
      length(unique(as.character(ortholog_tree_nodes$family_id[root_duplication_rows])))
    )
  }
  root_x_raw <- as.numeric(root_row$x[[1]])
  root_y <- as.numeric(root_row$y[[1]])
}

x_left <- 0
if (nrow(ci_plot_df) > 0) {
  x_left <- min(0, min(ci_plot_df$xmin, na.rm = TRUE))
}
x_right <- max(plot_data$x, na.rm = TRUE)
if (root_stem_needed) {
  root_stem_fraction <- max(0.08, 0.055 * root_mapped_family_count)
  root_stem_width <- max(
    (x_right - x_left) * root_stem_fraction,
    tree_span * root_stem_fraction,
    0.002
  )
  root_stem_start_raw <- root_x_raw - root_stem_width
  x_left <- min(x_left, root_stem_start_raw)
}
x_pad <- max((x_right - x_left) * 0.025, tree_height * 0.01, 0.001)
x_left <- x_left - x_pad
x_right <- x_right + x_pad

dated_ruler_mode <- tree_scale_mode == "ruler"
if (dated_ruler_mode && !ape::is.ultrametric(tree_for_plot)) {
  stop("--tree_scale=ruler requires an ultrametric species tree")
}
tree_axis_label <- "Branch length"
branch_scale_start <- NA_real_
branch_scale_end <- NA_real_
branch_scale_label <- character(0)
if (dated_ruler_mode) {
  tree_axis_label <- "Million years ago"
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
} else {
  tree_x_min <- min(plot_data$x, na.rm = TRUE)
  tree_x_max <- max(plot_data$x, na.rm = TRUE)
  branch_scale_width <- nice_scale_width(tree_x_max - tree_x_min)
  branch_scale_start <- tree_x_min
  branch_scale_end <- min(tree_x_max, branch_scale_start + branch_scale_width)
  branch_scale_label <- format_scale_bar_label(branch_scale_end - branch_scale_start)
  ruler_breaks <- numeric(0)
  ruler_labels <- character(0)
}
ruler_order <- order(ruler_breaks)
ruler_breaks <- ruler_breaks[ruler_order]
ruler_labels <- ruler_labels[ruler_order]

label_df <- unique(df[, c("species", "species_display"), drop = FALSE])
label_df$species <- factor(as.character(label_df$species), levels = tip_levels)
label_df$y <- unname(tip_y_by_species[as.character(label_df$species)])
label_df <- label_df[order(label_df$y), , drop = FALSE]

if (glyph_mode) {
  df$plot_value <- NA
  df$tile_fill <- "#ffffff"
} else if (value_mode == "copy_number") {
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
max_query_width <- max(nchar(as.character(query_plot_labels), type = "width"))
busco_df <- read_busco_summary_table(busco_table_path, tip_y_by_species)

tree_width <- 10.5
tree_left <- 0
tree_map <- function(x) tree_left + (x - x_left) / (x_right - x_left) * tree_width
tree_right <- tree_left + tree_width
root_stem_df <- data.frame(
  x = numeric(0),
  xend = numeric(0),
  y = numeric(0)
)
if (root_stem_needed) {
  root_stem_df <- data.frame(
    x = tree_map(root_stem_start_raw),
    xend = tree_map(root_x_raw),
    y = root_y,
    stringsAsFactors = FALSE
  )
}
label_right <- tree_right + max(4.0, max_label_width * 0.33)
species_label_matrix_gap <- 0.18
heatmap_left <- label_right + species_label_matrix_gap
heatmap_right <- heatmap_left + length(query_levels) * heatmap_cell_pitch
busco_left <- heatmap_right + ifelse(nrow(busco_df) > 0, 0.95, 0.15)
busco_width <- ifelse(nrow(busco_df) > 0, 3.6, 0)
busco_right <- busco_left + busco_width
x_max <- busco_right + ifelse(nrow(busco_df) > 0, 0.55, 0.15)
x_min <- tree_left - 0.25
layout_data_units_per_inch <- (x_max - x_min) / plot_width

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
branch_scale_df <- data.frame()
if (!dated_ruler_mode) {
  branch_scale_df <- data.frame(
    xmin = tree_map(branch_scale_start),
    xmax = tree_map(branch_scale_end),
    xmid = mean(c(tree_map(branch_scale_start), tree_map(branch_scale_end))),
    label = branch_scale_label,
    stringsAsFactors = FALSE
  )
}
tree_axis_label_x <- mean(c(tree_left, tree_right))
tree_axis_label_hjust <- 0.5
if (!dated_ruler_mode) {
  tree_axis_label_x <- branch_scale_df$xmin[[1]]
  tree_axis_label_hjust <- 0
}

heatmap_df <- df
heatmap_df$x <- heatmap_left + (heatmap_df$query_index - 0.5) * heatmap_cell_pitch
heatmap_cell_half <- heatmap_cell_pitch * 0.46
heatmap_df$xmin <- heatmap_df$x - heatmap_cell_half
heatmap_df$xmax <- heatmap_df$x + heatmap_cell_half
heatmap_df$ymin <- heatmap_df$y - heatmap_cell_half
heatmap_df$ymax <- heatmap_df$y + heatmap_cell_half
query_label_df <- data.frame(
  query = query_plot_labels,
  x = heatmap_left + (seq_along(query_levels) - 0.5) * heatmap_cell_pitch,
  stringsAsFactors = FALSE
)
family_boundary_df <- data.frame()
if (glyph_mode && length(unique(as.character(ortholog_columns$family_id))) > 1) {
  family_last_columns <- tapply(
    as.numeric(ortholog_columns$column_order),
    as.character(ortholog_columns$family_id),
    max
  )
  family_last_columns <- sort(as.numeric(family_last_columns))
  family_last_columns <- family_last_columns[family_last_columns < length(query_levels)]
  family_boundary_df <- data.frame(
    x = heatmap_left + family_last_columns * heatmap_cell_pitch,
    ymin = min(y_values, na.rm = TRUE) - heatmap_cell_half,
    ymax = max(y_values, na.rm = TRUE) + heatmap_cell_half,
    stringsAsFactors = FALSE
  )
}

glyph_rect_df <- data.frame()
if (glyph_mode && nrow(ortholog_glyphs) > 0) {
  glyph_rect_df <- ortholog_glyphs[as.character(ortholog_glyphs$species) %in% tip_levels, , drop = FALSE]
  glyph_rect_df$start_order <- as.numeric(glyph_rect_df$start_order)
  glyph_rect_df$end_order <- as.numeric(glyph_rect_df$end_order)
  glyph_rect_df$lane_index <- pmax(1, as.numeric(glyph_rect_df$lane_index))
  glyph_rect_df$lane_count <- pmax(glyph_rect_df$lane_index, as.numeric(glyph_rect_df$lane_count))
  glyph_rect_df$copy_number <- as.numeric(glyph_rect_df$copy_number)
  glyph_rect_df$y <- unname(tip_y_by_species[as.character(glyph_rect_df$species)])
  lane_height <- (2 * heatmap_cell_half) / glyph_rect_df$lane_count
  glyph_rect_df$xmin <- heatmap_left + (glyph_rect_df$start_order - 0.5) * heatmap_cell_pitch - heatmap_cell_half
  glyph_rect_df$xmax <- heatmap_left + (glyph_rect_df$end_order - 0.5) * heatmap_cell_pitch + heatmap_cell_half
  glyph_rect_df$ymin <- glyph_rect_df$y - heatmap_cell_half + (glyph_rect_df$lane_index - 1) * lane_height + 0.015
  glyph_rect_df$ymax <- glyph_rect_df$y - heatmap_cell_half + glyph_rect_df$lane_index * lane_height - 0.015
  glyph_rect_df$x <- (glyph_rect_df$xmin + glyph_rect_df$xmax) / 2
  glyph_rect_df$text_y <- (glyph_rect_df$ymin + glyph_rect_df$ymax) / 2
  glyph_rect_df$fill <- ifelse(
    glyph_rect_df$relation == "specific",
    "#2166ac",
    ifelse(glyph_rect_df$relation == "shared_ancestral", "#9ecae1", "#fdd49e")
  )
  glyph_rect_df$border <- ifelse(
    glyph_rect_df$relation == "specific",
    "#2166ac",
    ifelse(glyph_rect_df$relation == "shared_ancestral", "#6baed6", "#d95f0e")
  )
  glyph_rect_df$text_color <- ifelse(glyph_rect_df$relation == "specific", "white", "black")
  glyph_rect_df$copy_label <- ifelse(is.finite(glyph_rect_df$copy_number), format(glyph_rect_df$copy_number, trim = TRUE), "")
  evidence_rail_width <- heatmap_cell_half * 0.56 - 0.014
  if (identical(evidence_layout, "rail")) {
    glyph_rect_df$x <- glyph_rect_df$x - evidence_rail_width / 2
  }
}

evidence_group_fields <- c(
  "species", "family_id", "family_order", "relation", "column_order",
  "glyph_copy_number", "glyph_start_order", "glyph_end_order",
  "glyph_lane_index", "glyph_lane_count"
)
add_evidence_geometry <- function(
  evidence_df,
  evidence_half = c("top", "bottom"),
  span_glyph = FALSE
) {
  evidence_half <- match.arg(evidence_half)
  if (nrow(evidence_df) == 0) {
    return(evidence_df)
  }
  evidence_df$y <- unname(tip_y_by_species[as.character(evidence_df$species)])
  lane_height <- (2 * heatmap_cell_half) /
    as.numeric(evidence_df$glyph_lane_count)
  lane_ymin <- evidence_df$y - heatmap_cell_half +
    (as.numeric(evidence_df$glyph_lane_index) - 1) * lane_height + 0.015
  lane_ymax <- evidence_df$y - heatmap_cell_half +
    as.numeric(evidence_df$glyph_lane_index) * lane_height - 0.015
  evidence_df$cell_x <- heatmap_left +
    (as.numeric(evidence_df$column_order) - 0.5) * heatmap_cell_pitch
  if (identical(evidence_layout, "band")) {
    band_height <- (lane_ymax - lane_ymin) * 0.18
    if (identical(evidence_half, "top")) {
      evidence_df$evidence_ymin <- lane_ymax - band_height
      evidence_df$evidence_ymax <- lane_ymax
    } else {
      evidence_df$evidence_ymin <- lane_ymin
      evidence_df$evidence_ymax <- lane_ymin + band_height
    }
    glyph_xmin <- heatmap_left +
      (as.numeric(evidence_df$glyph_start_order) - 0.5) * heatmap_cell_pitch -
      heatmap_cell_half
    glyph_xmax <- heatmap_left +
      (as.numeric(evidence_df$glyph_end_order) - 0.5) * heatmap_cell_pitch +
      heatmap_cell_half
    logical_cell_xmin <- heatmap_left +
      (as.numeric(evidence_df$column_order) - 1) * heatmap_cell_pitch
    logical_cell_xmax <- heatmap_left +
      as.numeric(evidence_df$column_order) * heatmap_cell_pitch
    if (isTRUE(span_glyph)) {
      evidence_df$evidence_xmin <- glyph_xmin + 0.014
      evidence_df$evidence_xmax <- glyph_xmax - 0.014
    } else {
      evidence_df$evidence_xmin <- ifelse(
        evidence_df$column_order == evidence_df$glyph_start_order,
        glyph_xmin + 0.014,
        logical_cell_xmin
      )
      evidence_df$evidence_xmax <- ifelse(
        evidence_df$column_order == evidence_df$glyph_end_order,
        glyph_xmax - 0.014,
        logical_cell_xmax
      )
    }
  } else {
    lane_mid <- (lane_ymin + lane_ymax) / 2
    rail_gap <- pmin(0.006, lane_height * 0.03)
    if (identical(evidence_half, "top")) {
      evidence_df$evidence_ymin <- lane_mid + rail_gap
      evidence_df$evidence_ymax <- lane_ymax
    } else {
      evidence_df$evidence_ymin <- lane_ymin
      evidence_df$evidence_ymax <- lane_mid - rail_gap
    }
    evidence_df$evidence_xmin <- evidence_df$cell_x + heatmap_cell_half * 0.44
    evidence_df$evidence_xmax <- evidence_df$cell_x + heatmap_cell_half - 0.014
  }
  evidence_width <- evidence_df$evidence_xmax - evidence_df$evidence_xmin
  evidence_height <- evidence_df$evidence_ymax - evidence_df$evidence_ymin
  segment_y_inset <- pmin(
    0.014,
    evidence_height * 0.18
  )
  is_unavailable <- evidence_df$evidence_status == "unavailable"
  if (identical(evidence_layout, "band")) {
    segment_half_width <- pmin(0.060, evidence_width * 0.20)
    segment_mid_x <- (evidence_df$evidence_xmin + evidence_df$evidence_xmax) / 2
    evidence_df$status_x <- segment_mid_x - segment_half_width
    evidence_df$status_xend <- segment_mid_x + segment_half_width
  } else {
    segment_x_inset <- pmin(0.018, evidence_width * 0.16)
    evidence_df$status_x <- evidence_df$evidence_xmin + segment_x_inset
    evidence_df$status_xend <- evidence_df$evidence_xmax - segment_x_inset
  }
  evidence_df$status_y <- ifelse(
    is_unavailable,
    evidence_df$evidence_ymin + segment_y_inset,
    (evidence_df$evidence_ymin + evidence_df$evidence_ymax) / 2
  )
  evidence_df$status_yend <- ifelse(
    is_unavailable,
    evidence_df$evidence_ymax - segment_y_inset,
    (evidence_df$evidence_ymin + evidence_df$evidence_ymax) / 2
  )
  evidence_df
}

synteny_evidence_df <- data.frame()
synteny_marker_df <- data.frame()
synteny_legend_visible <- FALSE
synteny_anchor_scale_max <- 1L
synteny_window_label <- ""
if (glyph_mode && nrow(ortholog_synteny) > 0) {
  synteny_plot_rows <- ortholog_synteny[
    as.character(ortholog_synteny$species) %in% tip_levels,
    ,
    drop = FALSE
  ]
  if (nrow(synteny_plot_rows) > 0) {
    evaluated_synteny_status <- c("supported", "single_anchor", "no_support")
    synteny_window_radii <- sort(unique(as.integer(
      synteny_plot_rows$synteny_window_radius[
        synteny_plot_rows$synteny_status %in% evaluated_synteny_status &
          is.finite(synteny_plot_rows$synteny_window_radius)
      ]
    )))
    if (length(synteny_window_radii) == 1) {
      synteny_window_label <- paste0(
        "Window: ", synteny_window_radii[[1]], " upstream + ",
        synteny_window_radii[[1]], " downstream gene models"
      )
    } else if (length(synteny_window_radii) > 1) {
      synteny_window_label <- paste0(
        "Window: ", min(synteny_window_radii), "-", max(synteny_window_radii),
        " gene models per side (family-specific)"
      )
    }
    synteny_plot_rows$pair_count <- 1L
    synteny_plot_rows$evaluated_count <- as.integer(
      synteny_plot_rows$synteny_status %in% evaluated_synteny_status
    )
    synteny_plot_rows$reference_self_count <- as.integer(
      synteny_plot_rows$synteny_status == "reference_self"
    )
    synteny_plot_rows$anchor_for_max <- ifelse(
      synteny_plot_rows$synteny_status %in% evaluated_synteny_status,
      synteny_plot_rows$shared_anchor_count,
      -Inf
    )
    synteny_count_df <- stats::aggregate(
      synteny_plot_rows[
        , c("pair_count", "evaluated_count", "reference_self_count"),
        drop = FALSE
      ],
      by = synteny_plot_rows[, evidence_group_fields, drop = FALSE],
      FUN = sum
    )
    synteny_anchor_df <- stats::aggregate(
      synteny_plot_rows[, "anchor_for_max", drop = FALSE],
      by = synteny_plot_rows[, evidence_group_fields, drop = FALSE],
      FUN = max
    )
    colnames(synteny_anchor_df)[colnames(synteny_anchor_df) == "anchor_for_max"] <-
      "anchor_count"
    synteny_evidence_df <- merge(
      synteny_count_df,
      synteny_anchor_df,
      by = evidence_group_fields,
      all = TRUE,
      sort = FALSE
    )
    if (any(synteny_evidence_df$pair_count != synteny_evidence_df$glyph_copy_number)) {
      stop("Ortholog synteny pair counts must equal the corresponding glyph copy number")
    }
    synteny_evidence_df$evidence_status <- ifelse(
      synteny_evidence_df$reference_self_count == synteny_evidence_df$pair_count,
      "reference_self",
      ifelse(
        synteny_evidence_df$evaluated_count == synteny_evidence_df$pair_count,
        "evaluated",
        "unavailable"
      )
    )
    synteny_evidence_df$anchor_count[
      synteny_evidence_df$evidence_status != "evaluated"
    ] <- NA_real_
    evaluated_anchor_count <- synteny_evidence_df$anchor_count[
      synteny_evidence_df$evidence_status == "evaluated" &
        is.finite(synteny_evidence_df$anchor_count)
    ]
    synteny_anchor_scale_max <- if (length(evaluated_anchor_count) > 0) {
      max(1L, as.integer(max(evaluated_anchor_count)))
    } else {
      1L
    }
    synteny_evidence_df$marker_fill <- ifelse(
      synteny_evidence_df$evidence_status == "reference_self",
      "#ffffff",
      ifelse(
        synteny_evidence_df$evidence_status == "unavailable",
        "#e5e7eb",
        viridis_anchor_rail_colors(
          synteny_evidence_df$anchor_count,
          synteny_anchor_scale_max
        )
      )
    )
    synteny_evidence_df <- add_evidence_geometry(
      synteny_evidence_df,
      "top"
    )
    synteny_legend_visible <- any(
      synteny_evidence_df$evidence_status == "evaluated" &
        is.finite(synteny_evidence_df$anchor_count)
    )
    synteny_marker_df <- synteny_evidence_df[
      synteny_evidence_df$evidence_status == "evaluated" &
        is.finite(synteny_evidence_df$anchor_count) &
        synteny_evidence_df$anchor_count > 0,
      ,
      drop = FALSE
    ]
    if (nrow(synteny_marker_df) > 0) {
      synteny_marker_df$marker_y <- synteny_marker_df$evidence_ymax -
        pmin(0.020, synteny_marker_df$evidence_ymax - synteny_marker_df$evidence_ymin)
      synteny_marker_df$marker_x <- synteny_marker_df$cell_x +
        heatmap_cell_half * 0.78
      synteny_marker_df$marker_fill <- viridis_anchor_colors(
        synteny_marker_df$anchor_count,
        synteny_anchor_scale_max
      )
    }
  }
}

ufboot_evidence_df <- data.frame()
ufboot_marker_df <- data.frame()
ufboot_legend_visible <- FALSE
if (glyph_mode && nrow(ortholog_ufboot) > 0) {
  ufboot_plot_rows <- ortholog_ufboot[
    as.character(ortholog_ufboot$species) %in% tip_levels,
    ,
    drop = FALSE
  ]
  if (nrow(ufboot_plot_rows) > 0) {
    ufboot_glyph_group_fields <- setdiff(
      evidence_group_fields,
      "column_order"
    )
    ufboot_plot_glyph_keys <- do.call(
      paste,
      c(ufboot_plot_rows[, ufboot_glyph_group_fields, drop = FALSE], sep = "::")
    )
    ufboot_plot_row_groups <- split(
      seq_len(nrow(ufboot_plot_rows)),
      factor(ufboot_plot_glyph_keys, levels = unique(ufboot_plot_glyph_keys))
    )
    ufboot_evidence_rows <- lapply(ufboot_plot_row_groups, function(row_indices) {
      glyph_rows <- ufboot_plot_rows[row_indices, , drop = FALSE]
      candidate_count <- length(unique(glyph_rows$candidate_cds_fasta_id))
      reference_count <- length(unique(glyph_rows$reference_cds_fasta_id))
      expected_pair_count <- candidate_count * reference_count
      if (nrow(glyph_rows) != expected_pair_count) {
        stop(
          "Ortholog UFBoot candidate/reference rows must form a complete ",
          "pair set within each glyph"
        )
      }
      if (candidate_count != glyph_rows$glyph_copy_number[[1]]) {
        stop(
          "Ortholog UFBoot candidate counts must equal the corresponding ",
          "glyph copy number"
        )
      }

      summary_row <- glyph_rows[1, ufboot_glyph_group_fields, drop = FALSE]
      summary_row$column_order <- as.integer(summary_row$glyph_end_order)
      summary_row$pair_count <- nrow(glyph_rows)
      summary_row$evaluated_count <- sum(
        glyph_rows$orthology_ufboot_status == "evaluated"
      )
      summary_row$reference_self_count <- sum(
        glyph_rows$orthology_ufboot_status == "reference_self"
      )
      if (summary_row$reference_self_count == summary_row$pair_count) {
        summary_row$evidence_status <- "reference_self"
        summary_row$ufboot <- NA_real_
      } else if (summary_row$evaluated_count == summary_row$pair_count) {
        summary_row$evidence_status <- "evaluated"
        summary_row$ufboot <- glyph_rows$decisive_branch_ufboot[[1]]
      } else {
        summary_row$evidence_status <- "unavailable"
        summary_row$ufboot <- NA_real_
      }
      summary_row
    })
    ufboot_evidence_df <- do.call(rbind, ufboot_evidence_rows)
    rownames(ufboot_evidence_df) <- NULL
    ufboot_evidence_df$marker_fill <- ifelse(
      ufboot_evidence_df$evidence_status == "reference_self",
      "#ffffff",
      ifelse(
        ufboot_evidence_df$evidence_status == "unavailable",
        "#e5e7eb",
        viridis_ufboot_colors(ufboot_evidence_df$ufboot)
      )
    )
    ufboot_evidence_df <- add_evidence_geometry(
      ufboot_evidence_df,
      "bottom",
      span_glyph = TRUE
    )
    ufboot_legend_visible <- any(
      ufboot_evidence_df$evidence_status == "evaluated" &
        is.finite(ufboot_evidence_df$ufboot)
    )
    ufboot_marker_df <- ufboot_evidence_df[
      ufboot_evidence_df$evidence_status == "evaluated" &
        is.finite(ufboot_evidence_df$ufboot),
      ,
      drop = FALSE
    ]
    if (nrow(ufboot_marker_df) > 0) {
      ufboot_marker_df$marker_y <- ufboot_marker_df$evidence_ymin +
        pmin(0.020, ufboot_marker_df$evidence_ymax - ufboot_marker_df$evidence_ymin)
      ufboot_marker_df$marker_x <- ufboot_marker_df$cell_x +
        heatmap_cell_half * 0.78
      ufboot_marker_df$marker_fill <- inferno_ufboot_colors(ufboot_marker_df$ufboot)
    }
  }
}

query_tree_nodes_df <- data.frame()
query_tree_edges_df <- data.frame()
query_tree_title_df <- data.frame()
if (glyph_mode && nrow(ortholog_tree_nodes) > 0) {
  query_tree_nodes_df <- ortholog_tree_nodes[
    suppressWarnings(as.numeric(ortholog_tree_nodes$in_reference_tree)) == 1,
    ,
    drop = FALSE
  ]
  query_tree_nodes_df$family_id <- as.character(query_tree_nodes_df$family_id)
  query_tree_nodes_df$node_id <- as.character(query_tree_nodes_df$node_id)
  query_tree_nodes_df$parent_node_id <- as.character(query_tree_nodes_df$parent_node_id)
  query_tree_nodes_df$parent_node_id[is.na(query_tree_nodes_df$parent_node_id)] <- ""
  query_tree_nodes_df$is_tip <- as.numeric(query_tree_nodes_df$is_tip)
  query_tree_nodes_df$family_order <- as.numeric(query_tree_nodes_df$family_order)
  query_tree_nodes_df$node_height <- as.numeric(query_tree_nodes_df$node_height)
  query_tree_nodes_df$plot_order <- as.numeric(query_tree_nodes_df$plot_order)
  query_tree_nodes_df$duplication_index <- suppressWarnings(
    as.numeric(query_tree_nodes_df$duplication_index)
  )
  query_tree_nodes_df$mapped_species_node <- as.character(
    query_tree_nodes_df$mapped_species_node
  )
  query_tree_nodes_df$mapped_species_node[is.na(query_tree_nodes_df$mapped_species_node)] <- ""
  query_tree_nodes_df$node_key <- paste(query_tree_nodes_df$family_id, query_tree_nodes_df$node_id, sep = "::")
  query_tree_nodes_df$parent_key <- ifelse(
    nzchar(query_tree_nodes_df$parent_node_id),
    paste(query_tree_nodes_df$family_id, query_tree_nodes_df$parent_node_id, sep = "::"),
    ""
  )
  query_tree_nodes_df$event_label <- ifelse(
    query_tree_nodes_df$is_tip == 0,
    as.character(query_tree_nodes_df$event),
    ""
  )
  duplication_rows <- which(
    query_tree_nodes_df$is_tip == 0 & query_tree_nodes_df$event == "D"
  )
  if (length(duplication_rows) > 0) {
    duplication_rows <- duplication_rows[order(
      query_tree_nodes_df$family_order[duplication_rows],
      query_tree_nodes_df$duplication_index[duplication_rows],
      query_tree_nodes_df$node_id[duplication_rows]
    )]
    query_tree_nodes_df$event_label[duplication_rows] <- paste0(
      "D", seq_along(duplication_rows)
    )
  }
  query_tree_nodes_df$x <- heatmap_left + (query_tree_nodes_df$plot_order - 0.5) * heatmap_cell_pitch
  family_max_height <- ave(query_tree_nodes_df$node_height, query_tree_nodes_df$family_id, FUN = max)
  normalized_height <- ifelse(family_max_height > 0, query_tree_nodes_df$node_height / family_max_height, 0)
  query_tree_tip_y <- max(y_values, na.rm = TRUE) + heatmap_cell_half + 0.24
  query_tree_nodes_df$y <- query_tree_tip_y + normalized_height * 2.0
  query_tree_nodes_df$family_color <- family_colors[
    ((query_tree_nodes_df$family_order - 1) %% length(family_colors)) + 1
  ]
  duplication_node_fill <- rep("#d95f02", nrow(query_tree_nodes_df))
  if (duplication_family_mode) {
    duplication_node_fill <- query_tree_nodes_df$family_color
  }
  query_tree_nodes_df$node_fill <- ifelse(
    query_tree_nodes_df$event == "D",
    duplication_node_fill,
    ifelse(query_tree_nodes_df$event == "S", "#1b9e77", "#ffffff")
  )
  query_tree_nodes_df$node_label <- ifelse(
    duplication_family_mode & query_tree_nodes_df$event == "D",
    "",
    query_tree_nodes_df$event_label
  )
  query_tree_nodes_df$text_color <- ifelse(query_tree_nodes_df$event %in% c("D", "S"), "white", "black")

  parent_positions <- query_tree_nodes_df[, c("node_key", "x", "y"), drop = FALSE]
  colnames(parent_positions) <- c("parent_key", "parent_x", "parent_y")
  query_tree_edges_df <- merge(
    query_tree_nodes_df[nzchar(query_tree_nodes_df$parent_key), , drop = FALSE],
    parent_positions,
    by = "parent_key",
    all.x = FALSE,
    all.y = FALSE
  )
  root_nodes <- query_tree_nodes_df[!nzchar(query_tree_nodes_df$parent_key), , drop = FALSE]
  query_tree_title_df <- root_nodes
  single_family_plot <- length(unique(query_tree_nodes_df$family_id)) == 1
  query_tree_title_df$label <- if (single_family_plot) {
    paste0(reference_species_display, " gene tree")
  } else {
    query_tree_title_df$family_id
  }
  family_column_count <- table(as.character(ortholog_columns$family_id))
  family_column_center <- tapply(
    as.numeric(ortholog_columns$column_order),
    as.character(ortholog_columns$family_id),
    mean
  )
  query_tree_title_df$x <- heatmap_left +
    (unname(family_column_center[query_tree_title_df$family_id]) - 0.5) * heatmap_cell_pitch
  query_tree_title_df$title_angle <- ifelse(
    single_family_plot,
    0,
    90
  )
  query_tree_title_df$title_hjust <- ifelse(query_tree_title_df$title_angle == 90, 0, 0.5)
  query_tree_title_df$title_y <- max(query_tree_nodes_df$y, na.rm = TRUE) + 0.48
  query_tree_title_df$title_top_y <- query_tree_title_df$title_y + ifelse(
    query_tree_title_df$title_angle == 90,
    pmax(
      0.65,
      nchar(query_tree_title_df$label, type = "width") * 0.115 * layout_data_units_per_inch + 0.30
    ),
    0.35
  )
}

duplication_map_df <- data.frame()
duplication_family_key_df <- data.frame()
duplication_count_breaks <- c(0, 1)
duplication_count_scale_max <- 1
duplication_bar_height_max <- heatmap_cell_pitch * 0.62
if (glyph_mode && nrow(ortholog_tree_nodes) > 0) {
  all_mapped_species_nodes <- as.character(ortholog_tree_nodes$mapped_species_node)
  all_mapped_species_nodes[is.na(all_mapped_species_nodes)] <- ""
  duplication_event_df <- ortholog_tree_nodes[
    as.character(ortholog_tree_nodes$event) == "D" & nzchar(all_mapped_species_nodes),
    c("family_id", "family_order", "mapped_species_node"),
    drop = FALSE
  ]
  duplication_event_df$family_id <- as.character(duplication_event_df$family_id)
  duplication_event_df$family_order <- as.numeric(duplication_event_df$family_order)
  duplication_event_df$mapped_species_node <- as.character(
    duplication_event_df$mapped_species_node
  )
  duplication_event_df$duplication_count <- 1L
  species_duplication_df <- data.frame()
  if (nrow(duplication_event_df) > 0) {
    species_duplication_df <- stats::aggregate(
      duplication_count ~ family_id + family_order + mapped_species_node,
      data = duplication_event_df,
      FUN = sum
    )
  }
  species_node_df <- plot_data[
    nzchar(plot_data$species_mapping_label),
    c("node", "parent", "x", "y", "species_mapping_label"),
    drop = FALSE
  ]
  colnames(species_node_df)[colnames(species_node_df) == "x"] <- "species_x"
  colnames(species_node_df)[colnames(species_node_df) == "y"] <- "species_y"
  species_parent_df <- plot_data[, c("node", "x"), drop = FALSE]
  colnames(species_parent_df) <- c("parent", "species_parent_x")
  species_node_df <- merge(
    species_node_df,
    species_parent_df,
    by = "parent",
    all.x = TRUE,
    all.y = FALSE,
    sort = FALSE
  )
  if (nrow(species_duplication_df) > 0) {
    missing_species_nodes <- setdiff(
      unique(species_duplication_df$mapped_species_node),
      unique(species_node_df$species_mapping_label)
    )
    if (length(missing_species_nodes) > 0) {
      stop(
        "Mapped duplication species nodes were not found in the species mapping tree: ",
        paste(missing_species_nodes, collapse = ", ")
      )
    }
    duplication_map_df <- merge(
      species_duplication_df,
      species_node_df,
      by.x = "mapped_species_node",
      by.y = "species_mapping_label",
      all.x = FALSE,
      all.y = FALSE,
      sort = FALSE
    )
  }
  if (nrow(duplication_map_df) > 0) {
    root_mapping <- duplication_map_df$node == duplication_map_df$parent
    duplication_map_df$branch_start_raw <- ifelse(
      root_mapping,
      root_stem_start_raw,
      duplication_map_df$species_parent_x
    )
    duplication_map_df$branch_end_raw <- ifelse(
      root_mapping,
      root_x_raw,
      duplication_map_df$species_x
    )
    invalid_branch <- !is.finite(duplication_map_df$branch_start_raw) |
      !is.finite(duplication_map_df$branch_end_raw)
    duplication_map_df$branch_start_raw[invalid_branch] <-
      duplication_map_df$species_x[invalid_branch]
    duplication_map_df$branch_end_raw[invalid_branch] <-
      duplication_map_df$species_x[invalid_branch]
    duplication_map_df$marker_y <- duplication_map_df$species_y
    duplication_map_df$bar_width <- 0.10
    duplication_map_df$branch_start <- tree_map(duplication_map_df$branch_start_raw)
    duplication_map_df$branch_end <- tree_map(duplication_map_df$branch_end_raw)
    duplication_map_df$marker_x <-
      duplication_map_df$branch_start +
      (duplication_map_df$branch_end - duplication_map_df$branch_start) * 0.58
    minimum_bar_center_spacing <- 0.16
    collision_key <- duplication_map_df$mapped_species_node
    for (collision_rows in split(seq_len(nrow(duplication_map_df)), collision_key)) {
      collision_rows <- collision_rows[order(
        duplication_map_df$family_order[collision_rows]
      )]
      count <- length(collision_rows)
      branch_start <- duplication_map_df$branch_start[collision_rows[[1]]]
      branch_end <- duplication_map_df$branch_end[collision_rows[[1]]]
      preferred_centers <- if (count == 1) {
        branch_start + (branch_end - branch_start) * 0.58
      } else {
        seq(
          branch_start + (branch_end - branch_start) * 0.16,
          branch_start + (branch_end - branch_start) * 0.84,
          length.out = count
        )
      }
      preferred_spacing <- if (count > 1) {
        min(abs(diff(preferred_centers)))
      } else {
        Inf
      }
      if (!is.finite(preferred_spacing) || preferred_spacing < minimum_bar_center_spacing) {
        group_center <- branch_start + (branch_end - branch_start) * 0.58
        group_half_extent <-
          (count - 1) / 2 * minimum_bar_center_spacing + 0.06
        if (group_half_extent <= tree_width / 2) {
          group_center <- min(
            max(group_center, tree_left + group_half_extent),
            tree_right - group_half_extent
          )
        }
        preferred_centers <- group_center +
          (seq_len(count) - (count + 1) / 2) * minimum_bar_center_spacing
        preferred_spacing <- minimum_bar_center_spacing
      }
      duplication_map_df$marker_x[collision_rows] <- preferred_centers
      duplication_map_df$bar_width[collision_rows] <- min(
        0.28,
        max(0.10, preferred_spacing * 0.70)
      )
    }
    duplication_map_df$family_color <- family_colors[
      ((as.numeric(duplication_map_df$family_order) - 1) %% length(family_colors)) + 1
    ]
    duplication_count_breaks <- nice_count_breaks(
      max(duplication_map_df$duplication_count, na.rm = TRUE)
    )
    duplication_count_scale_max <- max(duplication_count_breaks)
    duplication_map_df$xmin <- duplication_map_df$marker_x - duplication_map_df$bar_width / 2
    duplication_map_df$xmax <- duplication_map_df$marker_x + duplication_map_df$bar_width / 2
    duplication_map_df$ymin <- duplication_map_df$marker_y
    duplication_map_df$ymax <- duplication_map_df$marker_y +
      duplication_map_df$duplication_count / duplication_count_scale_max *
      duplication_bar_height_max
    x_min <- min(x_min, min(duplication_map_df$xmin, na.rm = TRUE) - 0.05)
    x_max <- max(x_max, max(duplication_map_df$xmax, na.rm = TRUE) + 0.05)
  }
  if (duplication_family_mode && nrow(duplication_map_df) > 0) {
    family_headers <- unique(
      ortholog_columns[, c("family_id", "family_order"), drop = FALSE]
    )
    family_headers$family_order <- as.numeric(family_headers$family_order)
    family_headers <- family_headers[order(family_headers$family_order), , drop = FALSE]
    family_headers$family_color <- family_colors[
      ((family_headers$family_order - 1) %% length(family_colors)) + 1
    ]
    duplication_family_key_df <- data.frame(
      family_order = family_headers$family_order,
      family_id = as.character(family_headers$family_id),
      label = as.character(family_headers$family_id),
      family_color = family_headers$family_color,
      stringsAsFactors = FALSE
    )
  }
}

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
y_axis_label <- ifelse(dated_ruler_mode, row_y_min - 2.35, row_y_min - 1.95)
query_label_gap <- 0.06
y_query_label <- row_y_min - heatmap_cell_half - query_label_gap
query_label_depth <- max(
  2.6,
  max_query_width * 0.115 * layout_data_units_per_inch + 0.40
)
y_legend <- min(row_y_min - 4.15, y_query_label - query_label_depth - 0.80)
y_busco_axis <- row_y_min - 0.85
y_busco_label <- row_y_min - 2.35
y_min <- row_y_min - 6.70
y_max <- row_y_max + 0.75
duplication_family_key_title_df <- data.frame()
duplication_count_key_df <- data.frame()
duplication_count_key_title_df <- data.frame()
if (nrow(duplication_family_key_df) > 0) {
  family_key_label_width <- max(
    nchar(duplication_family_key_df$label, type = "width") *
      0.115 * layout_data_units_per_inch,
    na.rm = TRUE
  )
  family_key_two_column_width <- 0.42 + family_key_label_width + 0.35
  family_key_column_count <- if (
    nrow(duplication_family_key_df) > 6 &&
      family_key_two_column_width <= tree_width / 2
  ) {
    2
  } else {
    1
  }
  family_key_rows_per_column <- ceiling(
    nrow(duplication_family_key_df) / family_key_column_count
  )
  family_key_index <- seq_len(nrow(duplication_family_key_df)) - 1
  duplication_family_key_df$key_column <- floor(
    family_key_index / family_key_rows_per_column
  )
  duplication_family_key_df$key_row <- family_key_index %% family_key_rows_per_column
  family_key_top_y <- y_legend - ifelse(nrow(support_df) > 0, 1.25, 0)
  duplication_family_key_df$x <- tree_left +
    duplication_family_key_df$key_column * (tree_width / family_key_column_count)
  duplication_family_key_df$y <- family_key_top_y -
    duplication_family_key_df$key_row * 0.68
  duplication_family_key_title_df <- data.frame(
    x = tree_left,
    y = family_key_top_y + 0.75,
    label = "Duplication bar colors",
    stringsAsFactors = FALSE
  )
  count_key_values <- duplication_count_breaks[duplication_count_breaks > 0]
  count_key_title_y <- min(duplication_family_key_df$y, na.rm = TRUE) - 0.90
  count_key_y <- count_key_title_y - 0.85
  count_key_spacing <- max(1.35, tree_width / max(2, length(count_key_values)))
  duplication_count_key_df <- data.frame(
    x = tree_left + (seq_along(count_key_values) - 1) * count_key_spacing,
    y = count_key_y,
    duplication_count = count_key_values,
    label = as.character(count_key_values),
    stringsAsFactors = FALSE
  )
  duplication_count_key_df$xmin <- duplication_count_key_df$x
  duplication_count_key_df$xmax <- duplication_count_key_df$x + 0.24
  duplication_count_key_df$ymin <- duplication_count_key_df$y
  duplication_count_key_df$ymax <- duplication_count_key_df$y +
    duplication_count_key_df$duplication_count / duplication_count_scale_max *
    duplication_bar_height_max
  duplication_count_key_title_df <- data.frame(
    x = tree_left,
    y = count_key_title_y,
    label = "Bar height = duplication count",
    stringsAsFactors = FALSE
  )
  y_min <- min(y_min, count_key_y - 0.75)
}
if (nrow(query_tree_title_df) > 0) {
  y_max <- max(y_max, max(query_tree_title_df$title_top_y, na.rm = TRUE) + 0.20)
}

legend_df <- data.frame()
if (glyph_mode) {
  tree_legend_labels <- character(0)
  tree_legend_fills <- character(0)
  if (
    nrow(query_tree_nodes_df) > 0 &&
      any(query_tree_nodes_df$event == "D") &&
      !duplication_family_mode
  ) {
    tree_legend_labels <- c(tree_legend_labels, "D#: mapped duplication")
    tree_legend_fills <- c(tree_legend_fills, "#d95f02")
  }
  if (nrow(query_tree_nodes_df) > 0 && any(query_tree_nodes_df$event == "S")) {
    tree_legend_labels <- c(tree_legend_labels, "S: speciation node")
    tree_legend_fills <- c(tree_legend_fills, "#1b9e77")
  }
  glyph_legend_labels <- c(
    "undetected",
    "reference-gene-specific",
    "pre-duplication copy"
  )
  glyph_legend_fills <- c("#ffffff", "#2166ac", "#9ecae1")
  if (nrow(glyph_rect_df) > 0 && any(glyph_rect_df$relation == "ambiguous")) {
    glyph_legend_labels <- c(glyph_legend_labels, "non-contiguous orthology")
    glyph_legend_fills <- c(glyph_legend_fills, "#fdd49e")
  }
  legend_df <- data.frame(
    label = c(glyph_legend_labels, tree_legend_labels),
    fill = c(glyph_legend_fills, tree_legend_fills),
    x = heatmap_left + 0.05,
    y = y_legend - seq(
      0,
      by = 0.75,
      length.out = length(glyph_legend_labels) + length(tree_legend_labels)
    ),
    stringsAsFactors = FALSE
  )
} else if (value_mode == "presence") {
  legend_df <- data.frame(
    label = c("undetected", "detected"),
    fill = unname(ortholog_colors[c("undetected", "detected")]),
    x = heatmap_left + c(2.95, 7.35),
    stringsAsFactors = FALSE
  )
}
synteny_legend_df <- data.frame()
synteny_legend_title_df <- data.frame()
if (glyph_mode && synteny_legend_visible && !identical(evidence_layout, "off")) {
  synteny_legend_marker_y <- if (nrow(legend_df) > 0 && "y" %in% colnames(legend_df)) {
    min(legend_df$y, na.rm = TRUE) - 2.55
  } else {
    y_legend - 1.30
  }
  synteny_anchor_breaks <- nice_anchor_breaks(synteny_anchor_scale_max)
  if (evidence_strip_mode) {
    synteny_anchor_breaks <- sort(unique(c(0L, synteny_anchor_breaks)))
  }
  synteny_legend_df <- data.frame(
    anchor_count = synteny_anchor_breaks,
    label = as.character(synteny_anchor_breaks),
    fill = if (evidence_strip_mode) {
      viridis_anchor_rail_colors(
        synteny_anchor_breaks,
        synteny_anchor_scale_max
      )
    } else {
      viridis_anchor_colors(
        synteny_anchor_breaks,
        synteny_anchor_scale_max
      )
    },
    x = heatmap_left + 0.28 +
      (seq_along(synteny_anchor_breaks) - 1) * 0.95,
    y = synteny_legend_marker_y,
    stringsAsFactors = FALSE
  )
  synteny_legend_title_df <- data.frame(
    x = heatmap_left + 0.05,
    y = synteny_legend_marker_y + c(1.80, 1.15, 0.55),
    label = c(
      if (evidence_strip_mode) {
        paste0(evidence_layout_label, " (top): local synteny anchors")
      } else {
        "Local synteny anchors"
      },
      if (evidence_strip_mode) {
        "Fill: highest per-copy anchor count (0 = none)"
      } else {
        "Color: highest per-copy anchor count"
      },
      synteny_window_label
    ),
    stringsAsFactors = FALSE
  )
}
ufboot_legend_df <- data.frame()
ufboot_legend_title_df <- data.frame()
if (glyph_mode && ufboot_legend_visible && !identical(evidence_layout, "off")) {
  ufboot_legend_marker_y <- if (nrow(synteny_legend_df) > 0) {
    min(synteny_legend_df$y, na.rm = TRUE) - 2.55
  } else if (nrow(legend_df) > 0 && "y" %in% colnames(legend_df)) {
    min(legend_df$y, na.rm = TRUE) - 2.55
  } else {
    y_legend - 1.30
  }
  ufboot_breaks <- c(0, 50, 80, 95, 100)
  ufboot_legend_df <- data.frame(
    ufboot = ufboot_breaks,
    label = as.character(ufboot_breaks),
    fill = if (evidence_strip_mode) {
      viridis_ufboot_colors(ufboot_breaks)
    } else {
      inferno_ufboot_colors(ufboot_breaks)
    },
    x = heatmap_left + 0.28 +
      (seq_along(ufboot_breaks) - 1) * 1.15,
    y = ufboot_legend_marker_y,
    stringsAsFactors = FALSE
  )
  ufboot_legend_title_df <- data.frame(
    x = heatmap_left + 0.05,
    y = ufboot_legend_marker_y + c(1.80, 1.15, 0.55),
    label = c(
      paste0(
        if (evidence_strip_mode) paste0(evidence_layout_label, " (bottom): ") else "",
        "Gene tree UFBoot (%)"
      ),
      if (evidence_strip_mode) {
        "Fill: orthology-defining speciation-branch support"
      } else {
        "Color: orthology-defining speciation-branch support"
      },
      if (evidence_strip_mode) {
        "One value spans each ortholog glyph"
      } else {
        "No circle: unavailable or reference self"
      }
    ),
    stringsAsFactors = FALSE
  )
}
evidence_state_legend_df <- data.frame()
evidence_state_swatch_df <- data.frame()
evidence_state_legend_title_df <- data.frame()
observed_evidence_statuses <- unique(c(
  if (nrow(synteny_evidence_df) > 0) {
    as.character(synteny_evidence_df$evidence_status)
  } else {
    character(0)
  },
  if (nrow(ufboot_evidence_df) > 0) {
    as.character(ufboot_evidence_df$evidence_status)
  } else {
    character(0)
  }
))
displayed_evidence_statuses <- c("unavailable", "reference_self")
displayed_evidence_statuses <- displayed_evidence_statuses[
  displayed_evidence_statuses %in% observed_evidence_statuses
]
if (
  glyph_mode && evidence_strip_mode &&
    length(displayed_evidence_statuses) > 0
) {
  state_legend_y <- if (nrow(ufboot_legend_df) > 0) {
    min(ufboot_legend_df$y, na.rm = TRUE) - 1.65
  } else if (nrow(synteny_legend_df) > 0) {
    min(synteny_legend_df$y, na.rm = TRUE) - 1.65
  } else if (nrow(legend_df) > 0 && "y" %in% colnames(legend_df)) {
    min(legend_df$y, na.rm = TRUE) - 1.65
  } else {
    y_legend - 1.30
  }
  evidence_state_legend_df <- data.frame(
    evidence_status = displayed_evidence_statuses,
    label = ifelse(
      displayed_evidence_statuses == "reference_self",
      "reference self",
      "unavailable"
    ),
    fill = ifelse(
      displayed_evidence_statuses == "reference_self",
      "#ffffff",
      "#e5e7eb"
    ),
    x = heatmap_left + 0.28,
    y = state_legend_y - seq(
      0,
      by = 0.58,
      length.out = length(displayed_evidence_statuses)
    ),
    stringsAsFactors = FALSE
  )
  evidence_state_legend_df$xmin <- evidence_state_legend_df$x - 0.13
  evidence_state_legend_df$xmax <- evidence_state_legend_df$x + 0.13
  evidence_state_legend_df$ymin <- evidence_state_legend_df$y - 0.20
  evidence_state_legend_df$ymax <- evidence_state_legend_df$y + 0.20
  is_unavailable <- evidence_state_legend_df$evidence_status == "unavailable"
  evidence_state_legend_df$draw_swatch <- rep(
    TRUE,
    nrow(evidence_state_legend_df)
  )
  evidence_state_legend_df$label_x <- evidence_state_legend_df$xmax + 0.16
  evidence_state_legend_df$status_x <- evidence_state_legend_df$xmin + 0.035
  evidence_state_legend_df$status_xend <- evidence_state_legend_df$xmax - 0.035
  evidence_state_legend_df$status_y <- ifelse(
    is_unavailable,
    evidence_state_legend_df$ymin + 0.045,
    evidence_state_legend_df$y
  )
  evidence_state_legend_df$status_yend <- ifelse(
    is_unavailable,
    evidence_state_legend_df$ymax - 0.045,
    evidence_state_legend_df$y
  )
  evidence_state_swatch_df <- evidence_state_legend_df[
    evidence_state_legend_df$draw_swatch,
    ,
    drop = FALSE
  ]
  evidence_state_legend_title_df <- data.frame(
    x = heatmap_left + 0.05,
    y = state_legend_y + 0.65,
    label = paste0(evidence_layout_label, " states"),
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
if (nrow(legend_df) > 0) {
  ortholog_legend_right <- max(
    heatmap_left + 0.05 + nchar(ifelse(glyph_mode, paste0(reference_species_display, " orthologs"), "Ortholog"), type = "width") * 0.33,
    legend_df$x + 0.62 + nchar(legend_df$label, type = "width") * 0.33
  )
  x_max <- max(x_max, ortholog_legend_right + 0.15)
  if ("y" %in% colnames(legend_df) && any(is.finite(legend_df$y))) {
    y_min <- min(y_min, min(legend_df$y, na.rm = TRUE) - 0.75)
  }
}
if (nrow(synteny_legend_df) > 0) {
  synteny_legend_right <- max(
    synteny_legend_df$x + 0.24 +
      nchar(synteny_legend_df$label, type = "width") * 0.33,
    synteny_legend_title_df$x +
      nchar(synteny_legend_title_df$label, type = "width") * 0.33
  )
  x_max <- max(x_max, synteny_legend_right + 0.15)
  y_min <- min(y_min, min(synteny_legend_df$y, na.rm = TRUE) - 0.55)
}
if (nrow(ufboot_legend_df) > 0) {
  ufboot_legend_right <- max(
    ufboot_legend_df$x + 0.24 +
      nchar(ufboot_legend_df$label, type = "width") * 0.33,
    ufboot_legend_title_df$x +
      nchar(ufboot_legend_title_df$label, type = "width") * 0.33
  )
  x_max <- max(x_max, ufboot_legend_right + 0.15)
  y_min <- min(y_min, min(ufboot_legend_df$y, na.rm = TRUE) - 0.55)
}
if (nrow(evidence_state_legend_df) > 0) {
  evidence_state_legend_right <- max(
    evidence_state_legend_df$x + 0.25 +
      nchar(evidence_state_legend_df$label, type = "width") * 0.33,
    evidence_state_legend_title_df$x +
      nchar(evidence_state_legend_title_df$label, type = "width") * 0.33
  )
  x_max <- max(x_max, evidence_state_legend_right + 0.15)
  y_min <- min(y_min, min(evidence_state_legend_df$ymin, na.rm = TRUE) - 0.35)
}

combined <- ggplot()
if (nrow(ci_plot_df) > 0) {
  combined <- combined +
    geom_segment(data = ci_plot_df, aes(x = xmin_plot, xend = xmax_plot, y = y, yend = y), linewidth = 0.8, color = hpd_color, alpha = 0.75, lineend = "butt")
}
combined <- combined +
  geom_segment(data = root_stem_df, aes(x = x, xend = xend, y = y, yend = y), linewidth = 0.28, color = "black", lineend = "square") +
  geom_segment(data = edge_df, aes(x = parent_x, xend = x, y = y, yend = y), linewidth = 0.28, color = "black", lineend = "square") +
  geom_segment(data = edge_df, aes(x = parent_x, xend = parent_x, y = parent_y, yend = y), linewidth = 0.28, color = "black", lineend = "square")
if (nrow(duplication_map_df) > 0 && duplication_family_mode) {
  combined <- combined +
    geom_rect(
      data = duplication_map_df,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        fill = family_color
      ),
      linewidth = 0.16,
      color = "black"
    )
}
if (nrow(duplication_count_key_df) > 0) {
  combined <- combined +
    geom_text(
      data = duplication_count_key_title_df,
      aes(x = x, y = y, label = label),
      hjust = 0,
      size = font_size_mm,
      color = "black"
    ) +
    geom_rect(
      data = duplication_count_key_df,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      fill = "#bdbdbd",
      linewidth = 0.16,
      color = "black"
    ) +
    geom_text(
      data = duplication_count_key_df,
      aes(x = xmax + 0.18, y = (ymin + ymax) / 2, label = label),
      hjust = 0,
      size = font_size_mm * 0.88,
      color = "black"
    )
}
if (nrow(support_df) > 0) {
  combined <- combined +
    geom_text(data = support_df, aes(x = x_plot, y = y + 0.28, label = support_label), size = font_size_mm, color = support_color)
}
combined <- combined +
  geom_text(data = label_df, aes(x = label_right, y = y, label = species_display), hjust = 1, size = font_size_mm, color = "black", fontface = "italic") +
  geom_rect(data = heatmap_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = tile_fill), color = ifelse(glyph_mode, "#d9d9d9", "white"), linewidth = 0.18)
if (nrow(glyph_rect_df) > 0) {
  combined <- combined +
    geom_rect(data = glyph_rect_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill, color = border), linewidth = 0.18)
}
if (evidence_strip_mode && nrow(synteny_evidence_df) > 0) {
  synteny_display_df <- synteny_evidence_df
  synteny_status_df <- synteny_display_df[
    synteny_display_df$evidence_status != "evaluated",
    ,
    drop = FALSE
  ]
  combined <- combined +
    geom_rect(
      data = synteny_display_df,
      aes(
        xmin = evidence_xmin,
        xmax = evidence_xmax,
        ymin = evidence_ymin,
        ymax = evidence_ymax,
        fill = marker_fill
      ),
      linewidth = if (identical(evidence_layout, "band")) 0 else 0.12,
      color = if (identical(evidence_layout, "band")) NA else "#111827"
    )
  if (nrow(synteny_status_df) > 0) {
    combined <- combined +
      geom_segment(
        data = synteny_status_df,
        aes(x = status_x, xend = status_xend, y = status_y, yend = status_yend),
        linewidth = 0.22,
        color = "#4b5563",
        lineend = "round"
      )
  }
}
if (evidence_strip_mode && nrow(ufboot_evidence_df) > 0) {
  ufboot_display_df <- ufboot_evidence_df
  ufboot_status_df <- ufboot_display_df[
    ufboot_display_df$evidence_status != "evaluated",
    ,
    drop = FALSE
  ]
  combined <- combined +
    geom_rect(
      data = ufboot_display_df,
      aes(
        xmin = evidence_xmin,
        xmax = evidence_xmax,
        ymin = evidence_ymin,
        ymax = evidence_ymax,
        fill = marker_fill
      ),
      linewidth = if (identical(evidence_layout, "band")) 0 else 0.12,
      color = if (identical(evidence_layout, "band")) NA else "#111827"
    )
  if (nrow(ufboot_status_df) > 0) {
    combined <- combined +
      geom_segment(
        data = ufboot_status_df,
        aes(x = status_x, xend = status_xend, y = status_y, yend = status_yend),
        linewidth = 0.22,
        color = "#4b5563",
        lineend = "round"
      )
  }
}
if (identical(evidence_layout, "glyph") && nrow(synteny_marker_df) > 0) {
  combined <- combined +
    geom_point(
      data = synteny_marker_df,
      aes(x = marker_x, y = marker_y),
      shape = 23,
      size = 1.25,
      stroke = 0.22,
      fill = "white",
      color = "white"
    ) +
    geom_point(
      data = synteny_marker_df,
      aes(x = marker_x, y = marker_y, fill = marker_fill),
      shape = 23,
      size = 0.92,
      stroke = 0.16,
      color = "#111827"
    )
}
if (identical(evidence_layout, "glyph") && nrow(ufboot_marker_df) > 0) {
  combined <- combined +
    geom_point(
      data = ufboot_marker_df,
      aes(x = marker_x, y = marker_y),
      shape = 21,
      size = 1.15,
      stroke = 0.20,
      fill = "white",
      color = "white"
    ) +
    geom_point(
      data = ufboot_marker_df,
      aes(x = marker_x, y = marker_y, fill = marker_fill),
      shape = 21,
      size = 0.82,
      stroke = 0.15,
      color = "#111827"
    )
}
# Keep copy-number labels above every evidence representation. In particular,
# edge bands occupy part of the presence/absence cell and must never obscure
# the primary numeric value.
if (nrow(glyph_rect_df) > 0) {
  combined <- combined +
    geom_text(
      data = glyph_rect_df,
      aes(x = x, y = text_y, label = copy_label, color = text_color),
      size = font_size_mm,
      fontface = "bold"
    )
}
if (nrow(family_boundary_df) > 0) {
  combined <- combined +
    geom_segment(
      data = family_boundary_df,
      aes(x = x, xend = x, y = ymin, yend = ymax),
      linewidth = 0.45,
      color = "#7f7f7f"
    )
}
if (nrow(query_tree_edges_df) > 0) {
  combined <- combined +
    geom_segment(data = query_tree_edges_df, aes(x = x, xend = x, y = y, yend = parent_y), linewidth = 0.32, color = "black") +
    geom_segment(data = query_tree_edges_df, aes(x = x, xend = parent_x, y = parent_y, yend = parent_y), linewidth = 0.32, color = "black")
}
if (nrow(query_tree_nodes_df) > 0) {
  internal_query_tree_nodes <- query_tree_nodes_df[query_tree_nodes_df$is_tip == 0, , drop = FALSE]
  internal_duplication_nodes <- internal_query_tree_nodes[
    internal_query_tree_nodes$event == "D", , drop = FALSE
  ]
  internal_other_nodes <- internal_query_tree_nodes[
    internal_query_tree_nodes$event != "D", , drop = FALSE
  ]
  if (nrow(internal_other_nodes) > 0) {
    combined <- combined +
      geom_point(data = internal_other_nodes, aes(x = x, y = y, fill = node_fill), shape = 21, size = 2.8, stroke = 0.25, color = "black") +
      geom_text(data = internal_other_nodes, aes(x = x, y = y, label = node_label, color = text_color), size = font_size_mm * 0.62, fontface = "bold")
  }
  if (nrow(internal_duplication_nodes) > 0) {
    if (duplication_family_mode) {
      combined <- combined +
        geom_point(
          data = internal_duplication_nodes,
          aes(x = x, y = y, fill = family_color),
          shape = 21,
          size = 2.0,
          stroke = 0.25,
          color = "black"
        )
    } else {
      combined <- combined +
        geom_label(
        data = internal_duplication_nodes,
        aes(x = x, y = y, label = node_label),
        size = font_size_mm * 0.70,
        fontface = "bold",
        fill = "#d95f02",
        color = "white",
        linewidth = 0.16,
        label.padding = grid::unit(0.08, "lines"),
        label.r = grid::unit(0.08, "lines")
      )
    }
  }
  combined <- combined +
    geom_point(data = query_tree_nodes_df[query_tree_nodes_df$is_tip == 1, , drop = FALSE], aes(x = x, y = y), shape = 21, size = 1.5, stroke = 0.25, fill = "white", color = "black")
}
if (nrow(query_tree_title_df) > 0) {
  combined <- combined +
    geom_text(
      data = query_tree_title_df,
      aes(x = x, y = title_y, label = label, angle = title_angle, hjust = title_hjust),
      vjust = 0.5,
      size = font_size_mm,
      fontface = "plain",
      color = "black"
    )
}
combined <- combined +
  geom_text(data = query_label_df, aes(x = x, y = y_query_label, label = query), angle = 90, hjust = 1, vjust = 0.5, size = font_size_mm, color = "black")
if (dated_ruler_mode) {
  combined <- combined +
    geom_segment(aes(x = tree_map(x_left + x_pad), xend = tree_map(x_right - x_pad), y = y_ruler, yend = y_ruler), linewidth = 0.28, color = "black") +
    geom_segment(data = ruler_df, aes(x = x, xend = x, y = y_ruler, yend = y_ruler - 0.15), linewidth = 0.28, color = "black") +
    geom_text(data = ruler_df, aes(x = x, y = y_ruler - 0.30, label = label), size = font_size_mm, vjust = 1, color = "black")
} else {
  combined <- combined +
    geom_segment(data = branch_scale_df, aes(x = xmin, xend = xmax, y = y_ruler, yend = y_ruler), linewidth = 0.75, color = "black", lineend = "butt") +
    geom_text(data = branch_scale_df, aes(x = xmid, y = y_ruler - 0.22, label = label), size = font_size_mm, vjust = 1, color = "black")
}
combined <- combined +
  annotate("text", x = tree_axis_label_x, y = y_axis_label, label = tree_axis_label, hjust = tree_axis_label_hjust, size = font_size_mm, color = "black") +
  scale_fill_identity() +
  scale_color_identity() +
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

if (glyph_mode) {
  combined <- combined +
    annotate("text", x = heatmap_left, y = y_legend + 0.75, label = paste0(reference_species_display, " orthologs"), hjust = 0, size = font_size_mm, color = "black") +
    geom_rect(data = legend_df, aes(xmin = x, xmax = x + 0.36, ymin = y - 0.18, ymax = y + 0.18, fill = fill), color = "#808080", linewidth = 0.14) +
    geom_text(data = legend_df, aes(x = x + 0.62, y = y, label = label), hjust = 0, size = font_size_mm, color = "black")
} else if (value_mode == "presence") {
  combined <- combined +
    annotate("text", x = heatmap_left + 0.05, y = y_legend, label = "Ortholog", hjust = 0, size = font_size_mm, color = "black") +
    geom_rect(data = legend_df, aes(xmin = x, xmax = x + 0.36, ymin = y_legend - 0.18, ymax = y_legend + 0.18, fill = fill), color = "black", linewidth = 0.14) +
    geom_text(data = legend_df, aes(x = x + 0.62, y = y_legend, label = label), hjust = 0, size = font_size_mm, color = "black")
}
if (nrow(synteny_legend_df) > 0) {
  combined <- combined +
    geom_text(
      data = synteny_legend_title_df,
      aes(x = x, y = y, label = label),
      hjust = 0,
      size = font_size_mm,
      color = "black"
    )
  if (evidence_strip_mode) {
    combined <- combined +
      geom_rect(
        data = synteny_legend_df,
        aes(
          xmin = x - 0.12,
          xmax = x + 0.12,
          ymin = y - 0.18,
          ymax = y + 0.18,
          fill = fill
        ),
        linewidth = 0.16,
        color = "#111827"
      )
  } else {
    combined <- combined +
      geom_point(
        data = synteny_legend_df,
        aes(x = x, y = y, fill = fill),
        shape = 23,
        size = 1.55,
        stroke = 0.22,
        color = "#111827"
      )
  }
  combined <- combined +
    geom_text(
      data = synteny_legend_df,
      aes(x = x + 0.22, y = y, label = label),
      hjust = 0,
      size = font_size_mm,
      color = "black"
    )
}
if (nrow(ufboot_legend_df) > 0) {
  combined <- combined +
    geom_text(
      data = ufboot_legend_title_df,
      aes(x = x, y = y, label = label),
      hjust = 0,
      size = font_size_mm,
      color = "black"
    )
  if (evidence_strip_mode) {
    combined <- combined +
      geom_rect(
        data = ufboot_legend_df,
        aes(
          xmin = x - 0.12,
          xmax = x + 0.12,
          ymin = y - 0.18,
          ymax = y + 0.18,
          fill = fill
        ),
        linewidth = 0.16,
        color = "#111827"
      )
  } else {
    combined <- combined +
      geom_point(
        data = ufboot_legend_df,
        aes(x = x, y = y, fill = fill),
        shape = 21,
        size = 1.45,
        stroke = 0.20,
        color = "#111827"
      )
  }
  combined <- combined +
    geom_text(
      data = ufboot_legend_df,
      aes(x = x + 0.25, y = y, label = label),
      hjust = 0,
      size = font_size_mm * 0.88,
      color = "black"
    )
}
if (nrow(evidence_state_legend_df) > 0) {
  combined <- combined +
    geom_text(
      data = evidence_state_legend_title_df,
      aes(x = x, y = y, label = label),
      hjust = 0,
      size = font_size_mm,
      color = "black"
    ) +
    geom_rect(
      data = evidence_state_swatch_df,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
      linewidth = 0.16,
      color = "#111827"
    ) +
    geom_segment(
      data = evidence_state_swatch_df,
      aes(x = status_x, xend = status_xend, y = status_y, yend = status_yend),
      linewidth = 0.28,
      color = "#4b5563",
      lineend = "round"
    ) +
    geom_text(
      data = evidence_state_legend_df,
      aes(x = label_x, y = y, label = label),
      hjust = 0,
      size = font_size_mm * 0.88,
      color = "black"
    )
}
if (nrow(duplication_family_key_df) > 0) {
  combined <- combined +
    geom_text(
      data = duplication_family_key_title_df,
      aes(x = x, y = y, label = label),
      hjust = 0,
      size = font_size_mm,
      color = "black"
    ) +
    geom_rect(
      data = duplication_family_key_df,
      aes(xmin = x, xmax = x + 0.24, ymin = y - 0.22, ymax = y + 0.22, fill = family_color),
      linewidth = 0.16,
      color = "black"
    ) +
    geom_text(
      data = duplication_family_key_df,
      aes(x = x + 0.42, y = y, label = label),
      hjust = 0,
      size = font_size_mm * 0.88,
      color = "black"
    )
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
  ggsave(
    out_pdf,
    plot = combined,
    width = plot_width,
    height = plot_height,
    units = "in",
    dpi = 300,
    limitsize = FALSE,
    bg = "transparent",
    device = grDevices::cairo_pdf
  )
  cat("Wrote PDF:", out_pdf, "\n")
}
if (nzchar(out_svg)) {
  dir.create(dirname(out_svg), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_svg, plot = combined, width = plot_width, height = plot_height, units = "in", dpi = 300, limitsize = FALSE, bg = "transparent")
  cat("Wrote SVG:", out_svg, "\n")
}
