#!/usr/bin/env Rscript

# orthogroup_copy_number_trait_pgls.r
# Usage:
#   Rscript orthogroup_copy_number_trait_pgls.r \
#     --file_orthogroup_copy_number=orthogroup_copy_number/orthogroup_copy_number.tsv \
#     --file_sptree=dated_species_tree.nwk \
#     --file_trait=species_trait.tsv \
#     --outdir=trait_pgls
#
# Outputs:
#   outdir/orthogroup_copy_number_matrix.tsv
#   outdir/orthogroup_copy_number_trait_pgls.tsv
#   outdir/orthogroup_copy_number_trait_pgls.significant.tsv
#   outdir/orthogroup_copy_number_trait_pgls.summary.(pdf|svg)

cat(as.character(Sys.time()), "Starting orthogroup_copy_number_trait_pgls.r\n")

suppressPackageStartupMessages({
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
source_path <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
support_dir <- if (!is.null(source_path)) dirname(normalizePath(source_path, mustWork = TRUE)) else script_dir

parse_args <- function(argv) {
  out <- list()
  for (arg in argv) {
    if (!startsWith(arg, "--")) {
      stop("Arguments must use --name=value syntax: ", arg)
    }
    body <- substring(arg, 3)
    eq_pos <- regexpr("=", body, fixed = TRUE)[[1]]
    if (eq_pos < 1) {
      out[[body]] <- "1"
    } else {
      key <- substring(body, 1, eq_pos - 1)
      value <- substring(body, eq_pos + 1)
      out[[key]] <- value
    }
  }
  out
}

parse_string <- function(args, key, default = "") {
  value <- args[[key]]
  if (is.null(value) || length(value) == 0L || is.na(value) || value == "") {
    return(default)
  }
  as.character(value[[1]])
}

parse_integer <- function(args, key, default) {
  value <- parse_string(args, key, as.character(default))
  out <- suppressWarnings(as.integer(value))
  if (is.na(out)) {
    stop("Invalid integer for --", key, ": ", value)
  }
  out
}

parse_numeric <- function(args, key, default) {
  value <- parse_string(args, key, as.character(default))
  out <- suppressWarnings(as.numeric(value))
  if (!is.finite(out)) {
    stop("Invalid numeric value for --", key, ": ", value)
  }
  out
}

parse_bool <- function(args, key, default = TRUE) {
  value <- tolower(parse_string(args, key, ifelse(default, "1", "0")))
  if (value %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (value %in% c("0", "false", "f", "no", "n")) return(FALSE)
  stop("Invalid boolean value for --", key, ": ", value)
}

read_tsv_base <- function(path, na = character(), col_classes = NA) {
  fields <- count.fields(path, sep = "\t", quote = "\"", comment.char = "", blank.lines.skip = FALSE)
  # NA entries represent continued lines within quoted fields.
  if (!length(fields) || any(fields[!is.na(fields)] != fields[[1]])) stop("TSV rows have inconsistent field counts: ", path)
  header <- read.delim(path, nrows = 0L, sep = "\t", quote = "\"", comment.char = "", check.names = FALSE)
  labels <- names(header)
  if (anyDuplicated(labels) || any(!nzchar(trimws(labels))) || any(labels != trimws(labels)) || any(grepl("[\t\r\n]", labels))) {
    stop("TSV headers must be unique, non-empty and unpadded: ", path)
  }
  read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "\"",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fill = FALSE,
    colClasses = col_classes,
    na.strings = na
  )
}

write_tsv_base <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "NA")
}

safe_as_num <- function(x) {
  out <- suppressWarnings(as.numeric(x))
  if (length(out) == 0L) return(NA_real_)
  out[[1]]
}

normalize_species_label <- function(x) {
  gsub(" ", "_", as.character(x), fixed = TRUE)
}

split_tokens <- function(text) {
  if (is.null(text) || is.na(text) || !nzchar(text)) {
    return(character(0))
  }
  tokens <- unlist(strsplit(as.character(text), "[,[:space:]]+"))
  tokens[nzchar(tokens)]
}

deduplicate_ordered <- function(values) {
  out <- character(0)
  seen <- character(0)
  for (value in values) {
    if (!value %in% seen) {
      out <- c(out, value)
      seen <- c(seen, value)
    }
  }
  out
}

read_family_file <- function(path) {
  if (is.null(path) || is.na(path) || !nzchar(path)) {
    return(character(0))
  }
  if (!file.exists(path) || file.info(path)$size == 0) {
    stop("Family file was not found or empty: ", path)
  }
  raw <- readLines(path, warn = FALSE)
  raw <- raw[nzchar(trimws(raw))]
  raw <- raw[!startsWith(trimws(raw), "#")]
  if (!length(raw)) {
    return(character(0))
  }
  tokens <- vapply(strsplit(raw, "[\t, ]+"), function(x) x[[1]], character(1))
  header_tokens <- c("family_id", "query", "orthogroup", "gene_family_id", "id")
  if (length(tokens) > 0 && tolower(tokens[[1]]) %in% header_tokens) {
    tokens <- tokens[-1]
  }
  tokens[nzchar(tokens)]
}

parse_max_families <- function(value) {
  value <- tolower(trimws(as.character(value)))
  if (value %in% c("", "auto", "all", "none", "no", "unlimited")) {
    return(0L)
  }
  out <- suppressWarnings(as.integer(value))
  if (is.na(out) || out < 0L) {
    stop("--max_families must be all, auto, or a non-negative integer: ", value)
  }
  out
}

select_orthogroup_copy_number_families <- function(available_ids, family_ids = "", family_file = "", max_families = "all") {
  explicit <- deduplicate_ordered(c(split_tokens(family_ids), read_family_file(family_file)))
  available_ids <- as.character(available_ids)
  if (length(explicit)) {
    missing <- explicit[!explicit %in% available_ids]
    if (length(missing)) {
      stop("Requested orthogroup family ID(s) were not found: ", paste(head(missing, 20), collapse = ", "))
    }
    return(explicit)
  }
  max_count <- parse_max_families(max_families)
  if (max_count == 0L) {
    return(available_ids)
  }
  head(available_ids, max_count)
}

detect_family_col <- function(df) {
  candidates <- c("Orthogroup", "FamilyID", "#FamilyID", "HOG", "OG", "family_id", "query")
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit)) {
    return(hit[[1]])
  }
  if (ncol(df) >= 1L) {
    return(colnames(df)[[1]])
  }
  stop("Orthogroup copy-number table has no columns.")
}

assert_unique_labels <- function(labels, context) {
  dup <- unique(labels[duplicated(labels)])
  if (length(dup)) {
    stop(context, " has duplicate labels after normalization: ", paste(head(dup, 20), collapse = ", "))
  }
}

load_tree_normalized <- function(path) {
  tree <- ape::read.tree(path)
  tree$tip.label <- normalize_species_label(tree$tip.label)
  assert_unique_labels(tree$tip.label, "Species tree")
  tree
}

load_trait_table <- function(path) {
  trait <- read_tsv_base(path, col_classes = "character")
  if (ncol(trait) < 2L) {
    stop("Trait table must contain a species column and at least one trait column: ", path)
  }
  if ("species" %in% names(trait)[-1]) stop("Species key conflicts with a species trait column.")
  colnames(trait)[1] <- "species"
  trait$species <- normalize_species_label(trait$species)
  trait <- trait[!is.na(trait$species) & nzchar(trait$species), , drop = FALSE]
  assert_unique_labels(trait$species, "Trait table")
  attr(trait, "source_path") <- normalizePath(path, mustWork = TRUE)
  trait
}

resolve_trait_selection <- function(trait, trait_arg = "all") {
  path <- attr(trait, "source_path")
  if (is.null(path)) {
    cols <- resolve_trait_cols(trait, trait_arg)
    return(data.frame(trait = cols, value_type = "unspecified", status = "selected", reason = ""))
  }
  report_path <- tempfile("trait-selection-", fileext = ".tsv")
  on.exit(unlink(report_path), add = TRUE)
  command <- c(file.path(support_dir, "species_trait_schema.py"), "--table", path,
               "--trait", trait_arg, "--report", report_path)
  log <- suppressWarnings(system2("python", shQuote(command), stdout = TRUE, stderr = TRUE))
  status <- attr(log, "status")
  if (!is.null(status) && status != 0L) stop(paste(c("Trait selection failed:", log), collapse = "\n"))
  read_tsv_base(report_path)
}

resolve_trait_cols <- function(trait, trait_arg = "all") {
  if (!is.null(attr(trait, "source_path"))) {
    report <- resolve_trait_selection(trait, trait_arg)
    return(report$trait[report$status == "selected"])
  }
  if (anyDuplicated(names(trait))) stop("Trait headers must be unique.")
  available <- setdiff(colnames(trait), "species")
  requested <- split_tokens(trait_arg)
  if (anyDuplicated(requested)) stop("Requested traits must be unique.")
  if (!length(requested) || identical(tolower(trait_arg), "all")) {
    return(available)
  }
  missing <- requested[!requested %in% available]
  if (length(missing)) {
    stop("Requested trait column(s) were not found: ", paste(missing, collapse = ", "))
  }
  requested
}

load_orthogroup_copy_number_matrix <- function(file_orthogroup_copy_number, tree, family_ids = "", family_file = "", max_families = "all") {
  copy_number_df <- read_tsv_base(file_orthogroup_copy_number, col_classes = "character")
  family_col <- detect_family_col(copy_number_df)
  family_all <- as.character(copy_number_df[[family_col]])
  if (any(!nzchar(family_all))) {
    stop("Orthogroup copy-number table contains empty family IDs in column: ", family_col)
  }
  if (anyDuplicated(family_all)) {
    duplicated_ids <- unique(family_all[duplicated(family_all)])
    stop("Orthogroup copy-number table contains duplicate family IDs: ", paste(head(duplicated_ids, 20), collapse = ", "))
  }

  raw_cols <- colnames(copy_number_df)
  norm_cols <- normalize_species_label(raw_cols)
  assert_unique_labels(norm_cols, "Orthogroup copy-number table column names")
  col_map <- setNames(raw_cols, norm_cols)
  species <- tree$tip.label
  missing_species <- species[!species %in% names(col_map)]
  if (length(missing_species)) {
    stop("Orthogroup copy-number table is missing species column(s) from the species tree: ", paste(head(missing_species, 20), collapse = ", "))
  }

  selected_families <- select_orthogroup_copy_number_families(family_all, family_ids = family_ids, family_file = family_file, max_families = max_families)
  selected_rows <- match(selected_families, family_all)
  species_raw_cols <- unname(col_map[species])
  count_df <- copy_number_df[selected_rows, species_raw_cols, drop = FALSE]
  count_df[] <- lapply(count_df, function(x) {
    numeric <- suppressWarnings(as.numeric(x))
    if (any(!is.finite(numeric) | numeric < 0 | numeric != floor(numeric))) {
      stop("Predictor copy numbers must be finite non-negative integers.")
    }
    numeric
  })
  count_mat <- t(as.matrix(count_df))
  rownames(count_mat) <- species
  colnames(count_mat) <- selected_families
  storage.mode(count_mat) <- "numeric"
  count_mat
}

copy_matrix_to_table <- function(copy_matrix) {
  data.frame(species = rownames(copy_matrix), as.data.frame(copy_matrix, check.names = FALSE), check.names = FALSE)
}

subset_tree_to_species <- function(tree, species) {
  keep <- unique(as.character(species))
  drop_tips <- tree$tip.label[!tree$tip.label %in% keep]
  if (length(drop_tips)) {
    # Preserve the original root and shared Brownian history after missing tips
    # are removed, including when all retained tips lie in one descendant clade.
    tree <- ape::drop.tip(tree, drop_tips, collapse.singles = FALSE)
  }
  tree
}

make_empty_result_row <- function(family_id, trait_col, n_species = 0L, status = "skipped", skip_reason = "") {
  data.frame(
    Orthogroup = family_id,
    trait = trait_col,
    n_species = as.integer(n_species),
    R2 = NA_real_,
    pval = NA_real_,
    logLik = NA_real_,
    coefficient = NA_real_,
    standard_error = NA_real_,
    confidence_interval_lower = NA_real_,
    confidence_interval_upper = NA_real_,
    confidence_level = 0.95,
    statistic = NA_real_,
    degrees_of_freedom = NA_real_,
    evolutionary_rate = NA_real_,
    response_dispersion = NA_real_,
    response_reference = "",
    separation_warning = "",
    boundary_warning = "",
    optimizer_converged = "",
    response_family = "gaussian",
    link_function = "identity",
    coefficient_penalty = "none",
    predictor_transform = "log1p",
    covariance_estimator = "gaussian-REML",
    evolution_model = "brownian",
    measurement_error_model = "none",
    inference_method = "wald",
    inference_status = "",
    optimizer_message = "",
    small_sample_warning = "",
    PCC = NA_real_,
    OLS_slope = NA_real_,
    p.adj.global = NA_real_,
    p.adj.by_trait = NA_real_,
    fit_mode = status,
    status = status,
    skip_reason = skip_reason,
    error_message = "",
    stringsAsFactors = FALSE
  )
}

empty_orthogroup_copy_number_trait_result <- function() {
  out <- make_empty_result_row("", "")
  out[0, , drop = FALSE]
}

# Brownian GLS with intercept and REML residual rate preserves the old
# Rphylopars coefficient/SE test for one complete observation per species.
# Rphylopars disables phenotype error for those unreplicated inputs.
resolve_response_families <- function(spec, trait_cols) {
  out <- setNames(rep("gaussian", length(trait_cols)), trait_cols)
  seen <- character()
  for (token in split_tokens(spec)) {
    pair <- strsplit(token, "=", fixed = TRUE)[[1]]
    if (length(pair) != 2L || !pair[[1]] %in% trait_cols ||
        !pair[[2]] %in% c("gaussian", "binomial", "poisson", "negative-binomial")) {
      stop("Invalid response family mapping: ", token,
           ". Use trait=gaussian|binomial|poisson|negative-binomial.")
    }
    if (pair[[1]] %in% seen) stop("Duplicate response family: ", pair[[1]])
    seen <- c(seen, pair[[1]])
    out[[pair[[1]]]] <- pair[[2]]
  }
  out
}

parse_response_values <- function(values, family) {
  missing <- is.na(values) | as.character(values) %in% c("", "NA", "NaN", "nan")
  numeric <- suppressWarnings(as.numeric(values))
  if (any(!missing & !is.finite(numeric))) stop("Trait values must be numeric and finite, or explicitly missing.")
  numeric[missing] <- NA_real_
  validate_response_values(numeric, family)
  numeric
}

validate_response_values <- function(values, family) {
  values <- values[is.finite(values)]
  if (family == "binomial" && any(!values %in% c(0, 1))) stop("Binomial traits must contain only 0/1 or missing values.")
  if (family %in% c("poisson", "negative-binomial") && any(values < 0 | values != floor(values))) {
    stop("Count traits must contain non-negative integers or missing values.")
  }
}

fit_nwkit_copy_number_model <- function(model_df, tree, response_family = "gaussian") {
  executable <- Sys.which("nwkit")
  if (!nzchar(executable)) stop("nwkit regress is required for copy-number PGLS.")
  work <- tempfile("nwkit-copy-number-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  tree_path <- file.path(work, "tree.nwk")
  data_path <- file.path(work, "data.tsv")
  result_path <- file.path(work, "result.tsv")
  ape::write.tree(tree, file = tree_path, digits = 17)
  data <- model_df[, c("species", "trait_value", "copy_number"), drop = FALSE]
  names(data)[[1]] <- "leaf_name"
  write_tsv_base(data, data_path)
  command <- c("regress", "--tree", tree_path, "--data", data_path,
               "--responses", "trait_value", "--predictors", "copy_number",
               "--evolution-model", "brownian", "--intercept", "yes",
               "--reml", if (response_family == "gaussian") "yes" else "no",
               "--response-family", paste0("trait_value=", response_family),
               "--coefficient-penalty", "none", "--inference", "wald",
               "--confidence-level", "0.95", "--outfile", result_path)
  if (response_family == "binomial") {
    command <- c(command, "--categorical-responses", "trait_value", "--response-reference", "trait_value=0")
  }
  log <- suppressWarnings(system2(executable, shQuote(command), stdout = TRUE, stderr = TRUE))
  status <- attr(log, "status")
  if (!is.null(status) && status != 0L) {
    stop(paste(c("nwkit regress failed:", log), collapse = "\n"))
  }
  result <- read_tsv_base(result_path, na = c("NA", ""))
  result <- result[result$term == "copy_number", , drop = FALSE]
  if (nrow(result) != 1L) stop("NWKIT must return one copy_number coefficient.")
  result
}

fit_one_orthogroup_copy_number_trait <- function(model_df, tree, family_id, trait_col, min_species = 4L,
                               verbose = FALSE, response_family = "gaussian") {
  model_df$trait_value <- parse_response_values(model_df$trait_value, response_family)
  validate_response_values(model_df$trait_value, response_family)
  model_df$copy_number <- suppressWarnings(as.numeric(model_df$copy_number))
  if (any(!is.finite(model_df$copy_number) | model_df$copy_number < 0 |
          model_df$copy_number != floor(model_df$copy_number))) {
    stop("Predictor copy numbers must be finite non-negative integers.")
  }
  model_df <- model_df[
    !is.na(model_df$species) & nzchar(model_df$species) &
      is.finite(model_df$trait_value) & is.finite(model_df$copy_number),
    ,
    drop = FALSE
  ]
  model_df$copy_number <- log1p(model_df$copy_number)
  n_species <- length(unique(model_df$species))
  if (n_species < min_species) {
    return(make_empty_result_row(family_id, trait_col, n_species, "skipped", "too_few_species"))
  }
  if (length(unique(model_df$trait_value)) < 2L) {
    return(make_empty_result_row(family_id, trait_col, n_species, "skipped", "invariant_trait"))
  }
  if (length(unique(model_df$copy_number)) < 2L) {
    return(make_empty_result_row(family_id, trait_col, n_species, "skipped", "invariant_copy_number"))
  }

  tree_use <- subset_tree_to_species(tree, model_df$species)
  model_df <- model_df[model_df$species %in% tree_use$tip.label, , drop = FALSE]
  model_df <- model_df[match(tree_use$tip.label, model_df$species), , drop = FALSE]
  n_species <- nrow(model_df)
  if (n_species < min_species) {
    return(make_empty_result_row(family_id, trait_col, n_species, "skipped", "too_few_tree_matched_species"))
  }

  pcc <- suppressWarnings(cor(model_df$copy_number, model_df$trait_value, method = "pearson", use = "complete.obs"))
  ols_slope <- suppressWarnings(stats::coef(stats::lm(trait_value ~ copy_number, data = model_df))[["copy_number"]])

  fitted <- tryCatch(fit_nwkit_copy_number_model(model_df, tree_use, response_family), error = identity)
  if (inherits(fitted, "error")) {
    out <- make_empty_result_row(family_id, trait_col, n_species, "error", "fit_failed")
    out$PCC <- safe_as_num(pcc)
    out$OLS_slope <- safe_as_num(ols_slope)
    out$fit_mode <- if (response_family == "gaussian") "nwkit_brownian_reml" else paste0("nwkit_", response_family, "_ml")
    out$error_message <- conditionMessage(fitted)
    return(out)
  }
  usable <- identical(as.character(fitted$inference_status[[1]]), "ok")
  out <- make_empty_result_row(family_id, trait_col, n_species,
                               if (usable) "ok" else "not_estimable",
                               if (usable) "" else "nwkit_inference_unavailable")
  for (field in intersect(names(out), names(fitted))) out[[field]] <- fitted[[field]][[1]]
  out$R2 <- safe_as_num(fitted$r_squared)
  out$pval <- safe_as_num(fitted$p_value)
  out$logLik <- safe_as_num(fitted$log_likelihood)
  out$PCC <- safe_as_num(pcc)
  out$OLS_slope <- safe_as_num(ols_slope)
  out$fit_mode <- if (response_family == "gaussian") "nwkit_brownian_reml" else paste0("nwkit_", response_family, "_ml")
  out
}

run_orthogroup_copy_number_trait_associations <- function(copy_matrix, trait, tree, trait_cols, min_species = 4L,
                                        p_adjust_method = "BH", verbose = FALSE, response_families = "") {
  family_map <- resolve_response_families(response_families, trait_cols)
  rows <- list()
  idx <- 0L
  for (trait_col in trait_cols) {
    if (verbose) {
      cat("Processing trait:", trait_col, "\n")
    }
    trait_values <- parse_response_values(trait[[trait_col]], family_map[[trait_col]])
    validate_response_values(trait_values, family_map[[trait_col]])
    trait_df <- data.frame(species = trait$species, trait_value = trait_values, stringsAsFactors = FALSE)
    for (family_id in colnames(copy_matrix)) {
      idx <- idx + 1L
      copy_df <- data.frame(species = rownames(copy_matrix), copy_number = copy_matrix[, family_id], stringsAsFactors = FALSE)
      model_df <- merge(copy_df, trait_df, by = "species", all = FALSE, sort = FALSE)
      rows[[idx]] <- fit_one_orthogroup_copy_number_trait(
        model_df = model_df,
        tree = tree,
        family_id = family_id,
        trait_col = trait_col,
        min_species = min_species,
        verbose = verbose, response_family = family_map[[trait_col]]
      )
      rows[[idx]]$covariance_estimator <- if (family_map[[trait_col]] == "gaussian") "gaussian-REML" else "laplace-ML"
      rows[[idx]]$response_family <- family_map[[trait_col]]
      rows[[idx]]$link_function <- switch(family_map[[trait_col]], gaussian = "identity", binomial = "logit", "log")
    }
  }
  if (!length(rows)) {
    return(empty_orthogroup_copy_number_trait_result())
  }
  out <- do.call(rbind, rows)
  ok <- out$status == "ok" & is.finite(out$pval)
  if (any(ok)) {
    out$p.adj.global[ok] <- p.adjust(out$pval[ok], method = p_adjust_method)
    for (trait_col in unique(out$trait)) {
      trait_ok <- ok & out$trait == trait_col
      if (any(trait_ok)) {
        out$p.adj.by_trait[trait_ok] <- p.adjust(out$pval[trait_ok], method = p_adjust_method)
      }
    }
  }
  out <- out[order(out$p.adj.global, out$pval, out$trait, out$Orthogroup, na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}

save_summary_plot <- function(df_stat, outdir, alpha = 0.05, top_n = 50L) {
  plot_file_pdf <- file.path(outdir, "orthogroup_copy_number_trait_pgls.summary.pdf")
  plot_file_svg <- file.path(outdir, "orthogroup_copy_number_trait_pgls.summary.svg")
  plot_df <- df_stat[df_stat$status == "ok" & is.finite(df_stat$pval), , drop = FALSE]
  if (nrow(plot_df) == 0L) {
    p <- ggplot() +
      annotate("text", x = 0, y = 0, label = "No fitted orthogroup copy-number trait associations") +
      theme_void()
  } else {
    plot_df <- plot_df[order(plot_df$p.adj.global, plot_df$pval), , drop = FALSE]
    plot_df <- head(plot_df, max(1L, as.integer(top_n)))
    plot_df$label <- paste(plot_df$Orthogroup, plot_df$trait, sep = " / ")
    plot_df$label <- factor(plot_df$label, levels = rev(plot_df$label))
    plot_df$neg_log10_padj <- -log10(pmax(plot_df$p.adj.global, .Machine$double.xmin))
    p <- ggplot(plot_df, aes(x = label, y = neg_log10_padj, fill = trait)) +
      geom_col(width = 0.75) +
      geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "blue", linewidth = 0.25) +
      coord_flip() +
      xlab("") +
      ylab("-log10 adjusted P value") +
      theme_bw() +
      theme(
        text = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        legend.position = "bottom"
      )
  }
  height <- max(4, min(14, 2 + 0.18 * max(1, nrow(plot_df))))
  ggsave(plot_file_pdf, p, width = 7.2, height = height, dpi = 300)
  ggsave(plot_file_svg, p, width = 7.2, height = height, dpi = 300)
}

publish_copy_number_results <- function(work, outdir, inputs) {
  # Use the same NWKIT rollback primitive as the Python workflow adapters.
  code <- paste(c(
    "import pathlib, shutil, sys",
    "from nwkit.file_paths import validate_outputs_do_not_replace_inputs",
    "from nwkit.output_transaction import output_transaction",
    "work, target = map(pathlib.Path, sys.argv[1:3])",
    "sources = sorted(work.iterdir())",
    "outputs = [target / source.name for source in sources]",
    "validate_outputs_do_not_replace_inputs([(str(i), p) for i, p in enumerate(sys.argv[3:])], [(str(i), p) for i, p in enumerate(outputs)])",
    "with output_transaction(outputs, create_parents=True) as staged:",
    "    for source, output in zip(sources, outputs, strict=True):",
    "        shutil.copyfile(source, staged[output])"
  ), collapse = "\n")
  status <- system2("python", shQuote(c("-c", code, work, outdir, inputs)))
  if (status != 0L) stop("Could not publish the complete copy-number result bundle.")
}

run_orthogroup_copy_number_trait_pgls <- function(file_orthogroup_copy_number, file_sptree, file_trait, outdir,
                                trait_arg = "all", min_species = 4L,
                                family_ids = "", family_file = "", max_families = "all",
                                p_adjust_method = "BH", alpha = 0.05, plot_top_n = 50L,
                                verbose = FALSE, response_families = "") {
  if (!nzchar(Sys.which("nwkit"))) stop("nwkit regress is required for copy-number PGLS.")
  work <- tempfile("nwkit-copy-number-results-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  tree <- load_tree_normalized(file_sptree)
  trait <- load_trait_table(file_trait)
  trait_selection <- resolve_trait_selection(trait, trait_arg)
  trait_cols <- trait_selection$trait[trait_selection$status == "selected"]
  write_tsv_base(trait_selection, file.path(work, "trait_selection.tsv"))
  copy_matrix <- load_orthogroup_copy_number_matrix(
    file_orthogroup_copy_number = file_orthogroup_copy_number,
    tree = tree,
    family_ids = family_ids,
    family_file = family_file,
    max_families = max_families
  )

  write_tsv_base(copy_matrix_to_table(copy_matrix), file.path(work, "orthogroup_copy_number_matrix.tsv"))

  df_stat <- run_orthogroup_copy_number_trait_associations(
    copy_matrix = copy_matrix,
    trait = trait,
    tree = tree,
    trait_cols = trait_cols,
    min_species = min_species,
    p_adjust_method = p_adjust_method,
    response_families = response_families,
    verbose = verbose
  )
  write_tsv_base(df_stat, file.path(work, "orthogroup_copy_number_trait_pgls.tsv"))
  df_sig <- df_stat[df_stat$status == "ok" & is.finite(df_stat$p.adj.global) & df_stat$p.adj.global < alpha, , drop = FALSE]
  write_tsv_base(df_sig, file.path(work, "orthogroup_copy_number_trait_pgls.significant.tsv"))
  save_summary_plot(df_stat, outdir = work, alpha = alpha, top_n = plot_top_n)
  publish_copy_number_results(work, outdir, c(file_orthogroup_copy_number, file_sptree, file_trait,
                                            if (file.exists(paste0(file_trait, ".schema.json"))) paste0(file_trait, ".schema.json"),
                                            if (nzchar(family_file)) family_file))
  invisible(list(copy_matrix = copy_matrix, stats = df_stat, significant = df_sig))
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  required <- c("file_orthogroup_copy_number", "file_sptree", "file_trait", "outdir")
  missing <- required[vapply(required, function(key) {
    is.null(args[[key]]) || !nzchar(args[[key]])
  }, logical(1))]
  if (length(missing)) {
    stop("Missing required argument(s): ", paste(missing, collapse = ", "))
  }
  run_orthogroup_copy_number_trait_pgls(
    file_orthogroup_copy_number = parse_string(args, "file_orthogroup_copy_number"),
    file_sptree = parse_string(args, "file_sptree"),
    file_trait = parse_string(args, "file_trait"),
    outdir = parse_string(args, "outdir"),
    trait_arg = parse_string(args, "trait", "all"),
    response_families = parse_string(args, "response_families", ""),
    min_species = parse_integer(args, "min_species", 4L),
    family_ids = parse_string(args, "family_ids", ""),
    family_file = parse_string(args, "family_file", ""),
    max_families = parse_string(args, "max_families", "all"),
    p_adjust_method = parse_string(args, "p_adjust_method", "BH"),
    alpha = parse_numeric(args, "alpha", 0.05),
    plot_top_n = parse_integer(args, "plot_top_n", 50L),
    verbose = parse_bool(args, "verbose", TRUE)
  )
  cat(as.character(Sys.time()), "orthogroup_copy_number_trait_pgls.r completed successfully. Exiting\n")
}

if (!identical(Sys.getenv("GG_ORTHOGROUP_COPY_NUMBER_TRAIT_PGLS_NO_MAIN"), "1")) {
  main()
}
