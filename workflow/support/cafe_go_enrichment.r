# cafe_go_enrichment.r
# Usage:
# Rscript ${dir_myscript}/cafe_go_enrichment.r \
# "${dir_cafe_output}/Gamma_change.tab" \
# "${dir_cafe_output}/Gamma_branch_probabilities.tab" \
# "${file_gene_id}" \
# "${file_go_annotation}" \
# "${dir_go_enrichment}" \
# "${target_branch_go}" \
# "${change_direction_go}" \
# "${go_category}"

cat(as.character(Sys.time()), "Starting cafe_go_enrichment.r\n")

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

write_tsv_base <- function(df, path) {
  write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "NA")
}

resolve_target_branch_id <- function(target_branch, branch_columns) {
  exact_matches <- branch_columns[branch_columns == target_branch]
  if (length(exact_matches) == 1) {
    return(exact_matches)
  }
  if (length(exact_matches) > 1) {
    stop("Multiple exact branches found for target branch: ", paste(exact_matches, collapse = ", "))
  }
  fixed_matches <- grep(target_branch, branch_columns, value = TRUE, fixed = TRUE)
  suffix <- ""
  if (length(fixed_matches) > 0) {
    suffix <- paste0(". Fixed substring candidate(s): ", paste(fixed_matches, collapse = ", "))
  }
  stop("Target branch not found in the CAFE tree: ", target_branch, suffix)
}

split_tokens <- function(x, sep) {
  if (is.na(x) || x == "") {
    return(NA_character_)
  }
  strsplit(as.character(x), sep, fixed = TRUE)[[1]]
}

explode_rows <- function(df, col, sep) {
  if (nrow(df) == 0) {
    return(df)
  }
  pieces <- lapply(df[[col]], split_tokens, sep = sep)
  idx <- rep(seq_len(nrow(df)), lengths(pieces))
  out <- df[idx, , drop = FALSE]
  out[[col]] <- unlist(pieces, use.names = FALSE)
  rownames(out) <- NULL
  out
}

explode_rows_multi <- function(df, cols, sep) {
  if (nrow(df) == 0) {
    return(df)
  }
  out_list <- vector("list", nrow(df))
  for (i in seq_len(nrow(df))) {
    pieces <- lapply(cols, function(col) split_tokens(df[[col]][i], sep = sep))
    max_len <- max(lengths(pieces))
    row_df <- df[rep(i, max_len), , drop = FALSE]
    for (j in seq_along(cols)) {
      vals <- pieces[[j]]
      if (length(vals) < max_len) {
        vals <- c(vals, rep(NA_character_, max_len - length(vals)))
      }
      row_df[[cols[j]]] <- vals
    }
    out_list[[i]] <- row_df
  }
  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  out
}

pivot_longer_base <- function(df, id_col) {
  if (nrow(df) == 0) {
    return(data.frame(FamilyID = character(), branch = character(), is_change = numeric(), stringsAsFactors = FALSE))
  }
  value_cols <- setdiff(colnames(df), id_col)
  long_df <- reshape(
    df,
    direction = "long",
    varying = value_cols,
    v.names = "is_change",
    timevar = "branch",
    times = value_cols
  )
  long_df <- long_df[, c(id_col, "branch", "is_change"), drop = FALSE]
  rownames(long_df) <- NULL
  long_df
}

summarise_orthogroup_targets <- function(df) {
  if (nrow(df) == 0) {
    return(data.frame(
      FamilyID = character(),
      gene_ids = character(),
      sprot_bests = character(),
      sprot_recnames = character(),
      stringsAsFactors = FALSE
    ))
  }
  family_ids <- unique(df$FamilyID)
  out <- data.frame(
    FamilyID = family_ids,
    gene_ids = NA_character_,
    sprot_bests = NA_character_,
    sprot_recnames = NA_character_,
    stringsAsFactors = FALSE
  )
  for (i in seq_along(family_ids)) {
    fam <- family_ids[i]
    tmp <- df[df$FamilyID == fam, , drop = FALSE]
    out$gene_ids[i] <- paste0(tmp$gene_id, collapse = "; ")
    out$sprot_bests[i] <- paste0(tmp$sprot_best, collapse = "; ")
    out$sprot_recnames[i] <- paste0(tmp$sprot_recname, collapse = "; ")
  }
  out
}

summarise_go_enrichment <- function(event_go_df, total_target_events, total_other_events, p_value_threshold) {
  out_cols <- c(
    "go_ids", "go_aspects", "go_terms",
    "n_change_in_target_in_go", "n_change_out_target_in_go",
    "n_change_in_target_out_go", "n_change_out_target_out_go",
    "odds_ratio", "p_value", "p_value_adjusted", "orthogroup_in_target"
  )

  if (nrow(event_go_df) == 0) {
    empty <- as.data.frame(setNames(replicate(length(out_cols), character(), simplify = FALSE), out_cols), stringsAsFactors = FALSE)
    empty$n_change_in_target_in_go <- numeric()
    empty$n_change_out_target_in_go <- numeric()
    empty$n_change_in_target_out_go <- numeric()
    empty$n_change_out_target_out_go <- numeric()
    empty$odds_ratio <- numeric()
    empty$p_value <- numeric()
    empty$p_value_adjusted <- numeric()
    return(list(all = empty, significant = empty))
  }

  key <- paste(event_go_df$go_ids, event_go_df$go_aspects, event_go_df$go_terms, sep = "\r")
  uniq_key <- unique(key)
  rows <- vector("list", length(uniq_key))
  row_idx <- 0L

  for (k in uniq_key) {
    tmp <- event_go_df[key == k, , drop = FALSE]
    n_target <- sum(tmp$is_target, na.rm = TRUE)
    n_other <- sum(!tmp$is_target, na.rm = TRUE)
    if (n_target <= 0) {
      next
    }
    key_parts <- strsplit(k, "\r", fixed = TRUE)[[1]]
    fisher_mat <- matrix(
      c(
        n_target,
        n_other,
        total_target_events - n_target,
        total_other_events - n_other
      ),
      nrow = 2
    )
    fisher_res <- fisher.test(fisher_mat, alternative = "greater")

    row_idx <- row_idx + 1L
    rows[[row_idx]] <- data.frame(
      go_ids = key_parts[1],
      go_aspects = key_parts[2],
      go_terms = key_parts[3],
      n_change_in_target_in_go = n_target,
      n_change_out_target_in_go = n_other,
      n_change_in_target_out_go = total_target_events - n_target,
      n_change_out_target_out_go = total_other_events - n_other,
      odds_ratio = if (!is.null(fisher_res$estimate)) unname(fisher_res$estimate) else NA_real_,
      p_value = fisher_res$p.value,
      orthogroup_in_target = paste0(unique(tmp$FamilyID[tmp$is_target]), collapse = ", "),
      stringsAsFactors = FALSE
    )
  }

  if (row_idx == 0L) {
    empty <- as.data.frame(setNames(replicate(length(out_cols), character(), simplify = FALSE), out_cols), stringsAsFactors = FALSE)
    empty$n_change_in_target_in_go <- numeric()
    empty$n_change_out_target_in_go <- numeric()
    empty$n_change_in_target_out_go <- numeric()
    empty$n_change_out_target_out_go <- numeric()
    empty$odds_ratio <- numeric()
    empty$p_value <- numeric()
    empty$p_value_adjusted <- numeric()
    return(list(all = empty, significant = empty))
  }

  go_enrich_df <- do.call(rbind, rows[seq_len(row_idx)])
  go_enrich_df$p_value_adjusted <- p.adjust(go_enrich_df$p_value, method = "BH")
  go_enrich_df <- go_enrich_df[order(go_enrich_df$p_value_adjusted, go_enrich_df$p_value), , drop = FALSE]
  go_enrich_df <- go_enrich_df[, out_cols, drop = FALSE]

  go_enrich_sig_df <- go_enrich_df[go_enrich_df$p_value_adjusted < p_value_threshold, , drop = FALSE]
  go_enrich_sig_df <- go_enrich_sig_df[order(go_enrich_sig_df$p_value_adjusted), , drop = FALSE]
  list(all = go_enrich_df, significant = go_enrich_sig_df)
}


# CAFE tests a change in turnover rate; the input reconstruction supplies
# increase/decrease. Family BH is computed once, before direction selection.
select_cafe_families <- function(families, alpha, direction) {
  required <- c("FamilyID", "status", "p_value", "target_change", "lambda_target", "lambda_background")
  if (!all(required %in% names(families)) || anyNA(families$FamilyID) ||
      anyDuplicated(families$FamilyID) || any(!nzchar(families$FamilyID))) {
    stop("Malformed native CAFE family results.")
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1 ||
      !direction %in% c("increase", "decrease", "both")) stop("Invalid family alpha or direction.")
  if (anyNA(families$status) || any(families$status != "tested")) {
    stop("CAFE model comparison is incomplete; refusing to redefine the GO background.")
  }
  numeric_cols <- c("p_value", "target_change", "lambda_target", "lambda_background")
  if (!all(vapply(families[numeric_cols], is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(families[numeric_cols]))) ||
      any(families$p_value < 0 | families$p_value > 1) ||
      any(families$lambda_target < 0 | families$lambda_background < 0) ||
      any(families$target_change != trunc(families$target_change))) stop("Invalid native CAFE numerical results.")
  families$p_value_adjusted <- p.adjust(families$p_value, "BH")
  families$rate_shift <- ifelse(families$lambda_target > families$lambda_background, "faster",
                               ifelse(families$lambda_target < families$lambda_background, "slower", "equal"))
  families$selected <- families$p_value_adjusted < alpha & families$rate_shift == "faster" &
    if (direction == "increase") families$target_change > 0 else
      if (direction == "decrease") families$target_change < 0 else families$target_change != 0
  families
}

summarise_specific_go <- function(families, annotations, candidate_go, alpha = 0.05) {
  tested <- families[families$status == "tested", , drop = FALSE]
  selected <- tested$FamilyID[tested$selected]
  annotations <- unique(annotations[annotations$FamilyID %in% tested$FamilyID, , drop = FALSE])
  out <- candidate_go[, c("go_ids", "go_aspects", "go_terms"), drop = FALSE]
  for (col in c("n_specific_in_go", "n_other_in_go", "n_specific_out_go", "n_other_out_go")) out[[col]] <- integer(nrow(out))
  out$odds_ratio <- rep(NA_real_, nrow(out))
  out$p_value <- rep(1, nrow(out))
  out$orthogroup_in_target <- rep("", nrow(out))
  for (i in seq_len(nrow(out))) {
    ids <- unique(annotations$FamilyID[annotations$go_ids == out$go_ids[i]])
    a <- sum(selected %in% ids)
    b <- sum(!tested$selected & tested$FamilyID %in% ids)
    c <- length(selected) - a
    d <- nrow(tested) - length(selected) - b
    out[i, c("n_specific_in_go", "n_other_in_go", "n_specific_out_go", "n_other_out_go")] <- list(a, b, c, d)
    if (length(selected) > 0 && length(selected) < nrow(tested) && a + b > 0 && c + d > 0) {
      fit <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
      out$p_value[i] <- fit$p.value
      out$odds_ratio[i] <- unname(fit$estimate)
    }
    out$orthogroup_in_target[i] <- paste(selected[selected %in% ids], collapse = ", ")
  }
  out$p_value_adjusted <- p.adjust(out$p_value, "BH")
  out <- out[order(out$p_value_adjusted, out$p_value, out$go_ids), , drop = FALSE]
  list(all = out, significant = out[out$p_value_adjusted < alpha, , drop = FALSE])
}

# Input
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  stop("Usage: Rscript cafe_go_enrichment.r Gamma_change.tab Gamma_branch_probabilities.tab gene_id_file go_annotation_file outdir target_branch change_direction go_category [event|cafe_lrt] [family_alpha] [bootstrap_replicates] [fit_restarts] [max_iterations] [cores] [dated_tree]")
}
outdir <- args[5]
target_branch <- args[6]
direction <- args[7]
go_category <- strsplit(args[8], ",")[[1]]

go_method <- if (length(args) >= 9) args[9] else "event"
if (!go_method %in% c("event", "cafe_lrt")) stop("Invalid GO method: ", go_method)
if (go_method == "cafe_lrt") {
  if (!direction %in% c("increase", "decrease", "both") || grepl("[/\\\\]", target_branch)) {
    stop("Invalid native GO direction or target branch.")
  }
  # Invalidate published summaries before reading or validating new inputs.
  # Native run caches remain available for a subsequent corrected request.
  result_names <- c("family_specificity.tsv", "specificity_metadata.tsv",
    paste0("enrichment_significant_", direction, "_", target_branch, c("_all_go.tsv", "_significant_go.tsv")))
  unlink(file.path(outdir, result_names))
}
change_df <- read_tsv_base(args[1])
branch_probabilities_df <- read_tsv_base(args[2], na = c("N/A"))
orthogroup_df <- read_tsv_base(args[3])
ref_annotation_df <- read_tsv_base(args[4])
family_alpha <- if (length(args) >= 10) suppressWarnings(as.numeric(args[10])) else 0.05
bootstrap_replicates <- if (length(args) >= 11) args[11] else "999"
fit_restarts <- if (length(args) >= 12) args[12] else "5"
max_iterations <- if (length(args) >= 13) args[13] else "1000"
native_cores <- if (length(args) >= 14) args[14] else "1"
dated_tree <- if (length(args) >= 15) args[15] else ""
if (go_method == "cafe_lrt" && (length(family_alpha) != 1 || !is.finite(family_alpha) || family_alpha <= 0 || family_alpha >= 1)) {
  stop("family_alpha must be between 0 and 1.")
}
p_value_threshold <- 0.05

target_branch_id <- resolve_target_branch_id(target_branch, colnames(change_df))
go_ref_sp <- sub("\\.annotation\\.tsv$", "", basename(args[4]))
if (!all(go_category %in% c("BP", "MF", "CC"))) {
  stop("Invalid GO category specified. Valid categories are 'BP', 'MF', or 'CC'. Use comma to separate multiple categories.")
}
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

# Identify significant increase/decrease events from CAFE output
if (colnames(branch_probabilities_df)[1] == "#FamilyID") {
  colnames(branch_probabilities_df)[1] <- "FamilyID"
}
branch_probabilities_df <- branch_probabilities_df[, colnames(change_df), drop = FALSE]
original_colnames <- colnames(branch_probabilities_df)
value_cols <- setdiff(original_colnames, "FamilyID")

branch_probabilities_df$.ord <- seq_len(nrow(branch_probabilities_df))
merged_df <- merge(
  branch_probabilities_df,
  change_df,
  by = "FamilyID",
  all.x = TRUE,
  sort = FALSE,
  suffixes = c("_prob", "_change")
)
merged_df <- merged_df[order(merged_df$.ord), , drop = FALSE]

significant_change_df <- data.frame(FamilyID = merged_df$FamilyID, stringsAsFactors = FALSE)
for (col in value_cols) {
  prob_col <- paste0(col, "_prob")
  change_col <- paste0(col, "_change")
  prob_vec <- merged_df[[prob_col]]
  change_vec <- merged_df[[change_col]]
  significant_change_df[[col]] <- ifelse(
    prob_vec >= 0.05,
    0,
    ifelse(change_vec > 0, 1, ifelse(change_vec < 0, -1, 0))
  )
}

keep_cols <- c(
  "FamilyID",
  value_cols[!vapply(value_cols, function(col) all(is.na(significant_change_df[[col]])), logical(1))]
)
significant_change_df <- significant_change_df[, keep_cols, drop = FALSE]
value_cols <- setdiff(colnames(significant_change_df), "FamilyID")

if (direction == "increase") {
  for (col in value_cols) {
    significant_change_df[[col]] <- ifelse(significant_change_df[[col]] == 1, 1, 0)
  }
} else if (direction == "decrease") {
  for (col in value_cols) {
    significant_change_df[[col]] <- ifelse(significant_change_df[[col]] == -1, 1, 0)
  }
} else if (direction == "both" && go_method == "cafe_lrt") {
  for (col in value_cols) significant_change_df[[col]] <- abs(significant_change_df[[col]])
} else {
  stop("Invalid direction. Use 'increase' or 'decrease' (or 'both' with cafe_lrt).")
}

if (length(value_cols) > 0) {
  row_totals <- rowSums(significant_change_df[, value_cols, drop = FALSE])
  significant_change_df <- significant_change_df[!is.na(row_totals) & row_totals > 0, , drop = FALSE]
}

write_tsv_base(significant_change_df, file.path(outdir, paste0("orthogroup_table_significant_", direction, ".tsv")))


# Output orthogroups with significant increase/decrease in target branch
if (!target_branch_id %in% colnames(significant_change_df) && go_method == "cafe_lrt") {
  # No CAFE-significant target events means an empty legacy GO candidate set;
  # it must not prevent the independent native family model comparison.
  significant_change_df[[target_branch_id]] <- rep(0, nrow(significant_change_df))
}
if (!target_branch_id %in% colnames(significant_change_df)) {
  stop("Target branch column was removed after filtering NA-only columns: ", target_branch_id)
}
target_significant_df <- significant_change_df[significant_change_df[[target_branch_id]] == 1, c("FamilyID"), drop = FALSE]
target_significant_df <- merge(target_significant_df, orthogroup_df, by.x = "FamilyID", by.y = "Orthogroup", all.x = TRUE, sort = FALSE)
if (!go_ref_sp %in% colnames(target_significant_df)) {
  target_significant_df[[go_ref_sp]] <- rep(NA_character_, nrow(target_significant_df))
}
target_significant_df <- target_significant_df[, c("FamilyID", go_ref_sp), drop = FALSE]
target_significant_df <- explode_rows(target_significant_df, go_ref_sp, sep = ", ")
colnames(target_significant_df)[colnames(target_significant_df) == go_ref_sp] <- "gene_id"

orthogroup_significant_df <- merge(target_significant_df, ref_annotation_df, by = "gene_id", all.x = TRUE, sort = FALSE)
for (col in c("sprot_best", "sprot_recname")) {
  if (!col %in% colnames(orthogroup_significant_df)) {
    orthogroup_significant_df[[col]] <- rep(NA_character_, nrow(orthogroup_significant_df))
  }
}
orthogroup_significant_df <- orthogroup_significant_df[, c("FamilyID", "gene_id", "sprot_best", "sprot_recname"), drop = FALSE]
orthogroup_significant_df <- summarise_orthogroup_targets(orthogroup_significant_df)
write_tsv_base(orthogroup_significant_df, file.path(outdir, paste0("orthogroup_significant_", direction, "_", target_branch, ".tsv")))


# Link orthogroup with GO terms
if (!go_ref_sp %in% colnames(orthogroup_df)) {
  stop("Reference species column was not found in orthogroup table: ", go_ref_sp)
}
orthogroup_map_df <- orthogroup_df[, c("Orthogroup", go_ref_sp), drop = FALSE]
orthogroup_map_df <- explode_rows(orthogroup_map_df, go_ref_sp, sep = ", ")
colnames(orthogroup_map_df)[colnames(orthogroup_map_df) == "Orthogroup"] <- "FamilyID"
colnames(orthogroup_map_df)[colnames(orthogroup_map_df) == go_ref_sp] <- "gene_id"

for (col in c("gene_id", "sprot_recname", "go_ids", "go_aspects", "go_terms")) {
  if (!col %in% colnames(ref_annotation_df)) {
    ref_annotation_df[[col]] <- rep(NA_character_, nrow(ref_annotation_df))
  }
}
ref_annotation_df <- ref_annotation_df[, c("gene_id", "sprot_recname", "go_ids", "go_aspects", "go_terms"), drop = FALSE]
ref_annotation_df <- explode_rows_multi(ref_annotation_df, c("go_ids", "go_aspects", "go_terms"), sep = "; ")
ref_annotation_df <- ref_annotation_df[ref_annotation_df$go_aspects %in% go_category, , drop = FALSE]

orthogroup_go_df <- merge(orthogroup_map_df, ref_annotation_df, by = "gene_id", all.x = TRUE, sort = FALSE)
orthogroup_go_df <- orthogroup_go_df[!is.na(orthogroup_go_df$go_ids), c("FamilyID", "go_ids", "go_aspects", "go_terms"), drop = FALSE]
orthogroup_go_df <- unique(orthogroup_go_df)
if (nrow(orthogroup_go_df) > 0) {
  orthogroup_go_df <- orthogroup_go_df[order(orthogroup_go_df$FamilyID, orthogroup_go_df$go_ids), , drop = FALSE]
}


# GO enrichment analysis
if (go_method == "cafe_lrt") {
  if (!grepl("_change\\.tab$", args[1])) stop("cafe_lrt requires a native *_change.tab input path.")
  cafe_prefix <- sub("_change\\.tab$", "", args[1])
  native_inputs <- paste0(cafe_prefix, c("_count.tab", "_asr.tre"))
  if (!all(file.exists(native_inputs))) stop("cafe_lrt requires matching native *_count.tab and *_asr.tre files.")
  metadata <- unique(orthogroup_go_df[, c("go_ids", "go_aspects", "go_terms"), drop = FALSE])
  if (anyNA(metadata) || anyDuplicated(metadata$go_ids)) stop("GO IDs require consistent, nonmissing aspect and term metadata.")
  if (anyNA(change_df$FamilyID) || anyDuplicated(change_df$FamilyID)) stop("FamilyID must be unique and nonmissing.")
  family_ids <- intersect(change_df$FamilyID, orthogroup_go_df$FamilyID)
  if (!length(family_ids)) stop("No GO-annotated families available for CAFE model comparison.")
  family_file <- file.path(outdir, "specificity_family_ids.tsv")
  write_tsv_base(data.frame(FamilyID = family_ids), family_file)
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg) != 1) stop("Cannot locate native CAFE adapter.")
  native_script <- file.path(dirname(normalizePath(sub("^--file=", "", script_arg))), "cafe_branch_specificity.py")
  native_dir <- file.path(outdir, "native_cafe")
  native_args <- c(native_script, "--counts", native_inputs[1], "--asr-tree", native_inputs[2],
    "--changes", args[1], "--family-ids", family_file, "--target-branch", target_branch_id,
    "--output-dir", native_dir, "--bootstrap-replicates", bootstrap_replicates,
    "--fit-restarts", fit_restarts, "--max-iterations", max_iterations, "--cores", native_cores)
  if (nzchar(dated_tree)) native_args <- c(native_args, "--tree", dated_tree)
  error_file <- paste0(cafe_prefix, "_error_model.txt")
  if (file.exists(error_file)) native_args <- c(native_args, "--error-model", error_file)
  python <- Sys.which("python")
  if (!nzchar(python)) stop("Python is required for the native CAFE adapter.")
  status <- system2(python, shQuote(native_args))
  if (status != 0) stop("Native CAFE comparison failed; inspect native_cafe/metadata.json and preserved logs.")
  native <- read_tsv_base(file.path(native_dir, "family_lrt.tsv"), na = "NA")
  if (!setequal(native$FamilyID, family_ids)) stop("Native CAFE returned a different family universe.")
  families <- select_cafe_families(native, family_alpha, direction)
  write_tsv_base(families, file.path(outdir, "family_specificity.tsv"))
  directions <- if (direction == "both") c("increase", "decrease") else direction
  go_parts <- lapply(directions, function(d) {
    # Preserve the legacy tested GO IDs for each direction, before LRT selection.
    signed_ids <- change_df$FamilyID[if (d == "increase") change_df[[target_branch_id]] > 0 else change_df[[target_branch_id]] < 0]
    candidate_ids <- intersect(target_significant_df$FamilyID, signed_ids)
    candidate_go <- unique(orthogroup_go_df[orthogroup_go_df$FamilyID %in% candidate_ids,
      c("go_ids", "go_aspects", "go_terms"), drop = FALSE])
    selected <- select_cafe_families(native, family_alpha, d)
    out <- summarise_specific_go(selected, orthogroup_go_df, candidate_go, p_value_threshold)
    out$all$direction <- rep(d, nrow(out$all))
    out$significant$direction <- rep(d, nrow(out$significant))
    out
  })
  go_out <- list(all = do.call(rbind, lapply(go_parts, `[[`, "all")),
                 significant = do.call(rbind, lapply(go_parts, `[[`, "significant")))
  write_tsv_base(data.frame(
    method = go_method, family_alpha = family_alpha, family_adjustment = "BH_before_direction_selection",
    n_tested_families = nrow(families), n_selected_families = sum(families$selected),
    n_go_tests = nrow(go_out$all), go_scope = "legacy_target_observed_go_per_direction",
    interpretation = "native_CAFE_turnover_LRT_with_reconstructed_change_direction",
    target_branch = target_branch_id, direction = direction,
    go_category = paste(go_category, collapse = ","), go_adjustment = "BH_per_direction",
    bootstrap_replicates = as.integer(bootstrap_replicates),
    minimum_bootstrap_p = 1 / (as.integer(bootstrap_replicates) + 1)
  ), file.path(outdir, "specificity_metadata.tsv"))
  write_tsv_base(go_out$all, file.path(outdir, paste0("enrichment_significant_", direction, "_", target_branch, "_all_go.tsv")))
  write_tsv_base(go_out$significant, file.path(outdir, paste0("enrichment_significant_", direction, "_", target_branch, "_significant_go.tsv")))
  cat("cafe_lrt: selected", sum(families$selected), "families; see native_cafe/metadata.json and docs/go-enrichment.md.\n")
  quit(save = "no", status = 0)
}
all_families_with_go <- intersect(significant_change_df$FamilyID, orthogroup_go_df$FamilyID)
target_families_with_go <- intersect(significant_change_df$FamilyID[significant_change_df[[target_branch_id]] == 1], all_families_with_go)

if (length(target_families_with_go) == 0) {
  stop("No significant families with GO annotation in the target branch: ", target_branch_id)
}

event_df <- significant_change_df[significant_change_df$FamilyID %in% all_families_with_go, , drop = FALSE]
event_df <- pivot_longer_base(event_df, "FamilyID")
event_df <- event_df[event_df$is_change == 1, , drop = FALSE]
event_df$is_target <- (event_df$branch == target_branch_id)

target_events <- event_df[event_df$is_target, , drop = FALSE]
other_events <- event_df[!event_df$is_target, , drop = FALSE]
total_target_events <- nrow(target_events)
total_other_events <- nrow(other_events)

if (total_target_events == 0) {
  stop("No change events in the target branch after filtering.")
}
if (total_other_events == 0) {
  stop("No change events in other branches after filtering.")
}

event_go_df <- merge(event_df, orthogroup_go_df, by = "FamilyID", all.x = TRUE, sort = FALSE)
event_go_df <- event_go_df[!is.na(event_go_df$go_ids), , drop = FALSE]

go_out <- summarise_go_enrichment(
  event_go_df = event_go_df,
  total_target_events = total_target_events,
  total_other_events = total_other_events,
  p_value_threshold = p_value_threshold
)

write_tsv_base(go_out$all, file.path(outdir, paste0("enrichment_significant_", direction, "_", target_branch, "_all_go.tsv")))
write_tsv_base(go_out$significant, file.path(outdir, paste0("enrichment_significant_", direction, "_", target_branch, "_significant_go.tsv")))

cat(as.character(Sys.time()), "cafe_go_enrichment.r completed successfully. Exiting\n")
