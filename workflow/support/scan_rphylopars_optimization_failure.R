#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[[1]] else file.path("workspace", "output", "rphylopars_failure_scan")
n_random <- if (length(args) >= 2) suppressWarnings(as.integer(args[[2]])) else 1800L
n_edge <- if (length(args) >= 3) suppressWarnings(as.integer(args[[3]])) else 700L
seed <- if (length(args) >= 4) suppressWarnings(as.integer(args[[4]])) else 20260220L
ncpu <- if (length(args) >= 5) suppressWarnings(as.integer(args[[5]])) else NA_integer_
if (is.na(n_random) || n_random < 0) n_random <- 1800L
if (is.na(n_edge) || n_edge < 0) n_edge <- 700L
if (is.na(seed)) seed <- 20260220L
if (is.na(ncpu) || ncpu < 1) {
  ncpu <- suppressWarnings(as.integer(Sys.getenv("NSLOTS", unset = "0")))
}
if (is.na(ncpu) || ncpu < 1) {
  ncpu <- max(1L, min(8L, suppressWarnings(parallel::detectCores(logical = TRUE))))
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ape)
  library(Rphylopars)
  library(parallel)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_path <- sub("^--file=", "", script_arg[1])
  script_path <- gsub("~\\+~", " ", script_path)
  script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
} else {
  script_dir <- getwd()
}
source(file.path(script_dir, "pgls_common.R"))

normalize_error_type <- function(msg) {
  if (is.null(msg) || length(msg) == 0 || is.na(msg) || !nzchar(msg)) return("none")
  m <- tolower(as.character(msg))
  if (grepl("non[- ]positive definite|not positive definite|leading minor|lapack|chol|eigen", m)) return("matrix_non_pd")
  if (grepl("na/nan/inf|nan|infinite|non-finite|missing value", m)) return("non_finite")
  if (grepl("singular|computationally singular|system is exactly singular", m)) return("singular_matrix")
  if (grepl("subscript out of bounds|replacement has length zero|argument is of length zero", m)) return("index_or_length")
  if (grepl("cannot allocate|memory|vector memory exhausted", m)) return("memory")
  if (grepl("did not converge|maxit|maximum number of", m)) return("non_convergence")
  return("other")
}

safe_as_num <- function(x) {
  out <- suppressWarnings(as.numeric(x))
  if (length(out) == 0) return(NA_real_)
  out[1]
}

log_uniform <- function(min_v, max_v) {
  exp(runif(1, log(min_v), log(max_v)))
}

build_one_config_random <- function() {
  n_species <- sample(c(4L, 5L, 6L, 8L, 10L, 12L, 16L, 20L), 1)
  rep_min <- sample(1:3, 1)
  rep_max <- sample(rep_min:7, 1)
  list(
    scenario = "random",
    n_species = n_species,
    rep_min = rep_min,
    rep_max = rep_max,
    miss_prob = runif(1, 0, 0.92),
    drop_species_prob = runif(1, 0, 0.45),
    sigma_phylo_x = log_uniform(1e-7, 3),
    sigma_phylo_y = log_uniform(1e-7, 3),
    sigma_within_x = log_uniform(1e-8, 30),
    beta = runif(1, -2.5, 2.5),
    branch_scale = log_uniform(1e-8, 2),
    near_zero_edge_frac = runif(1, 0, 0.85),
    outlier_prob = runif(1, 0, 0.2),
    outlier_scale = log_uniform(1, 1e6)
  )
}

build_edge_config_pool <- function() {
  edge_df <- expand.grid(
    n_species = c(4L, 6L, 12L, 20L),
    rep_mode = c("tight", "wide"),
    miss_prob = c(0, 0.5, 0.85),
    drop_species_prob = c(0, 0.2, 0.4),
    sigma_within_x = c(0, 1e-8, 1e-4, 1, 20),
    sigma_phylo_x = c(1e-7, 0.05, 1),
    sigma_phylo_y = c(1e-7, 0.05, 1),
    beta = c(0, 1.5),
    branch_scale = c(1, 1e-4, 1e-8),
    near_zero_edge_frac = c(0, 0.8),
    outlier_prob = c(0, 0.15),
    outlier_scale = c(1, 1e4),
    stringsAsFactors = FALSE
  )
  configs <- lapply(seq_len(nrow(edge_df)), function(i) {
    row <- edge_df[i, , drop = FALSE]
    if (row$rep_mode == "tight") {
      rep_min <- 1L
      rep_max <- 2L
    } else {
      rep_min <- 2L
      rep_max <- 6L
    }
    list(
      scenario = "edge",
      n_species = as.integer(row$n_species),
      rep_min = rep_min,
      rep_max = rep_max,
      miss_prob = safe_as_num(row$miss_prob),
      drop_species_prob = safe_as_num(row$drop_species_prob),
      sigma_phylo_x = safe_as_num(row$sigma_phylo_x),
      sigma_phylo_y = safe_as_num(row$sigma_phylo_y),
      sigma_within_x = safe_as_num(row$sigma_within_x),
      beta = safe_as_num(row$beta),
      branch_scale = safe_as_num(row$branch_scale),
      near_zero_edge_frac = safe_as_num(row$near_zero_edge_frac),
      outlier_prob = safe_as_num(row$outlier_prob),
      outlier_scale = safe_as_num(row$outlier_scale)
    )
  })
  configs
}

simulate_dataset <- function(cfg) {
  species <- sprintf("sp%03d", seq_len(cfg$n_species))
  tree <- ape::rtree(cfg$n_species)
  tree$tip.label <- species
  tree$edge.length <- suppressWarnings(as.numeric(tree$edge.length))
  tree$edge.length[!is.finite(tree$edge.length)] <- 1
  tree$edge.length <- tree$edge.length * cfg$branch_scale
  tree$edge.length[tree$edge.length < 0] <- 0
  if (cfg$near_zero_edge_frac > 0) {
    n_e <- length(tree$edge.length)
    n_z <- max(1L, floor(n_e * cfg$near_zero_edge_frac))
    idx <- sample.int(n_e, n_z, replace = FALSE)
    tree$edge.length[idx] <- tree$edge.length[idx] * 1e-12
  }
  tree$edge.length[!is.finite(tree$edge.length)] <- 0

  x_sp <- try(as.numeric(ape::rTraitCont(tree, model = "BM", sigma = cfg$sigma_phylo_x)), silent = TRUE)
  if (inherits(x_sp, "try-error") || length(x_sp) != cfg$n_species) x_sp <- stats::rnorm(cfg$n_species, 0, 1e-6)
  y_noise <- try(as.numeric(ape::rTraitCont(tree, model = "BM", sigma = cfg$sigma_phylo_y)), silent = TRUE)
  if (inherits(y_noise, "try-error") || length(y_noise) != cfg$n_species) y_noise <- stats::rnorm(cfg$n_species, 0, 1e-6)
  y_sp <- cfg$beta * x_sp + y_noise

  max_rep <- cfg$rep_max
  wide <- data.frame(species = species, y = as.numeric(y_sp), stringsAsFactors = FALSE)
  for (r in seq_len(max_rep)) {
    wide[[sprintf("x_%d", r)]] <- NA_real_
  }
  rep_counts <- integer(cfg$n_species)

  for (i in seq_len(cfg$n_species)) {
    sp <- species[i]
    k <- sample(cfg$rep_min:cfg$rep_max, 1)
    rep_counts[i] <- k
    x_rep <- x_sp[i] + stats::rnorm(k, mean = 0, sd = cfg$sigma_within_x)
    if (cfg$outlier_prob > 0 && stats::runif(1) < cfg$outlier_prob) {
      outlier <- stats::rnorm(1, mean = 0, sd = cfg$outlier_scale)
      x_rep[sample.int(k, 1)] <- x_rep[sample.int(k, 1)] + outlier
    }
    keep <- stats::runif(k) > cfg$miss_prob
    if (stats::runif(1) < cfg$drop_species_prob) {
      keep[] <- FALSE
    }
    x_rep[!keep] <- NA_real_
    wide[i, sprintf("x_%d", seq_len(k))] <- x_rep
  }

  x_cols <- grep("^x_[0-9]+$", colnames(wide), value = TRUE)
  wide$mean_x <- apply(wide[, x_cols, drop = FALSE], 1, function(v) {
    vv <- suppressWarnings(as.numeric(v))
    if (all(is.na(vv))) return(NA_real_)
    mean(vv, na.rm = TRUE)
  })

  list(tree = tree, wide = wide, x_cols = x_cols, rep_counts = rep_counts)
}

drop_tree_to_species <- function(tree, species_vec) {
  keep <- unique(as.character(species_vec))
  t_use <- tree
  drop_tips <- t_use$tip.label[!(t_use$tip.label %in% keep)]
  if (length(drop_tips) > 0) t_use <- ape::drop.tip(t_use, drop_tips)
  t_use
}

run_raw_phenocov <- function(td, t_use, pcl) {
  fit <- suppressWarnings(try(do.call(Rphylopars::phylopars.lm, list(
    formula = y ~ x,
    trait_data = td,
    tree = t_use,
    model = "BM",
    pheno_error = FALSE,
    phylo_correlated = TRUE,
    pheno_correlated = FALSE,
    phenocov_list = pcl
  )), silent = TRUE))
  if (!inherits(fit, "try-error")) {
    return(list(status = "ok", error_message = "", pval = safe_as_num(fit$pval)))
  }
  list(status = "error", error_message = extract_try_error_message(fit), pval = NA_real_)
}

run_mean_fallback <- function(wide, tree) {
  td <- wide[, c("species", "y", "mean_x"), drop = FALSE]
  colnames(td)[3] <- "x"
  td$y <- suppressWarnings(as.numeric(td$y))
  td$x <- suppressWarnings(as.numeric(td$x))
  td <- td[is.finite(td$y) & is.finite(td$x) & !is.na(td$species) & nzchar(td$species), , drop = FALSE]
  if (nrow(td) < 3 || length(unique(td$species)) < 3) {
    return(list(status = "no_data", fit_mode = "mean_fallback_no_data", error_message = "", pval = NA_real_))
  }
  t_use <- drop_tree_to_species(tree, td$species)
  td <- td[td$species %in% t_use$tip.label, , drop = FALSE]
  fit <- suppressWarnings(phylopars.lm.safe(
    formula_obj = y ~ x,
    trait_data = td,
    tree = t_use,
    model = "BM",
    pheno_error = TRUE,
    phylo_correlated = TRUE,
    pheno_correlated = TRUE,
    phenocov_list = list(),
    trait_col = "y",
    expression_col = "x",
    fit_mode_label = "mean_fallback",
    verbose = FALSE
  ))
  if (is.null(fit$fit)) {
    return(list(status = "error", fit_mode = fit$fit_mode, error_message = fit$error_message, pval = NA_real_))
  }
  list(status = "ok", fit_mode = fit$fit_mode, error_message = "", pval = safe_as_num(fit$fit$pval))
}

evaluate_case <- function(cfg, case_id) {
  dat <- simulate_dataset(cfg)
  tree <- dat$tree
  wide <- dat$wide
  x_cols <- dat$x_cols
  rep_counts <- dat$rep_counts

  ph <- build_phenocov_input(wide, trait_col = "y", expression_cols = x_cols, expression_base = "x")
  td <- ph[["trait_data"]]
  pcl <- ph[["phenocov_list"]]
  has_var <- isTRUE(ph[["has_within_species_variation"]])

  se2 <- rep(NA_real_, length(pcl))
  if (length(pcl) > 0) {
    se2 <- sapply(pcl, function(m) {
      if (is.matrix(m) && all(dim(m) >= c(2, 2))) safe_as_num(m[2, 2]) else NA_real_
    })
  }

  if (nrow(td) > 0) {
    t_use <- drop_tree_to_species(tree, td$species)
    td <- td[td$species %in% t_use$tip.label, , drop = FALSE]
    pcl <- pcl[td$species]
    se2 <- se2[td$species]
  } else {
    t_use <- tree
  }

  analyzable <- (nrow(td) >= 3 && length(unique(td$species)) >= 3)
  vcv_kappa <- NA_real_
  if (length(t_use$tip.label) >= 3) {
    vcv_kappa <- suppressWarnings(tryCatch(kappa(ape::vcv.phylo(t_use), exact = FALSE), error = function(e) NA_real_))
  }

  raw_status <- "no_data"
  raw_error <- ""
  raw_pval <- NA_real_
  safe_status <- "no_data"
  safe_error <- ""
  safe_fit_mode <- "no_data"
  safe_pval <- NA_real_
  pipeline_status <- "no_data"
  pipeline_fit_mode <- "no_data"
  pipeline_error <- ""
  pipeline_pval <- NA_real_

  if (analyzable) {
    raw <- run_raw_phenocov(td, t_use, pcl)
    raw_status <- raw$status
    raw_error <- raw$error_message
    raw_pval <- raw$pval

    safe <- suppressWarnings(phylopars.lm.safe(
      formula_obj = y ~ x,
      trait_data = td,
      tree = t_use,
      model = "BM",
      pheno_error = FALSE,
      phylo_correlated = TRUE,
      pheno_correlated = FALSE,
      phenocov_list = pcl,
      trait_col = "y",
      expression_col = "x",
      fit_mode_label = "phenocov",
      verbose = FALSE
    ))
    if (is.null(safe$fit)) {
      safe_status <- "error"
      safe_error <- safe$error_message
      safe_fit_mode <- safe$fit_mode
    } else {
      safe_status <- "ok"
      safe_fit_mode <- safe$fit_mode
      safe_pval <- safe_as_num(safe$fit$pval)
    }

    if (!has_var) {
      mean_out <- run_mean_fallback(wide, tree)
      pipeline_status <- mean_out$status
      pipeline_fit_mode <- "phenocov_no_variation"
      if (mean_out$status == "ok") {
        pipeline_fit_mode <- mean_out$fit_mode
      }
      pipeline_error <- mean_out$error_message
      pipeline_pval <- mean_out$pval
    } else if (safe_status == "ok") {
      pipeline_status <- "ok"
      pipeline_fit_mode <- safe_fit_mode
      pipeline_pval <- safe_pval
    } else {
      mean_out <- run_mean_fallback(wide, tree)
      pipeline_status <- mean_out$status
      pipeline_fit_mode <- mean_out$fit_mode
      pipeline_error <- mean_out$error_message
      pipeline_pval <- mean_out$pval
    }
  }

  data.frame(
    case_id = case_id,
    scenario = cfg$scenario,
    n_species_input = cfg$n_species,
    rep_min = cfg$rep_min,
    rep_max = cfg$rep_max,
    rep_mean_observed = mean(rep_counts),
    missing_rate_input = cfg$miss_prob,
    drop_species_prob_input = cfg$drop_species_prob,
    missing_rate_actual = mean(is.na(as.matrix(wide[, x_cols, drop = FALSE]))),
    sigma_phylo_x = cfg$sigma_phylo_x,
    sigma_phylo_y = cfg$sigma_phylo_y,
    sigma_within_x = cfg$sigma_within_x,
    beta = cfg$beta,
    branch_scale = cfg$branch_scale,
    near_zero_edge_frac = cfg$near_zero_edge_frac,
    outlier_prob = cfg$outlier_prob,
    outlier_scale = cfg$outlier_scale,
    n_species_usable = nrow(td),
    n_species_unique_usable = length(unique(td$species)),
    has_within_species_variation = has_var,
    se2_min = suppressWarnings(min(se2, na.rm = TRUE)),
    se2_median = suppressWarnings(stats::median(se2, na.rm = TRUE)),
    se2_max = suppressWarnings(max(se2, na.rm = TRUE)),
    min_edge = suppressWarnings(min(tree$edge.length, na.rm = TRUE)),
    max_edge = suppressWarnings(max(tree$edge.length, na.rm = TRUE)),
    tree_vcv_kappa = vcv_kappa,
    analyzable = analyzable,
    raw_status = raw_status,
    raw_fit_error_type = normalize_error_type(raw_error),
    raw_fit_error_message = substr(as.character(raw_error), 1, 220),
    raw_pval = raw_pval,
    safe_status = safe_status,
    safe_fit_mode = safe_fit_mode,
    safe_fit_error_type = normalize_error_type(safe_error),
    safe_fit_error_message = substr(as.character(safe_error), 1, 220),
    safe_pval = safe_pval,
    pipeline_status = pipeline_status,
    pipeline_fit_mode = pipeline_fit_mode,
    pipeline_error_type = normalize_error_type(pipeline_error),
    pipeline_error_message = substr(as.character(pipeline_error), 1, 220),
    pipeline_pval = pipeline_pval,
    stringsAsFactors = FALSE
  )
}

set.seed(seed)
edge_pool <- build_edge_config_pool()
if (length(edge_pool) > 0 && n_edge > 0) {
  edge_idx <- sample.int(length(edge_pool), size = min(n_edge, length(edge_pool)), replace = FALSE)
  edge_cfg <- edge_pool[edge_idx]
} else {
  edge_cfg <- list()
}
random_cfg <- if (n_random > 0) replicate(n_random, build_one_config_random(), simplify = FALSE) else list()
all_cfg <- c(edge_cfg, random_cfg)

cat("scan_rphylopars_optimization_failure.R\n")
cat("seed:", seed, "\n")
cat("ncpu:", ncpu, "\n")
cat("cases_total:", length(all_cfg), "\n")
cat("cases_edge:", length(edge_cfg), "\n")
cat("cases_random:", length(random_cfg), "\n")

if (length(all_cfg) == 0) {
  stop("No configurations were generated.")
}

runner <- function(i) evaluate_case(all_cfg[[i]], i)
if (ncpu > 1) {
  rows <- parallel::mclapply(seq_along(all_cfg), runner, mc.cores = ncpu, mc.preschedule = TRUE)
} else {
  rows <- lapply(seq_along(all_cfg), runner)
}
df <- do.call(rbind, rows)

num_cols <- c("se2_min", "se2_median", "se2_max", "tree_vcv_kappa")
for (cn in num_cols) {
  bad <- !is.finite(df[[cn]])
  if (any(bad)) df[[cn]][bad] <- NA_real_
}

results_tsv <- file.path(out_dir, "rphylopars_failure_scan.results.tsv.gz")
write.table(df, gzfile(results_tsv), sep = "\t", row.names = FALSE, quote = FALSE)

df_a <- df[df$analyzable, , drop = FALSE]

overall <- data.frame(
  metric = c(
    "cases_total",
    "cases_analyzable",
    "raw_fail_rate",
    "safe_fail_rate",
    "pipeline_fail_rate",
    "pipeline_mean_fallback_rate",
    "pipeline_no_variation_rate"
  ),
  value = c(
    nrow(df),
    nrow(df_a),
    if (nrow(df_a) > 0) mean(df_a$raw_status == "error") else NA_real_,
    if (nrow(df_a) > 0) mean(df_a$safe_status == "error") else NA_real_,
    if (nrow(df_a) > 0) mean(df_a$pipeline_status == "error") else NA_real_,
    if (nrow(df_a) > 0) mean(grepl("^mean_fallback", df_a$pipeline_fit_mode)) else NA_real_,
    if (nrow(df_a) > 0) mean(df_a$pipeline_fit_mode == "phenocov_no_variation") else NA_real_
  ),
  stringsAsFactors = FALSE
)
write.table(overall, file.path(out_dir, "rphylopars_failure_scan.summary_overall.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

tab_raw_err <- if (nrow(df_a) > 0) as.data.frame(table(df_a$raw_fit_error_type), stringsAsFactors = FALSE) else data.frame()
if (nrow(tab_raw_err) > 0) {
  colnames(tab_raw_err) <- c("error_type", "count")
  tab_raw_err <- tab_raw_err[order(tab_raw_err$count, decreasing = TRUE), , drop = FALSE]
}
write.table(tab_raw_err, file.path(out_dir, "rphylopars_failure_scan.raw_error_type.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

tab_safe_err <- if (nrow(df_a) > 0) as.data.frame(table(df_a$safe_fit_error_type), stringsAsFactors = FALSE) else data.frame()
if (nrow(tab_safe_err) > 0) {
  colnames(tab_safe_err) <- c("error_type", "count")
  tab_safe_err <- tab_safe_err[order(tab_safe_err$count, decreasing = TRUE), , drop = FALSE]
}
write.table(tab_safe_err, file.path(out_dir, "rphylopars_failure_scan.safe_error_type.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

df_a$miss_bin <- cut(df_a$missing_rate_actual, breaks = c(-Inf, 0.2, 0.5, 0.8, Inf), labels = c("<=0.2", "0.2-0.5", "0.5-0.8", ">0.8"))
df_a$species_bin <- cut(df_a$n_species_usable, breaks = c(-Inf, 5, 10, Inf), labels = c("<=5", "6-10", ">10"))
df_a$within_bin <- cut(df_a$sigma_within_x, breaks = c(-Inf, 1e-4, 1, Inf), labels = c("<1e-4", "1e-4-1", ">1"))
df_a$branch_bin <- cut(df_a$branch_scale, breaks = c(-Inf, 1e-6, 1e-3, Inf), labels = c("<1e-6", "1e-6-1e-3", ">1e-3"))

summarize_fail_by_bin <- function(df_in, bin_col) {
  if (nrow(df_in) == 0) return(data.frame())
  bins <- unique(as.character(df_in[[bin_col]]))
  bins <- bins[!is.na(bins)]
  do.call(rbind, lapply(bins, function(b) {
    dd <- df_in[as.character(df_in[[bin_col]]) == b, , drop = FALSE]
    data.frame(
      variable = bin_col,
      bin = b,
      n = nrow(dd),
      raw_fail_rate = mean(dd$raw_status == "error"),
      safe_fail_rate = mean(dd$safe_status == "error"),
      pipeline_fail_rate = mean(dd$pipeline_status == "error"),
      mean_fallback_rate = mean(grepl("^mean_fallback", dd$pipeline_fit_mode)),
      stringsAsFactors = FALSE
    )
  }))
}

bin_summaries <- rbind(
  summarize_fail_by_bin(df_a, "species_bin"),
  summarize_fail_by_bin(df_a, "miss_bin"),
  summarize_fail_by_bin(df_a, "within_bin"),
  summarize_fail_by_bin(df_a, "branch_bin")
)
write.table(bin_summaries, file.path(out_dir, "rphylopars_failure_scan.fail_by_bin.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

hot_raw_fail <- df_a[df_a$raw_status == "error", c(
  "case_id", "scenario", "n_species_input", "rep_min", "rep_max", "missing_rate_actual",
  "sigma_phylo_x", "sigma_phylo_y", "sigma_within_x", "branch_scale", "near_zero_edge_frac",
  "outlier_prob", "outlier_scale", "n_species_usable", "has_within_species_variation",
  "se2_min", "se2_median", "se2_max", "tree_vcv_kappa", "raw_fit_error_type", "raw_fit_error_message",
  "safe_status", "safe_fit_mode", "pipeline_status", "pipeline_fit_mode"
), drop = FALSE]
if (nrow(hot_raw_fail) > 0) {
  ord <- order(hot_raw_fail$tree_vcv_kappa, decreasing = TRUE, na.last = TRUE)
  hot_raw_fail <- hot_raw_fail[ord, , drop = FALSE]
}
write.table(hot_raw_fail, gzfile(file.path(out_dir, "rphylopars_failure_scan.raw_fail_cases.tsv.gz")), sep = "\t", row.names = FALSE, quote = FALSE)

report_file <- file.path(out_dir, "rphylopars_failure_scan.report.md")
con <- file(report_file, open = "wt")
writeLines(c(
  "# Rphylopars Optimization Failure Scan",
  "",
  paste0("- Seed: ", seed),
  paste0("- CPU cores: ", ncpu),
  paste0("- Cases total: ", nrow(df)),
  paste0("- Cases analyzable: ", nrow(df_a)),
  "",
  "## Overall",
  paste0("- Raw phenocov failure rate: ", sprintf("%.4f", overall$value[overall$metric == "raw_fail_rate"])),
  paste0("- Safe wrapper failure rate: ", sprintf("%.4f", overall$value[overall$metric == "safe_fail_rate"])),
  paste0("- Pipeline final failure rate: ", sprintf("%.4f", overall$value[overall$metric == "pipeline_fail_rate"])),
  paste0("- Mean fallback rate: ", sprintf("%.4f", overall$value[overall$metric == "pipeline_mean_fallback_rate"])),
  paste0("- No-variation fallback rate: ", sprintf("%.4f", overall$value[overall$metric == "pipeline_no_variation_rate"])),
  "",
  "## Files",
  paste0("- Results: ", basename(results_tsv)),
  "- Overall summary: rphylopars_failure_scan.summary_overall.tsv",
  "- Raw error types: rphylopars_failure_scan.raw_error_type.tsv",
  "- Safe error types: rphylopars_failure_scan.safe_error_type.tsv",
  "- Fail rates by bin: rphylopars_failure_scan.fail_by_bin.tsv",
  "- Raw fail cases: rphylopars_failure_scan.raw_fail_cases.tsv.gz"
), con = con)
close(con)

cat("Wrote:", results_tsv, "\n")
cat("Wrote:", report_file, "\n")
