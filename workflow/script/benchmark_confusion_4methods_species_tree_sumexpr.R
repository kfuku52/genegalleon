#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[[1]] else file.path("workspace", "output", "benchmark_epgls")
n_iter <- if (length(args) >= 2) suppressWarnings(as.integer(args[[2]])) else 80L
rrpp_iter <- if (length(args) >= 3) suppressWarnings(as.integer(args[[3]])) else 49L
alpha <- if (length(args) >= 4) suppressWarnings(as.numeric(args[[4]])) else 0.05
rphylopars_mode <- if (length(args) >= 5) tolower(args[[5]]) else "mean"
if (is.na(n_iter) || n_iter < 10) n_iter <- 80L
if (is.na(rrpp_iter) || rrpp_iter < 9) rrpp_iter <- 49L
if (is.na(alpha) || alpha <= 0 || alpha >= 1) alpha <- 0.05
if (!(rphylopars_mode %in% c("mean", "phenocov"))) rphylopars_mode <- "mean"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ape)
  library(Rphylopars)
  library(ggplot2)
  library(nlme)
})
.libPaths(c("/tmp/epgls_lib", .libPaths()))
suppressPackageStartupMessages(library(RRPP))

patch_rrpp_detect_cores <- function() {
  dc <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (is.finite(dc) && !is.na(dc) && dc >= 1) return(invisible(NULL))

  rrpp_imports <- parent.env(asNamespace("RRPP"))
  if (!exists("detectCores", envir = rrpp_imports, inherits = FALSE)) {
    return(invisible(NULL))
  }

  orig_detect <- get("detectCores", envir = rrpp_imports, inherits = FALSE)
  if (bindingIsLocked("detectCores", rrpp_imports)) {
    unlockBinding("detectCores", rrpp_imports)
  }
  assign("detectCores", function(all.tests = FALSE, logical = TRUE) {
    out <- suppressWarnings(orig_detect(all.tests = all.tests, logical = logical))
    if (is.na(out) || !is.finite(out) || out < 1) return(1L)
    as.integer(out)
  }, envir = rrpp_imports)
  lockBinding("detectCores", rrpp_imports)
  invisible(NULL)
}
patch_rrpp_detect_cores()

build_obs_cov_from_species_tree <- function(dat, species_tree) {
  if (nrow(dat) == 0) return(list(data = dat, Cov = NULL))
  cov_species <- ape::vcv.phylo(species_tree)
  sp_ids <- rownames(cov_species)
  dat <- dat[dat[["species"]] %in% sp_ids, , drop = FALSE]
  if (nrow(dat) == 0) return(list(data = dat, Cov = NULL))
  idx <- match(dat[["species"]], sp_ids)
  cov_obs <- cov_species[idx, idx, drop = FALSE]
  eps <- suppressWarnings(max(diag(cov_obs), na.rm = TRUE)) * 1e-8
  if (!is.finite(eps) || eps <= 0) eps <- 1e-8
  diag(cov_obs) <- diag(cov_obs) + eps
  obs_ids <- paste0("obs_", seq_len(nrow(dat)))
  rownames(cov_obs) <- obs_ids
  colnames(cov_obs) <- obs_ids
  rownames(dat) <- obs_ids
  list(data = dat, Cov = cov_obs)
}

simulate_species_tree_sum_dataset <- function(n_species = 12,
                                              dup_fraction = 0.45,
                                              max_copy = 4,
                                              min_rep = 2,
                                              max_rep = 5,
                                              miss_prob = 0.25,
                                              beta = 0.7,
                                              sigma_species_latent = 1.0,
                                              sigma_species_noise = 1.0,
                                              sigma_copy = 0.6,
                                              sigma_rep = 0.5) {
  species <- sprintf("sp%02d", seq_len(n_species))
  species_tree <- ape::rtree(n_species)
  species_tree$tip.label <- species

  copy_count <- rep(1L, n_species)
  names(copy_count) <- species
  n_dup <- max(1L, floor(n_species * dup_fraction))
  dup_species <- sample(species, size = n_dup, replace = FALSE)
  for (sp in dup_species) {
    copy_count[[sp]] <- sample(2:max_copy, 1)
  }

  z_sp <- as.numeric(ape::rTraitCont(species_tree, model = "BM", sigma = sigma_species_latent))
  names(z_sp) <- species_tree$tip.label

  gene_base <- setNames(vector("list", length(species)), species)
  for (sp in species) {
    k <- copy_count[[sp]]
    gene_base[[sp]] <- z_sp[[sp]] + stats::rnorm(k, mean = 0, sd = sigma_copy)
  }
  x_latent_sum <- sapply(species, function(sp) sum(gene_base[[sp]]))
  x_latent_scaled <- as.numeric(scale(x_latent_sum))
  names(x_latent_scaled) <- species

  y_noise <- as.numeric(ape::rTraitCont(species_tree, model = "BM", sigma = sigma_species_noise))
  names(y_noise) <- species_tree$tip.label
  y_sp <- beta * x_latent_scaled + as.numeric(y_noise[species])
  names(y_sp) <- species

  long <- data.frame(
    species = character(0),
    replicate = integer(0),
    y = numeric(0),
    x_sum = numeric(0),
    n_copy_total = integer(0),
    n_copy_obs = integer(0),
    stringsAsFactors = FALSE
  )

  for (sp in species) {
    k_rep <- sample(min_rep:max_rep, 1)
    k_copy <- length(gene_base[[sp]])
    for (r in seq_len(k_rep)) {
      vals <- gene_base[[sp]] + stats::rnorm(k_copy, mean = 0, sd = sigma_rep)
      keep <- stats::runif(k_copy) > miss_prob
      if (!any(keep)) keep[sample.int(k_copy, 1)] <- TRUE
      x_sum <- sum(vals[keep])
      long <- rbind(long, data.frame(
        species = sp,
        replicate = r,
        y = y_sp[[sp]],
        x_sum = x_sum,
        n_copy_total = k_copy,
        n_copy_obs = sum(keep),
        stringsAsFactors = FALSE
      ))
    }
  }

  wide <- aggregate(
    cbind(y, x_sum) ~ species,
    data = long,
    FUN = function(v) mean(as.numeric(v), na.rm = TRUE)
  )
  colnames(wide)[colnames(wide) == "x_sum"] <- "mean_x_sum"
  wide$n_copy_total <- as.integer(copy_count[wide$species])
  rownames(wide) <- NULL

  list(
    species_tree = species_tree,
    long = long,
    wide_species = wide,
    copy_count = copy_count
  )
}

fit_rphylopars_species_tree <- function(wide_species, species_tree) {
  out <- list(status = "error", pval = NA_real_)
  td <- wide_species[, c("species", "y", "mean_x_sum"), drop = FALSE]
  colnames(td)[3] <- "x"
  td$y <- suppressWarnings(as.numeric(td$y))
  td$x <- suppressWarnings(as.numeric(td$x))
  td <- td[is.finite(td$y) & is.finite(td$x) & nzchar(td$species), , drop = FALSE]
  if (nrow(td) < 4 || length(unique(td$species)) < 4) {
    out$status <- "no_data"
    return(out)
  }
  keep_tips <- intersect(species_tree$tip.label, td$species)
  if (length(keep_tips) < 4) {
    out$status <- "no_data"
    return(out)
  }
  t_use <- species_tree
  drop_tips <- t_use$tip.label[!(t_use$tip.label %in% keep_tips)]
  if (length(drop_tips) > 0) t_use <- ape::drop.tip(t_use, drop_tips)
  td <- td[td$species %in% t_use$tip.label, , drop = FALSE]
  fit <- try(do.call(Rphylopars::phylopars.lm, list(
    formula = y ~ x,
    trait_data = td,
    tree = t_use,
    model = "BM",
    pheno_error = TRUE,
    phylo_correlated = TRUE,
    pheno_correlated = TRUE
  )), silent = TRUE)
  if (inherits(fit, "try-error")) {
    out$status <- "error"
    return(out)
  }
  out$status <- "ok"
  out$pval <- suppressWarnings(as.numeric(fit$pval))
  out
}

fit_rphylopars_species_tree_phenocov <- function(long, species_tree) {
  out <- list(status = "error", pval = NA_real_)
  dat <- long[, c("species", "y", "x_sum"), drop = FALSE]
  colnames(dat)[3] <- "x"
  dat <- dat[is.finite(dat$y) & is.finite(dat$x) & nzchar(dat$species), , drop = FALSE]
  if (nrow(dat) < 8 || length(unique(dat$species)) < 4) {
    out$status <- "no_data"
    return(out)
  }

  split_x <- split(dat$x, dat$species)
  sp <- names(split_x)
  x_mean <- sapply(split_x, function(v) mean(v, na.rm = TRUE))
  x_var <- sapply(split_x, function(v) {
    vv <- v[is.finite(v)]
    if (length(vv) > 1) stats::var(vv) else 0
  })
  x_n <- sapply(split_x, function(v) sum(is.finite(v)))
  x_se2 <- x_var / pmax(1, x_n)
  y_by <- tapply(dat$y, dat$species, function(v) mean(v, na.rm = TRUE))

  td <- data.frame(
    species = sp,
    y = as.numeric(y_by[sp]),
    x = as.numeric(x_mean),
    stringsAsFactors = FALSE
  )
  td <- td[is.finite(td$y) & is.finite(td$x) & nzchar(td$species), , drop = FALSE]
  if (nrow(td) < 4 || length(unique(td$species)) < 4) {
    out$status <- "no_data"
    return(out)
  }

  keep_tips <- intersect(species_tree$tip.label, td$species)
  if (length(keep_tips) < 4) {
    out$status <- "no_data"
    return(out)
  }
  t_use <- species_tree
  drop_tips <- t_use$tip.label[!(t_use$tip.label %in% keep_tips)]
  if (length(drop_tips) > 0) t_use <- ape::drop.tip(t_use, drop_tips)
  td <- td[td$species %in% t_use$tip.label, , drop = FALSE]
  if (nrow(td) < 4 || length(unique(td$species)) < 4) {
    out$status <- "no_data"
    return(out)
  }

  x_se2 <- x_se2[td$species]
  x_se2[!is.finite(x_se2)] <- 0
  x_se2[x_se2 < 0] <- 0

  phenocov_list <- vector("list", length = nrow(td))
  names(phenocov_list) <- td$species
  for (j in seq_len(nrow(td))) {
    m <- matrix(c(0, 0, 0, x_se2[[j]]), nrow = 2, byrow = TRUE)
    rownames(m) <- c("y", "x")
    colnames(m) <- c("y", "x")
    phenocov_list[[j]] <- m
  }

  fit <- try(do.call(Rphylopars::phylopars.lm, list(
    formula = y ~ x,
    trait_data = td,
    tree = t_use,
    model = "BM",
    pheno_error = FALSE,
    phylo_correlated = TRUE,
    pheno_correlated = FALSE,
    phenocov_list = phenocov_list
  )), silent = TRUE)
  if (inherits(fit, "try-error")) {
    out$status <- "error"
    return(out)
  }
  out$status <- "ok"
  out$pval <- suppressWarnings(as.numeric(fit$pval))
  out
}

fit_epgls_ws_species <- function(long, species_tree, iter = 49L) {
  out <- list(status = "error", pval = NA_real_)
  dat <- long[, c("species", "y", "x_sum"), drop = FALSE]
  colnames(dat)[3] <- "x"
  dat <- dat[is.finite(dat$y) & is.finite(dat$x) & nzchar(dat$species), , drop = FALSE]
  if (nrow(dat) < 8 || length(unique(dat$species)) < 4) {
    out$status <- "no_data"
    return(out)
  }

  cov_obj <- build_obs_cov_from_species_tree(dat, species_tree)
  dat <- cov_obj[["data"]]
  cov_mat <- cov_obj[["Cov"]]
  if (is.null(cov_mat) || nrow(dat) < 8) {
    out$status <- "no_data"
    return(out)
  }
  dat$species_block <- factor(dat$species)
  if (length(unique(dat$species_block)) < 4) {
    out$status <- "no_data"
    return(out)
  }

  fit <- try({
    suppressWarnings(utils::capture.output({
      fit_obj <- RRPP::lm.rrpp.ws(
        y ~ x,
        data = dat,
        subjects = "species_block",
        Cov = cov_mat,
        iter = iter,
        print.progress = FALSE,
        Parallel = FALSE
      )
    }))
    fit_obj
  }, silent = TRUE)
  if (inherits(fit, "try-error")) {
    out$status <- "error"
    return(out)
  }
  sm <- summary(fit)
  tb <- sm[["table"]]
  if (is.null(tb) || nrow(tb) == 0) {
    out$status <- "error"
    return(out)
  }
  row_name <- if ("x" %in% rownames(tb)) "x" else rownames(tb)[1]
  out$pval <- suppressWarnings(as.numeric(tb[row_name, "Pr(>F)"]))
  out$status <- "ok"
  out
}

fit_pgls_species_mean <- function(wide_species, species_tree) {
  out <- list(status = "error", pval = NA_real_)
  dat <- wide_species[, c("species", "y", "mean_x_sum"), drop = FALSE]
  colnames(dat)[3] <- "x"
  dat <- dat[is.finite(dat$y) & is.finite(dat$x) & nzchar(dat$species), , drop = FALSE]
  if (nrow(dat) < 4 || length(unique(dat$species)) < 4) {
    out$status <- "no_data"
    return(out)
  }

  cov_species <- ape::vcv.phylo(species_tree)
  dat <- dat[dat$species %in% rownames(cov_species), , drop = FALSE]
  if (nrow(dat) < 4 || length(unique(dat$species)) < 4) {
    out$status <- "no_data"
    return(out)
  }
  dat <- dat[!duplicated(dat$species), , drop = FALSE]
  keep_tips <- intersect(species_tree$tip.label, dat$species)
  if (length(keep_tips) < 4) {
    out$status <- "no_data"
    return(out)
  }
  t_use <- species_tree
  drop_tips <- t_use$tip.label[!(t_use$tip.label %in% keep_tips)]
  if (length(drop_tips) > 0) t_use <- ape::drop.tip(t_use, drop_tips)
  dat <- dat[dat$species %in% t_use$tip.label, , drop = FALSE]
  if (nrow(dat) < 4 || length(unique(dat$species)) < 4) {
    out$status <- "no_data"
    return(out)
  }
  dat$species <- factor(dat$species, levels = t_use$tip.label)

  fit <- try(nlme::gls(
    y ~ x,
    data = dat,
    correlation = ape::corBrownian(phy = t_use, form = ~species),
    method = "ML"
  ), silent = TRUE)
  if (inherits(fit, "try-error")) {
    out$status <- "error"
    return(out)
  }
  sm <- summary(fit)
  tb <- sm[["tTable"]]
  if (is.null(tb) || nrow(tb) == 0 || !("x" %in% rownames(tb))) {
    out$status <- "error"
    return(out)
  }
  out$pval <- suppressWarnings(as.numeric(tb["x", "p-value"]))
  out$status <- "ok"
  out
}

fit_hierarchical_species_prototype <- function(long, species_tree) {
  out <- list(status = "error", pval = NA_real_)
  dat <- long[, c("species", "y", "x_sum"), drop = FALSE]
  colnames(dat)[3] <- "x"
  dat <- dat[is.finite(dat$y) & is.finite(dat$x) & nzchar(dat$species), , drop = FALSE]
  if (nrow(dat) < 8 || length(unique(dat$species)) < 4) {
    out$status <- "no_data"
    return(out)
  }
  dat$species <- factor(dat$species)
  species_levels <- levels(dat$species)

  cov_species <- ape::vcv.phylo(species_tree)
  if (!all(species_levels %in% rownames(cov_species))) {
    keep <- dat$species %in% rownames(cov_species)
    dat <- dat[keep, , drop = FALSE]
    dat$species <- droplevels(dat$species)
    species_levels <- levels(dat$species)
  }
  if (nrow(dat) < 8 || length(species_levels) < 4) {
    out$status <- "no_data"
    return(out)
  }

  y <- as.numeric(scale(dat$y))
  x <- as.numeric(scale(dat$x))
  if (!all(is.finite(y)) || !all(is.finite(x)) || stats::sd(x) <= 0) {
    out$status <- "no_data"
    return(out)
  }

  X <- cbind(Intercept = 1, x = x)
  Zs <- model.matrix(~0 + dat$species)
  colnames(Zs) <- species_levels

  Cs <- cov_species[species_levels, species_levels, drop = FALSE]
  eps <- suppressWarnings(max(diag(Cs), na.rm = TRUE)) * 1e-8
  if (!is.finite(eps) || eps <= 0) eps <- 1e-8
  diag(Cs) <- diag(Cs) + eps

  A_phy <- Zs %*% Cs %*% t(Zs)
  A_species <- tcrossprod(Zs)
  I_n <- diag(nrow(dat))

  solve_chol <- function(R, b) {
    backsolve(R, forwardsolve(t(R), b))
  }

  nll_fun <- function(par) {
    v_phy <- exp(par[1])
    v_species <- exp(par[2])
    v_resid <- exp(par[3])
    V <- v_phy * A_phy + v_species * A_species + v_resid * I_n
    R <- try(chol(V), silent = TRUE)
    if (inherits(R, "try-error")) return(1e30)
    VinvX <- solve_chol(R, X)
    XtVinvX <- crossprod(X, VinvX)
    bad <- suppressWarnings(kappa(XtVinvX, exact = FALSE))
    if (!is.finite(bad) || bad > 1e12) return(1e30)
    Vinvy <- solve_chol(R, y)
    beta_hat <- try(solve(XtVinvX, crossprod(X, Vinvy)), silent = TRUE)
    if (inherits(beta_hat, "try-error")) return(1e30)
    resid <- y - as.numeric(X %*% beta_hat)
    Vinv_resid <- solve_chol(R, resid)
    quad <- sum(resid * Vinv_resid)
    if (!is.finite(quad)) return(1e30)
    logdet <- 2 * sum(log(diag(R)))
    0.5 * (length(y) * log(2 * pi) + logdet + quad)
  }

  init_scale <- stats::var(y)
  if (!is.finite(init_scale) || init_scale <= 0) init_scale <- 1
  start_list <- list(
    log(c(init_scale * 0.3, init_scale * 0.3, init_scale * 0.4)),
    log(c(0.5, 0.5, 0.5)),
    log(c(1.0, 0.2, 0.8))
  )

  best <- NULL
  best_value <- Inf
  for (st in start_list) {
    fit <- try(optim(
      par = st,
      fn = nll_fun,
      method = "L-BFGS-B",
      lower = rep(log(1e-8), 3),
      upper = rep(log(1e3), 3),
      control = list(maxit = 400)
    ), silent = TRUE)
    if (!inherits(fit, "try-error") && is.finite(fit$value) && fit$value < best_value) {
      best <- fit
      best_value <- fit$value
    }
  }
  if (is.null(best)) {
    out$status <- "error"
    return(out)
  }

  v_phy <- exp(best$par[1])
  v_species <- exp(best$par[2])
  v_resid <- exp(best$par[3])
  V <- v_phy * A_phy + v_species * A_species + v_resid * I_n
  R <- try(chol(V), silent = TRUE)
  if (inherits(R, "try-error")) {
    out$status <- "error"
    return(out)
  }
  VinvX <- solve_chol(R, X)
  XtVinvX <- crossprod(X, VinvX)
  Vinvy <- solve_chol(R, y)
  beta_hat <- try(solve(XtVinvX, crossprod(X, Vinvy)), silent = TRUE)
  cov_beta <- try(solve(XtVinvX), silent = TRUE)
  if (inherits(beta_hat, "try-error") || inherits(cov_beta, "try-error")) {
    out$status <- "error"
    return(out)
  }
  se <- suppressWarnings(sqrt(diag(cov_beta)))
  if (length(se) < 2 || !is.finite(se[2]) || se[2] <= 0) {
    out$status <- "error"
    return(out)
  }
  z <- as.numeric(beta_hat[2] / se[2])
  out$pval <- 2 * stats::pnorm(-abs(z))
  out$status <- "ok"
  out
}

scenarios <- data.frame(
  scenario = c("A_balanced", "B_missing_mid", "C_missing_high", "D_few_species"),
  n_species = c(12, 12, 12, 8),
  dup_fraction = c(0.35, 0.45, 0.55, 0.50),
  max_copy = c(4, 4, 5, 3),
  min_rep = c(3, 2, 2, 2),
  max_rep = c(5, 5, 4, 3),
  miss_prob = c(0.10, 0.30, 0.55, 0.40),
  sigma_rep = c(0.40, 0.55, 0.75, 0.65),
  stringsAsFactors = FALSE
)

truth_states <- data.frame(
  truth = c("null", "alt"),
  beta = c(0.0, 0.7),
  stringsAsFactors = FALSE
)

set.seed(42)
rows <- list()
idx <- 1L
rphy_method_label <- if (rphylopars_mode == "phenocov") "rphylopars_species_tree_phenocov" else "rphylopars_species_tree"
for (s in seq_len(nrow(scenarios))) {
  sc <- scenarios[s, ]
  for (t in seq_len(nrow(truth_states))) {
    tr <- truth_states[t, ]
    cat("Running scenario:", sc$scenario, "| truth:", tr$truth, "\n")
    for (i in seq_len(n_iter)) {
      d <- simulate_species_tree_sum_dataset(
        n_species = sc$n_species,
        dup_fraction = sc$dup_fraction,
        max_copy = sc$max_copy,
        min_rep = sc$min_rep,
        max_rep = sc$max_rep,
        miss_prob = sc$miss_prob,
        beta = tr$beta,
        sigma_rep = sc$sigma_rep
      )

      t0 <- proc.time()[["elapsed"]]
      m_rphy <- if (rphylopars_mode == "phenocov") {
        fit_rphylopars_species_tree_phenocov(d$long, d$species_tree)
      } else {
        fit_rphylopars_species_tree(d$wide_species, d$species_tree)
      }
      rt_rphy <- proc.time()[["elapsed"]] - t0

      t0 <- proc.time()[["elapsed"]]
      m_epws <- fit_epgls_ws_species(d$long, d$species_tree, iter = rrpp_iter)
      rt_epws <- proc.time()[["elapsed"]] - t0

      t0 <- proc.time()[["elapsed"]]
      m_pgls <- fit_pgls_species_mean(d$wide_species, d$species_tree)
      rt_pgls <- proc.time()[["elapsed"]] - t0

      t0 <- proc.time()[["elapsed"]]
      m_hie <- fit_hierarchical_species_prototype(d$long, d$species_tree)
      rt_hie <- proc.time()[["elapsed"]] - t0

      rows[[idx]] <- data.frame(
        scenario = sc$scenario,
        truth = tr$truth,
        truth_positive = as.integer(tr$truth == "alt"),
        iter = i,
        method = rphy_method_label,
        pval = m_rphy$pval,
        status = m_rphy$status,
        runtime_sec = rt_rphy,
        stringsAsFactors = FALSE
      ); idx <- idx + 1L
      rows[[idx]] <- data.frame(
        scenario = sc$scenario,
        truth = tr$truth,
        truth_positive = as.integer(tr$truth == "alt"),
        iter = i,
        method = "epgls_ws_species_sum",
        pval = m_epws$pval,
        status = m_epws$status,
        runtime_sec = rt_epws,
        stringsAsFactors = FALSE
      ); idx <- idx + 1L
      rows[[idx]] <- data.frame(
        scenario = sc$scenario,
        truth = tr$truth,
        truth_positive = as.integer(tr$truth == "alt"),
        iter = i,
        method = "pgls_species_mean_sum",
        pval = m_pgls$pval,
        status = m_pgls$status,
        runtime_sec = rt_pgls,
        stringsAsFactors = FALSE
      ); idx <- idx + 1L
      rows[[idx]] <- data.frame(
        scenario = sc$scenario,
        truth = tr$truth,
        truth_positive = as.integer(tr$truth == "alt"),
        iter = i,
        method = "hierarchical_species_tree",
        pval = m_hie$pval,
        status = m_hie$status,
        runtime_sec = rt_hie,
        stringsAsFactors = FALSE
      ); idx <- idx + 1L
    }
  }
}

df <- do.call(rbind, rows)
df$pred_positive <- as.integer(is.finite(df$pval) & (df$pval < alpha))
df$fit_ok <- as.integer(df$status == "ok")

out_runs <- file.path(out_dir, "benchmark_4methods_speciesTree_sumexpr_confusion_runs.tsv")
write.table(df, out_runs, sep = "\t", row.names = FALSE, quote = FALSE)

calc_confusion <- function(dat) {
  tp <- sum(dat$truth_positive == 1 & dat$pred_positive == 1)
  fn <- sum(dat$truth_positive == 1 & dat$pred_positive == 0)
  tn <- sum(dat$truth_positive == 0 & dat$pred_positive == 0)
  fp <- sum(dat$truth_positive == 0 & dat$pred_positive == 1)
  tpr <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  fnr <- if ((tp + fn) > 0) fn / (tp + fn) else NA_real_
  tnr <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  fpr <- if ((tn + fp) > 0) fp / (tn + fp) else NA_real_
  data.frame(
    TP = tp, FP = fp, TN = tn, FN = fn,
    TPR = tpr, FPR = fpr, TNR = tnr, FNR = fnr,
    n = nrow(dat),
    ok_rate = mean(dat$fit_ok == 1),
    stringsAsFactors = FALSE
  )
}

make_summary <- function(dat, scope_name) {
  groups <- split(dat, list(dat$scenario, dat$method), drop = TRUE)
  rows <- lapply(groups, function(gd) {
    cbind(
      data.frame(
        scenario = unique(gd$scenario),
        method = unique(gd$method),
        scope = scope_name,
        stringsAsFactors = FALSE
      ),
      calc_confusion(gd)
    )
  })
  do.call(rbind, rows)
}

summary_all <- make_summary(df, "all_runs")
summary_ok <- make_summary(df[df$fit_ok == 1, , drop = FALSE], "ok_only")
summary_tab <- rbind(summary_all, summary_ok)
out_summary <- file.path(out_dir, "benchmark_4methods_speciesTree_sumexpr_confusion_summary.tsv")
write.table(summary_tab, out_summary, sep = "\t", row.names = FALSE, quote = FALSE)

make_overall <- function(dat, scope_name) {
  groups <- split(dat, dat$method, drop = TRUE)
  rows <- lapply(groups, function(gd) {
    cbind(
      data.frame(method = unique(gd$method), scope = scope_name, stringsAsFactors = FALSE),
      calc_confusion(gd)
    )
  })
  do.call(rbind, rows)
}

overall_tab <- rbind(
  make_overall(df, "overall_all_runs"),
  make_overall(df[df$fit_ok == 1, , drop = FALSE], "overall_ok_only")
)
out_overall <- file.path(out_dir, "benchmark_4methods_speciesTree_sumexpr_confusion_overall.tsv")
write.table(overall_tab, out_overall, sep = "\t", row.names = FALSE, quote = FALSE)

plot_df <- overall_tab
plot_df <- plot_df[plot_df$scope == "overall_all_runs", , drop = FALSE]
plot_long <- rbind(
  data.frame(method = plot_df$method, metric = "TPR", value = plot_df$TPR, stringsAsFactors = FALSE),
  data.frame(method = plot_df$method, metric = "FPR", value = plot_df$FPR, stringsAsFactors = FALSE),
  data.frame(method = plot_df$method, metric = "TNR", value = plot_df$TNR, stringsAsFactors = FALSE),
  data.frame(method = plot_df$method, metric = "FNR", value = plot_df$FNR, stringsAsFactors = FALSE)
)
g <- ggplot(plot_long, aes(x = method, y = value, fill = metric)) +
  geom_col(position = "dodge") +
  ylim(0, 1) +
  labs(
    x = "Method",
    y = "Rate",
    title = sprintf("Species tree (sum expression): confusion rates, alpha=%.3f", alpha)
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))
out_plot_pdf <- file.path(out_dir, "benchmark_4methods_speciesTree_sumexpr_confusion_rates.pdf")
out_plot_png <- file.path(out_dir, "benchmark_4methods_speciesTree_sumexpr_confusion_rates.png")
ggsave(out_plot_pdf, g, width = 11.5, height = 5.8)
ggsave(out_plot_png, g, width = 11.5, height = 5.8, dpi = 200)

cat("runs:", out_runs, "\n")
cat("summary:", out_summary, "\n")
cat("overall:", out_overall, "\n")
cat("rates plot pdf:", out_plot_pdf, "\n")
cat("rates plot png:", out_plot_png, "\n")
cat("n_iter_per_truth:", n_iter, "\n")
cat("rrpp_iter:", rrpp_iter, "\n")
cat("alpha:", alpha, "\n")
cat("rphylopars_mode:", rphylopars_mode, "\n")
