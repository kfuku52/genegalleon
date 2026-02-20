#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[[1]] else file.path("workspace", "output", "benchmark_epgls")
n_iter <- if (length(args) >= 2) suppressWarnings(as.integer(args[[2]])) else 60L
rrpp_iter <- if (length(args) >= 3) suppressWarnings(as.integer(args[[3]])) else 49L
alpha <- if (length(args) >= 4) suppressWarnings(as.numeric(args[[4]])) else 0.05
epgls_formula_mode <- if (length(args) >= 5) tolower(args[[5]]) else "subject_x"
cov_diag_scale <- if (length(args) >= 6) suppressWarnings(as.numeric(args[[6]])) else 1e-8
epgls_delta <- if (length(args) >= 7) suppressWarnings(as.numeric(args[[7]])) else 0.001
epgls_gamma <- if (length(args) >= 8) tolower(args[[8]]) else "sample"
epgls_cov_mode <- if (length(args) >= 9) tolower(args[[9]]) else "obs"
if (is.na(n_iter) || n_iter < 10) n_iter <- 60L
if (is.na(rrpp_iter) || rrpp_iter < 9) rrpp_iter <- 49L
if (is.na(alpha) || alpha <= 0 || alpha >= 1) alpha <- 0.05
if (!(epgls_formula_mode %in% c("subject_x", "x_only"))) epgls_formula_mode <- "subject_x"
if (is.na(cov_diag_scale) || cov_diag_scale <= 0) cov_diag_scale <- 1e-8
if (is.na(epgls_delta) || epgls_delta < 0 || epgls_delta > 1) epgls_delta <- 0.001
if (!(epgls_gamma %in% c("sample", "equal"))) epgls_gamma <- "sample"
if (!(epgls_cov_mode %in% c("obs", "subject"))) epgls_cov_mode <- "obs"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ape)
  library(Rphylopars)
  library(ggplot2)
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

build_obs_cov_from_gene_tree <- function(dat, gene_tree, diag_scale = 1e-8) {
  if (nrow(dat) == 0) return(list(data = dat, Cov = NULL))
  cov_gene <- ape::vcv.phylo(gene_tree)
  tip_ids <- rownames(cov_gene)
  dat <- dat[dat[["gene_id"]] %in% tip_ids, , drop = FALSE]
  if (nrow(dat) == 0) return(list(data = dat, Cov = NULL))
  idx <- match(dat[["gene_id"]], tip_ids)
  cov_obs <- cov_gene[idx, idx, drop = FALSE]
  eps <- suppressWarnings(max(diag(cov_obs), na.rm = TRUE)) * diag_scale
  if (!is.finite(eps) || eps <= 0) eps <- 1e-8
  diag(cov_obs) <- diag(cov_obs) + eps
  obs_ids <- paste0("obs_", seq_len(nrow(dat)))
  rownames(cov_obs) <- obs_ids
  colnames(cov_obs) <- obs_ids
  list(data = dat, Cov = cov_obs)
}

build_subject_cov_from_gene_tree <- function(dat, gene_tree, ws_subject = c("species", "gene_id"), diag_scale = 1e-8) {
  ws_subject <- match.arg(ws_subject)
  cov_gene <- ape::vcv.phylo(gene_tree)
  if (ws_subject == "gene_id") {
    subjects <- unique(dat$gene_id)
    subjects <- subjects[subjects %in% rownames(cov_gene)]
    if (length(subjects) < 3) return(NULL)
    cov_sub <- cov_gene[subjects, subjects, drop = FALSE]
  } else {
    tip_species <- sub("_copy[0-9]+$", "", rownames(cov_gene))
    subjects <- unique(dat$species)
    subjects <- subjects[subjects %in% unique(tip_species)]
    if (length(subjects) < 3) return(NULL)
    cov_sub <- matrix(
      0,
      nrow = length(subjects),
      ncol = length(subjects),
      dimnames = list(subjects, subjects)
    )
    for (i in seq_along(subjects)) {
      tips_i <- names(tip_species)[tip_species == subjects[i]]
      if (length(tips_i) == 0) next
      for (j in i:length(subjects)) {
        tips_j <- names(tip_species)[tip_species == subjects[j]]
        if (length(tips_j) == 0) next
        v <- mean(cov_gene[tips_i, tips_j, drop = FALSE], na.rm = TRUE)
        cov_sub[i, j] <- v
        cov_sub[j, i] <- v
      }
    }
    if (any(!is.finite(cov_sub))) return(NULL)
  }
  eps <- suppressWarnings(max(diag(cov_sub), na.rm = TRUE)) * diag_scale
  if (!is.finite(eps) || eps <= 0) eps <- 1e-8
  diag(cov_sub) <- diag(cov_sub) + eps
  cov_sub
}

build_gene_tree_with_lineage_dup <- function(base_tree, max_copy, dup_fraction) {
  gt <- base_tree
  sp_names <- base_tree$tip.label
  dup_count <- rep(1L, length(sp_names))
  names(dup_count) <- sp_names
  n_dup <- max(1L, floor(length(sp_names) * dup_fraction))
  dup_species <- sample(sp_names, size = n_dup, replace = FALSE)
  for (sp in dup_species) {
    dup_count[sp] <- sample(2:max_copy, size = 1)
  }

  for (sp in sp_names) {
    idx <- which(gt$tip.label == sp)
    if (length(idx) != 1) next
    k <- dup_count[[sp]]
    if (k <= 1) {
      gt$tip.label[idx] <- sprintf("%s_copy01", sp)
      next
    }
    parent_edge <- which(gt$edge[, 2] == idx)
    tip_bl <- if (length(parent_edge) == 1) gt$edge.length[parent_edge] else NA_real_
    if (!is.finite(tip_bl) || tip_bl <= 0) {
      vv <- gt$edge.length[is.finite(gt$edge.length) & gt$edge.length > 0]
      tip_bl <- if (length(vv) > 0) stats::median(vv) else 1
    }
    sub <- ape::rtree(k, tip.label = sprintf("%s_copy%02d", sp, seq_len(k)))
    sub$edge.length <- pmax(sub$edge.length, 1e-6) * tip_bl * 0.08
    gt <- ape::bind.tree(gt, sub, where = idx, position = 0)
    if (sp %in% gt$tip.label) gt <- ape::drop.tip(gt, sp)
  }
  gt
}

simulate_gene_tree_dataset <- function(n_species = 12,
                                       dup_fraction = 0.45,
                                       max_copy = 4,
                                       min_rep = 2,
                                       max_rep = 5,
                                       miss_prob = 0.25,
                                       beta = 0.7,
                                       sigma_species_latent = 1.0,
                                       sigma_species_noise = 1.0,
                                       sigma_gene_phy = 0.8,
                                       sigma_rep = 0.5) {
  species <- sprintf("sp%02d", seq_len(n_species))
  species_tree <- ape::rtree(n_species)
  species_tree$tip.label <- species

  z_sp <- as.numeric(ape::rTraitCont(species_tree, model = "BM", sigma = sigma_species_latent))
  names(z_sp) <- species_tree$tip.label
  y_noise <- as.numeric(ape::rTraitCont(species_tree, model = "BM", sigma = sigma_species_noise))
  names(y_noise) <- species_tree$tip.label
  y_sp <- beta * z_sp + y_noise
  y_sp <- as.numeric(y_sp[species])
  names(y_sp) <- species

  gene_tree <- build_gene_tree_with_lineage_dup(
    base_tree = species_tree,
    max_copy = max_copy,
    dup_fraction = dup_fraction
  )
  gene_ids <- gene_tree$tip.label
  gene_species <- sub("_copy[0-9]+$", "", gene_ids)

  gene_phy <- as.numeric(ape::rTraitCont(gene_tree, model = "BM", sigma = sigma_gene_phy))
  names(gene_phy) <- gene_tree$tip.label
  gene_phy <- as.numeric(scale(gene_phy))
  names(gene_phy) <- gene_tree$tip.label

  max_r <- max_rep
  wide <- data.frame(
    species = gene_ids,
    y = as.numeric(y_sp[gene_species]),
    stringsAsFactors = FALSE
  )
  for (r in seq_len(max_r)) wide[[sprintf("x_%d", r)]] <- NA_real_

  long <- data.frame(
    gene_id = character(0),
    species = character(0),
    y = numeric(0),
    x = numeric(0),
    stringsAsFactors = FALSE
  )
  for (g in seq_along(gene_ids)) {
    gid <- gene_ids[g]
    sp <- gene_species[g]
    x_gene <- z_sp[sp] + gene_phy[gid]
    k_rep <- sample(min_rep:max_rep, 1)
    x_rep <- x_gene + stats::rnorm(k_rep, mean = 0, sd = sigma_rep)
    keep <- stats::runif(k_rep) > miss_prob
    x_rep[!keep] <- NA_real_
    if (all(is.na(x_rep))) x_rep[1] <- x_gene + stats::rnorm(1, 0, sigma_rep)

    wide[g, sprintf("x_%d", seq_len(k_rep))] <- as.numeric(x_rep)

    long_g <- data.frame(
      gene_id = rep(gid, k_rep),
      species = rep(sp, k_rep),
      y = rep(y_sp[sp], k_rep),
      x = as.numeric(x_rep),
      stringsAsFactors = FALSE
    )
    long <- rbind(long, long_g)
  }

  x_cols <- grep("^x_[0-9]+$", colnames(wide), value = TRUE)
  wide$mean_x <- apply(wide[, x_cols, drop = FALSE], 1, function(v) {
    vv <- suppressWarnings(as.numeric(v))
    if (all(is.na(vv))) return(NA_real_)
    mean(vv, na.rm = TRUE)
  })

  long <- long[is.finite(long$x) & is.finite(long$y) & nzchar(long$species) & nzchar(long$gene_id), , drop = FALSE]
  rownames(long) <- NULL

  list(gene_tree = gene_tree, wide_gene = wide, long = long)
}

fit_rphylopars_gene_tree <- function(wide_gene, gene_tree) {
  out <- list(status = "error", pval = NA_real_)
  if (!("mean_x" %in% colnames(wide_gene))) {
    out$status <- "no_data"
    return(out)
  }
  td <- wide_gene[, c("species", "y", "mean_x"), drop = FALSE]
  colnames(td)[3] <- "x"
  td$y <- suppressWarnings(as.numeric(td$y))
  td$x <- suppressWarnings(as.numeric(td$x))
  td <- td[is.finite(td$y) & is.finite(td$x) & nzchar(td$species), , drop = FALSE]
  if (nrow(td) < 4 || length(unique(td$species)) < 4) {
    out$status <- "no_data"
    return(out)
  }
  keep_tips <- intersect(gene_tree$tip.label, td$species)
  if (length(keep_tips) < 4) {
    out$status <- "no_data"
    return(out)
  }
  t_use <- gene_tree
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

fit_epgls <- function(long, gene_tree, ws_subject = c("species", "gene_id"), iter = 49L, formula_mode = c("subject_x", "x_only"), cov_diag_scale = 1e-8, delta = 0.001, gamma = c("sample", "equal"), cov_mode = c("obs", "subject")) {
  ws_subject <- match.arg(ws_subject)
  formula_mode <- match.arg(formula_mode)
  gamma <- match.arg(gamma)
  cov_mode <- match.arg(cov_mode)
  out <- list(status = "error", pval = NA_real_)
  dat <- long[, c("gene_id", "species", "y", "x"), drop = FALSE]
  dat <- dat[is.finite(dat$y) & is.finite(dat$x) & nzchar(dat$species) & nzchar(dat$gene_id), , drop = FALSE]
  if (nrow(dat) < 6 || length(unique(dat$species)) < 3 || length(unique(dat$gene_id)) < 3 || length(unique(dat[[ws_subject]])) < 3) {
    out$status <- "no_data"
    return(out)
  }

  cov_mat <- NULL
  if (cov_mode == "subject") {
    cov_mat <- build_subject_cov_from_gene_tree(dat, gene_tree, ws_subject = ws_subject, diag_scale = cov_diag_scale)
    if (!is.null(cov_mat)) {
      keep <- dat[[ws_subject]] %in% rownames(cov_mat)
      dat <- dat[keep, , drop = FALSE]
    }
  }
  if (is.null(cov_mat)) {
    cov_obj <- build_obs_cov_from_gene_tree(dat, gene_tree, diag_scale = cov_diag_scale)
    dat <- cov_obj[["data"]]
    cov_mat <- cov_obj[["Cov"]]
  }
  if (is.null(cov_mat) || nrow(dat) < 6) {
    out$status <- "no_data"
    return(out)
  }
  dat$subject_block <- factor(dat[[ws_subject]])
  if (length(unique(dat$subject_block)) < 3) {
    out$status <- "no_data"
    return(out)
  }

  fit <- try({
    suppressWarnings(utils::capture.output({
      if (formula_mode == "subject_x") {
        fit_obj <- RRPP::lm.rrpp.ws(
          y ~ subject_block + x,
          data = dat,
          subjects = "subject_block",
          Cov = cov_mat,
          delta = delta,
          gamma = gamma,
          iter = iter,
          print.progress = FALSE,
          Parallel = FALSE
        )
      } else {
        fit_obj <- RRPP::lm.rrpp.ws(
          y ~ x,
          data = dat,
          subjects = "subject_block",
          Cov = cov_mat,
          delta = delta,
          gamma = gamma,
          iter = iter,
          print.progress = FALSE,
          Parallel = FALSE
        )
      }
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

fit_hierarchical_prototype <- function(long, gene_tree) {
  out <- list(status = "error", pval = NA_real_)
  dat <- long[, c("gene_id", "species", "y", "x"), drop = FALSE]
  dat <- dat[is.finite(dat$y) & is.finite(dat$x) & nzchar(dat$species) & nzchar(dat$gene_id), , drop = FALSE]
  if (nrow(dat) < 8 || length(unique(dat$species)) < 3 || length(unique(dat$gene_id)) < 3) {
    out$status <- "no_data"
    return(out)
  }
  dat$species <- factor(dat$species)
  dat$gene_id <- factor(dat$gene_id)

  cov_gene_full <- ape::vcv.phylo(gene_tree)
  keep <- dat$gene_id %in% rownames(cov_gene_full)
  dat <- dat[keep, , drop = FALSE]
  dat$species <- droplevels(dat$species)
  dat$gene_id <- droplevels(dat$gene_id)
  if (nrow(dat) < 8 || nlevels(dat$species) < 3 || nlevels(dat$gene_id) < 3) {
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
  species_levels <- levels(dat$species)
  gene_levels <- levels(dat$gene_id)
  Zs <- model.matrix(~0 + dat$species)
  Zg <- model.matrix(~0 + dat$gene_id)
  colnames(Zs) <- species_levels
  colnames(Zg) <- gene_levels
  Cg <- cov_gene_full[gene_levels, gene_levels, drop = FALSE]
  eps <- suppressWarnings(max(diag(Cg), na.rm = TRUE)) * 1e-8
  if (!is.finite(eps) || eps <= 0) eps <- 1e-8
  diag(Cg) <- diag(Cg) + eps

  A_species <- tcrossprod(Zs)
  A_gene <- Zg %*% Cg %*% t(Zg)
  I_n <- diag(nrow(dat))

  solve_chol <- function(R, b) {
    backsolve(R, forwardsolve(t(R), b))
  }

  nll_fun <- function(par) {
    v_species <- exp(par[1])
    v_gene <- exp(par[2])
    v_resid <- exp(par[3])
    V <- v_species * A_species + v_gene * A_gene + v_resid * I_n
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
    log(c(init_scale * 0.4, init_scale * 0.4, init_scale * 0.2)),
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

  v_species <- exp(best$par[1])
  v_gene <- exp(best$par[2])
  v_resid <- exp(best$par[3])
  V <- v_species * A_species + v_gene * A_gene + v_resid * I_n
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
for (s in seq_len(nrow(scenarios))) {
  sc <- scenarios[s, ]
  for (t in seq_len(nrow(truth_states))) {
    tr <- truth_states[t, ]
    cat("Running scenario:", sc$scenario, "| truth:", tr$truth, "\n")
    for (i in seq_len(n_iter)) {
      d <- simulate_gene_tree_dataset(
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
      m_rphy <- fit_rphylopars_gene_tree(d$wide_gene, d$gene_tree)
      rt_rphy <- proc.time()[["elapsed"]] - t0

      t0 <- proc.time()[["elapsed"]]
      m_eps <- fit_epgls(
        d$long, d$gene_tree,
        ws_subject = "species",
        iter = rrpp_iter,
        formula_mode = epgls_formula_mode,
        cov_diag_scale = cov_diag_scale,
        delta = epgls_delta,
        gamma = epgls_gamma,
        cov_mode = epgls_cov_mode
      )
      rt_eps <- proc.time()[["elapsed"]] - t0

      t0 <- proc.time()[["elapsed"]]
      m_epg <- fit_epgls(
        d$long, d$gene_tree,
        ws_subject = "gene_id",
        iter = rrpp_iter,
        formula_mode = epgls_formula_mode,
        cov_diag_scale = cov_diag_scale,
        delta = epgls_delta,
        gamma = epgls_gamma,
        cov_mode = epgls_cov_mode
      )
      rt_epg <- proc.time()[["elapsed"]] - t0

      t0 <- proc.time()[["elapsed"]]
      m_hie <- fit_hierarchical_prototype(d$long, d$gene_tree)
      rt_hie <- proc.time()[["elapsed"]] - t0

      rows[[idx]] <- data.frame(
        scenario = sc$scenario,
        truth = tr$truth,
        truth_positive = as.integer(tr$truth == "alt"),
        iter = i,
        method = "rphylopars_gene_tree",
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
        method = "epgls_ws_species",
        pval = m_eps$pval,
        status = m_eps$status,
        runtime_sec = rt_eps,
        stringsAsFactors = FALSE
      ); idx <- idx + 1L
      rows[[idx]] <- data.frame(
        scenario = sc$scenario,
        truth = tr$truth,
        truth_positive = as.integer(tr$truth == "alt"),
        iter = i,
        method = "epgls_ws_gene_id",
        pval = m_epg$pval,
        status = m_epg$status,
        runtime_sec = rt_epg,
        stringsAsFactors = FALSE
      ); idx <- idx + 1L
      rows[[idx]] <- data.frame(
        scenario = sc$scenario,
        truth = tr$truth,
        truth_positive = as.integer(tr$truth == "alt"),
        iter = i,
        method = "hierarchical_prototype",
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

out_runs <- file.path(out_dir, "benchmark_4methods_confusion_runs.tsv")
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
out_summary <- file.path(out_dir, "benchmark_4methods_confusion_summary.tsv")
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
out_overall <- file.path(out_dir, "benchmark_4methods_confusion_overall.tsv")
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
    title = sprintf("Confusion rates (lineage-specific duplication, alpha=%.3f)", alpha)
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))
out_plot_pdf <- file.path(out_dir, "benchmark_4methods_confusion_rates.pdf")
out_plot_png <- file.path(out_dir, "benchmark_4methods_confusion_rates.png")
ggsave(out_plot_pdf, g, width = 11, height = 5.5)
ggsave(out_plot_png, g, width = 11, height = 5.5, dpi = 200)

cat("runs:", out_runs, "\n")
cat("summary:", out_summary, "\n")
cat("overall:", out_overall, "\n")
cat("rates plot pdf:", out_plot_pdf, "\n")
cat("rates plot png:", out_plot_png, "\n")
cat("n_iter_per_truth:", n_iter, "\n")
cat("rrpp_iter:", rrpp_iter, "\n")
cat("alpha:", alpha, "\n")
cat("epgls_formula_mode:", epgls_formula_mode, "\n")
cat("cov_diag_scale:", cov_diag_scale, "\n")
cat("epgls_delta:", epgls_delta, "\n")
cat("epgls_gamma:", epgls_gamma, "\n")
cat("epgls_cov_mode:", epgls_cov_mode, "\n")
