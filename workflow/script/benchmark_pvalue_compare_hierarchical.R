#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[[1]] else file.path("workspace", "output", "benchmark_epgls")
n_iter <- if (length(args) >= 2) suppressWarnings(as.integer(args[[2]])) else 120L
rrpp_iter <- if (length(args) >= 3) suppressWarnings(as.integer(args[[3]])) else 49L
if (is.na(n_iter) || n_iter < 10) n_iter <- 120L
if (is.na(rrpp_iter) || rrpp_iter < 9) rrpp_iter <- 49L
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

build_obs_cov_from_gene_tree <- function(dat, gene_tree) {
  if (nrow(dat) == 0) return(list(data = dat, Cov = NULL))
  cov_gene <- ape::vcv.phylo(gene_tree)
  tip_ids <- rownames(cov_gene)
  dat <- dat[dat[["gene_id"]] %in% tip_ids, , drop = FALSE]
  if (nrow(dat) == 0) return(list(data = dat, Cov = NULL))
  idx <- match(dat[["gene_id"]], tip_ids)
  cov_obs <- cov_gene[idx, idx, drop = FALSE]
  eps <- suppressWarnings(max(diag(cov_obs), na.rm = TRUE)) * 1e-8
  if (!is.finite(eps) || eps <= 0) eps <- 1e-8
  diag(cov_obs) <- diag(cov_obs) + eps
  obs_ids <- paste0("obs_", seq_len(nrow(dat)))
  rownames(cov_obs) <- obs_ids
  colnames(cov_obs) <- obs_ids
  list(data = dat, Cov = cov_obs)
}

simulate_gene_tree_dataset <- function(n_species = 12,
                                       dup_fraction = 0.4,
                                       max_copy = 4,
                                       min_rep = 2,
                                       max_rep = 5,
                                       miss_prob = 0.2,
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
      } else {
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
    }
    gt
  }

  gene_tree <- build_gene_tree_with_lineage_dup(species_tree, max_copy = max_copy, dup_fraction = dup_fraction)
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

  list(
    gene_tree = gene_tree,
    species_tree = species_tree,
    wide_gene = wide,
    long = long
  )
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

fit_epgls <- function(long, gene_tree, ws_subject = c("species", "gene_id"), iter = 49L) {
  ws_subject <- match.arg(ws_subject)
  out <- list(status = "error", pval = NA_real_)
  dat <- long[, c("gene_id", "species", "y", "x"), drop = FALSE]
  dat <- dat[is.finite(dat$y) & is.finite(dat$x) & nzchar(dat$species) & nzchar(dat$gene_id), , drop = FALSE]
  if (nrow(dat) < 6 || length(unique(dat$species)) < 3 || length(unique(dat$gene_id)) < 3 || length(unique(dat[[ws_subject]])) < 3) {
    out$status <- "no_data"
    return(out)
  }

  cov_obj <- build_obs_cov_from_gene_tree(dat, gene_tree)
  dat <- cov_obj[["data"]]
  cov_mat <- cov_obj[["Cov"]]
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
      fit_obj <- RRPP::lm.rrpp.ws(
        y ~ subject_block + x,
        data = dat,
        subjects = "subject_block",
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

fit_hierarchical_prototype <- function(long, gene_tree) {
  out <- list(status = "error", pval = NA_real_, beta = NA_real_, se = NA_real_)
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
    rcond_x <- suppressWarnings(kappa(XtVinvX, exact = FALSE))
    if (!is.finite(rcond_x) || rcond_x > 1e12) return(1e30)

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
  pval <- 2 * stats::pnorm(-abs(z))
  out$status <- "ok"
  out$pval <- pval
  out$beta <- as.numeric(beta_hat[2])
  out$se <- as.numeric(se[2])
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

set.seed(42)
rows <- list()
idx <- 1L

for (s in seq_len(nrow(scenarios))) {
  sc <- scenarios[s, ]
  cat("Running scenario:", sc$scenario, "\n")
  for (i in seq_len(n_iter)) {
    d <- simulate_gene_tree_dataset(
      n_species = sc$n_species,
      dup_fraction = sc$dup_fraction,
      max_copy = sc$max_copy,
      min_rep = sc$min_rep,
      max_rep = sc$max_rep,
      miss_prob = sc$miss_prob,
      sigma_rep = sc$sigma_rep
    )

    t0 <- proc.time()[["elapsed"]]
    m_rphy <- fit_rphylopars_gene_tree(d$wide_gene, d$gene_tree)
    rt_rphy <- proc.time()[["elapsed"]] - t0

    t0 <- proc.time()[["elapsed"]]
    m_epsp <- fit_epgls(d$long, d$gene_tree, ws_subject = "species", iter = rrpp_iter)
    rt_epsp <- proc.time()[["elapsed"]] - t0

    t0 <- proc.time()[["elapsed"]]
    m_epgid <- fit_epgls(d$long, d$gene_tree, ws_subject = "gene_id", iter = rrpp_iter)
    rt_epgid <- proc.time()[["elapsed"]] - t0

    t0 <- proc.time()[["elapsed"]]
    m_hier <- fit_hierarchical_prototype(d$long, d$gene_tree)
    rt_hier <- proc.time()[["elapsed"]] - t0

    rows[[idx]] <- data.frame(
      scenario = sc$scenario,
      iter = i,
      p_rphylopars_gene_tree = m_rphy$pval,
      status_rphylopars_gene_tree = m_rphy$status,
      runtime_rphylopars_gene_tree = rt_rphy,
      p_epgls_species = m_epsp$pval,
      status_epgls_species = m_epsp$status,
      runtime_epgls_species = rt_epsp,
      p_epgls_gene_id = m_epgid$pval,
      status_epgls_gene_id = m_epgid$status,
      runtime_epgls_gene_id = rt_epgid,
      p_hierarchical = m_hier$pval,
      status_hierarchical = m_hier$status,
      runtime_hierarchical = rt_hier,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }
}

df <- do.call(rbind, rows)

out_pval <- file.path(out_dir, "benchmark_4methods_pvalues.tsv")
write.table(df, out_pval, sep = "\t", row.names = FALSE, quote = FALSE)

status_long <- rbind(
  data.frame(scenario = df$scenario, iter = df$iter, method = "rphylopars_gene_tree", status = df$status_rphylopars_gene_tree, stringsAsFactors = FALSE),
  data.frame(scenario = df$scenario, iter = df$iter, method = "epgls_ws_species", status = df$status_epgls_species, stringsAsFactors = FALSE),
  data.frame(scenario = df$scenario, iter = df$iter, method = "epgls_ws_gene_id", status = df$status_epgls_gene_id, stringsAsFactors = FALSE),
  data.frame(scenario = df$scenario, iter = df$iter, method = "hierarchical_prototype", status = df$status_hierarchical, stringsAsFactors = FALSE)
)
status_summary <- aggregate(
  status ~ scenario + method,
  data = status_long,
  FUN = function(x) sprintf("%d/%d (%.1f%%)", sum(x == "ok"), length(x), 100 * mean(x == "ok"))
)
colnames(status_summary)[3] <- "success_rate"
out_summary <- file.path(out_dir, "benchmark_4methods_summary.tsv")
write.table(status_summary, out_summary, sep = "\t", row.names = FALSE, quote = FALSE)

status_count <- as.data.frame(table(status_long$scenario, status_long$method, status_long$status), stringsAsFactors = FALSE)
colnames(status_count) <- c("scenario", "method", "status", "count")
status_count <- status_count[status_count$count > 0, , drop = FALSE]
out_status <- file.path(out_dir, "benchmark_4methods_status_counts.tsv")
write.table(status_count, out_status, sep = "\t", row.names = FALSE, quote = FALSE)

runtime_long <- rbind(
  data.frame(
    scenario = df$scenario,
    iter = df$iter,
    method = "rphylopars_gene_tree",
    status = df$status_rphylopars_gene_tree,
    runtime_sec = suppressWarnings(as.numeric(df$runtime_rphylopars_gene_tree)),
    stringsAsFactors = FALSE
  ),
  data.frame(
    scenario = df$scenario,
    iter = df$iter,
    method = "epgls_ws_species",
    status = df$status_epgls_species,
    runtime_sec = suppressWarnings(as.numeric(df$runtime_epgls_species)),
    stringsAsFactors = FALSE
  ),
  data.frame(
    scenario = df$scenario,
    iter = df$iter,
    method = "epgls_ws_gene_id",
    status = df$status_epgls_gene_id,
    runtime_sec = suppressWarnings(as.numeric(df$runtime_epgls_gene_id)),
    stringsAsFactors = FALSE
  ),
  data.frame(
    scenario = df$scenario,
    iter = df$iter,
    method = "hierarchical_prototype",
    status = df$status_hierarchical,
    runtime_sec = suppressWarnings(as.numeric(df$runtime_hierarchical)),
    stringsAsFactors = FALSE
  )
)
runtime_long <- runtime_long[is.finite(runtime_long$runtime_sec) & runtime_long$runtime_sec >= 0, , drop = FALSE]
out_runtime <- file.path(out_dir, "benchmark_4methods_runtime.tsv")
write.table(runtime_long, out_runtime, sep = "\t", row.names = FALSE, quote = FALSE)

summarize_runtime <- function(v) {
  v <- v[is.finite(v) & v >= 0]
  if (length(v) == 0) return(c(mean = NA_real_, median = NA_real_, p90 = NA_real_, n = 0))
  c(
    mean = mean(v),
    median = stats::median(v),
    p90 = as.numeric(stats::quantile(v, probs = 0.9, na.rm = TRUE)),
    n = length(v)
  )
}

runtime_ok <- runtime_long[runtime_long$status == "ok", , drop = FALSE]
make_runtime_summary <- function(dat, type_label) {
  if (nrow(dat) == 0) {
    return(data.frame(
      scenario = character(0),
      method = character(0),
      mean = numeric(0),
      median = numeric(0),
      p90 = numeric(0),
      n = numeric(0),
      type = character(0),
      stringsAsFactors = FALSE
    ))
  }
  groups <- split(dat, list(dat$scenario, dat$method), drop = TRUE)
  rows <- lapply(groups, function(gd) {
    ss <- summarize_runtime(gd$runtime_sec)
    data.frame(
      scenario = unique(gd$scenario),
      method = unique(gd$method),
      mean = as.numeric(ss["mean"]),
      median = as.numeric(ss["median"]),
      p90 = as.numeric(ss["p90"]),
      n = as.numeric(ss["n"]),
      type = type_label,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

runtime_summary <- rbind(
  make_runtime_summary(runtime_long, "all_runs"),
  make_runtime_summary(runtime_ok, "ok_only")
)

make_runtime_overall <- function(dat, type_label) {
  if (nrow(dat) == 0) {
    return(data.frame(
      method = character(0),
      mean = numeric(0),
      median = numeric(0),
      p90 = numeric(0),
      n = numeric(0),
      type = character(0),
      stringsAsFactors = FALSE
    ))
  }
  groups <- split(dat, dat$method, drop = TRUE)
  rows <- lapply(groups, function(gd) {
    ss <- summarize_runtime(gd$runtime_sec)
    data.frame(
      method = unique(gd$method),
      mean = as.numeric(ss["mean"]),
      median = as.numeric(ss["median"]),
      p90 = as.numeric(ss["p90"]),
      n = as.numeric(ss["n"]),
      type = type_label,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

runtime_overall <- rbind(
  make_runtime_overall(runtime_long, "overall_all_runs"),
  make_runtime_overall(runtime_ok, "overall_ok_only")
)

out_runtime_summary <- file.path(out_dir, "benchmark_4methods_runtime_summary.tsv")
write.table(runtime_summary, out_runtime_summary, sep = "\t", row.names = FALSE, quote = FALSE)
out_runtime_overall <- file.path(out_dir, "benchmark_4methods_runtime_overall.tsv")
write.table(runtime_overall, out_runtime_overall, sep = "\t", row.names = FALSE, quote = FALSE)

method_cols <- c(
  p_rphylopars_gene_tree = "Rphylopars (gene tree)",
  p_epgls_species = "E-PGLS (ws=species)",
  p_epgls_gene_id = "E-PGLS (ws=gene_id)",
  p_hierarchical = "Hierarchical prototype"
)
pair_idx <- combn(names(method_cols), 2, simplify = FALSE)
plot_rows <- list()
pid <- 1L
for (pr in pair_idx) {
  tmp <- df[, c("scenario", "iter", pr), drop = FALSE]
  colnames(tmp) <- c("scenario", "iter", "p1", "p2")
  tmp <- tmp[is.finite(tmp$p1) & is.finite(tmp$p2) & tmp$p1 > 0 & tmp$p2 > 0, , drop = FALSE]
  if (nrow(tmp) == 0) next
  tmp$comparison <- paste0(method_cols[[pr[1]]], " vs ", method_cols[[pr[2]]])
  plot_rows[[pid]] <- tmp
  pid <- pid + 1L
}
if (length(plot_rows) == 0) {
  stop("No paired p-values were available for scatter plots.")
}
df_plot <- do.call(rbind, plot_rows)

metrics <- do.call(rbind, lapply(split(df_plot, df_plot$comparison), function(d) {
  rho <- suppressWarnings(cor(d$p1, d$p2, method = "spearman"))
  data.frame(
    comparison = unique(d$comparison),
    n = nrow(d),
    rho = rho,
    label_x = 0.03,
    label_y = 0.97,
    stringsAsFactors = FALSE
  )
}))
metrics$label <- sprintf("n = %d, rho = %.3f", metrics$n, metrics$rho)

g <- ggplot(df_plot, aes(x = p1, y = p2, color = scenario)) +
  geom_point(alpha = 0.55, size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.4, color = "gray35") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  facet_wrap(~comparison, ncol = 3, scales = "fixed") +
  labs(
    x = "p value (method 1)",
    y = "p value (method 2)",
    color = "Scenario",
    title = "Synthetic benchmark: p value comparisons across four methods"
  ) +
  theme_bw(base_size = 11)
g <- g + geom_text(
  data = metrics,
  aes(x = label_x, y = label_y, label = label),
  inherit.aes = FALSE,
  hjust = 0,
  vjust = 1,
  size = 3.1,
  color = "black"
)

out_pdf <- file.path(out_dir, "benchmark_4methods_pvalue_scatter.pdf")
out_png <- file.path(out_dir, "benchmark_4methods_pvalue_scatter.png")
ggsave(out_pdf, g, width = 14, height = 8.2)
ggsave(out_png, g, width = 14, height = 8.2, dpi = 200)

g_rt <- ggplot(runtime_long, aes(x = method, y = runtime_sec, fill = method)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  facet_wrap(~scenario, ncol = 2, scales = "free_y") +
  labs(
    x = "Method",
    y = "Runtime (sec/run)",
    title = "Synthetic benchmark: runtime comparison"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "none")

out_rt_pdf <- file.path(out_dir, "benchmark_4methods_runtime_boxplot.pdf")
out_rt_png <- file.path(out_dir, "benchmark_4methods_runtime_boxplot.png")
ggsave(out_rt_pdf, g_rt, width = 11, height = 7)
ggsave(out_rt_png, g_rt, width = 11, height = 7, dpi = 200)

cat("p-value table:", out_pval, "\n")
cat("summary:", out_summary, "\n")
cat("status counts:", out_status, "\n")
cat("runtime table:", out_runtime, "\n")
cat("runtime summary:", out_runtime_summary, "\n")
cat("runtime overall:", out_runtime_overall, "\n")
cat("scatter pdf:", out_pdf, "\n")
cat("scatter png:", out_png, "\n")
cat("runtime boxplot pdf:", out_rt_pdf, "\n")
cat("runtime boxplot png:", out_rt_png, "\n")
cat("n_iter:", n_iter, "\n")
cat("rrpp_iter:", rrpp_iter, "\n")
