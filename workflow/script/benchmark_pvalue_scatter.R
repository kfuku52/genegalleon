#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[[1]] else file.path("workspace", "output", "benchmark_epgls")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ape)
  library(Rphylopars)
  library(ggplot2)
})

# Use separate RRPP install prepared for E-PGLS verification.
.libPaths(c("/tmp/epgls_lib", .libPaths()))
suppressPackageStartupMessages(library(RRPP))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_path <- sub("^--file=", "", script_arg[1])
  script_path <- gsub("~\\+~", " ", script_path)
  script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
} else {
  script_dir <- getwd()
}
source(file.path(script_dir, "pgls_common.R"))

simulate_one <- function(n_species = 12,
                         min_rep = 2,
                         max_rep = 5,
                         miss_prob = 0.2,
                         sigma_phylo_x = 1.0,
                         sigma_phylo_y = 1.0,
                         sigma_within_x = 0.6,
                         beta = 0.7) {
  tree <- ape::rtree(n_species)
  tree$tip.label <- sprintf("sp%02d", seq_len(n_species))

  x_sp <- as.numeric(ape::rTraitCont(tree, model = "BM", sigma = sigma_phylo_x))
  names(x_sp) <- tree$tip.label
  y_phy <- as.numeric(ape::rTraitCont(tree, model = "BM", sigma = sigma_phylo_y))
  names(y_phy) <- tree$tip.label
  y_sp <- beta * x_sp + y_phy

  reps <- sample(min_rep:max_rep, n_species, replace = TRUE)
  max_r <- max_rep

  wide <- data.frame(species = tree$tip.label, y = as.numeric(y_sp[tree$tip.label]), stringsAsFactors = FALSE)
  long <- data.frame(species = character(0), y = numeric(0), x = numeric(0), stringsAsFactors = FALSE)
  for (r in seq_len(max_r)) wide[[sprintf("x_%d", r)]] <- NA_real_

  for (i in seq_len(n_species)) {
    sp <- tree$tip.label[i]
    k <- reps[i]
    x_rep <- x_sp[sp] + stats::rnorm(k, mean = 0, sd = sigma_within_x)
    keep <- stats::runif(k) > miss_prob
    x_rep[!keep] <- NA_real_
    if (all(is.na(x_rep))) x_rep[1] <- x_sp[sp] + stats::rnorm(1, 0, sigma_within_x)
    wide[i, sprintf("x_%d", seq_len(k))] <- x_rep
    long_i <- data.frame(species = rep(sp, k), y = rep(y_sp[sp], k), x = as.numeric(x_rep), stringsAsFactors = FALSE)
    long <- rbind(long, long_i)
  }

  x_cols <- grep("^x_[0-9]+$", colnames(wide), value = TRUE)
  wide$mean_x <- apply(wide[, x_cols, drop = FALSE], 1, function(v) {
    vv <- suppressWarnings(as.numeric(v))
    if (all(is.na(vv))) return(NA_real_)
    mean(vv, na.rm = TRUE)
  })

  long <- long[!is.na(long$x) & !is.na(long$y) & !is.na(long$species) & nzchar(long$species), , drop = FALSE]
  rownames(long) <- NULL

  list(tree = tree, wide = wide, long = long)
}

drop_tree_to_species <- function(tree, species_vec) {
  keep <- unique(as.character(species_vec))
  t_use <- tree
  drop_tips <- t_use$tip.label[!(t_use$tip.label %in% keep)]
  if (length(drop_tips) > 0) t_use <- ape::drop.tip(t_use, drop_tips)
  t_use
}

fit_rphylopars_mean <- function(wide, tree) {
  out <- list(status = "error", pval = NA_real_)
  if (!("mean_x" %in% colnames(wide))) {
    out$status <- "no_data"
    return(out)
  }
  td <- wide[, c("species", "y", "mean_x"), drop = FALSE]
  colnames(td)[3] <- "x"
  td$y <- suppressWarnings(as.numeric(td$y))
  td$x <- suppressWarnings(as.numeric(td$x))
  td <- td[!is.na(td$species) & !is.na(td$y) & !is.na(td$x), , drop = FALSE]
  if (nrow(td) < 3 || length(unique(td$species)) < 3) {
    out$status <- "no_data"
    return(out)
  }
  t_use <- drop_tree_to_species(tree, td$species)
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

fit_rphylopars_phenocov <- function(wide, tree) {
  out <- list(status = "error", pval = NA_real_)
  x_cols <- grep("^x(_[0-9]+)?$", colnames(wide), value = TRUE)
  if (length(x_cols) == 0) {
    out$status <- "no_data"
    return(out)
  }
  ph <- build_phenocov_input(wide, trait_col = "y", expression_cols = x_cols, expression_base = "x")
  td <- ph[["trait_data"]]
  pcl <- ph[["phenocov_list"]]
  if (is.null(td) || nrow(td) < 3) {
    out$status <- "no_data"
    return(out)
  }
  t_use <- drop_tree_to_species(tree, td$species)
  td <- td[td$species %in% t_use$tip.label, , drop = FALSE]
  fit <- try(do.call(Rphylopars::phylopars.lm, list(
    formula = y ~ x,
    trait_data = td,
    tree = t_use,
    model = "BM",
    pheno_error = FALSE,
    phylo_correlated = TRUE,
    pheno_correlated = FALSE,
    phenocov_list = pcl
  )), silent = TRUE)
  if (inherits(fit, "try-error")) {
    out$status <- "error"
    return(out)
  }
  out$status <- "ok"
  out$pval <- suppressWarnings(as.numeric(fit$pval))
  out
}

fit_epgls <- function(long, tree, iter = 19) {
  out <- list(status = "error", pval = NA_real_)
  if (nrow(long) < 4) {
    out$status <- "no_data"
    return(out)
  }
  t_use <- drop_tree_to_species(tree, long$species)
  dat <- long[long$species %in% t_use$tip.label, , drop = FALSE]
  dat <- dat[!is.na(dat$x) & !is.na(dat$y), , drop = FALSE]
  if (nrow(dat) < 4 || length(unique(dat$species)) < 3) {
    out$status <- "no_data"
    return(out)
  }
  fit <- try(RRPP::lm.rrpp.ws(
    y ~ x,
    data = dat,
    subjects = "species",
    Cov = ape::vcv.phylo(t_use),
    iter = iter,
    print.progress = FALSE,
    Parallel = FALSE
  ), silent = TRUE)
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

scenarios <- data.frame(
  scenario = c("A_balanced", "B_missing_mid", "C_missing_high", "D_few_species"),
  n_species = c(12, 12, 12, 8),
  min_rep = c(3, 2, 2, 2),
  max_rep = c(5, 5, 4, 3),
  miss_prob = c(0.10, 0.35, 0.60, 0.45),
  sigma_within_x = c(0.4, 0.6, 0.8, 0.7),
  stringsAsFactors = FALSE
)

set.seed(42)
n_iter <- 250
rows <- list()
idx <- 1L

for (s in seq_len(nrow(scenarios))) {
  sc <- scenarios[s, ]
  for (i in seq_len(n_iter)) {
    d <- simulate_one(
      n_species = sc$n_species,
      min_rep = sc$min_rep,
      max_rep = sc$max_rep,
      miss_prob = sc$miss_prob,
      sigma_within_x = sc$sigma_within_x
    )
    m_mean <- fit_rphylopars_mean(d$wide, d$tree)
    m_pheno <- fit_rphylopars_phenocov(d$wide, d$tree)
    m_ep <- fit_epgls(d$long, d$tree)
    rows[[idx]] <- data.frame(
      scenario = sc$scenario,
      iter = i,
      p_epgls = m_ep$pval,
      status_epgls = m_ep$status,
      p_rphylopars_mean = m_mean$pval,
      status_rphylopars_mean = m_mean$status,
      p_rphylopars_phenocov = m_pheno$pval,
      status_rphylopars_phenocov = m_pheno$status,
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }
}

df <- do.call(rbind, rows)
out_tsv <- file.path(out_dir, "benchmark_rphylopars_epgls_pvalues.tsv")
write.table(df, out_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

make_scatter_df <- function(dat, y_col, comparison_name) {
  tmp <- dat[, c("scenario", "iter", "p_epgls", y_col), drop = FALSE]
  colnames(tmp)[4] <- "p_other"
  tmp <- tmp[is.finite(tmp$p_epgls) & is.finite(tmp$p_other), , drop = FALSE]
  tmp <- tmp[tmp$p_epgls > 0 & tmp$p_other > 0, , drop = FALSE]
  tmp$comparison <- comparison_name
  tmp
}

df_mean <- make_scatter_df(df, "p_rphylopars_mean", "E-PGLS vs Rphylopars(mean)")
df_pheno <- make_scatter_df(df, "p_rphylopars_phenocov", "E-PGLS vs Rphylopars(phenocov)")
df_plot <- rbind(df_mean, df_pheno)

if (nrow(df_plot) == 0) {
  stop("No paired p-values are available for plotting.")
}

metrics <- do.call(rbind, lapply(split(df_plot, df_plot$comparison), function(d) {
  rho <- suppressWarnings(cor(d$p_epgls, d$p_other, method = "spearman"))
  data.frame(
    comparison = unique(d$comparison),
    n = nrow(d),
    rho = rho,
    label_x = min(d$p_epgls, na.rm = TRUE),
    label_y = max(d$p_other, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
metrics$label <- sprintf("n = %d, rho = %.3f", metrics$n, metrics$rho)

g <- ggplot(df_plot, aes(x = p_epgls, y = p_other, color = scenario)) +
  geom_point(alpha = 0.6, size = 1.4) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.4, color = "gray35") +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~comparison, ncol = 2, scales = "free") +
  labs(
    x = "E-PGLS p value",
    y = "Rphylopars p value",
    color = "Scenario",
    title = "Synthetic benchmark: p value comparison"
  ) +
  theme_bw(base_size = 11)
g <- g + geom_text(
  data = metrics,
  aes(x = label_x, y = label_y, label = label),
  inherit.aes = FALSE,
  hjust = 0,
  vjust = 1,
  size = 3.2,
  color = "black"
)

out_pdf <- file.path(out_dir, "benchmark_rphylopars_epgls_pvalue_scatter.pdf")
out_png <- file.path(out_dir, "benchmark_rphylopars_epgls_pvalue_scatter.png")
ggsave(out_pdf, g, width = 10, height = 5.6)
ggsave(out_png, g, width = 10, height = 5.6, dpi = 200)

cat("p-value table:", out_tsv, "\n")
cat("scatter pdf:", out_pdf, "\n")
cat("scatter png:", out_png, "\n")
cat("paired points (mean):", nrow(df_mean), "\n")
cat("paired points (phenocov):", nrow(df_pheno), "\n")
