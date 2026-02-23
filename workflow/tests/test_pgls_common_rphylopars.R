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
  if (dir.exists(file.path(repo_root_candidate, "workflow", "support"))) {
    repo_root <- normalizePath(repo_root_candidate, winslash = "/", mustWork = TRUE)
  }
}

if (is.na(repo_root)) {
  cwd_candidate <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (dir.exists(file.path(cwd_candidate, "workflow", "support"))) {
    repo_root <- cwd_candidate
  }
}

if (is.na(repo_root)) {
  stop(sprintf(
    "Could not resolve repository root from script path (%s) or current working directory (%s).",
    script_path,
    getwd()
  ))
}

source(file.path(repo_root, "workflow", "support", "pgls_common.R"))

# build_phenocov_input: no within-species variation
df_no_var <- data.frame(
  species = c("sp1", "sp2"),
  expr_1 = c(1, 2),
  expr_2 = c(1, 2),
  trait = c(0.1, 0.2),
  stringsAsFactors = FALSE
)
res_no_var <- build_phenocov_input(df_no_var, "trait", c("expr_1", "expr_2"), "expr")
stopifnot(isFALSE(res_no_var$has_within_species_variation))

# build_phenocov_input: within-species variation exists
df_with_var <- df_no_var
df_with_var$expr_2 <- c(1.6, 2.4)
res_with_var <- build_phenocov_input(df_with_var, "trait", c("expr_1", "expr_2"), "expr")
stopifnot(isTRUE(res_with_var$has_within_species_variation))

# sanitize_phenocov_list
raw_phenocov <- list(
  sp1 = matrix(c(NA, Inf, -Inf, -1), nrow = 2, byrow = TRUE),
  sp2 = matrix(c(0, 3, 5, 0), nrow = 2, byrow = TRUE)
)
names(raw_phenocov) <- c("sp1", "sp2")
sanitized <- sanitize_phenocov_list(raw_phenocov, "trait", "expr", ridge = 1e-6)
for (nm in names(sanitized)) {
  m <- sanitized[[nm]]
  stopifnot(all(is.finite(m)))
  stopifnot(isTRUE(all.equal(m, t(m))))
  stopifnot(all(diag(m) >= 1e-6))
}

# sanitize_phenocov_list_generic + sanitize_cov_matrix
generic_raw <- list(
  sp1 = matrix(c(NA, Inf, -Inf, -2), nrow = 2, byrow = TRUE),
  sp2 = matrix(c(0, 7, 3, 0), nrow = 2, byrow = TRUE)
)
generic_sanitized <- sanitize_phenocov_list_generic(generic_raw, ridge = 1e-6)
for (nm in names(generic_sanitized)) {
  m <- generic_sanitized[[nm]]
  stopifnot(all(is.finite(m)))
  stopifnot(isTRUE(all.equal(m, t(m))))
  stopifnot(all(diag(m) >= 1e-6))
}

# extract_try_error_message
mock_err <- try(stop("extract error"), silent = TRUE)
stopifnot(grepl("extract error", extract_try_error_message(mock_err), fixed = TRUE))

# fit_phylopars_lm_with_retries: succeeds after retries
.mock_counter <- 0L
phylopars.lm <- function(...) {
  .mock_counter <<- .mock_counter + 1L
  if (.mock_counter < 3L) {
    stop("mock failure")
  }
  structure(
    list(R2 = 0.5, R2adj = 0.4, sigma = 1.0, Fstat = 3.0, pval = 0.1, logLik = -10),
    class = "mock_phylopars"
  )
}
AIC.mock_phylopars <- function(object, ...) 123
BIC.mock_phylopars <- function(object, ...) 234

dummy_trait <- data.frame(
  species = c("sp1", "sp2"),
  trait = c(1.0, 1.5),
  expr = c(2.0, 2.5),
  stringsAsFactors = FALSE
)
dummy_tree <- structure(
  list(
    tip.label = c("sp1", "sp2"),
    edge = matrix(c(3, 1, 3, 2), nrow = 2, byrow = TRUE),
    edge.length = c(1, 1),
    Nnode = 1
  ),
  class = "phylo"
)
fit_ok <- fit_phylopars_lm_with_retries(
  formula_obj = trait ~ expr,
  trait_data = dummy_trait,
  tree = dummy_tree,
  model = "BM",
  pheno_error = TRUE,
  phylo_correlated = TRUE,
  pheno_correlated = TRUE,
  phenocov_list = list(),
  fit_mode_label = "wide",
  fit_fun = phylopars.lm
)
stopifnot(!is.null(fit_ok$fit))
stopifnot(identical(fit_ok$fit_mode, "wide_npd_noz"))

# fit_phylopars_lm_with_retries: all attempts fail
phylopars.lm <- function(...) {
  stop("always fail")
}
fit_fail <- fit_phylopars_lm_with_retries(
  formula_obj = trait ~ expr,
  trait_data = dummy_trait,
  tree = dummy_tree,
  model = "BM",
  pheno_error = TRUE,
  phylo_correlated = TRUE,
  pheno_correlated = TRUE,
  phenocov_list = list(),
  fit_mode_label = "wide",
  fit_fun = phylopars.lm
)
stopifnot(is.null(fit_fail$fit))
stopifnot(grepl("always fail", fit_fail$error_message, fixed = TRUE))

cat("test_pgls_common_rphylopars.R: OK\n")
