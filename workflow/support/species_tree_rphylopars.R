#!/usr/bin/env Rscript

# Strict Rphylopars comparator for GeneGalleon expression-trait PGLS.
#
# Input preparation, paralog aggregation, and sampling-variance estimation are
# performed once by species_tree_pgls.py.  This script does not retry with a
# different model or silently discard uncertainty when a fit fails.

parse_args <- function(values) {
  out <- list()
  for (value in values) {
    if (!startsWith(value, "--") || !grepl("=", value, fixed = TRUE)) {
      stop(sprintf("Arguments must use --name=value syntax: %s", value), call. = FALSE)
    }
    pieces <- strsplit(sub("^--", "", value), "=", fixed = TRUE)[[1]]
    name <- pieces[[1]]
    setting <- paste(pieces[-1], collapse = "=")
    if (!nzchar(name) || name %in% names(out)) {
      stop(sprintf("Invalid or duplicate argument: %s", value), call. = FALSE)
    }
    out[[name]] <- setting
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c(
  "tree", "summary", "plan", "responses", "tree_id", "model", "parameter",
  "predictor_model", "predictor_parameter", "branch_length", "predictor_branch_length",
  "reml", "confidence_level", "sampling_covariance", "outfile", "status_out"
)
missing <- required[!required %in% names(args)]
if (length(missing) > 0) {
  stop(sprintf("Missing required arguments: %s", paste(missing, collapse = ", ")), call. = FALSE)
}

if (!requireNamespace("ape", quietly = TRUE)) {
  stop("ape is required for the Rphylopars comparator.", call. = FALSE)
}
if (!requireNamespace("Rphylopars", quietly = TRUE)) {
  stop("Rphylopars is required when species-rphylopars is selected.", call. = FALSE)
}

read_tsv <- function(path) {
  read.delim(
    path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("NA", "NaN", "nan", "?", "missing", "unknown", ".", "")
  )
}

summary_long <- read_tsv(args$summary)
plan <- read_tsv(args$plan)
tree <- ape::read.tree(args$tree)
branch_length_mode <- tolower(args$branch_length)
predictor_branch_length_mode <- tolower(args$predictor_branch_length)
if (!branch_length_mode %in% c("original", "unit") ||
    !predictor_branch_length_mode %in% c("original", "unit")) {
  stop("branch-length modes must be original or unit.", call. = FALSE)
}
if (identical(branch_length_mode, "unit")) {
  tree$edge.length <- rep(1, nrow(tree$edge))
}
responses <- Filter(nzchar, trimws(strsplit(args$responses, ",", fixed = TRUE)[[1]]))
confidence_level <- as.numeric(args$confidence_level)
if (!is.finite(confidence_level) || confidence_level <= 0 || confidence_level >= 1) {
  stop("confidence_level must lie strictly between zero and one.", call. = FALSE)
}
reml <- identical(tolower(args$reml), "yes")
package_version <- as.character(utils::packageVersion("Rphylopars"))

model_map <- c(
  brownian = "BM",
  ou = "OU",
  lambda = "lambda",
  kappa = "kappa",
  delta = "delta",
  eb = "EB",
  independent = "star"
)
model_key <- tolower(args$model)
model_supported <- model_key %in% names(model_map)
rphylopars_model <- if (model_supported) unname(model_map[[model_key]]) else NA_character_
parameter <- if (tolower(args$parameter) == "auto" || !nzchar(args$parameter)) NA_real_ else as.numeric(args$parameter)
if (!is.na(parameter) && !is.finite(parameter)) {
  stop("parameter must be auto or a finite number.", call. = FALSE)
}

shared_model_reason <- ""
if (!identical(tolower(args$model), tolower(args$predictor_model))) {
  shared_model_reason <- paste0(
    "Rphylopars uses one joint evolutionary model and cannot match distinct response/predictor models (",
    args$model, " versus ", args$predictor_model, ")"
  )
} else if (!identical(branch_length_mode, predictor_branch_length_mode)) {
  shared_model_reason <- "Rphylopars cannot match distinct response/predictor branch-length modes"
} else {
  response_parameter <- tolower(trimws(args$parameter))
  predictor_parameter <- tolower(trimws(args$predictor_parameter))
  both_auto <- response_parameter %in% c("", "auto") && predictor_parameter %in% c("", "auto")
  both_numeric_equal <- suppressWarnings({
    left <- as.numeric(response_parameter)
    right <- as.numeric(predictor_parameter)
    is.finite(left) && is.finite(right) && isTRUE(all.equal(left, right, tolerance = 1e-12))
  })
  if (!both_auto && !both_numeric_equal) {
    shared_model_reason <- "Rphylopars cannot match distinct response/predictor evolution parameters"
  }
}

split_setting <- function(value) {
  if (is.na(value) || !nzchar(value) || identical(value, ".")) character(0) else Filter(nzchar, trimws(strsplit(value, ",", fixed = TRUE)[[1]]))
}

empty_result <- data.frame(
  analysis_method = character(), aggregation = character(), estimand = character(),
  analysis_id = character(), model_id = character(), tree_id = character(),
  response = character(), response_type = character(), response_family = character(),
  link_function = character(), term = character(), source_term = character(),
  predictor_type = character(), coefficient = numeric(), standard_error = numeric(),
  statistic = numeric(), degrees_of_freedom = numeric(), p_value = numeric(),
  confidence_level = numeric(), confidence_interval_lower = numeric(),
  confidence_interval_upper = numeric(), n_species = integer(),
  evolution_model = character(), evolution_parameter_name = character(),
  evolution_parameter = numeric(), evolution_parameter_status = character(),
  evolutionary_rate = numeric(), r_squared = numeric(), adjusted_r_squared = numeric(),
  log_likelihood = numeric(), aic = numeric(), bic = numeric(),
  inference_method = character(), inference_status = character(),
  optimizer_converged = character(), optimizer_message = character(),
  sampling_covariance_mode = character(), engine_version = character(),
  stringsAsFactors = FALSE
)

empty_status <- data.frame(
  tree_id = character(), analysis_method = character(), aggregation = character(),
  analysis_id = character(), response = character(), status = character(),
  reason = character(), n_species = integer(), n_result_rows = integer(),
  engine_version = character(), stringsAsFactors = FALSE
)

result_rows <- list()
status_rows <- list()
result_index <- 1L
status_index <- 1L

add_status <- function(aggregation, analysis_id, response, status, reason, n_species = NA_integer_, n_result_rows = 0L) {
  status_rows[[status_index]] <<- data.frame(
    tree_id = args$tree_id,
    analysis_method = "species_rphylopars",
    aggregation = aggregation,
    analysis_id = analysis_id,
    response = response,
    status = status,
    reason = reason,
    n_species = n_species,
    n_result_rows = n_result_rows,
    engine_version = package_version,
    stringsAsFactors = FALSE
  )
  status_index <<- status_index + 1L
}

trait_slice <- function(aggregation, source, trait) {
  summary_long[
    summary_long$aggregation == aggregation &
      summary_long$source == source &
      summary_long$trait == trait,
    c("leaf_name", "value", "sampling_variance", "has_offdiagonal_sampling_covariance"),
    drop = FALSE
  ]
}

for (aggregation in unique(summary_long$aggregation)) {
  for (plan_index in seq_len(nrow(plan))) {
    analysis_id <- as.character(plan$analysis_id[[plan_index]])
    predictors <- split_setting(as.character(plan$predictors[[plan_index]]))
    categorical <- split_setting(as.character(plan$categorical_predictors[[plan_index]]))
    ordered <- split_setting(as.character(plan$ordered_predictors[[plan_index]]))
    for (response in responses) {
      unsupported_reason <- ""
      if (!model_supported) {
        unsupported_reason <- sprintf("Rphylopars does not support evolution model '%s'", args$model)
      } else if (nzchar(shared_model_reason)) {
        unsupported_reason <- shared_model_reason
      } else if (length(categorical) > 0 || length(ordered) > 0) {
        unsupported_reason <- "Rphylopars comparison is restricted to continuous predictors"
      }
      if (nzchar(unsupported_reason)) {
        add_status(aggregation, analysis_id, response, "not_estimable", unsupported_reason)
        next
      }

      selected <- list(trait_slice(aggregation, "response", response))
      names(selected) <- response
      for (predictor in predictors) {
        selected[[predictor]] <- trait_slice(
          aggregation,
          paste0("predictor:", analysis_id),
          predictor
        )
      }
      if (any(vapply(selected, nrow, integer(1)) == 0L)) {
        add_status(aggregation, analysis_id, response, "not_estimable", "prepared species summary lacks a selected trait")
        next
      }
      offdiagonal <- any(vapply(selected, function(frame) {
        any(as.character(frame$has_offdiagonal_sampling_covariance) == "yes")
      }, logical(1)))
      if (offdiagonal && identical(args$sampling_covariance, "require-diagonal")) {
        add_status(
          aggregation,
          analysis_id,
          response,
          "not_estimable",
          "Rphylopars phenocov_list cannot represent cross-species sampling covariance"
        )
        next
      }

      species <- tree$tip.label
      trait_data <- data.frame(species = factor(species, levels = species), check.names = FALSE)
      variances <- matrix(0, nrow = length(species), ncol = length(selected))
      colnames(variances) <- names(selected)
      rownames(variances) <- species
      valid <- TRUE
      for (trait in names(selected)) {
        frame <- selected[[trait]]
        indices <- match(species, as.character(frame$leaf_name))
        if (anyNA(indices)) {
          valid <- FALSE
          break
        }
        values <- suppressWarnings(as.numeric(frame$value[indices]))
        sampling_variance <- suppressWarnings(as.numeric(frame$sampling_variance[indices]))
        if (any(!is.finite(values)) || any(!is.finite(sampling_variance)) || any(sampling_variance < 0)) {
          valid <- FALSE
          break
        }
        trait_data[[trait]] <- values
        variances[, trait] <- sampling_variance
      }
      if (!valid) {
        add_status(aggregation, analysis_id, response, "not_estimable", "species summary contains missing or invalid values")
        next
      }

      phenocov_list <- list()
      if (any(variances > 0)) {
        phenocov_list <- lapply(seq_along(species), function(index) {
          matrix_value <- diag(variances[index, ], nrow = ncol(variances), ncol = ncol(variances))
          dimnames(matrix_value) <- list(colnames(variances), colnames(variances))
          matrix_value
        })
        names(phenocov_list) <- species
      }
      formula <- stats::reformulate(predictors, response = response)
      fit_args <- list(
        formula = formula,
        trait_data = trait_data,
        tree = tree,
        model = rphylopars_model,
        pheno_error = FALSE,
        phylo_correlated = TRUE,
        pheno_correlated = FALSE,
        REML = reml,
        phenocov_list = phenocov_list
      )
      if (!is.na(parameter)) {
        fit_args$model_par_fixed <- parameter
      }
      fit <- tryCatch(
        suppressWarnings(do.call(Rphylopars::phylopars.lm, fit_args)),
        error = function(error) error
      )
      if (inherits(fit, "error")) {
        add_status(
          aggregation,
          analysis_id,
          response,
          "not_estimable",
          paste("Rphylopars fit failed:", conditionMessage(fit)),
          n_species = length(species)
        )
        next
      }

      coefficients <- as.numeric(fit$coefficients)
      coefficient_names <- names(fit$coefficients)
      standard_errors <- rep(NA_real_, length(coefficients))
      statistics <- rep(NA_real_, length(coefficients))
      p_values <- rep(NA_real_, length(coefficients))
      if (length(fit$SEs) == length(coefficients)) standard_errors <- as.numeric(fit$SEs)
      if (length(fit$ts) == length(coefficients)) statistics <- as.numeric(fit$ts)
      if (length(fit$ps) == length(coefficients)) p_values <- as.numeric(fit$ps)
      degrees_of_freedom <- if (!is.null(fit$df2)) as.numeric(fit$df2) else length(species) - length(coefficients)
      critical <- stats::qt((1 + confidence_level) / 2, df = degrees_of_freedom)
      lower <- coefficients - critical * standard_errors
      upper <- coefficients + critical * standard_errors
      parameter_name <- switch(
        model_key,
        lambda = "lambda", ou = "alpha", kappa = "kappa", delta = "delta",
        eb = "rate_change", ""
      )
      parameter_field <- if (identical(model_key, "eb")) "rate" else parameter_name
      estimated_parameter <- NA_real_
      if (nzchar(parameter_field) && !is.null(fit$PPE$model[[parameter_field]])) {
        estimated_parameter <- suppressWarnings(as.numeric(fit$PPE$model[[parameter_field]])[1])
      } else if (!is.na(parameter)) {
        estimated_parameter <- parameter
      }
      parameter_status <- if (!nzchar(parameter_name)) {
        "not_applicable"
      } else if (!is.na(parameter)) {
        "fixed"
      } else if (is.finite(estimated_parameter)) {
        "estimated"
      } else {
        "not_reported"
      }
      fit_rows <- lapply(seq_along(coefficients), function(index) {
        term <- coefficient_names[[index]]
        data.frame(
          analysis_method = "species_rphylopars",
          aggregation = aggregation,
          estimand = sprintf("species-level %s paralog expression", aggregation),
          analysis_id = analysis_id,
          model_id = sprintf("species_rphylopars:%s:%s:%s", aggregation, analysis_id, response),
          tree_id = args$tree_id,
          response = response,
          response_type = "continuous",
          response_family = "gaussian",
          link_function = "identity",
          term = term,
          source_term = if (term %in% predictors) term else term,
          predictor_type = if (term == "(Intercept)") "intercept" else "continuous",
          coefficient = coefficients[[index]],
          standard_error = standard_errors[[index]],
          statistic = statistics[[index]],
          degrees_of_freedom = degrees_of_freedom,
          p_value = p_values[[index]],
          confidence_level = confidence_level,
          confidence_interval_lower = lower[[index]],
          confidence_interval_upper = upper[[index]],
          n_species = length(species),
          evolution_model = model_key,
          evolution_parameter_name = parameter_name,
          evolution_parameter = estimated_parameter,
          evolution_parameter_status = parameter_status,
          evolutionary_rate = if (!is.null(fit$sigma)) as.numeric(fit$sigma)^2 else NA_real_,
          r_squared = if (!is.null(fit$R2)) as.numeric(fit$R2) else NA_real_,
          adjusted_r_squared = if (!is.null(fit$R2adj)) as.numeric(fit$R2adj) else NA_real_,
          log_likelihood = if (!is.null(fit$logLik)) as.numeric(fit$logLik) else NA_real_,
          aic = suppressWarnings(as.numeric(stats::AIC(fit))),
          bic = suppressWarnings(as.numeric(stats::BIC(fit))),
          inference_method = "wald",
          inference_status = if (is.finite(p_values[[index]])) "ok" else "not_reported",
          optimizer_converged = "not_reported",
          optimizer_message = "Rphylopars does not expose optimizer convergence in phylopars.lm output",
          sampling_covariance_mode = if (offdiagonal) "diagonalized" else "exact-diagonal",
          engine_version = package_version,
          stringsAsFactors = FALSE
        )
      })
      result_rows[[result_index]] <- do.call(rbind, fit_rows)
      result_index <- result_index + 1L
      finite_inference <- any(is.finite(p_values))
      add_status(
        aggregation,
        analysis_id,
        response,
        if (finite_inference) "ok" else "not_estimable",
        if (!finite_inference) {
          "Rphylopars returned no finite coefficient p-values"
        } else if (offdiagonal) {
          "cross-species sampling covariance was diagonalized"
        } else {
          ""
        },
        n_species = length(species),
        n_result_rows = length(coefficients)
      )
    }
  }
}

results <- if (length(result_rows) > 0) do.call(rbind, result_rows) else empty_result
statuses <- if (length(status_rows) > 0) do.call(rbind, status_rows) else empty_status
write.table(results, args$outfile, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
write.table(statuses, args$status_out, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
