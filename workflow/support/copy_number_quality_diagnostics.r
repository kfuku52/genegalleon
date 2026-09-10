#!/usr/bin/env Rscript

# BUSCO is a quality diagnostic, never an automatic false-positive classifier.
.quality_source <- if (sys.nframe() > 0 && !is.null(sys.frame(1)$ofile)) sys.frame(1)$ofile else
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
.quality_dir <- dirname(normalizePath(.quality_source, mustWork = TRUE))
.quality_pgls <- new.env(parent = globalenv())
.quality_old_guard <- Sys.getenv("GG_ORTHOGROUP_COPY_NUMBER_TRAIT_PGLS_NO_MAIN", unset = NA_character_)
Sys.setenv(GG_ORTHOGROUP_COPY_NUMBER_TRAIT_PGLS_NO_MAIN = "1")
sys.source(file.path(.quality_dir, "orthogroup_copy_number_trait_pgls.r"), envir = .quality_pgls)
.quality_pgls$support_dir <- .quality_dir
.quality_pgls$script_dir <- .quality_dir
if (is.na(.quality_old_guard)) Sys.unsetenv("GG_ORTHOGROUP_COPY_NUMBER_TRAIT_PGLS_NO_MAIN") else
  Sys.setenv(GG_ORTHOGROUP_COPY_NUMBER_TRAIT_PGLS_NO_MAIN = .quality_old_guard)

quality_read_table <- function(path) {
  df <- .quality_pgls$read_tsv_base(path, na = c("", "NA", "NaN", "nan"), col_classes = "character")
  if (!ncol(df) || any(!nzchar(names(df))) || anyDuplicated(names(df))) stop("Invalid table headers: ", path)
  names(df)[[1]] <- "species"
  if (anyDuplicated(names(df))) stop("Duplicate species header: ", path)
  df$species <- .quality_pgls$normalize_species_label(df$species)
  if (anyNA(df$species) || any(!nzchar(df$species))) stop("Missing species labels: ", path)
  .quality_pgls$assert_unique_labels(df$species, path)
  df
}

quality_numeric <- function(values, label, strict = TRUE) {
  out <- suppressWarnings(as.numeric(values))
  bad <- !is.na(values) & !is.finite(out)
  if (any(bad)) {
    if (strict) stop("Non-numeric or infinite values in ", label)
    return(NULL)
  }
  out
}

quality_load_busco <- function(file_busco, busco_short_dir, species, work) {
  if (nzchar(file_busco) && nzchar(busco_short_dir)) stop("Choose a BUSCO table or short-summary directory, not both.")
  if (nzchar(busco_short_dir)) {
    file_busco <- file.path(work, "busco_input.tsv")
    status <- system2("python", shQuote(c(file.path(.quality_dir, "busco_quality_metadata.py"),
                                        "--directory", busco_short_dir, "--output", file_busco)))
    if (status != 0L) stop("Could not extract BUSCO quality metadata.")
  }
  if (nzchar(file_busco)) {
    busco <- quality_read_table(file_busco)
    if (!"busco_complete_pct" %in% names(busco)) stop("BUSCO table requires busco_complete_pct (0–100).")
    busco$busco_complete_pct <- quality_numeric(busco$busco_complete_pct, "busco_complete_pct")
    if (any(busco$busco_complete_pct < 0 | busco$busco_complete_pct > 100, na.rm = TRUE)) {
      stop("busco_complete_pct must be in percent units, between 0 and 100.")
    }
  } else busco <- data.frame(species = character(), busco_complete_pct = numeric())
  busco <- busco[match(species, busco$species), , drop = FALSE]
  busco$species <- species
  for (column in c("lineage", "mode", "busco_version", "source")) {
    if (!column %in% names(busco)) busco[[column]] <- NA_character_
  }
  # Known incompatible protocols must not be pooled into an apparent quality gradient.
  for (column in c("lineage", "mode", "busco_version")) {
    observed <- busco[[column]][is.finite(busco$busco_complete_pct)]
    observed <- unique(observed[!is.na(observed) & nzchar(observed)])
    if (length(observed) > 1L) stop("Mixed BUSCO ", column, ": use a comparable quality table/subset.")
  }
  busco$quality_status <- ifelse(is.finite(busco$busco_complete_pct), "available", "missing")
  metadata_known <- !is.na(busco$lineage) & nzchar(busco$lineage) & !is.na(busco$mode) & nzchar(busco$mode)
  busco$comparability_status <- ifelse(metadata_known, "recorded_protocol", "metadata_incomplete")
  busco
}

quality_numeric_traits <- function(trait, busco, trait_arg = "all") {
  trait <- trait[match(busco$species, trait$species), , drop = FALSE]
  trait$species <- busco$species
  # The diagnostic is included even when the biological trait list is explicit.
  cols <- .quality_pgls$resolve_trait_cols(trait, trait_arg)
  if ("busco_complete_pct" %in% names(trait)) {
    existing <- quality_numeric(trait$busco_complete_pct, "trait busco_complete_pct")
    supplied <- busco$busco_complete_pct
    if (any(is.finite(existing) & (!is.finite(supplied) | abs(existing - supplied) > 1e-8))) {
      stop("Trait busco_complete_pct conflicts with the designated BUSCO quality source.")
    }
  }
  out <- data.frame(species = busco$species, busco_complete_pct = busco$busco_complete_pct)
  excluded <- data.frame(trait = character(), reason = character())
  for (column in setdiff(cols, "busco_complete_pct")) {
    values <- quality_numeric(trait[[column]], column, strict = FALSE)
    if (is.null(values)) excluded <- rbind(excluded, data.frame(trait = column, reason = "non_numeric")) else
      out[[column]] <- values
  }
  list(values = out, excluded = excluded)
}

quality_correlations <- function(traits, method = "spearman") {
  if (!method %in% c("pearson", "spearman")) stop("correlation_method must be pearson or spearman.")
  cols <- setdiff(names(traits), "species")
  pairs <- expand.grid(trait_x = cols, trait_y = cols, stringsAsFactors = FALSE)
  pairs$n_species <- 0L
  pairs$correlation <- NA_real_
  pairs$status <- "too_few_species"
  pairs$method <- method
  for (i in seq_len(nrow(pairs))) {
    x <- traits[[pairs$trait_x[[i]]]]
    y <- traits[[pairs$trait_y[[i]]]]
    usable <- is.finite(x) & is.finite(y)
    pairs$n_species[[i]] <- sum(usable)
    if (sum(usable) < 3L) next
    if (length(unique(x[usable])) < 2L || length(unique(y[usable])) < 2L) {
      pairs$status[[i]] <- "invariant_trait"
      next
    }
    pairs$correlation[[i]] <- cor(x[usable], y[usable], method = method)
    pairs$status[[i]] <- "ok"
  }
  pairs
}

quality_plot_correlations <- function(pairs, outdir, observation_caption = "") {
  cols <- unique(pairs$trait_x)
  pairs$trait_x <- factor(pairs$trait_x, levels = cols)
  pairs$trait_y <- factor(pairs$trait_y, levels = rev(cols))
  pairs$label <- ifelse(is.finite(pairs$correlation),
                        sprintf("%.2f\nn=%d", pairs$correlation, pairs$n_species),
                        sprintf("NA\nn=%d", pairs$n_species))
  pairs$text_color <- ifelse(is.finite(pairs$correlation) & abs(pairs$correlation) >= 0.65, "white", "#172033")
  labels <- function(x) ifelse(x == "busco_complete_pct", "BUSCO completeness (%)", x)
  p <- ggplot(pairs, aes(trait_x, trait_y, fill = correlation)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = label, color = text_color), size = 2.6) + scale_color_identity() +
    scale_x_discrete(labels = labels) + scale_y_discrete(labels = labels) +
    scale_fill_gradient2(low = "#2563A6", mid = "#FFFFFF", high = "#B63C45", midpoint = 0,
                         limits = c(-1, 1), na.value = "#E5E7EB", name = "Correlation") +
    coord_fixed() + labs(x = NULL, y = NULL, title = "Traits and BUSCO completeness",
                        subtitle = paste0(pairs$method[[1]], " correlation | pairwise complete species"),
                        caption = paste0("Descriptive correlations; not phylogenetically corrected.\nNA: insufficient data or no variation.", observation_caption)) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 55, hjust = 1),
          legend.position = "bottom", plot.title = element_text(face = "bold"),
          plot.title.position = "plot", plot.caption.position = "plot",
          plot.caption = element_text(hjust = 0, size = 8), plot.margin = margin(10, 12, 10, 12))
  size <- max(6, 3 + 0.6 * length(cols))
  # Embed PDF fonts so exported figures render consistently outside the container.
  ggsave(file.path(outdir, "trait_correlations.pdf"), p, device = grDevices::cairo_pdf,
         width = size + 1, height = size, limitsize = FALSE)
  ggsave(file.path(outdir, "trait_correlations.svg"), p, width = size + 1, height = size, limitsize = FALSE)
}

quality_fit <- function(data, tree, family, trait, min_species, adjusted = FALSE, response_family = "gaussian") {
  # The shared adapter validates raw integer counts and applies log1p exactly once.
  out <- .quality_pgls$fit_one_orthogroup_copy_number_trait(data, tree, family, trait,
    min_species = min_species, response_family = response_family,
    covariates = if (adjusted) "busco_quality" else character())
  out$response_family <- response_family
  out$link_function <- switch(response_family, gaussian = "identity", binomial = "logit", "log")
  out$covariance_estimator <- if (response_family == "gaussian") "gaussian-REML" else "laplace-ML"
  out
}

quality_adjust <- function(df, group = NULL) {
  df$quality_qvalue <- rep(NA_real_, nrow(df))
  groups <- if (is.null(group)) rep("all", nrow(df)) else df[[group]]
  for (key in unique(groups)) {
    planned <- which(groups == key)
    usable <- planned[df$status[planned] == "ok" & is.finite(df$pval[planned])]
    if (length(usable)) df$quality_qvalue[usable] <- p.adjust(df$pval[usable], "BH", n = length(planned))
  }
  df
}

quality_family_flags <- function(copy_matrix, busco, tree, min_species = 4L, alpha = 0.05) {
  rows <- lapply(colnames(copy_matrix), function(family) {
    data <- data.frame(species = rownames(copy_matrix), copy_number = copy_matrix[, family],
                       trait_value = busco$busco_complete_pct[match(rownames(copy_matrix), busco$species)])
    quality_fit(data, tree, family, "busco_complete_pct", min_species)
  })
  out <- if (length(rows)) do.call(rbind, rows) else .quality_pgls$empty_orthogroup_copy_number_trait_result()
  out <- quality_adjust(out)
  out$quality_flag <- ifelse(!is.finite(out$quality_qvalue), "not_assessable",
                            ifelse(out$quality_qvalue < alpha, "quality_associated", "no_quality_association_detected"))
  out$predictor_transform <- rep("log1p", nrow(out))
  out$interpretation <- rep("Quality diagnostic only; not a false-positive classification", nrow(out))
  out
}

quality_sensitivity <- function(copy_matrix, traits, busco, tree, min_species, high_threshold, response_families = "") {
  family_map <- .quality_pgls$resolve_response_families(response_families, setdiff(names(traits), c("species", "busco_complete_pct")))
  rows <- list()
  membership <- list()
  for (trait in setdiff(names(traits), c("species", "busco_complete_pct"))) {
    for (family in colnames(copy_matrix)) {
      species <- rownames(copy_matrix)
      data <- data.frame(species = species, copy_number = copy_matrix[, family],
                         trait_value = traits[[trait]][match(species, traits$species)],
                         busco_quality = busco$busco_complete_pct[match(species, busco$species)])
      complete <- is.finite(data$trait_value) & is.finite(data$copy_number)
      masks <- list(baseline_all = complete,
                    baseline_busco_observed = complete & is.finite(data$busco_quality),
                    busco_adjusted = complete & is.finite(data$busco_quality),
                    high_completeness = complete & is.finite(data$busco_quality) & data$busco_quality >= high_threshold)
      for (variant in names(masks)) {
        row <- quality_fit(data[masks[[variant]], , drop = FALSE], tree, family, trait, min_species,
                           adjusted = variant == "busco_adjusted", response_family = family_map[[trait]])
        row$variant <- variant
        rows[[length(rows) + 1L]] <- row
      }
      # Counts are required to be complete for all families; cohort membership is trait-specific.
      if (identical(family, colnames(copy_matrix)[[1]])) {
        membership[[length(membership) + 1L]] <- data.frame(species = species, trait = trait,
          baseline_all = masks$baseline_all, busco_observed = masks$baseline_busco_observed,
          high_completeness = masks$high_completeness)
      }
    }
  }
  out <- if (length(rows)) do.call(rbind, rows) else .quality_pgls$empty_orthogroup_copy_number_trait_result()
  if (!"variant" %in% names(out)) out$variant <- character()
  out <- quality_adjust(out, "variant")
  out$predictor_transform <- rep("log1p", nrow(out))
  out$high_completeness_threshold <- rep(high_threshold, nrow(out))
  out$interpretation <- rep("Sensitivity analysis; loss of significance is not evidence of false detection", nrow(out))
  cohorts <- if (length(membership)) do.call(rbind, membership) else
    data.frame(species = character(), trait = character(), baseline_all = logical(), busco_observed = logical(), high_completeness = logical())
  list(results = out, cohorts = cohorts)
}

quality_annotate_cafe <- function(path, flags) {
  if (!nzchar(path)) return(data.frame())
  cafe <- .quality_pgls$read_tsv_base(path, col_classes = "character")
  if (!"#FamilyID" %in% names(cafe)) stop("Expected CAFE family_results with #FamilyID column.")
  if (anyDuplicated(cafe[["#FamilyID"]])) stop("Duplicate CAFE family IDs.")
  index <- match(cafe[["#FamilyID"]], flags$Orthogroup)
  cafe$busco_quality_flag <- flags$quality_flag[index]
  cafe$busco_quality_flag[is.na(index)] <- "not_assessed"
  cafe$busco_quality_qvalue <- flags$quality_qvalue[index]
  cafe$quality_diagnostic_scope <- rep("Extant family counts; not a CAFE branch test or false-positive classification", nrow(cafe))
  cafe
}

quality_annotate_pgls <- function(path, flags) {
  result <- .quality_pgls$read_tsv_base(path, col_classes = "character")
  if (!all(c("Orthogroup", "trait") %in% names(result))) stop("Expected PGLS Orthogroup and trait columns.")
  if (anyDuplicated(result[, c("Orthogroup", "trait"), drop = FALSE])) stop("Duplicate PGLS family/trait rows.")
  index <- match(result$Orthogroup, flags$Orthogroup)
  result$busco_quality_flag <- flags$quality_flag[index]
  result$busco_quality_flag[is.na(index)] <- "not_assessed"
  result$busco_quality_qvalue <- flags$quality_qvalue[index]
  result$quality_diagnostic_scope <- rep("BUSCO versus log1p extant counts; not a false-positive classification", nrow(result))
  result
}

run_copy_number_quality_diagnostics <- function(file_sptree, outdir, file_trait = "", file_busco = "",
    busco_short_dir = "", file_copy_number = "", file_cafe_results = "", trait_arg = "all",
    family_ids = "", family_file = "", max_families = "all", min_species = 4L, alpha = 0.05,
    high_threshold = 95, correlation_method = "spearman", sensitivity = TRUE, file_pgls_results = "", response_families = "") {
  if (nzchar(file_copy_number) && !nzchar(Sys.which("nwkit"))) stop("nwkit regress is required for copy-number quality diagnostics.")
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) stop("alpha must be between 0 and 1.")
  if (!is.finite(high_threshold) || high_threshold < 0 || high_threshold > 100) stop("high_threshold must be between 0 and 100.")
  if (length(min_species) != 1L || !is.finite(min_species) || min_species < 4L || min_species != floor(min_species)) stop("min_species must be an integer of at least 4.")
  work <- tempfile("gg-quality-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  tree <- .quality_pgls$load_tree_normalized(file_sptree)
  busco <- quality_load_busco(file_busco, busco_short_dir, tree$tip.label, work)
  trait <- if (nzchar(file_trait)) quality_read_table(file_trait) else data.frame(species = tree$tip.label)
  # Validate original labels and any reserved BUSCO column before selecting traits.
  if ("busco_complete_pct" %in% names(trait)) quality_numeric_traits(trait, busco, "busco_complete_pct")
  selection <- NULL
  observation_caption <- ""
  if (nzchar(file_trait)) {
    selected <- file.path(work, "selected_species_traits.tsv")
    status <- system2("python", shQuote(c(file.path(.quality_dir, "species_trait_contract.py"),
      "--input", file_trait, "--select", trait_arg, "--numeric-only", "--output", selected,
      "--report", file.path(work, "species_trait_input.json"),
      "--selection-report", file.path(work, "trait_selection.tsv"))))
    if (status != 0L) stop("Species-trait input contract failed; retain valid type/provenance sidecars and explicitly select observation metrics.")
    caption_flag <- system2("python", shQuote(c("-c",
      "import json,sys; print(int('gbif' in json.load(open(sys.argv[1]))))",
      file.path(work, "species_trait_input.json"))), stdout = TRUE)
    if (!is.null(attr(caption_flag, "status"))) stop("Could not read the species-trait interpretation report.")
    if (identical(caption_flag, "1")) observation_caption <- "\nGBIF values describe retained occurrence records, not true biological ranges."
    trait <- quality_read_table(selected)
    selection <- .quality_pgls$read_tsv_base(file.path(work, "trait_selection.tsv"), col_classes = "character")
  } else if (!tolower(trait_arg) %in% c("", "all")) stop("Explicit trait selection requires a trait table.")
  selected_arg <- if (nzchar(file_trait)) paste(setdiff(names(trait), "species"), collapse = ",") else "all"
  numeric_traits <- quality_numeric_traits(trait, busco, selected_arg)
  if (!is.null(selection)) {
    numeric_traits$excluded <- selection[selection$status != "selected", c("trait", "reason"), drop = FALSE]
  }
  pairs <- quality_correlations(numeric_traits$values, correlation_method)
  write <- .quality_pgls$write_tsv_base
  write(busco, file.path(work, "busco_quality.tsv"))
  write(numeric_traits$values, file.path(work, "traits_with_busco.tsv"))
  write(numeric_traits$excluded, file.path(work, "excluded_traits.tsv"))
  write(pairs, file.path(work, "trait_correlations.tsv"))
  quality_plot_correlations(pairs, work, observation_caption)
  flags <- .quality_pgls$empty_orthogroup_copy_number_trait_result()
  if (nzchar(file_copy_number)) {
    copies <- .quality_pgls$load_orthogroup_copy_number_matrix(file_copy_number, tree, family_ids, family_file, max_families)
    if (any(!is.finite(copies) | copies < 0 | copies != floor(copies))) stop("Copy numbers must be finite non-negative integers.")
    flags <- quality_family_flags(copies, busco, tree, min_species, alpha)
    write(flags, file.path(work, "family_busco_associations.tsv"))
    if (sensitivity) {
      sensitivities <- quality_sensitivity(copies, numeric_traits$values, busco, tree, min_species, high_threshold, response_families)
      match_flags <- match(sensitivities$results$Orthogroup, flags$Orthogroup)
      sensitivities$results$busco_quality_flag <- flags$quality_flag[match_flags]
      sensitivities$results$busco_quality_qvalue <- flags$quality_qvalue[match_flags]
      write(sensitivities$results, file.path(work, "trait_quality_sensitivity.tsv"))
      write(sensitivities$cohorts, file.path(work, "sensitivity_cohorts.tsv"))
    }
    if (nzchar(file_cafe_results)) write(quality_annotate_cafe(file_cafe_results, flags), file.path(work, "cafe_family_quality.tsv"))
    if (nzchar(file_pgls_results)) write(quality_annotate_pgls(file_pgls_results, flags), file.path(work, "pgls_quality.tsv"))
  } else if (nzchar(file_cafe_results) || nzchar(file_pgls_results)) stop("Result annotation requires copy-number data.")
  # Only manifest-listed members belong to this run; unrelated/previous files are preserved.
  unlink(file.path(work, "busco_input.tsv"))
  write(data.frame(parameter = c("alpha", "min_species", "high_threshold", "correlation_method", "predictor_transform", "sensitivity", "response_families", "classification"),
                   value = c(alpha, min_species, high_threshold, correlation_method, "log1p", sensitivity, response_families, "diagnostic_only")),
        file.path(work, "parameters.tsv"))
  write(data.frame(file = sort(list.files(work))), file.path(work, "manifest.tsv"))
  inputs <- c(file_sptree, file_trait, file_busco, busco_short_dir, file_copy_number, file_cafe_results, file_pgls_results, family_file)
  if (nzchar(file_trait)) {
    sidecars <- paste0(file_trait, c(".schema.json", ".metadata.json", ".gbif-quality.tsv", ".gbif-observations.tsv"))
    inputs <- c(inputs, sidecars[file.exists(sidecars)])
  }
  .quality_pgls$publish_copy_number_results(work, outdir, inputs[nzchar(inputs)])
  invisible(list(flags = flags, correlations = pairs))
}

quality_main <- function() {
  args <- .quality_pgls$parse_args(commandArgs(TRUE))
  str <- function(key, default = "") .quality_pgls$parse_string(args, key, default)
  num <- function(key, default) .quality_pgls$parse_numeric(args, key, default)
  if (!nzchar(str("file_sptree")) || !nzchar(str("outdir"))) stop("file_sptree and outdir are required.")
  run_copy_number_quality_diagnostics(str("file_sptree"), str("outdir"), str("file_trait"),
    str("file_busco"), str("busco_short_dir"), str("file_copy_number"), str("file_cafe_results"),
    str("trait", "all"), str("family_ids"), str("family_file"), str("max_families", "all"),
    num("min_species", 4L), num("alpha", 0.05), num("high_threshold", 95),
    str("correlation_method", "spearman"), .quality_pgls$parse_bool(args, "sensitivity", TRUE), str("file_pgls_results"), str("response_families"))
}

if (!identical(Sys.getenv("GG_COPY_NUMBER_QUALITY_NO_MAIN"), "1")) quality_main()
