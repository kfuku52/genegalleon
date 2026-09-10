script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
repo <- normalizePath(file.path(dirname(script), "..", ".."))
Sys.setenv(GG_COPY_NUMBER_QUALITY_NO_MAIN = "1")
source(file.path(repo, "workflow/support/copy_number_quality_diagnostics.r"))
tmp <- tempfile("gg-quality-test-")
dir.create(tmp)

expect_error <- function(expr, pattern) {
  error <- tryCatch({force(expr); NULL}, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}
write <- .quality_pgls$write_tsv_base
tree_path <- file.path(tmp, "tree.nwk")
writeLines("(((sp1:1,sp2:1):1,(sp3:1,sp4:1):1):1,((sp5:1,sp6:1):1,(sp7:1,sp8:1):1):1);", tree_path)
tree <- ape::read.tree(tree_path)
species <- tree$tip.label
copies <- cbind(OG1 = 1:8, OG_constant = rep(2, 8))
rownames(copies) <- species
busco <- data.frame(species = species, busco_complete_pct = 80 + 8 * log1p(1:8) + c(.1,-.2,.1,.2,-.1,-.2,.1,.1),
                    lineage = "eukaryota_odb12", mode = "proteins", busco_version = "6.0.0")
busco_path <- file.path(tmp, "busco.tsv")
write(busco[8:1, ], busco_path)
loaded <- quality_load_busco(busco_path, "", species, tmp)
stopifnot(identical(loaded$species, species), isTRUE(all.equal(loaded$busco_complete_pct, busco$busco_complete_pct)))
empty <- quality_load_busco("", "", species, tmp)
stopifnot(all(is.na(empty$busco_complete_pct)), all(empty$quality_status == "missing"))
mixed <- busco
mixed$lineage[[1]] <- "metazoa_odb12"
write(mixed, file.path(tmp, "mixed.tsv"))
expect_error(quality_load_busco(file.path(tmp, "mixed.tsv"), "", species, tmp), "Mixed BUSCO lineage")

traits <- data.frame(species = species,
  size = 1 + 1.2 * log1p(1:8) + .2 * busco$busco_complete_pct + c(.2,-.1,.3,-.1,.1,-.3,.1,-.2),
  binary = rep(0:1, 4), constant = 1, habitat = "forest")
numeric_traits <- quality_numeric_traits(traits, loaded)
stopifnot(identical(numeric_traits$excluded$trait, "habitat"))
pairs <- quality_correlations(numeric_traits$values)
pair <- pairs[pairs$trait_x == "size" & pairs$trait_y == "busco_complete_pct", ]
stopifnot(pair$n_species == 8, isTRUE(all.equal(pair$correlation, cor(traits$size, busco$busco_complete_pct, method = "spearman"))))
stopifnot(all(is.na(pairs$correlation[pairs$trait_x == "constant"])))
missing_traits <- numeric_traits$values
missing_traits$size[1:6] <- NA_real_
missing_pairs <- quality_correlations(missing_traits)
stopifnot(all(missing_pairs$status[missing_pairs$trait_x == "size"] == "too_few_species"))

flags <- quality_family_flags(copies, loaded, tree)
stopifnot(flags$quality_flag[flags$Orthogroup == "OG1"] == "quality_associated",
          flags$quality_flag[flags$Orthogroup == "OG_constant"] == "not_assessable")
stopifnot(isTRUE(all.equal(flags$quality_qvalue[[1]], p.adjust(flags$pval[[1]], "BH", n = 2))))
missing_flags <- quality_family_flags(copies, empty, tree)
stopifnot(all(missing_flags$quality_flag == "not_assessable"))

# Verify adjusted coefficients/SE against independent matrix GLS, not another adapter call.
data <- data.frame(species = species, trait_value = traits$size, copy_number = 1:8, busco_quality = busco$busco_complete_pct)
fit <- quality_fit(data, tree, "OG1", "size", 4, adjusted = TRUE)
X <- cbind(1, log1p(1:8), busco$busco_complete_pct)
W <- solve(ape::vcv.phylo(tree))
beta <- solve(t(X) %*% W %*% X, t(X) %*% W %*% traits$size)
resid <- traits$size - as.vector(X %*% beta)
variance <- as.numeric(t(resid) %*% W %*% resid) / 5
se <- sqrt(diag(variance * solve(t(X) %*% W %*% X)))[[2]]
stopifnot(fit$status == "ok", isTRUE(all.equal(fit$coefficient, beta[[2]], tolerance = 1e-7)),
          isTRUE(all.equal(fit$standard_error, se, tolerance = 1e-7)))
data$busco_quality <- log1p(data$copy_number)
rank_fit <- quality_fit(data, tree, "OG1", "size", 4, adjusted = TRUE)
stopifnot(rank_fit$status == "not_estimable", rank_fit$skip_reason == "rank_deficient_or_no_residual_df")

# Restricting to one clade must retain Brownian covariance from the original root.
subset <- .quality_pgls$subset_tree_to_species(tree, species[1:4])
stopifnot(isTRUE(all.equal(ape::vcv.phylo(subset), ape::vcv.phylo(tree)[species[1:4], species[1:4]])))
busco_missing <- loaded
busco_missing$busco_complete_pct[[8]] <- NA_real_
sensitivity <- quality_sensitivity(copies[, "OG1", drop = FALSE], numeric_traits$values[, c("species", "size")],
                                   busco_missing, tree, 4, 95)
stopifnot(sensitivity$results$n_species[sensitivity$results$variant == "baseline_all"] == 8,
          sensitivity$results$n_species[sensitivity$results$variant == "baseline_busco_observed"] == 7,
          sensitivity$results$n_species[sensitivity$results$variant == "busco_adjusted"] == 7,
          sensitivity$results$status[sensitivity$results$variant == "high_completeness"] == "skipped")

cafe_path <- file.path(tmp, "Gamma_family_results.txt")
writeLines(c("#FamilyID\tpvalue\tSignificant at 0.05", "OG1\t0.001\ty", "OG_missing\t0.04\ty"), cafe_path)
annotated <- quality_annotate_cafe(cafe_path, flags)
stopifnot(nrow(annotated) == 2, identical(annotated[["Significant at 0.05"]], c("y", "y")),
          identical(annotated$busco_quality_flag, c("quality_associated", "not_assessed")))
pgls_path <- file.path(tmp, "original_pgls.tsv")
write(data.frame(Orthogroup = c("OG1", "OG1"), trait = c("size", "binary"), p.adj.global = c(.01, .2)), pgls_path)
annotated_pgls <- quality_annotate_pgls(pgls_path, flags)
stopifnot(nrow(annotated_pgls) == 2, identical(annotated_pgls$p.adj.global, c("0.01", "0.2")),
          all(annotated_pgls$busco_quality_flag == "quality_associated"))

copy_path <- file.path(tmp, "copies.tsv")
write(data.frame(Description = "test", Orthogroup = colnames(copies), t(copies), check.names = FALSE), copy_path)
trait_path <- file.path(tmp, "traits.tsv")
write(traits, trait_path)
# Formal type schema excludes category text and numeric-looking category codes.
status <- system2("python", shQuote(c("-c", paste0(
  "import sys; sys.path.insert(0,sys.argv[1]); from pathlib import Path; ",
  "from species_trait_schema import schema_payload; p=Path(sys.argv[2]); ",
  "Path(str(p)+'.schema.json').write_bytes(schema_payload(p.read_bytes(), ",
  "{'size':'numeric','binary':'binary','constant':'numeric','habitat':'categorical'}))"),
  file.path(repo, "workflow/support"), trait_path)))
stopifnot(status == 0L)
# Preserve leading zeros in identifiers and original annotation fields.
writeLines(c("species\tbusco_complete_pct", "001\t98.5"), file.path(tmp, "leading.tsv"))
stopifnot(quality_read_table(file.path(tmp, "leading.tsv"))$species == "001")
writeLines(c("#FamilyID\tpvalue", "001\t1.00e-03"), file.path(tmp, "leading_cafe.tsv"))
leading <- quality_annotate_cafe(file.path(tmp, "leading_cafe.tsv"), flags)
stopifnot(leading[["#FamilyID"]] == "001", leading$pvalue == "1.00e-03")
# Empty planned sets produce valid empty tables rather than a scalar assignment error.
stopifnot(nrow(quality_family_flags(copies[, FALSE, drop = FALSE], loaded, tree)) == 0)
stopifnot(nrow(quality_sensitivity(copies, numeric_traits$values[, "species", drop = FALSE], loaded, tree, 4, 95)$results) == 0)
stopifnot(nrow(quality_sensitivity(copies[, FALSE, drop = FALSE], numeric_traits$values, loaded, tree, 4, 95)$results) == 0)
# Each sensitivity variant carries the requested response family, including skipped fits.
binary_sensitivity <- quality_sensitivity(copies[, "OG1", drop = FALSE],
  numeric_traits$values[, c("species", "binary")], busco_missing, tree, 4, 95, "binary=binomial")
stopifnot(all(binary_sensitivity$results$response_family == "binomial"),
          all(binary_sensitivity$results$link_function == "logit"),
          !any(binary_sensitivity$results$status == "error"))
# A missing fitting executable is a workflow failure, not an unassessable biological result.
saved_path <- Sys.getenv("PATH")
Sys.setenv(PATH = file.path(tmp, "no-executables"))
expect_error(run_copy_number_quality_diagnostics(tree_path, file.path(tmp, "missing-nwkit"),
  file_copy_number = copy_path), "nwkit regress is required")
Sys.setenv(PATH = saved_path)
out <- file.path(tmp, "out")
run_copy_number_quality_diagnostics(tree_path, out, trait_path, busco_path,
  file_copy_number = copy_path, file_cafe_results = cafe_path, trait_arg = "size", sensitivity = TRUE,
  file_pgls_results = pgls_path)
manifest <- read.delim(file.path(out, "manifest.tsv"))
stopifnot(all(file.exists(file.path(out, manifest$file))), "parameters.tsv" %in% manifest$file,
          file.info(file.path(out, "trait_correlations.svg"))$size > 1000)
# Late failure must not replace the previously published analysis.
before <- tools::md5sum(file.path(out, manifest$file))
saved_plot <- quality_plot_correlations
quality_plot_correlations <- function(...) stop("injected plot failure")
expect_error(run_copy_number_quality_diagnostics(tree_path, out, trait_path, busco_path), "injected plot failure")
quality_plot_correlations <- saved_plot
stopifnot(identical(before, tools::md5sum(file.path(out, manifest$file))))
preview <- Sys.getenv("GG_QUALITY_TEST_OUTPUT")
if (nzchar(preview)) {
  dir.create(preview, recursive = TRUE, showWarnings = FALSE)
  file.copy(list.files(out, full.names = TRUE), preview, overwrite = TRUE)
}
unlink(tmp, recursive = TRUE)
cat("BUSCO quality diagnostics tests passed.\n")
