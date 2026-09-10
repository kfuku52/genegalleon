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

resolved_script_path <- normalizePath(script_path, winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(resolved_script_path), "..", ".."), winslash = "/", mustWork = TRUE)

Sys.setenv(GG_ORTHOGROUP_COPY_NUMBER_TRAIT_PGLS_NO_MAIN = "1")
source(file.path(repo_root, "workflow", "support", "orthogroup_copy_number_trait_pgls.r"))

tmp <- tempfile("gg_orthogroup_copy_number_trait_pgls_")
dir.create(tmp, recursive = TRUE)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

tree_file <- file.path(tmp, "tree.nwk")
writeLines("((sp1:1,sp2:1):1,(sp3:1,sp4:1):1);", tree_file)

genecount_file <- file.path(tmp, "Orthogroups.GeneCount.selected.tsv")
writeLines(
  c(
    "besthit_0.95\tOrthogroup\tsp1\tsp2\tsp3\tsp4",
    "hit1\tOG1\t1\t2\t3\t4",
    "hit2\tOG10\t9\t9\t9\t9",
    "hit3\tOG2\t4\t3\t2\t1",
    "hit4\tOG_TOO_WIDE\t1\t1\t1\t99"
  ),
  genecount_file
)

python_bin <- Sys.which("python")
if (!nzchar(python_bin)) {
  python_bin <- Sys.which("python3")
}
if (!nzchar(python_bin)) {
  stop("python or python3 is required for the orthogroup copy-number preparation integration test.")
}
prepare_dir <- file.path(tmp, "prepared_copy_number")
prepare_script <- file.path(repo_root, "workflow", "support", "prepare_orthogroup_copy_number.py")
prepare_out <- system2(
  python_bin,
  c(
    prepare_script,
    "--genecount", genecount_file,
    "--dated_species_tree", tree_file,
    "--output_dir", prepare_dir,
    "--max_size_differential", "10"
  ),
  stdout = TRUE,
  stderr = TRUE
)
prepare_status <- attr(prepare_out, "status")
if (!is.null(prepare_status) && prepare_status != 0) {
  stop(paste(c("prepare_orthogroup_copy_number.py failed:", prepare_out), collapse = "\n"))
}

copy_number_file <- file.path(prepare_dir, "orthogroup_copy_number.tsv")
removed_file <- file.path(prepare_dir, "removed_orthogroups.tsv")
stopifnot(file.exists(copy_number_file))
stopifnot(file.exists(removed_file))
prepared_copy_number <- read.delim(copy_number_file, sep = "\t", check.names = FALSE)
removed_orthogroups <- read.delim(removed_file, sep = "\t", check.names = FALSE)
stopifnot(identical(prepared_copy_number$Orthogroup, c("OG1", "OG10", "OG2")))
stopifnot(identical(removed_orthogroups$Orthogroup, "OG_TOO_WIDE"))

trait_file <- file.path(tmp, "species_trait.tsv")
writeLines(
  c(
    "species\theight\tbinary_trait\tconstant_trait",
    "sp1\t1\t0\t5",
    "sp2\t3\t0\t5",
    "sp3\t2\t1\t5",
    "sp4\t6\t1\t5"
  ),
  trait_file
)

tree <- load_tree_normalized(tree_file)
copy_matrix <- load_orthogroup_copy_number_matrix(copy_number_file, tree, family_ids = "OG1 OG10", max_families = "all")
stopifnot(identical(colnames(copy_matrix), c("OG1", "OG10")))
stopifnot(identical(rownames(copy_matrix), c("sp1", "sp2", "sp3", "sp4")))
stopifnot(identical(as.numeric(copy_matrix[, "OG1"]), c(1, 2, 3, 4)))
stopifnot(identical(as.numeric(copy_matrix[, "OG10"]), c(9, 9, 9, 9)))

family_file <- file.path(tmp, "families.tsv")
writeLines(c("family_id", "OG2"), family_file)
copy_matrix_file <- load_orthogroup_copy_number_matrix(copy_number_file, tree, family_file = family_file, max_families = "all")
stopifnot(identical(colnames(copy_matrix_file), "OG2"))

trait <- load_trait_table(trait_file)
stopifnot(identical(resolve_trait_cols(trait, "all"), c("height", "binary_trait", "constant_trait")))
stopifnot(identical(resolve_trait_cols(trait, "height,binary_trait"), c("height", "binary_trait")))

df_stat <- run_orthogroup_copy_number_trait_associations(
  copy_matrix = copy_matrix,
  trait = trait,
  tree = tree,
  trait_cols = c("height", "constant_trait"),
  min_species = 4,
  p_adjust_method = "BH",
  verbose = FALSE
)

og1_height <- df_stat[df_stat$Orthogroup == "OG1" & df_stat$trait == "height", , drop = FALSE]
stopifnot(nrow(og1_height) == 1)
stopifnot(identical(og1_height$status, "ok"),
          identical(og1_height$predictor_transform, "log1p"))
# Independent GLS calculation using ape's Brownian covariance.
X <- cbind(1, log1p(1:4))
y <- c(1, 3, 2, 6)
V <- ape::vcv.phylo(tree)
precision <- solve(V)
beta <- solve(t(X) %*% precision %*% X, t(X) %*% precision %*% y)
residual <- y - as.vector(X %*% beta)
rate <- as.numeric(t(residual) %*% precision %*% residual) / (length(y) - ncol(X))
standard_error <- sqrt(diag(rate * solve(t(X) %*% precision %*% X)))[[2]]
p_value <- 2 * pt(-abs(beta[[2]] / standard_error), df = length(y) - ncol(X))
interval <- beta[[2]] + c(-1, 1) * qt(0.975, df = length(y) - ncol(X)) * standard_error
stopifnot(isTRUE(all.equal(og1_height$coefficient, beta[[2]], tolerance = 1e-10)),
          isTRUE(all.equal(og1_height$standard_error, standard_error, tolerance = 1e-10)),
          isTRUE(all.equal(og1_height$pval, p_value, tolerance = 1e-10)),
          isTRUE(all.equal(og1_height$confidence_interval_lower, interval[[1]], tolerance = 1e-10)),
          isTRUE(all.equal(og1_height$confidence_interval_upper, interval[[2]], tolerance = 1e-10)),
          identical(og1_height$confidence_level, 0.95),
          identical(og1_height$fit_mode, "nwkit_brownian_reml"),
          identical(og1_height$covariance_estimator, "gaussian-REML"))

og10_height <- df_stat[df_stat$Orthogroup == "OG10" & df_stat$trait == "height", , drop = FALSE]
stopifnot(nrow(og10_height) == 1)
stopifnot(identical(og10_height$status, "skipped"))
stopifnot(identical(og10_height$skip_reason, "invariant_copy_number"))

og1_constant <- df_stat[df_stat$Orthogroup == "OG1" & df_stat$trait == "constant_trait", , drop = FALSE]
stopifnot(nrow(og1_constant) == 1)
stopifnot(identical(og1_constant$status, "skipped"))
stopifnot(is.na(og1_constant$confidence_interval_lower), is.na(og1_constant$confidence_interval_upper))
stopifnot(identical(og1_constant$skip_reason, "invariant_trait"))

outdir <- file.path(tmp, "out")
run_orthogroup_copy_number_trait_pgls(
  file_orthogroup_copy_number = copy_number_file,
  file_sptree = tree_file,
  file_trait = trait_file,
  outdir = outdir,
  trait_arg = "height",
  min_species = 4,
  family_ids = "OG1",
  max_families = "all",
  verbose = FALSE
)

stopifnot(file.exists(file.path(outdir, "orthogroup_copy_number_matrix.tsv")))
stopifnot(file.exists(file.path(outdir, "orthogroup_copy_number_trait_pgls.tsv")))
stopifnot(file.exists(file.path(outdir, "orthogroup_copy_number_trait_pgls.significant.tsv")))
stopifnot(file.exists(file.path(outdir, "orthogroup_copy_number_trait_pgls.summary.pdf")))

written_stats <- read.delim(file.path(outdir, "orthogroup_copy_number_trait_pgls.tsv"), sep = "\t", check.names = FALSE)
written_significant <- read.delim(file.path(outdir, "orthogroup_copy_number_trait_pgls.significant.tsv"), sep = "\t", check.names = FALSE)
stopifnot(identical(names(written_stats), names(written_significant)))
stopifnot(identical(written_stats$Orthogroup, "OG1"))
stopifnot(identical(written_stats$trait, "height"))
stopifnot(identical(written_stats$status, "ok"))
stopifnot(isTRUE(all.equal(written_stats$confidence_interval_lower, interval[[1]], tolerance = 1e-10)),
          isTRUE(all.equal(written_stats$confidence_interval_upper, interval[[2]], tolerance = 1e-10)),
          identical(written_stats$confidence_level, 0.95))

# A failed late plot must leave every member of the previous bundle intact.
output_names <- c("orthogroup_copy_number_matrix.tsv", "orthogroup_copy_number_trait_pgls.tsv",
                  "orthogroup_copy_number_trait_pgls.significant.tsv",
                  "orthogroup_copy_number_trait_pgls.summary.pdf", "orthogroup_copy_number_trait_pgls.summary.svg")
for (name in output_names) writeLines(paste("previous", name), file.path(outdir, name))
previous <- lapply(file.path(outdir, output_names), readBin, what = "raw", n = 10000)
original_save_summary_plot <- save_summary_plot
save_summary_plot <- function(...) stop("injected plot failure")
failed <- tryCatch(run_orthogroup_copy_number_trait_pgls(
  file_orthogroup_copy_number = copy_number_file, file_sptree = tree_file,
  file_trait = trait_file, outdir = outdir, trait_arg = "height", family_ids = "OG1"
), error = identity)
save_summary_plot <- original_save_summary_plot
stopifnot(inherits(failed, "error"), grepl("injected plot failure", conditionMessage(failed)))
stopifnot(identical(previous, lapply(file.path(outdir, output_names), readBin, what = "raw", n = 10000)))

# Output names must not overwrite a supplied input, even through an alias.
alias_out <- file.path(tmp, "alias output")
dir.create(alias_out)
alias_input <- file.path(alias_out, "orthogroup_copy_number_trait_pgls.tsv")
file.copy(copy_number_file, alias_input)
alias_before <- readBin(alias_input, "raw", n = 100000)
alias_failure <- tryCatch(run_orthogroup_copy_number_trait_pgls(
  file_orthogroup_copy_number = alias_input, file_sptree = tree_file,
  file_trait = trait_file, outdir = alias_out, trait_arg = "height", family_ids = "OG1"
), error = identity)
stopifnot(inherits(alias_failure, "error"),
          identical(alias_before, readBin(alias_input, "raw", n = 100000)),
          length(list.files(alias_out)) == 1L)

# Destination errors must also preserve the other members of the old bundle.
blocked_name <- file.path(outdir, "orthogroup_copy_number_trait_pgls.summary.svg")
unlink(blocked_name)
dir.create(blocked_name)
blocked_before <- lapply(file.path(outdir, head(output_names, -1)), readBin, what = "raw", n = 10000)
blocked_failure <- tryCatch(run_orthogroup_copy_number_trait_pgls(
  file_orthogroup_copy_number = copy_number_file, file_sptree = tree_file,
  file_trait = trait_file, outdir = outdir, trait_arg = "height", family_ids = "OG1"
), error = identity)
stopifnot(inherits(blocked_failure, "error"), dir.exists(blocked_name),
          identical(blocked_before, lapply(file.path(outdir, head(output_names, -1)), readBin, what = "raw", n = 10000)))

cat("test_orthogroup_copy_number_trait_pgls.R: OK\n")

# Family routing uses the actual NWKIT engine, including numeric 0/1 coding.
stopifnot(identical(unname(resolve_response_families("binary_trait=binomial,height=poisson", c("height", "binary_trait"))), c("poisson", "binomial")))
for (bad in c("unknown=binomial", "height=bad", "height=poisson,height=gaussian")) {
  stopifnot(inherits(tryCatch(resolve_response_families(bad, "height"), error = identity), "error"))
}
stopifnot(inherits(tryCatch(validate_response_values(c(0, 2), "binomial"), error = identity), "error"))
stopifnot(inherits(tryCatch(validate_response_values(c(0, 1.5), "poisson"), error = identity), "error"))
mixed <- run_orthogroup_copy_number_trait_associations(
  copy_matrix, load_trait_table(trait_file), tree, c("height", "binary_trait"),
  response_families = "height=poisson,binary_trait=binomial"
)
stopifnot(all(mixed$response_family[mixed$trait == "height"] == "poisson"))
stopifnot(all(mixed$response_family[mixed$trait == "binary_trait"] == "binomial"))
stopifnot(all(mixed$link_function[mixed$trait == "binary_trait"] == "logit"))
stopifnot(!any(mixed$status == "error"))
usable_intervals <- mixed$status == "ok" & is.finite(mixed$standard_error)
stopifnot(any(usable_intervals),
          isTRUE(all.equal(mixed$confidence_interval_lower[usable_intervals],
                           mixed$coefficient[usable_intervals] - qnorm(0.975) * mixed$standard_error[usable_intervals],
                           tolerance = 1e-10)),
          isTRUE(all.equal(mixed$confidence_interval_upper[usable_intervals],
                           mixed$coefficient[usable_intervals] + qnorm(0.975) * mixed$standard_error[usable_intervals],
                           tolerance = 1e-10)))
cat("mixed-family regression: OK\n")

stopifnot(inherits(tryCatch(parse_response_values(c("0", "oops"), "binomial"), error = identity), "error"))
stopifnot(inherits(tryCatch(parse_response_values(c("1", "Inf"), "poisson"), error = identity), "error"))
stopifnot(identical(parse_response_values(c("0", "1", "NA"), "binomial"), c(0, 1, NA_real_)))

# Missing responses confined to one root clade must preserve shared history.
rooted_subset_source <- ape::read.tree(text = "(((a:1.123456789012345,b:1):1,(c:1,d:1):1):2,(e:1,f:1):3);")
rooted_subset_names <- c("a", "b", "c", "d")
rooted_subset <- subset_tree_to_species(rooted_subset_source, rooted_subset_names)
stopifnot(isTRUE(all.equal(
  ape::vcv(rooted_subset),
  ape::vcv(rooted_subset_source)[rooted_subset_names, rooted_subset_names],
  tolerance = 1e-14
)))
