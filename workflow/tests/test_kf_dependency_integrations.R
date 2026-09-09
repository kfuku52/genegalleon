args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_arg <- args[startsWith(args, file_arg)]
if (length(script_arg) != 1L) {
  stop("Could not determine test script path from commandArgs().")
}

test_path <- normalizePath(sub(file_arg, "", script_arg[[1]], fixed = TRUE), mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(test_path), "..", ".."), mustWork = TRUE)
support_dir <- file.path(repo_root, "workflow", "support")

run_r_script <- function(script_name, script_args, work_dir) {
  old_wd <- setwd(work_dir)
  on.exit(setwd(old_wd), add = TRUE)
  output <- suppressWarnings(system2(
    "Rscript",
    c(file.path(support_dir, script_name), script_args),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  if (status != 0L) {
    stop(script_name, " failed:\n", paste(output, collapse = "\n"))
  }
  output
}

tree_text <- "(((sp1:1,sp2:1):1,(sp3:1,sp4:1):1):1,((sp5:1,sp6:1):1,(sp7:1,sp8:1):1):1);"

ou_dir <- tempfile("gg-kfl1ou-")
dir.create(ou_dir)
writeLines(tree_text, file.path(ou_dir, "tree.nwk"))
ou_traits <- data.frame(
  root_1 = c(1.0, 1.2, 1.1, 1.3, 3.8, 4.0, 4.1, 3.9),
  root_2 = c(1.1, 1.3, 1.0, 1.4, 3.9, 4.2, 4.0, 3.8),
  leaf_1 = c(2.0, 2.1, 2.2, 2.0, 4.0, 4.2, 4.1, 4.3),
  leaf_2 = c(2.1, 2.2, 2.3, 1.9, 4.1, 4.1, 4.2, 4.4),
  row.names = paste0("sp", 1:8),
  check.names = FALSE
)
write.table(
  ou_traits,
  file.path(ou_dir, "traits.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
ou_output <- run_r_script(
  "detect_OU_shift_kfl1ou.r",
  c(
    "--tree_file=tree.nwk",
    "--trait_file=traits.tsv",
    "--max_nshift=1",
    "--nbootstrap=0",
    "--nslots=1",
    "--criterion=pBIC",
    "--detect_convergence=0",
    "--fit_ind_file=",
    "--fit_conv_file="
  ),
  ou_dir
)
if (!any(grepl("kfl1ou completed!", ou_output, fixed = TRUE))) {
  stop("The direct kfl1ou wrapper did not complete.")
}
ou_files <- c("fit_ind.RData", "l1ou_tree.tsv", "l1ou_regime.tsv", "l1ou_leaf.tsv", "l1ou_plot.pdf")
if (!all(file.exists(file.path(ou_dir, ou_files)))) {
  stop("The direct kfl1ou wrapper did not create every required output.")
}

ou_empty_dir <- tempfile("gg-kfl1ou-empty-")
dir.create(ou_empty_dir)
writeLines(tree_text, file.path(ou_empty_dir, "tree.nwk"))
writeLines(c("gene_id", paste0("sp", 1:8)), file.path(ou_empty_dir, "traits.tsv"))
ou_empty_output <- run_r_script(
  "detect_OU_shift_kfl1ou.r",
  c(
    "--tree_file=tree.nwk",
    "--trait_file=traits.tsv",
    "--max_nshift=1",
    "--nbootstrap=1",
    "--nslots=1",
    "--criterion=pBIC",
    "--detect_convergence=0",
    "--fit_ind_file=",
    "--fit_conv_file="
  ),
  ou_empty_dir
)
if (!any(grepl("No trait columns were found", ou_empty_output, fixed = TRUE))) {
  stop("The empty-trait kfl1ou path did not report placeholder generation.")
}
if (!all(file.exists(file.path(ou_empty_dir, ou_files)))) {
  stop("The empty-trait kfl1ou path did not create every required placeholder output.")
}
ou_empty_bootstrap <- read.delim(file.path(ou_empty_dir, "l1ou_bootstrap.tsv"), check.names = FALSE)
if (nrow(ou_empty_bootstrap) != 15L || sum(is.na(ou_empty_bootstrap$bootstrap_support)) != 1L) {
  stop("The empty-trait kfl1ou path did not create a valid placeholder bootstrap table.")
}

# PGLS is exercised by the NWKIT runtime and copy-number integration tests.
unlink(c(ou_dir, ou_empty_dir), recursive = TRUE)
cat("test_kf_dependency_integrations.R: OK\n")
