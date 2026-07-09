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
    "sp2\t2\t0\t5",
    "sp3\t3\t1\t5",
    "sp4\t4\t1\t5"
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

mock_phylopars <- function(...) {
  structure(
    list(R2 = 0.8, R2adj = 0.7, sigma = 0.1, Fstat = 10, pval = 0.01, logLik = -2),
    class = "mock_copy_number_phylopars"
  )
}
AIC.mock_copy_number_phylopars <- function(object, ...) 12
BIC.mock_copy_number_phylopars <- function(object, ...) 14

df_stat <- run_orthogroup_copy_number_trait_associations(
  copy_matrix = copy_matrix,
  trait = trait,
  tree = tree,
  trait_cols = c("height", "constant_trait"),
  min_species = 4,
  p_adjust_method = "BH",
  fit_fun = mock_phylopars,
  verbose = FALSE
)

og1_height <- df_stat[df_stat$Orthogroup == "OG1" & df_stat$trait == "height", , drop = FALSE]
stopifnot(nrow(og1_height) == 1)
stopifnot(identical(og1_height$status, "ok"))
stopifnot(isTRUE(all.equal(og1_height$R2, 0.8)))
stopifnot(isTRUE(all.equal(og1_height$pval, 0.01)))
stopifnot(isTRUE(all.equal(og1_height$PCC, 1)))

og10_height <- df_stat[df_stat$Orthogroup == "OG10" & df_stat$trait == "height", , drop = FALSE]
stopifnot(nrow(og10_height) == 1)
stopifnot(identical(og10_height$status, "skipped"))
stopifnot(identical(og10_height$skip_reason, "invariant_copy_number"))

og1_constant <- df_stat[df_stat$Orthogroup == "OG1" & df_stat$trait == "constant_trait", , drop = FALSE]
stopifnot(nrow(og1_constant) == 1)
stopifnot(identical(og1_constant$status, "skipped"))
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
  fit_fun = mock_phylopars,
  verbose = FALSE
)

stopifnot(file.exists(file.path(outdir, "orthogroup_copy_number_matrix.tsv")))
stopifnot(file.exists(file.path(outdir, "orthogroup_copy_number_trait_pgls.tsv")))
stopifnot(file.exists(file.path(outdir, "orthogroup_copy_number_trait_pgls.significant.tsv")))
stopifnot(file.exists(file.path(outdir, "orthogroup_copy_number_trait_pgls.summary.pdf")))

written_stats <- read.delim(file.path(outdir, "orthogroup_copy_number_trait_pgls.tsv"), sep = "\t", check.names = FALSE)
stopifnot(identical(written_stats$Orthogroup, "OG1"))
stopifnot(identical(written_stats$trait, "height"))
stopifnot(identical(written_stats$status, "ok"))

cat("test_orthogroup_copy_number_trait_pgls.R: OK\n")
