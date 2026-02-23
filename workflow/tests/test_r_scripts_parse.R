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

workflow_script_dir <- file.path(repo_root, "workflow", "support")

if (!dir.exists(workflow_script_dir)) {
  stop(sprintf("Workflow script directory not found: %s", workflow_script_dir))
}

r_files <- sort(list.files(
  workflow_script_dir,
  pattern = "[.][rR]$",
  recursive = TRUE,
  full.names = TRUE
))

if (length(r_files) == 0) {
  stop(sprintf("No R scripts found under: %s", workflow_script_dir))
}

parse_failures <- character(0)
for (r_file in r_files) {
  tryCatch(
    parse(file = r_file),
    error = function(err) {
      message(sprintf("PARSE_FAIL: %s :: %s", r_file, conditionMessage(err)))
      parse_failures <<- c(parse_failures, r_file)
    }
  )
}

if (length(parse_failures) > 0) {
  stop(sprintf("R parse failed for %d script(s).", length(parse_failures)))
}

cat(sprintf("Parsed %d R scripts successfully.\n", length(r_files)))
