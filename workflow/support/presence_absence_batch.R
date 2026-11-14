# Render related views in one R process; retain only shared parsed inputs.
gg_render_presence_batch <- function(arguments, renderer) {
  separator <- which(arguments == "--next-plot")
  boundaries <- c(0L, separator, length(arguments) + 1L)
  jobs <- lapply(seq_len(length(boundaries) - 1L), function(i) {
    start <- boundaries[[i]] + 1L
    end <- boundaries[[i + 1L]] - 1L
    if (start > end) stop("Empty presence/absence plot job")
    arguments[start:end]
  })
  if (length(jobs) > 3L) stop("Expected at most three presence/absence views")
  cache <- new.env(parent = emptyenv())
  observed <- new.env(parent = emptyenv())
  signature <- function(path, verify_bytes = TRUE) {
    info <- file.info(path)[, c("size", "mtime", "ctime"), drop = FALSE]
    if (dir.exists(path)) return(list(info = info, entries = list.files(path, all.files = TRUE, no.. = TRUE)))
    result <- list(info = info)
    if (verify_bytes) result$digest <- unname(tools::md5sum(path))
    result
  }
  observe <- function(path, verify_bytes = FALSE) {
    path <- normalizePath(path, mustWork = TRUE)
    if (exists(path, observed, inherits = FALSE)) {
      expected <- observed[[path]]
      if (!verify_bytes) expected$digest <- NULL
      if (!identical(expected, signature(path, verify_bytes))) stop("Plot input changed during batch: ", path)
    } else {
      observed[[path]] <- signature(path)
    }
    path
  }
  check <- function(verify_bytes = FALSE) {
    for (path in ls(observed, all.names = TRUE)) observe(path, verify_bytes)
  }
  read_shared <- function(path, kind, reader, retain = TRUE) {
    path <- observe(path)
    key <- paste(kind, path, sep = "\n")
    value <- if (retain && exists(key, cache, inherits = FALSE)) cache[[key]] else reader(path)
    if (retain) cache[[key]] <- value
    observe(path)
    value
  }
  renderer <- observe(renderer)
  expressions <- parse(renderer)
  outputs <- character()
  staged <- character()
  backups <- character()
  installed <- character()
  render_messages <- character()
  committed <- FALSE
  on.exit({
    retained <- character()
    if (!committed) for (path in installed) {
      if (path %in% names(backups)) {
        if (!file.rename(backups[[path]], path)) {
          retained <- c(retained, backups[[path]])
          warning("Could not restore plot; previous bytes retained at: ", backups[[path]])
        }
      } else unlink(path)
    }
    unlink(c(staged, setdiff(backups, retained)))
  }, add = TRUE)
  for (job in jobs) {
    inputs <- grepl("^--(species_tree|species_tree_ci|support_tree|species_mapping_tree|busco_table|long_table|label_map|taxonomy_lineage|ortholog_[a-z_]+_table)=", job)
    for (argument in job[inputs]) {
      path <- sub("^--[^=]+=", "", argument)
      if (file.exists(path) && !dir.exists(path)) observe(path)
    }
    for (index in which(grepl("^--out_(pdf|svg)=.", job))) {
      output <- sub("^--[^=]+=", "", job[[index]])
      dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
      output <- file.path(normalizePath(dirname(output), mustWork = TRUE), basename(output))
      if (output %in% outputs || dir.exists(output) || output %in% ls(observed, all.names = TRUE)) {
        stop("Duplicate output or plot output replaces an input: ", output)
      }
      stage <- tempfile(".gg-presence-", tmpdir = dirname(output), fileext = paste0(".", tools::file_ext(output)))
      outputs <- c(outputs, output)
      staged <- c(staged, stage)
      job[[index]] <- paste0(sub("=.*$", "=", job[[index]]), stage)
    }
    messages <- local({
      before_options <- options()
      before_theme <- ggplot2::theme_get()
      on.exit({
        while (grDevices::dev.cur() > 1L) grDevices::dev.off()
        added <- setdiff(names(options()), names(before_options))
        if (length(added)) options(setNames(rep(list(NULL), length(added)), added))
        options(before_options)
        ggplot2::theme_set(before_theme)
      })
      environment <- new.env(parent = globalenv())
      environment$.gg_presence_args <- job
      environment$.gg_presence_read <- read_shared
      environment$.gg_presence_script_path <- renderer
      messages <- capture.output(eval(expressions, envir = environment))
      for (index in seq_along(staged)) messages <- gsub(staged[[index]], outputs[[index]], messages, fixed = TRUE)
      messages
    })
    render_messages <- c(render_messages, messages)
    check()
    gc(verbose = FALSE)
  }
  check()
  if (any(outputs %in% ls(observed, all.names = TRUE))) stop("Plot output replaces a later input")
  if (any(!file.exists(staged)) || any(file.info(staged)$size <= 0)) stop("Plot batch produced an empty output")
  for (path in outputs[file.exists(outputs)]) {
    backup <- tempfile(".gg-presence-backup-", tmpdir = dirname(path))
    backups[[path]] <- backup
    if (!file.copy(path, backup, copy.date = TRUE)) stop("Could not preserve plot: ", path)
  }
  check(verify_bytes = TRUE)
  for (index in seq_along(outputs)) {
    if (!file.rename(staged[[index]], outputs[[index]])) stop("Could not publish plot: ", outputs[[index]])
    installed <- c(installed, outputs[[index]])
  }
  committed <- TRUE
  if (length(render_messages)) cat(paste(render_messages, collapse = "\n"), "\n")
  invisible(outputs)
}
