.gg_topology_helper_dir <- local({
  paths <- Filter(function(x) is.character(x) && length(x) == 1L && nzchar(x),
                  lapply(sys.frames(), function(frame) frame$ofile))
  dirname(normalizePath(tail(paths, 1L)[[1L]], mustWork = TRUE))
})

# Shared topology-based support transfer; callers choose the allowed tip scope.
gg_descendant_tips <- function(tree, node) {
  n_tip <- length(tree[["tip.label"]])
  child_map <- split(tree[["edge"]][, 2], tree[["edge"]][, 1])
  tips <- integer()
  pending <- node
  while (length(pending)) {
    current <- pending[[length(pending)]]
    pending <- pending[-length(pending)]
    children <- child_map[[as.character(current)]]
    tips <- c(tips, children[children <= n_tip])
    pending <- c(pending, rev(children[children > n_tip]))
  }
  tips
}

gg_descendant_tip_labels <- function(tree, node, tip_set = NULL) {
  labels <- tree[["tip.label"]][gg_descendant_tips(tree, node)]
  if (!is.null(tip_set)) {
    labels <- labels[labels %in% tip_set]
  }
  labels
}

gg_node_clade_key <- local({
  cache_tree <- NULL
  cache_tips <- NULL
  cache_keys <- NULL
  function(tree, node, tip_set = NULL) {
    if (!identical(cache_tree, tree) || !identical(cache_tips, tip_set)) {
      request <- tempfile(fileext = ".json")
      response <- tempfile(fileext = ".json")
      on.exit(unlink(c(request, response)), add = TRUE)
      data <- list(operation = jsonlite::unbox("clade_keys"), tree = list(
        tip.label = tree[["tip.label"]], Nnode = jsonlite::unbox(as.integer(tree[["Nnode"]])),
        edge = tree[["edge"]], node.label = tree[["node.label"]]))
      if (!is.null(tip_set)) data$tip_set <- tip_set
      jsonlite::write_json(data, request)
      status <- system2("python", c(shQuote(file.path(.gg_topology_helper_dir, "transfer_tree_labels.py")),
                                   shQuote(request), shQuote(response)))
      if (status != 0L) stop("NWKIT clade indexing failed")
      cache_keys <<- jsonlite::read_json(response, simplifyVector = FALSE)
      cache_tree <<- tree
      cache_tips <<- tip_set
    }
    value <- cache_keys[[node]]
    if (is.null(value)) NA_character_ else as.character(value)
  }
})

gg_transfer_internal_labels <- function(target_tree, source_tree, allow_subset = FALSE) {
  target_tips <- target_tree[["tip.label"]]
  compatible <- if (allow_subset) {
    length(setdiff(target_tips, source_tree[["tip.label"]])) == 0
  } else {
    setequal(target_tips, source_tree[["tip.label"]])
  }
  if (!compatible) {
    warning(if (allow_subset) {
      "Cannot transfer support values: support tree does not contain all target tips."
    } else {
      "Cannot transfer support values: target and source trees have different tip labels."
    })
    return(target_tree)
  }

  source_labels <- source_tree[["node.label"]]
  if (is.null(source_labels) || length(source_labels) == 0 || target_tree[["Nnode"]] < 1L) {
    return(target_tree)
  }
  numbered_tree <- function(tree) list(
    tip.label = tree[["tip.label"]], Nnode = jsonlite::unbox(as.integer(tree[["Nnode"]])),
    edge = tree[["edge"]], node.label = tree[["node.label"]]
  )
  request <- tempfile("gg-support-", fileext = ".json")
  response <- tempfile("gg-support-", fileext = ".json")
  diagnostics <- tempfile("gg-support-")
  on.exit(unlink(c(request, response, diagnostics)), add = TRUE)
  jsonlite::write_json(list(target = numbered_tree(target_tree), source = numbered_tree(source_tree),
                           allow_subset = jsonlite::unbox(allow_subset)), request)
  status <- system2("python", c(shQuote(file.path(.gg_topology_helper_dir, "transfer_tree_labels.py")),
                                shQuote(request), shQuote(response)), stderr = diagnostics)
  if (status != 0L) stop("NWKIT support-label transfer failed: ", paste(readLines(diagnostics, warn = FALSE), collapse = "\n"))
  target_tree[["node.label"]] <- as.character(jsonlite::read_json(response, simplifyVector = TRUE))
  target_tree
}
