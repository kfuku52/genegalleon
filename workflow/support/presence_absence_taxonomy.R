read_presence_taxonomy_lineage <- function(path) {
  fields <- c("species", "lineage_index", "taxid", "name", "rank")
  lineage <- read.delim(
    path, colClasses = "character", quote = "", comment.char = "",
    check.names = FALSE, na.strings = NULL, fill = FALSE
  )
  if (!all(fields %in% names(lineage)) || anyDuplicated(names(lineage))) {
    stop("Taxonomy lineage must have species, lineage_index, taxid, name and rank columns: ", path)
  }
  lineage <- lineage[, fields]
  if (anyNA(lineage) || any(vapply(lineage, function(column) any(!nzchar(trimws(column))), logical(1))) ||
      any(!grepl("^[0-9]+$", lineage$lineage_index)) || any(!grepl("^[1-9][0-9]*$", lineage$taxid))) {
    stop("Taxonomy lineage has an empty field or invalid lineage_index/TaxID: ", path)
  }
  # Accept spaces and underscores interchangeably, retaining strain qualifiers.
  original_species <- unique(lineage$species)
  if (anyDuplicated(presence_taxonomy_species_key(original_species))) {
    stop("Taxonomy lineage uses multiple spellings for the same species key: ", path)
  }
  lineage$species <- presence_taxonomy_species_key(lineage$species)
  lineage$lineage_index <- suppressWarnings(as.integer(lineage$lineage_index))
  if (anyNA(lineage$lineage_index) || anyDuplicated(lineage[c("species", "lineage_index")])) {
    stop("Taxonomy lineage has an invalid or repeated lineage_index for a species: ", path)
  }
  lineage[order(lineage$lineage_index, decreasing = TRUE), , drop = FALSE]
}

presence_taxonomy_species_key <- function(species) {
  gsub("[[:space:]]+", "_", trimws(species))
}

presence_taxonomy_ranks <- function(lineage, rank_arg, tree_species, include_upper = TRUE) {
  major_ranks <- c("domain", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus")
  rank_arg <- tolower(trimws(rank_arg))
  if (!identical(rank_arg, "auto")) {
    ranks <- trimws(strsplit(rank_arg, ",", fixed = TRUE)[[1]])
    valid <- c(major_ranks, "species", unique(lineage$rank))
    if (!length(ranks) || endsWith(rank_arg, ",") || anyDuplicated(ranks) || any(!ranks %in% valid)) {
      stop("--taxonomy_ranks must be auto or a comma-separated list of distinct NCBI ranks: ", rank_arg)
    }
    return(ranks)
  }

  tree_lineage <- lineage[lineage$species %in% presence_taxonomy_species_key(tree_species), ]
  resolved_count <- length(unique(tree_lineage$species))
  rank_stats <- vapply(major_ranks, function(rank) {
    rows <- tree_lineage[tree_lineage$rank == rank, ]
    sizes <- table(rows$taxid[!duplicated(rows$species)])
    shared <- sizes[sizes > 1L]
    # Reward repeated groups outside the largest one, so a rank with one
    # dominant category and many singletons has little effect on the choice.
    c(coverage = sum(sizes), shared = sum(shared), score = sum(shared) - max(c(0L, shared)))
  }, numeric(3))
  available <- which(rank_stats["coverage", ] > 0)
  if (!length(available)) return("phylum")
  candidates <- which(rank_stats["coverage", ] >= resolved_count / 2 & rank_stats["score", ] > 0)
  if (!length(candidates)) {
    candidates <- available
    shared <- rank_stats["shared", candidates]
    best <- order(-rank_stats["coverage", candidates], -shared, ifelse(shared > 0, -candidates, candidates))
  } else {
    # Prefer coverage on equal scores, and a coarser rank on an exact tie.
    best <- order(-rank_stats["score", candidates], -rank_stats["coverage", candidates], -rank_stats["shared", candidates])
  }
  selected <- candidates[best[1L]]
  if (!include_upper) return(major_ranks[selected])

  upper <- available[available < selected]
  covered_upper <- upper[rank_stats["coverage", upper] >= resolved_count / 2]
  if (length(covered_upper)) upper <- covered_upper
  major_ranks[c(tail(upper, 1L), selected)]
}

presence_taxonomy_tracks <- function(lineage, group_arg, rank_arg, species, tree_species) {
  group_arg <- trimws(group_arg)
  if (nzchar(group_arg) && !grepl("^[1-9][0-9]*([[:space:]]*,[[:space:]]*[1-9][0-9]*)*$", group_arg)) {
    stop("--taxonomy_groups must be a comma-separated list of positive NCBI TaxIDs")
  }
  groups <- trimws(strsplit(group_arg, ",", fixed = TRUE)[[1]])
  keys <- presence_taxonomy_species_key(species)
  resolved <- keys %in% lineage$species
  track <- function(rows, missing_label) {
    # Rows are ordered deepest first, so the first match wins even if a
    # broader group appears first in the user's argument.
    index <- match(keys, rows$species)
    taxid <- rows$taxid[index]
    label <- rows$name[index]
    label[is.na(index)] <- ifelse(resolved[is.na(index)], missing_label, "Unresolved")
    # Derive the palette from the full lineage input so views that omit taxa
    # still assign the same color to a TaxID.
    ids <- sort(unique(rows$taxid), method = "radix")
    palette <- if (length(ids)) stats::setNames(grDevices::hcl.colors(length(ids), "Pastel 1"), ids) else character()
    fill <- unname(palette[taxid])
    fill[is.na(index)] <- ifelse(resolved[is.na(index)], "#f4f4f4", "#d1d5db")
    data.frame(species = species, key = ifelse(is.na(taxid), label, taxid), label = label, fill = fill)
  }
  ranks <- presence_taxonomy_ranks(lineage, rank_arg, tree_species, include_upper = length(groups) == 0L)
  headings <- paste0(toupper(substr(ranks, 1L, 1L)), substring(ranks, 2L))
  tracks <- stats::setNames(lapply(ranks, function(rank) {
    track(lineage[lineage$rank == rank, ], paste("No", rank))
  }), headings)
  if (!length(groups)) return(tracks)
  c(list("NCBI group" = track(lineage[lineage$taxid %in% groups, ], "Other")), tracks)
}

presence_taxonomy_blocks <- function(track, tip_y, row_half_height) {
  tree_order <- names(sort(tip_y))
  tip_index <- match(track$species, tree_order)
  starts <- which(c(TRUE, diff(tip_index) != 1L | head(track$key, -1L) != tail(track$key, -1L)))
  ends <- c(tail(starts, -1L) - 1L, nrow(track))
  first <- track$species[starts]
  last <- track$species[ends]
  blocks <- track[starts, c("label", "fill"), drop = FALSE]
  blocks$ymin <- unname(tip_y[first] - row_half_height[first])
  blocks$ymax <- unname(tip_y[last] + row_half_height[last])
  blocks
}
