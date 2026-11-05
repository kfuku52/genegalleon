# Keep the public plotting helpers on the common species-label contract.
gg_tree_label_sources <- Filter(function(path) is.character(path) && length(path) == 1L && nzchar(path),
                              lapply(sys.frames(), function(frame) frame$ofile))
gg_tree_label_source <- tail(gg_tree_label_sources, 1L)[[1L]]
source(file.path(dirname(normalizePath(gg_tree_label_source, mustWork = TRUE)),
                 "species_label_utils.r"), local = TRUE)
rm(gg_tree_label_source, gg_tree_label_sources)

taxonomic_genus_only_placeholders <- GG_TAXONOMIC_GENUS_ONLY_PLACEHOLDERS
taxonomic_proximity_qualifiers <- GG_TAXONOMIC_PROXIMITY_QUALIFIERS
taxonomic_hybrid_connectors <- GG_TAXONOMIC_HYBRID_CONNECTORS
taxonomic_rank_aliases <- GG_TAXONOMIC_RANK_ALIASES
taxonomic_display_ranks <- GG_TAXONOMIC_DISPLAY_RANKS
taxonomic_token_key <- gg_taxonomic_token_key
is_hybrid_connector_token <- gg_is_hybrid_connector_token
is_hybrid_binomial_connector <- gg_is_hybrid_binomial_connector
species_prefix_token_count <- function(parts) as.integer(gg_species_prefix_token_count(parts))
extract_species_label <- function(x) gg_species_label_from_value_one(x)

strip_species_label <- function(x) {
  text <- as.character(x)
  species_label <- extract_species_label(text)
  prefix <- paste0(species_label, "_")
  if (nzchar(species_label) && startsWith(text, prefix)) {
    return(substr(text, nchar(prefix) + 1, nchar(text)))
  }
  text
}

scientific_name_from_label <- function(x) {
  species_label <- extract_species_label(x)
  if (!nzchar(species_label)) {
    species_label <- as.character(x)
  }
  parts <- strsplit(species_label, "_", fixed = TRUE)[[1]]
  parts <- parts[nzchar(parts)]
  if (length(parts) >= 5 && is_hybrid_connector_token(parts[[3]])) {
    return(sprintf("%s %s x %s %s", parts[[1]], parts[[2]], parts[[4]], parts[[5]]))
  }
  if (length(parts) >= 3 && taxonomic_token_key(parts[[2]]) %in% taxonomic_proximity_qualifiers) {
    return(sprintf("%s %s. %s", parts[[1]], taxonomic_token_key(parts[[2]]), parts[[3]]))
  }
  if (length(parts) >= 3 && taxonomic_token_key(parts[[3]]) %in% taxonomic_proximity_qualifiers) {
    return(sprintf("%s %s. %s", parts[[1]], taxonomic_token_key(parts[[3]]), parts[[2]]))
  }
  if (length(parts) >= 3 && taxonomic_token_key(parts[[2]]) == "sp") {
    return(sprintf("%s sp. %s", parts[[1]], parts[[3]]))
  }
  if (length(parts) >= 4 && taxonomic_token_key(parts[[3]]) %in% unname(taxonomic_rank_aliases)) {
    rank <- taxonomic_token_key(parts[[3]])
    return(sprintf("%s %s %s %s", parts[[1]], parts[[2]], taxonomic_display_ranks[[rank]], parts[[4]]))
  }
  gsub("_", " ", species_label)
}
