taxonomic_genus_only_placeholders <- c("sp", "spp")
taxonomic_proximity_qualifiers <- c("cf", "aff", "nr")
taxonomic_hybrid_connectors <- c("x", "\u00d7", "hybrid")
taxonomic_rank_aliases <- c(
  "subsp" = "subsp",
  "ssp" = "subsp",
  "subspecies" = "subsp",
  "var" = "var",
  "variety" = "var",
  "forma" = "forma",
  "form" = "forma",
  "f" = "forma",
  "strain" = "strain",
  "substrain" = "substrain",
  "serovar" = "serovar",
  "serotype" = "serotype",
  "serogroup" = "serogroup",
  "pathovar" = "pathovar",
  "pv" = "pathovar",
  "biovar" = "biovar",
  "biotype" = "biotype",
  "chemovar" = "chemovar",
  "morphovar" = "morphovar",
  "cultivar" = "cultivar",
  "cv" = "cultivar",
  "isolate" = "isolate",
  "group" = "group",
  "subgroup" = "subgroup",
  "complex" = "complex",
  "clade" = "clade",
  "lineage" = "lineage",
  "section" = "section",
  "series" = "series",
  "ecotype" = "ecotype",
  "breed" = "breed"
)
taxonomic_display_ranks <- c(
  "subsp" = "subsp.",
  "var" = "var.",
  "forma" = "f.",
  "strain" = "strain",
  "substrain" = "substrain",
  "serovar" = "serovar",
  "serotype" = "serotype",
  "serogroup" = "serogroup",
  "pathovar" = "pathovar",
  "biovar" = "biovar",
  "biotype" = "biotype",
  "chemovar" = "chemovar",
  "morphovar" = "morphovar",
  "cultivar" = "cultivar",
  "isolate" = "isolate",
  "group" = "group",
  "subgroup" = "subgroup",
  "complex" = "complex",
  "clade" = "clade",
  "lineage" = "lineage",
  "section" = "section",
  "series" = "series",
  "ecotype" = "ecotype",
  "breed" = "breed"
)

taxonomic_token_key <- function(token) {
  key <- tolower(sub("[.]$", "", as.character(token)))
  if (key %in% taxonomic_hybrid_connectors) {
    return("x")
  }
  if (key %in% taxonomic_genus_only_placeholders) {
    return("sp")
  }
  if (key %in% names(taxonomic_rank_aliases)) {
    return(unname(taxonomic_rank_aliases[[key]]))
  }
  key
}

is_hybrid_connector_token <- function(token) {
  taxonomic_token_key(token) == "x"
}

is_hybrid_binomial_connector <- function(parts, index = 3L) {
  if (length(parts) < index + 2L) {
    return(FALSE)
  }
  if (!is_hybrid_connector_token(parts[[index]])) {
    return(FALSE)
  }
  next_genus <- trimws(as.character(parts[[index + 1L]]))
  nzchar(next_genus) && grepl("^[A-Z]", next_genus)
}

species_prefix_token_count <- function(parts) {
  parts <- parts[nzchar(parts)]
  if (length(parts) < 2) {
    return(0L)
  }
  second <- taxonomic_token_key(parts[[2]])
  third <- if (length(parts) >= 3) taxonomic_token_key(parts[[3]]) else ""
  if (second == "x") {
    return(if (length(parts) >= 3) 3L else 2L)
  }
  if (second %in% taxonomic_genus_only_placeholders) {
    return(if (length(parts) >= 3) 3L else 2L)
  }
  if (second %in% taxonomic_proximity_qualifiers) {
    return(if (length(parts) >= 3) 3L else 2L)
  }
  if (is_hybrid_binomial_connector(parts, 3L)) {
    return(5L)
  }
  if (third %in% taxonomic_proximity_qualifiers) {
    return(3L)
  }
  if (third %in% unname(taxonomic_rank_aliases)) {
    return(if (length(parts) >= 4) 4L else 3L)
  }
  return(2L)
}

extract_species_label <- function(x) {
  text <- as.character(x)
  if (!nzchar(text)) {
    return("")
  }
  parts <- strsplit(text, "_", fixed = TRUE)[[1]]
  parts <- parts[nzchar(parts)]
  prefix_count <- species_prefix_token_count(parts)
  if (prefix_count == 0) {
    return("")
  }
  selected <- parts[seq_len(prefix_count)]
  is_hybrid_connector <- vapply(selected, is_hybrid_connector_token, logical(1), USE.NAMES = FALSE)
  if (any(is_hybrid_connector)) {
    selected[is_hybrid_connector] <- "x"
  }
  paste(selected, collapse = "_")
}

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
