# CAFE labels combine a species/branch key with an ancestral copy-number suffix.
# The tree fortifier passes whole label vectors, including internal branch IDs.
norm_key <- function(lbl) {
  # The table key includes CAFE's <branch> suffix. Taxonomic prefix parsing
  # would truncate longer species keys before that suffix.
  sub("_[0-9]+$", "", gsub("\\*", "", lbl))
}

get_species <- function(lbl) {
  vapply(norm_key(lbl), function(label) {
    x <- sub("<.*$", "", label)
    if (is.na(x)) return(NA_character_)
    if (grepl("^[A-Za-z][A-Za-z0-9]*_", x)) {
      species_label <- extract_species_label(x)
      if (identical(species_label, x)) return(scientific_name_from_label(species_label))
    }
    gsub("_", " ", x)
  }, character(1), USE.NAMES = FALSE)
}
