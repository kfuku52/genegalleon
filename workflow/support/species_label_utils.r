GG_TAXONOMIC_PROXIMITY_QUALIFIERS <- c("cf", "aff", "nr")
GG_TAXONOMIC_GENUS_ONLY_PLACEHOLDERS <- c("sp", "spp")
GG_TAXONOMIC_INFRASPECIFIC_RANKS <- c(
    "subsp",
    "var",
    "forma",
    "strain",
    "substrain",
    "serovar",
    "serotype",
    "serogroup",
    "pathovar",
    "biovar",
    "biotype",
    "chemovar",
    "morphovar",
    "cultivar",
    "isolate",
    "group",
    "subgroup",
    "complex",
    "clade",
    "lineage",
    "section",
    "series",
    "ecotype",
    "breed"
)
GG_TAXONOMIC_RANK_ALIASES <- c(
    subsp = "subsp",
    ssp = "subsp",
    subspecies = "subsp",
    var = "var",
    variety = "var",
    forma = "forma",
    form = "forma",
    f = "forma",
    strain = "strain",
    substrain = "substrain",
    serovar = "serovar",
    serotype = "serotype",
    serogroup = "serogroup",
    pathovar = "pathovar",
    pv = "pathovar",
    biovar = "biovar",
    biotype = "biotype",
    chemovar = "chemovar",
    morphovar = "morphovar",
    cultivar = "cultivar",
    cv = "cultivar",
    isolate = "isolate",
    group = "group",
    subgroup = "subgroup",
    complex = "complex",
    clade = "clade",
    lineage = "lineage",
    section = "section",
    series = "series",
    ecotype = "ecotype",
    breed = "breed"
)
GG_TAXONOMIC_DISPLAY_RANKS <- c(
    subsp = "subsp.",
    var = "var.",
    forma = "f.",
    strain = "strain",
    substrain = "substrain",
    serovar = "serovar",
    serotype = "serotype",
    serogroup = "serogroup",
    pathovar = "pathovar",
    biovar = "biovar",
    biotype = "biotype",
    chemovar = "chemovar",
    morphovar = "morphovar",
    cultivar = "cultivar",
    isolate = "isolate",
    group = "group",
    subgroup = "subgroup",
    complex = "complex",
    clade = "clade",
    lineage = "lineage",
    section = "section",
    series = "series",
    ecotype = "ecotype",
    breed = "breed"
)
GG_SPECIES_LABEL_TERMINAL_SUFFIXES <- c(
    ".fasta.busco.full.tsv",
    ".fa.busco.full.tsv",
    ".faa.busco.full.tsv",
    ".fna.busco.full.tsv",
    ".ffn.busco.full.tsv",
    ".fasta.busco.short.txt",
    ".fa.busco.short.txt",
    ".faa.busco.short.txt",
    ".fna.busco.short.txt",
    ".ffn.busco.short.txt",
    ".busco.full.tsv",
    "_busco.full.tsv",
    ".busco.short.txt",
    "_busco.short.txt",
    ".busco.full",
    "_busco.full",
    ".busco.short",
    "_busco.short",
    ".busco",
    "_busco",
    "_longestCDS.transcript.fa.gz",
    ".longestCDS.transcript.fa.gz",
    "_longestCDS_contamination_removal.fa.gz",
    ".longestCDS_contamination_removal.fa.gz",
    "_longestCDS.fx2tab_cds.tsv",
    ".longestCDS.fx2tab_cds.tsv",
    "_longestCDS.mmseqs2taxonomy.tsv",
    ".longestCDS.mmseqs2taxonomy.tsv",
    "_longestCDS.fa.gz",
    ".longestCDS.fa.gz",
    "_isoform.fa.gz",
    ".isoform.fa.gz",
    "_fx2tab_cds.tsv",
    ".fx2tab_cds.tsv",
    "_fx2tab_genome.tsv",
    ".fx2tab_genome.tsv",
    "_annotation.tsv",
    ".annotation.tsv",
    "_raw_counts.tsv",
    "_percentages.tsv",
    "_uniprot.tsv",
    ".uniprot.tsv",
    "_assembly_stat.tsv",
    ".fastq.gz",
    ".fq.gz",
    ".fasta.gz",
    ".fa.gz",
    ".faa.gz",
    ".fna.gz",
    ".ffn.gz",
    ".gff3.gz",
    ".gff.gz",
    ".gtf.gz",
    ".tsv.gz",
    ".txt.gz",
    ".csv.gz",
    ".fastq",
    ".fq",
    ".fasta",
    ".fa",
    ".faa",
    ".fna",
    ".ffn",
    ".gff3",
    ".gff",
    ".gtf",
    ".tsv",
    ".txt",
    ".csv"
)

gg_species_strip_terminal_suffixes_one <- function(value) {
    stem <- basename(as.character(value))
    stem_lower <- tolower(stem)
    for (suffix in GG_SPECIES_LABEL_TERMINAL_SUFFIXES) {
        suffix_lower <- tolower(suffix)
        if (endsWith(stem_lower, suffix_lower)) {
            return(substr(stem, 1, nchar(stem) - nchar(suffix)))
        }
    }
    return(stem)
}

gg_species_strip_terminal_suffixes <- function(values) {
    vapply(values, gg_species_strip_terminal_suffixes_one, character(1), USE.NAMES = FALSE)
}

gg_taxonomic_token_key <- function(token) {
    lowered <- tolower(sub("[.]$", "", trimws(as.character(token))))
    if (lowered %in% GG_TAXONOMIC_GENUS_ONLY_PLACEHOLDERS) {
        return("sp")
    }
    if (lowered %in% GG_TAXONOMIC_PROXIMITY_QUALIFIERS) {
        return(lowered)
    }
    if (lowered %in% names(GG_TAXONOMIC_RANK_ALIASES)) {
        return(unname(GG_TAXONOMIC_RANK_ALIASES[[lowered]]))
    }
    return(lowered)
}

gg_species_prefix_token_count <- function(parts) {
    normalized <- parts[nzchar(parts)]
    if (length(normalized) < 2) {
        return(length(normalized))
    }
    second <- gg_taxonomic_token_key(normalized[[2]])
    third <- if (length(normalized) >= 3) gg_taxonomic_token_key(normalized[[3]]) else ""
    if (second == "sp") {
        return(if (length(normalized) >= 3) 3 else 2)
    }
    if (second %in% GG_TAXONOMIC_PROXIMITY_QUALIFIERS) {
        return(if (length(normalized) >= 3) 3 else 2)
    }
    if (third %in% GG_TAXONOMIC_PROXIMITY_QUALIFIERS) {
        return(3)
    }
    if (third %in% GG_TAXONOMIC_INFRASPECIFIC_RANKS) {
        return(if (length(normalized) >= 4) 4 else 3)
    }
    return(2)
}

gg_species_label_from_value_one <- function(value, strip_extension = FALSE) {
    text <- as.character(value)
    if (is.na(text) || !nzchar(text)) {
        return("")
    }
    text <- basename(text)
    if (strip_extension) {
        text <- gg_species_strip_terminal_suffixes_one(text)
    }
    text <- gsub("[[:space:]]+", "_", trimws(text), perl = TRUE)
    parts <- strsplit(text, "_", fixed = TRUE)[[1]]
    parts <- parts[nzchar(parts)]
    count <- gg_species_prefix_token_count(parts)
    if (count == 0) {
        return("")
    }
    return(paste(parts[seq_len(count)], collapse = "_"))
}

gg_species_label_from_text <- function(values) {
    vapply(values, gg_species_label_from_value_one, character(1), strip_extension = FALSE, USE.NAMES = FALSE)
}

gg_species_label_from_filename <- function(values) {
    vapply(values, gg_species_label_from_value_one, character(1), strip_extension = TRUE, USE.NAMES = FALSE)
}

gg_species_display_rank <- function(rank) {
    if (rank == "sp") {
        return("sp.")
    }
    if (rank %in% GG_TAXONOMIC_PROXIMITY_QUALIFIERS) {
        return(paste0(rank, "."))
    }
    if (rank %in% names(GG_TAXONOMIC_DISPLAY_RANKS)) {
        return(unname(GG_TAXONOMIC_DISPLAY_RANKS[[rank]]))
    }
    return(rank)
}

gg_species_display_from_key_one <- function(value) {
    species_label <- as.character(value)
    if (is.na(species_label) || !nzchar(species_label)) {
        return("")
    }
    parts <- strsplit(species_label, "_", fixed = TRUE)[[1]]
    parts <- parts[nzchar(parts)]
    if (length(parts) >= 3 && gg_taxonomic_token_key(parts[[2]]) %in% GG_TAXONOMIC_PROXIMITY_QUALIFIERS) {
        return(sprintf("%s %s %s", parts[[1]], gg_species_display_rank(gg_taxonomic_token_key(parts[[2]])), parts[[3]]))
    }
    if (length(parts) >= 3 && gg_taxonomic_token_key(parts[[3]]) %in% GG_TAXONOMIC_PROXIMITY_QUALIFIERS) {
        return(sprintf("%s %s %s", parts[[1]], gg_species_display_rank(gg_taxonomic_token_key(parts[[3]])), parts[[2]]))
    }
    if (length(parts) >= 3 && gg_taxonomic_token_key(parts[[2]]) == "sp") {
        return(sprintf("%s sp. %s", parts[[1]], parts[[3]]))
    }
    if (length(parts) >= 2 && gg_taxonomic_token_key(parts[[2]]) == "sp") {
        return(sprintf("%s sp.", parts[[1]]))
    }
    if (length(parts) >= 4 && gg_taxonomic_token_key(parts[[3]]) %in% GG_TAXONOMIC_INFRASPECIFIC_RANKS) {
        rank <- gg_taxonomic_token_key(parts[[3]])
        return(sprintf("%s %s %s %s", parts[[1]], parts[[2]], gg_species_display_rank(rank), parts[[4]]))
    }
    return(gsub("_", " ", species_label, fixed = TRUE))
}

gg_species_display_from_key <- function(values) {
    vapply(values, gg_species_display_from_key_one, character(1), USE.NAMES = FALSE)
}

gg_species_display_from_label <- function(values, strip_extension = FALSE) {
    keys <- vapply(values, gg_species_label_from_value_one, character(1), strip_extension = strip_extension, USE.NAMES = FALSE)
    displays <- gg_species_display_from_key(keys)
    fallback <- !nzchar(displays)
    if (any(fallback)) {
        values_for_display <- as.character(values[fallback])
        if (strip_extension) {
            values_for_display <- gg_species_strip_terminal_suffixes(values_for_display)
        }
        displays[fallback] <- gsub("_", " ", values_for_display, fixed = TRUE)
    }
    return(displays)
}

gg_stop_on_duplicate_species_keys <- function(keys, labels, context) {
    if (length(keys) == 0) {
        return(invisible(NULL))
    }
    labels <- as.character(labels)
    keys <- as.character(keys)
    valid <- !is.na(keys) & nzchar(keys)
    duplicate_keys <- unique(keys[valid][duplicated(keys[valid])])
    if (length(duplicate_keys) == 0) {
        return(invisible(NULL))
    }
    lines <- vapply(
        duplicate_keys,
        function(key) {
            sprintf("  %s: %s", key, paste(labels[which(keys == key)], collapse = ", "))
        },
        character(1),
        USE.NAMES = FALSE
    )
    stop(
        sprintf(
            "Duplicate species labels in %s after parsing:\n%s",
            context,
            paste(lines, collapse = "\n")
        ),
        call. = FALSE
    )
}
