# Dataset identity comes from BUSCO's own header, never from the gene count.
gg_busco_lineage_from_header <- function(lines) {
    matches <- grep('lineage dataset is:[[:blank:]]*[^[:space:]]+', lines, value=TRUE, ignore.case=TRUE)
    values <- unique(sub('.*lineage dataset is:[[:blank:]]*([^[:space:]]+).*', '\\1',
                         matches, ignore.case=TRUE))
    if (length(values) > 1L) stop('Conflicting BUSCO lineage datasets in one result.')
    if (length(values) == 0L) return(NA_character_)
    values[[1L]]
}

gg_busco_axis_label <- function(lineages) {
    lineages <- trimws(lineages)
    known <- unique(lineages[!is.na(lineages) & nzchar(lineages)])
    if (length(known) > 1L) stop('Mixed BUSCO lineage datasets cannot share a completeness axis.')
    if (length(known) == 1L && all(!is.na(lineages) & nzchar(lineages))) {
        return(paste0('Number of BUSCO genes\n(', known[[1L]], ')'))
    }
    'Number of BUSCO genes'
}
