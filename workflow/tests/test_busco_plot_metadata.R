args <- commandArgs(trailingOnly=FALSE)
script <- sub('^--file=', '', args[grepl('^--file=', args)][[1L]])
script <- normalizePath(script)
source(file.path(dirname(script), '..', 'support', 'busco_plot_metadata.r'))

lineage <- gg_busco_lineage_from_header(c('# BUSCO version is: 5.8.2',
  '# The lineage dataset is: embryophyta_odb12 (Creation date: 2025-04-11, number of BUSCOs: 2026)'))
stopifnot(identical(lineage, 'embryophyta_odb12'))
stopifnot(identical(gg_busco_axis_label(rep(lineage, 3)), 'Number of BUSCO genes\n(embryophyta_odb12)'))
stopifnot(is.na(gg_busco_lineage_from_header('# Busco id\tStatus')))
stopifnot(is.na(gg_busco_lineage_from_header(c('# lineage dataset is:  ', '# Busco id\tStatus'))))
stopifnot(identical(gg_busco_axis_label(c(lineage, '  ')), 'Number of BUSCO genes'))
stopifnot(identical(gg_busco_axis_label(c(lineage, NA_character_)), 'Number of BUSCO genes'))
mixed <- try(gg_busco_axis_label(c('embryophyta_odb10', 'embryophyta_odb12')), silent=TRUE)
stopifnot(inherits(mixed, 'try-error'))
conflict <- try(gg_busco_lineage_from_header(c('# lineage dataset is: a_odb12', '# lineage dataset is: b_odb12')), silent=TRUE)
stopifnot(inherits(conflict, 'try-error'))

# Exercise the real producer, including propagation into the saved summary
# and the two-line label in its SVG, without writing curated workflow data.
check_summary_producer <- function(with_lineage=TRUE) {
    fixture <- tempfile('busco-plot-metadata-')
    dir.create(fixture)
    on.exit(unlink(fixture, recursive=TRUE), add=TRUE)
    input <- file.path(fixture, 'full')
    dir.create(input)
    lines <- c('# The lineage dataset is: embryophyta_odb12 (number of BUSCOs: 4)',
               '# Busco id\tStatus\tSequence', 'id1\tComplete\tseq1',
               'id2\tDuplicated\tseq2', 'id3\tFragmented\tseq3', 'id4\tMissing\t')
    if (!with_lineage) lines[[1L]] <- '# The lineage dataset is:'
    for (species in c('A_a', 'B_b')) {
        writeLines(lines, file.path(input, paste0(species, '.busco.full.tsv')))
    }
    # Native BUSCO exports keep reusable single-copy records beside the tables.
    # A directory must never enter the table reader, even with a table suffix.
    dir.create(file.path(input, 'single_copy'))
    writeLines('{}', file.path(input, 'single_copy', 'A_a.json.gz'))
    dir.create(file.path(input, 'C_c.busco.full.tsv'))
    previous <- getwd()
    on.exit(setwd(previous), add=TRUE)
    setwd(fixture)
    output <- system2(file.path(R.home('bin'), 'Rscript'),
        c(shQuote(file.path(dirname(script), '..', 'support', 'annotation_summary.r')),
          paste0('--dir_species_cds_busco=', shQuote(input)), '--min_og_species=auto'),
        stdout=TRUE, stderr=TRUE)
    if (!is.null(attr(output, 'status'))) stop(paste(output, collapse='\n'))
    stopifnot(any(grepl('Number of BUSCO full tables: 2', output, fixed=TRUE)))
    summary <- read.delim('annotation_summary.tsv', check.names=FALSE)
    stopifnot(nrow(summary) == 2L, all(summary$busco_cds_total == 4L), all(summary$busco_cds_single == 1L))
    if (with_lineage) stopifnot(all(summary$busco_cds_lineage == 'embryophyta_odb12'))
    else stopifnot(all(is.na(summary$busco_cds_lineage)))
    svg <- readLines('busco_cds.svg')
    stopifnot(any(grepl('Number of BUSCO genes', svg, fixed=TRUE)),
              any(grepl('(embryophyta_odb12)', svg, fixed=TRUE)) == with_lineage)
}
check_summary_producer()
check_summary_producer(FALSE)
cat('BUSCO plot metadata checks passed.\n')
