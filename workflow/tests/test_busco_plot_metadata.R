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

# Compare the production annotation reader with the original full-table counts,
# including missing columns, comments, quoted text and numeric/NA inference.
check_annotation_counts <- function() {
    fixture <- tempfile('annotation-column-counts-')
    dir.create(fixture)
    on.exit(unlink(fixture, recursive=TRUE), add=TRUE)
    env <- new.env(parent=globalenv())
    source(file.path(dirname(script), '..', 'support', 'species_label_utils.r'), local=env)
    env$path_exists_or_null <- function(path) dir.exists(path)
    expressions <- parse(file.path(dirname(script), '..', 'support', 'annotation_summary.r'))
    assignment <- Filter(function(expr) is.call(expr) && length(expr)>=3 &&
        identical(expr[[2]], as.name('generate_annotation_summary')), expressions)[[1]]
    eval(assignment, envir=env)
    fn <- env$generate_annotation_summary
    expressions <- as.list(body(fn))
    cut <- which(vapply(expressions, function(expr) is.call(expr) &&
        identical(expr[[1]], as.name('write.table')), logical(1)))[1]
    body(fn) <- as.call(c(expressions[seq_len(cut-1)], list(quote(df))))
    cases <- list(
        c('gene_id\torthogroup\tsprot_best\tbusco_id\textra',
          'g1\t001\tNA\t"id1"\tunused', 'g2\t\tP1\t2\t# comment'),
        c('gene_id\textra', 'g1\tunused', 'g2\tNA'),
        c('extra', 'unused', 'NA'),
        c('gene_id\tbusco_id', 'g1\t1', 'g2\tNA'))
    columns <- c('gene_id', 'orthogroup', 'sprot_best', 'cdskit_localize_predicted_class', 'busco_id')
    for (lines in cases) {
        path <- file.path(fixture, 'Species_a.tsv')
        writeLines(lines, path)
        original <- read.table(path, sep='\t', header=TRUE, comment.char='#', fill=TRUE, quote='')
        expected <- apply(original, 2, function(x) sum(!is.na(x) & x!=''))[columns]
        expected[is.na(expected)] <- 0
        invisible(capture.output(actual <- fn(fixture, 'unused')))
        stopifnot(identical(as.integer(actual[1, -1]), as.integer(expected)))
    }
}
check_annotation_counts()
cat('Selective annotation count checks passed.\n')
