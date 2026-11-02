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
    # Annotation tables are generated TSVs: hashes inside their identifiers
    # and annotation values are data, not end-of-line comments.
    writeLines(c('gene_id\tbusco_id', 'g#1\tb#1', 'g2\tb2'), path)
    invisible(capture.output(actual <- fn(fixture, 'unused')))
    stopifnot(identical(as.integer(actual[1, -1]), c(2L, 0L, 0L, 0L, 2L)))
}
check_annotation_counts()
cat('Selective annotation count checks passed.\n')

# Load the production trait transformation without starting its plotting CLI.
# Literal expectations protect value coercion, row names and trait-column order.
check_trait_tidy <- function() {
    expressions <- parse(file.path(dirname(script), '..', 'support', 'annotation_summary.r'))
    branch <- Filter(function(expr) is.call(expr) &&
        identical(expr[[1]], as.name('if')) && 'trait_tidy' %in% all.names(expr), expressions)[[1]]
    while (is.call(branch) && identical(branch[[1]], as.name('if'))) branch <- branch[[4L]]
    expressions <- as.list(branch)
    assigns <- function(expr, name) is.call(expr) && length(expr) >= 3L &&
        as.character(expr[[1]]) %in% c('=', '<-') && identical(expr[[2]], as.name(name))
    begin <- which(vapply(expressions, assigns, logical(1), name='trait_tidy'))[[1L]]
    end <- which(vapply(expressions, assigns, logical(1), name='g'))[[1L]] - 1L
    transform <- function(trait) NULL
    body(transform) <- as.call(c(list(as.name('{')), expressions[begin:end], list(quote(trait_tidy))))

    ordinary <- data.frame(label=c('A_a','B_b'), first=c(1,NA), second=c(2,3))
    expected <- data.frame(label=c('A_a','B_b','A_a','B_b'), value=c(1,NA,2,3),
                           trait=c('first','first','second','second'))
    stopifnot(identical(transform(ordinary), expected))
    one <- ordinary[,c('label','second')]
    stopifnot(identical(transform(one), data.frame(label=c('A_a','B_b'), value=c(2,3), trait=c('second','second'))))
    mixed <- data.frame(label=c('A_a','B_b'), first=c(1,NA), second=c('up','down'))
    expected$value <- c('1',NA,'up','down')
    stopifnot(identical(transform(mixed), expected))
    factors <- data.frame(label=c('A_a','B_b'), first=ordered(c('a','b')), second=ordered(c('b','c')))
    expected$value <- ordered(c('a','b','b','c'), levels=c('a','b','c'))
    stopifnot(identical(transform(factors), expected))
    dates <- data.frame(label=c('A_a','B_b'), first=as.Date(c('2020-01-01',NA)),
                        second=as.Date(c('2021-01-01','2022-01-01')))
    expected$value <- as.Date(c('2020-01-01',NA,'2021-01-01','2022-01-01'))
    stopifnot(identical(transform(dates), expected))
    empty <- ordinary[FALSE,]
    empty_error <- tryCatch(transform(empty), error=identity)
    stopifnot(inherits(empty_error,'error'),
              grepl('replacement has 1 row, data has 0',conditionMessage(empty_error),fixed=TRUE))
    no_columns <- tryCatch(transform(data.frame()), error=identity)
    stopifnot(inherits(no_columns,'error'), grepl('undefined columns selected',conditionMessage(no_columns),fixed=TRUE))
    rownames(ordinary) <- c('first_row','second_row')
    expected$value <- c(1,NA,2,3)
    rownames(expected) <- c('first_row','second_row','first_row1','second_row1')
    stopifnot(identical(transform(ordinary), expected))
    duplicate <- data.frame(label=c('A_a','B_b'), first=1:2, second=3:4)
    names(duplicate) <- c('label','same','same')
    expected <- data.frame(label=c('A_a','B_b','A_a','B_b'),value=c(1L,2L,1L,2L),trait=rep('same',4))
    stopifnot(identical(transform(duplicate), expected))
    missing <- tryCatch(transform(data.frame(label=c('A_a','B_b'))), error=identity)
    stopifnot(inherits(missing,'error'), grepl('undefined columns selected',conditionMessage(missing),fixed=TRUE))
}
check_trait_tidy()
cat('Trait tidy value and type checks passed.\n')
