# Benchmark real preprocessing blocks; validate plotting separately.
args <- commandArgs(TRUE)
support <- args[[1]]
case <- args[[2]]
fixture <- args[[3]]
scratch <- args[[4]]
env <- new.env(parent=globalenv())
source(file.path(support, 'species_label_utils.r'), local=env)
source(file.path(support, 'busco_plot_metadata.r'), local=env)
env$path_exists_or_null <- function(path) !is.null(path) && length(path) && file.exists(path)
name <- if (case == 'r-annotation') 'generate_annotation_summary' else 'generate_busco_summary'
assignment <- Filter(function(expr) is.call(expr) && length(expr) >= 3 &&
    identical(expr[[2]], as.name(name)), parse(file.path(support, 'annotation_summary.r')))[[1]]
eval(assignment, envir=env)
fn <- env[[name]]
expressions <- as.list(body(fn))
if (case == 'r-annotation') {
    cutoff <- which(vapply(expressions, function(expr) is.call(expr) &&
        identical(expr[[1]], as.name('write.table')), logical(1)))[1]
    input <- file.path(fixture, 'annotation')
    counts <- list(species=3, genes_per_species=30000, columns=80, counted_columns=5)
} else {
    cutoff <- which(vapply(expressions, function(expr) is.call(expr) &&
        identical(expr[[1]], as.name('if')) && grepl('nrow(df)',
            paste(deparse(expr[[2]]), collapse=''), fixed=TRUE), logical(1)))[1]
    input <- file.path(scratch, 'busco')
    dir.create(input)
    for (i in seq_len(400)) {
        writeLines(c('# Busco id\tStatus', paste0('id', seq_len(1000), '\tComplete')),
                   file.path(input, paste0('Species_', sprintf('%03d', i), '.tsv')))
    }
    counts <- list(species=400, records_per_species=1000)
}
body(fn) <- as.call(c(expressions[seq_len(cutoff - 1)], list(quote(df))))
invisible(gc())
started <- proc.time()[['elapsed']]
if (case == 'r-annotation') {
    invisible(capture.output(result <- fn(input, file.path(scratch, 'result'))))
} else {
    invisible(capture.output(result <- fn(input, 'unused', data.frame())))
}
elapsed <- proc.time()[['elapsed']] - started
path <- file.path(scratch, 'counts.tsv')
write.table(result, path, sep='\t', quote=FALSE, row.names=FALSE)
cat(jsonlite::toJSON(list(seconds=elapsed, output_sha256=digest::digest(file=path, algo='sha256'),
    counts=counts), auto_unbox=TRUE), '\n')
