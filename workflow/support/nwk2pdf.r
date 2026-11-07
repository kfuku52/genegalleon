args <- commandArgs(trailingOnly=TRUE)
status <- system2("python", c("-m", "nwkit.r_render", vapply(args, shQuote, character(1))))
if (status != 0L) stop("NWKIT basic tree rendering failed")
