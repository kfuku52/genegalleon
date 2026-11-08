#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
status <- system2("python", c("-m", "nwkit.posterior_diagnostics", "--profile", "chain", vapply(args, shQuote, character(1))))
if (status != 0L) stop("NWKIT posterior diagnostics failed")
