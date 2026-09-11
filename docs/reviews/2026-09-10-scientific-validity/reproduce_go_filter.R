# Run in a GeneGalleon container; argument 1 is the repository root.
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1]] else "/review"
e <- new.env()
exprs <- parse(file.path(root, "workflow/support/cafe_go_enrichment.r"))
# Load only function definitions, without executing the helper's CLI.
for (x in exprs) {
  if (is.call(x) && as.character(x[[1]]) %in% c("<-", "=") &&
      is.call(x[[3]]) && identical(x[[3]][[1]], as.name("function"))) eval(x, e)
}
events <- data.frame(FamilyID = paste0("F", 1:100),
                     is_target = c(rep(TRUE, 10), rep(FALSE, 90)))
rows <- events[c(1:4, 11:16), ]
rows$go_ids <- "GO:target"
rows$go_aspects <- "BP"
rows$go_terms <- "target"
common <- events
common$go_ids <- "GO:common"
common$go_aspects <- "BP"
common$go_terms <- "common"
rows <- rbind(rows, common)
for (i in 1:8) {
  more <- events[17 + i, ]
  more$go_ids <- paste0("GO:other", i)
  more$go_aspects <- "BP"
  more$go_terms <- paste0("other", i)
  rows <- rbind(rows, more)
}
result <- e$summarise_go_enrichment(rows, 10, 90, 0.05)
print(result$all[, c("go_ids", "p_value", "p_value_adjusted")])
cat("Number of eligible GO terms:", length(unique(rows$go_ids)), "\n")
cat("Number retained for BH:", nrow(result$all), "\n")
cat("BH with all 10 tested terms:",
    p.adjust(c(result$all$p_value, rep(1, 8)), method = "BH")[1], "\n")
cat("Number declared significant by current implementation:",
    nrow(result$significant), "\n")
