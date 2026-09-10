# Run from the repository root in a GeneGalleon container.
e <- new.env()
for (x in parse('workflow/support/cafe_go_enrichment.r')) {
  if (is.call(x) && as.character(x[[1]]) %in% c('<-', '=') &&
      is.call(x[[3]]) && identical(x[[3]][[1]], as.name('function'))) eval(x, e)
}
expect_error <- function(expr) stopifnot(inherits(tryCatch(force(expr), error = identity), 'error'))

# Native model P values are adjusted across all families before selecting sign.
# An accelerated lambda is turnover, not necessarily an increase in copies.
families <- data.frame(FamilyID = c('gain', 'loss', 'slow', 'flat', 'background'),
  status = 'tested', p_value = c(.001, .002, .001, .001, 1),
  target_change = c(10, -10, 2, 0, 1),
  lambda_target = c(.2, .2, .01, .2, .1), lambda_background = .1)
inc <- e$select_cafe_families(families, .05, 'increase')
dec <- e$select_cafe_families(families, .05, 'decrease')
both <- e$select_cafe_families(families, .05, 'both')
stopifnot(identical(inc$FamilyID[inc$selected], 'gain'),
          identical(dec$FamilyID[dec$selected], 'loss'),
          setequal(both$FamilyID[both$selected], c('gain', 'loss')),
          identical(inc$p_value_adjusted, dec$p_value_adjusted),
          identical(both$p_value_adjusted, p.adjust(families$p_value, 'BH')))
bad <- families; bad$status[1] <- 'failed'
expect_error(e$select_cafe_families(bad, .05, 'both'))
bad <- families; bad$p_value[1] <- NA
expect_error(e$select_cafe_families(bad, .05, 'both'))
bad <- families; bad$FamilyID[1] <- bad$FamilyID[2]
expect_error(e$select_cafe_families(bad, .05, 'both'))
expect_error(e$select_cafe_families(families, NA_real_, 'both'))

# One observation per tested family; zero selected GO terms stay in the set.
ann <- data.frame(FamilyID = families$FamilyID,
  go_ids = c('GO:G', 'GO:L', 'GO:U', 'GO:U', 'GO:U'), go_aspects = 'BP',
  go_terms = c('gain', 'loss', 'ubiquitous', 'ubiquitous', 'ubiquitous'))
candidate <- unique(ann[, c('go_ids', 'go_aspects', 'go_terms')])
g <- e$summarise_specific_go(inc, rbind(ann, ann), candidate)
stopifnot(setequal(g$all$go_ids, candidate$go_ids),
          g$all$n_specific_in_go[g$all$go_ids == 'GO:G'] == 1,
          all(g$all$p_value[g$all$go_ids != 'GO:G'] == 1))
none <- inc; none$selected <- FALSE
g0 <- e$summarise_specific_go(none, ann, candidate)
stopifnot(nrow(g0$all) == 3, all(g0$all$p_value == 1), nrow(g0$significant) == 0)
stopifnot(nrow(e$summarise_specific_go(none, ann, candidate[FALSE, ])$all) == 0)

# Legacy default remains unchanged, including its original GO filtering rule.
events <- data.frame(FamilyID = paste0('F', 1:100), is_target = c(rep(TRUE, 10), rep(FALSE, 90)))
rows <- events[c(1:4, 11:16), ]; rows$go_ids <- 'GO:T'; rows$go_aspects <- 'BP'; rows$go_terms <- 'target'
common <- events; common$go_ids <- 'GO:C'; common$go_aspects <- 'BP'; common$go_terms <- 'common'
rows <- rbind(rows, common)
for (i in 1:8) {
  z <- events[17+i, ]; z$go_ids <- paste0('GO:O', i); z$go_aspects <- 'BP'; z$go_terms <- paste0('other', i)
  rows <- rbind(rows, z)
}
legacy <- e$summarise_go_enrichment(rows, 10, 90, .05)
stopifnot(nrow(legacy$all) == 2, abs(legacy$all$p_value_adjusted[legacy$all$go_ids == 'GO:T'] - .01644975288516) < 1e-12)

cat('CAFE GO native selection and legacy numerical regressions passed\n')

# Eight-argument legacy CLI and explicit event mode must remain equivalent,
# without any of the new native-comparison input files being present.
root <- normalizePath('.')
tmp <- tempfile('cafe-go-default-'); dir.create(tmp)
tryCatch({
  put <- function(d, name) e$write_tsv_base(d, file.path(tmp, name))
  put(data.frame(FamilyID = events$FamilyID, A = as.integer(events$is_target),
                 B = as.integer(!events$is_target)), 'change.tsv')
  put(data.frame(FamilyID = events$FamilyID, A = .001, B = .001), 'prob.tsv')
  put(data.frame(Orthogroup = events$FamilyID, ref = paste0('g', 1:100)), 'ids.tsv')
  annotation <- rows[, c('go_ids', 'go_aspects', 'go_terms')]
  annotation$gene_id <- paste0('g', match(rows$FamilyID, events$FamilyID))
  put(annotation, 'ref.annotation.tsv')
  run <- function(out, extra = character()) {
    argv <- c(file.path(root, 'workflow/support/cafe_go_enrichment.r'),
              file.path(tmp, c('change.tsv', 'prob.tsv', 'ids.tsv', 'ref.annotation.tsv')),
              file.path(tmp, out), 'A', 'increase', 'BP', extra)
    result <- suppressWarnings(system2(file.path(R.home('bin'), 'Rscript'), shQuote(argv), stdout = TRUE, stderr = TRUE))
    if (!is.null(attr(result, 'status'))) stop(paste(result, collapse = '\n'))
  }
  run('default'); run('explicit', 'event')
  for (name in list.files(file.path(tmp, 'default'))) {
    stopifnot(identical(readLines(file.path(tmp, 'default', name)), readLines(file.path(tmp, 'explicit', name))))
  }
  actual <- e$read_tsv_base(file.path(tmp, 'default', 'enrichment_significant_increase_A_all_go.tsv'))
  stopifnot(nrow(actual) == 2, abs(actual$p_value_adjusted[actual$go_ids == 'GO:T'] - .01644975288516) < 1e-12)
}, finally = unlink(tmp, recursive = TRUE))
cat('CAFE GO default CLI equivalence passed\n')

# A bad new input must not leave old optional GO summaries published.
tmp <- tempfile('cafe-go-failure-'); dir.create(tmp)
tryCatch({
  names <- c('family_specificity.tsv', 'specificity_metadata.tsv',
    'enrichment_significant_both_A<1>_all_go.tsv',
    'enrichment_significant_both_A<1>_significant_go.tsv')
  for (name in names) writeLines('old success', file.path(tmp, name))
  dir.create(file.path(tmp, 'native_cafe'))
  writeLines('retained native evidence', file.path(tmp, 'native_cafe', 'cache.txt'))
  argv <- c(file.path(root, 'workflow/support/cafe_go_enrichment.r'),
    rep(file.path(tmp, 'missing.tsv'), 4), tmp, 'A<1>', 'both', 'BP', 'cafe_lrt')
  result <- suppressWarnings(system2(file.path(R.home('bin'), 'Rscript'), shQuote(argv), stdout=TRUE, stderr=TRUE))
  stopifnot(!is.null(attr(result, 'status')), !any(file.exists(file.path(tmp, names))),
    file.exists(file.path(tmp, 'native_cafe', 'cache.txt')))
}, finally=unlink(tmp, recursive=TRUE))
cat('CAFE GO failed-input publication regression passed\n')
