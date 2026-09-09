# Exercise the input boundary without running the plotting program's CLI.
expressions = parse('workflow/support/plot_query2family_presence_absence.R')
for (expression in expressions) {
  if (is.call(expression) && as.character(expression[[1]]) %in% c('=', '<-') &&
      identical(expression[[2]], as.name('extract_ci_table'))) {
    eval(expression)
  }
}
stopifnot(exists('extract_ci_table'))
check_input = function(text, valid) {
  path = tempfile(fileext = '.nwk')
  on.exit(unlink(path))
  writeLines(text, path)
  result = tryCatch(extract_ci_table(paste(readLines(path), collapse = "\n")), error = function(e) e)
  stopifnot(inherits(result, 'error') != valid)
  if (valid) stopifnot(inherits(result$tree_ci, "phylo"), length(result$tree_ci$tip.label) >= 2L)
}
check_input('(A:2,(B:1,C:1):1);', TRUE)
check_input('((B:1,C:1):1,A:2);', TRUE)
check_input('(A:1,B:1);', TRUE)
check_input(c('#NEXUS', 'BEGIN TREES;', 'END;'), FALSE)
check_input('(A:1,B:1)', FALSE)
check_input('(A:1,B:1);garbage', FALSE)
check_input('(A:1,B:1);(C:1,D:1);', FALSE)
check_input('();', FALSE)
check_input('(A:1,B:1;', FALSE)
check_input('(A:1,B:1));', FALSE)
check_input(c('((A:1,B:1):1,C:2);', '((A:1,B:1)[&95%HPD={1,2}]:1,C:2);'), FALSE)
check_input(c('UTREE 1 = (A:1,B:1);', 'UTREE 2 = (A:1,B:1);'), FALSE)
cat('dated-tree interval input regression passed\n')

check_input(c('#NEXUS', 'BEGIN TREES;', 'TREE dated = [&R] (A:1,', 'B:1)[&95%HPD={0.5,1.5}];', 'END;'), TRUE)

# The same named clade must receive its interval in the query-presence plotting consumer.
for (script in c('workflow/support/plot_query2family_presence_absence.R')) {
  for (expression in parse(script)) {
    if (is.call(expression) && as.character(expression[[1]]) %in% c('=', '<-') &&
        identical(expression[[2]], as.name('extract_ci_table'))) eval(expression)
  }
  for (text in c(
      '((A:1,B:1)named:1[&95%HPD={0.5,1.5}],C:2);',
      '((A:1,B:1)named:1[&&NHX:age_ci_low=0.5:age_ci_high=1.5:age_ci_kind=HPD:age_ci_level=0.95],C:2);')) {
    parsed = extract_ci_table(text)
    stopifnot(nrow(parsed$ci_table) == 1L,
              parsed$ci_table$lower == 0.5, parsed$ci_table$upper == 1.5)
    clade = ape::extract.clade(parsed$tree_ci, parsed$ci_table$node_ci)
    stopifnot(setequal(clade$tip.label, c('A', 'B')))
  }
  stopifnot(nrow(extract_ci_table('(A:1,B:1)95;')$ci_table) == 0L)
}
cat('NHX node-age mapping regression passed\n')

# Explicit age attributes can accompany CIs without being dropped to Newick.
rich = extract_ci_table('(A:1,B:1)[&&NHX:age=1:age_ci_low=0.5:age_ci_high=1.5:age_ci_kind=HPD:age_ci_level=0.95];')
stopifnot(nrow(rich$ci_table) == 1L, all(rich$tree_ci$edge.length == 1))
