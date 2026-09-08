# Exercise the input boundary without running the plotting program's CLI.
expressions = parse('workflow/support/plot_mcmctreer.r')
for (expression in expressions) {
  if (is.call(expression) && as.character(expression[[1]]) %in% c('=', '<-') &&
      identical(expression[[2]], as.name('extract_tree_strings'))) {
    eval(expression)
  }
}
stopifnot(exists('extract_tree_strings'))
check_input = function(text, valid) {
  path = tempfile(fileext = '.nwk')
  on.exit(unlink(path))
  writeLines(text, path)
  result = tryCatch(extract_tree_strings(path), error = function(e) e)
  stopifnot(inherits(result, 'error') != valid)
  if (valid) stopifnot(nzchar(result$mean_tree), nzchar(result$ci_tree))
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
check_input(c('((A:1,B:1):1,C:2);', '((A:1,B:1)[&95%HPD={1,2}]:1,C:2);'), TRUE)
check_input(c('UTREE 1 = (A:1,B:1);', 'UTREE 2 = (A:1,B:1);'), TRUE)
cat('plot_mcmctreer input regression passed\n')
