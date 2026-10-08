# Bounded renderer. Results attest rendering only; callers still record and
# verify each family's declared provenance and exact workflow attempt.
argv = commandArgs(trailingOnly=TRUE)
if (length(argv) != 1) stop('Usage: Rscript tree_plot_batch.r PLAN.json')
if (file.info(argv[[1]])$size > 1024 * 1024) stop('Plot plan exceeds 1 MiB')
jobs = jsonlite::fromJSON(argv[[1]], simplifyVector=FALSE)
if (!is.list(jobs) || length(jobs) < 1 || length(jobs) > 32) stop('Expected 1..32 plot jobs')
valid_text = function(x) is.character(x) && length(x) == 1 && !is.na(x) && nzchar(x)
for (job in jobs) {
  if (!is.list(job) || !all(c('id','cwd','args','output') %in% names(job)) ||
      !all(vapply(job[c('id','cwd','output')], valid_text, logical(1))) ||
      !is.list(job$args) || length(job$args) > 256 ||
      !all(vapply(job$args, valid_text, logical(1))) ||
      (!is.null(job$check_ggimage) && (!is.logical(job$check_ggimage) || length(job$check_ggimage) != 1 || is.na(job$check_ggimage))) ||
      !startsWith(job$cwd, '/') || !dir.exists(job$cwd) || !startsWith(job$output, '/'))
    stop('Invalid plot job; cwd/output and all input file arguments must be absolute')
}
if (anyDuplicated(vapply(jobs, function(job) job$id, '')) ||
    anyDuplicated(vapply(jobs, function(job) job$output, ''))) stop('Duplicate plot job/output')
script_dir = dirname(normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(), value=TRUE)[[1]])))
renderer_path = file.path(script_dir, 'stat_branch2tree_plot.r')
renderer_bytes = function() {
  size = file.info(renderer_path)$size
  if (is.na(size) || size > 4 * 1024 * 1024) stop('Invalid renderer source size')
  readBin(renderer_path, 'raw', n=size + 1)
}
source_bytes = renderer_bytes()
expressions = parse(text=rawToChar(source_bytes))
render = function(job) {
  env = NULL
  before_dir = getwd()
  before_options = options()
  before_theme = ggplot2::theme_get()
  had_seed = exists('.Random.seed', envir=globalenv(), inherits=FALSE)
  before_seed = if (had_seed) get('.Random.seed', envir=globalenv()) else NULL
  before_env = Sys.getenv(c('TREEVIS_SPECIES_PARSER','GG_TREE_PLOT_CHECK_GGIMAGE'), unset=NA)
  scratch = tempfile('gg-plot-', tmpdir=job$cwd)
  dir.create(scratch)
  on.exit({
    while (grDevices::dev.cur() > 1) grDevices::dev.off()
    added_options = setdiff(names(options()), names(before_options))
    if (length(added_options)) options(setNames(rep(list(NULL), length(added_options)), added_options))
    options(before_options)
    ggplot2::theme_set(before_theme)
    if (had_seed) assign('.Random.seed', before_seed, envir=globalenv()) else if (exists('.Random.seed', envir=globalenv(), inherits=FALSE)) rm('.Random.seed', envir=globalenv())
    for (name in names(before_env)) {
      if (is.na(before_env[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, setNames(list(before_env[[name]]), name))
    }
    setwd(before_dir)
    unlink(scratch, recursive=TRUE)
    env = NULL
    gc(verbose=FALSE)
  })
  result = tryCatch({
    if (!identical(source_bytes, renderer_bytes())) stop('Renderer source changed during batch')
    values = sub('^[^=]*=', '', unlist(job$args, use.names=FALSE))
    tokens = unlist(lapply(values, function(value) {
      # A complete existing file is one argument, even when its name contains
      # commas. Split panel/convergence lists only after that exact-path check.
      if (startsWith(value, '/') && file.exists(value) && !dir.exists(value)) value
      else strsplit(value, ',', fixed=TRUE)[[1]]
    }), use.names=FALSE)
    paths = unique(tokens[startsWith(tokens, '/') & file.exists(tokens) & !dir.exists(tokens)])
    stat_arg = grep('^--stat_branch=', unlist(job$args, use.names=FALSE), value=TRUE)
    if (length(stat_arg) != 1 || !startsWith(sub('^--stat_branch=', '', stat_arg), '/'))
      stop('Batch stat_branch must use one absolute path')
    identities = function() list(info=file.info(paths)[,c('size','mtime','ctime'),drop=FALSE], md5=tools::md5sum(paths))
    before_inputs = identities()
    setwd(scratch)
    if (!is.null(job$species_parser)) Sys.setenv(TREEVIS_SPECIES_PARSER=job$species_parser)
    # The optional dependency decision belongs to this request, never a stale PDF.
    if (!identical(job$check_ggimage, FALSE) && !requireNamespace('ggimage', quietly=TRUE))
      return(list(id=job$id, exit_code=42L, output=job$output))
    Sys.setenv(GG_TREE_PLOT_CHECK_GGIMAGE='0')
    env = new.env(parent=globalenv())
    env$.gg_tree_plot_args = unlist(job$args, use.names=FALSE)
    eval(expressions, envir=env)
    if (!identical(source_bytes, renderer_bytes())) stop('Renderer source changed before publication')
    if (!identical(before_inputs, identities())) stop('Plot input changed before publication')
    output = file.path(scratch, 'stat_branch2tree_plot.pdf')
    if (!file.exists(output) || file.info(output)$size <= 0) stop('Renderer produced no PDF')
    dir.create(dirname(job$output), recursive=TRUE, showWarnings=FALSE)
    staged = tempfile('.gg-plot-', tmpdir=dirname(job$output))
    on.exit(unlink(staged), add=TRUE)
    if (!file.copy(output, staged) || !file.rename(staged, job$output)) stop('PDF publication failed')
    list(id=job$id, exit_code=0L, output=job$output, size_bytes=unname(file.info(job$output)$size))
  }, error=function(e) list(id=job$id, exit_code=1L, output=job$output, detail=conditionMessage(e)))
  result
}
results = lapply(jobs, render)
result_file = paste0(argv[[1]], '.results.json')
staged = tempfile('.gg-plot-results-', tmpdir=dirname(result_file))
jsonlite::write_json(list(schema='genegalleon-plot-batch-result-v1', completion_evidence=FALSE, results=results), staged, auto_unbox=TRUE)
if (!file.rename(staged, result_file)) stop('Result publication failed')
if (any(vapply(results, function(result) result$exit_code != 0L, logical(1)))) quit(status=1)
