suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(rkftools))

cli_args = commandArgs(trailingOnly = TRUE)
mode = ifelse('--debug' %in% cli_args, 'debug', 'batch')
cli_args = cli_args[cli_args != '--debug']

if (mode == 'debug') {
  dir_pg = '/Users/kf/Dropbox/data/evolutionary_transcriptomics/20231016_pg_pipeline/workspace'
  dir_sptree = file.path(dir_pg, 'species_tree')
  args = c()
  args = c(args, paste0('--infile=', file.path(dir_sptree, 'constrained_tree', 'constrained.nwk')))
  args = c(args, paste0('--outfile=', file.path(dir_sptree, 'constrained_tree', 'constrained_tree_constraints.pdf')))
} else {
  args = cli_args
}

cat('Arguments:\n')
args = rkftools::get_parsed_args(args, print = TRUE)

infile = args[['infile']]
if (is.null(infile) || nchar(infile) == 0) {
  stop('Missing required argument: --infile')
}

if (is.null(args[['outfile']]) || nchar(args[['outfile']]) == 0) {
  outfile = file.path(dirname(infile), 'constrained_tree_constraints.pdf')
} else {
  outfile = args[['outfile']]
}

parse_constraint_label = function(label_text) {
  if (is.null(label_text) || is.na(label_text) || nchar(trimws(label_text)) == 0) {
    return(NULL)
  }
  s = trimws(gsub("^'+|'+$", '', label_text))
  m = regexec(
    'B\\(([+-]?[0-9]*\\.?[0-9]+(?:[eE][+-]?[0-9]+)?),([+-]?[0-9]*\\.?[0-9]+(?:[eE][+-]?[0-9]+)?)(?:,([+-]?[0-9]*\\.?[0-9]+(?:[eE][+-]?[0-9]+)?),([+-]?[0-9]*\\.?[0-9]+(?:[eE][+-]?[0-9]+)?))?\\)',
    s,
    perl = TRUE
  )
  hit = regmatches(s, m)[[1]]
  if (length(hit) < 3) {
    return(NULL)
  }
  lo = suppressWarnings(as.numeric(hit[2]))
  hi = suppressWarnings(as.numeric(hit[3]))
  if (!is.finite(lo) || !is.finite(hi)) {
    return(NULL)
  }
  if (lo > hi) {
    tmp = lo
    lo = hi
    hi = tmp
  }
  list(lower = lo, upper = hi)
}

extract_constraints = function(tree) {
  out = data.frame()
  if (is.null(tree$node.label) || length(tree$node.label) == 0) {
    return(out)
  }
  n_tip = length(tree$tip.label)
  rows = list()
  k = 0
  for (i in seq_along(tree$node.label)) {
    parsed = parse_constraint_label(tree$node.label[i])
    if (is.null(parsed)) {
      next
    }
    k = k + 1
    rows[[k]] = data.frame(
      node = n_tip + i,
      lower = parsed$lower,
      upper = parsed$upper,
      stringsAsFactors = FALSE
    )
  }
  if (k == 0) {
    return(data.frame())
  }
  do.call(rbind, rows)
}

format_age = function(x) {
  if (!is.finite(x)) {
    return('NA')
  }
  if (abs(x - round(x)) < 1e-6) {
    return(as.character(as.integer(round(x))))
  }
  sprintf('%.1f', x)
}

tree = read.tree(infile)
tree = ladderize(tree, right = TRUE)
constraints = extract_constraints(tree)
if (nrow(constraints) > 0) {
  age_labels = sprintf(
    '%s-%s Ma',
    vapply(constraints$lower, format_age, character(1)),
    vapply(constraints$upper, format_age, character(1))
  )
  max_ci_chars = max(nchar(age_labels))
} else {
  age_labels = character(0)
  max_ci_chars = 0
}

n_tip = length(tree$tip.label)
max_tip_chars = max(nchar(tree$tip.label))
max_depth = max(node.depth(tree)[1:n_tip])
x_left = -max(0.8, max_ci_chars * 0.12)
x_right = max_depth + max(5.0, max_tip_chars * 0.40)

pdf_height_inch = max(2.9, 1.1 + (n_tip * 0.24))
pdf(outfile, width = 7.0, height = pdf_height_inch)
par(mar = c(1.6, 1.0, 0.4, 2.2), xpd = NA)

plot(
  tree,
  use.edge.length = FALSE,
  show.tip.label = TRUE,
  cex = 0.7,
  no.margin = FALSE,
  x.lim = c(x_left, x_right),
  y.lim = c(0.35, n_tip + 0.5),
  label.offset = 0.03
)

plot_info = get('last_plot.phylo', envir = .PlotPhyloEnv)

if (nrow(constraints) > 0) {
  y = plot_info$yy[constraints$node]
  x_node = plot_info$xx[constraints$node]
  node_mark_size = 0.33
  points(x_node, y, pch = 16, cex = node_mark_size, col = 'navy')

  ord = order(y)
  y_jitter = rep(0, length(y))
  if (length(ord) > 1) {
    for (k in 2:length(ord)) {
      i_prev = ord[k - 1]
      i_curr = ord[k]
      if (abs((y[i_curr] + y_jitter[i_curr]) - (y[i_prev] + y_jitter[i_prev])) < 0.35) {
        y_jitter[i_curr] = y_jitter[i_prev] + 0.18
      }
    }
  }

  x_span = diff(range(plot_info$xx))
  x_gap = max(0.02 * x_span, 0.06)

  text(
    x = x_node - x_gap,
    y = y + y_jitter,
    labels = age_labels,
    adj = c(1, 0.5),
    cex = 0.62,
    col = 'navy'
  )
  mtext('Blue labels: constrained node age ranges (Ma)', side = 1, line = 0.2, cex = 0.73, col = 'navy')
} else {
  warning('No parsable B(lower,upper,...) constraints found in node labels.')
  mtext('No parsable constrained ranges found in node labels.', side = 1, line = 0.2, cex = 0.73)
}

dev.off()
cat('plot_constrained_tree.r done!\n')
