suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(rkftools))

cli_args = commandArgs(trailingOnly = TRUE)
args = cli_args
cat('Arguments:\n')
args = rkftools::get_parsed_args(args, print = TRUE)
mtree_file = args[['infile']]
if (is.null(args[['outfile']]) || nchar(args[['outfile']]) == 0) {
  infile_prefix = sub('\\.[^.]*$', '', mtree_file)
  outfile = paste0(infile_prefix, '.pdf')
} else {
  outfile = args[['outfile']]
}

extract_tree_strings = function(infile) {
  lines = readLines(infile, warn = FALSE)
  lines = trimws(lines)
  lines = lines[nchar(lines) > 0]
  if (length(lines) == 0) {
    stop(paste0('Input file is empty: ', infile))
  }

  utree_lines = grep('UTREE', lines, value = TRUE)
  if (length(utree_lines) > 0) {
    to_newick = function(x) {
      y = sub('^.*UTREE[[:space:]]+[0-9]+[[:space:]]*=[[:space:]]*', '', x)
      semicolon_pos = regexpr(';', y, fixed = TRUE)[1]
      if (semicolon_pos > 0) {
        y = substr(y, 1, semicolon_pos)
      }
      y
    }
    ci_tree = to_newick(utree_lines[1])
    mean_tree = ifelse(length(utree_lines) >= 2, to_newick(utree_lines[2]), ci_tree)
    return(list(ci_tree = ci_tree, mean_tree = mean_tree, source = 'UTREE'))
  }

  tree_lines = grep('^\\(\\(', lines, value = TRUE)
  if (length(tree_lines) == 0) {
    stop(paste0('No Newick tree line found in: ', infile))
  }
  if (length(tree_lines) == 1) {
    return(list(ci_tree = tree_lines[1], mean_tree = tree_lines[1], source = 'FIGTREE_SINGLE'))
  }

  hpd_lines = tree_lines[grepl('95%HPD|95%=', tree_lines)]
  non_hpd_lines = tree_lines[!grepl('95%HPD|95%=', tree_lines)]
  colon_count = vapply(non_hpd_lines, function(x) lengths(regmatches(x, gregexpr(':', x, fixed = TRUE))), integer(1))
  if (length(non_hpd_lines) > 0) {
    mean_tree = non_hpd_lines[which.max(colon_count)]
  } else {
    mean_tree = tree_lines[which.max(vapply(tree_lines, function(x) lengths(regmatches(x, gregexpr(':', x, fixed = TRUE))), integer(1)))]
  }
  if (length(hpd_lines) > 0) {
    ci_tree = hpd_lines[1]
    return(list(ci_tree = ci_tree, mean_tree = mean_tree, source = 'FIGTREE_HPD'))
  }
  list(ci_tree = tree_lines[1], mean_tree = mean_tree, source = 'FIGTREE_DOUBLE')
}

extract_numeric_ci_from_label = function(label_text) {
  if (is.na(label_text) || nchar(trimws(label_text)) == 0) {
    return(NULL)
  }
  nums = regmatches(label_text, gregexpr('[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?', label_text, perl = TRUE))[[1]]
  nums = suppressWarnings(as.numeric(nums))
  nums = nums[is.finite(nums)]
  if (length(nums) < 2) {
    return(NULL)
  }
  c(min(nums[1], nums[2]), max(nums[1], nums[2]))
}

extract_ci_table = function(ci_tree_text) {
  ci_pattern = '\\[&95%(?:HPD)?=\\{[[:space:]]*[-+0-9.eE]+[[:space:]]*,[[:space:]]*[-+0-9.eE]+[[:space:]]*\\}\\]'
  ci_hits = gregexpr(ci_pattern, ci_tree_text, perl = TRUE)
  hit_texts = regmatches(ci_tree_text, ci_hits)[[1]]
  ci_tree_tagged = ci_tree_text

  ci_table = data.frame(
    tag = character(0),
    lower = numeric(0),
    upper = numeric(0),
    stringsAsFactors = FALSE
  )

  if (length(hit_texts) == 0) {
    cleaned = gsub('\\[&[^]]*\\]', '', ci_tree_text, perl = TRUE)
    tree_ci = read.tree(text = cleaned)
  } else {
    for (i in seq_along(hit_texts)) {
      h = hit_texts[i]
      mm = regexec('\\{[[:space:]]*([-+0-9.eE]+)[[:space:]]*,[[:space:]]*([-+0-9.eE]+)[[:space:]]*\\}', h, perl = TRUE)
      cap = regmatches(h, mm)[[1]]
      if (length(cap) < 3) {
        next
      }
      nums = as.numeric(cap[2:3])
      if (any(is.na(nums))) {
        next
      }
      lo = min(nums[1], nums[2])
      hi = max(nums[1], nums[2])
      tag = paste0('CI_TAG_', i)
      ci_tree_tagged = sub(ci_pattern, tag, ci_tree_tagged, perl = TRUE)
      ci_table = rbind(
        ci_table,
        data.frame(tag = tag, lower = lo, upper = hi, stringsAsFactors = FALSE)
      )
    }
    ci_tree_tagged = gsub('\\[&[^]]*\\]', '', ci_tree_tagged, perl = TRUE)
    tree_ci = read.tree(text = ci_tree_tagged)
  }

  # Some FigTree outputs store CIs as internal-node labels (e.g. "81.3,96.9")
  node_labels = tree_ci$node.label
  if (!is.null(node_labels) && length(node_labels) > 0) {
    n_tip = length(tree_ci$tip.label)
    for (i in seq_along(node_labels)) {
      ci_vals = extract_numeric_ci_from_label(node_labels[i])
      if (is.null(ci_vals)) {
        next
      }
      node_id = n_tip + i
      ci_table = rbind(
        ci_table,
        data.frame(
          tag = paste0('NODE_LABEL_', i),
          lower = ci_vals[1],
          upper = ci_vals[2],
          node_ci = node_id,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  if (nrow(ci_table) == 0) {
    return(list(tree_ci = tree_ci, ci_table = data.frame()))
  }

  if (!'node_ci' %in% colnames(ci_table)) {
    ci_table$node_ci = NA_integer_
  }
  idx_missing_node = is.na(ci_table$node_ci)
  if (any(idx_missing_node)) {
    n_tip = length(tree_ci$tip.label)
    node_labels = tree_ci$node.label
    ci_table$node_ci[idx_missing_node] = n_tip + match(ci_table$tag[idx_missing_node], node_labels)
  }
  ci_table = ci_table[!is.na(ci_table$node_ci), c('tag', 'lower', 'upper', 'node_ci'), drop = FALSE]
  list(tree_ci = tree_ci, ci_table = ci_table)
}

get_clade_signature = function(tree, node_id) {
  clade = extract.clade(tree, node = node_id)
  paste(sort(clade$tip.label), collapse = '|')
}

map_ci_nodes_to_mean_tree = function(tree_ci, tree_mean, ci_table) {
  if (nrow(ci_table) == 0) {
    return(ci_table)
  }
  n_tip_mean = length(tree_mean$tip.label)
  n_node_mean = tree_mean$Nnode
  sig_to_node_mean = integer(0)
  for (i in seq_len(n_node_mean)) {
    node_id = n_tip_mean + i
    sig = get_clade_signature(tree_mean, node_id)
    sig_to_node_mean[sig] = node_id
  }

  ci_table$signature = vapply(ci_table$node_ci, function(x) get_clade_signature(tree_ci, x), character(1))
  ci_table$node_mean = unname(sig_to_node_mean[ci_table$signature])
  ci_table = ci_table[!is.na(ci_table$node_mean), , drop = FALSE]
  ci_table
}

plot_dated_tree = function(tree_mean, ci_table, out_pdf) {
  tree_mean = ladderize(tree_mean, right = TRUE)
  n_tip = length(tree_mean$tip.label)
  root_age = max(node.depth.edgelength(tree_mean))
  if (!is.finite(root_age) || root_age <= 0) {
    root_age = 1
  }

  max_time = root_age
  if (nrow(ci_table) > 0 && any(is.finite(ci_table$upper))) {
    max_time = max(max_time, max(ci_table$upper, na.rm = TRUE))
  }

  left_extra = max(0, max_time - root_age)
  left_pad = max(max_time * 0.02, root_age * 0.01)
  x_left = -left_extra - left_pad
  max_tip_chars = max(nchar(tree_mean$tip.label))
  x_right = root_age + max(root_age * 0.03, root_age * (max_tip_chars * 0.012))
  axis_y = 0.35

  pdf_height_inch = max(3.0, 1.2 + (n_tip * 0.24))
  pdf(file = out_pdf, width = 7.0, height = pdf_height_inch)
  par(mar = c(3.9, 0.4, 0.3, 0.4), xpd = NA)
  plot(
    tree_mean,
    show.tip.label = TRUE,
    cex = 0.7,
    no.margin = FALSE,
    x.lim = c(x_left, x_right),
    y.lim = c(0.6, n_tip + 0.4),
    label.offset = max(0.001, root_age * 0.008)
  )
  plot_info = get('last_plot.phylo', envir = .PlotPhyloEnv)

  ticks_time = pretty(c(0, max_time), n = 6)
  ticks_time = ticks_time[ticks_time >= 0 & ticks_time <= max_time]
  if (!any(abs(ticks_time - 0) < 1e-12)) {
    ticks_time = c(0, ticks_time)
  }
  ticks_time = sort(unique(ticks_time))
  ticks_pos = root_age - ticks_time
  max_tick = max(ticks_time)
  tick_decimals = if (max_tick >= 20) 0 else if (max_tick >= 2) 1 else 2
  tick_labels = sprintf(paste0('%.', tick_decimals, 'f'), ticks_time)
  axis(1, at = ticks_pos, labels = tick_labels, cex.axis = 0.75, lwd = 1, lwd.ticks = 1, pos = axis_y)
  mtext('Time before present (Ma)', side = 1, line = 1.9, cex = 0.85)

  if (nrow(ci_table) > 0) {
    y = plot_info$yy[ci_table$node_mean]
    x_old = root_age - ci_table$upper
    x_young = root_age - ci_table$lower
    x_old = pmax(x_left, pmin(x_right, x_old))
    x_young = pmax(x_left, pmin(x_right, x_young))
    segments(x0 = x_old, y0 = y, x1 = x_young, y1 = y, col = 'navy', lwd = 2)
    node_x = plot_info$xx[ci_table$node_mean]
    points(node_x, y, pch = 16, cex = 0.4, col = 'navy')
  }

  dev.off()
}

tree_strings = extract_tree_strings(mtree_file)
mean_tree = read.tree(text = tree_strings$mean_tree)
ci_parsed = extract_ci_table(tree_strings$ci_tree)
ci_table = map_ci_nodes_to_mean_tree(ci_parsed$tree_ci, mean_tree, ci_parsed$ci_table)
if (nrow(ci_table) == 0) {
  warning('No parsable time CI annotations were found in the input tree. Plotting the mean dated tree without CI bars.')
}

plot_dated_tree(mean_tree, ci_table, outfile)
cat('plot_mcmctreer.r done!\n')
