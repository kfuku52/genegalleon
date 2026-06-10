suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(rkftools))

cli_args = commandArgs(trailingOnly = TRUE)
cat('Arguments:\n')
args = cli_args
args = rkftools::get_parsed_args(args, print = TRUE)

font_size = 8
font_size_factor = 0.352777778 # https://stackoverflow.com/questions/17311917/ggplot2-the-unit-of-size

has_nonempty_file <- function(file_path) {
  if (is.null(file_path) || length(file_path) == 0 || is.na(file_path) || !nzchar(file_path)) {
    return(FALSE)
  }
  if (!file.exists(file_path)) {
    return(FALSE)
  }
  file_info <- file.info(file_path)
  !is.na(file_info$size[[1]]) && file_info$size[[1]] > 0
}

has_value <- function(x) {
  !(is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x)))
}


extract_param_from_iqtree_log <- function(file_path) {
  if (!has_nonempty_file(file_path)) {
    return(list(model = NA_character_, nsite = NA_integer_, loglik = NA_real_, version = NA_character_))
  }
  # Read the contents of the file
  log_text <- readLines(file_path, warn = FALSE)
  full_text <- paste(log_text, collapse = "\n")

  # Regular expressions for the desired values
  model_pattern <- "-m\\s+([A-Za-z+\\d]+)"
  nsite_pattern <- "Alignment has \\d+ sequences with (\\d+) columns"
  loglik_pattern <- "BEST SCORE FOUND :\\s+(-?\\d+\\.\\d+)"
  version_pattern <- "version (\\d+(?:\\.\\d+)*(?:[a-zA-Z])?)"

  # Using regexpr and regmatches to extract the matches
  model_match <- regexpr(model_pattern, full_text)
  nsite_match <- regexpr(nsite_pattern, full_text)
  loglik_match <- regexpr(loglik_pattern, full_text)
  version_match <- regexpr(version_pattern, full_text)

  model <- ifelse(model_match[1] != -1, regmatches(full_text, model_match)[[1]], NA)
  nsite <- ifelse(nsite_match[1] != -1, regmatches(full_text, nsite_match)[[1]], NA)
  loglik <- ifelse(loglik_match[1] != -1, regmatches(full_text, loglik_match)[[1]], NA)
  version <- ifelse(version_match[1] != -1, regmatches(full_text, version_match)[[1]], NA)

  # Extracting the actual values
  model_value <- ifelse(!is.na(model), sub(model_pattern, "\\1", model), NA)
  nsite_value <- ifelse(!is.na(nsite), as.integer(sub(nsite_pattern, "\\1", nsite)), NA)
  loglik_value <- ifelse(!is.na(loglik), as.numeric(sub(loglik_pattern, "\\1", loglik)), NA)
  version_value <- ifelse(!is.na(version), sub(version_pattern, "\\1", version), NA)

  return(list(model = model_value, nsite = nsite_value, loglik = loglik_value, version = version_value))
}

extract_first_group <- function(full_text, patterns) {
  for (pattern in patterns) {
    reg_match <- regexec(pattern, full_text, perl = TRUE)
    extracted <- regmatches(full_text, reg_match)[[1]]
    if (length(extracted) >= 2) {
      return(extracted[2])
    }
  }
  return(NA_character_)
}

extract_param_from_astral_log <- function(file_path) {
  if (!has_nonempty_file(file_path)) {
    return(list(ntree = NA_integer_, version = NA_character_))
  }
  log_text <- readLines(file_path, warn = FALSE)
  full_text <- paste(log_text, collapse = "\n")
  ntree_str <- extract_first_group(full_text, c(
    "#Genetrees: (\\d+)",
    "Number of gene trees: (\\d+)"
  ))
  version_value <- extract_first_group(full_text, c(
    "Version: v(\\d+(?:\\.\\d+)*(?:[a-zA-Z])?)",
    "This is ASTRAL version (\\d+(?:\\.\\d+)*(?:[a-zA-Z])?)"
  ))
  ntree_value <- suppressWarnings(as.integer(ntree_str))
  return(list(ntree = ntree_value, version = version_value))
}

fmt_num <- function(x, digits = 0) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !is.finite(x)) {
    return("NA")
  }
  formatC(as.numeric(x), format = "f", digits = digits, big.mark = ",")
}

build_iqtree_label <- function(sequence_label, params) {
  if (has_value(params[['version']]) && has_value(params[['nsite']]) && has_value(params[['model']]) && has_value(params[['loglik']]) && has_value(params[['ngene']])) {
    return(sprintf(
      'IQ-TREE v%s with concatenated %s alignment\n(%s genes, %s sites, %s model, loglik = %s)',
      params[['version']],
      sequence_label,
      fmt_num(params[['ngene']], digits = 0),
      fmt_num(params[['nsite']], digits = 0),
      params[['model']],
      fmt_num(params[['loglik']], digits = 1)
    ))
  }
  sprintf('IQ-TREE concatenated %s alignment', sequence_label)
}

build_astral_label <- function(sequence_label, astral_params, iqtree_params) {
  if (has_value(astral_params[['version']]) && has_value(astral_params[['ntree']]) && has_value(iqtree_params[['model']])) {
    return(sprintf(
      'ASTRAL v%s with %s %s ML gene trees (%s)',
      astral_params[['version']],
      fmt_num(astral_params[['ntree']], digits = 0),
      sequence_label,
      iqtree_params[['model']]
    ))
  }
  sprintf('ASTRAL %s gene trees', sequence_label)
}

iqtree_dna_params = extract_param_from_iqtree_log(file_path = args[['iqtree_dna_log']])
iqtree_pep_params = extract_param_from_iqtree_log(file_path = args[['iqtree_pep_log']])
astral_dna_params = extract_param_from_astral_log(file_path = args[['astral_dna_log']])
astral_pep_params = extract_param_from_astral_log(file_path = args[['astral_pep_log']])
iqtree_dna_params[['ngene']] = astral_dna_params[['ntree']]
iqtree_pep_params[['ngene']] = astral_pep_params[['ntree']]

infer_iqtree_support_nwk <- function(rooted_nwk) {
  if (!has_value(rooted_nwk)) {
    return(NA_character_)
  }
  candidate <- sub("\\.rooted\\.nwk$", ".iqtree.nwk", rooted_nwk)
  if (has_nonempty_file(candidate)) {
    return(candidate)
  }
  return(NA_character_)
}

descendant_tips <- function(tree, node) {
  n_tip <- length(tree[['tip.label']])
  children <- tree[['edge']][tree[['edge']][, 1] == node, 2]
  tips <- children[children <= n_tip]
  internal_children <- children[children > n_tip]
  for (child in internal_children) {
    tips <- c(tips, descendant_tips(tree, child))
  }
  tips
}

node_split_key <- function(tree, node) {
  n_tip <- length(tree[['tip.label']])
  tips <- sort(tree[['tip.label']][descendant_tips(tree, node)])
  if (length(tips) == 0 || length(tips) == n_tip) {
    return(NA_character_)
  }
  other <- sort(setdiff(tree[['tip.label']], tips))
  if (length(other) < length(tips)) {
    tips <- other
  } else if (length(other) == length(tips)) {
    tip_key <- paste(tips, collapse = "\t")
    other_key <- paste(other, collapse = "\t")
    if (other_key < tip_key) {
      tips <- other
    }
  }
  paste(tips, collapse = "\t")
}

transfer_internal_labels <- function(target_tree, source_tree) {
  n_tip <- length(target_tree[['tip.label']])
  if (!setequal(target_tree[['tip.label']], source_tree[['tip.label']])) {
    warning('Cannot transfer support values: target and source trees have different tip labels.')
    return(target_tree)
  }
  source_nodes <- seq.int(length(source_tree[['tip.label']]) + 1, length(source_tree[['tip.label']]) + source_tree[['Nnode']])
  source_labels <- source_tree[['node.label']]
  label_by_split <- list()
  for (idx in seq_along(source_nodes)) {
    label <- source_labels[[idx]]
    if (is.na(label) || !nzchar(label)) {
      next
    }
    split_key <- node_split_key(source_tree, source_nodes[[idx]])
    if (!is.na(split_key) && nzchar(split_key)) {
      label_by_split[[split_key]] <- label
    }
  }

  target_nodes <- seq.int(n_tip + 1, n_tip + target_tree[['Nnode']])
  target_labels <- rep('', target_tree[['Nnode']])
  for (idx in seq_along(target_nodes)) {
    split_key <- node_split_key(target_tree, target_nodes[[idx]])
    if (!is.na(split_key) && nzchar(split_key) && !is.null(label_by_split[[split_key]])) {
      target_labels[[idx]] <- label_by_split[[split_key]]
    }
  }
  target_tree[['node.label']] <- target_labels
  target_tree
}

tree_specs = list(
  list(label = build_iqtree_label('DNA', iqtree_dna_params), nwk = args[['iqtree_dna_nwk']], support_nwk = infer_iqtree_support_nwk(args[['iqtree_dna_nwk']]), support_note = 'Support: ultrafast bootstrap (%).\nBranch length: nucleotide substitutions per site.'),
  list(label = build_iqtree_label('protein', iqtree_pep_params), nwk = args[['iqtree_pep_nwk']], support_nwk = infer_iqtree_support_nwk(args[['iqtree_pep_nwk']]), support_note = 'Support: ultrafast bootstrap (%).\nBranch length: amino acid substitutions per site.'),
  list(label = build_astral_label('DNA', astral_dna_params, iqtree_dna_params), nwk = args[['astral_dna_nwk']], support_note = 'Support: ASTRAL local posterior probability-like support.\nBranch length: nucleotide substitutions per site.'),
  list(label = build_astral_label('protein', astral_pep_params, iqtree_pep_params), nwk = args[['astral_pep_nwk']], support_note = 'Support: ASTRAL local posterior probability-like support.\nBranch length: amino acid substitutions per site.')
)
tree_specs = Filter(function(spec) has_nonempty_file(spec[['nwk']]), tree_specs)

if (length(tree_specs) == 0) {
  stop('No species tree files were provided.')
}

trees = list()
for (spec in tree_specs) {
  key = spec[['label']]
  trees[[key]] = read.tree(spec[['nwk']])
  if (has_value(spec[['support_nwk']]) && has_nonempty_file(spec[['support_nwk']])) {
    trees[[key]] = transfer_internal_labels(trees[[key]], read.tree(spec[['support_nwk']]))
  }
  trees[[key]] = ladderize(trees[[key]])
}

plots = list()
for (spec in tree_specs) {
  key = spec[['label']]
  cat('Processing:', key, '\n')
  plots[[key]] = ggtree(trees[[key]])
  isTip = plots[[key]][['data']][['isTip']]
  if (grepl('concatenated', key)) {
    support_color = 'red'
    support_values = suppressWarnings(as.numeric(plots[[key]][['data']][!isTip,][['label']]))
    support_values[!is.finite(support_values)] = NA_real_
    support_labels = ifelse(
      is.na(support_values),
      '',
      format(round(support_values, digits = 0), nsmall = 0)
    )
    plots[[key]][['data']][!isTip, 'label'] = support_labels
  } else if (grepl('ASTRAL', key)) {
    support_color = 'blue'
    support_values = suppressWarnings(as.numeric(plots[[key]][['data']][!isTip,][['label']]))
    support_values[!is.finite(support_values)] = NA_real_
    support_labels = ifelse(
      is.na(support_values),
      '',
      format(round(support_values, digits = 2), nsmall = 2)
    )
    plots[[key]][['data']][!isTip, 'label'] = support_labels
  }
  plots[[key]][['data']][, 'label'] = gsub('_', ' ', plots[[key]][['data']][['label']])
  is_root = (plots[[key]][['data']][['parent']] == plots[[key]][['data']][['node']])
  plots[[key]][['data']][is_root, 'label'] = ''
  root_num = c(plots[[key]][['data']][is_root, 'node'])
  is_subroot = (plots[[key]][['data']][['parent']] == root_num)
  plots[[key]][['data']][is_subroot & (!isTip), 'label'] = ''
  plots[[key]] = plots[[key]] + geom_tiplab(fontface = 'italic', size = font_size * font_size_factor)
  plots[[key]] = plots[[key]] + hexpand(1, direction = 1)
  plots[[key]] = plots[[key]] + geom_text2(aes(subset = !isTip, label = label, x = branch, y = y + 0.5), size = font_size * font_size_factor, color = support_color)
  plots[[key]] = plots[[key]] + geom_treescale(x = max(plots[[key]][['data']][['x']]) * 1.5, y = 1, fontsize = font_size * font_size_factor)
  plots[[key]] = plots[[key]] + labs(caption = spec[['support_note']])
  plots[[key]] = plots[[key]] + theme(
    plot.caption = element_text(size = font_size * 0.85, hjust = 0),
    strip.text.x = element_text(size = font_size),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
  )
}

max_tip_count = max(vapply(trees, function(tree) length(tree[['tip.label']]), integer(1)))
plot_ncol = min(2, length(plots))
plot_nrow = ceiling(length(plots) / plot_ncol)
panel_height = max(2.6, 0.95 + max_tip_count * font_size / 72 * 1.03)
pdf_height = 0.5 + plot_nrow * panel_height
out = plot_list(gglist = plots, ncol = plot_ncol)
ggsave('species_trees.pdf', plot = out, width = 7.2, height = pdf_height, units = 'in')
out
