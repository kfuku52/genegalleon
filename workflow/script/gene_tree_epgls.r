#!/usr/bin/env Rscript

cli_args = commandArgs(trailingOnly = TRUE)
run_mode = ifelse('--debug' %in% cli_args, 'debug', 'batch')
cli_args = cli_args[cli_args != '--debug']

load_package = function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf('Required package is missing: %s', pkg), call. = FALSE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

load_package('ape')
load_package('ggplot2')
load_package('RRPP')
load_package('rkftools')

script_arg = grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) > 0) {
  script_path = sub('^--file=', '', script_arg[1])
  script_path = gsub('~\\+~', ' ', script_path)
  script_dir = dirname(normalizePath(script_path, mustWork = FALSE))
} else {
  script_dir = getwd()
}
source(file.path(script_dir, 'pgls_common.R'))

if (run_mode == 'debug') {
  dir_pg = '/Users/kf/Dropbox/collaborators/Yudai_Okuyama/20240606_Asarum_PGLS/workspace_debug'
  og_id = 'HOG0000272'
  args = c()
  args = c(args, paste0('--file_tree=', file.path(dir_pg, 'orthogroup', 'dated_tree', paste0(og_id, '.dated.nwk'))))
  args = c(args, paste0('--file_exp=', file.path(dir_pg, 'orthogroup', 'character.expression', paste0(og_id, '.expression.tsv'))))
  args = c(args, paste0('--file_trait=', file.path(dir_pg, 'species_trait', 'species_trait.tsv')))
  args = c(args, paste0('--replicate_sep=', '_'))
  args = c(args, paste0('--merge_replicates=', 'no'))
  args = c(args, paste0('--iter=', 199))
  setwd(dir_pg)
} else {
  args = cli_args
}

cat('arguments:\n')
args = rkftools::get_parsed_args(args, print = TRUE)
if (is.null(args[['replicate_sep']]) || nchar(args[['replicate_sep']]) == 0) {
  args[['replicate_sep']] = '_'
}
if (is.null(args[['merge_replicates']]) || !(args[['merge_replicates']] %in% c('yes', 'no'))) {
  args[['merge_replicates']] = 'yes'
}
if (is.null(args[['ws_subject']]) || !(args[['ws_subject']] %in% c('auto', 'species', 'gene_id'))) {
  args[['ws_subject']] = 'auto'
}
iter = suppressWarnings(as.integer(args[['iter']]))
if (is.na(iter) || iter < 9) {
  iter = 199
}

leaf_to_species = function(leaf_names) {
  leaf_names = as.character(leaf_names)
  leaf_names = sub('_', '|', leaf_names)
  leaf_names = sub('_.*', '', leaf_names)
  leaf_names = sub('\\|', '_', leaf_names)
  return(leaf_names)
}

sort_exp_local = function(exp, tree) {
  if (!('gene_id' %in% colnames(exp))) {
    return(exp)
  }
  ids = as.character(exp[['gene_id']])
  exp = exp[ids %in% tree[['tip.label']], , drop = FALSE]
  ord = match(tree[['tip.label']], as.character(exp[['gene_id']]))
  ord = ord[!is.na(ord)]
  exp = exp[ord, , drop = FALSE]
  rownames(exp) = NULL
  return(exp)
}

get_expression_bases_local = function(exp, replicate_sep = '_') {
  cols = colnames(exp)
  cols = cols[cols != 'gene_id']
  if (length(cols) == 0) {
    return(character(0))
  }
  suffix_pat = paste0(replicate_sep, '[0-9]+$')
  base = sub(suffix_pat, '', cols)
  base = base[!duplicated(base)]
  return(base)
}

sort_exp_fn = if (exists('sort_exp', mode = 'function')) get('sort_exp') else sort_exp_local
get_expression_bases_fn = if (exists('get_expression_bases', mode = 'function')) get('get_expression_bases') else get_expression_bases_local

build_observation_cov_from_gene_tree = function(dat, gene_tree) {
  if (nrow(dat) == 0) {
    return(list(data = dat, Cov = NULL))
  }
  cov_gene = ape::vcv.phylo(gene_tree)
  tip_ids = rownames(cov_gene)
  dat = dat[dat[['gene_id']] %in% tip_ids, , drop = FALSE]
  if (nrow(dat) == 0) {
    return(list(data = dat, Cov = NULL))
  }
  idx = match(dat[['gene_id']], tip_ids)
  cov_obs = cov_gene[idx, idx, drop = FALSE]
  eps = suppressWarnings(max(diag(cov_obs), na.rm = TRUE)) * 1e-8
  if (!is.finite(eps) || eps <= 0) {
    eps = 1e-8
  }
  diag(cov_obs) = diag(cov_obs) + eps
  obs_ids = paste0('obs_', seq_len(nrow(dat)))
  rownames(cov_obs) = obs_ids
  colnames(cov_obs) = obs_ids
  return(list(data = dat, Cov = cov_obs))
}

empty_stat_row = function(output_cols, trait_col, expression_base, fit_mode = 'no_data') {
  out = as.list(rep(NA_real_, length(output_cols)))
  names(out) = output_cols
  out[['trait']] = trait_col
  out[['variable']] = expression_base
  out[['fit_mode']] = fit_mode
  return(as.data.frame(out, stringsAsFactors = FALSE))
}

prepare_epgls_data = function(df_trait_exp, trait_col, expression_base, expression_cols, merge_replicates) {
  if (merge_replicates == 'yes') {
    mean_col = paste0('mean_', expression_base)
    x = if (mean_col %in% colnames(df_trait_exp)) {
      suppressWarnings(as.numeric(df_trait_exp[[mean_col]]))
    } else if (length(expression_cols) > 0) {
      apply(df_trait_exp[, expression_cols, drop = FALSE], 1, function(v) {
        vv = suppressWarnings(as.numeric(v))
        if (all(is.na(vv))) {
          return(NA_real_)
        }
        return(mean(vv, na.rm = TRUE))
      })
    } else {
      rep(NA_real_, nrow(df_trait_exp))
    }
    dat = data.frame(
      gene_id = as.character(df_trait_exp[['gene_id']]),
      species = as.character(df_trait_exp[['species']]),
      y = suppressWarnings(as.numeric(df_trait_exp[[trait_col]])),
      x = x,
      stringsAsFactors = FALSE
    )
  } else {
    if (length(expression_cols) == 0 && (expression_base %in% colnames(df_trait_exp))) {
      expression_cols = c(expression_base)
    }
    if (length(expression_cols) == 0) {
      return(data.frame())
    }
    dat_list = lapply(expression_cols, function(expr_col) {
      data.frame(
        gene_id = as.character(df_trait_exp[['gene_id']]),
        species = as.character(df_trait_exp[['species']]),
        y = suppressWarnings(as.numeric(df_trait_exp[[trait_col]])),
        x = suppressWarnings(as.numeric(df_trait_exp[[expr_col]])),
        stringsAsFactors = FALSE
      )
    })
    dat = do.call(rbind, dat_list)
  }
  dat = dat[
    !is.na(dat[['gene_id']]) & nzchar(dat[['gene_id']]) &
      !is.na(dat[['species']]) & nzchar(dat[['species']]) &
      !is.na(dat[['y']]) & !is.na(dat[['x']]),
    ,
    drop = FALSE
  ]
  rownames(dat) = NULL
  return(dat)
}

resolve_ws_subject = function(ws_subject_arg, merge_replicates) {
  if (ws_subject_arg == 'auto') {
    if (merge_replicates == 'no') {
      return('gene_id')
    }
    return('species')
  }
  return(ws_subject_arg)
}

fit_epgls_one = function(df_trait_exp, gene_tree, trait_col, expression_base, expression_cols, merge_replicates, ws_subject, iter, output_cols) {
  dat = prepare_epgls_data(df_trait_exp, trait_col, expression_base, expression_cols, merge_replicates)
  if (nrow(dat) == 0) {
    return(empty_stat_row(output_cols, trait_col, expression_base, 'no_data'))
  }

  cov_obj = build_observation_cov_from_gene_tree(dat, gene_tree)
  dat = cov_obj[['data']]
  cov_mat = cov_obj[['Cov']]
  if (!(ws_subject %in% colnames(dat))) {
    return(empty_stat_row(output_cols, trait_col, expression_base, 'no_data'))
  }
  dat[['subject_block']] = factor(dat[[ws_subject]])
  if (nrow(dat) < 4 || length(unique(dat[['species']])) < 3 || length(unique(dat[['gene_id']])) < 3 || length(unique(dat[['subject_block']])) < 3 || is.null(cov_mat)) {
    return(empty_stat_row(output_cols, trait_col, expression_base, 'no_data'))
  }

  fit = try(
    RRPP::lm.rrpp.ws(
      y ~ subject_block + x,
      data = dat,
      subjects = 'subject_block',
      Cov = cov_mat,
      iter = iter,
      print.progress = FALSE,
      Parallel = FALSE
    ),
    silent = TRUE
  )

  if (inherits(fit, 'try-error')) {
    return(empty_stat_row(output_cols, trait_col, expression_base, 'epgls_failed'))
  }

  sm = summary(fit)
  tbl = sm[['table']]
  if (is.null(tbl) || nrow(tbl) == 0) {
    return(empty_stat_row(output_cols, trait_col, expression_base, 'epgls_failed'))
  }

  stat_row = if ('x' %in% rownames(tbl)) 'x' else rownames(tbl)[1]
  tt = tbl[stat_row, , drop = FALSE]
  r2 = suppressWarnings(as.numeric(tt[['Rsq']]))
  n = nrow(dat)
  p = 1
  r2adj = NA_real_
  if (!is.na(r2) && n > (p + 1)) {
    r2adj = 1 - ((1 - r2) * (n - 1) / (n - p - 1))
  }

  rss = suppressWarnings(as.numeric(tt[['Residual SS']]))
  rdf = suppressWarnings(as.numeric(tt[['Residual Df']]))
  sigma = if (!is.na(rss) && !is.na(rdf) && rdf > 0) sqrt(rss / rdf) else NA_real_
  fstat = suppressWarnings(as.numeric(tt[['F']]))
  pval = suppressWarnings(as.numeric(tt[['Pr(>F)']]))
  pcc = suppressWarnings(cor(dat[['x']], dat[['y']], method = 'pearson', use = 'complete.obs'))
  fit_mode = if (merge_replicates == 'yes') {
    'epgls_mean'
  } else if (ws_subject == 'gene_id') {
    'epgls_within_copy'
  } else {
    'epgls_within_species'
  }

  out = empty_stat_row(output_cols, trait_col, expression_base, fit_mode)
  out[['R2']] = r2
  out[['R2adj']] = r2adj
  out[['sigma']] = sigma
  out[['Fstat']] = fstat
  out[['pval']] = pval
  out[['PCC']] = pcc
  return(out)
}

tree = ape::read.tree(args[['file_tree']])
trait = read.table(args[['file_trait']], header = TRUE, sep = '\t')
trait[, 'species'] = sub(' ', '_', trait[, 'species'])
exp = read.table(args[['file_exp']], header = TRUE, sep = '\t')

if (args[['merge_replicates']] == 'yes') {
  cat('Expression replicates will be merged before E-PGLS.\n')
  rownames(exp) = exp[, 1]
  exp = exp[, 2:ncol(exp), drop = FALSE]
  exp = merge_replicates(exp, args[['replicate_sep']])
  exp = cbind(data.frame(gene_id = rownames(exp), stringsAsFactors = FALSE), exp)
  rownames(exp) = NULL
} else {
  cat('Expression replicates will not be merged (raw replicate rows retained for RRPP).\n')
}

exp = sort_exp_fn(exp, tree)
expression_bases = get_expression_bases_fn(exp, args[['replicate_sep']])
exp = add_expression_mean_cols(exp, expression_bases)

exp[, 'species'] = leaf_to_species(exp[, 'gene_id'])
exp = sort_exp_fn(exp, tree)
df_trait_exp = merge(exp, trait, by = 'species', all.x = TRUE, sort = FALSE)
trait_cols = colnames(trait[2:length(trait)])

ws_subject = resolve_ws_subject(args[['ws_subject']], args[['merge_replicates']])
cat('E-PGLS covariance mode: observation-level covariance derived directly from gene tree tips (gene_id).\n')
cat('E-PGLS within-subject blocking variable:', ws_subject, '\n')
output_cols = c('R2', 'R2adj', 'sigma', 'Fstat', 'pval', 'logLik', 'AIC', 'BIC', 'PCC', 'p.adj', 'trait', 'variable', 'fit_mode')
res_list = list()

i = 1
for (trait_col in trait_cols) {
  for (expression_base in expression_bases) {
    cat('Working with', trait_col, 'vs', expression_base, '\n')
    expression_cols = colnames(df_trait_exp)[startsWith(colnames(df_trait_exp), expression_base)]
    expression_cols = expression_cols[(expression_cols == expression_base) | (grepl('.*[0-9]$', expression_cols))]
    res_list[[i]] = fit_epgls_one(
      df_trait_exp = df_trait_exp,
      gene_tree = tree,
      trait_col = trait_col,
      expression_base = expression_base,
      expression_cols = expression_cols,
      merge_replicates = args[['merge_replicates']],
      ws_subject = ws_subject,
      iter = iter,
      output_cols = output_cols
    )
    i = i + 1
  }
}

if (length(res_list) > 0) {
  df_stat = do.call(rbind, res_list)
} else {
  df_stat = data.frame(matrix(nrow = 0, ncol = length(output_cols)))
  colnames(df_stat) = output_cols
}

df_stat = finalize_pgls_stats(df_stat)
if (nrow(df_stat) > 0) {
  pvals = suppressWarnings(as.numeric(df_stat[, 'pval']))
  df_stat[, 'p.adj'] = p.adjust(pvals, length(pvals), method = 'fdr')
}
write.table(df_stat, 'gene_tree_EPGLS.tsv', row.names = FALSE, sep = '\t')

plot_df = df_stat
plot_df[, 'trait'] = factor(plot_df[, 'trait'], levels = trait_cols)
plot_df[, 'logp'] = suppressWarnings(-log(plot_df[, 'pval']))

if (sum(is.finite(plot_df[, 'logp'])) > 0) {
  g = ggplot(plot_df, aes(x = trait, y = logp, fill = variable)) +
    geom_bar(position = 'dodge', stat = 'identity') +
    geom_hline(yintercept = -log(0.05), linetype = 'dashed', color = 'red', linewidth = 0.5)
  ggsave('gene_EPGLS.barplot.pdf', g)
} else {
  grDevices::pdf('gene_EPGLS.barplot.pdf')
  plot.new()
  text(0.5, 0.5, 'No valid E-PGLS result to plot')
  dev.off()
}
