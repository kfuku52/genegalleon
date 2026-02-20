build_phenocov_input = function(df_trait_exp, trait_col, expression_cols, expression_base) {
  if (length(expression_cols) == 0) {
    return(list(trait_data = data.frame(), phenocov_list = list()))
  }

  expr_df = df_trait_exp[, expression_cols, drop = FALSE]
  expr_df = as.data.frame(lapply(expr_df, function(x) suppressWarnings(as.numeric(x))))
  expr_mat = as.matrix(expr_df)

  mean_col = paste0('mean_', expression_base)
  if (mean_col %in% colnames(df_trait_exp)) {
    x_mean = suppressWarnings(as.numeric(df_trait_exp[[mean_col]]))
  } else {
    x_mean = apply(expr_mat, 1, function(v) {
      if (all(is.na(v))) {
        return(NA_real_)
      }
      return(mean(v, na.rm = TRUE))
    })
  }
  x_var = apply(expr_mat, 1, function(v) {
    vv = v[!is.na(v)]
    if (length(vv) > 1) {
      return(stats::var(vv))
    }
    return(0)
  })
  n_obs = rowSums(!is.na(expr_mat))
  x_se2 = rep(NA_real_, length(n_obs))
  x_se2[n_obs > 0] = x_var[n_obs > 0] / n_obs[n_obs > 0]

  out = data.frame(species = as.character(df_trait_exp[['species']]), stringsAsFactors = FALSE)
  out[[trait_col]] = suppressWarnings(as.numeric(df_trait_exp[[trait_col]]))
  out[[expression_base]] = x_mean
  is_valid = !is.na(out[['species']]) & nzchar(out[['species']]) & !is.na(out[[trait_col]]) & !is.na(out[[expression_base]])
  if (!any(is_valid)) {
    return(list(trait_data = data.frame(), phenocov_list = list()))
  }
  out = out[is_valid, c('species', trait_col, expression_base), drop = FALSE]
  x_se2 = x_se2[is_valid]

  if (anyDuplicated(out[['species']])) {
    sp_raw = out[['species']]
    trait_df = aggregate(out[, c(trait_col, expression_base), drop = FALSE], by = list(species = sp_raw), function(v) mean(v, na.rm = TRUE))
    se_df = aggregate(data.frame(se2 = x_se2), by = list(species = sp_raw), function(v) {
      vv = v[is.finite(v)]
      if (length(vv) == 0) {
        return(0)
      }
      return(mean(vv))
    })
    out = merge(trait_df, se_df, by = 'species', all.x = TRUE, sort = FALSE)
  } else {
    out[['se2']] = x_se2
  }

  out[['se2']][!is.finite(out[['se2']])] = 0
  out[['se2']][out[['se2']] < 0] = 0

  phenocov_list = vector('list', length = nrow(out))
  names(phenocov_list) = out[['species']]
  for (j in seq_len(nrow(out))) {
    this_se2 = out[['se2']][j]
    m = matrix(c(0, 0, 0, this_se2), nrow = 2, byrow = TRUE)
    rownames(m) = c(trait_col, expression_base)
    colnames(m) = c(trait_col, expression_base)
    phenocov_list[[j]] = m
  }

  out = out[, c('species', trait_col, expression_base), drop = FALSE]
  return(list(trait_data = out, phenocov_list = phenocov_list))
}

run_phylopars_regression = function(df_trait_exp, tree, trait_cols, expression_bases, output_cols, include_foreground_lineage = FALSE, verbose_working = FALSE, use_phenocov = FALSE) {
  fg_nums = NULL
  if (include_foreground_lineage) {
    fg_nums = rkftools::count_foreground_lineage(tree, df_trait_exp[, c('species', trait_cols), drop = FALSE])
  }

  nrow_stat = length(trait_cols) * length(expression_bases)
  df_stat = data.frame(matrix(NA, nrow_stat, length(output_cols)))
  colnames(df_stat) = output_cols

  i = 1
  for (trait_col in trait_cols) {
    for (expression_base in expression_bases) {
      if (verbose_working) {
        cat('Working with', trait_col, 'vs', expression_base, '\n')
      }
      expression_cols = colnames(df_trait_exp)[startsWith(colnames(df_trait_exp), expression_base)]
      expression_cols = expression_cols[(expression_cols == expression_base) | (grepl('.*[0-9]$', expression_cols))]
      explanatory_variables = paste(expression_cols, collapse = ' + ')
      explained_variable = trait_col
      formula_string = paste(explained_variable, '~', explanatory_variables)
      df_model = df_trait_exp
      run_mode = 'wide'
      fit_mode = run_mode
      pheno_error = TRUE
      pheno_correlated = TRUE
      phenocov_list = list()
      if (use_phenocov) {
        pheno_input = build_phenocov_input(df_trait_exp, trait_col, expression_cols, expression_base)
        df_model = pheno_input[['trait_data']]
        phenocov_list = pheno_input[['phenocov_list']]
        formula_string = paste(explained_variable, '~', expression_base)
        run_mode = 'phenocov'
        fit_mode = run_mode
        pheno_error = FALSE
        pheno_correlated = FALSE
      } else {
        if (length(expression_cols) > 0) {
          for (expr_col in expression_cols) {
            df_model[[expr_col]] = suppressWarnings(as.numeric(df_model[[expr_col]]))
          }
        }
        df_model[[trait_col]] = suppressWarnings(as.numeric(df_model[[trait_col]]))
      }
      cat('PGLS formula:', formula_string, '\n')

      has_signal = TRUE
      if (use_phenocov) {
        has_signal = (expression_base %in% colnames(df_model)) && any(!is.na(df_model[[expression_base]]))
      } else {
        has_signal = any(rowSums(!is.na(df_model[, expression_cols, drop = FALSE])) > 0)
      }
      if (nrow(df_model) == 0 || !has_signal) {
        cat('no valid observations for formula:', formula_string, '\n')
        df_stat[i, 'trait'] = trait_col
        df_stat[i, 'variable'] = expression_base
        if ('fit_mode' %in% colnames(df_stat)) {
          df_stat[i, 'fit_mode'] = 'no_data'
        }
        if (include_foreground_lineage) {
          df_stat[i, 'num_foreground_lineage'] = fg_nums[[trait_col]]
        }
        i = i + 1
        next
      }

      out_phylopars = try(phylopars.lm(
        as.formula(formula_string),
        trait_data = df_model,
        tree = tree,
        model = "BM",
        pheno_error = pheno_error,
        phylo_correlated = TRUE,
        pheno_correlated = pheno_correlated,
        phenocov_list = phenocov_list
      ))
      if (class(out_phylopars) == "try-error" && use_phenocov) {
        mean_col = paste0('mean_', expression_base)
        if (mean_col %in% colnames(df_trait_exp)) {
          cat('phenocov model failed; fallback to mean-based model:', formula_string, '\n')
          df_model = df_trait_exp[, c('species', trait_col, mean_col), drop = FALSE]
          colnames(df_model)[3] = expression_base
          df_model[[trait_col]] = suppressWarnings(as.numeric(df_model[[trait_col]]))
          df_model[[expression_base]] = suppressWarnings(as.numeric(df_model[[expression_base]]))
          df_model = df_model[!is.na(df_model[['species']]) & !is.na(df_model[[trait_col]]) & !is.na(df_model[[expression_base]]), , drop = FALSE]
          formula_string = paste(explained_variable, '~', expression_base)
          run_mode = 'mean_fallback'
          fit_mode = run_mode
          if (nrow(df_model) > 0) {
            out_phylopars = try(phylopars.lm(
              as.formula(formula_string),
              trait_data = df_model,
              tree = tree,
              model = "BM",
              pheno_error = TRUE,
              phylo_correlated = TRUE,
              pheno_correlated = TRUE
            ))
          }
        }
      }

      if (class(out_phylopars) != "try-error") {
        cat('PGLS fit mode:', fit_mode, '\n')
        for (stat in output_cols[1:6]) {
          df_stat[i, stat] = out_phylopars[stat]
        }
        df_stat[i, 'AIC'] = AIC(out_phylopars)
        df_stat[i, 'BIC'] = BIC(out_phylopars)
        mean_col = paste0('mean_', expression_base)
        if (mean_col %in% colnames(df_trait_exp)) {
          pcc = suppressWarnings(cor(df_trait_exp[[mean_col]], df_trait_exp[[trait_col]], method = 'pearson', use = 'complete.obs'))
          df_stat[i, 'PCC'] = pcc
        }
      } else {
        if (fit_mode == 'mean_fallback') {
          fit_mode = 'mean_fallback_failed'
        } else if (fit_mode == 'phenocov') {
          fit_mode = 'phenocov_failed'
        } else if (fit_mode == 'wide') {
          fit_mode = 'wide_failed'
        }
        cat('error during phylopars process:', formula_string, '\n')
      }

      df_stat[i, 'trait'] = trait_col
      df_stat[i, 'variable'] = expression_base
      if ('fit_mode' %in% colnames(df_stat)) {
        df_stat[i, 'fit_mode'] = fit_mode
      }
      if (include_foreground_lineage) {
        df_stat[i, 'num_foreground_lineage'] = fg_nums[[trait_col]]
      }
      i = i + 1
    }
  }

  cat('adjusted pvalues are calculated:\n')
  cat('method = fdr\n')
  cat('length =', length(df_stat[, 'pval']), '\n')
  df_stat[, 'p.adj'] = p.adjust(df_stat[, 'pval'], length(df_stat[, 'pval']), method = 'fdr')
  return(df_stat)
}

add_expression_mean_cols = function(exp, expression_bases) {
  for (col in expression_bases) {
    is_col = grepl(col, colnames(exp))
    mean_col = paste('mean_', col, sep = '')
    if (sum(is_col) > 1) {
      exp[, mean_col] = apply(exp[, is_col], 1, function(x) { mean(x, na.rm = TRUE) })
    } else {
      exp[, mean_col] = exp[, is_col]
    }
  }
  return(exp)
}

finalize_pgls_stats = function(df_stat) {
  if (nrow(df_stat) == 0) {
    return(df_stat)
  }
  if ('variable' %in% colnames(df_stat)) {
    df_stat = df_stat[order(df_stat[, 'variable']), , drop = FALSE]
  }
  df_stat = df_stat[apply(df_stat, 1, function(x) !all(is.na(x))), , drop = FALSE]
  rownames(df_stat) = NULL
  return(df_stat)
}
