add_ortholog_column <- function(g,
                                args,
                                tree,
                                ortholog_prefix,
                                path_ortho_nwk = 'no',
                                show_all_ortholog = FALSE) {
  cat(as.character(Sys.time()), 'Adding ortholog column.\n')
  
  font_size <- args[['font_size']]
  df_tip <- get_df_tip(g[['tree']])
  
  # Decide which tree to use for ortholog inference
  if (path_ortho_nwk != 'no') {
    if (file.exists(path_ortho_nwk)) {
      cat('The user-provided tree will be used for ortholog inference.\n')
      tree2 <- ape::read.tree(path_ortho_nwk)
    } else {
      path_ortho_nwk <- 'no'
    }
  }
  if (path_ortho_nwk == 'no') {
    cat('The input tree will be used for ortholog inference.\n')
    tree2 <- tree
  }
  
  # Identify which tips are queries vs. subjects
  is_query <- tree2[['tip.label']] %in% tree[['tip.label']]
  queries  <- tree2[['tip.label']][is_query]
  subjects <- tree2[['tip.label']][startsWith(tree2[['tip.label']], ortholog_prefix)]
  
  # If no subjects found, bail out
  if (length(subjects) == 0) {
    cat(ortholog_prefix, ' not found in the tree specified for ortholog inference.\n')
    return(g)
  }
  
  # Find nearest orthologs for each query
  # We use lapply to build the rows in parallel, then rbind them.
  cat('Finding nearest ortholog tips...\n')
  mrca_matrix_2 <- ape::mrca(tree2)
  nearest_info_list <- lapply(queries, function(q) {
    tmp <- rkftools::get_nearest_tips(phy = tree2, query = q, subjects = subjects, mrca_matrix = mrca_matrix_2)
    data.frame(
      query    = q,
      nearests = paste(sort(tmp[['nearests']]), collapse = '\n'),
      mrca     = tmp[['mrca']],
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, nearest_info_list)
  
  # If a separate tree file was provided, we need to recalc MRCA in the "main" tree
  if (path_ortho_nwk != 'no') {
    mrca_matrix_main <- ape::mrca(tree)
    uniq_vals <- unique(df[['nearests']])
    for (val in uniq_vals) {
      qs <- df[['query']][df[['nearests']] == val]
      
      # Subset the matrix to only these queries
      sub_mrcas <- mrca_matrix_main[qs, qs]
      # Flatten into a vector (diagonal repeated but usually not a big deal)
      sub_mrcas_vec <- as.vector(sub_mrcas)
      
      new_mrca <- sub_mrcas_vec[1]
      for (nm in sub_mrcas_vec) {
        tmp_mrcas <- rkftools::get_ancestor_num(tree, nm)
        if (any(sub_mrcas_vec %in% tmp_mrcas)) {
          # The first match we find becomes the new MRCA
          new_mrca <- sub_mrcas_vec[sub_mrcas_vec %in% tmp_mrcas][1]
          break
        }
      }
      df[df[['nearests']] == val, 'mrca'] <- new_mrca
    }
  } else {
    cat(paste0('Clade orthologs will be inferred from ', args[['stat_branch']], '\n'))
  }
  
  # Adjust queries in case they accidentally matched themselves
  df <- correct_self_hit_query(df, tree2, mrca_matrix_2, ortholog_prefix)
  
  # Merge the 'y' position from the g[['tree']] data so we know where to plot
  df <- merge(df, g[['tree']][['data']][, c('label', 'y')], 
              by.x = 'query', by.y = 'label', sort = FALSE)
  
  # Build a data frame describing the clade structure (ymin, ymax, etc.)
  df_clade <- get_df_clade(df, ortholog_prefix)
  if (!show_all_ortholog) {
    df_clade <- abbreviate_clade_labels(df_clade, num_line_offset = 1)
  }
  df_clade <- split_nested_clade(df_clade)
  
  # Axis label
  xlabel <- ortholog_prefix
  xlabel <- sub('_$', '', xlabel)
  xlabel <- sub('_', ' ', xlabel)
  xlabel <- paste0(xlabel, '\nclosest gene')
  
  gname <- paste0('ortholog,', ortholog_prefix)
  
  # Build the ggplot layer for orthologs
  g[[gname]] <- ggplot() +
    # "Invisible" points to define y-coordinates
    geom_point(data = df_tip, mapping = aes(y = label), x = 0, size = 0, colour = rgb(1,1,1,0)) +
    xlim(0, 1) +
    labs(x = xlabel, y = NULL) +
    theme_minimal(base_size = font_size) +
    coord_cartesian(clip = "off") +
    theme(
      axis.text       = element_blank(),
      axis.ticks      = element_blank(),
      panel.grid      = element_blank(),
      plot.margin     = unit(args[['margins']], "cm")
    ) +
    # Segment lines showing each clade's vertical span
    geom_segment(
      aes(y = ymin, yend = ymax, x = 0, xend = 0),
      data = df_clade,
      linewidth = 0.25,
      color='black'
    ) +
    # Text labels for each clade's nearest ortholog(s)
    geom_text(
      aes(x = 0.05, y = ymid, label = nearests, hjust = 0, vjust = 0.5),
      data = df_clade,
      size = args[['font_size']] * args[['font_size_factor']],
      lineheight = 0.8
    )
  
  return(g)
}

get_line_coordinate = function(g, branch_ids, jitter=FALSE, pairwise=FALSE) {
    dat = g[['data']]
    if (!all(c('branch_id', 'branch', 'y') %in% colnames(dat))) {
        return(data.frame())
    }
    x_map = setNames(dat[['branch']], dat[['branch_id']])
    y_map = setNames(dat[['y']], dat[['branch_id']])
    if (pairwise) {
        if (length(branch_ids) < 2) {
            return(data.frame())
        }
        nl_combinations = t(combn(branch_ids, 2))
        out = data.frame(
            nl_start = nl_combinations[, 1],
            nl_end = nl_combinations[, 2],
            stringsAsFactors = FALSE
        )
        out[['x_start']] = unname(x_map[as.character(out[['nl_start']])])
        out[['y_start']] = unname(y_map[as.character(out[['nl_start']])])
        out[['x_end']] = unname(x_map[as.character(out[['nl_end']])])
        out[['y_end']] = unname(y_map[as.character(out[['nl_end']])])
        out = out[is.finite(out[['x_start']]) & is.finite(out[['x_end']]) &
                  is.finite(out[['y_start']]) & is.finite(out[['y_end']]), , drop = FALSE]
        return(out)
    }
    idx = match(as.character(branch_ids), as.character(dat[['branch_id']]))
    idx = idx[!is.na(idx)]
    if (length(idx) < 2) {
        return(data.frame())
    }
    out = dat[idx, c('branch_id', 'branch', 'y'), drop = FALSE]
    colnames(out) = sub('^y$', 'y_start', colnames(out))
    out = out[order(out[['y_start']], decreasing = TRUE), , drop = FALSE]
    rownames(out) = NULL
    out[,'x_start'] = out[['branch']]
    if (jitter) {
        out[,'x_start'] = jitter(out[['x_start']], amount = jitter)
    }
    n = nrow(out)
    out[,'x_end'] = c(out[['x_start']][-1], NA)
    out[,'y_end'] = c(out[['y_start']][-1], NA)
    if (n > 1) {
        out = out[seq_len(n - 1), , drop = FALSE]
    } else {
        out = out[0, , drop = FALSE]
    }
    rownames(out) = NULL
    return(out)
}

overlay_convergence = function(g, top_percent_to_show, max_num_to_show, is_target_only, arity_range, cb_path, cutoff_stat, node_label='highest_arity') {
  if (max_num_to_show == 0) {
    cat('Maximum number of convergence connection is set to 0. Skipping.\n')
    return(g)
  }
  min_arity = as.integer(strsplit(arity_range, '-')[[1]][1])
  max_arity = as.integer(strsplit(arity_range, '-')[[1]][2])
  df_label_coords = list()
  df_line_coords = list()
  reserved_combinations = list()
  for (arity in max_arity:min_arity) {
    cat(as.character(Sys.time()), 'Adding convergence connection for arity (K) =', arity, '\n')
    arity_txt = paste0('arity_', arity)
    cb_file = gsub('ARITY', arity, cb_path)
    if (!file.exists(cb_file)) {
      cat('Skipping. CSUBST cb file was not found for arity (K) =', arity, '\n')
      next
    }
    cb = read.table(cb_file, sep = '\t', header = TRUE, stringsAsFactors = FALSE, check.name = FALSE)
    bid_cols = colnames(cb)[grep('branch_id_', colnames(cb))]
    fg_cols = colnames(cb)[grep('is_fg', colnames(cb))]
    if (is_target_only) {
      is_fg = apply(cb[, fg_cols, drop=FALSE], 1, function(x) { any(x == 'Y') })
      cb_target = cb[is_fg,]
    } else {
      cb_target = cb
    }
    max_num_show = floor(nrow(cb) * top_percent_to_show * 0.01)
    max_num_show = max(max_num_show, 1)
    cat(paste0(top_percent_to_show, '% (', max_num_show, '+tiers/', nrow(cb), ') of top protein convergence will be analyzed.\n'))
    cutoff_stat_list = list()
    cutoff_stat_names = c()
    if (!is.null(cutoff_stat) && nzchar(cutoff_stat)) {
      for (split in strsplit(strsplit(cutoff_stat, '\\|')[[1]], ',')) {
        if (length(split) < 2) {
          next
        }
        stat_name = split[[1]]
        stat_value = as.numeric(split[[2]])
        if (!(stat_name %in% colnames(cb_target)) || is.na(stat_value)) {
          cat('Skipping invalid cutoff stat: ', paste(split, collapse=','), '\n')
          next
        }
        cutoff_stat_list[[stat_name]] = stat_value
        cutoff_stat_names = c(cutoff_stat_names, stat_name)
      }
    }
    if (length(cutoff_stat_names) > 0) {
      cb_target = cb_target[rev(do.call(order, cb_target[,cutoff_stat_names, drop=FALSE])),] # sort by the first cutoff_stat
      for (cutoff_stat_name in cutoff_stat_names) {
        cb_target = cb_target[(cb_target[[cutoff_stat_name]] >= cutoff_stat_list[[cutoff_stat_name]]),]
      }
    }
    rownames(cb_target) = NULL
    cat(nrow(cb_target), 'branch combinations satisfied cutoff stats:', cutoff_stat, '\n')
    if (nrow(cb_target) > max_num_to_show) {
      cat('For visualization, top', max_num_to_show, 'protein convergence will be visualized.\n')
      cb_target = cb_target[1:max_num_to_show,]
    }
    counts_list = list()
    counts_idx = 1
    if (nrow(cb_target) > 0) {
      line_coords_list = list()
      line_coords_idx = 1
      for (i in 1:nrow(cb_target)) {
        branch_ids = as.vector(unlist(cb_target[i, bid_cols]))
        skip_flag = FALSE
        for (rc in reserved_combinations) {
          if (all(branch_ids %in% rc)) {
            skip_flag = TRUE
            break
          }
        }
        if (skip_flag) {
          next
        } else {
          reserved_combinations[[length(reserved_combinations)+1]] = branch_ids
        }
        line_coords = get_line_coordinate(g, branch_ids, jitter=0.01, pairwise=FALSE)
        if (nrow(line_coords) == 0) {
          next
        }
        line_coords[1, 'OCNany2spe'] = round(cb_target[i, 'OCNany2spe'], digits = 1)
        line_coords[1, 'omegaCany2spe'] = round(cb_target[i, 'omegaCany2spe'], digits = 1)
        g = g + geom_curve(aes(x = x_start, y = y_start, xend = x_end, yend = y_end, alpha = 0.1),
                           linewidth = 0.1, curvature = jitter(0, amount = 0.05), colour = 'firebrick',
                           data = line_coords, show.legend = FALSE)
        g$layers = c(g$layers[[length(g$layers)]], g$layers[-length(g$layers)]) # Move it to the bottommost layer
        counts_list[[counts_idx]] = as.numeric(unname(branch_ids))
        counts_idx = counts_idx + 1
        line_coords_list[[line_coords_idx]] = line_coords
        line_coords_idx = line_coords_idx + 1
      }
      if (length(line_coords_list) == 0) {
        df_line_coords[[arity_txt]] = data.frame()
      } else {
        df_line_coords[[arity_txt]] = do.call(rbind, line_coords_list)
      }
      labels = paste0(df_line_coords[[arity_txt]][['OCNany2spe']], '/', df_line_coords[[arity_txt]][['omegaCany2spe']])
      labels = ifelse(is.na(df_line_coords[[arity_txt]][['OCNany2spe']]), '', labels)
      df_label_coords[[arity_txt]] = data.frame(
        'x' = df_line_coords[[arity_txt]][['x_start']] + ((df_line_coords[[arity_txt]][['x_end']] - df_line_coords[[arity_txt]][['x_start']]) / 2),
        'y' = df_line_coords[[arity_txt]][['y_start']] + ((df_line_coords[[arity_txt]][['y_end']] - df_line_coords[[arity_txt]][['y_start']]) / 2),
        label = labels
      )
      g = g + ggrepel::geom_text_repel(mapping = aes(x = x, y = y, label = label), data = df_label_coords[[arity_txt]],
                                       color = 'firebrick', size = 2, hjust = 0.5, vjust = 0.5, force = 0.02)
    }
    count_col = paste0('count_', arity)
    if (length(counts_list) > 0) {
      counts = table(unlist(counts_list))
      circle_coords = data.frame(branch_id = as.numeric(names(counts)))
      circle_coords[,count_col] = c(unname(counts))
      g[['data']] = merge(g[['data']], circle_coords, by = 'branch_id', sort = FALSE, all.x = TRUE)
    } else {
      g[['data']][,count_col] = 0
    }
    g[['data']][(is.na(g[['data']][[count_col]])), count_col] = 0
  }
  g[['data']][,'count'] = 0
  for (count_col in colnames(g[['data']])[startsWith(colnames(g[['data']]), 'count_')]) {
    count_col_arity = as.integer(sub('^count_', '', count_col))
    if (node_label=='pairwise_count') {
      g[['data']][,'count'] = g[['data']][['count']] + (g[['data']][[count_col]] * (count_col_arity-1))
    } else if (node_label=='highest_arity') {
      is_detected = as.logical(g[['data']][[count_col]])
      if (any(is_detected)) {
        g[['data']][is_detected,'count'] = sapply(unlist(g[['data']][is_detected,'count']), function(x){max(x, count_col_arity)})
      }
    }
  }
  g[['data']][, 'show_circle'] = as.logical(g[['data']][['count']])
  g = g + geom_point2(aes(x = branch, y = y, subset = show_circle), color = 'firebrick', shape = 16, size = 2.5, show.legend = FALSE)
  g$layers = c(g$layers[[length(g$layers)]], g$layers[-length(g$layers)]) # Move it to the bottommost layer
  g = g + geom_text2(aes(x = branch, y = y, label = count, subset = show_circle), color = 'white', size = 2, show.legend = FALSE)
  g = g + guides(alpha='none')
  return(g)
}

overlay_branch_combination = function(g, branch_combination, jitter=FALSE) {
    if (length(branch_combination)==0) {
        cat('Branch combination is empty. Skipping.\n')
        return(g)
    }
    cat(as.character(Sys.time()), 'Adding specified branch combination:', branch_combination, '\n')
    line_coords = get_line_coordinate(g, branch_combination, jitter=jitter, pairwise=FALSE)
    g = g + geom_curve(aes(x = x_start, y = y_start, xend = x_end, yend = y_end, alpha = 0.1),
                       linewidth = 1.0, curvature = jitter(0, amount = 0.05), colour = 'firebrick',
                       data = line_coords, show.legend = FALSE)
    g$layers = c(g$layers[[length(g$layers)]], g$layers[-length(g$layers)]) # Move it to the bottommost layer
    g[['data']][, 'show_circle'] = (g[['data']][['branch_id']] %in% branch_combination)
    g = g + geom_point2(aes(x = branch, y = y, subset = show_circle), color = 'firebrick', shape = 16, size = 2.5, show.legend = FALSE)
    g$layers = c(g$layers[[length(g$layers)]], g$layers[-length(g$layers)]) # Move it to the bottommost layer
    return(g)
}

propagate_tiplab_colors_to_internal_branches = function(gtree, tree) {
    dat = gtree[['data']]
    required_cols = c('node', 'parent', 'isTip', 'branch_color', 'x')
    if (!all(required_cols %in% colnames(dat))) {
        return(gtree)
    }
    idx_by_node = setNames(seq_len(nrow(dat)), as.character(dat[['node']]))
    children_by_parent = split(as.character(dat[['node']]), as.character(dat[['parent']]), drop = TRUE)
    is_black = as.character(dat[['branch_color']]) %in% c('#000000', 'black')
    desc_color = rep(NA_character_, nrow(dat))
    tip_idx = which(dat[['isTip']])
    desc_color[tip_idx] = as.character(dat[['branch_color']][tip_idx])
    desc_color[tip_idx][is_black[tip_idx]] = NA_character_
    ord = order(dat[['x']], decreasing = TRUE, na.last = TRUE)
    for (ri in ord) {
        if (isTRUE(dat[['isTip']][ri])) {
            next
        }
        node = as.character(dat[['node']][ri])
        child_nodes = children_by_parent[[node]]
        if (is.null(child_nodes) || length(child_nodes) == 0) {
            next
        }
        child_idx = idx_by_node[child_nodes]
        child_idx = child_idx[!is.na(child_idx)]
        if (length(child_idx) == 0) {
            next
        }
        colors = unique(desc_color[child_idx])
        colors = colors[!is.na(colors) & nzchar(colors) & !(colors %in% c('#000000', 'black'))]
        if (length(colors) == 1) {
            desc_color[ri] = colors[1]
        }
    }
    can_update = (!dat[['isTip']]) & is_black & !is.na(desc_color)
    dat[['branch_color']][can_update] = desc_color[can_update]
    gtree[['data']] = dat
    return(gtree)
}

parse_optional_numeric = function(value, default_value) {
    if (is.null(value)) {
        return(default_value)
    }
    parsed = suppressWarnings(as.numeric(value))
    if (!is.finite(parsed)) {
        return(default_value)
    }
    return(parsed)
}

get_long_branch_display_settings = function(args) {
    mode = 'auto'
    if ('long_branch_display' %in% names(args)) {
        mode = tolower(as.character(args[['long_branch_display']]))
    }
    if (!(mode %in% c('auto', 'no'))) {
        mode = 'auto'
    }
    settings = list(
        mode = mode,
        ref_quantile = parse_optional_numeric(args[['long_branch_ref_quantile']], 0.95),
        detect_ratio = parse_optional_numeric(args[['long_branch_detect_ratio']], 5),
        cap_ratio = parse_optional_numeric(args[['long_branch_cap_ratio']], 2.5),
        tail_shrink = parse_optional_numeric(args[['long_branch_tail_shrink']], 0.02),
        max_fraction = parse_optional_numeric(args[['long_branch_max_fraction']], 0.1)
    )
    settings[['ref_quantile']] = max(0.5, min(0.99, settings[['ref_quantile']]))
    settings[['detect_ratio']] = max(1, settings[['detect_ratio']])
    settings[['cap_ratio']] = max(1, settings[['cap_ratio']])
    settings[['tail_shrink']] = max(0, min(1, settings[['tail_shrink']]))
    settings[['max_fraction']] = max(0.01, min(1, settings[['max_fraction']]))
    return(settings)
}

compress_long_branches_for_display = function(tree, args) {
    out = list(
        tree = tree,
        compressed = FALSE,
        note = '',
        all_edge = data.frame(
            node = numeric(0),
            original_edge_length = numeric(0),
            display_edge_length = numeric(0),
            is_long_branch_compressed = logical(0),
            stringsAsFactors = FALSE
        ),
        summary = data.frame(
            node = numeric(0),
            original_edge_length = numeric(0),
            display_edge_length = numeric(0),
            stringsAsFactors = FALSE
        )
    )
    if (is.null(tree[['edge']]) || is.null(tree[['edge.length']])) {
        return(out)
    }
    edge_length = suppressWarnings(as.numeric(tree[['edge.length']]))
    if (length(edge_length) != nrow(tree[['edge']])) {
        return(out)
    }
    out[['all_edge']] = data.frame(
        node = tree[['edge']][,2],
        original_edge_length = edge_length,
        display_edge_length = edge_length,
        is_long_branch_compressed = FALSE,
        stringsAsFactors = FALSE
    )
    settings = get_long_branch_display_settings(args)
    if (settings[['mode']] == 'no') {
        return(out)
    }
    finite_positive = is.finite(edge_length) & (edge_length > 0)
    if (sum(finite_positive) < 10) {
        return(out)
    }
    ref_value = suppressWarnings(as.numeric(stats::quantile(
        edge_length[finite_positive],
        probs = settings[['ref_quantile']],
        names = FALSE,
        na.rm = TRUE,
        type = 7
    )))
    if (!is.finite(ref_value) || ref_value <= 0) {
        return(out)
    }
    detect_threshold = ref_value * settings[['detect_ratio']]
    if (!is.finite(detect_threshold) || detect_threshold <= 0) {
        return(out)
    }
    long_idx = which(finite_positive & (edge_length > detect_threshold))
    if (length(long_idx) == 0) {
        return(out)
    }
    long_fraction = length(long_idx) / sum(finite_positive)
    if (long_fraction > settings[['max_fraction']]) {
        cat(
            'Long-branch display compression skipped: too many candidates.',
            'num_long =', length(long_idx),
            'num_edge =', sum(finite_positive),
            'ratio =', signif(long_fraction, 4), '\n'
        )
        return(out)
    }
    cap = ref_value * settings[['cap_ratio']]
    if (!is.finite(cap) || cap <= 0) {
        cap = detect_threshold
    }
    display_edge_length = edge_length
    display_edge_length[long_idx] = cap + (pmax(edge_length[long_idx] - cap, 0) * settings[['tail_shrink']])
    changed_idx = which(is.finite(edge_length) & is.finite(display_edge_length) & (abs(display_edge_length - edge_length) > .Machine$double.eps))
    if (length(changed_idx) == 0) {
        return(out)
    }
    out[['tree']][['edge.length']] = display_edge_length
    out[['compressed']] = TRUE
    out[['all_edge']][['display_edge_length']] = display_edge_length
    out[['all_edge']][['is_long_branch_compressed']][changed_idx] = TRUE
    out[['summary']] = out[['all_edge']][out[['all_edge']][['is_long_branch_compressed']], c('node', 'original_edge_length', 'display_edge_length'), drop = FALSE]
    max_original = max(out[['summary']][['original_edge_length']], na.rm = TRUE)
    max_display = max(out[['summary']][['display_edge_length']], na.rm = TRUE)
    out[['note']] = ''
    cat(
        'Long-branch display compression applied.',
        'num_long =', nrow(out[['summary']]),
        'detect_threshold =', signif(detect_threshold, 4),
        'max_original =', signif(max_original, 4),
        'max_display =', signif(max_display, 4), '\n'
    )
    return(out)
}

add_long_branch_compression_marks = function(g, args, compression_info) {
    if (is.null(compression_info) || !isTRUE(compression_info[['compressed']])) {
        return(g)
    }
    if (!all(c('is_long_branch_compressed', 'node', 'parent', 'x', 'y', 'branch_thickness', 'branch_color') %in% colnames(g[['data']]))) {
        return(g)
    }
    is_mark = (
        as.logical(g[['data']][['is_long_branch_compressed']]) &
        is.finite(g[['data']][['x']]) &
        is.finite(g[['data']][['y']])
    )
    is_mark[is.na(is_mark)] = FALSE
    df_mark = g[['data']][is_mark, c('node', 'parent', 'x', 'y', 'branch_thickness', 'branch_color'), drop = FALSE]
    if (nrow(df_mark) == 0) {
        return(g)
    }
    node_idx = setNames(seq_len(nrow(g[['data']])), as.character(g[['data']][['node']]))
    parent_idx = node_idx[as.character(df_mark[['parent']])]
    parent_x = rep(NA_real_, nrow(df_mark))
    has_parent = !is.na(parent_idx)
    parent_x[has_parent] = suppressWarnings(as.numeric(g[['data']][['x']][parent_idx[has_parent]]))
    df_mark[['parent_x']] = parent_x
    # If parent x is unavailable (e.g., root edge), skip those rows.
    is_valid_parent = is.finite(df_mark[['parent_x']])
    df_mark = df_mark[is_valid_parent, , drop = FALSE]
    if (nrow(df_mark) == 0) {
        return(g)
    }
    df_mark[['x_start']] = pmin(df_mark[['parent_x']], df_mark[['x']])
    df_mark[['x_end']] = pmax(df_mark[['parent_x']], df_mark[['x']])
    df_mark[['edge_len']] = df_mark[['x_end']] - df_mark[['x_start']]
    df_mark = df_mark[is.finite(df_mark[['edge_len']]) & (df_mark[['edge_len']] > 0), , drop = FALSE]
    if (nrow(df_mark) == 0) {
        return(g)
    }
    df_mark[['left_end']] = df_mark[['x_start']] + (df_mark[['edge_len']] / 3)
    df_mark[['mid_start']] = df_mark[['left_end']]
    df_mark[['mid_end']] = df_mark[['x_start']] + (2 * df_mark[['edge_len']] / 3)
    df_mark[['right_start']] = df_mark[['mid_end']]
    df_mark[['branch_color']] = as.character(df_mark[['branch_color']])
    bad_color = is.na(df_mark[['branch_color']]) | !nzchar(df_mark[['branch_color']])
    df_mark[['branch_color']][bad_color] = 'black'
    df_mark = df_mark[df_mark[['mid_end']] > df_mark[['mid_start']], , drop = FALSE]
    if (nrow(df_mark) == 0) {
        return(g)
    }
    erase_width = pmax(suppressWarnings(as.numeric(df_mark[['branch_thickness']])), 0.25) * 1.5
    branch_width = pmax(suppressWarnings(as.numeric(df_mark[['branch_thickness']])), 0.25)
    g = g + geom_segment(
        data = df_mark,
        mapping = aes(x = x_start, xend = x_end, y = y, yend = y),
        inherit.aes = FALSE,
        color = 'white',
        linewidth = erase_width,
        lineend = 'butt'
    )
    g = g + geom_segment(
        data = df_mark,
        mapping = aes(x = x_start, xend = left_end, y = y, yend = y),
        inherit.aes = FALSE,
        color = df_mark[['branch_color']],
        linewidth = branch_width,
        lineend = 'butt'
    )
    g = g + geom_segment(
        data = df_mark,
        mapping = aes(x = mid_start, xend = mid_end, y = y, yend = y),
        inherit.aes = FALSE,
        color = df_mark[['branch_color']],
        linewidth = branch_width,
        linetype = 'dotted',
        lineend = 'butt'
    )
    g = g + geom_segment(
        data = df_mark,
        mapping = aes(x = right_start, xend = x_end, y = y, yend = y),
        inherit.aes = FALSE,
        color = df_mark[['branch_color']],
        linewidth = branch_width,
        lineend = 'butt'
    )
    return(g)
}

add_tree_column = function(g, args, b, dist_col, nodelabel_col, branch_color, orientation, 
                           pie_chart_value_transformation, path_species_color_table) {
    cat(as.character(Sys.time()), 'Adding tree column.\n')
    stopifnot(dist_col %in% colnames(b))
    tree = table2phylo(df=b, name_col='node_name', dist_col=dist_col)
    if (dist_col=='bl_dated') {
        try_out = try(rkftools::force_ultrametric(tree, stop_if_larger_change=0.1), silent = TRUE)
        if (class(try_out)!="try-error") {
            tree = try_out
        }
    }
    compression_info = compress_long_branches_for_display(tree, args)
    tree_display = compression_info[['tree']]
    gname = paste('tree', dist_col, nodelabel_col, branch_color, orientation, sep=',')
    g[[gname]] = ggtree(tree_display, size=0, layout='rectangular') + theme_void()
    g[[gname]][['data']] = data.frame(g[[gname]][['data']], stringsAsFactors=FALSE)
    g[[gname]] = append_branch_stats(g[[gname]], b)
    if (nrow(compression_info[['all_edge']]) > 0) {
        idx = match(g[[gname]][['data']][['node']], compression_info[['all_edge']][['node']])
        g[[gname]][['data']][['original_edge_length']] = compression_info[['all_edge']][['original_edge_length']][idx]
        g[[gname]][['data']][['display_edge_length']] = compression_info[['all_edge']][['display_edge_length']][idx]
        compressed_flag = rep(FALSE, nrow(g[[gname]][['data']]))
        has_idx = !is.na(idx)
        compressed_flag[has_idx] = compression_info[['all_edge']][['is_long_branch_compressed']][idx[has_idx]]
        g[[gname]][['data']][['is_long_branch_compressed']] = compressed_flag
    }
    g[[gname]] = append_branch_tiplab_colors(g[[gname]], branch_color, path_species_color_table)
    g[[gname]] = propagate_tiplab_colors_to_internal_branches(g[[gname]], tree_display)
    if (orientation=='R') {
        # reverse x-axis
        x = g[[gname]][['data']][['x']]
        br = g[[gname]][['data']][['branch']]
        g[[gname]][['data']][['x']] = max(x) - x
        g[[gname]][['data']][['branch']] = max(br) - br
    }
    if (grepl('_regime$', branch_color)) {
        if (branch_color %in% colnames(g[[gname]][['data']])) {
            g[[gname]] = add_pie_charts(g[[gname]], args, tree, branch_color, pie_chart_value_transformation)
        } else {
            cat(branch_color, 'column not found in --stat_branch\n')
        }
    }
    g[[gname]] = g[[gname]] + geom_tree(
        color=g[[gname]][['data']][['branch_color']],
        size=g[[gname]][['data']][['branch_thickness']]
    )
    g[[gname]] = add_node_points(g[[gname]], args)
    g[[gname]] = add_node_labels(g[[gname]], args, nodelabel_col)
    g[[gname]] = g[[gname]] + 
        coord_cartesian(clip = "off") + 
        theme(
            panel.border=element_blank(),
            panel.background=element_blank(),
            legend.background = element_blank(),
            #legend.key = element_blank(),
            plot.margin=unit(args[['margins']], "cm")
        )
    g[[gname]] = add_scale_bar(g[[gname]], args, tree_display)
    if (args[['show_branch_id']]=='yes') {
        g[[gname]] = g[[gname]] + 
        geom_text(aes(x=branch, y=y, label=paste0('nl', branch_id)), size=1, color='gray70', vjust=1.25)
    }
    g[[gname]] = add_long_branch_compression_marks(g[[gname]], args, compression_info)
    return(g)
}

append_branch_stats = function(g, b) {
    g[['data']] = merge(g[['data']], b, by.x='label', by.y='node_name', all.x=TRUE)
    colnames(g[['data']]) = sub('\\.x$', '', colnames(g[['data']]))
    g[['data']] = g[['data']][order(g[['data']][['node']]),]
    rownames(g[['data']]) = NULL
    return(g)
}

add_node_labels = function(g, args, nodelabel_col) {
    if (!nodelabel_col %in% c('no','')) {
        font_size = args[['font_size']] * args[['font_size_factor']]
        g = g + geom_nodelab(
            mapping=aes(x=branch, y=y, label=!!rlang::sym(nodelabel_col)),
            nudge_x=0,
            nudge_y=0, 
            geom='text', 
            colour=args[['nodelabel_color']],
            hjust=0.5, 
            vjust=-0.25,
            size=font_size
        )
    } else {
        cat('Node labels were not added: "', nodelabel_col, '"\n')
    }
    return(g)
}

add_signal_peptide_column = function(g, args) {
    cat(as.character(Sys.time()), 'Adding signal peptide column.\n')
    gname = 'signal_peptide'    
    cols = c('targetp_noTP','targetp_SP','targetp_mTP','targetp_cTP','targetp_luTP')
    df_tip = get_df_tip(g[['tree']])
    if (!all(cols %in% colnames(df_tip))) {
        cat('TargetP info were not found. Signal peptide column will not be shown.\n')
        return(g)
    }
    colnames(df_tip) = sub('targetp_', '', colnames(df_tip))
    newcols = sub('targetp_', '', cols)
    df_tip2 = reshape_long_base(
        df = df_tip[, c('y', newcols), drop = FALSE],
        id_col = 'y',
        value_cols = newcols,
        key_col = 'key',
        value_col = 'value'
    )
    df_tip2 = merge(df_tip2, df_tip[,c('label','y')], by='y', all.x=TRUE)
    g[[gname]] = ggplot(data=df_tip2) +
        geom_blank(aes(y=label)) + 
        geom_bar(mapping=aes(x=value, y=label, fill=key), data=df_tip2, position='fill', stat='identity') +
        coord_cartesian(clip="off") +
        ylab(NULL) +
        xlab(NULL) +
        theme_minimal(base_size=args[['font_size']]) + 
        labs(fill='Signal peptide') +
        guides(
            fill=guide_legend(title=NULL, nrow=6, byrow=FALSE)
        ) +
        theme(
            axis.title.x=element_text(size=args[['font_size']]),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid  = element_blank(),
            legend.position="bottom",
            legend.title=element_text(size=args[['font_size']]),
            legend.text=element_text(size=args[['font_size']]),
            legend.box.just='center',
            legend.key.size=unit(0.4, 'lines'), 
            rect=element_rect(fill="transparent"),
            plot.margin=unit(args[['margins']]/4, "cm")
        )
    return(g)
}

add_integer_column = function(g, args, gname, col, xlab) {
    cat(as.character(Sys.time()), 'Adding integer column.\n')
    df_tip = get_df_tip(g[['tree']])
    if (! col %in% colnames(df_tip)) {
        cat('Column not found: ', col, '. ', gname, ' will not be displayed.\n')
        return(g)
    }
    df_tip[,'x_dummy'] = 0.5
    df_tip[,'hjust'] = 0.5
    df_tip[is.na(df_tip[[col]]),col] = 0
    df_tip[,col] = as.integer(df_tip[[col]])
    g[[gname]] = ggplot(df_tip, aes(x=x_dummy, y=label, label=!!rlang::sym(col), hjust=hjust)) +
        geom_blank(aes(y=label)) + 
        geom_text(size=args[['font_size']] * args[['font_size_factor']], colour=df_tip[['tiplab_color']]) +
        xlim(0,1) +
        xlab(xlab) +
        theme(
            axis.title.y=element_blank(),
            axis.title.x=element_text(size=args[['font_size']]),
            axis.line=element_blank(), 
            axis.ticks=element_blank(), 
            axis.text.y=element_blank(),
            axis.text.x=element_blank(),
            panel.grid.major.y=element_blank(),
            rect=element_blank(),
            plot.margin=unit(args[['margins']]/4, "cm")
        )
    return(g)
}

