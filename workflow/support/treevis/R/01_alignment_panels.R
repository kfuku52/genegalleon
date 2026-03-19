add_alignment_column <- function(g, args, seqs = NULL, df_rpsblast = NULL, seqs_untrim = NULL) {
  cat(as.character(Sys.time()), 'Adding alignment column.\n')
  
  if (is.null(seqs)) {
    cat('Input alignment is empty. Alignment column will not be added.\n')
    return(g)
  }
  
  df_tip <- get_df_tip(g[['tree']])
  seq_names <- names(seqs)
  if (is.null(seq_names)) {
    seq_names <- character(0)
  }
  missing_labels <- setdiff(as.character(df_tip[['label']]), seq_names)
  if (length(missing_labels) > 0) {
    preview <- paste(utils::head(missing_labels, 10), collapse = ',')
    cat('Alignment is missing', length(missing_labels), 'tip sequence(s). Displaying backbone only. Examples:', preview, '\n')
  }
  df_rps = prepare_df_rps_for_plot(df_rpsblast, df_tip)
  domain_fill_colors = c()
  if (nrow(df_rps) > 0) {
    df_domain = get_df_domain(df_rps)
    domain_fill_colors = get_domain_fill_colors(df_domain)
  }
  non_domain_fill = '__non_domain__'
  fill_colors = c(setNames('gray30', non_domain_fill), domain_fill_colors)
  domain_rank = setNames(seq_along(names(domain_fill_colors)), names(domain_fill_colors))
  
  key_sep <- ':::DOMAINSEP:::'

  # We'll accumulate data frames for each tip in a list (faster than rbind in a loop)
  out_list <- vector("list", length(df_tip[['label']]))
  idx <- 1
  xmax_val <- 0
  
  for (seqname in df_tip[['label']]) {
    seqname <- as.character(seqname)
    # Get the sequence associated with this tip
    seq_vec <- as.character(seqs[[seqname]])
    if (length(seq_vec) == 0) next
    xmax_val <- max(xmax_val, length(seq_vec) - 1)
    
    # Convert to a logical vector: TRUE if not '04', FALSE if '04'
    is_atgc <- (seq_vec != '04')
    mapped_untrim_nt = rep(NA_integer_, length(seq_vec))
    has_untrim_map = FALSE
    if (!is.null(seqs_untrim) && !is.null(seqs_untrim[[seqname]])) {
      seq_untrim_vec = as.character(seqs_untrim[[seqname]])
      if (length(seq_untrim_vec) > 0) {
        mapped_untrim_nt = map_trimmed_to_untrimmed_nongap_index(seq_trim = seq_vec, seq_untrim = seq_untrim_vec, gap_code = '04')
        has_untrim_map = any(!is.na(mapped_untrim_nt[is_atgc]))
      }
    }

    active_domains <- vector("list", length(seq_vec))
    if (nrow(df_rps) > 0) {
      df_seq_rps <- df_rps[(df_rps[['qacc']] == seqname), c('sacc', 'qstart', 'qend', 'qlen'), drop = FALSE]
      if (nrow(df_seq_rps) > 0) {
        nt_pos_trim = cumsum(is_atgc)
        seq_aa_len <- suppressWarnings(as.integer(df_seq_rps[1, 'qlen']))
        if (!is.na(seq_aa_len) && seq_aa_len > 0) {
          ord <- order(
            ifelse(df_seq_rps[['sacc']] %in% names(domain_rank), domain_rank[df_seq_rps[['sacc']]], Inf),
            suppressWarnings(as.numeric(df_seq_rps[['qstart']])),
            method = 'radix',
            na.last = TRUE
          )
          for (k in ord) {
            qstart_aa <- suppressWarnings(as.integer(df_seq_rps[k, 'qstart']))
            qend_aa <- suppressWarnings(as.integer(df_seq_rps[k, 'qend']))
            if (!is.finite(qstart_aa) || !is.finite(qend_aa)) {
              next
            }
            nt_start <- (qstart_aa - 1L) * 3L + 1L
            nt_end <- qend_aa * 3L
            nt_start <- max(1L, nt_start)
            nt_end <- min(seq_aa_len * 3L, nt_end)
            if (nt_start > nt_end) {
              next
            }
            if (has_untrim_map) {
              is_domain_pos <- is_atgc & !is.na(mapped_untrim_nt) & (mapped_untrim_nt >= nt_start) & (mapped_untrim_nt <= nt_end)
            } else {
              is_domain_pos <- is_atgc & (nt_pos_trim >= nt_start) & (nt_pos_trim <= nt_end)
            }
            hit_idx <- which(is_domain_pos)
            if (length(hit_idx) > 0) {
              dlabel <- as.character(df_seq_rps[k, 'sacc'])
              for (h in hit_idx) {
                if (is.null(active_domains[[h]]) || length(active_domains[[h]]) == 0) {
                  active_domains[[h]] <- dlabel
                } else if (!(dlabel %in% active_domains[[h]])) {
                  active_domains[[h]] <- c(active_domains[[h]], dlabel)
                }
              }
            }
          }
        }
      }
    }
    key_vec <- rep(NA_character_, length(seq_vec))
    valid_idx <- which(is_atgc)
    for (i in valid_idx) {
      labels_i <- active_domains[[i]]
      if (is.null(labels_i) || length(labels_i) == 0) {
        key_vec[i] <- non_domain_fill
      } else {
        labels_i <- unique(as.character(labels_i))
        ord_i <- order(
          ifelse(labels_i %in% names(domain_rank), domain_rank[labels_i], Inf),
          labels_i,
          method = 'radix',
          na.last = TRUE
        )
        key_vec[i] <- paste(labels_i[ord_i], collapse = key_sep)
      }
    }

    # Run-length encode
    runs <- rle(key_vec)

    # Identify the start and end indices of each run
    run_ends   <- cumsum(runs$lengths)
    run_starts <- run_ends - runs$lengths + 1
    
    # Keep only runs in non-gap positions
    keep_idx <- which(!is.na(runs$value))
    if (length(keep_idx) > 0) {
      local_rows <- list()
      local_n <- 1
      for (rk in keep_idx) {
        xmin_val <- run_starts[rk] - 1
        xmax_val_run <- run_ends[rk] - 1
        key_val <- runs$value[rk]
        if (identical(key_val, non_domain_fill)) {
          local_rows[[local_n]] <- data.frame(
            xmin = xmin_val,
            xmax = xmax_val_run,
            label = seqname,
            fill = non_domain_fill,
            split_index = 1L,
            split_total = 1L,
            stringsAsFactors = FALSE
          )
          local_n <- local_n + 1
        } else {
          labels_run <- strsplit(as.character(key_val), key_sep, fixed = TRUE)[[1]]
          n_split <- length(labels_run)
          if (n_split == 0) {
            next
          }
          for (j in seq_len(n_split)) {
            local_rows[[local_n]] <- data.frame(
              xmin = xmin_val,
              xmax = xmax_val_run,
              label = seqname,
              fill = labels_run[j],
              split_index = as.integer(j),
              split_total = as.integer(n_split),
              stringsAsFactors = FALSE
            )
            local_n <- local_n + 1
          }
        }
      }
      if (length(local_rows) > 0) {
        df_local <- do.call(rbind, local_rows)
      } else {
        df_local <- data.frame(xmin = numeric(0), xmax = numeric(0), label = character(0), fill = character(0), split_index = integer(0), split_total = integer(0))
      }
    } else {
      # If there are no TRUE runs, just store an empty data frame
      df_local <- data.frame(xmin = numeric(0), xmax = numeric(0), label = character(0), fill = character(0), split_index = integer(0), split_total = integer(0))
    }
    
    out_list[[idx]] <- df_local
    idx <- idx + 1
  }
  
  # Combine all individual tip data frames
  df_aln <- do.call(rbind, out_list)
  
  # Factor tip labels in the order they appear in df_tip
  df_aln[['label']] <- factor(df_aln[['label']], levels = df_tip[['label']])
  df_aln[['fill']] <- factor(df_aln[['fill']], levels = names(fill_colors))
  
  # Add y-coordinates for geom_rect
  y_start <- as.numeric(df_aln[['label']]) - 0.45
  y_step <- 0.9 / pmax(1, as.numeric(df_aln[['split_total']]))
  df_aln[['ymin']] <- y_start + (as.numeric(df_aln[['split_index']]) - 1) * y_step
  df_aln[['ymax']] <- y_start + as.numeric(df_aln[['split_index']]) * y_step
  
  # Create a "backbone" segment spanning the alignment range for each tip
  df_backbone <- data.frame(
    label = df_tip[['label']],
    x     = 0,
    xend  = xmax_val
  )
  
  # Build the alignment ggplot object
  g[['alignment']] <- ggplot(data = df_tip, aes(y = label)) +
    geom_segment(
      data    = df_backbone,
      aes(x = x, xend = xend, y = label, yend = label),
      linewidth = 0.25, color = 'gray50', linetype = 3
    ) +
    geom_rect(
      data = df_aln,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
      alpha = 1
    ) +
    scale_fill_manual(values = fill_colors, guide = 'none') +
    xlim(0, xmax_val) +
    xlab('Alignment position (nt)') +
    theme_minimal(base_size = args[['font_size']]) +
    theme(
      axis.title.y        = element_blank(),
      axis.title.x        = element_text(size = args[['font_size']]),
      axis.line           = element_blank(), 
      axis.ticks          = element_blank(), 
      axis.text.y         = element_blank(),
      axis.text.x         = element_text(size = args[['font_size']], color = 'black'),
      panel.grid.major.y  = element_blank(),
      legend.position     = 'none',
      rect                = element_rect(fill = "transparent"),
      plot.margin         = unit(args[['margins']], "cm")
    )
  
  return(g)
}

get_df_trait = function(b, transform, scale, trait_prefix, negative2zero=TRUE) {
    is_leaf = (b[['so_event']]=='L')
    df_trait = b[is_leaf,startsWith(colnames(b), trait_prefix), drop=FALSE]
    rownames(df_trait) = b[is_leaf,'node_name']
    colnames(df_trait) = sub(paste0('^', trait_prefix), '', colnames(df_trait))
    if (ncol(df_trait)==0) {
        cat('Trait columns not found in the stat_branch table.\n')
        df_trait = data.frame()
        return(df_trait)
    }
    # Ensure all values are numeric; non-numeric values become NA.
    df_trait[] = lapply(df_trait, function(x) suppressWarnings(as.numeric(x)))
    if (negative2zero) {
        cat('Negative trait values will be replaced with 0.\n')
        df_trait[df_trait<0] = 0
    }
    if (transform=='log2') {
        cat('Trait data are transformed by log2(x).', '\n')
        df_trait = log2(df_trait)
    } else if (transform=='log2p1') {
        cat('Trait data are transformed by log2(x+1).', '\n')
        df_trait = log2(df_trait+1)
    } else if (transform=='log10') {
        cat('Trait data are transformed by log10(x).', '\n')
        df_trait = log10(df_trait)
    } else if (transform=='log10p1') {
        cat('Trait data are transformed by log10(x+1).', '\n')
        df_trait = log10(df_trait+1)
    } else {
        cat('Trait data are not transformed.', '\n')
    }
    if (scale=='abs') {
        cat('Absolute trait values will be displayed.\n')
    } else if (scale=='rel') {
        cat('Relative trait values (x/(max(xs_leaf))) will be displayed.\n')
        df_trait_fillzero = df_trait
        df_trait_fillzero[is.na(df_trait_fillzero)] = 0
        row_max = apply(df_trait_fillzero, 1, function(x){max(na.omit(x))})
        row_max[row_max == 0] = 1
        df_trait = df_trait / row_max
        is_all_na = apply(is.na(df_trait), 1, sum)==apply(df_trait, 1, length)
        df_trait[is_all_na,] = 0
    }
    return(df_trait)
}

add_heatmap_column = function(g, args, df_trait, fill_label='Expression') {
    cat(as.character(Sys.time()), 'Adding heatmap column.\n')
    if ((is.null(ncol(df_trait)))|(ncol(df_trait)==0)) {
        cat('df_trait is emply. Heatmap panel will not be added.\n')
        return(g)
    }
    font_size = args[['font_size']]
    if (any(grepl('^pointplot', unlist(args[grep("^panel", names(args))])))) {
        trait_colors = args[['trait_colors']]
    } else {
        trait_colors = 'black'
    }
    df_tip = get_df_tip(g[['tree']])
    df_tip = merge_tip_trait(df_tip, df_trait)
    df_tip_tidy = tidy_df_tip(df_tip, df_trait)
    df_tip_tidy[,'group'] = as.factor(df_tip_tidy[['group']])
    finite_values = df_tip_tidy[['value']][is.finite(df_tip_tidy[['value']])]
    if (length(finite_values) == 0) {
        cat('All heatmap values are NA/NaN/Inf. Heatmap panel will not be added.\n')
        return(g)
    }
    max_val = max(finite_values)
    if (max_val <= 0) {
        max_val = 1
    }
    grid_color = ifelse(any(is.na(df_tip_tidy[['value']])), rgb(0,0,0,0), rgb(0,0,0,1))
    g[['heatmap']] = ggplot(data=df_tip_tidy)
    if (any(apply(df_tip[, colnames(df_trait), drop = FALSE], 1, function(x){all(is.na(x))}))) {
        df_tip_tidy_seg = unique(df_tip_tidy[,c('label','y')])
        df_tip_tidy_seg[,'xmin'] = df_tip_tidy[which.min(df_tip_tidy[['group']]),'group']
        df_tip_tidy_seg[,'xmax'] = df_tip_tidy[which.max(df_tip_tidy[['group']]),'group']
        g[['heatmap']] = g[['heatmap']] + geom_segment(
            data=df_tip_tidy_seg, 
            mapping=aes(y=label, yend=label, x=xmin, xend=xmax), 
            linewidth=0.25, 
            color='gray50', 
            linetype=3
        )
    }
    g[['heatmap']] = g[['heatmap']] + 
        geom_tile(mapping=aes(x=group, y=label, fill=value), color=grid_color) +
        scale_fill_gradient2(
            low="#00204c", 
            mid="#7b7b78", 
            high="#ffe945", 
            midpoint=max_val/2, 
            limits=c(0, max_val), 
            na.value=rgb(1,1,1,0)
        ) +
        ggplot2::labs(fill=fill_label) +
        theme_classic(base_size=font_size) +
        theme(
            axis.title=element_blank(), 
            axis.line=element_blank(), 
            axis.ticks=element_blank(), 
            axis.text.y=element_blank(),
            axis.text.x=element_text(angle=45, hjust=1, colour=trait_colors, size=font_size),
            legend.position="bottom",
            legend.title=element_text(size=font_size),
            legend.text=element_text(size=font_size),
            legend.box.just='center',
            rect=element_rect(fill="transparent"),
            plot.margin=unit(args[['margins']], "cm")
        ) +
        guides(
            fill=guide_colourbar(
                title.position="bottom", 
                title.hjust=0.25, 
                label.position="bottom", 
                barheight=0.3,
                barwidth=ncol(df_trait)*0.7
            )
        )
    return(g)
}

add_tiplabel_column = function(g, args) {
    cat(as.character(Sys.time()), 'Adding tip label column.\n')
    df_tip = get_df_tip(g[['tree']])
    df_tip[,'x_dummy'] = 0
    df_tip[,'hjust'] = 0
    g[['tiplabel']] = ggplot(df_tip, aes(x=x_dummy, y=label, label=label, hjust=hjust)) +
        geom_text(size=args[['font_size']] * args[['font_size_factor']], colour=df_tip[['tiplab_color']]) +
        xlim(0,1) +
        theme_void() +
        theme(
            plot.margin=unit(args[['margins']], "cm")
        )
    return(g)
}

truncate_panel_text = function(values, max_nchar = 18) {
    max_nchar = suppressWarnings(as.integer(max_nchar))
    if (is.na(max_nchar)) {
        max_nchar = 18
    }
    values = as.character(values)
    values[is.na(values)] = ''
    values = trimws(values)
    values = gsub('[[:space:]]+', ' ', values)
    if (max_nchar < 1) {
        return(values)
    }
    too_long = (nchar(values, type = 'width') > max_nchar)
    if (max_nchar <= 3) {
        values[too_long] = substr(values[too_long], 1, max_nchar)
    } else {
        values[too_long] = paste0(substr(values[too_long], 1, max_nchar - 3), '...')
    }
    return(values)
}

add_text_column = function(g, args, gname, col, xlab, max_nchar = 18) {
    cat(as.character(Sys.time()), 'Adding text column.\n')
    df_tip = get_df_tip(g[['tree']])
    if (! col %in% colnames(df_tip)) {
        cat('Column not found: ', col, '. ', gname, ' will not be displayed.\n')
        return(g)
    }
    values = truncate_panel_text(df_tip[[col]], max_nchar = max_nchar)
    if (all(values == '')) {
        cat('All values were empty in column: ', col, '. ', gname, ' will not be displayed.\n')
        return(g)
    }
    df_tip[['x_dummy']] = 0
    df_tip[['hjust']] = 0
    df_tip[['plot_value']] = values
    text_margins = args[['margins']] / 4
    text_margins[2] = text_margins[2] + 0.12
    g[[gname]] = ggplot(df_tip, aes(x=x_dummy, y=label, label=plot_value, hjust=hjust)) +
        geom_text(size=args[['font_size']] * args[['font_size_factor']], colour=df_tip[['tiplab_color']], lineheight=0.92) +
        coord_cartesian(xlim = c(0, 1), clip = 'off') +
        xlab(xlab) +
        theme_void() +
        theme(
            axis.title.x=element_text(size=args[['font_size']]),
            plot.margin=unit(text_margins, "cm")
        )
    return(g)
}

normalize_category_text = function(values, missing_label = '-') {
    values = as.character(values)
    values[is.na(values)] = ''
    values = trimws(values)
    values = gsub('_', ' ', values, fixed = TRUE)
    values = gsub('[[:space:]]+', ' ', values)
    values[values == ''] = missing_label
    return(values)
}

get_categorical_levels = function(values, missing_label = '-') {
    values = unique(as.character(values))
    values = values[values != missing_label]
    preferred = c(
        'species',
        'genus',
        'family',
        'order',
        'class',
        'phylum',
        'kingdom',
        'superkingdom',
        'domain',
        'genus mismatch'
    )
    preferred_present = preferred[preferred %in% values]
    other_values = sort(setdiff(values, preferred_present))
    c(preferred_present, other_values, missing_label)
}

add_categorical_column = function(g, args, gname, col, xlab, missing_label = '-') {
    cat(as.character(Sys.time()), 'Adding categorical column.\n')
    df_tip = get_df_tip(g[['tree']])
    if (! col %in% colnames(df_tip)) {
        cat('Column not found: ', col, '. ', gname, ' will not be displayed.\n')
        return(g)
    }
    values = normalize_category_text(df_tip[[col]], missing_label = missing_label)
    levels = get_categorical_levels(values, missing_label = missing_label)
    palette = grDevices::hcl.colors(max(1, length(levels) - 1), palette = 'Dark 3')
    names(palette) = levels[seq_len(max(1, length(levels) - 1))]
    if (missing_label %in% levels) {
        palette[missing_label] = '#e6e6e6'
    }
    df_tip[['plot_value']] = factor(values, levels = levels)
    df_tip[['panel_x']] = 1
    g[[gname]] = ggplot(df_tip, aes(x = panel_x, y = label, fill = plot_value)) +
        geom_tile(width = 0.9, height = 0.9, color = 'white', linewidth = 0.2) +
        scale_fill_manual(values = palette, drop = FALSE) +
        coord_cartesian(xlim = c(0.5, 1.5), clip = 'off') +
        xlab(xlab) +
        theme_classic(base_size = args[['font_size']]) +
        theme(
            axis.title.y = element_blank(),
            axis.title.x = element_text(size = args[['font_size']]),
            axis.text.y = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            legend.position = 'bottom',
            legend.title = element_blank(),
            legend.text = element_text(size = args[['font_size']]),
            legend.key.size = unit(0.35, 'lines'),
            rect = element_rect(fill = 'transparent'),
            plot.margin = unit(args[['margins']] / 4, 'cm')
        ) +
        guides(fill = guide_legend(nrow = min(4, length(levels)), byrow = TRUE))
    return(g)
}

extract_xml_attr = function(tag, candidates) {
    for (attr in candidates) {
        pattern = paste0(attr, '="([^"]*)"')
        hit = regexec(pattern, tag, perl = TRUE)
        regmatch = regmatches(tag, hit)[[1]]
        if (length(regmatch) >= 2) {
            return(regmatch[2])
        }
    }
    return('')
}

read_meme_motifs = function(path_meme, max_motif = 8) {
    if (is.na(path_meme) || !nzchar(path_meme) || !file.exists(path_meme)) {
        return(data.frame())
    }
    xml_lines = tryCatch(readLines(path_meme, warn = FALSE), error = function(e) character(0))
    if (length(xml_lines) == 0) {
        return(data.frame())
    }
    xml_text = paste(xml_lines, collapse = ' ')
    motif_match = gregexpr('<motif\\b[^>]*>', xml_text, perl = TRUE)
    motif_tags = regmatches(xml_text, motif_match)[[1]]
    if (length(motif_tags) == 0) {
        return(data.frame())
    }
    motif_tags = motif_tags[seq_len(min(length(motif_tags), max_motif))]
    motif_rows = lapply(seq_along(motif_tags), function(i) {
        tag = motif_tags[[i]]
        motif_id = extract_xml_attr(tag, c('alt', 'name', 'id'))
        if (!nzchar(motif_id)) {
            motif_id = paste0('motif_', i)
        }
        width = extract_xml_attr(tag, c('width', 'w'))
        sites = extract_xml_attr(tag, c('sites', 'nsites'))
        evalue = extract_xml_attr(tag, c('e_value', 'evalue'))
        data.frame(
            motif_id = motif_id,
            width = width,
            sites = sites,
            evalue = evalue,
            stringsAsFactors = FALSE
        )
    })
    out = do.call(rbind, motif_rows)
    rownames(out) = NULL
    return(out)
}

add_meme_column = function(g, args, path_meme, max_motif = 8) {
    cat(as.character(Sys.time()), 'Adding meme column.\n')
    df_meme = read_meme_motifs(path_meme, max_motif = max_motif)
    if (nrow(df_meme) == 0) {
        cat('No MEME motifs were detected. MEME panel will not be added.\n')
        return(g)
    }
    motif_text = ifelse(
        nzchar(df_meme[['width']]) | nzchar(df_meme[['sites']]) | nzchar(df_meme[['evalue']]),
        paste0(
            truncate_panel_text(df_meme[['motif_id']], 14),
            ' (w=', ifelse(nzchar(df_meme[['width']]), df_meme[['width']], '?'),
            ', n=', ifelse(nzchar(df_meme[['sites']]), df_meme[['sites']], '?'),
            ', E=', ifelse(nzchar(df_meme[['evalue']]), df_meme[['evalue']], '?'),
            ')'
        ),
        truncate_panel_text(df_meme[['motif_id']], 24)
    )
    df_plot = data.frame(
        x = 0,
        y = rev(seq_len(nrow(df_meme))),
        label = truncate_panel_text(motif_text, 32),
        stringsAsFactors = FALSE
    )
    g[['meme']] = ggplot(df_plot, aes(x = x, y = y, label = label)) +
        geom_text(hjust = 0, size = args[['font_size']] * args[['font_size_factor']], colour = 'black') +
        xlim(0, 1) +
        ylim(0.5, nrow(df_plot) + 0.5) +
        xlab('MEME motif') +
        theme_void() +
        theme(
            axis.title.x = element_text(size = args[['font_size']]),
            plot.margin = unit(args[['margins']]/4, "cm")
        )
    return(g)
}

