get_df_domain <- function(df_rps) {
  domain_values <- as.character(df_rps[['sacc']])
  domain_values <- domain_values[!is.na(domain_values) & nzchar(domain_values)]
  domain_order <- character(0)
  if (length(domain_values) > 0) {
    domain_counts <- table(domain_values)
    domain_order <- names(domain_counts)[order(-as.numeric(domain_counts), names(domain_counts), method = 'radix')]
  }
  domain_rank <- setNames(seq_along(domain_order), domain_order)
  
  # Prepare output
  out_rows <- list()
  out_idx <- 1
  
  # Process each 'qacc' separately
  for (qacc in unique(df_rps$qacc)) {
    
    tmp <- df_rps[df_rps$qacc == qacc, ]
    if (nrow(tmp) == 0) next
    
    qlen <- suppressWarnings(as.integer(tmp[1, "qlen"]))
    if (is.na(qlen) || qlen <= 0) {
      next
    }
    ymin <- tmp[1, "ymin"]
    ymax <- tmp[1, "ymax"]
    if (!is.finite(ymin) || !is.finite(ymax)) {
      next
    }
    tmp$qstart <- suppressWarnings(as.integer(tmp$qstart))
    tmp$qend <- suppressWarnings(as.integer(tmp$qend))
    tmp <- tmp[is.finite(tmp$qstart) & is.finite(tmp$qend), , drop = FALSE]
    if (nrow(tmp) == 0) next
    tmp$qstart <- pmax(1L, pmin(qlen, tmp$qstart))
    tmp$qend <- pmax(1L, pmin(qlen, tmp$qend))
    tmp <- tmp[tmp$qstart <= tmp$qend, , drop = FALSE]
    if (nrow(tmp) == 0) next
    
    # Create "events" for each domain start/end
    #   event = +1 when domain starts
    #   event = -1 when domain ends
    events <- data.frame(
      pos   = c(tmp$qstart, tmp$qend + 1),    # end at qend => next position is qend+1
      sacc  = c(tmp$sacc,   tmp$sacc),
      event = c(rep(+1, nrow(tmp)), 
                rep(-1, nrow(tmp))),
      stringsAsFactors = FALSE
    )
    events <- na.omit(events)
    
    # Also ensure we don't go beyond qlen + 1
    # (in case any domain ends exactly at qlen; qend+1 would be qlen+1)
    # and at least one event at pos=1 if needed
    # -- If you want explicit guaranteed boundaries, you might do something like:
    # events <- rbind(events, data.frame(pos=1, sacc="", event=0))
    # events <- rbind(events, data.frame(pos=qlen+1, sacc="", event=0))
    # but this depends on how strictly you'd like to define bounding events.
    
    # Sort events by position
    events <- events[order(events$pos), ]
    
    # Track current domain set, start of that set's range
    current_set   <- character()  # empty vector
    current_start <- 1
    
    # Traverse all events
    for (i in seq_len(nrow(events))) {
      new_pos <- events$pos[i]
      
      # Build the interval [current_start, new_pos - 1] for the old set
      # only if the old set is not empty AND we haven't gone past qlen
      if (length(current_set) > 0 && current_start <= qlen && new_pos - 1 >= current_start) {
        current_set_sorted <- current_set
        if (length(current_set_sorted) > 1) {
          ord_set <- order(
            ifelse(current_set_sorted %in% names(domain_rank), domain_rank[current_set_sorted], Inf),
            current_set_sorted,
            method = 'radix',
            na.last = TRUE
          )
          current_set_sorted <- current_set_sorted[ord_set]
        }
        step_value <- (ymax - ymin) / length(current_set_sorted)
        borders    <- seq(ymin, ymax, step_value)
        
        # One row per active domain
        for (j in seq_along(current_set_sorted)) {
          out_rows[[out_idx]] <- data.frame(
            xmin  = current_start - 1,
            xmax  = min(new_pos - 1, qlen - 1),
            ymin  = borders[j],
            ymax  = borders[j + 1],
            sacc  = current_set_sorted[j],
            label = qacc,
            stringsAsFactors = FALSE
          )
          out_idx <- out_idx + 1
        }
      }
      
      # Update the domain set based on this event
      if (events$event[i] == +1) {
        # Domain is starting
        current_set <- union(current_set, events$sacc[i])
      } else {
        # Domain is ending
        current_set <- setdiff(current_set, events$sacc[i])
      }
      
      # Move our start pointer to the new position
      current_start <- new_pos
    }
    
    # If anything is still active after the last event (and within qlen)
    if (length(current_set) > 0 && current_start <= qlen) {
      current_set_sorted <- current_set
      if (length(current_set_sorted) > 1) {
        ord_set <- order(
          ifelse(current_set_sorted %in% names(domain_rank), domain_rank[current_set_sorted], Inf),
          current_set_sorted,
          method = 'radix',
          na.last = TRUE
        )
        current_set_sorted <- current_set_sorted[ord_set]
      }
      step_value <- (ymax - ymin) / length(current_set_sorted)
      borders    <- seq(ymin, ymax, step_value)
      for (j in seq_along(current_set_sorted)) {
        out_rows[[out_idx]] <- data.frame(
          xmin  = current_start - 1,
          xmax  = qlen - 1,
          ymin  = borders[j],
          ymax  = borders[j + 1],
          sacc  = current_set_sorted[j],
          label = qacc,
          stringsAsFactors = FALSE
        )
        out_idx <- out_idx + 1
      }
    }
  }
  if (length(out_rows) == 0) {
    out <- data.frame(
      xmin = numeric(0),
      xmax = numeric(0),
      ymin = numeric(0),
      ymax = numeric(0),
      sacc = character(0),
      label = character(0),
      stringsAsFactors = FALSE
    )
    return(out)
  }
  out <- do.call(rbind, out_rows)
  
  detected_levels <- unique(as.character(out$sacc))
  if (length(domain_order) == 0) {
    sacc_levels <- detected_levels[order(detected_levels, method = 'radix')]
  } else {
    missing_levels <- detected_levels[!(detected_levels %in% domain_order)]
    missing_levels <- missing_levels[order(missing_levels, method = 'radix')]
    sacc_levels <- c(domain_order, missing_levels)
  }
  out$sacc <- factor(out$sacc, levels = sacc_levels)
  
  out
}

get_df_intron = function(df_tip) {
    intron_rows = list()
    row_idx = 1
    for (i in 1:nrow(df_tip)) {
        if (is.na(df_tip[i,'intron_positions'])) {
            next
        }
        if (df_tip[i,'intron_positions']=='') {
            next
        }
        intron_val = as.character(df_tip[i,'intron_positions'])
        # gff2genestat now emits ';' delimited positions. Keep ',' for backward compatibility.
        intron_pos = strsplit(intron_val, '[;,]')[[1]]
        intron_pos = trimws(intron_pos)
        intron_pos = intron_pos[intron_pos != '']
        if (length(intron_pos) == 0) {
            next
        }
        intron_site = as.numeric(intron_pos)/3
        df_tmp = data.frame(x=intron_site)
        df_tmp[,'label'] = df_tip[i,'label']
        intron_rows[[row_idx]] = df_tmp
        row_idx = row_idx + 1
    }
    if (length(intron_rows) == 0) {
        df_intron = data.frame(x = numeric(0), label = character(0), stringsAsFactors = FALSE)
        df_intron[,'label'] = factor(df_intron[['label']], levels=df_tip[['label']])
        return(df_intron)
    }
    df_intron = do.call(rbind, intron_rows)
    df_intron[,'label'] = factor(df_intron[['label']], levels=df_tip[['label']])
    return(df_intron)
}

add_protein_domain_column = function(g, args, df_rpsblast=NULL) {
    cat(as.character(Sys.time()), 'Adding protein domain column.\n')
    if (is.null(df_rpsblast)) {
        cat('df_rpsblast is empty. Protein domain column will not be added.\n')
        return(g)
    }
    df_tip = get_df_tip(g[['tree']])
    df_rps = prepare_df_rps_for_plot(df_rpsblast, df_tip)
    if (nrow(df_rps)==0) {
        cat('No protein domain was detected. The protein domain column will not be added.\n')
        return(g)
    }
    df_pb = get_df_protein_backbone(df_rps)
    df_pb[,'label'] = factor(df_pb[['label']], levels=df_tip[['label']])
    df_domain = get_df_domain(df_rps)
    if ('label' %in% colnames(df_domain)) {
        df_domain[,'label'] = factor(df_domain[,'label'], levels=df_tip[,'label'])
    }
    if ('intron_positions' %in% colnames(df_tip)) {
        cat('Intron position info is detected. Adding it to the protein domain column.\n')
        df_intron = get_df_intron(df_tip)
    }

    g[['domain']] = ggplot(data=df_tip, aes(y=label)) +
        geom_segment(mapping=aes(y=label, yend=label), x=0, xend=max(df_rps[['qlen']]), linewidth=0.25, color='gray50', linetype=3) +
        geom_segment(data=df_pb, aes(x=x, xend=xend, y=label, yend=label), linewidth=0.5, color='black') +
        xlim(0, max(df_rps[['qlen']])) +
        xlab('Amino acid position (aa)') +
        theme_minimal(base_size=args[['font_size']]) +
        guides(
            fill=guide_legend(title=NULL, nrow=6, byrow=FALSE)
        ) +
        theme(
            axis.title.y=element_blank(),
            axis.title.x=element_text(size=args[['font_size']]),
            axis.line=element_blank(), 
            axis.ticks=element_blank(), 
            axis.text.y=element_blank(),
            axis.text.x=element_text(size=args[['font_size']], color='black'),
            panel.grid.major.y=element_blank(),
            legend.position="bottom",
            legend.title=element_text(size=args[['font_size']]),
            legend.text=element_text(size=args[['font_size']]),
            legend.box.just='center',
            legend.key.size=unit(0.4, 'lines'), 
            rect=element_rect(fill="transparent"),
            plot.margin=unit(args[['margins']], "cm")
        )
    if (nrow(df_domain)>0) {
        mycolors = get_domain_fill_colors(df_domain)
        g[['domain']] = g[['domain']] + 
            geom_rect(data=df_domain, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=sacc), alpha=1, inherit.aes=FALSE) +
            scale_fill_manual(values=mycolors)
    }
    if ('intron_positions' %in% colnames(df_tip)) {
        if (!all(is.na(df_tip[['intron_positions']]))) {
            g[['domain']] = g[['domain']] + geom_point(data=df_intron, aes(x=x, y=label), shape='|', color='gray50')
        }
    }
    return(g)
}

correct_self_hit_query = function(df, tree2, mrca_matrix, ortholog_prefix) {
    self_queries = df[(df[['query']]==df[['nearests']]),'query']            
    for (sq in self_queries) {
        df_sq_ortho = df[(df[['nearests']]==sq),]
        if (nrow(df_sq_ortho)>1) {
            ph_mrca = df_sq_ortho[(df_sq_ortho[['query']]!=sq),'mrca'][1]
            new_nearests = df_sq_ortho[(df_sq_ortho[['query']]!=sq),'nearests'][1]
            df[(df[['query']]==sq),'mrca'] = ph_mrca
            df[(df[['query']]==sq),'nearests'] = new_nearests
        #} else {
        #    new_nearests = df[grepl(sq, df[['nearests']]),'nearests']
        #    if (length(new_nearests)>1) {                    
        #        all_except_sq = tree2[['tip.label']][tree2[['tip.label']]!=sq]
        #        tmp = rkftools::get_nearest_tips(phy=tree2, sq, all_except_sq, mrca_matrix)
        #        #if (all(startsWith(tmp[['nearests']], ortholog_prefix))) {
        #        if (!all(tmp[['nearests']]==sq)) {
        #            new_nearests = new_nearests[new_nearests!=sq]
        #        }
        #        new_nearests = new_nearests[nchar(new_nearests)==min(nchar(new_nearests))][1]
        #        ph_mrca = df[(df[['nearests']]==new_nearests),'mrca'][1]
        #    } else {
        #        new_nearests = sq
        #        ph_mrca = df_sq_ortho[['mrca']]
        #    }            
        }
        #df[(df[['query']]==sq),'mrca'] = ph_mrca
        #df[(df[['query']]==sq),'nearests'] = new_nearests
    }
    return(df)
}

get_label_num = function(nearests) {
    match_positions = gregexpr('\n', nearests)[[1]]
    label_num = ifelse(match_positions[1]==-1, 1, length(match_positions)+1)
    return(label_num)
}

get_df_clade = function(df, ortholog_prefix) {
    clade_bar_offset = 0.3 # should be < 0.5
    if (nrow(df) == 0) {
        return(data.frame(
            nearests = character(0),
            ymin = numeric(0),
            ymid = numeric(0),
            ymax = numeric(0),
            is_self = logical(0),
            num_line = numeric(0),
            num_label = numeric(0),
            stringsAsFactors = FALSE
        ))
    }
    df = df[(order(df[['y']])),]
    rownames(df) = NULL
    runs = rle(as.character(df[['nearests']]))
    end_idx = cumsum(runs[['lengths']])
    start_idx = end_idx - runs[['lengths']] + 1
    nearests = runs[['values']]
    ymin = df[['y']][start_idx] - clade_bar_offset
    ymax = df[['y']][end_idx] + clade_bar_offset
    ymid = (ymin + ymax) / 2
    is_self = startsWith(as.character(df[['query']][start_idx]), ortholog_prefix)
    df_clade = data.frame(
        nearests = nearests,
        ymin = as.numeric(ymin),
        ymid = as.numeric(ymid),
        ymax = as.numeric(ymax),
        is_self = as.logical(is_self),
        stringsAsFactors = FALSE
    )
    df_clade[['nearests']] = gsub(ortholog_prefix, '', df_clade[['nearests']])
    df_clade[['num_line']] = ceiling(df_clade[['ymax']] - df_clade[['ymin']])
    df_clade[['num_label']] = sapply(df_clade[['nearests']], get_label_num)
    return(df_clade)
}

abbreviate_clade_labels = function(df_clade, num_line_offset=1) {
    for (i in 1:nrow(df_clade)) {
        num_line = max(c(1, df_clade[i,'num_line'] - num_line_offset))
        num_label = df_clade[i,'num_label']
        if (num_line < num_label) {
            nearests = df_clade[i,'nearests']
            match_positions = gregexpr('\n', nearests)[[1]]
            abbreviated_nearests = substr(nearests, 1, match_positions[num_line]-1)
            num_rest = num_label - num_line
            abbreviated_nearests = paste0(abbreviated_nearests, ' + ', num_rest)
            df_clade[i,'nearests'] = abbreviated_nearests
        }
    }
    return(df_clade)    
}

split_nested_clade = function(df_clade, margin=0.3) {
    if (nrow(df_clade) <= 1) {
        return(df_clade)
    }
    min_span = (1 - (margin * 2))
    repeat {
        changed = FALSE
        n = nrow(df_clade)
        for (i in seq_len(n)) {
            for (j in seq_len(n)) {
                if (i == j) {
                    next
                }
                is_nested = (df_clade[i,'ymin'] < df_clade[j,'ymin']) && (df_clade[i,'ymax'] > df_clade[j,'ymax'])
                if (!is_nested) {
                    next
                }
                cat('Splitting nested clade orthologs. Row numbers', i, 'and', j, '\n')
                base_row = df_clade[i, , drop = FALSE]
                upper = base_row
                lower = base_row
                upper[,'ymax'] = df_clade[j,'ymin'] - min_span
                upper[,'ymid'] = upper[,'ymin'] + ((upper[,'ymax'] - upper[,'ymin'])/2)
                lower[,'ymin'] = df_clade[j,'ymax'] + min_span
                lower[,'ymid'] = lower[,'ymin'] + ((lower[,'ymax'] - lower[,'ymin'])/2)
                df_clade = df_clade[-i, , drop = FALSE]
                candidates = rbind(df_clade, upper, lower)
                span = abs(candidates[['ymax']] - candidates[['ymin']])
                df_clade = candidates[span > min_span, , drop = FALSE]
                rownames(df_clade) = NULL
                changed = TRUE
                break
            }
            if (changed) {
                break
            }
        }
        if (!changed) {
            break
        }
    }
    return(df_clade)
}

