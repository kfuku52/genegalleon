attach_neighbor_stats = function(df, neighbor, columns=c()){
    columns = c('branch_id', columns)
    df_tmp = df[,columns]
    colnames(df_tmp) = paste0(neighbor, '_', colnames(df_tmp))
    df = merge(df, df_tmp, by.x=neighbor, by.y=paste0(neighbor,'_branch_id'), all.x=TRUE)
    return(df)
}

get_df_tip = function(g) {
    df_tip = subset(g[['data']], isTip)
    tip_order = with(df_tip, label[order(y, decreasing=FALSE)])
    tip_order_index = match(tip_order, df_tip[['label']])
    df_tip = df_tip[tip_order_index,]
    df_tip[['label']] = factor(df_tip[['label']], levels=df_tip[['label']])
    rownames(df_tip) = NULL
    return(df_tip)
}

tidy_df_tip = function(df_tip, df_trait) {
    trait_cols = colnames(df_trait)
    if (length(trait_cols) == 0) {
        return(data.frame())
    }
    trait_cols = sort(trait_cols)
    # Keep only trait columns to avoid reshaping unrelated columns in df_tip.
    keep_cols = c('label', intersect(trait_cols, colnames(df_tip)))
    df_tmp = tidyr::pivot_longer(
        data = df_tip[, keep_cols, drop = FALSE],
        cols = -label,
        names_to = "group",
        values_to = "value"
    )
    df_tmp[['group']] = factor(df_tmp[['group']], levels = trait_cols)
    # Preserve tree tip order (small y at the top).
    tip_order = with(df_tip, as.character(label[order(y, decreasing = FALSE)]))
    df_tmp[['label']] = factor(as.character(df_tmp[['label']]), levels = tip_order)
    df_tmp = df_tmp[order(df_tmp[['group']], df_tmp[['label']]), , drop = FALSE]
    df_tmp[['value']] = suppressWarnings(as.numeric(df_tmp[['value']]))
    df_tmp = merge(df_tmp, df_tip[, c('label', 'y')], by = 'label', sort = FALSE, all.x = TRUE)
    return(df_tmp)
}

merge_tip_trait = function(df_tip, df_trait) {
    if (is.null(df_trait) || ncol(df_trait) == 0) {
        return(df_tip)
    }
    trait_df = data.frame(label = rownames(df_trait), df_trait, stringsAsFactors = FALSE, check.names = FALSE)
    idx = match(as.character(df_tip[['label']]), trait_df[['label']])
    add_cols = trait_df[idx, setdiff(colnames(trait_df), 'label'), drop = FALSE]
    out = cbind(df_tip, add_cols)
    return(out)
}

get_trait_colors = function(num_trait, method) {
    if (is.null(num_trait)) {
        return(NULL)
    }
    if (num_trait==0) {
        trait_colors = c()
    } else {
        if (method=='FP2020') {
            if (num_trait<=6) {
                trait_colors = c('#1B9E77','#D95F02','#7570B3','#E7298A','#66A61E','#E6AB02')[1:num_trait]
            } else {
                stop('FP2020 works only with <=6 traits.')
            }
        } else if (method=='continuous') {
            trait_colors = viridis::viridis_pal(alpha=1, begin=0, end=0.8, direction=1, option='plasma')(num_trait)
        } else if (method=='discrete') {
            trait_colors = scales::hue_pal()(num_trait)            
        }
    }
    return(trait_colors)
}

geom_divtime = function(g, y=0, linewidth=0.5, fontsize=8, font_size_factor=0.352777778, label='') {
    options(scipen=18)
    xmax = max(g[['data']][['x']])
    if (xmax>1) {
        ndigit = nchar(round(xmax, digits=0)) - 1
    } else {
        ndigit = -nchar(sub('[1-9].*', '', xmax)) - 1
    }
    xunit = 10^ndigit
    if (xmax%/%(10^ndigit)==1) {
        ndigit = ndigit - 1
        xunit = xunit / 2
    }
    yunit = 0.25
    out = list()
    lineend = 'square'
    linejoin = 'round'
    out[[1]] = ggplot2::geom_segment(x=0, xend=xmax, y=y, yend=y, linewidth=linewidth, lineend=lineend, linejoin=linejoin)
    count = 2
    for (i in 1:(xmax%/%xunit+1)) {
        x = max(0, xmax-(xunit*i))
        xend = xmax-(xunit*(i-1))
        out[[count+0]] = ggplot2::geom_segment(x=xend, xend=xend, y=y-yunit, yend=y, linewidth=linewidth, lineend=lineend, linejoin=linejoin)
        out[[count+1]] = ggplot2::annotation_custom(
            grob=grid::textGrob(
                label=xunit*(i-1), 
                hjust=0.5,
                vjust=1.25,
                gp=grid::gpar(cex=1/font_size_factor, fontsize=fontsize)
            ),
            ymin=y-yunit,
            ymax=y-yunit,
            xmin=xend,
            xmax=xend
        )
        count = count+2
    }
    out[[count+1]] = ggplot2::annotation_custom(
        grob=grid::textGrob(
            label=label, 
            hjust=0.5,
            vjust=3,
            gp=grid::gpar(cex=1/font_size_factor, fontsize=fontsize)
        ),
        ymin=y-yunit,
        ymax=y-yunit,
        xmin=xmax/2,
        xmax=xmax/2
    )    
    return(out)
}

add_scale_bar = function(g, args, tree) {
    if (ape::is.ultrametric(tree)) {
        g = g + geom_divtime(
            g, 
            y=0, 
            linewidth=args[['scale_thickness']], 
            fontsize=args[['font_size']]*args[['font_size_factor']], 
            label=args[['scale_label']])
    } else {
        g = g + geom_treescale(
            x=max(g[['data']][['x']])*0.8, 
            y=1.5, 
            fontsize=args[['font_size']]*args[['font_size_factor']], 
            linesize=args[['scale_thickness']], 
            offset=0.25)
    }
    return(g)
}

add_node_points = function(g, args) {
    g = g + geom_nodepoint(
        mapping=aes(subset=!isTip, color=node_category), 
        size = args[['nodepoint_size']],
        position='identity', 
        show.legend=TRUE
    ) + 
    scale_colour_manual(values=args[['node_colors']]) +
    theme(
        legend.title = element_text(size=args[['font_size']]),
        legend.text = element_text(size=args[['font_size']]),
        legend.position = c(0, 1), 
        legend.justification = c(0, 1), 
        legend.background = element_rect(colour=NA, fill=NA),
        legend.key.size=unit(0.4, 'lines')
    ) +
    guides(color=guide_legend(title="Branching event"))
    return(g)
}

get_rel_widths = function(g, args_rel_widths) {
    rel_widths = rep(1, length(g))
    names(rel_widths) = names(g)
    for (gname in names(rel_widths)) {
        if (grepl('^tree', gname)) {
            rel_widths[gname] = 1.5
        } else if (grepl('^heatmap$', gname)) {
            rel_widths[gname] = 0.5
        } else if (grepl('^pointplot$', gname)) {
            rel_widths[gname] = 0.5
        } else if (grepl('^signal_peptide$', gname)) {
            rel_widths[gname] = 0.1
        } else if (grepl('^tm$', gname)) {
            rel_widths[gname] = 0.12
        } else if (grepl('^intron$', gname)) {
            rel_widths[gname] = 0.12
        } else if (grepl('^cluster_membership$', gname)) {
            rel_widths[gname] = 0.5
        } else if (grepl('^domain$', gname)) {
            rel_widths[gname] = 1.0
        } else if (grepl('^alignment$', gname)) {
            rel_widths[gname] = 1.0
        } else if (grepl('^fimo$', gname)) {
            rel_widths[gname] = 1.0
        } else if (grepl('^amino_acid_site$', gname)) {
            built_plot <- ggplot_build(g[["amino_acid_site"]])
            xmax = built_plot$layout$panel_scales_x[[1]]$range$range[2]
            num_amino_acid_site = floor(xmax) - 1
            rel_widths[gname] = 0.09 * num_amino_acid_site
        } else if (grepl('^ortholog,', gname)) {
            rel_widths[gname] = 0.7
        } else {
            rel_widths[gname] = 1
        }
    }
    for (gname in names(rel_widths)) { # tiplabel should be processed in the very last
        if (grepl('^tiplabel$', gname)) {
            rel_widths[gname] = 0 # initialize
            max_nchar = max(nchar(as.character(g[['tiplabel']][['data']][['label']])))
            rel_widths[gname] = min(0.5, max_nchar * 0.006)
            is_other = !grepl('^tiplabel$', names(rel_widths))
            rel_widths[is_other] = rel_widths[is_other] / (sum(rel_widths[is_other])/(1-rel_widths[gname]))
        }
    }
    rel_widths = rel_widths / sum(rel_widths)
    if (!is.null(args_rel_widths) && !is.na(args_rel_widths) && nzchar(args_rel_widths)) {
        splitted = strsplit(args_rel_widths, ',')[[1]]
        keys = splitted[seq(1, length(splitted), 2)]
        cat('Width of folloowing columns will be updated according to --rel_widths:', keys, '\n')
        values = as.numeric(splitted[seq(2, length(splitted), 2)])
        for (i in 1:length(keys)) {
            rel_widths[startsWith(names(rel_widths), keys[i])] = values[i]
        }        
    }
    return(rel_widths)
}

add_pointplot_column = function(g, args, df_trait, replicate_sep) {
    cat(as.character(Sys.time()), 'Adding point plot column.\n')
    if ((is.null(ncol(df_trait)))|(ncol(df_trait)==0)) {
        cat('df_trait is emply. Pointplot panel will not be added.\n')
        return(g)
    }
    font_size = args[['font_size']]
    xlabel = args[['trait_axis_label']]
    df_tip = get_df_tip(g[['tree']])
    df_tip = merge_tip_trait(df_tip, df_trait)
    df_tip_tidy = tidy_df_tip(df_tip, df_trait)
    # Remove replicate suffix (e.g. "_1") in a vectorized way.
    if (!is.null(replicate_sep) && nzchar(replicate_sep)) {
        group_chr = as.character(df_tip_tidy[['group']])
        parts = strsplit(group_chr, replicate_sep, fixed = TRUE)
        group_chr = vapply(parts, function(x) {
            if (length(x) <= 1) {
                return(x[1])
            }
            paste(x[1:(length(x)-1)], collapse = replicate_sep)
        }, character(1))
        df_tip_tidy[,'group'] = group_chr
    } else {
        df_tip_tidy[,'group'] = as.character(df_tip_tidy[['group']])
    }
    df_tip_tidy[,'group'] = factor(df_tip_tidy[['group']], levels=unique(sort(df_tip_tidy[['group']])))
    finite_values = df_tip_tidy[['value']][is.finite(df_tip_tidy[['value']])]
    if (length(finite_values) == 0) {
        cat('All pointplot values are NA/NaN/Inf. Pointplot panel will not be added.\n')
        return(g)
    }
    xtickmin = floor(min(finite_values))
    xtickmax = ceiling(max(finite_values))
    if ((xtickmax==1)&(xtickmin>0)) {
        xtickmin = 0
    }
    shape_numbers = c(0:25, 33:127)
    g[['pointplot']] = ggplot(df_tip_tidy, aes(x=value, y=label, colour=group, shape=group)) +
        geom_jitter(width=0, height=0.05, size=args[['trait_point_size']], alpha=0.7) +
        scale_shape_manual(values=shape_numbers[1:nlevels(df_tip_tidy[['group']])]) +
        scale_color_manual(values=args[['trait_colors']]) +
        scale_x_continuous(breaks=c(xtickmin,(xtickmax-xtickmin)/2, xtickmax)) +
        xlab(xlabel) + 
        theme_linedraw(base_size=font_size) +
        guides(colour=guide_legend(nrow=6, byrow=FALSE)) +
        theme(
            axis.text=element_text(size=font_size),
            axis.title=element_text(size=font_size),
            axis.title.y=element_blank(), 
            axis.text.y=element_blank(), 
            axis.ticks.y=element_blank(),
            legend.title=element_blank(),
            legend.text=element_text(size=font_size),
            legend.position='bottom',
            legend.key.size=unit(0.4, 'lines'), 
            legend.box.just='center',
            rect=element_rect(fill="transparent"),
            plot.margin=unit(args[['margins']], "cm")
        )
    return(g)
}

add_pie_charts = function(gtree, args, tree, branch_color, pie_chart_value_transformation) {
    pcm_prefix = sub('regime$', '', branch_color)
    trait_colors = args[['trait_colors']]
    is_mu_col = startsWith(colnames(gtree[['data']]), paste0(pcm_prefix, 'mu_'))
    regime_col = paste0(pcm_prefix, 'regime')
    shift_col = paste0(pcm_prefix,'is_shift')
    if (any(is_mu_col)) {
        mu_cols = colnames(gtree[['data']])[is_mu_col]
        mu_cols = mu_cols[mu_cols!=paste0(pcm_prefix, 'mu_complementarity')]
        trait_colors = args[['trait_colors']][1:length(mu_cols)]
        names(trait_colors) = mu_cols
        is_shift = (gtree[['data']][[shift_col]]==1)
        shift_node_nums = unlist(c(gtree[['data']][is_shift,'node']))
        pie_node_nums = c(rkftools::get_root_num(tree), shift_node_nums)
        is_pie_node = (gtree[['data']][['node']]%in%pie_node_nums)
        pie_data = gtree[['data']][is_pie_node, c(mu_cols,'node')]
        if (pie_chart_value_transformation=='identity') {
            cat('Pie charts on the tree will be relative values.\n')
            pie_data = pie_data
        } else if (pie_chart_value_transformation=='delog2') {
            cat('Pie charts on the tree will be delog2 relative values.\n')
            pie_data[,mu_cols] = 2**pie_data[,mu_cols]
        } else if (pie_chart_value_transformation=='delog2p1') {
            cat('Pie charts on the tree will be delog2p1 relative values.\n')
            pie_data[,mu_cols] = 2**pie_data[,mu_cols]-1
        } else if (pie_chart_value_transformation=='delog10') {
            cat('Pie charts on the tree will be delog10 relative values.\n')
            pie_data[,mu_cols] = 10**pie_data[,mu_cols]
        } else if (pie_chart_value_transformation=='delog10p1') {
            cat('Pie charts on the tree will be delog10p1 relative values.\n')
            pie_data[,mu_cols] = 10**pie_data[,mu_cols]-1
        }
        pie_data[,mu_cols][pie_data[,mu_cols]<0] = 0
        pie_data[,mu_cols] = t(apply(pie_data[,mu_cols], 1, function(x){x/sum(x)}))
        pie_data = pie_data[order(pie_data[['node']]),]
        rownames(pie_data) = NULL
        if (regime_col %in% colnames(gtree[['data']])) {
            pies = ggtree::nodepie(data=pie_data, cols=1:length(trait_colors), color=trait_colors, alpha=0.7)
            pie_size = 0.25
            gtree = ggtree::inset(tree_view=gtree, insets=pies, width=pie_size, height=pie_size, x='branch', hjust=0, vjust=0)
        }
    }
    return(gtree)
}

append_branch_tiplab_colors = function(g, branch_color, path_species_color_table) {
    tiplab_colors = 'black'
    branch_colors = 'black'
    if (grepl('_regime$', branch_color)) {
        regime_col = branch_color
        if (regime_col %in% colnames(g[['data']])) {
            cat('The regime column was detected. Branches are painted according to the character regimes.\n')
            regime_nos = sort(unique(g[['data']][[regime_col]]))
            regime_colors = c('#000000', colorspace::rainbow_hcl(length(regime_nos)-1, c=100))
            tmp = data.frame(regime=regime_nos, regime_color=regime_colors, stringsAsFactors=FALSE)
            g[['data']] = merge(g[['data']], tmp, by.x=regime_col, by.y='regime')
            tiplab_colors = g[['data']][['regime_color']]
            branch_colors = g[['data']][['regime_color']]
        } else {
            cat('The regime column was not detected. Branches are not painted.\n')
        }
    } else if (branch_color=='species') {
        is_file_available = ifelse(is.null(path_species_color_table), FALSE, file.exists(path_species_color_table))
        if (is_file_available) {
            cat('Species color table was found at:', path_species_color_table, '\n')
            df_spcolor = read.table(path_species_color_table, header=TRUE, sep='\t', stringsAsFactors=FALSE, comment.char='')
        } else {
            cat('Generating species color. File was not found at:', path_species_color_table, '\n')
            has_ubar = grepl('_',g[['data']][['label']])
            spp = g[['data']][['label']][has_ubar]
            spp = sort(unique(get_species_name(spp)))
            num_spp = length(spp)
            # https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
            #qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info[['category']] == 'qual',]
            #col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals[['maxcolors']], rownames(qual_col_pals)))
            #col_vector = col_vector[c(5:length(col_vector), 1:4)]
            col_vector = scales::hue_pal(h=c(0, 360)+15, c=100, l=65, h.start=0, direction = 1)(num_spp)
            colors = col_vector[1:num_spp]
            df_spcolor = data.frame(species=spp, color=col_vector, stringsAsFactors=FALSE)
            write.table(df_spcolor, 'species_color.tsv', sep='\t', row.names=FALSE, quote=FALSE)
        }
        g[['data']][,'species_color'] = '#000000'
        for (i in 1:nrow(df_spcolor)) {
            species_ub = sub(' ', '_', df_spcolor[i,'species'])
            conditions = (g[['data']][['spnode_coverage']]==species_ub)
            g[['data']][conditions,'species_color'] = df_spcolor[i,'color']
        }
        tiplab_colors = g[['data']][['species_color']]
        branch_colors = g[['data']][['species_color']]
    } else if (branch_color=='no') {
        tiplab_colors = 'black'
        branch_colors = 'black'
    } else if (is.numeric(g[['data']][[branch_color]])) {
        values = g[['data']][[branch_color]]
        omega_max = 1 # min(1, max(values))
        values[values>omega_max] = omega_max
        rel_values = values / omega_max
        branch_colors =rgb(0,rel_values,0)
    }
    g[['data']][,'tiplab_color'] = tiplab_colors
    g[['data']][,'branch_color'] = branch_colors
    g[['data']] = g[['data']][order(g[['data']][['node']]),] # If pie charts are messed up, this line might have to be commented out in ggtree 3.x
    return(g)
}

get_df_protein_backbone = function(df_rps, mode='one_line') {
    rows = list()
    idx = 1L
    qacc_vals = as.character(df_rps[['qacc']])
    by_qacc_index = split(seq_len(nrow(df_rps)), qacc_vals, drop = TRUE)
    if (mode=='segmented') {
        for (qacc in names(by_qacc_index)) {
            qidx = by_qacc_index[[qacc]]
            qlen = suppressWarnings(as.integer(df_rps[qidx[1], 'qlen']))
            if (is.na(qlen) || qlen <= 0) {
                next
            }
            qstart = suppressWarnings(as.integer(df_rps[qidx, 'qstart']))
            qend = suppressWarnings(as.integer(df_rps[qidx, 'qend']))
            is_valid = is.finite(qstart) & is.finite(qend)
            if (!any(is_valid)) {
                rows[[idx]] = data.frame(x = 0, xend = qlen, label = qacc)
                idx = idx + 1L
                next
            }
            qstart = qstart[is_valid]
            qend = qend[is_valid]
            qstart = pmax(1L, pmin(qlen, qstart))
            qend = pmax(1L, pmin(qlen, qend))
            keep = (qstart <= qend)
            if (!any(keep)) {
                rows[[idx]] = data.frame(x = 0, xend = qlen, label = qacc)
                idx = idx + 1L
                next
            }
            qstart = qstart[keep]
            qend = qend[keep]
            ord = order(qstart, qend, method = 'radix')
            qstart = qstart[ord]
            qend = qend[ord]

            merged_start = qstart[1]
            merged_end = qend[1]
            merged_rows = list()
            m_idx = 1L
            if (length(qstart) >= 2L) {
                for (i in 2:length(qstart)) {
                    if (qstart[i] <= (merged_end + 1L)) {
                        if (qend[i] > merged_end) {
                            merged_end = qend[i]
                        }
                    } else {
                        merged_rows[[m_idx]] = c(merged_start, merged_end)
                        m_idx = m_idx + 1L
                        merged_start = qstart[i]
                        merged_end = qend[i]
                    }
                }
            }
            merged_rows[[m_idx]] = c(merged_start, merged_end)

            prev_end = 0L
            for (m in seq_len(length(merged_rows))) {
                seg = merged_rows[[m]]
                seg_start = seg[1]
                seg_end = seg[2]
                if (seg_start > (prev_end + 1L)) {
                    rows[[idx]] = data.frame(x = prev_end, xend = seg_start, label = qacc)
                    idx = idx + 1L
                }
                prev_end = seg_end
            }
            if (prev_end < qlen) {
                rows[[idx]] = data.frame(x = prev_end, xend = qlen, label = qacc)
                idx = idx + 1L
            }
        }
    } else if (mode=='one_line') {
        for (qacc in names(by_qacc_index)) {
            qidx = by_qacc_index[[qacc]]
            xend = suppressWarnings(as.numeric(df_rps[qidx[1], 'qlen']))
            if (!is.finite(xend) || xend <= 0) {
                next
            }
            rows[[idx]] = data.frame(x=0, xend=xend, label=qacc)
            idx = idx + 1L
        }
    }
    if (length(rows) == 0) {
        df_pb = data.frame(x = numeric(0), xend = numeric(0), label = character(0))
    } else {
        df_pb = do.call(rbind, rows)
    }
    df_pb[,'xend'] = as.numeric(df_pb[['xend']])
    return(df_pb)
}

prepare_df_rps_for_plot = function(df_rpsblast, df_tip) {
    if (is.null(df_rpsblast) || nrow(df_rpsblast)==0) {
        return(data.frame())
    }
    if (!('qacc' %in% colnames(df_rpsblast))) {
        cat('The qacc column was not found in the rpsblast table.\n')
        return(data.frame())
    }
    if (!('sacc' %in% colnames(df_rpsblast))) {
        df_rpsblast[,'sacc'] = NA_character_
    }
    if (!('stitle' %in% colnames(df_rpsblast))) {
        df_rpsblast[,'stitle'] = as.character(df_rpsblast[['sacc']])
    }
    required_num_cols = c('qlen', 'slen', 'qstart', 'qend')
    missing_num_cols = required_num_cols[!(required_num_cols %in% colnames(df_rpsblast))]
    if (length(missing_num_cols)>0) {
        cat('Required rpsblast columns are missing:', paste(missing_num_cols, collapse=', '), '\n')
        return(data.frame())
    }
    df_rps = df_rpsblast[(df_rpsblast[['qacc']] %in% df_tip[['label']]),,drop=FALSE]
    if (nrow(df_rps)==0) {
        return(df_rps)
    }
    df_rps[,'sacc_original'] = df_rps[['sacc']]
    df_rps[,'stitle_original'] = df_rps[['stitle']]
    df_rps[,'stitle'] = sub(',', '---PLACEHOLDER---', df_rps[['stitle']])
    df_rps[,'stitle'] = sub(',.*', '', df_rps[['stitle']])
    df_rps[,'stitle'] = sub('.*---PLACEHOLDER---', '', df_rps[['stitle']])
    df_rps[,'sacc'] = df_rps[['stitle']]
    df_rps[,'qlen'] = suppressWarnings(as.numeric(sub(',', '', as.character(df_rps[['qlen']]))))
    df_rps[,'slen'] = suppressWarnings(as.numeric(sub(',', '', as.character(df_rps[['slen']]))))
    df_rps[,'qstart'] = suppressWarnings(as.numeric(sub(',', '', as.character(df_rps[['qstart']]))))
    df_rps[,'qend'] = suppressWarnings(as.numeric(sub(',', '', as.character(df_rps[['qend']]))))
    df_rps[,'label'] = factor(df_rps[['qacc']], levels=df_tip[['label']])
    df_rps[,'ymin'] = as.numeric(df_rps[['label']]) - 0.375
    df_rps[,'ymax'] = as.numeric(df_rps[['label']]) + 0.375
    return(df_rps)
}

get_domain_fill_colors = function(df_domain) {
    if (is.null(df_domain) || nrow(df_domain)==0 || !('sacc' %in% colnames(df_domain))) {
        return(c())
    }
    domain_levels = levels(df_domain[['sacc']])
    if (is.null(domain_levels) || length(domain_levels)==0) {
        domain_levels = unique(as.character(df_domain[['sacc']]))
    }
    num_domain = length(domain_levels)
    if (num_domain==0) {
        return(c())
    }
    defaultW <- getOption("warn")
    options(warn = -1)
    mycolors = colorRampPalette(RColorBrewer::brewer.pal(12, "Set1"))(num_domain)
    options(warn = defaultW)
    names(mycolors) = domain_levels
    return(mycolors)
}

map_trimmed_to_untrimmed_nongap_index = function(seq_trim, seq_untrim, gap_code='04') {
  seq_trim = as.character(seq_trim)
  seq_untrim = as.character(seq_untrim)
  out = rep(NA_integer_, length(seq_trim))
  if (length(seq_trim)==0 || length(seq_untrim)==0) {
    return(out)
  }
  trim_non_gap_idx = which(seq_trim != gap_code)
  untrim_non_gap_idx = which(seq_untrim != gap_code)
  if (length(trim_non_gap_idx)==0 || length(untrim_non_gap_idx)==0) {
    return(out)
  }
  trim_tokens = seq_trim[trim_non_gap_idx]
  untrim_tokens = seq_untrim[untrim_non_gap_idx]
  u_ptr = 1L
  for (i in seq_along(trim_tokens)) {
    tok = trim_tokens[i]
    while (u_ptr <= length(untrim_tokens) && untrim_tokens[u_ptr] != tok) {
      u_ptr = u_ptr + 1L
    }
    if (u_ptr > length(untrim_tokens)) {
      return(rep(NA_integer_, length(seq_trim)))
    }
    out[trim_non_gap_idx[i]] = u_ptr
    u_ptr = u_ptr + 1L
  }
  return(out)
}

add_alignment_column <- function(g, args, seqs = NULL, df_rpsblast = NULL, seqs_untrim = NULL) {
  cat(as.character(Sys.time()), 'Adding alignment column.\n')
  
  if (is.null(seqs)) {
    cat('Input alignment is empty. Alignment column will not be added.\n')
    return(g)
  }
  
  df_tip <- get_df_tip(g[['tree']])
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

add_heatmap_column = function(g, args, df_trait) {
    cat(as.character(Sys.time()), 'Adding heatmap column.\n')
    if ((is.null(ncol(df_trait)))|(ncol(df_trait)==0)) {
        cat('df_trait is emply. Heatmap panel will not be added.\n')
        return(g)
    }
    font_size = args[['font_size']]
    if (any(grepl('^pointplot', unlist(args[grep("^panel", names(args))])))) {
        trait_colors = args[['trait_colors']] # TODO Not matching to pointplot colors.
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
    if (any(apply(df_tip[,colnames(df_trait)], 1, function(x){all(is.na(x))}))) {
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
        ggplot2::labs(fill="Expression") +
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
        if (grepl(',', df_tip[i,'intron_positions'])) {
            intron_val = df_tip[i,'intron_positions']
            intron_pos = strsplit(as.character(intron_val), ',')[[1]]
        } else {
            intron_pos = df_tip[i,'intron_positions']
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
    df_tip2 = df_tip[,c('y',newcols)] %>% tidyr::gather('key','value',-y)
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

enhance_branch_table = function(b, args, event_method=c('auto','generax','species_overlap')) {
    if ('mapdnds_omega' %in% colnames(b)) {
        b[is.na(b[['mapdnds_omega']]),'mapdnds_omega'] = 0
    }
    b[['branch_category']] = 'S'
    if ((event_method=='generax')|((event_method=='auto')&('generax_event' %in% colnames(b)))) {
        cat('generax_event_parent will be used to define branch events.\n')
        b[['is_parent_dup']] = (b[['generax_event_parent']]=='D')
        b[['is_parent_transfer']] = (b[['generax_event_parent']]=='H')
    } else if (event_method=='species_overlap') {
        cat('so_event_parent will be used to define branch events.\n')
        b[['is_parent_dup']] = (b[['so_event_parent']]=='D')
        b[['is_parent_transfer']] = FALSE
    } else {
        cat('Non-supported --event_method specified. Exiting.\n')
        quit()
    }
    if ('delta_intron_present' %in% colnames(b)) {
        cat('delta_intron_present will be used to define R branches.\n')
        cat('Nodes will be categorized into: S, D, R\n')
        b = attach_neighbor_stats(df=b, neighbor='sister', columns='delta_intron_present')
        b[['is_retrotransposition']] = (b[['delta_intron_present']] <= args[['max_delta_intron_present']])
        b[['is_retrotransposition']][is.na(b[['is_retrotransposition']])] = FALSE
        b[['is_sister_retrotransposition']] = (b[['sister_delta_intron_present']] <= args[['max_delta_intron_present']])
        b[['is_sister_retrotransposition']][is.na(b[['is_sister_retrotransposition']])] = FALSE
        b[['is_lower_delta_intron_present']] = (b[['delta_intron_present']] <= b[['sister_delta_intron_present']])    
        b[(b[['is_parent_dup']]),'branch_category'] = 'D'
        b[(b[['is_parent_dup']])&(b[['is_retrotransposition']])&(!b[['is_sister_retrotransposition']]),'branch_category'] = 'R'
    } else {
        cat('Nodes will be categorized into: S, D', '\n')
        b[(b[['is_parent_dup']]),'branch_category'] = 'D'
    }
    if (sum(b[['is_parent_transfer']])>0) {
        b[(b[['is_parent_transfer']]),'branch_category'] = 'H'
    }
    
    b = attach_neighbor_stats(df=b, neighbor='child1', columns='branch_category')
    b = attach_neighbor_stats(df=b, neighbor='child2', columns='branch_category')
    b[is.na(b[['child1_branch_category']]),'child1_branch_category'] = ''
    b[is.na(b[['child2_branch_category']]),'child2_branch_category'] = ''
    b[['node_category']] = ''
    child_cat = b[,c('child1_branch_category','child2_branch_category'), drop = FALSE]
    num_children_S = rowSums(child_cat == 'S')
    num_children_D = rowSums(child_cat == 'D')
    num_children_R = rowSums(child_cat == 'R')
    num_children_H = rowSums(child_cat == 'H')
    b[(num_children_S>=1),'node_category'] = 'S'
    b[(num_children_D>=1),'node_category'] = 'D'
    b[((num_children_D==1)&(num_children_R==1)),'node_category'] = 'R'
    b[(num_children_H>=1),'node_category'] = 'H'
    cat('detected node categories:', unique(b[['node_category']][b[['node_category']]!='NaN']), '\n')

    # Support values
    support_cols = colnames(b)[grepl('^support_', colnames(b))]
    root_nlabel = max(b[['branch_id']])
    subroot_nlabels = b[(b[['parent']]==root_nlabel),'branch_id']
    for (col in support_cols) {
        root_support = b[(b[['parent']]==-999),col]
        if (!is.na(root_support)) {
            cat('Warning:', col, ': support value at root node is not NA (', root_support, '). Masking.\n')
            cat('Non-NA root support is no problem if this is a pruned tree. If not, check the data carefully.\n')
            b[(b[['parent']]==-999),col] = NA
        }
        is_reconciled_branch = ((is.na(b[[col]]))&(b[['parent']]!=-999))
        b[is_reconciled_branch,col] = '-'
        # Suppress the support value of one of the subroot branches
        subroot_idx = which(b[['branch_id']] %in% subroot_nlabels)
        if (length(subroot_idx) > 0) {
            b[subroot_idx[length(subroot_idx)], col] = NA
        }
    }

    # dNdS
    omega_method = 'mapdnds'
    omega_column = paste0(omega_method, '_omega')
    if (omega_column %in% colnames(b)) {
        b = attach_neighbor_stats(df=b, neighbor='sister', columns=omega_column)
        b[['is_higher_dnds']] = (b[[omega_column]]>=b[[paste0('sister_', omega_column)]])
    }

    # csubst
    if ('csubst_is_foreground' %in% colnames(b)) {
        b[,'csubst_is_foreground'] = ifelse(b[['csubst_is_foreground']] %in% c('Y','yes'), TRUE, FALSE)
        b[,'csubst_is_marginal'] = ifelse(b[['csubst_is_marginal']] %in% c('Y','yes'), TRUE, FALSE)
    }

    # Branch width
    b[['branch_thickness']] = 0.25
    #b[(b[['is_higher_dnds']]&b[['is_parent_dup']]),'branch_thickness'] = 1

    #img_info = get_img_info(b, args[['dir_phylopic']])

    # sorting
    b = b[order(b[['branch_id']]),]
    rownames(b) = NULL
    return(b)
}

add_gene_cluster_membership = function(df_tip, max_bp_membership) {
    df_tip[,'cluster_membership'] = ''
    df_tip[,'species'] = ''
    is_annotated_leaf = (df_tip[['so_event']]=='L')
    is_annotated_leaf = is_annotated_leaf & (!is.na(df_tip[['start']]))
    is_annotated_leaf = is_annotated_leaf & (!is.na(df_tip[['end']]))
    spp = unique(df_tip[is_annotated_leaf, 'taxon', drop = FALSE][['taxon']])
    for (sp in spp) {
        sp_ub = sub(' ', '_', sp)
        gene_membership_counter = 0
        is_sp = ((df_tip[['so_event']]=='L')&(df_tip[['taxon']]==sp))
        df_tip[is_sp,'species'] = sp_ub
        chromosome_names = unique(df_tip[is_sp, 'chromosome', drop = FALSE][['chromosome']])
        chromosome_names = chromosome_names[chromosome_names!='']
        for (chromosome_name in chromosome_names) {
            is_chromosome = ((is_sp)&(df_tip[['chromosome']]==chromosome_name))
            is_chromosome[is.na(is_chromosome)] = FALSE
            n_chr = sum(is_chromosome)
            if (n_chr == 0) {
                next
            }
            b_cr = data.frame(df_tip[is_chromosome,], stringsAsFactors=FALSE)
            b_cr_start_new = apply(b_cr[,c('start','end')], 1, min)
            b_cr_end_new = apply(b_cr[,c('start','end')], 1, max)
            b_cr[,'start'] = b_cr_start_new
            b_cr[,'end'] = b_cr_end_new
            b_cr = b_cr[order(b_cr[['start']]), , drop=FALSE]
            rownames(b_cr) = NULL

            # Start a new cluster at each chromosome, then split by intergenic gaps.
            gene_membership_counter = gene_membership_counter + 1
            if (nrow(b_cr) == 1) {
                cluster_ids = gene_membership_counter
            } else {
                intergenic = b_cr[['start']][-1] - b_cr[['end']][-nrow(b_cr)]
                is_new_cluster = c(TRUE, intergenic > max_bp_membership)
                cluster_ids = (gene_membership_counter - 1L) + cumsum(is_new_cluster)
                gene_membership_counter = max(cluster_ids)
            }
            b_cr[['cluster_membership']] = paste0(sp_ub, '_', cluster_ids)
            df_tip[match(as.character(b_cr[['label']]), as.character(df_tip[['label']])), 'cluster_membership'] = b_cr[['cluster_membership']]
        }
    }
    return(df_tip)
}

add_cluster_membership_column = function(g, args, gname, max_bp_membership){
    cat(as.character(Sys.time()), 'Adding gene cluster membership column.\n')
    if (!'start'%in%colnames(g[['tree']][['data']])) {
        cat('Scaffold positions are not available. cluster_membership plot will not be displayed.\n')
        return(g)
    }   
    
    df_tip = get_df_tip(g[['tree']])
    df_tip = add_gene_cluster_membership(df_tip=df_tip, max_bp_membership=max_bp_membership)
    font_size = args[['font_size']]

    duplicated_membership = duplicated(df_tip[['cluster_membership']])
    duplicated_membership = duplicated_membership & (!is.na(df_tip[['cluster_membership']]))
    duplicated_membership = duplicated_membership & (df_tip[['cluster_membership']] != '')
    multimember_clusters = unique(df_tip[['cluster_membership']][duplicated_membership])
    is_multimember = (df_tip[['cluster_membership']]%in%multimember_clusters)
    is_na = (is.na(df_tip[['cluster_membership']]))|(df_tip[['cluster_membership']]=='')
    df_single = df_tip[(!is_multimember)&(!is_na),]
    df_multi = df_tip[(is_multimember)&(!is_na),]

    g[[gname]] = ggplot(mapping=aes(x=species, y=y)) +
        geom_blank(data=df_tip, aes(y=label)) + 
        geom_point(data=df_single, color='gray90', alpha=1) +
        geom_line(data=df_multi, mapping=aes(group=cluster_membership, color=cluster_membership), alpha=0.5) +
        geom_point(data=df_multi, mapping=aes(color=cluster_membership), alpha=1) +
        #aplot::ylim2(gg=g[['tree']]) +
        #ylim(c(1, max(df_tip[['y']]))) +
        xlab(paste0('Gene cluster\nmembership\n(max dist =\n', max_bp_membership, ' bp)')) +
        theme_void() +
        guides(
            fill=guide_legend(title=NULL, nrow=6, byrow=FALSE)
        ) +
        coord_cartesian(clip = "off") +
        theme(
            axis.text=element_blank(),
            axis.title.x=element_text(size=args[['font_size']]),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.text.x=element_blank(),
            legend.title=element_blank(),
            legend.text=element_blank(),
            legend.position="none",
            legend.key.size=unit(0.4, 'lines'),
            legend.box.just='center',
            plot.margin=unit(args[['margins']], "cm")
        )
  return(g)
}

gg_color_hue = function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

add_motif_colors = function(df_fimo, min_count=2) {
    motif_keys = as.character(df_fimo[['motif_altid']])
    motif_tab = sort(table(motif_keys), decreasing=TRUE)
    df_count = data.frame(
        altid = names(motif_tab),
        count = as.integer(motif_tab),
        stringsAsFactors = FALSE
    )
    if (nrow(df_count) == 0) {
        df_fimo[,'count'] = rep(NA_integer_, nrow(df_fimo))
        df_fimo[,'color_id'] = factor(rep(NA_integer_, nrow(df_fimo)))
        return(df_fimo)
    }
    df_count[,'color_id'] = seq_len(nrow(df_count))
    df_count[(df_count[['count']]<min_count),'color_id'] = 0
    color_ids = unname(df_count[['color_id']][match(motif_keys, df_count[['altid']])])
    df_fimo[,'count'] = unname(df_count[['count']][match(motif_keys, df_count[['altid']])])
    df_fimo[,'color_id'] = factor(color_ids, levels=sort(unique(df_count[,'color_id'])))
    return(df_fimo)
}

add_motif_counts = function(df_fimo) {
    motif_keys = as.character(df_fimo[['motif_altid']])
    motif_tab = table(motif_keys)
    df_fimo[,'count'] = as.integer(unname(motif_tab[motif_keys]))
    return(df_fimo)
}

get_df_polygon = function(df_fimo, min_count, max_count, multiple_connection='align_from_tss', alpha_max=0.5, ncpu=1, verbose=TRUE) {
    empty_polygon_df = function() {
        data.frame(
            x = numeric(0),
            y = numeric(0),
            polygon_id = factor(),
            motif_altid = factor(),
            alpha_value = numeric(0)
        )
    }
    n = nrow(df_fimo)
    if (n == 0L) {
        if (verbose) {
            cat('Number of motifs with displayed connections: 0 \n')
        }
        return(empty_polygon_df())
    }
    motif_raw = df_fimo[['motif_altid']]
    motif_vals = as.character(motif_raw)
    valid_motif = !is.na(motif_raw) & nzchar(motif_vals)
    motif_tab = table(motif_vals[valid_motif])
    keep_altids = names(motif_tab)[(motif_tab >= min_count) & (motif_tab <= max_count)]
    if (verbose) {
        cat('Number of motifs with displayed connections:', length(keep_altids), '\n')
    }
    if (length(keep_altids) == 0) {
        return(empty_polygon_df())
    }

    y_num = suppressWarnings(as.numeric(df_fimo[['y']]))
    start_num = suppressWarnings(as.numeric(df_fimo[['start']]))
    end_num = suppressWarnings(as.numeric(df_fimo[['end']]))
    ymin_num = suppressWarnings(as.numeric(df_fimo[['ymin']]))
    ymax_num = suppressWarnings(as.numeric(df_fimo[['ymax']]))
    keep_map = match(motif_vals, keep_altids, nomatch = 0L)
    valid_row = valid_motif & is.finite(y_num) & is.finite(start_num) & is.finite(end_num) &
        is.finite(ymin_num) & is.finite(ymax_num) & (keep_map > 0L)
    if (!any(valid_row)) {
        return(empty_polygon_df())
    }
    num_max = length(unique(y_num[valid_row]))
    if (num_max == 0L) {
        return(empty_polygon_df())
    }
    ncpu = suppressWarnings(as.integer(ncpu))
    if (!is.finite(ncpu) || is.na(ncpu) || (ncpu < 1L)) {
        ncpu = 1L
    }

    idx_keep = which(valid_row)
    idx_by_altid = split(idx_keep, motif_vals[idx_keep], drop = TRUE)
    build_altid_polygon_parts = function(altid) {
        idx_altid = idx_by_altid[[altid]]
        if (is.null(idx_altid) || (length(idx_altid) < 2L)) {
            return(NULL)
        }
        ord = order(y_num[idx_altid], start_num[idx_altid], end_num[idx_altid], method = 'radix')
        idx_altid = idx_altid[ord]
        n_altid = length(idx_altid)
        if (n_altid < 2L) {
            return(NULL)
        }
        yvals = y_num[idx_altid]
        starts = start_num[idx_altid]
        ends = end_num[idx_altid]
        ymins = ymin_num[idx_altid]
        ymaxs = ymax_num[idx_altid]
        y_change = c(TRUE, yvals[-1L] != yvals[-n_altid])
        group_starts = which(y_change)
        n_group = length(group_starts)
        if (n_group < 2L) {
            return(NULL)
        }
        group_ends = c(group_starts[-1L] - 1L, n_altid)
        alpha_value = (n_group / num_max) * alpha_max
        x_parts_local = list()
        y_parts_local = list()
        id_parts_local = list()
        motif_parts_local = list()
        alpha_parts_local = list()
        part_idx_local = 1L
        counter_local = 0L
        append_local = function(x1, x2, x3, x4, y1, y2, y3, y4) {
            n_poly = length(x1)
            if (n_poly == 0L) {
                return()
            }
            total_len = n_poly * 4L
            pos = seq.int(1L, total_len, by = 4L)
            x_vec = numeric(total_len)
            y_vec = numeric(total_len)
            x_vec[pos] = x1
            x_vec[pos + 1L] = x2
            x_vec[pos + 2L] = x3
            x_vec[pos + 3L] = x4
            y_vec[pos] = y1
            y_vec[pos + 1L] = y2
            y_vec[pos + 2L] = y3
            y_vec[pos + 3L] = y4
            poly_ids = counter_local + seq_len(n_poly) - 1L
            x_parts_local[[part_idx_local]] <<- x_vec
            y_parts_local[[part_idx_local]] <<- y_vec
            id_parts_local[[part_idx_local]] <<- rep(poly_ids, each = 4L)
            motif_parts_local[[part_idx_local]] <<- rep(altid, total_len)
            alpha_parts_local[[part_idx_local]] <<- rep(alpha_value, total_len)
            part_idx_local <<- part_idx_local + 1L
            counter_local <<- counter_local + n_poly
        }
        for (k in seq_len(n_group - 1L)) {
            current_ind = seq.int(group_starts[k], group_ends[k])
            next_ind = seq.int(group_starts[k + 1L], group_ends[k + 1L])
            n_current = length(current_ind)
            n_next = length(next_ind)
            if ((n_current == 0L) || (n_next == 0L)) {
                next
            }
            if ((n_current == 1L) || (n_next == 1L)) {
                if (multiple_connection == 'align_from_tss') {
                    current_ind = tail(current_ind, 1L)
                    next_ind = tail(next_ind, 1L)
                    n_current = length(current_ind)
                    n_next = length(next_ind)
                }
                if (n_current == 1L) {
                    ci_vec = rep(current_ind, n_next)
                    ni_vec = next_ind
                } else if (n_next == 1L) {
                    ci_vec = current_ind
                    ni_vec = rep(next_ind, n_current)
                } else {
                    ci_vec = rep(current_ind, each = n_next)
                    ni_vec = rep(next_ind, times = n_current)
                }
                if (length(ci_vec) == 0L) {
                    next
                }
                append_local(
                    starts[ci_vec], ends[ci_vec], ends[ni_vec], starts[ni_vec],
                    ymaxs[ci_vec], ymaxs[ci_vec], ymins[ni_vec], ymins[ni_vec]
                )
            } else {
                if (multiple_connection == 'align_from_tss') {
                    min_num_ind = min(n_current, n_next)
                    if (min_num_ind == 0L) {
                        next
                    }
                    current_ind2 = tail(current_ind, min_num_ind)
                    next_ind2 = tail(next_ind, min_num_ind)
                    append_local(
                        starts[current_ind2], ends[current_ind2], ends[next_ind2], starts[next_ind2],
                        ymaxs[current_ind2], ymaxs[current_ind2], ymins[next_ind2], ymins[next_ind2]
                    )
                } else if (multiple_connection == 'all') {
                    current_mid = (starts[current_ind] + ends[current_ind]) / 2
                    next_mid = (starts[next_ind] + ends[next_ind]) / 2
                    node_x = mean(c(current_mid, next_mid))
                    current_size = ends[current_ind] - starts[current_ind]
                    next_size = ends[next_ind] - starts[next_ind]
                    node_width = mean(c(current_size, next_size))
                    node_y = mean(c(mean(yvals[current_ind]), mean(yvals[next_ind])))
                    if (is.finite(node_y) && (node_y %% 1 == 0)) {
                        node_y = node_y + 0.5
                    }
                    node_x_start = node_x - (node_width / 2)
                    node_x_end = node_x + (node_width / 2)
                    append_local(
                        starts[current_ind], ends[current_ind], rep(node_x_end, n_current), rep(node_x_start, n_current),
                        ymaxs[current_ind], ymaxs[current_ind], rep(node_y, n_current), rep(node_y, n_current)
                    )
                    append_local(
                        starts[next_ind], ends[next_ind], rep(node_x_end, n_next), rep(node_x_start, n_next),
                        ymins[next_ind], ymins[next_ind], rep(node_y, n_next), rep(node_y, n_next)
                    )
                }
            }
        }
        if (part_idx_local == 1L) {
            return(NULL)
        }
        return(list(
            x = unlist(x_parts_local, use.names = FALSE),
            y = unlist(y_parts_local, use.names = FALSE),
            id = unlist(id_parts_local, use.names = FALSE),
            motif = unlist(motif_parts_local, use.names = FALSE),
            alpha = unlist(alpha_parts_local, use.names = FALSE),
            n_poly = as.integer(counter_local)
        ))
    }
    use_parallel = (ncpu > 1L) && (length(keep_altids) > 1L) && (.Platform$OS.type != 'windows')
    if (use_parallel) {
        motif_results = parallel::mclapply(keep_altids, build_altid_polygon_parts, mc.cores = ncpu)
        keep_result = !vapply(motif_results, is.null, logical(1))
        motif_results = motif_results[keep_result]
        if (length(motif_results) == 0L) {
            return(empty_polygon_df())
        }
        seg_len = vapply(motif_results, function(x) length(x[['x']]), integer(1))
        total_len = sum(seg_len)
        x_all = numeric(total_len)
        y_all = numeric(total_len)
        id_all = integer(total_len)
        motif_all = character(total_len)
        alpha_all = numeric(total_len)
        write_pos = 1L
        poly_offset = 0L
        for (ri in seq_along(motif_results)) {
            seg_n = seg_len[ri]
            idx = write_pos:(write_pos + seg_n - 1L)
            res = motif_results[[ri]]
            x_all[idx] = res[['x']]
            y_all[idx] = res[['y']]
            id_all[idx] = as.integer(res[['id']]) + poly_offset
            motif_all[idx] = res[['motif']]
            alpha_all[idx] = res[['alpha']]
            poly_offset = poly_offset + as.integer(res[['n_poly']])
            write_pos = write_pos + seg_n
        }
        df_polygon = data.frame(
            x = x_all,
            y = y_all,
            polygon_id = factor(id_all),
            motif_altid = factor(motif_all),
            alpha_value = alpha_all,
            stringsAsFactors = FALSE
        )
        return(df_polygon)
    }

    x_parts = list()
    y_parts = list()
    id_parts = list()
    motif_parts = list()
    alpha_parts = list()
    part_idx = 1L
    counter = 0L

    append_polygon_batch = function(x1, x2, x3, x4, y1, y2, y3, y4, altid, alpha_value) {
        n_poly = length(x1)
        if (n_poly == 0L) {
            return()
        }
        total_len = n_poly * 4L
        pos = seq.int(1L, total_len, by = 4L)
        x_vec = numeric(total_len)
        y_vec = numeric(total_len)
        x_vec[pos] = x1
        x_vec[pos + 1L] = x2
        x_vec[pos + 2L] = x3
        x_vec[pos + 3L] = x4
        y_vec[pos] = y1
        y_vec[pos + 1L] = y2
        y_vec[pos + 2L] = y3
        y_vec[pos + 3L] = y4
        poly_ids = counter + seq_len(n_poly) - 1L
        x_parts[[part_idx]] <<- x_vec
        y_parts[[part_idx]] <<- y_vec
        id_parts[[part_idx]] <<- rep(poly_ids, each = 4L)
        motif_parts[[part_idx]] <<- rep(altid, total_len)
        alpha_parts[[part_idx]] <<- rep(alpha_value, total_len)
        part_idx <<- part_idx + 1L
        counter <<- counter + n_poly
    }

    for (altid in keep_altids) {
        idx_altid = idx_by_altid[[altid]]
        if (is.null(idx_altid) || (length(idx_altid) < 2L)) {
            next
        }
        ord = order(y_num[idx_altid], start_num[idx_altid], end_num[idx_altid], method = 'radix')
        idx_altid = idx_altid[ord]
        n_altid = length(idx_altid)
        if (n_altid < 2L) {
            next
        }
        yvals = y_num[idx_altid]
        starts = start_num[idx_altid]
        ends = end_num[idx_altid]
        ymins = ymin_num[idx_altid]
        ymaxs = ymax_num[idx_altid]
        y_change = c(TRUE, yvals[-1L] != yvals[-n_altid])
        group_starts = which(y_change)
        n_group = length(group_starts)
        if (n_group < 2L) {
            next
        }
        group_ends = c(group_starts[-1L] - 1L, n_altid)
        alpha_value = (n_group / num_max) * alpha_max

        for (k in seq_len(n_group - 1L)) {
            current_ind = seq.int(group_starts[k], group_ends[k])
            next_ind = seq.int(group_starts[k + 1L], group_ends[k + 1L])
            n_current = length(current_ind)
            n_next = length(next_ind)
            if ((n_current == 0L) || (n_next == 0L)) {
                next
            }
            if ((n_current == 1L) || (n_next == 1L)) {
                if (multiple_connection == 'align_from_tss') {
                    current_ind = tail(current_ind, 1L)
                    next_ind = tail(next_ind, 1L)
                    n_current = length(current_ind)
                    n_next = length(next_ind)
                }
                if (n_current == 1L) {
                    ci_vec = rep(current_ind, n_next)
                    ni_vec = next_ind
                } else if (n_next == 1L) {
                    ci_vec = current_ind
                    ni_vec = rep(next_ind, n_current)
                } else {
                    ci_vec = rep(current_ind, each = n_next)
                    ni_vec = rep(next_ind, times = n_current)
                }
                if (length(ci_vec) == 0L) {
                    next
                }
                append_polygon_batch(
                    starts[ci_vec], ends[ci_vec], ends[ni_vec], starts[ni_vec],
                    ymaxs[ci_vec], ymaxs[ci_vec], ymins[ni_vec], ymins[ni_vec],
                    altid, alpha_value
                )
            } else {
                if (multiple_connection == 'align_from_tss') {
                    min_num_ind = min(n_current, n_next)
                    if (min_num_ind == 0L) {
                        next
                    }
                    current_ind2 = tail(current_ind, min_num_ind)
                    next_ind2 = tail(next_ind, min_num_ind)
                    append_polygon_batch(
                        starts[current_ind2], ends[current_ind2], ends[next_ind2], starts[next_ind2],
                        ymaxs[current_ind2], ymaxs[current_ind2], ymins[next_ind2], ymins[next_ind2],
                        altid, alpha_value
                    )
                } else if (multiple_connection == 'all') {
                    current_mid = (starts[current_ind] + ends[current_ind]) / 2
                    next_mid = (starts[next_ind] + ends[next_ind]) / 2
                    node_x = mean(c(current_mid, next_mid))
                    current_size = ends[current_ind] - starts[current_ind]
                    next_size = ends[next_ind] - starts[next_ind]
                    node_width = mean(c(current_size, next_size))
                    node_y = mean(c(mean(yvals[current_ind]), mean(yvals[next_ind])))
                    if (is.finite(node_y) && (node_y %% 1 == 0)) {
                        node_y = node_y + 0.5
                    }
                    node_x_start = node_x - (node_width / 2)
                    node_x_end = node_x + (node_width / 2)
                    append_polygon_batch(
                        starts[current_ind], ends[current_ind], rep(node_x_end, n_current), rep(node_x_start, n_current),
                        ymaxs[current_ind], ymaxs[current_ind], rep(node_y, n_current), rep(node_y, n_current),
                        altid, alpha_value
                    )
                    append_polygon_batch(
                        starts[next_ind], ends[next_ind], rep(node_x_end, n_next), rep(node_x_start, n_next),
                        ymins[next_ind], ymins[next_ind], rep(node_y, n_next), rep(node_y, n_next),
                        altid, alpha_value
                    )
                }
            }
        }
    }

    if (part_idx == 1L) {
        return(empty_polygon_df())
    }
    df_polygon = data.frame(
        x = unlist(x_parts, use.names = FALSE),
        y = unlist(y_parts, use.names = FALSE),
        polygon_id = factor(unlist(id_parts, use.names = FALSE)),
        motif_altid = factor(unlist(motif_parts, use.names = FALSE)),
        alpha_value = unlist(alpha_parts, use.names = FALSE),
        stringsAsFactors = FALSE
    )
    return(df_polygon)
}

merge_overlap_group = function(df_fimo, merge_level='TF', ncpu=1, verbose=TRUE) {
    if (merge_level=='TF') {
        altid_col = 'motif_altid'
    } else if (merge_level=='TFfamily') {
        altid_col = 'merged_altid'
        motif_raw = df_fimo[['motif_altid']]
        motif_chr = as.character(motif_raw)
        uniq_motif = unique(motif_chr)
        uniq_family = sub('([a-zA-Z])[-_]{0,1}[0-9]{1,3}$', '\\1', uniq_motif)
        mapped_family = uniq_family[match(motif_chr, uniq_motif)]
        mapped_family[is.na(motif_raw)] = NA_character_
        df_fimo[,altid_col] = mapped_family
    }
    df_fimo = add_complete_overlap_groups(df_fimo, altid_col, ncpu = ncpu, verbose = verbose)
    n = nrow(df_fimo)
    if (n == 0) {
        return(df_fimo[0, , drop = FALSE])
    }
    labels_chr = as.character(df_fimo[['label']])
    altid_chr = as.character(df_fimo[[altid_col]])
    starts = suppressWarnings(as.numeric(df_fimo[['start']]))
    ends = suppressWarnings(as.numeric(df_fimo[['end']]))
    swap_idx = is.finite(starts) & is.finite(ends) & (starts > ends)
    if (any(swap_idx)) {
        tmp_swap = starts[swap_idx]
        starts[swap_idx] = ends[swap_idx]
        ends[swap_idx] = tmp_swap
    }
    sort_order = order(labels_chr, altid_chr, starts, ends, method='radix', na.last=TRUE)
    labels_chr = labels_chr[sort_order]
    altid_chr = altid_chr[sort_order]
    starts = starts[sort_order]
    ends = ends[sort_order]
    y_vals = suppressWarnings(as.numeric(df_fimo[['y']]))[sort_order]
    strand_vals = as.character(df_fimo[['strand']])
    strand_vals = strand_vals[sort_order]
    finite_start = is.finite(starts)
    finite_end = is.finite(ends)

    label_levels = unique(labels_chr)
    altid_levels = unique(altid_chr)
    label_id = match(labels_chr, label_levels)
    altid_id = match(altid_chr, altid_levels)
    label_key = label_id
    altid_key = altid_id
    label_key[is.na(label_key)] = 0L
    altid_key[is.na(altid_key)] = 0L
    same_key = c(FALSE, (label_key[-1] == label_key[-n]) & (altid_key[-1] == altid_key[-n]))
    key_new = c(TRUE, !same_key[-1])
    key_starts = which(key_new)
    key_ends = c(key_starts[-1] - 1L, n)
    end_scan = ends
    end_scan[!finite_end] = -Inf
    prev_max_end = rep(-Inf, n)
    for (k in seq_len(length(key_starts))) {
        s = key_starts[k]
        e = key_ends[k]
        if (s < e) {
            prev_max_end[(s + 1L):e] = cummax(end_scan[s:(e - 1L)])
        }
    }
    new_group = key_new | !finite_start | (starts > prev_max_end)
    num_groups = sum(new_group)

    y_clean = y_vals
    y_clean[!is.finite(y_clean)] = Inf
    start_clean = starts
    start_clean[!finite_start] = Inf
    end_clean = ends
    end_clean[!finite_end] = -Inf
    strand_mask = integer(n)
    strand_not_na = !is.na(strand_vals)
    strand_mask[strand_not_na & (strand_vals == '+')] = 1L
    strand_mask[strand_not_na & (strand_vals == '-')] = 2L
    strand_mask[strand_not_na & ((strand_vals == '+-') | (strand_vals == '-+'))] = 3L
    # Fast bitwise-union lookup for masks in {0,1,2,3}.
    mask_or = matrix(c(
        0L, 1L, 2L, 3L,
        1L, 1L, 3L, 3L,
        2L, 3L, 2L, 3L,
        3L, 3L, 3L, 3L
    ), nrow = 4L, byrow = TRUE)
    out_label_id = rep(NA_integer_, num_groups)
    out_altid_id = rep(NA_integer_, num_groups)
    out_group = seq_len(num_groups)
    out_y = rep(NA_real_, num_groups)
    out_start = rep(NA_real_, num_groups)
    out_end = rep(NA_real_, num_groups)
    out_strand = rep(NA_character_, num_groups)
    current_group = 1L
    current_label_id = label_id[1L]
    current_altid_id = altid_id[1L]
    current_y_min = y_clean[1L]
    current_start_min = start_clean[1L]
    current_end_max = end_clean[1L]
    current_strand_mask = strand_mask[1L]
    current_first_raw_idx = if (strand_not_na[1L]) 1L else 0L

    if (n > 1L) {
        for (i in 2:n) {
            if (new_group[i]) {
                out_label_id[current_group] = current_label_id
                out_altid_id[current_group] = current_altid_id
                if (is.finite(current_y_min)) {
                    out_y[current_group] = current_y_min
                }
                if (is.finite(current_start_min)) {
                    out_start[current_group] = current_start_min
                }
                if (is.finite(current_end_max)) {
                    out_end[current_group] = current_end_max
                }
                if (current_strand_mask == 3L) {
                    out_strand[current_group] = '+-'
                } else if (current_strand_mask == 1L) {
                    out_strand[current_group] = '+'
                } else if (current_strand_mask == 2L) {
                    out_strand[current_group] = '-'
                } else if (current_first_raw_idx > 0L) {
                    out_strand[current_group] = strand_vals[current_first_raw_idx]
                }
                current_group = current_group + 1L
                current_label_id = label_id[i]
                current_altid_id = altid_id[i]
                current_y_min = y_clean[i]
                current_start_min = start_clean[i]
                current_end_max = end_clean[i]
                current_strand_mask = strand_mask[i]
                current_first_raw_idx = if (strand_not_na[i]) i else 0L
            } else {
                if (y_clean[i] < current_y_min) {
                    current_y_min = y_clean[i]
                }
                if (start_clean[i] < current_start_min) {
                    current_start_min = start_clean[i]
                }
                if (end_clean[i] > current_end_max) {
                    current_end_max = end_clean[i]
                }
                current_strand_mask = mask_or[current_strand_mask + 1L, strand_mask[i] + 1L]
                if ((current_first_raw_idx == 0L) && strand_not_na[i]) {
                    current_first_raw_idx = i
                }
            }
        }
    }
    out_label_id[current_group] = current_label_id
    out_altid_id[current_group] = current_altid_id
    if (is.finite(current_y_min)) {
        out_y[current_group] = current_y_min
    }
    if (is.finite(current_start_min)) {
        out_start[current_group] = current_start_min
    }
    if (is.finite(current_end_max)) {
        out_end[current_group] = current_end_max
    }
    if (current_strand_mask == 3L) {
        out_strand[current_group] = '+-'
    } else if (current_strand_mask == 1L) {
        out_strand[current_group] = '+'
    } else if (current_strand_mask == 2L) {
        out_strand[current_group] = '-'
    } else if (current_first_raw_idx > 0L) {
        out_strand[current_group] = strand_vals[current_first_raw_idx]
    }
    out_label = label_levels[out_label_id]
    out_altid = altid_levels[out_altid_id]
    df_out = data.frame(
        label = out_label,
        motif_altid = out_altid,
        overlap_group = out_group,
        y = out_y,
        start = out_start,
        end = out_end,
        strand = out_strand,
        stringsAsFactors = FALSE
    )
    rownames(df_out) = NULL
    return(df_out)
}

add_complete_overlap_groups = function(df_fimo, altid_col, ncpu=1, verbose=TRUE) {
    all_intervals_overlap = function(iv_a, iv_b) {
        n_a = length(iv_a[['start']])
        n_b = length(iv_b[['start']])
        if (n_a == 0L || n_b == 0L) {
            return(FALSE)
        }
        idx = findInterval(iv_a[['end']], iv_b[['start']])
        if (any(idx == 0L)) {
            return(FALSE)
        }
        return(all(iv_b[['prefix_end_max']][idx] >= iv_a[['start']]))
    }
    altid_raw = df_fimo[[altid_col]]
    altid_vals = as.character(altid_raw)
    valid_altid = !is.na(altid_raw) & nzchar(altid_vals)
    df_fimo[,'complete_overlap_group'] = rep(NA_integer_, nrow(df_fimo))
    if (sum(valid_altid) < 2L) {
        return(df_fimo)
    }
    by_altid_index = split(which(valid_altid), altid_vals[valid_altid], drop = TRUE)
    altids = names(by_altid_index)
    if (length(altids) < 2L) {
        return(df_fimo)
    }
    label_vals = as.character(df_fimo[['label']])
    label_codes = as.integer(factor(label_vals))
    starts_num_all = suppressWarnings(as.numeric(df_fimo[['start']]))
    ends_num_all = suppressWarnings(as.numeric(df_fimo[['end']]))
    genes_by_altid = lapply(by_altid_index, function(idx) {
        codes = unique(label_codes[idx])
        codes = codes[!is.na(codes)]
        if (length(codes) <= 1L) {
            return(codes)
        }
        sort.int(codes, method = 'quick')
    })
    gene_signature = vapply(genes_by_altid, function(g) paste(g, collapse = ","), character(1))
    candidate_groups = split(names(gene_signature), gene_signature, drop = TRUE)
    candidate_groups = candidate_groups[vapply(candidate_groups, length, integer(1)) > 1L]
    if (length(candidate_groups) == 0) {
        return(df_fimo)
    }
    candidate_altids = unique(unlist(candidate_groups, use.names=FALSE))
    intervals_by_altid_gene = setNames(vector('list', length(candidate_altids)), candidate_altids)
    for (altid in candidate_altids) {
        idx_altid = by_altid_index[[altid]]
        label_index = split(idx_altid, label_codes[idx_altid], drop = TRUE)
        genes = genes_by_altid[[altid]]
        n_genes = length(genes)
        if (n_genes == 0L) {
            intervals_by_altid_gene[[altid]] = list()
            next
        }
        gene_intervals = vector('list', n_genes)
        for (gi in seq_len(n_genes)) {
            idx = label_index[[as.character(genes[gi])]]
            if (is.null(idx) || (length(idx) == 0L)) {
                gene_intervals[[gi]] = NULL
                next
            }
            starts = starts_num_all[idx]
            ends = ends_num_all[idx]
            is_valid = is.finite(starts) & is.finite(ends)
            if (!any(is_valid)) {
                gene_intervals[[gi]] = NULL
                next
            }
            starts = starts[is_valid]
            ends = ends[is_valid]
            start_sorted = pmin(starts, ends)
            end_sorted = pmax(starts, ends)
            ord = order(start_sorted, end_sorted, method = 'radix')
            start_sorted = start_sorted[ord]
            end_sorted = end_sorted[ord]
            gene_intervals[[gi]] = list(
                start = start_sorted,
                end = end_sorted,
                prefix_end_max = cummax(end_sorted),
                min_start = start_sorted[1L],
                max_end = end_sorted[length(end_sorted)]
            )
        }
        intervals_by_altid_gene[[altid]] = gene_intervals
    }

    ncpu = suppressWarnings(as.integer(ncpu))
    if (!is.finite(ncpu) || is.na(ncpu) || (ncpu < 1L)) {
        ncpu = 1L
    }
    get_components_in_group = function(grp) {
        n_grp = length(grp)
        if (n_grp < 2L) {
            return(list())
        }
        parent = seq_len(n_grp)
        uf_find_local = function(x) {
            while (parent[x] != x) {
                parent[x] <<- parent[parent[x]]
                x <- parent[x]
            }
            x
        }
        uf_union_local = function(a, b) {
            ra = uf_find_local(a)
            rb = uf_find_local(b)
            if (ra != rb) {
                parent[rb] <<- ra
            }
        }
        for (ii in seq_len(n_grp - 1L)) {
            iv1 = intervals_by_altid_gene[[grp[ii]]]
            n_gene = length(iv1)
            for (jj in seq.int(ii + 1L, n_grp)) {
                if (uf_find_local(ii) == uf_find_local(jj)) {
                    next
                }
                iv2 = intervals_by_altid_gene[[grp[jj]]]
                if ((n_gene == 0L) || (length(iv2) != n_gene)) {
                    next
                }
                is_overlap = TRUE
                for (gi in seq_len(n_gene)) {
                    g1 = iv1[[gi]]
                    g2 = iv2[[gi]]
                    if (is.null(g1) || is.null(g2)) {
                        is_overlap = FALSE
                        break
                    }
                    if ((g1[['min_start']] > g2[['max_end']]) || (g2[['min_start']] > g1[['max_end']])) {
                        is_overlap = FALSE
                        break
                    }
                    if (!all_intervals_overlap(g1, g2) || !all_intervals_overlap(g2, g1)) {
                        is_overlap = FALSE
                        break
                    }
                }
                if (is_overlap) {
                    uf_union_local(ii, jj)
                }
            }
        }
        roots = vapply(seq_len(n_grp), uf_find_local, integer(1))
        comps = split(seq_len(n_grp), roots)
        comps = comps[vapply(comps, length, integer(1)) > 1L]
        if (length(comps) == 0L) {
            return(list())
        }
        lapply(comps, function(idx) grp[idx])
    }
    use_parallel = (ncpu > 1L) && (length(candidate_groups) > 1L) && (.Platform$OS.type != 'windows')
    if (use_parallel) {
        component_lists = parallel::mclapply(candidate_groups, get_components_in_group, mc.cores = ncpu)
    } else {
        component_lists = lapply(candidate_groups, get_components_in_group)
    }
    complete_overlap_groups = unlist(component_lists, recursive = FALSE, use.names = FALSE)

    n_altid = length(altids)
    altid_to_idx = setNames(seq_len(n_altid), altids)
    if (length(complete_overlap_groups)>0) {
        if (verbose) {
            for (i in seq_len(length(complete_overlap_groups))) {
                members = complete_overlap_groups[[i]]
                preview_n = min(6L, length(members))
                preview = paste(members[seq_len(preview_n)], collapse=' ')
                if (length(members) > preview_n) {
                    preview = paste0(preview, ' ... (', length(members), ' motifs)')
                }
                cat('Detected universal overlap (motifs will be merged):', preview, '\n')
            }
        }
        mapped_altid = rep(NA_character_, n_altid)
        mapped_group = rep(NA_integer_, n_altid)
        counter = 1L
        for (i in seq_len(length(complete_overlap_groups))) {
            member_altids = complete_overlap_groups[[i]]
            member_idx = as.integer(altid_to_idx[member_altids])
            my_order = order(nchar(member_altids))
            new_altid = paste(member_altids[my_order], collapse=',')
            if (nchar(new_altid)>20) {
                new_altid = paste0(substr(new_altid, 1, 20), '...')
            }
            mapped_altid[member_idx] = new_altid
            mapped_group[member_idx] = counter
            counter = counter + 1
        }
        idx = as.integer(altid_to_idx[altid_vals])
        is_mapped = valid_altid & !is.na(idx) & !is.na(mapped_group[idx])
        if (any(is_mapped)) {
            df_fimo[is_mapped, 'complete_overlap_group'] = mapped_group[idx[is_mapped]]
            altid_vals[is_mapped] = mapped_altid[idx[is_mapped]]
            df_fimo[[altid_col]] = altid_vals
        }
    }
    return(df_fimo)
}

get_df_N = function(df_tip) {
    df_tip = data.frame(df_tip, stringsAsFactors=FALSE)
    cols = c('x_start', 'x_end', 'y_start', 'y_end')
    x_start_parts = list()
    x_end_parts = list()
    y_start_parts = list()
    y_end_parts = list()
    part_idx = 1L
    for (i in seq_len(nrow(df_tip))) {
        N_slice = df_tip[i,'promoter_N']
        if (is.na(N_slice)) {
            next
        }
        if (N_slice=='') {
            next
        }
        y_val = suppressWarnings(as.numeric(df_tip[i,'y']))
        if (!is.finite(y_val)) {
            next
        }
        slice_vec = strsplit(as.character(N_slice), ',', fixed = TRUE)[[1]]
        if (length(slice_vec) == 0L) {
            next
        }
        for (j in seq_along(slice_vec)) {
            token = trimws(slice_vec[j])
            if (token == '') {
                next
            }
            if (!grepl(':', token, fixed = TRUE)) {
                x = suppressWarnings(as.numeric(token))
                if (!is.finite(x)) {
                    next
                }
                x_start_parts[[part_idx]] = x - 0.5
                x_end_parts[[part_idx]] = x + 0.5
                y_start_parts[[part_idx]] = y_val
                y_end_parts[[part_idx]] = y_val
                part_idx = part_idx + 1L
            } else {
                start_stop = strsplit(token, ':', fixed = TRUE)[[1]]
                if (length(start_stop) < 2L) {
                    next
                }
                x1 = suppressWarnings(as.numeric(start_stop[1]))
                x2 = suppressWarnings(as.numeric(start_stop[2]))
                if (!is.finite(x1) || !is.finite(x2)) {
                    next
                }
                lo = min(x1, x2)
                hi = max(x1, x2)
                x_start_parts[[part_idx]] = lo - 0.5
                x_end_parts[[part_idx]] = hi + 0.5
                y_start_parts[[part_idx]] = y_val
                y_end_parts[[part_idx]] = y_val
                part_idx = part_idx + 1L
            }
        }
    }
    if (part_idx == 1L) {
        df_N = data.frame(matrix(numeric(0), nrow = 0, ncol = length(cols)))
    } else {
        df_N = data.frame(
            x_start = as.numeric(unlist(x_start_parts, use.names = FALSE)),
            x_end = as.numeric(unlist(x_end_parts, use.names = FALSE)),
            y_start = as.numeric(unlist(y_start_parts, use.names = FALSE)),
            y_end = as.numeric(unlist(y_end_parts, use.names = FALSE)),
            stringsAsFactors = FALSE
        )
    }
    colnames(df_N) = cols
    return(df_N)
}

add_fimo_column = function(g, args, qname, xmax=2000, qvalue=0.01, multiple_connection='align_from_tss', ncpu=1) {
    cat(as.character(Sys.time()), 'Adding fimo column.\n')
    if (!'fimo_start'%in%colnames(g[['tree']][['data']])) {
        cat('Scaffold positions are not available. FIMO plot will not be displayed.\n')
        return(g)
    }
    ncpu = suppressWarnings(as.integer(ncpu))
    if (!is.finite(ncpu) || is.na(ncpu) || (ncpu < 1L)) {
        ncpu = 1L
    }
    df_tip = get_df_tip(g[['tree']])
    df_tip = data.frame(df_tip, stringsAsFactors=FALSE)
    to_char_no_na = function(x) {
        v = as.character(x)
        v[is.na(v)] = ''
        v
    }
    split_column_tokens = function(v) {
        out = strsplit(v, ';', fixed = TRUE)
        if (length(out) == 0L) {
            return(out)
        }
        empty_idx = (v == '')
        if (any(empty_idx)) {
            out[empty_idx] = replicate(sum(empty_idx), character(0), simplify = FALSE)
        }
        has_space = grepl(' ', v, fixed = TRUE)
        if (any(has_space)) {
            out[has_space] = lapply(out[has_space], trimws)
            empty_idx2 = has_space & (v == '')
            if (any(empty_idx2)) {
                out[empty_idx2] = replicate(sum(empty_idx2), character(0), simplify = FALSE)
            }
        }
        out
    }
    fimo_start_list = split_column_tokens(to_char_no_na(df_tip[['fimo_start']]))
    fimo_end_list = split_column_tokens(to_char_no_na(df_tip[['fimo_end']]))
    fimo_strand_list = split_column_tokens(to_char_no_na(df_tip[['fimo_strand']]))
    fimo_altid_list = split_column_tokens(to_char_no_na(df_tip[['fimo_alt_id']]))
    fimo_qvalue_list = split_column_tokens(to_char_no_na(df_tip[['fimo_qvalue']]))
    labels_vec = as.character(df_tip[['label']])
    y_vec = suppressWarnings(as.numeric(df_tip[['y']]))
    parse_tip_fimo = function(i) {
        raw_start = fimo_start_list[[i]]
        if (length(raw_start) == 0L) {
            return(list(core = 0L, pass = 0L))
        }
        fimo_starts = suppressWarnings(as.integer(raw_start))
        fimo_ends = suppressWarnings(as.integer(fimo_end_list[[i]]))
        fimo_strands = fimo_strand_list[[i]]
        fimo_altids = fimo_altid_list[[i]]
        fimo_qvalues = suppressWarnings(as.numeric(fimo_qvalue_list[[i]]))
        n = min(length(fimo_starts), length(fimo_ends), length(fimo_strands), length(fimo_altids), length(fimo_qvalues))
        if (n == 0L) {
            return(list(core = 0L, pass = 0L))
        }
        if (length(fimo_starts) != n) {
            fimo_starts = fimo_starts[seq_len(n)]
        }
        if (length(fimo_ends) != n) {
            fimo_ends = fimo_ends[seq_len(n)]
        }
        if (length(fimo_strands) != n) {
            fimo_strands = fimo_strands[seq_len(n)]
        }
        if (length(fimo_altids) != n) {
            fimo_altids = fimo_altids[seq_len(n)]
        }
        if (length(fimo_qvalues) != n) {
            fimo_qvalues = fimo_qvalues[seq_len(n)]
        }
        keep_core = is.finite(fimo_starts) & is.finite(fimo_ends) & is.finite(fimo_qvalues) &
            !is.na(fimo_altids) & (fimo_altids != '')
        n_core = sum(keep_core)
        if (n_core == 0L) {
            return(list(core = 0L, pass = 0L))
        }
        keep_q = keep_core & (fimo_qvalues <= qvalue)
        n_pass = sum(keep_q)
        if (n_pass == 0L) {
            return(list(core = n_core, pass = 0L))
        }
        return(list(
            core = n_core,
            pass = n_pass,
            row_idx = rep.int(i, n_pass),
            start = fimo_starts[keep_q],
            end = fimo_ends[keep_q],
            strand = fimo_strands[keep_q],
            qvalue = fimo_qvalues[keep_q],
            motif_altid = fimo_altids[keep_q]
        ))
    }
    index_vec = seq_len(nrow(df_tip))
    use_parallel = (ncpu > 1L) && (length(index_vec) > 1L) && (.Platform$OS.type != 'windows')
    if (use_parallel) {
        parsed_parts = parallel::mclapply(index_vec, parse_tip_fimo, mc.cores = ncpu)
    } else {
        if ((ncpu > 1L) && (.Platform$OS.type == 'windows')) {
            cat('ncpu > 1 was requested, but mclapply is unsupported on Windows. Falling back to ncpu=1.\n')
        }
        parsed_parts = vector('list', length(index_vec))
        for (i in index_vec) {
            parsed_parts[[i]] = parse_tip_fimo(i)
        }
    }
    core_counts = vapply(parsed_parts, function(x) x[['core']], integer(1))
    pass_counts = vapply(parsed_parts, function(x) x[['pass']], integer(1))
    num_core_records = sum(core_counts)
    num_pass_records = sum(pass_counts)
    if (num_core_records == 0L) {
        cat('No FIMO records are available. FIMO plot will not be displayed.\n')
        return(g)
    }
    if (num_pass_records == 0L) {
        cat('No FIMO records passed the qvalue cutoff. FIMO plot will not be displayed.\n')
        return(g)
    }
    row_idx_out = integer(num_pass_records)
    start_out = numeric(num_pass_records)
    end_out = numeric(num_pass_records)
    strand_out = character(num_pass_records)
    qvalue_out = numeric(num_pass_records)
    motif_out = character(num_pass_records)
    write_pos = 1L
    for (i in seq_along(parsed_parts)) {
        n_pass = pass_counts[i]
        if (n_pass == 0L) {
            next
        }
        idx = write_pos:(write_pos + n_pass - 1L)
        part = parsed_parts[[i]]
        row_idx_out[idx] = part[['row_idx']]
        start_out[idx] = part[['start']]
        end_out[idx] = part[['end']]
        strand_out[idx] = part[['strand']]
        qvalue_out[idx] = part[['qvalue']]
        motif_out[idx] = part[['motif_altid']]
        write_pos = write_pos + n_pass
    }
    df_fimo = data.frame(
        y = y_vec[row_idx_out],
        label = labels_vec[row_idx_out],
        start = start_out,
        end = end_out,
        strand = strand_out,
        qvalue = qvalue_out,
        motif_altid = motif_out,
        stringsAsFactors = FALSE
    )
    df_fimo = merge_overlap_group(df_fimo, merge_level='TFfamily', ncpu=ncpu, verbose=FALSE)
    strand_chr = as.character(df_fimo[['strand']])
    delta_up = rep(0, nrow(df_fimo))
    delta_down = rep(0, nrow(df_fimo))
    hjust_vals = rep(0.5, nrow(df_fimo))
    is_plus = (strand_chr == '+')
    is_minus = (strand_chr == '-')
    is_both = (strand_chr == '+-') | (strand_chr == '-+')
    delta_up[is_plus] = 0.3
    delta_down[is_minus] = 0.3
    delta_up[is_both] = 0.15
    delta_down[is_both] = 0.15
    hjust_vals[is_plus] = -0.05
    hjust_vals[is_minus] = 1.05
    df_fimo[,'ymin'] = df_fimo[['y']] - delta_down
    df_fimo[,'ymax'] = df_fimo[['y']] + delta_up
    df_fimo[,'hjust'] = hjust_vals
    df_fimo[,'start'] = xmax - df_fimo[['start']]
    df_fimo[,'end'] = xmax - df_fimo[['end']]
    df_fimo[,'xmid'] = (df_fimo[,'start']+df_fimo[,'end'])/2
    df_fimo = add_motif_counts(df_fimo)
    df_fimo = df_fimo[order(df_fimo[['y']], df_fimo[['start']]),]
    df_polygon = get_df_polygon(df_fimo, min_count=2, max_count=Inf, multiple_connection=multiple_connection, ncpu=ncpu, verbose=FALSE)
    legend_texts = unique(df_fimo[['motif_altid']])
    unique_altid = unique(df_fimo[(df_fimo[['count']]==1),'motif_altid'])
    my_palette = rep('gray50', length(legend_texts))
    num_colored = length(legend_texts) - length(unique_altid)
    if (num_colored > 0L) {
        my_palette[!legend_texts%in%unique_altid] = gg_color_hue(num_colored)
    }
    df_fimo[,'motif_altid'] = factor(df_fimo[,'motif_altid'], levels=legend_texts)
    df_N = get_df_N(df_tip)
    df_N[,'x_start'] = xmax - df_N[['x_start']]
    df_N[,'x_end'] = xmax - df_N[['x_end']]
    df_backbone = df_tip[(df_tip[['promoter_available']]=='Y'),c('label','y')]
    for (power in 0:10) {
        if ((10**(power+1))>xmax) {
            break
        }
    }
    xbreaks = (0:(xmax%/%(10**power))) * (10**power)
    xlabels = xbreaks - xmax
    polygon_layer = NULL
    alpha_scale = NULL
    if (nrow(df_polygon) > 0L) {
        alpha_vals = df_polygon[['alpha_value']]
        alpha_vals = alpha_vals[is.finite(alpha_vals)]
        alpha_range = c(0, 0.5)
        if (length(alpha_vals) > 0L) {
            alpha_range = range(alpha_vals)
            if (alpha_range[1] == alpha_range[2]) {
                alpha_range = c(max(0, alpha_range[1] - 1e-6), alpha_range[2] + 1e-6)
            }
        }
        polygon_layer = geom_polygon(data=df_polygon, mapping=aes(x=x, y=y, fill=motif_altid, group=polygon_id, alpha=alpha_value))
        alpha_scale = scale_alpha_continuous(range=alpha_range)
    }
    g[[qname]] = ggplot(data=df_fimo) +
        geom_blank(data=df_tip, aes(y=label)) + 
        geom_segment(data=df_backbone, mapping=aes(y=y, yend=y), x=1, xend=xmax, linewidth=0.1, color='black') +
        geom_segment(data=df_N, mapping=aes(x=x_start, xend=x_end, y=y_start, yend=y_end), linewidth=0.1, color='gray70') +
        geom_rect(mapping=aes(xmin=start, xmax=end, ymin=ymin, ymax=ymax, fill=motif_altid), alpha=0.5) +
        geom_text(mapping=aes(x=xmid, y=y, label=motif_altid, color=motif_altid, hjust=hjust), size=0.2, angle=90, alpha=0.1) +
        #geom_text_repel(mapping=aes(x=xmid, y=y, label=motif_altid), size=0.15, angle=90, max.overlaps=1000, segment.size=0.1) +
        polygon_layer +
        scale_x_continuous(breaks=xbreaks, labels=xlabels, limits=c(-0.5, xmax+0.5)) +
        scale_fill_manual(values=my_palette) +
        scale_color_manual(values=my_palette) +
        alpha_scale +
        #aplot::ylim2(gg=g[['tree']]) +
        xlab('Motif position') +
        theme_minimal(base_size=args[['font_size']]) +
        guides(
            fill=guide_legend(title="", nrow=6, byrow=FALSE),
            color='none',
            alpha='none'
        ) +
        coord_cartesian(clip = "off") +
        theme(
            axis.title.y=element_blank(),
            axis.title.x=element_text(size=args[['font_size']]),
            axis.line=element_blank(),
            axis.ticks=element_blank(),
            axis.text.y=element_blank(),
            axis.text.x=element_text(size=args[['font_size']], color='black'),
            panel.grid.major.y=element_blank(),
            panel.grid.minor.y=element_blank(),
            panel.grid.minor.x=element_blank(),
            legend.position="bottom",
            legend.justification = "right",
            legend.title=element_text(size=args[['font_size']]),
            legend.text=element_text(size=args[['font_size']]),
            legend.box.just='center',
            legend.key.size=unit(0.4, 'lines'),
            rect=element_rect(fill="transparent"),
            plot.margin=unit(args[['margins']], "cm")
        )
    return(g)
}

extract_alignment_sites <- function(tidy_aln, selected_amino_acid_sites) {
  # Filter for only the selected positions
  tidy_aln <- tidy_aln[tidy_aln[["position"]] %in% selected_amino_acid_sites, ]
  
  # Keep the original position
  tidy_aln[["original_position"]] <- tidy_aln[["position"]]
  
  # Re-index positions using match() rather than a loop
  unique_positions <- unique(tidy_aln[["position"]])
  tidy_aln[["position"]] <- match(tidy_aln[["position"]], unique_positions)
  
  rownames(tidy_aln) <- NULL
  tidy_aln
}

add_amino_acid_site_column <- function(g, args, tidy_aln, selected_amino_acid_sites) {
  cat(as.character(Sys.time()), "Adding amino acid site column.\n")
  df_tip <- get_df_tip(g[["tree"]])[,c('label','branch_id')]
  tidy_site <- extract_alignment_sites(tidy_aln, selected_amino_acid_sites)
  tidy_site <- merge(tidy_site, df_tip, by.x = "name", by.y = "label", all.x = TRUE)
  tidy_site[["name"]] <- factor(tidy_site[["name"]], levels = df_tip[["label"]])
  qname <- "amino_acid_site"
    g[[qname]] <- ggplot() + 
    geom_blank(data = df_tip, aes(y = label)) + 
    ggmsa::geom_msa(data = tidy_site, color = "LETTER") +
    scale_x_continuous(
        breaks = 1:length(selected_amino_acid_sites),
        labels = selected_amino_acid_sites
    ) +
    xlab("Amino acid position (aa)") +
    theme(
        axis.title.y        = element_blank(),
        axis.title.x        = element_text(size = args[["font_size"]]),
        axis.text.x         = element_text(angle = 90, hjust = 1, vjust = 0.5, color = 'black', size = args[["font_size"]]),
        axis.ticks.x        = element_blank(),
        axis.line.x         = element_blank(),
        axis.text.y         = element_blank(),
        axis.ticks.y        = element_blank(),
        axis.line.y         = element_blank(), 
        panel.grid.major.y  = element_blank(),
        rect                = element_blank(),
        plot.margin         = unit(args[["margins"]] / 4, "cm")
    )
  
  g
}
