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

reshape_long_base = function(df, id_col, value_cols, key_col, value_col) {
    if (!(id_col %in% colnames(df))) {
        stop(paste0("reshape_long_base: id column not found: ", id_col))
    }
    if (length(value_cols) == 0 || nrow(df) == 0) {
        out = data.frame(id = character(0), key = character(0), value = character(0), stringsAsFactors = FALSE, check.names = FALSE)
        colnames(out) = c(id_col, key_col, value_col)
        return(out)
    }
    out = data.frame(
        id = rep(df[[id_col]], times = length(value_cols)),
        key = rep(value_cols, each = nrow(df)),
        value = as.vector(as.matrix(df[, value_cols, drop = FALSE])),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    colnames(out) = c(id_col, key_col, value_col)
    return(out)
}

tidy_df_tip = function(df_tip, df_trait) {
    trait_cols = colnames(df_trait)
    if (length(trait_cols) == 0) {
        return(data.frame())
    }
    trait_cols = sort(trait_cols)
    # Keep only trait columns to avoid reshaping unrelated columns in df_tip.
    keep_cols = c('label', intersect(trait_cols, colnames(df_tip)))
    value_cols = setdiff(keep_cols, 'label')
    df_tmp = reshape_long_base(
        df = df_tip[, keep_cols, drop = FALSE],
        id_col = 'label',
        value_cols = value_cols,
        key_col = 'group',
        value_col = 'value'
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
        legend.position = "inside",
        legend.position.inside = c(0, 1),
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
        } else if (grepl('^categorical,', gname)) {
            rel_widths[gname] = 0.55
        } else if (grepl('^text,', gname)) {
            rel_widths[gname] = 0.7
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
        } else if (grepl('^meme$', gname)) {
            rel_widths[gname] = 0.9
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
    treevis_species_parser = Sys.getenv('TREEVIS_SPECIES_PARSER', unset='taxonomic')
    if (!nzchar(treevis_species_parser)) {
        treevis_species_parser = 'taxonomic'
    }
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
            spp = sort(unique(get_species_name(spp, species_parser=treevis_species_parser)))
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
