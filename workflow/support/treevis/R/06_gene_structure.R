# Coordinates are 1-based inclusive GFF blocks, serialized in transcript order.
treevis_gene_structure_data = function(tips, mode='compressed') {
    if (!(mode %in% c('compressed', 'linear'))) stop('Unknown gene structure mode: ', mode)
    boxes = list()
    lines = list()
    trans_splices = list()
    for (i in seq_len(nrow(tips))) {
        value = as.character(tips[['feature_blocks']][i])
        if (is.na(value) || !nzchar(value)) next
        if (is.na(tips[['feature_type']][i]) || tips[['feature_type']][i] != 'CDS') next
        fields = strsplit(value, ';', fixed=TRUE)[[1]]
        if (any(!grepl('^[0-9]+-[0-9]+$', fields))) stop('Invalid GFF blocks for ', tips[['label']][i])
        blocks = do.call(rbind, lapply(strsplit(fields, '-', fixed=TRUE), as.numeric))
        if (any(!is.finite(blocks)) || any(blocks[,1] < 1 | blocks[,2] < blocks[,1])) {
            stop('Invalid GFF coordinates for ', tips[['label']][i])
        }
        trans = 'splice_mode' %in% names(tips) &&
            !is.na(tips[['splice_mode']][i]) && tips[['splice_mode']][i] == 'trans-splicing'
        if (trans) {
            required = c('feature_block_sequences','feature_block_strands','transcript_junction_positions')
            if (!all(required %in% names(tips)) || any(is.na(tips[i,required]))) {
                stop('Missing trans-splicing coordinates for ', tips[['label']][i])
            }
            sequences = strsplit(as.character(tips[['feature_block_sequences']][i]), ';', fixed=TRUE)[[1]]
            strands = strsplit(as.character(tips[['feature_block_strands']][i]), ';', fixed=TRUE)[[1]]
            if (length(sequences) != nrow(blocks) || any(!nzchar(sequences)) ||
                length(strands) != nrow(blocks) || any(!strands %in% c('+','-'))) {
                stop('Invalid trans-splicing coordinate systems')
            }
            lens = blocks[,2] - blocks[,1] + 1
            boundary_text = as.character(tips[['transcript_junction_positions']][i])
            boundaries = if (nzchar(boundary_text)) suppressWarnings(as.numeric(strsplit(boundary_text, ';', fixed=TRUE)[[1]])) else numeric()
            expected = utils::head(cumsum(lens), -1)
            if (!identical(boundaries, expected)) stop('Invalid trans-splicing junction positions')
            if ('num_intron' %in% names(tips) && !is.na(tips[['num_intron']][i])) {
                stop('Trans-splicing must not be reported as a cis-intron count')
            }
            if ('intron_feature_size' %in% names(tips) && !is.na(tips[['intron_feature_size']][i]) &&
                tips[['intron_feature_size']][i] != sum(lens)) stop('CDS length disagrees with GFF blocks')
            if ('utr_blocks' %in% names(tips) && !is.na(tips[['utr_blocks']][i]) && nzchar(tips[['utr_blocks']][i])) {
                stop('Trans-spliced UTR order is not represented')
            }
            starts = c(0, expected)
            color = if ('tiplab_color' %in% names(tips)) as.character(tips[['tiplab_color']][i]) else 'black'
            boxes[[length(boxes)+1]] = data.frame(label=as.character(tips[['label']][i]),
                y=tips[['y']][i], start=starts, end=starts+lens, colour=color,
                feature='CDS', fill=color, half_height=0.27)
            if (length(boundaries)) trans_splices[[length(trans_splices)+1]] = data.frame(
                label=as.character(tips[['label']][i]), y=tips[['y']][i], position=boundaries)
            next
        }
        strand = as.character(tips[['strand']][i])
        if (is.na(strand) || !(strand %in% c('+', '-'))) stop('Unknown GFF strand for ', tips[['label']][i])
        lens = blocks[,2] - blocks[,1] + 1
        gaps = if (nrow(blocks) > 1) {
            if (strand == '+') blocks[-1,1] - blocks[-nrow(blocks),2] - 1 else
                blocks[-nrow(blocks),1] - blocks[-1,2] - 1
        } else numeric()
        if (any(gaps < 0)) stop('Overlapping or unordered GFF blocks for ', tips[['label']][i])
        if ('num_intron' %in% names(tips) && !is.na(tips[['num_intron']][i]) &&
            tips[['num_intron']][i] != sum(gaps > 0)) stop('Intron count disagrees with GFF blocks')
        if ('intron_feature_size' %in% names(tips) && !is.na(tips[['intron_feature_size']][i]) &&
            tips[['intron_feature_size']][i] != sum(lens)) stop('CDS length disagrees with GFF blocks')
        # UTRs are explicit annotations from the selected CDS transcript only.
        feature = rep('CDS', nrow(blocks))
        if ('utr_blocks' %in% names(tips)) {
            utr = as.character(tips[['utr_blocks']][i])
            if (!is.na(utr) && nzchar(utr)) {
                parts = strsplit(utr, ';', fixed=TRUE)[[1]]
                if (any(!grepl('^[0-9]+-[0-9]+$', parts))) stop('Invalid UTR blocks')
                utr = do.call(rbind, lapply(strsplit(parts, '-', fixed=TRUE), as.numeric))
                if (any(!is.finite(utr)) || any(utr[,1] < 1 | utr[,2] < utr[,1])) stop('Invalid UTR coordinates')
                feature = c(feature, rep('UTR', nrow(utr)))
                blocks = rbind(blocks, utr)
                order = order(blocks[,1], decreasing=(strand == '-'))
                blocks = blocks[order,,drop=FALSE]
                feature = feature[order]
            }
        }
        lens = blocks[,2] - blocks[,1] + 1
        gaps = if (nrow(blocks) > 1) {
            if (strand == '+') blocks[-1,1] - blocks[-nrow(blocks),2] - 1 else
                blocks[-nrow(blocks),1] - blocks[-1,2] - 1
        } else numeric()
        if (any(gaps < 0)) stop('Overlapping CDS/UTR blocks')
        displayed_gaps = if (mode == 'compressed') 100 * log1p(gaps / 100) else gaps
        starts = c(0, cumsum(lens[-length(lens)] + displayed_gaps))
        color = if ('tiplab_color' %in% names(tips)) as.character(tips[['tiplab_color']][i]) else 'black'
        boxes[[length(boxes)+1]] = data.frame(
            label=as.character(tips[['label']][i]), y=tips[['y']][i],
            start=starts, end=starts+lens, colour=color, feature=feature,
            fill=ifelse(feature == 'CDS', color,
                grDevices::rgb(t(grDevices::col2rgb(color) / 255 * 0.3 + 0.7))),
            half_height=ifelse(feature == 'CDS', 0.27, 0.18))
        if (any(gaps > 0)) lines[[length(lines)+1]] = data.frame(
            label=as.character(tips[['label']][i]), y=tips[['y']][i],
            start=(starts+lens)[seq_along(gaps)], end=starts[-1],
            length_bp=gaps, colour=color,
            cds_offset=cumsum(ifelse(feature == 'CDS', lens, 0))[seq_along(gaps)])[gaps > 0,,drop=FALSE]
    }
    list(boxes=if (length(boxes)) do.call(rbind, boxes) else data.frame(),
         introns=if (length(lines)) do.call(rbind, lines) else data.frame(),
         trans_splices=if (length(trans_splices)) do.call(rbind, trans_splices) else data.frame())
}

add_gene_structure_column = function(g, args, mode='compressed', width_mm=23, path_alignment=NULL) {
    if (length(width_mm) != 1 || !is.finite(width_mm) || width_mm <= 0) stop('Invalid gene structure width')
    tips = get_df_tip(g[['tree']])
    if (!all(c('feature_blocks','feature_type','strand') %in% names(tips))) return(g)
    data = treevis_gene_structure_data(tips, mode)
    if (!nrow(data$boxes)) return(g)
    correspondence = NULL
    labels = data.frame()
    if (!is.null(path_alignment) && !is.na(path_alignment) && file.exists(path_alignment) &&
        file.info(path_alignment)$size > 0 &&
        all(c('num_intron','intron_positions','intron_feature_size','cds_first_phase') %in% names(tips))) {
        correspondence = treevis_intron_site_data(tips, treevis_similarity_fasta(path_alignment))
        if (!is.null(correspondence)) {
            mismatches = correspondence$diagnostics$node_name[correspondence$diagnostics$reason == 'CDS_length_mismatch']
            if (length(mismatches)) stop('GFF/CDS length mismatch; repair the input annotation before plotting: ',
                                        paste(mismatches, collapse=', '))
            labels = treevis_structure_intron_labels(data$introns, correspondence$events)
            if (nrow(correspondence$diagnostics)) cat('Intron correspondence excluded tips:',
                paste(correspondence$diagnostics$node_name, correspondence$diagnostics$reason, sep=':',collapse=', '), '\n')
        }
    }
    axis_divisor = if (max(data$boxes$end) >= 1000) 1000 else 1
    p = ggplot(data$boxes, aes(y=.data[['y']])) +
        geom_rect(aes(xmin=.data[['start']], xmax=.data[['end']],
                      ymin=.data[['y']]- .data[['half_height']], ymax=.data[['y']]+ .data[['half_height']],
                      fill=.data[['fill']], colour=.data[['colour']]), linewidth=0.12) +
        scale_fill_identity() + scale_colour_identity() +
        scale_x_continuous(expand=expansion(mult=c(0.01,0.01)),
            breaks=if (mode == 'linear') scales::breaks_pretty(n=3) else NULL,
            labels=scales::label_number(scale=1/axis_divisor)) +
        scale_y_continuous(limits=ggplot_build(g[['tree']])$layout$panel_params[[1]]$y.range,
            expand=c(0,0)) +
        labs(x=if (mode == 'linear') paste0('Exon / intron (', if (axis_divisor == 1000) 'kb' else 'bp', ')') else 'Exon / intron\n(compressed)',
             caption=if (mode == 'compressed') 'CDS: solid; UTR: pale\nIntrons: 100 ln(1 + bp/100)' else 'CDS: solid; UTR: pale') +
        theme_void(base_size=args[['font_size']]) +
        theme(axis.title.x=element_text(size=args[['font_size']]),
              axis.text.x=element_text(size=args[['font_size']]),
              plot.caption=element_text(size=args[['font_size']], hjust=0.5),
              plot.margin=unit(args[['margins']], 'cm'))
    if (nrow(data$introns)) p = p +
        geom_segment(data=data$introns,
            aes(x=.data[['start']], xend=.data[['end']], yend=.data[['y']],
                colour=.data[['colour']]), linewidth=0.25)
    if (nrow(data$trans_splices)) {
        p = p + geom_segment(data=data$trans_splices,
            aes(x=.data[['position']], xend=.data[['position']],
                y=.data[['y']]-0.33, yend=.data[['y']]+0.33),
            inherit.aes=FALSE, linetype='dashed', linewidth=0.3) +
            labs(caption=paste0(p$labels$caption,
                '\nDashed vertical: trans-spliced joins\n(no genomic distance inferred)'))
    }
    if (nrow(labels)) {
        units_per_mm = max(data$boxes$end) / width_mm
        polygons = treevis_intron_connection_polygons(labels, half_width=0.3 * units_per_mm)
        if (nrow(polygons)) {
            p = p + geom_polygon(data=polygons,
                aes(x=.data[['x']],y=.data[['y']],group=.data[['connection']],alpha=.data[['opacity']]),
                inherit.aes=FALSE, fill='grey40',colour=NA) + scale_alpha_identity()
            # The connection bands sit behind CDS/UTR boxes and intron lines.
            p$layers = p$layers[c(length(p$layers),seq_len(length(p$layers)-1))]
        }
        attr(p, 'treevis_intron_connections') = polygons
        p = p + labs(caption=paste0(p$labels$caption,
            '\nDark: exact position\n+ phase\nLight: nearby position\n(<=3 nt)'))
    }
    attr(p, 'treevis_intron_sites') = correspondence
    attr(p, 'treevis_intron_labels') = labels
    attr(p, 'treevis_width_mm') = width_mm
    g[['gene_structure']] = p
    g
}

# Join by CDS offset, so UTR-only introns do not shift intron numbering.
treevis_structure_intron_labels = function(introns, events) {
    if (!nrow(introns) || !nrow(events)) return(data.frame())
    key = function(label, offset) paste(label, offset, sep='::')
    index = match(key(introns$label,introns$cds_offset), key(events$node_name,events$cds_offset))
    matched = !is.na(index)
    labels = introns[matched,,drop=FALSE]
    event = events[index[matched],,drop=FALSE]
    labels$site_id = event$site_id
    labels$status = event$status
    labels$alignment_left = event$alignment_left
    labels$phase = event$phase
    number = as.character(suppressWarnings(as.integer(sub('^I','',event$site_id))))
    labels$display_id = ifelse(is.na(number), '?',
        ifelse(event$status == 'mapped',number,paste0(number,'?')))
    labels$x = (labels$start+labels$end)/2
    labels
}

# Follow the next nearby occurrence in tip order, requiring reciprocal unique
# nearest boundaries. An uncertain occurrence interrupts a strong band.
treevis_intron_connection_polygons = function(labels, half_width, tolerance_nt=3) {
    if (!nrow(labels)) return(data.frame())
    if (length(tolerance_nt)!=1 || !is.finite(tolerance_nt) || tolerance_nt<0 ||
        tolerance_nt!=floor(tolerance_nt)) stop('Invalid intron position tolerance')
    labels = labels[labels$status %in% c('mapped','position_only'),,drop=FALSE]
    labels$connection_index = seq_len(nrow(labels))
    polygons = list()
    seen = character()
    add_band = function(a,b,kind) {
        pair = paste(sort(c(a$connection_index,b$connection_index)),collapse=':')
        if (pair %in% seen || a$y == b$y) return(invisible(NULL))
        seen <<- c(seen,pair)
        wa = min(half_width,(a$end-a$start)/2)
        wb = min(half_width,(b$end-b$start)/2)
        polygons[[length(polygons)+1]] <<- data.frame(
            connection=pair,site_id=if(kind=='exact_phase') a$site_id else NA_character_,
            kind=kind,opacity=if(kind=='exact_phase') 0.16 else 0.05,
            from_index=a$connection_index,to_index=b$connection_index,
            position_delta=abs(a$alignment_left-b$alignment_left),
            x=c(a$x-wa,a$x+wa,b$x+wb,b$x-wb),y=c(a$y,a$y,b$y,b$y))
    }
    rows = sort(unique(labels$y))
    if (length(rows)>1) for (r in seq_len(length(rows)-1)) {
        a = labels[labels$y == rows[r],,drop=FALSE]
        for (i in seq_len(nrow(a))) for (next_row in rows[(r+1):length(rows)]) {
            b = labels[labels$y == next_row,,drop=FALSE]
            distance = abs(outer(a$alignment_left,b$alignment_left,'-'))
            if (min(distance[i,])>tolerance_nt) next
            candidates = which(distance[i,] == min(distance[i,]))
            # An ambiguous nearest row is not bypassed to manufacture a link.
            if (length(candidates)!=1) break
            j = candidates[1]
            reverse = which(distance[,j] == min(distance[,j]))
            if (length(reverse)!=1 || reverse[1]!=i) break
            exact = a$status[i]=='mapped' && b$status[j]=='mapped' &&
                !is.na(a$site_id[i]) && !is.na(b$site_id[j]) && a$site_id[i]==b$site_id[j]
            add_band(a[i,],b[j,],if(exact) 'exact_phase' else 'position_candidate')
            break
        }
    }
    if (length(polygons)) do.call(rbind,polygons) else data.frame()
}
