treevis_tip_pairs <- function(n) {
    if (n < 2L) stop('n < m')
    # Keep combn order in two vectors, without assembling a pair matrix.
    first <- seq_len(n - 1L)
    remaining <- n - first
    list(first=rep.int(first, remaining), second=sequence(remaining, from=first + 1L))
}

treevis_tip_groups <- function(labels, node_names, group_ids) {
    if (anyNA(labels) || anyNA(node_names)) {
        return(lapply(labels, function(label) unique(group_ids[node_names == label])))
    }
    by_tip <- split(group_ids, node_names)
    lapply(by_tip[match(labels, names(by_tip))], unique)
}

treevis_jaccard <- function(groups, observed) {
    sizes <- as.double(lengths(groups))
    n <- length(groups)
    first_tips <- seq_len(n - 1L)
    offsets <- c(0, cumsum(as.double(n - first_tips)))
    shared <- integer(offsets[n])
    group_ids <- unlist(groups, use.names=FALSE)
    members <- split(rep.int(seq_len(n), sizes), match(group_ids, unique(group_ids)))
    # Within each group, add one to the triangular index for every pair of
    # member tips. The grouped input has at most one entry per tip and group.
    for (tips in members[lengths(members) > 1L]) {
        count <- length(tips)
        for (index in seq_len(count - 1L)) {
            first <- tips[index]
            second <- tips[seq.int(index + 1L, count)]
            positions <- offsets[first] + (second - first)
            shared[positions] <- shared[positions] + 1L
        }
    }
    result <- rep(NA_real_, length(shared))
    for (first in first_tips[observed[first_tips]]) {
        second <- seq.int(first + 1L, n)
        second <- second[observed[second]]
        positions <- offsets[first] + (second - first)
        common <- shared[positions]
        total <- sizes[first] + sizes[second] - common
        valid <- total > 0
        result[positions[valid]] <- common[valid] / total[valid]
    }
    result
}

# Pairwise neighborhood group similarity in the plotted tree's tip order.
add_synteny_similarity_column = function(g, args, path_synteny, search_window=20, width_mm=15) {
    if (!file.exists(path_synteny) || file.info(path_synteny)$size == 0) return(g)
    if (length(search_window) != 1 || !is.finite(search_window) ||
        search_window < 1 || search_window != as.integer(search_window)) {
        stop('Synteny similarity search window must be a positive integer.')
    }
    if (length(width_mm) != 1 || !is.finite(width_mm) || width_mm <= 0) {
        stop('Synteny similarity width must be positive (mm).')
    }
    tips = get_df_tip(g[['tree']])
    tips = tips[order(tips[['y']]), , drop=FALSE]
    if (nrow(tips) < 2) return(g)
    syn = read.table(path_synteny, header=TRUE, sep='\t', quote='', comment.char='',
                     colClasses='character', na.strings=character(0), stringsAsFactors=FALSE)
    if (nrow(syn) == 0) return(g)
    if (!all(c('node_name', 'offset', 'group_id') %in% names(syn))) {
        stop('Synteny similarity input requires node_name, offset, and group_id.')
    }
    offsets = suppressWarnings(as.numeric(syn[['offset']]))
    valid = is.finite(offsets) & offsets != 0 & abs(offsets) <= search_window &
        !is.na(syn[['group_id']]) & nzchar(syn[['group_id']]) &
        syn[['node_name']] %in% tips[['label']]
    syn = syn[valid, , drop=FALSE]
    groups = treevis_tip_groups(as.character(tips[['label']]), syn[['node_name']], syn[['group_id']])
    observed = lengths(groups) > 0L
    if (sum(observed) < 2) return(g)
    similarities = treevis_jaccard(groups, observed)
    treevis_pairwise_column(g, args, tips, similarities, 'synteny_similarity',
        paste0('Syntenic\nsimilarity\n(',search_window,' genes\nper side)'), width_mm)
}

treevis_pairwise_column = function(g, args, tips, similarities, gname, axis_label, width_mm=15) {
    if (length(width_mm) != 1 || !is.finite(width_mm) || width_mm <= 0) {
        stop('Pairwise similarity width must be positive (mm).')
    }
    pairs = treevis_tip_pairs(nrow(tips))
    # Use a common unit for tip spacing; each diamond is one pair, not a tip row.
    spacing = min(diff(sort(unique(tips[['y']]))))
    y1 = tips[['y']][pairs$first]; y2 = tips[['y']][pairs$second]
    x = (y2-y1)/2; y = (y1+y2)/2
    polygons = data.frame(
        pair=rep(seq_along(x), each=4),
        x=as.vector(rbind(x-spacing/2, x, x+spacing/2, x)),
        y=as.vector(rbind(y, y-spacing/2, y, y+spacing/2)),
        similarity=rep(similarities, each=4))
    p = ggplot(polygons, aes(x=.data[['x']], y=.data[['y']])) +
        geom_polygon(aes(group=.data[['pair']], fill=.data[['similarity']]), colour=NA) +
        scale_fill_viridis_c(limits=c(0,1), na.value='grey80', name=NULL,
                            breaks=c(0,0.5,1), guide=guide_colorbar(ticks=TRUE, barwidth=unit(10,'mm'), barheight=unit(2,'mm'))) +
        scale_x_continuous(limits=c(0, (max(tips[['y']])-min(tips[['y']])+spacing)/2), expand=c(0,0), breaks=NULL) +
        scale_y_continuous(limits=ggplot_build(g[['tree']])$layout$panel_params[[1]]$y.range, expand=c(0,0)) +
        xlab(axis_label) +
        theme_void(base_size=args[['font_size']]) +
        theme(axis.title.x=element_text(size=args[['font_size']], margin=margin(t=3)),
              legend.position='bottom', legend.margin=margin(0,0,0,0),
              legend.text=element_text(size=args[['font_size']]),
              legend.ticks=element_line(colour='black', linewidth=0.3, lineend='butt'),
              legend.ticks.length=unit(c(-0.6,0),'mm'),
              plot.margin=margin(0,0,0,0))
    attr(p, 'treevis_width_mm') = width_mm
    # Insert labels only when the first available triangle is actually added.
    # Keep the suffix literal (including transcript suffixes); the full IDs
    # remain in the main tip-label column, in exactly the same row order.
    if (!('pairwise_tip_suffix' %in% names(g))) {
        ids = as.character(tips[['label']])
        labels = data.frame(y=tips[['y']], suffix=paste0('...', substring(ids, pmax(1, nchar(ids)-2))))
        label_plot = ggplot(labels, aes(y=.data[['y']])) +
            geom_text(aes(label=.data[['suffix']]), x=0.95, hjust=1,
                      size=args[['font_size']] / (72.27 / 25.4),
                      colour=tips[['tiplab_color']]) +
            scale_x_continuous(limits=c(0,1), expand=c(0,0)) +
            scale_y_continuous(limits=ggplot_build(g[['tree']])$layout$panel_params[[1]]$y.range,
                               expand=c(0,0)) +
            theme_void() + theme(plot.margin=margin(0,0,0,0))
        attr(label_plot, 'treevis_width_mm') = 8
        g[['pairwise_tip_suffix']] = label_plot
    }
    g[[gname]] = p
    return(g)
}

# Read IDs exactly as FASTA tools do (first whitespace-delimited token).
treevis_similarity_fasta = function(path) {
    con = if (grepl('\\.gz$',path)) gzfile(path,'rt') else file(path,'rt')
    on.exit(close(con))
    lines = trimws(readLines(con,warn=FALSE))
    lines = lines[nzchar(lines)]
    starts = which(startsWith(lines,'>'))
    if (!length(starts)) {
        if (length(lines)) stop('Nonempty FASTA input has no record headers.')
        return(setNames(character(0),character(0)))
    }
    if (starts[1] != 1) stop('FASTA sequence appeared before its header.')
    ids = sub('[[:space:]].*$', '', substring(lines[starts],2))
    if (any(!nzchar(ids)) || anyDuplicated(ids)) stop('FASTA IDs must be nonempty and unique.')
    ends = c(starts[-1]-1,length(lines))
    seqs = vapply(seq_along(starts),function(i) {
        if (ends[i] == starts[i]) '' else toupper(paste0(lines[(starts[i]+1):ends[i]],collapse=''))
    },character(1))
    setNames(seqs,ids)
}

add_sequence_similarity_column = function(g, args, path_alignment, sequence_mode='cds', width_mm=15, genetic_code=1) {
    if (!file.exists(path_alignment) || file.info(path_alignment)$size == 0) return(g)
    if (!sequence_mode %in% c('cds','protein')) stop('Sequence similarity mode must be cds or protein.')
    tips = get_df_tip(g[['tree']]); tips = tips[order(tips[['y']]), ,drop=FALSE]
    if (nrow(tips) < 2) return(g)
    seqs = treevis_similarity_fasta(path_alignment)
    seqs = seqs[names(seqs) %in% as.character(tips[['label']])]
    if (length(seqs) < 2) return(g)
    if (any(nchar(seqs)==0) || length(unique(nchar(seqs))) != 1) {
        stop('Sequence similarity requires a nonempty, equal-length alignment; raw unaligned FASTA is not valid.')
    }
    if (sequence_mode == 'cds') {
        if (any(nchar(seqs) %% 3 != 0)) stop('Amino-acid identity requires a codon alignment with length divisible by three.')
        if (length(genetic_code)!=1 || !is.finite(genetic_code) || genetic_code<1 || genetic_code!=as.integer(genetic_code)) {
            stop('Genetic code must be a positive integer.')
        }
    }
    labels = as.character(tips[['label']])
    pairs = treevis_tip_pairs(nrow(tips))
    values <- treevis_cdskit_request(list(operation="identity",sequences=as.list(seqs),
        pairs=lapply(seq_along(pairs$first), function(k) list(labels[[pairs$first[[k]]]], labels[[pairs$second[[k]]]])),
        mode=sequence_mode,
        genetic_code=genetic_code))
    values <- vapply(values,function(x) if (is.null(x)) NA_real_ else as.numeric(x),numeric(1))
    treevis_pairwise_column(g,args,tips,values,'sequence_similarity',
        'Sequence\nidentity\n(amino acid)',width_mm)
}

add_cis_similarity_column = function(g, args, path_fimo, path_promoters, qvalue=0.01, width_mm=15) {
    if (!file.exists(path_fimo) || !file.exists(path_promoters) ||
        file.info(path_fimo)$size == 0 || file.info(path_promoters)$size == 0) return(g)
    if (length(qvalue)!=1 || !is.finite(qvalue) || qvalue<0 || qvalue>1) stop('FIMO q-value cutoff must be between 0 and 1.')
    tips = get_df_tip(g[['tree']]); tips = tips[order(tips[['y']]), ,drop=FALSE]
    if (nrow(tips)<2) return(g)
    promoters = treevis_similarity_fasta(path_promoters)
    observed = names(promoters)[nchar(promoters)>0]
    if (sum(as.character(tips[['label']]) %in% observed)<2) return(g)
    # Retain legacy commented headers; ignore FIMO's trailing command comments.
    lines = readLines(path_fimo,warn=FALSE)
    headers = which(grepl('sequence[ _]name',lines,ignore.case=TRUE) & grepl('motif|pattern',lines,ignore.case=TRUE))
    if (!length(headers)) stop('FIMO similarity input lacks a recognized table header.')
    lines = lines[seq.int(headers[1],length(lines))]
    lines[1] = sub('^#[[:space:]]*','',lines[1])
    lines = c(lines[1],lines[-1][!grepl('^[[:space:]]*#',lines[-1])])
    hits = read.table(text=paste(lines,collapse='\n'),header=TRUE,sep='\t',quote='',comment.char='',
                      check.names=FALSE,colClasses='character',na.strings=character(0),stringsAsFactors=FALSE)
    names(hits) = gsub('[ -]','_',tolower(names(hits)))
    if ('pattern_name' %in% names(hits)) names(hits)[names(hits)=='pattern_name']='motif_id'
    if (!all(c('sequence_name','motif_id','q_value') %in% names(hits))) {
        stop('Cis similarity requires FIMO sequence_name, motif_id and q-value columns (p-values are not interchangeable).')
    }
    q = suppressWarnings(as.numeric(hits[['q_value']]))
    if (any(!is.finite(q) | q<0 | q>1)) {
        stop('FIMO q-values must be finite numbers between 0 and 1; missing/invalid values are not absence of hits.')
    }
    if (any(!nzchar(hits[['sequence_name']]) | !nzchar(hits[['motif_id']]))) {
        stop('FIMO sequence_name and motif_id must be nonempty.')
    }
    hits = hits[is.finite(q) & q>=0 & q<=qvalue & !is.na(hits[['sequence_name']]) & nzchar(hits[['sequence_name']]) & !is.na(hits[['motif_id']]) & nzchar(hits[['motif_id']]),,drop=FALSE]
    labels = as.character(tips[['label']])
    groups = treevis_tip_groups(labels, hits[['sequence_name']], hits[['motif_id']])
    values = treevis_jaccard(groups, labels %in% observed)
    treevis_pairwise_column(g,args,tips,values,'cis_similarity',
        paste0('Promoter cis\nsimilarity\n(q <= ',format(qvalue,trim=TRUE),')'),width_mm)
}
