# %%
cat('Starting annotation_summary.r\n')
cli_args = commandArgs(trailingOnly=TRUE)
args = cli_args
cat(paste0('Working at: ', getwd(), '\n'))

library(ape, quietly=TRUE)
library(cowplot, quietly=TRUE)
library(ggplot2, quietly=TRUE)
library(phytools, quietly=TRUE)
library(ggtree, quietly=TRUE)
library(svglite, quietly=TRUE)
library(rkftools, quietly=TRUE)
options(stringsAsFactors=FALSE)
script_file_arg = grep('^--file=', commandArgs(), value=TRUE)
if (length(script_file_arg) > 0) {
    script_dir = dirname(normalizePath(sub('^--file=', '', script_file_arg[[1]]), winslash='/', mustWork=FALSE))
} else {
    script_dir = getwd()
}

cat('arguments:\n')
args = rkftools::get_parsed_args(args, print=TRUE)
Rfiles = list.files(file.path(args[['tree_annotation_dir']], 'R'))
for (Rfile in Rfiles) {
    source(file.path(args[['tree_annotation_dir']], 'R', Rfile))
}

font_size = 8
args[['font_size']] = font_size
args[['font_size_factor']] = 0.352777778 # https://stackoverflow.com/questions/17311917/ggplot2-the-unit-of-size
args[['scale_thickness']] = 1
args[['scale_label']] = 'Million years ago'
args[['margins']] = c(0.05,0.05,0.2,0.05) # cm
args[['width']] = 7.2

# %%
path_exists_or_null <- function(path) {
  if (is.null(path) || length(path) == 0) {
    FALSE
  } else {
    file.exists(path)
  }
}

# %%
cat('Reading species tree.\n')
dir_species_tree_summary = file.path(args[['dir_species_tree']], 'species_tree_summary')
dated_sptree_path = file.path(dir_species_tree_summary, 'dated_species_tree.nwk')
undated_sptree_path = file.path(dir_species_tree_summary, 'undated_species_tree.nwk')
if (path_exists_or_null(dated_sptree_path)|path_exists_or_null(undated_sptree_path)) {
    if (path_exists_or_null(dated_sptree_path)) {
        cat(paste0('Reading: ', dated_sptree_path, '\n'))
        tree = read.tree(dated_sptree_path)
    } else if (path_exists_or_null(undated_sptree_path)) {
        cat(paste0('Reading: ', undated_sptree_path, '\n'))
        tree = read.tree(undated_sptree_path)        
    }
    tree[['tip.label']] = sub('_', ' ', tree[['tip.label']])
    tree = ladderize(tree)
    df_out = data.frame(Species=tree[['tip.label']])
    cat(paste0('Number of species in the tree: ', length(tree[['tip.label']]), '\n'))
    tr = list()
    tr[['tree']] = ggtree(tree, size=0.5, layout='rectangular')
    tr[['tree']][['data']][,'tiplab_color'] = 'black'
    tr[['tree']] = tr[['tree']] + coord_cartesian(clip = "off")
    tr[['tree']] = tr[['tree']] + theme(
            rect=element_rect(fill="transparent"),
            plot.margin=unit(args[['margins']], "cm")
        )
    if (path_exists_or_null(dated_sptree_path)) {
        tr[['tree']] = tr[['tree']] + geom_divtime(
            tr[['tree']], 
            y=0.5, 
            linewidth=0.5, 
            fontsize=font_size*args[['font_size_factor']], 
            label='Million years ago'
        )
    } else if (path_exists_or_null(undated_sptree_path)) {
        tr[['tree']] = tr[['tree']] + geom_treescale(
            x=max(tr[['tree']][['data']][['x']])*0.8, 
            y=1.5, 
            fontsize=args[['font_size']]*args[['font_size_factor']], 
            linesize=args[['scale_thickness']], 
            offset=0.25
        )
    }
    tr = add_tiplabel_column(tr, args)
} else {
    message('No species tree is detected. Please make sure one of the following two PATHs is available.')
    message(dated_sptree_path)
    message(undated_sptree_path)
    message('Proceeding without species tree.')
    df_out = data.frame(matrix(NA, nrow=0, ncol=1))
    colnames(df_out) = 'Species'
    tr = NA
}

cat('\n')

# %%
cat('Starting the analysis of species traits.\n')
if (!path_exists_or_null(args[['file_species_trait']])) {
    cat('Trait file not found:', args[['file_species_trait']], '\n')
} else if (any(is.na(tr))) {
    cat('Tree file not found. Skipping the analysis of species traits.\n')
} else {
    trait = read.table(args[['file_species_trait']], sep='\t', header=TRUE, stringsAsFactors=FALSE)
    colnames(trait)[1] = 'label'
    trait_tidy = data.frame(matrix(NA,0,3))
    colnames(trait_tidy) = c('label', 'value', 'trait')
    for (col in colnames(trait)[2:length(colnames(trait))]) {
        tmp = trait[,c('label',col)]
        colnames(tmp) = c('label','value')
        tmp[,'trait'] = col
        trait_tidy = rbind(trait_tidy, tmp)
    }
    g = tr
    rel_widths = c(1.0 ,0.7, 0.5)
    df_tip = get_df_tip(g[['tree']])
    trait_tidy[['label']] = factor(sub('_', ' ', trait_tidy[['label']]), levels=levels(df_tip[['label']]))
    df_tip = merge(df_tip, trait_tidy, all.x=TRUE, by='label', sort=FALSE)
    p = ggplot(df_tip)
    p = p + theme(axis.text.y=element_blank())
    p = p + geom_tile(aes(x=trait, y=label, fill=value))
    p = p + theme(
        axis.text=element_text(size=font_size),
        axis.title=element_text(size=font_size),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.title=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        legend.justification=c(0,0),
        legend.text=element_text(size=font_size),
        legend.position=c(0.1,0.1),
        legend.key.size=unit(0.4, 'lines'),
        legend.box.just='center',
        #rect=element_rect(fill="transparent"),
        plot.margin=unit(args[['margins']], "cm")
    )
    g[['trait']] = p

    base_width = as.numeric(args[['width']])
    cp = cowplot::plot_grid(plotlist=g, nrow=1, align='h', axis='bt', rel_widths=rel_widths)
    cat('Writing the plot pdf and svg.\n')
    height = nrow(trait) / 8
    extensions = c('.pdf','.svg')
    for (extension in extensions) {
        cowplot::save_plot(
            filename=paste0('annotation_summary_trait', extension),
            plot=cp,
            nrow=1,
            base_height=max(1.5,height),
            base_width=base_width,
            units='in',
            dpi=300,
            limitsize=FALSE,
            bg = 'transparent'
        )
    }
}
cat('\n')

# %%
generate_busco_summary = function(dir_busco, outbase, df_out, tr = NA, font_size = 8, args = list()) {
    cat('Starting the analysis of BUSCO outputs:', dir_busco, '\n')
    if (!path_exists_or_null(dir_busco)) {
        cat('BUSCO directory not found:', dir_busco, '\n')
        return(df_out)
    }    
    files = list.files(dir_busco)
    cat(paste0('Number of BUSCO full tables: ', length(files), '\n'))
    if (length(files)==0) {
        cat('Skipping the analysis of BUSCO tables.\n')
        return(df_out)
    }
    df = data.frame()
    for (file in files) {
        file_path = file.path(dir_busco, file)
        sp = sub('_busco.*', '', sub('\\.busco.*', '', sub('_', ' ', file)))
        all_lines <- readLines(file_path)
        header_line <- grep("^# Busco id", all_lines, value = TRUE)
        header <- unlist(strsplit(header_line, "\t"))
        header = tolower(header)
        header = sub('# ', '', header)
        header = sub(' ', '_', header)
        data_lines <- grep("^[^#]", all_lines, value = TRUE)
        tmp = read.table(text=data_lines, sep='\t', header=FALSE, comment.char='#', col.names=header, fill=TRUE, quote='')
        tmp[,'label'] = sp
        df = rbind(df, tmp)
    }
    df2 = unique(df[,c('status','label','busco_id')])
    df2 = aggregate(df2, by=df2[,c('status','label')], FUN=length)
    df2 = df2[,c('status','label','busco_id')]
    df2[(df2[['status']]=='Complete'),'status'] = 'Single'
    colnames(df2) = c('status','label','count')
    num_busco_gene = aggregate(df2[,'count'], by=list(df2[['label']]), FUN=sum)[1,2]
    
    df3 = reshape(df2, v.names='count', timevar='status', idvar='label', direction='wide')
    df3[is.na(df3)] = 0
    for (col in colnames(df3)[2:length(colnames(df3))]) {
        df3[,col] = as.integer(df3[[col]])
    }
    colnames(df3) = sub('count.', '', colnames(df3))
    colnames(df3) = sub('label', 'Species', colnames(df3))

    expected_statuses = c('Single', 'Duplicated', 'Fragmented', 'Missing')
    for (status in expected_statuses) {
        if (!(status %in% colnames(df3))) {
            df3[,status] = 0L
        }
    }
    df3[,'Total'] = apply(df3[,expected_statuses], 1, sum)
    safe_percent <- function(numerator, denominator) {
        out <- rep(NA_real_, length(denominator))
        valid <- denominator > 0
        out[valid] <- round(numerator[valid] / denominator[valid] * 100, digits=1)
        return(out)
    }
    pc = as.character(safe_percent(df3[['Single']] + df3[['Duplicated']], df3[['Total']]))
    ps = as.character(safe_percent(df3[['Single']], df3[['Total']]))
    pd = as.character(safe_percent(df3[['Duplicated']], df3[['Total']]))
    pf = as.character(safe_percent(df3[['Fragmented']], df3[['Total']]))
    pm = as.character(safe_percent(df3[['Missing']], df3[['Total']]))
    df3[,'Summary'] = paste0('C:',pc,'%[S:',ps,'%,D:',pd,'%],F:',pf,'%,M:',pm,'%,n:',df3[['Total']])
    df3 = df3[,c('Species','Summary','Single','Duplicated','Fragmented','Missing','Total')]
    colnames(df3) = paste0(outbase, '_', tolower(colnames(df3)))
    colnames(df3)[1] = 'Species'
    df_out = merge(df_out, df3, by='Species', all=TRUE, sort=TRUE)
    levels = c('Missing','Fragmented','Duplicated','Single')
    colors = c('gray80','gray40','firebrick','black')
    names(colors) = levels
    df2[,'status'] = factor(df2[,'status'], levels=levels)
    
    if (any(is.na(tr))) {
        cat('The BUSCO stacked barplot will be generated without a species tree.\n')
        g = list()
        rel_widths = c(1.0)
        levels = rev(sort(unique(df2[['label']])))
        df2[,'label'] = factor(df2[,'label'], levels=levels)
        p = ggplot(df2, aes(fill=status, y=label, x=count))
    } else {
        cat('The BUSCO stacked barplot will be generated with a species tree.\n')
        g = tr
        rel_widths = c(1.0 ,0.7, 1.0)
        df_tip = get_df_tip(g[['tree']])
        df_tip = merge(df_tip, df2, all.x=TRUE, by='label', sort=FALSE)
        p = ggplot(df_tip, aes(fill=status, y=label, x=count))
    }
    
    p = p + geom_bar(position='stack', stat='identity')
    p = p + scale_fill_manual(values=colors)
    p = p + theme_linedraw(base_size=font_size)
    p = p + xlab('Number of BUSCO single-copy genes')
    p = p + guides(colour=guide_legend(nrow=6, byrow=FALSE))
    p = p + scale_x_continuous(limits=c(0, num_busco_gene), sec.axis=sec_axis(trans=~./num_busco_gene*100, name='%'))
    p = p + theme(
        axis.text.x=element_text(size=font_size),
        axis.title=element_text(size=font_size),
        axis.title.y=element_blank(), 
        axis.ticks.y=element_blank(),
        legend.title=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        legend.justification=c(0,0),
        legend.text=element_text(size=font_size),
        legend.position=c(0.1,0.1),
        legend.key.size=unit(0.4, 'lines'), 
        legend.box.just='center',
        #rect=element_rect(fill="transparent"),
        plot.margin=unit(args[['margins']], "cm")
    )
    if (any(is.na(tr))) {
        p = p + theme(axis.text.y=element_text(size=font_size))
    } else {
        p = p + theme(axis.text.y=element_blank())
    }
    g[['barplot']] = p
    
    base_width = as.numeric(args[['width']])
    cp = cowplot::plot_grid(plotlist=g, nrow=1, align='h', axis='bt', rel_widths=rel_widths)
    cat('Writing the plot pdf and svg.\n')
    height = nrow(df3) / 8
    extensions = c('.pdf','.svg')
    for (extension in extensions) {
        cowplot::save_plot(
            filename=paste0(outbase, extension), 
            plot=cp, 
            nrow=1,
            base_height=max(1.5,height),
            base_width=base_width,
            units='in', 
            dpi=300, 
            limitsize=FALSE,
            bg = 'transparent'
        )
    }
    return(df_out)
    cat('\n')
}

df_out = generate_busco_summary(args[['dir_species_cds_busco']], 'busco_cds', df_out, tr=tr, font_size=font_size, args=args)
df_out = generate_busco_summary(args[['dir_species_genome_busco']], 'busco_genome', df_out, tr=tr, font_size=font_size, args=args)

# %%
calculate_N50 <- function(data) {
  # Ensure the data frame contains a "length" column
  if (!"length" %in% colnames(data)) {
    stop("The input data frame must have a 'length' column.")
  }
  
  # Extract the lengths and sort them in descending order
  lengths <- sort(data$length, decreasing = TRUE)
  
  # Calculate the total length of all sequences
  total_length <- sum(lengths)
  
  # Determine the N50 value
  cumulative_sum <- cumsum(lengths)
  n50_index <- which(cumulative_sum >= total_length / 2)[1]
  n50 <- lengths[n50_index]
  
  return(n50)
}

generate_fx2tab_summary <- function(dir_fx2tab, outbase, df_out, tr = NA, font_size = 8, args = list()) {
    if (!path_exists_or_null(dir_fx2tab)) {
        cat('fx2tab directory not found:', dir_fx2tab, '\n')
        return(df_out)
    }
    files <- list.files(dir_fx2tab)
    cat(paste0('Number of fx2tab tables: ', length(files), '\n'))
    if (length(files) == 0) {
        cat('Skipping the analysis of fx2tab tables.\n')
        return(df_out)
    }
    cat('Starting the analysis of fx2tab tables.\n')
    df <- data.frame()
    for (file in files) {
        file_path <- file.path(dir_fx2tab, file)
        sp <- sub('_.*', '', sub('\\..*', '', sub('_', ' ', file)))
        tmp <- read.table(file_path, sep = '\t', header = TRUE, comment.char = '', fill = TRUE, quote = '')
        
        # Calculate weighted GC %
        weighted_gc <- sum(tmp[['GC']] * tmp[['length']], na.rm = TRUE) / sum(tmp[['length']], na.rm = TRUE)
        
        tmp2 <- data.frame(label = sp)
        tmp2[['num_seq']] <- nrow(tmp)
        tmp2[['total_length']] <- sum(tmp[['length']], na.rm = TRUE)
        tmp2[['mean_length']] <- mean(tmp[['length']], na.rm = TRUE)
        tmp2[['weighted_gc']] <- weighted_gc
        tmp2[['mean_gc_skew']] <- mean(tmp[['GC.Skew']], na.rm = TRUE)
        tmp2[['N50']] <- calculate_N50(tmp)
        rownames(tmp2) <- NULL
        df <- rbind(df, tmp2)
    }
    
    write.table(df, file = paste0(outbase, '.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
    df2 <- df
    colnames(df2) <- sub('label', 'Species', colnames(df2))
    df_out <- merge(df_out, df2, by = 'Species', all = TRUE, sort = TRUE)

    if (!any(is.na(tr))) {
        cat('Plotting with species tree.\n')
        g <- tr
        df_tip <- get_df_tip(g[['tree']])
        df_tip <- merge(df_tip, df, all.x = TRUE, by = 'label', sort = FALSE)
    } else {
        cat('Plotting without species tree.\n')
        df_tip <- df
        g <- list()
        g[['tree']] <- ggplot(df_tip, aes(y = label)) + 
            geom_blank() + 
            scale_y_discrete() + 
            theme_void() + 
            theme(
                axis.text.y = element_text(size = font_size),
                axis.ticks.y = element_blank(),
                axis.title.y = element_blank(),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank()
            )
    }
    
    xs <- c(
        'num_seq',
        'total_length',
        'N50',
        'weighted_gc'
    )
    xlabels <- c(
        '# of seqs',
        'Total length (bp)',
        'N50 (bp)',
        'GC content (%)'
    )
    
    for (i in 1:length(xs)) {
        x <- xs[i]
        if (startsWith(x, 'num_')) {
            xmax <- max(df_tip[[x]], na.rm = TRUE) * 1.1
            df_tip[['xtext']] <- format(df_tip[[x]], big.mark = ",", scientific = FALSE)
        } else if (x == 'weighted_gc') {
            xmax <- 100
            df_tip[['xtext']] <- format(round(df_tip[[x]], digits = 1), nsmall = 1)
        } else {
            xmax <- max(df_tip[[x]], na.rm = TRUE) * 1.1
            df_tip[['xtext']] <- format(df_tip[[x]], big.mark = ",", scientific = FALSE)
        }
        
        p <- ggplot(df_tip, mapping = aes(y = label)) +
            geom_bar(mapping = aes(x = !!rlang::sym(x)), position = 'stack', stat = 'identity', fill = 'black') +
            geom_text(mapping = aes(label = xtext), x = xmax / 50, size = font_size * args[['font_size_factor']], 
                      hjust = 0, color = 'gray50') +
            theme_linedraw(base_size = font_size) +
            xlab(xlabels[i]) +
            coord_cartesian(xlim = c(0, xmax)) +
            theme(
                axis.text = element_text(size = font_size),
                axis.title = element_text(size = font_size),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                axis.ticks.y = element_blank(),
                legend.title = element_blank(),
                panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank(),
                legend.justification = c(0, 0),
                legend.text = element_text(size = font_size),
                legend.position = c(0.1, 0.1),
                legend.key.size = unit(0.4, 'lines'),
                legend.box.just = 'center',
                plot.margin = unit(args[['margins']], "cm")
            )
        
        g[[x]] <- p
    }

    base_width <- as.numeric(args[['width']])
    cp <- cowplot::plot_grid(plotlist = g, nrow = 1, align = 'h', axis = 'bt')
    cat('Writing the plot pdf and svg.\n')
    height <- nrow(df) / 8
    extensions = c('.pdf', '.svg')
    for (extension in extensions) {
        cowplot::save_plot(
            filename = paste0(outbase, extension),
            plot = cp, 
            nrow = 1,
            base_height = max(1.5, height),
            base_width = base_width,
            units = 'in', 
            dpi = 300, 
            limitsize = FALSE,
            bg = 'transparent'
        )
    }
    return(df_out)
    cat('\n')
}

df_out = generate_fx2tab_summary(dir_fx2tab=args[['dir_species_cds_fx2tab']], outbase='fx2tab_cds', df_out, tr=tr, font_size=font_size, args=args)
df_out = generate_fx2tab_summary(dir_fx2tab=args[['dir_species_genome_fx2tab']], outbase='fx2tab_genome', df_out, tr=tr, font_size=font_size, args=args)

# %%
generate_annotation_summary <- function(dir_species_annotation, outbase, tr = NA, font_size = 8, args = list()) {
    if (!path_exists_or_null(dir_species_annotation)) {
        cat("Annotation directory not found:", dir_species_annotation, "\n")
        return()
    }
    files <- list.files(dir_species_annotation)
    cat(paste0("Number of annotation tables: ", length(files), "\n"))
    if (length(files) == 0) {
        cat("Skipping the analysis of annotation tables.\n")
        return()
    }
    
    cat("Starting the analysis of annotation tables.\n")
    df <- data.frame()
    old_cols <- c('gene_id', 'orthogroup', 'sprot_best', 'busco_id')
    new_cols <- c('num_gene', 'num_orthogroup', 'num_sprot', 'num_busco')
    
    for (file in files) {
        file_path <- file.path(dir_species_annotation, file)
        sp <- sub('_.*', '', sub('\\..*', '', sub('_', ' ', file)))
        tmp <- read.table(file_path, sep = '\t', header = TRUE, comment.char = '#', fill = TRUE, quote = '')
        tmp2 <- apply(tmp, 2, function(x) { sum((!is.na(x)) & (x != '')) })
        tmp2 <- tmp2[old_cols]
        tmp2[is.na(tmp2)] <- 0
        names(tmp2) <- new_cols
        tmp2['label'] <- sp
        tmp2 <- tmp2[c('label', new_cols)]
        tmp2 <- t(data.frame(tmp2))
        rownames(tmp2) <- NULL
        df <- rbind(df, tmp2)
    }
    
    # Save the raw counts to a separate file
    write.table(df, file = paste0(outbase, "_raw_counts.tsv"), sep = '\t', row.names = FALSE, quote = FALSE)
    
    # Add percentage columns
    for (col in colnames(df)) {
        if (col == 'label') next
        df[, col] <- as.integer(df[, col])
        if (col == 'num_gene') next
        new_col <- sub('num_', 'percent_', col)
        df[, new_col] <- df[[col]] / df[['num_gene']] * 100
    }
    
    # Save the percentage table to a separate file
    write.table(df, file = paste0(outbase, "_percentages.tsv"), sep = '\t', row.names = FALSE, quote = FALSE)
    
    if (length(files) > 0) {
        if (!any(is.na(tr))) {
            cat("Plotting with species tree.\n")
            raw_g <- tr
            percent_g <- tr
            df_tip <- get_df_tip(raw_g[['tree']])
            df_tip <- merge(df_tip, df, all.x = TRUE, by = 'label', sort = FALSE)
        } else {
            cat("Plotting without species tree.\n")
            df_tip <- df
            raw_g <- list()
            raw_g[['tree']] <- ggplot(df_tip, aes(y = label)) + geom_blank() + scale_y_discrete() + theme_void() + 
                theme(
                    axis.text.y = element_text(size = font_size),
                    axis.ticks.y = element_blank(),
                    axis.title.y = element_blank(),
                    axis.text.x = element_blank(),
                    axis.ticks.x = element_blank()
                )
            percent_g <- raw_g
        }
        
        # Raw count plots
        raw_xs <- c(
            'num_gene',
            'num_orthogroup',
            'num_sprot',
            'num_busco'
        )
        raw_xlabels <- c(
            '# genes',
            '# OG',
            '# UNIPROT',
            '# BUSCO'
        )
        
        for (i in 1:length(raw_xs)) {
            x <- raw_xs[i]
            xmax <- max(df_tip[[x]], na.rm = TRUE) * 1.1
            df_tip[['xtext']] <- format(df_tip[[x]], big.mark = ",", scientific = FALSE)
            
            p <- ggplot(df_tip, mapping = aes(y = label)) +
                geom_bar(mapping = aes(x = !!rlang::sym(x)), position = 'stack', stat = 'identity', fill = 'black') +
                geom_text(mapping = aes(label = xtext), x = xmax / 50, size = font_size * args[['font_size_factor']], 
                          hjust = 0, color = 'gray50') +
                theme_linedraw(base_size = font_size) +
                xlab(raw_xlabels[i]) +
                coord_cartesian(xlim = c(0, xmax)) +
                theme(
                    axis.text = element_text(size = font_size),
                    axis.title = element_text(size = font_size),
                    axis.title.y = element_blank(),
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    legend.title = element_blank(),
                    panel.grid.major.y = element_blank(),
                    panel.grid.minor.y = element_blank(),
                    panel.grid.major.x = element_blank(),
                    panel.grid.minor.x = element_blank(),
                    legend.justification = c(0, 0),
                    legend.text = element_text(size = font_size),
                    legend.position = c(0.1, 0.1),
                    legend.key.size = unit(0.4, 'lines'),
                    legend.box.just = 'center',
                    plot.margin = unit(args[['margins']], "cm")
                )
            
            raw_g[[x]] <- p
        }
        
        # Combine and save raw count plots
        base_width <- as.numeric(args[['width']])
        raw_cp <- cowplot::plot_grid(plotlist = raw_g, nrow = 1, align = 'h', axis = 'bt')
        cat("Writing the raw count combined plot pdf and svg.\n")
        for (ext in c(".pdf", ".svg")) {
            cowplot::save_plot(
                filename = paste0(outbase, "_raw_counts", ext),
                plot = raw_cp,
                nrow = 1,
                base_height = max(1.5, nrow(df) / 8),
                base_width = base_width,
                units = 'in',
                dpi = 300,
                limitsize = FALSE,
                bg = 'transparent'
            )
        }
        
        # Percentage plots
        percent_xs <- c(
            'percent_orthogroup',
            'percent_sprot',
            'percent_busco'
        )
        percent_xlabels <- c(
            '% OG',
            '% UNIPROT',
            '% BUSCO'
        )
        
        for (i in 1:length(percent_xs)) {
            x <- percent_xs[i]
            xmax <- 100
            df_tip[['xtext']] <- format(round(df_tip[[x]], digits = 1), nsmall = 1)
            
            p <- ggplot(df_tip, mapping = aes(y = label)) +
                geom_bar(mapping = aes(x = !!rlang::sym(x)), position = 'stack', stat = 'identity', fill = 'black') +
                geom_text(mapping = aes(label = xtext), x = xmax / 50, size = font_size * args[['font_size_factor']], 
                          hjust = 0, color = 'gray50') +
                theme_linedraw(base_size = font_size) +
                xlab(percent_xlabels[i]) +
                coord_cartesian(xlim = c(0, xmax)) +
                theme(
                    axis.text = element_text(size = font_size),
                    axis.title = element_text(size = font_size),
                    axis.title.y = element_blank(),
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    legend.title = element_blank(),
                    panel.grid.major.y = element_blank(),
                    panel.grid.minor.y = element_blank(),
                    panel.grid.major.x = element_blank(),
                    panel.grid.minor.x = element_blank(),
                    legend.justification = c(0, 0),
                    legend.text = element_text(size = font_size),
                    legend.position = c(0.1, 0.1),
                    legend.key.size = unit(0.4, 'lines'),
                    legend.box.just = 'center',
                    plot.margin = unit(args[['margins']], "cm")
                )
            
            percent_g[[x]] <- p
        }
        
        # Combine and save percentage plots
        percent_cp <- cowplot::plot_grid(plotlist = percent_g, nrow = 1, align = 'h', axis = 'bt')
        cat("Writing the percentage combined plot pdf and svg.\n")
        for (ext in c(".pdf", ".svg")) {
            cowplot::save_plot(
                filename = paste0(outbase, "_percentages", ext),
                plot = percent_cp,
                nrow = 1,
                base_height = max(1.5, nrow(df) / 8),
                base_width = base_width,
                units = 'in',
                dpi = 300,
                limitsize = FALSE,
                bg = 'transparent'
            )
        }
    }
    cat("\n")
}


generate_annotation_summary(dir_species_annotation=args[['dir_species_annotation']], outbase='annotation', tr=tr, font_size=font_size, args=args)

# %%
if ((!any(is.na(tr)))&(path_exists_or_null(args[['file_orthogroup_gene_count']]))) {
    cat('Starting the analysis of orthogroup tables.\n')
    if (args[['min_og_species']]=='auto') {
        min_og_species = max(2, ceiling(length(tree[['tip.label']])/2))
        cat(paste0('--min_og_species=auto: orthogroups with genes from>=50% of species will be included (>=', min_og_species, ' species).\n'))
    } else {
        min_og_species = as.integer(args[['min_og_species']])
    }
    cat(paste('Orthogroups with genes from a minimum of', min_og_species, 'species are used for PCA.\n'))
    gc = read.table(args[['file_orthogroup_gene_count']], sep='\t', header=TRUE, fill=TRUE, quote='')
    rownames(gc) = gc[['Orthogroup']]
    gc = gc[,!(colnames(gc) %in% c('Orthogroup', 'Total'))]
    if (ncol(gc)<min_og_species) {
        cat(paste('There are only ', ncol(gc), 'species in the genecount table. ', min_og_species, 'species are required for NMDS. Skipping.\n'))
    } else {
        colnames(gc) = sub('_', ' ', colnames(gc))
        gc[,] = (gc[,] > 0)
        gc = gc[(apply(gc, 1, function(x){sum(x) >= min_og_species})), , drop=FALSE] # Remove single-species orthogroups
        cat(paste0('Number of analyzed orthogroups: ', nrow(gc), '\n'))
        if (nrow(gc) < 2) {
            cat('Too few orthogroups remained after filtering. Skipping orthogroup PCA.\n')
        } else {
            sd_by_species = apply(gc, 2, function(x){sd(as.numeric(x), na.rm=TRUE)})
            keep_species = names(sd_by_species)[(!is.na(sd_by_species)) & (sd_by_species > 0)]
            if (length(keep_species) < 2) {
                cat('Too few species with variable orthogroup content for PCA. Skipping orthogroup PCA.\n')
            } else {
                gc = gc[,keep_species,drop=FALSE]
                tc_dist_matrix = cor(gc, method='pearson', use='pairwise.complete.obs')
                tc_dist_matrix[is.na(tc_dist_matrix)] = 0
                diag(tc_dist_matrix) = 1
                pca = prcomp(tc_dist_matrix)
                xlabel = paste0("PC 1 (", round(summary(pca)[['importance']][2, 1] * 100, digits = 1), "%)")
                ylabel = paste0("PC 2 (", round(summary(pca)[['importance']][2, 2] * 100, digits = 1), "%)")
                pca_xy = pca[['x']][,c('PC1','PC2'),drop=FALSE]
                shared_species = intersect(tree[['tip.label']], rownames(pca_xy))
                if (length(shared_species) < 3) {
                    cat('Fewer than 3 shared species between tree and PCA matrix. Skipping orthogroup PCA.\n')
                } else {
                    tree_pca = keep.tip(tree, shared_species)
                    pca_xy = pca_xy[tree_pca[['tip.label']],,drop=FALSE]

                    pdf('annotation_summary_orthogroup_PCA.pdf', height = 3.2, width = 3.6, fonts = "Helvetica", pointsize = font_size)
                    phylomorphospace(tree_pca, pca_xy, label="horizontal", xlab=xlabel, ylab=ylabel, node.by.map=TRUE, las=1)
                    dev.off()
                }
            }
        }
    }
} else {
    cat('Skipping the analysis of orthogroup tables.\n')
}
cat('\n')

# %%
if (all(c('num_seq', 'num_gene') %in% colnames(df_out))) {
    df_out[['num_seq']] = NULL
}
cat(paste('Column names of output table:', paste(colnames(df_out), collapse=' '), '\n'))
write.table(df_out, file='annotation_summary.tsv', sep='\t', row.names=FALSE, quote=FALSE)
cat('Ending annotation_summary.r\n')

# %%
