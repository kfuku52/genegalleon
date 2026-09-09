# Regression coverage for tree-plot legend, label, and width customizations.
suppressPackageStartupMessages(library(genegalleon.treevis))
suppressPackageStartupMessages(library(ggplot2))
output <- tempfile('treevis-legends-')
dir.create(output)
# Keep any device opened during grob measurements out of the worktree.
pdf(file.path(output, 'measurements.pdf'))

phy <- ape::read.tree(text='(((g1:1,g2:1):1,g3:2):1,(g4:1,g5:1):2);')
base <- list(tree=ggtree::ggtree(phy))
args <- list(font_size=8, font_size_factor=0.352777778, margins=rep(0.1,4),
             nodepoint_size=2, node_colors=c(D='red',R='orange',S='blue',H='darkturquoise'))

# Legend labels must retain the original category-to-color mapping, even when
# only a subset is present. No inferred event absent from the tree is added.
for (events in list(c('D','S','H','R'), c('D','S','D','S'), rep('H',4))) {
    p <- base$tree
    p$data$node_category <- NA_character_
    p$data$node_category[!p$data$isTip] <- events
    p <- add_node_points(p,args)
    scale <- ggplot_build(p)$plot$scales$get_scales('colour')
    breaks <- scale$get_breaks()
    breaks <- breaks[!is.na(breaks)]
    stopifnot(setequal(breaks, unique(events)))
    expected <- c(D='Duplication', S='Speciation', H='Transfer', R='Retrotransposition')
    stopifnot(identical(unname(scale$get_labels(breaks)), unname(expected[breaks])),
              identical(unname(scale$map(breaks)), unname(args$node_colors[breaks])))
    stopifnot(inherits(cowplot::get_legend(p), 'gtable'))
}

# Query styling includes legacy marker spellings and all-missing/all-hit cases.
key <- 'categorical,query_marker,Query'
for (values in list(c('Best hit','-','best_blast_hit',NA,''), rep('-',5), rep('Best hit',5))) {
    input <- base
    input$tree$data$query_marker <- NA_character_
    input$tree$data$query_marker[input$tree$data$isTip] <- values
    g <- add_categorical_column(input,args,key,'query_marker','Query')
    p <- g[[key]]
    built <- ggplot_build(p)
    stopifnot(is.null(cowplot::get_legend(p)), p$labels$x == 'Query\nbest\nhit')
    expected_values <- normalize_category_text(values)
    expected_fill <- ifelse(expected_values == 'Best hit','black','#e6e6e6')
    # Compare by tip label, rather than assuming the tree's plotted row order.
    scale <- built$plot$scales$get_scales('fill')
    stopifnot(identical(unname(scale$map(expected_values)), unname(expected_fill)))
    custom <- add_categorical_column(input,args,key,'query_marker','User label')[[key]]
    stopifnot(custom$labels$x == 'User label')
}

# Ordinary categorical panels keep their own palettes and graphical legends.
input <- base
input$tree$data$other_category <- 'family'
g <- add_categorical_column(input,args,'categorical,other,Other','other_category','Other')
stopifnot(g[['categorical,other,Other']]$labels$x == 'Other',
          inherits(cowplot::get_legend(g[['categorical,other,Other']]), 'gtable'))
g[[key]] <- g[['categorical,other,Other']]
w <- get_rel_widths(g,'')
stopifnot(abs(w[key] / w['categorical,other,Other'] - 0.5) < 1e-12)
overridden <- get_rel_widths(g,'categorical,0.8')
stopifnot(all(overridden[grepl('^categorical,',names(overridden))] == 0.8))

# Rich species labels must render as text, including names with punctuation.
for (prefix in c('Arabidopsis_thaliana_', 'Species_a_subsp_b_', 'Species_a_x_Species_b_')) {
    label <- treevis_ortholog_axis_label(prefix)
    p <- ggplot(data.frame(x=1,y=1),aes(x,y)) + geom_blank() + xlab(label)
    gt <- ggplotGrob(p)
    stopifnot(is.expression(label), inherits(gt,'gtable'))
}

# A long domain legend is limited to the available width. Exercise the
# actual driver and confirm that fitting keys preserves panel coordinates.
file_arg <- grep('^--file=', commandArgs(), value=TRUE)[1]
repo <- normalizePath(file.path(dirname(sub('^--file=', '', file_arg)), '..','..'))
driver <- file.path(repo,'workflow/support/stat_branch2tree_plot.r')
branch <- data.frame(branch_id=c(4,2,0,1,3), parent=c(-999,4,2,2,4),
    sister=c(-999,3,1,0,2), child1=c(2,0,NA,NA,NA), child2=c(3,1,NA,NA,NA),
    node_name=c('root','n4','g1','g2','g3'), bl_rooted=c(0,1,1,1,2),
    so_event=c('S','S','L','L','L'), so_event_parent='S', expression_A_rep1=1:5)
write.table(branch,file.path(output,'branch.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
rps <- data.frame(qacc=rep('g1',7),sacc=paste0('PF',1:7),
    stitle=paste0('PF',1:7,',Long_protein_domain_name_',1:7),
    qlen=300,slen=300,qstart=seq(1,241,40),qend=seq(30,270,40))
write.table(rps,file.path(output,'rps.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
runner <- file.path(output,'render.R')
writeLines(c(sprintf('source(%s)',encodeString(driver,quote='"')),
    "stopifnot(length(cp$layers) == length(g))",
    "stopifnot(length(layout_mm$legend_entries$domain) < 7)",
    "domain_position <- match('domain',names(g))",
    "expected_x <- sum(rel_widths[seq_len(domain_position-1)]) / sum(rel_widths)",
    "last <- cp$layers[[length(cp$layers)]]$geom_params",
    "stopifnot(abs(last$xmin - expected_x) < 1e-12)",
    "stopifnot(abs((last$xmax-last$xmin) - rel_widths[domain_position]/sum(rel_widths)) < 1e-12)",
    "ggplot2::ggsave('domain-overflow.png',cp,width=9,height=4,dpi=120,bg='white')"), runner)
oldwd <- setwd(output)
logs <- system2('Rscript', c(shQuote(runner), '--stat_branch=branch.tsv',
    '--panel_widths_mm=tree:60,domain:45,pointplot:90', '--panel1=tree,bl_rooted,no,no,L',
    '--panel2=domain,rps.tsv','--panel3=pointplot,no,rel,_,expression_',
    '--show_branch_id=no','--event_method=species_overlap','--species_color_table=PLACEHOLDER',
    '--pie_chart_value_transformation=identity','--max_delta_intron_present=-0.5',
    '--long_branch_display=no'), stdout=TRUE, stderr=TRUE)
setwd(oldwd)
if (!is.null(attr(logs,'status'))) stop(paste(logs,collapse='\n'))
stopifnot(file.info(file.path(output,'stat_branch2tree_plot.pdf'))$size > 1000)
dev.off()
cat('Treevis event, query, species-label, and domain-overflow tests passed.\n')

# Intron overlays are opt-in even when annotation is available.
input <- base
input$tree$data$intron_positions <- '30;60'
default_domain <- add_protein_domain_column(input,args,rps)$domain
marked_domain <- add_protein_domain_column(input,args,rps,show_introns=TRUE)$domain
has_marks <- function(p) any(vapply(p$layers,function(layer)
    isTRUE(attr(layer,'treevis_intron_marks')),logical(1)))
stopifnot(!has_marks(default_domain),has_marks(marked_domain))
