suppressPackageStartupMessages(library(genegalleon.treevis))
suppressPackageStartupMessages(library(ggplot2))
tree = ape::read.tree(text='((a:1,b:1):1,c:2,d:2);')
g = list(tree=ggtree::ggtree(tree))
args = list(font_size=6, margins=rep(0.1,4))
f = tempfile(fileext='.tsv')
syn = data.frame(node_name=c('a','a','a','b','b','c'), offset=c(-1,1,20,-1,1,1),
                 group_id=c('A','A','B','A','C','Z'))
write.table(syn,f,sep='\t',quote=FALSE,row.names=FALSE)
p = add_synteny_similarity_column(g,args,f)$synteny_similarity
stopifnot(attr(p,'treevis_width_mm') == 15,
          grepl('20 genes', p$labels$x), nrow(p$data) == 6*4)
tips = get_df_tip(g$tree); tips=tips[order(tips$y),]
pairs=utils::combn(as.character(tips$label),2)
vals=p$data$similarity[seq(1,nrow(p$data),4)]
for (i in seq_along(vals)) {
    pair=pairs[,i]
    if ('d' %in% pair) stopifnot(is.na(vals[i]))
    else if (setequal(pair,c('a','b'))) stopifnot(abs(vals[i]-1/3)<1e-10)
    else stopifnot(vals[i]==0)
}
# Group IDs, duplicate copies, row order and irrelevant tips do not alter similarity.
syn2=rbind(syn,syn[1,],data.frame(node_name='absent',offset=1,group_id='A'))
syn2$group_id=paste0('renamed_',syn2$group_id)
write.table(syn2[nrow(syn2):1,],f,sep='\t',quote=FALSE,row.names=FALSE)
p2=add_synteny_similarity_column(g,args,f)$synteny_similarity
stopifnot(identical(vals,p2$data$similarity[seq(1,nrow(p2$data),4)]))
# Search range is independent of the neighboring-gene display range.
p5=add_synteny_similarity_column(g,args,f,search_window=5)$synteny_similarity
ab=which(apply(pairs,2,function(x) setequal(x,c('a','b'))))
stopifnot(p5$data$similarity[(ab-1)*4+1] == 1/2)
built=ggplot_build(p)
stopifnot(identical(built$layout$panel_params[[1]]$y.range,
                   ggplot_build(g$tree)$layout$panel_params[[1]]$y.range),
          any(built$data[[1]]$fill=='grey80'))
stopifnot(identical(add_synteny_similarity_column(g,args,paste0(f,'.missing')),g))
write.table(syn[FALSE,],f,sep='\t',quote=FALSE,row.names=FALSE)
stopifnot(identical(add_synteny_similarity_column(g,args,f),g))
unlink(f)
cat('Synteny similarity values, missing data, window and geometry passed.\n')

# Run the workflow's actual PDF driver; fixed column adds exactly 15 mm.
test_file=sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE)[1])
repo=normalizePath(file.path(dirname(test_file),'..','..'))
outdir=tempfile('synteny-similarity-driver-'); dir.create(outdir)
branch=data.frame(branch_id=c(4,2,0,1,3),parent=c(-999,4,2,2,4),
    sister=c(-999,3,1,0,2),child1=c(2,0,NA,NA,NA),child2=c(3,1,NA,NA,NA),
    node_name=c('root','n4','a','b','c'),bl_rooted=c(0,1,1,1,2),
    so_event=c('S','S','L','L','L'),so_event_parent='S')
write.table(branch,file.path(outdir,'branch.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(syn,file.path(outdir,'syn.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
writeLines(c('>a','GCTAAA','>b','GCCGAA','>c','GCCAAA'),file.path(outdir,'aligned.fa'))
writeLines(c('motif_id\tsequence_name\tq-value','M1\ta\t0.001','M1\tb\t0.001'),file.path(outdir,'fimo.tsv'))
oldwd=setwd(outdir)
for (panel_count in c(0,1,3)) {
    has_input=panel_count>=1
    result=system2('Rscript',c(file.path(repo,'workflow/support/stat_branch2tree_plot.r'),
        '--stat_branch=branch.tsv','--panel_widths_mm=tree:60',
        '--panel1=tree,bl_rooted,no,no,L','--panel2=tiplabel',
        '--panel3=ortholog,none_,missing.nwk',
        paste0('--panel4=synteny_similarity,',if(has_input) 'syn.tsv' else 'missing.tsv',',20,15'),
        paste0('--panel5=cis_similarity,',if(panel_count==3) 'fimo.tsv' else 'missing.tsv',',aligned.fa,0.01,15'),
        paste0('--panel6=sequence_similarity,',if(panel_count==3) 'aligned.fa' else 'missing.fa',',cds,15'),
        '--show_branch_id=no','--event_method=species_overlap',
        '--species_color_table=PLACEHOLDER','--pie_chart_value_transformation=identity',
        '--max_delta_intron_present=-0.5','--long_branch_display=no'),stdout=TRUE,stderr=TRUE)
    if (!is.null(attr(result,'status'))) stop(paste(result,collapse='\n'))
    read_width="import re; data=open('stat_branch2tree_plot.pdf','rb').read(); m=re.search(rb'/MediaBox\\s*\\[\\s*0\\s+0\\s+([0-9.]+)', data); print(float(m.group(1)))"
    width=as.numeric(system2('python',c('-c',shQuote(read_width)),stdout=TRUE))
    # R's PDF device truncates MediaBox dimensions to whole points.
    if (panel_count == 0) {
        baseline_width = width
    } else {
        # PDF rounds to whole points; additions retain their physical widths.
        stopifnot((width - baseline_width) >=
            (panel_count * 15 + 8) / 25.4 * 72 - 1)
    }
}
setwd(oldwd)
unlink(outdir,recursive=TRUE)
cat('Workflow PDF generation, absent input, and physical width passed.\n')

# Adding fixed-width columns must preserve every ordinary column's proportions,
# including tiplabel's special size rule and explicit relative-width overrides.
base=list(tree=g$tree,tiplabel=ggplot(data.frame(label=c('long_gene_a','long_gene_b'))),domain=ggplot())
for (override in c('', 'tree,2,domain,1')) {
    before=get_rel_widths(base,override)
    for (count in c(1,3)) {
        expanded=base
        for(i in seq_len(count)) {
            panel=ggplot();attr(panel,'treevis_width_mm')=15
            expanded[[paste0('triangle',i)]]=panel
        }
        after=get_rel_widths(expanded,override)[names(before)]
        stopifnot(isTRUE(all.equal(before/sum(before),after/sum(after))))
    }
}
f=tempfile()
write.table(data.frame(node_name=c('a','b'),offset=1,group_id=c('001','1')),f,
            sep='\t',quote=FALSE,row.names=FALSE)
p=add_synteny_similarity_column(g,args,f)$synteny_similarity
# All scored pairs are disjoint; numeric-looking group IDs stay distinct.
stopifnot(all(p$data$similarity[!is.na(p$data$similarity)]==0))
unlink(f)
cat('Fixed-width isolation and exact synteny group IDs passed.\n')
