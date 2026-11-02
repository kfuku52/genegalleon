# Synthetic 24-species/1500-family CAFE aggregation benchmark, three repetitions.

suppressPackageStartupMessages({library(ggtree);library(ape)})
args<-commandArgs(TRUE);if(length(args)!=2)stop('Usage: Rscript workflow/tests/benchmark_cafe_summary.R source.r counts.tsv');e<-new.env()
source('workflow/support/tree_label_formatting.R',local=e);source('workflow/support/cafe_plot_labels.R',local=e)
for(x in parse(args[1])) if(is.call(x)&&as.character(x[[1]])%in%c('<-','=')&&length(x)>=3&&is.call(x[[3]])&&identical(x[[3]][[1]],as.name('function')))eval(x,e)
set.seed(904);f<-1500L;tree<-ape::rtree(24);tree$tip.label<-paste0('Species_',seq_len(24),'<',seq_len(24),'>_3')
tree$node.label<-paste0('<',25:47,'>_3');keys<-e$norm_key(c(tree$tip.label,tree$node.label))
changes<-as.data.frame(matrix(sample(-3:3,f*47,replace=TRUE),nrow=f));names(changes)<-keys;changes$FamilyID<-paste0('OG',seq_len(f))
trees<-lapply(seq_len(f),function(i){t<-tree;t$tip.label[i%%24+1]<-paste0(t$tip.label[i%%24+1],'*');t$node.label[i%%23+1]<-paste0(t$node.label[i%%23+1],'*');t});names(trees)<-changes$FamilyID
invisible(suppressWarnings(e$accumulate_branch_counts(trees[1:30],changes[1:30,],changes$FamilyID[1:30],tree)))
times<-numeric(3)
for(i in seq_len(3))times[i]<-system.time(out<-suppressWarnings(e$accumulate_branch_counts(trees,changes,changes$FamilyID,tree)))[['elapsed']]
out<-out[order(out$key),];rownames(out)<-NULL
write.table(out,args[2],sep='\t',quote=FALSE,row.names=FALSE)
cat('synthetic families:',f,'species:',24,'repetitions:',3,'elapsed seconds:',paste(times,collapse=','),'median:',median(times),'\n')

system2('sha256sum', args[2])
