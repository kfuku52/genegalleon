"""CAFE summary aggregation must preserve native branch counts and reference geometry."""
import json
import subprocess
from pathlib import Path

SUPPORT = Path(__file__).resolve().parents[1] / "support"


def test_summary_counts_equal_independent_family_loop_and_validate_identity():
    code = r'''

suppressPackageStartupMessages({library(ggtree);library(ape)})
e<-new.env()
source(file.path(SUPPORT,'tree_label_formatting.R'),local=e)
source(file.path(SUPPORT,'cafe_plot_labels.R'),local=e)
for(x in parse(file.path(SUPPORT,'cafe_plot_summary.r')))if(is.call(x)&&as.character(x[[1]])%in%c('<-','=')&&length(x)>=3&&is.call(x[[3]])&&identical(x[[3]][[1]],as.name('function')))eval(x,e)
set.seed(47);ref<-ape::rtree(8);ref$tip.label<-paste0('Species_',1:8,'<',1:8,'>_2');ref$node.label<-paste0('<',9:15,'>_2')
keys<-e$norm_key(c(ref$tip.label,ref$node.label));ids<-paste0('OG',1:100)
changes<-as.data.frame(matrix(sample(c(-2:2,NA),1500,replace=TRUE),nrow=100));names(changes)<-keys;changes$FamilyID<-ids
# Put both the input rows and branch columns in a different order.
changes<-changes[sample(100),c('FamilyID',sample(keys))]
trees<-lapply(1:100,function(i){t<-ref;t$tip.label[i%%8+1]<-paste0(t$tip.label[i%%8+1],'*');if(i%%2)t$node.label[i%%7+1]<-paste0(t$node.label[i%%7+1],'*');t});names(trees)<-ids
out<-suppressWarnings(e$accumulate_branch_counts(trees,changes,ids,ref))
fort<-suppressWarnings(ggtree::fortify(ref));stopifnot(identical(out$key,e$norm_key(fort$label)),identical(out$x,fort$x),identical(out$y,fort$y),identical(out$isTip,fort$isTip))
# Reference uses the original per-family fortifier and an independent indexed sum.
expected<-matrix(0L,nrow=nrow(out),ncol=4)
for(id in ids){ann<-suppressWarnings(e$create_annotation_df(id,trees[[id]],changes));positions<-match(out$key,ann$key);v<-ann$change[positions];flag<-ann$has_star[positions];expected<-expected+cbind(as.integer(!is.na(v)&v>0),as.integer(!is.na(v)&v<0),as.integer(!is.na(v)&v>0&flag),as.integer(!is.na(v)&v<0&flag))}
stopifnot(identical(unname(as.matrix(out[c('inc','dec','inc_sig','dec_sig')])),unname(expected)))
empty<-suppressWarnings(e$accumulate_branch_counts(trees,changes,character(),ref));stopifnot(all(as.matrix(empty[c('inc','dec','inc_sig','dec_sig')])==0))
bad<-trees;bad[[1]]$tip.label[1]<-'Other_species<99>_2'
expect_error<-function(expr)stopifnot(inherits(tryCatch(force(expr),error=identity),'error'))
expect_error(suppressWarnings(e$accumulate_branch_counts(bad,changes,ids,ref)))
expect_error(suppressWarnings(e$accumulate_branch_counts(trees,rbind(changes,changes[1,]),ids,ref)))
expect_error(suppressWarnings(e$accumulate_branch_counts(trees,changes[,setdiff(names(changes),keys[1])],ids,ref)))
cat('All summary counts and geometry agree; malformed identities are rejected\n')
'''
    result = subprocess.run(["Rscript", "-e", "SUPPORT<-" + json.dumps(str(SUPPORT)) + "\n" + code],
                            text=True, capture_output=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr
