library(kfl1ou)
library(ape)
library(jsonlite)
args <- commandArgs(trailingOnly=TRUE)
directory <- args[1]; output <- args[2]; criterion <- args[3]
started <- proc.time()[['elapsed']]
result <- tryCatch({
  tree <- reorder.phylo(read.tree(file.path(directory,'tree.nwk')), 'postorder')
  tab <- read.delim(file.path(directory,'traits.tsv'),row.names=1,check.names=FALSE)
  Y <- as.matrix(tab[tree$tip.label,,drop=FALSE])
  fit <- estimate_shift_configuration(tree,Y,max.nShifts=10,criterion=criterion,root.model='OUfixedRoot',nCores=1,rescale=FALSE,measurement_error=FALSE,quietly=TRUE)
  descendants <- function(node) {
    if(node <= length(fit$tree$tip.label)) return(fit$tree$tip.label[node])
    sort(unlist(lapply(fit$tree$edge[fit$tree$edge[,1]==node,2],descendants)))
  }
  selected <- lapply(as.integer(fit$shift.configuration),function(e) descendants(fit$tree$edge[e,2]))
  predictions <- as.numeric(fitted(fit))
  list(status='complete',selection=criterion,score=fit$score,selected=list(shift_clades=selected,predicted=as.list(setNames(predictions,fit$tree$tip.label)),log_likelihood=as.numeric(logLik(fit))),alpha=as.numeric(fit$alpha),package_version=as.character(packageVersion('kfl1ou')))
},error=function(e) list(status='failed',error=conditionMessage(e)))
result$elapsed_seconds <- proc.time()[['elapsed']]-started
write_json(result,output,auto_unbox=TRUE,pretty=TRUE,digits=16)
if(result$status=='failed') quit(status=1)
