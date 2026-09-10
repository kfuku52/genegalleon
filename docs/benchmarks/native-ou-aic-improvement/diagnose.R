library(kfl1ou); library(ape); library(jsonlite)
root <- '/bench'; old <- '/old'; rows <- list()
for (directory in list.dirs(file.path(old,'data'), recursive=FALSE)) {
 tree <- reorder.phylo(read.tree(file.path(directory,'tree.nwk')),'postorder')
 tab <- read.delim(file.path(directory,'traits.tsv'), row.names=1)
 Y <- as.matrix(tab[tree$tip.label,,drop=FALSE])
 descendants <- function(node) {if(node<=length(tree$tip.label)) return(tree$tip.label[node]); sort(unlist(lapply(tree$edge[tree$edge[,1]==node,2],descendants)))}
 keys <- vapply(tree$edge[,2],function(n) paste(descendants(n),collapse=','),'')
 for (origin in c('nwkit-AIC','kfl1ou-AIC','truth')) {
  saved <- fromJSON(file.path(directory,paste0(origin,'.json')), simplifyVector=FALSE)
  clades <- if(origin=='truth') saved$shift_clades else saved$selected$shift_clades
  shifts <- vapply(clades,function(x) match(paste(sort(unlist(x)),collapse=','),keys),1L)
  fit <- fit_OU(tree,Y,shifts,criterion='AIC',root.model='OUfixedRoot',compute.hessian=FALSE)
  rows[[length(rows)+1]] <- list(dataset=basename(directory),origin=origin,alpha=as.numeric(fit$alpha),score=fit$score,log_likelihood=as.numeric(logLik(fit)),predicted=as.list(setNames(as.numeric(fitted(fit)),tree$tip.label)))
 }
}
write_json(rows,file.path(root,'crossfit-r.json'),auto_unbox=TRUE,pretty=TRUE,digits=16)
