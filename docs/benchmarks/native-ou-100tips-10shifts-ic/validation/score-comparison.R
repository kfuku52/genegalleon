library(kfl1ou)
library(ape)
library(jsonlite)
tree <- reorder.phylo(read.tree(text="(((a:1,b:1):1,(c:1,d:1):1):1,((e:1,f:1):1,(g:1,h:1):1):1);"), "postorder")
y <- c(1.5549029785536523,1.3137883823887941,-.21225628554971665,-.24002282307176648,1.7035058708492061,1.8521601652123756,.7166573577624096,.2547434554367455)
Y <- matrix(y,ncol=1,dimnames=list(letters[1:8],"x"))[tree$tip.label,,drop=FALSE]
shifts <- match(match(c("a","e"), tree$tip.label),tree$edge[,2])
rows <- list()
for(root in c("OUfixedRoot","OUrandomRoot")) for(shared in c(FALSE,TRUE)) for(criterion in c("AIC","BIC","pBIC")) {
 fit <- fit_OU(tree,Y,shifts,criterion=criterion,cr.regimes=if(shared) list(0L,shifts) else NULL,root.model=root,alpha.lower=.4,alpha.upper=.4,compute.hessian=FALSE)
 rows[[length(rows)+1]] <- list(root=root,shared=shared,criterion=criterion,score=fit$score,log_likelihood=sum(fit$logLik))
}
write_json(rows,commandArgs(TRUE)[1],auto_unbox=TRUE,pretty=TRUE,digits=16)
