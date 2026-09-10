"""Independent-backend fixed-alpha score agreement before benchmark timing."""
import json
import subprocess
from pathlib import Path
import numpy as np
from nwkit.shift_backend_probe import R_PBIC_PROBE, collect_pbic_validation
from nwkit.shift_native_fit import NativeFitOptions, fit_native_layout
from nwkit.shift_native_ic import native_information_criterion
from nwkit.shift_native_model import ShiftData, ShiftLayout
from nwkit.util import read_tree

ROOT=Path(__file__).resolve().parent
validation=ROOT/'validation'
validation.mkdir(exist_ok=True)
probe=validation/'pbic-probe.R'
probe.write_text(R_PBIC_PROBE+'\nnwkit_pbic_preflight(directory="'+str(validation)+'")\n')
subprocess.run(['Rscript',str(probe)],check=True)
attestation=collect_pbic_validation(validation)
(validation/'pbic-attestation.json').write_text(json.dumps(attestation,indent=2)+'\n')
r=validation/'score-comparison.R'
r.write_text('''library(kfl1ou)
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
''')
output=validation/'kfl-fixed-alpha-scores.json'
subprocess.run(['Rscript',str(r),str(output)],check=True)
tree=read_tree('(((a:1,b:1):1,(c:1,d:1):1):1,((e:1,f:1):1,(g:1,h:1):1):1);','auto',True,quiet=True)
y=np.array([1.5549029785536523,1.3137883823887941,-.21225628554971665,-.24002282307176648,1.7035058708492061,1.8521601652123756,.7166573577624096,.2547434554367455])
data=ShiftData.build(tree,y[:,None],['x'])
nodes=list(tree.traverse('levelorder'))
shifts=[i for i,n in enumerate(nodes) if n.name in ['a','e']]
checks=[]
for row in json.loads(output.read_text()):
 layout=ShiftLayout.build(data.tree,shifts,[[0],shifts] if row['shared'] else None)
 fit=fit_native_layout(data,layout,alpha_height=1.2,options=NativeFitOptions(root_model=row['root']))
 score=native_information_criterion(data,fit,row['criterion'])['score']
 errors={'score_error':abs(score-row['score']),'likelihood_error':abs(fit['log_likelihood']-row['log_likelihood'])}
 checks.append({**row,'native_score':score,**errors})
 assert max(errors.values())<2e-6, checks[-1]
(validation/'native-kfl-score-agreement.json').write_text(json.dumps({'status':'passed','tolerance':2e-6,'checks':checks},indent=2)+'\n')
print('pBIC probe:',len(attestation['checks']),'checks; native/kfl score agreement:',len(checks),'checks')
