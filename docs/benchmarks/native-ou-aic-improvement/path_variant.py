"""Prototype covariance-updated OU-effect lasso path, no final penalized fits."""
import math
import numpy as np
from nwkit.shift_native_fit import fit_native_layout
from nwkit.shift_native_model import ShiftLayout, covariance_geometry
from nwkit.shift_native_screen import descendant_design, _proximal_step, _kkt_residual
from nwkit.gaussian_whitening import TreeWhitening
from nwkit.shift_native_search import NativeLayoutEvaluator,observable_layout

def matrices(data,fit):
 design,branches=descendant_design(data.tree)
 age={}
 remaining=np.zeros(len(data.tree.branch_ids))
 for i in data.tree.compiled.postorder:
  if i:
   parent=data.tree.compiled.parents[i];remaining[parent]=max(remaining[parent],data.tree.times[i]+remaining[i])
 for i,b in enumerate(data.tree.branch_ids):
  if i:age[b]=remaining[data.tree.compiled.parents[i]]
 ages=np.array([age[b] for b in branches]);xs=[];ys=[]
 for j,f in enumerate(fit['fits']):
  alpha=f.alpha_height
  weights=ages if alpha==0 else np.ones(len(ages)) if math.isinf(alpha) else -np.expm1(-alpha*ages)/-np.expm1(-alpha)
  mask=np.isfinite(data.values[:,j]);obs=tuple(i for i,m in zip(data.tree.compiled.leaf_indices,mask) if m)
  slopes,innovations,rv=covariance_geometry(data.tree,alpha,f.process_variance,f.root_model)
  factor=TreeWhitening.build(data.tree.compiled,obs,slopes,innovations,data.variances[mask,j]+f.measurement_variance,root_variance=rv)
  white=factor.apply(np.column_stack((np.ones(sum(mask)),data.values[mask,j],design[mask]*weights)))
  intercept=white[:,0]/np.linalg.norm(white[:,0]);projected=white[:,1:]-intercept[:,None]*(intercept@white[:,1:])[None,:]
  xs.append(projected[:,1:]);ys.append(projected[:,0])
 return xs,ys,branches

def path(data,fit,max_shifts,iterations=300,paths=80):
 xs,ys,branches=matrices(data,fit)
 grad=np.column_stack([x.T@y for x,y in zip(xs,ys)])
 maximum=float(np.max(np.linalg.norm(grad,axis=1)))
 coefficients=np.zeros((len(branches),len(ys)));step=1/max(1,max(np.sum(x*x) for x in xs));layouts=[];records=[]
 for fraction in np.geomspace(.999,.005,paths):
  strength=maximum*fraction
  for iteration in range(iterations):
   coefficients,gradient,step=_proximal_step(xs,ys,coefficients,strength,step)
   residual=_kkt_residual(coefficients,gradient,strength)
   if residual<=1e-5*max(1,maximum):break
   step*=1.05
  magnitudes=np.linalg.norm(coefficients,axis=1)
  active=np.flatnonzero(magnitudes>1e-8*max(1,float(np.max(magnitudes))))
  records.append({'fraction':float(fraction),'active':len(active),'kkt':residual,'iterations':iteration+1})
  if len(active)>max_shifts:
   # Preserve configurations across sizes when a coarse path jumps over a size.
   active=sorted(active,key=lambda j:(-magnitudes[j],branches[j]))[:max_shifts]
  if not len(active):continue
  try:layout=ShiftLayout.build(data.tree,[branches[j] for j in active])
  except ValueError:continue
  if observable_layout(data,layout) and layout not in layouts:layouts.append(layout)
  if records[-1]['active']>2*max_shifts:break
 return layouts,records

def search(data,options,variant):
 evaluator=NativeLayoutEvaluator(data,{},'AIC');null=ShiftLayout.build(data.tree);evaluator.evaluate(null)
 fit=fit_native_layout(data,null,alpha_height=0.0)
 metadata=[]
 for round in range(2):
  layouts,records=path(data,fit,options.max_shifts,iterations=options.lasso_iterations)
  # Preserve diversity over shift counts within the same bounded refit budget.
  remaining=options.refit_budget-len(evaluator.records)
  limit=remaining//(2-round)
  fresh=[x for x in layouts if x not in evaluator.scores]
  if len(fresh)>limit:
   indices=np.linspace(0,len(fresh)-1,limit).round().astype(int);fresh=[fresh[i] for i in indices]
  for layout in fresh:evaluator.evaluate(layout)
  metadata.append({'path':records,'candidate_count':len(layouts),'refits':len(fresh)})
  fit=evaluator.best_information
 return evaluator.finish({'variant':variant,'paths':metadata,'refits':len(evaluator.records)})
