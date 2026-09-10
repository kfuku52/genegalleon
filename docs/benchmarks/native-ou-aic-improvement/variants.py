"""Development-only search ablations; not public NWKIT options."""
import itertools
import math
from nwkit.shift_native_heuristic import _HeuristicSearch, added_layouts, neighboring_layouts
from nwkit.shift_native_quick import NativeQuickProfile
from nwkit.shift_native_screen import group_lasso_screen
from nwkit.shift_native_model import ShiftLayout

class ProfiledQuick(NativeQuickProfile):
    def _trait_profile(self, trait, alpha, fit, design):
        # Whiten at unit process scale, then profile scale analytically in score.
        # This experimental adapter is only for no observation errors/free scale.
        from dataclasses import replace
        return super()._trait_profile(trait,alpha,replace(fit,process_variance=1.0),design)

# Keep the orthogonal fit/rank algorithm identical, changing only the scale profile.
import inspect
source=inspect.getsource(NativeQuickProfile.score)
source=source.replace('    def score','def score',1)
source='\n'.join(line[4:] if line.startswith('    ') else line for line in source.splitlines()) if False else source
# The method body already has eight-space indentation; valid for a free function.
source=source.replace('total += constant - 0.5 * float(residual @ residual)', 'quadratic = float(residual @ residual)\n                if quadratic <= 0:\n                    total = -math.inf\n                    break\n                total += constant - 0.5 * observed_count * (1 + math.log(quadratic / observed_count))')
import numpy as np
exec(source)
ProfiledQuick.score=score

class AdaptiveSearch(_HeuristicSearch):
    def __init__(self,data,options,fit_arguments,criterion,variant):
        super().__init__(data,options,fit_arguments,criterion)
        self.variant=variant
        self.refined=set()
        self.updates=[]
        if 'scale' in variant:
            self.profile=ProfiledQuick(data,self.evaluator.best[0],self.pool)
    def update(self,fit):
        pool,metadata=group_lasso_screen(self.data,fit,pool_size=self.options.candidate_pool,iterations=self.options.lasso_iterations,memory_limit=self.options.memory_limit)
        # Preserve selected branches and half of the original pool for diversity.
        pool=list(dict.fromkeys([*fit['layout'].shifts,*pool,*self.pool[:self.options.candidate_pool//2]]))
        evaluations=self.profile.evaluations
        cls=ProfiledQuick if 'scale' in self.variant else NativeQuickProfile
        self.profile=cls(self.data,fit,pool)
        self.profile.evaluations=evaluations
        self.pool=pool; self.quick_scores={}
        self.updates.append({'from':list(fit['layout'].shifts),'pool':pool})
    def forward(self):
        if 'adaptive' not in self.variant and 'refine' not in self.variant:
            return super().forward()
        frontier=[self.null]
        for size in range(1,self.options.max_shifts+1):
            if 'adaptive' in self.variant and size in (4,7,10):
                fit=self.evaluator.best[2*(size-1)]
                self.update(fit)
                # Retain other frontier branches when the covariance is updated.
                needed=list(dict.fromkeys([*self.pool,*(b for x in frontier for b in x.shifts)]))
                evaluations=self.profile.evaluations
                cls=ProfiledQuick if 'scale' in self.variant else NativeQuickProfile
                self.pool=needed; self.profile=cls(self.data,fit,needed); self.profile.evaluations=evaluations
            ranked=self.rank(itertools.chain.from_iterable(added_layouts(self.data,x,self.pool,False) for x in frontier))
            chosen=self.refit(ranked,self.options.beam_width)
            if not chosen:break
            frontier=chosen[:self.options.beam_width]
    def refine(self):
        if 'refine' not in self.variant and 'adaptive' not in self.variant:
            return super().refine()
        while len(self.evaluator.records)<self.options.refit_budget:
            fits=sorted(self.evaluator.best.values(),key=lambda f:(f['information_criterion']['score'],f['layout'].shifts))
            fit=next((f for f in fits if f['layout'] not in self.refined),None)
            if fit is None:break
            layout=fit['layout'];self.refined.add(layout)
            if 'adaptive' in self.variant:self.update(fit)
            proposed=neighboring_layouts(self.data,layout,self.pool,False)
            if len(layout.shifts)<self.options.max_shifts:
                proposed=itertools.chain(proposed,added_layouts(self.data,layout,self.pool,False))
            ranked=self.rank(proposed)
            # Fast AIC approximation: location and mean penalties, covariance common.
            ranked.sort(key=lambda x:(-2*self.quick_scores[x]+4*len(x.shifts),x.shifts))
            # Evaluate directly so complexity representatives do not override IC rank.
            for candidate in ranked[:min(self.options.beam_width,self.options.refit_budget-len(self.evaluator.records))]:
                self.evaluator.evaluate(candidate)
            if self.profile.evaluations>=self.options.screening_budget:break

def search(data,options,variant):
    s=AdaptiveSearch(data,options,{},'AIC',variant)
    s.forward();s.refine()
    return s.evaluator.finish({'variant':variant,'pool_updates':s.updates,'quick_evaluations':s.profile.evaluations,'refits':len(s.evaluator.records),'initial_screen':s.screen_metadata})
