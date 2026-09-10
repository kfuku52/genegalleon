"""Independent branch-recursion scope check with unknown shift count and tree shape."""
import csv,json,math,hashlib
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent/'extension'
for seed in range(29101,29104):
 for shape in ['balanced','random']:
  rng=np.random.default_rng(seed);nodes=[]
  def balanced(start,count):
   node={'name':f't{start}' if count==1 else '', 'children':[]}
   if count>1:node['children']=[balanced(start,count//2),balanced(start+count//2,count-count//2)]
   node['height']=max((c['height'] for c in node['children']),default=-1)+1
   node['tips']=sorted(sum((c['tips'] for c in node['children']),[])) if count>1 else [node['name']]
   nodes.append(node);return node
  if shape=='balanced':tree=balanced(0,100)
  else:
   active=[{'name':f't{i}','children':[],'height':0,'tips':[f't{i}']} for i in range(100)];nodes.extend(active);height=0
   while len(active)>1:
    chosen=sorted(rng.choice(len(active),2,replace=False),reverse=True);children=[active.pop(int(i)) for i in chosen]
    height+=float(rng.exponential(1/len(active))) if active else float(rng.exponential())
    node={'name':'','children':children,'height':height,'tips':sorted(sum((c['tips'] for c in children),[]))};nodes.append(node);active.append(node)
   tree=active[0]
  H=tree['height']
  candidates=[n for n in nodes if 3<=len(n['tips'])<=10];order=rng.permutation(len(candidates));chosen=[];occupied=set()
  for i in order:
   n=candidates[int(i)]
   if occupied.isdisjoint(n['tips']):chosen.append(n);occupied.update(n['tips'])
   if len(chosen)==5:break
  assert len(chosen)==5
  multipliers=rng.uniform(.8,1.2,5)*rng.choice([-1,1],5);noise=rng.normal(size=len(nodes));indices={id(n):i for i,n in enumerate(nodes)}
  def newick(n):return n['name'] if not n['children'] else '('+','.join(newick(c)+':'+format((n['height']-c['height'])/H,'.17g') for c in n['children'])+')'
  for strength in [2,6]:
   d=ROOT/'data'/f'{shape}-effect{strength}-seed{seed}';d.mkdir(parents=True,exist_ok=False);optima={id(n):strength*m for n,m in zip(chosen,multipliers)};rows=[];alpha=3.;variance=1/-math.expm1(-2*alpha)
   def simulate(n,state=0.,mean=0.,theta=0.,parent=None):
    theta=optima.get(id(n),theta)
    if parent is not None:
     length=(parent-n['height'])/H;decay=math.exp(-alpha*length);state=decay*state+math.sqrt(variance*(-math.expm1(-2*alpha*length)))*noise[indices[id(n)]];mean=decay*mean+(1-decay)*theta
    for c in n['children']:simulate(c,state,mean,theta,n['height'])
    if not n['children']:rows.append([n['name'],state+mean,mean])
   simulate(tree);(d/'tree.nwk').write_text(newick(tree)+';\n')
   with (d/'traits.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['taxon','x']);w.writerows(r[:2] for r in rows)
   truth={'seed':seed,'effect':strength,'shape':shape,'alpha_height':alpha,'process_tip_variance':1,'shift_clades':[n['tips'] for n in chosen],'tip_mean':{r[0]:r[2] for r in rows},'sha256':{p:hashlib.sha256((d/p).read_bytes()).hexdigest() for p in ['tree.nwk','traits.tsv']}}
   (d/'truth.json').write_text(json.dumps(truth,indent=2)+'\n')
