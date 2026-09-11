import ast, importlib, json, pathlib, re, subprocess
root=pathlib.Path('/review/workflow')
results=[]
for p in list((root/'core').glob('*.sh'))+list((root/'support').rglob('*.sh')):
 lines=p.read_text().splitlines()
 for i,line in enumerate(lines):
  m=re.search(r'(?:^|\s|\|)(csubst|cdskit|nwkit|amalgkit)\s+([a-z][a-z0-9-]*)\b',line)
  if not m or line.lstrip().startswith('#') or any(x in line[:m.start()] for x in ['echo ', 'task=', '"', "'"]): continue
  tool,sub=m.groups(); block=line[m.end():]; j=i
  while j+1<len(lines) and (lines[j].rstrip().endswith('\\') or (not block.strip() and lines[j+1].lstrip().startswith('--'))):
   j+=1
   if re.search(r'\|\s*(nwkit|cdskit|csubst|amalgkit)',lines[j]): break
   block+='\n'+lines[j]
  block=block.split('|')[0]
  opts=sorted(set(re.findall(r'(?<![\w-])--[a-zA-Z][a-zA-Z0-9_-]*',block)))
  results.append({'file':str(p.relative_to(root)), 'line':i+1,'tool':tool,'command':sub,'options':opts})
help_cache={}
for r in results:
 k=(r['tool'],r['command'])
 if k not in help_cache:
  q=subprocess.run([*k,'--help-advanced' if k[0]=='csubst' else '--help'],text=True,capture_output=True)
  help_cache[k]=(q.returncode,set(re.findall(r'--[a-zA-Z][a-zA-Z0-9_-]*',q.stdout+q.stderr)))
 rc,opts=help_cache[k]
 r['help_exit']=rc; r['unknown_options']=sorted(set(r['options'])-opts)
imports=[]
for p in (root/'support').glob('*.py'):
 for n in ast.walk(ast.parse(p.read_text())):
  if isinstance(n,ast.ImportFrom) and n.module and n.module.split('.')[0] in {'nwkit','csubst','amalgkit','cdskit'}:
   try:
    mod=importlib.import_module(n.module)
    missing=[]
    for a in n.names:
     if not hasattr(mod,a.name):
      try: importlib.import_module(n.module+'.'+a.name)
      except ModuleNotFoundError: missing.append(a.name)
    imports.append({'file':p.name,'line':n.lineno,'module':n.module,'missing':missing})
   except Exception as e: imports.append({'file':p.name,'line':n.lineno,'module':n.module,'error':str(e)})
pathlib.Path('/audit/cli-import-check.json').write_text(json.dumps({'shell_calls':results,'imports':imports},indent=2))
print('shell_calls',len(results),'unique_commands',len(help_cache),'imports',len(imports))
print(json.dumps([r for r in results if r['help_exit'] or r['unknown_options']],indent=2))
print(json.dumps([r for r in imports if r.get('missing') or r.get('error')],indent=2))
