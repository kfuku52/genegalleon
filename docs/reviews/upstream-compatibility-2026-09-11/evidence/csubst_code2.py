import pathlib,subprocess,sys,os
p=pathlib.Path('/tmp/gg-csubst-code2');p.mkdir(exist_ok=True)
fixture=pathlib.Path('/review/workflow/tests/data/csubst_scan_inference')
lines=(fixture/'input.fa').read_text().splitlines();out=[]
for line in lines:
 if not line.startswith('>'):line=''.join('GGA' if line[i:i+3] in {'AGA','AGG','TAA','TAG'} else line[i:i+3] for i in range(0,len(line),3))
 out.append(line)
(p/'csubst.fasta').write_text('\n'.join(out)+'\n');(p/'csubst.nwk').write_bytes((fixture/'tree.nwk').read_bytes())
log=open('/audit/sites-code2.log','w')
subprocess.run(['iqtree','-s','csubst.fasta','-te','csubst.nwk','-m','GY+F+R4','-T','1','--seqtype','CODON2','--prefix','csubst','--ancestral','--rate','--seed','12345','--redo'],cwd=p,stdout=log,stderr=log,check=True,timeout=120)
sys.path.insert(0,'/review/workflow/support');import csubst_site_wrapper as w
cmd=w.build_csubst_sites_command(str(p),str(p),'1,2',1,'no',pdb='none')
for name,extra in [('current',[]),('explicit-code2',['--genetic_code','2'])]:
 d=p/name;d.mkdir(exist_ok=True)
 f=open('/audit/sites-code2-'+name+'.log','w')
 q=subprocess.run(cmd+extra,cwd=d,stdout=f,stderr=f,timeout=120)
 print(name,q.returncode,file=log,flush=True)
