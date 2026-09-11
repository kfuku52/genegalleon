import pathlib, subprocess, os, re
p=pathlib.Path('/tmp/gg-csubst-search');p.mkdir(exist_ok=True)
repo=pathlib.Path('/review');core=(repo/'workflow/core/gg_gene_evolution_core.sh').read_text()
entry=(repo/'workflow/gg_gene_evolution_entrypoint.sh').read_text()
fixture=repo/'workflow/tests/data/csubst_scan_inference'
(p/'csubst.fasta').write_bytes((fixture/'input.fa').read_bytes());(p/'csubst.nwk').write_bytes((fixture/'tree.nwk').read_bytes())
(p/'foreground.tsv').write_text('name\ttrait\n'+''.join(f'{x}\t{1 if x=="a" else 2 if x=="e" else 0}\n' for x in 'abcdefgh'))
log=open('/audit/search-flow.log','w')
subprocess.run(['iqtree','-s','csubst.fasta','-te','csubst.nwk','-m','ECMK07+F+R4','-T','AUTO','--threads-max','1','--seqtype','CODON1','--prefix','csubst','--ancestral','--rate','--seed','12345','--redo'],cwd=p,stdout=log,stderr=log,check=True,timeout=120)
start=core.index('  csubst search \\');end=core.index('\n\n',start)
config='\n'.join(x for x in entry.splitlines() if re.match('csubst_(max_arity|exhaustive_until|cutoff_stat|max_combination|fg_exclude_wg|fg_stem_only|nonsyn_recode)=',x))
shell='set -euo pipefail\n'+config+'\nforeground_params=(--foreground foreground.tsv --fg_format 2)\n'+core[start:end]
proc=subprocess.run(['bash','-c',shell],cwd=p,env=dict(os.environ,genetic_code='1',csubst_input_base='./csubst',codon_model='ECMK07+F+R4',GG_TASK_CPUS='1'),stdout=log,stderr=log,timeout=180)
print('search exit',proc.returncode,file=log,flush=True)
raise SystemExit(proc.returncode)
