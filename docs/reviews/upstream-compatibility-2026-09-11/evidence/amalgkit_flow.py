import pathlib, subprocess, random, gzip, re, json, csv, sys
p=pathlib.Path('/tmp/gg-amalgkit-flow'); p.mkdir(exist_ok=True)
raw=p/'raw'/'Test_species';raw.mkdir(parents=True,exist_ok=True)
rng=random.Random(71);seq=''.join(rng.choice('ACGT') for _ in range(1000))
with gzip.open(raw/'R1.fastq.gz','wt') as f:
 for i in range(2000):
  start=rng.randrange(850); read=seq[start:start+100]; f.write(f'@r{i}\n{read}\n+\n'+('I'*100)+'\n')
fasta=p/'fasta';fasta.mkdir(exist_ok=True);(fasta/'Test_species.fa').write_text('>tx1\n'+seq+'\n')
log=open('/audit/amalgkit-flow.log','w')
def run(cmd,**kwargs):
 print('RUN',cmd,file=log,flush=True)
 x=subprocess.run(cmd,stdout=log,stderr=log,cwd=p,timeout=120,**kwargs)
 print('EXIT',x.returncode,file=log,flush=True)
 if x.returncode: raise SystemExit(x.returncode)
from amalgkit.metadata_utils import Metadata
import pandas
meta=p/'metadata_private_fastq.tsv'
Metadata.from_DataFrame(pandas.DataFrame([dict(run='R1', scientific_name='Test species', private_file='yes', read1_path=str(raw/'R1.fastq.gz'), read2_path='', lib_layout='single', total_spots=2000, total_bases=200000, spot_length=100, exclusion='no', is_sampled='yes', instrument='Illumina', taxid='1')])).df.to_csv(meta, sep='\t', index=False)
run(['amalgkit','getfastq','--out_dir',str(p),'--metadata',str(meta),'--threads','1','--rrna_filter','no','--contam_filter','no','--read_name','trinity','--remove_sra','yes','--remove_tmp','yes','--dump_print','yes'])
core=pathlib.Path('/review/workflow/core/gg_transcriptome_generation_core.sh').read_text()
b=core[core.index('bind_amalgkit_getfastq_completion_manifest()'):]
embedded=re.search(r"<<'PY'\n(.*?)\nPY",b,re.S).group(1)
run(['python','-',str(next((p/'getfastq').glob('*completion*.json'))),str(meta)],input=embedded,text=True)
run(['amalgkit','quant','--out_dir',str(p),'--metadata',str(meta),'--threads','1','--clean_fastq','no','--fasta_dir',str(fasta),'--build_index','yes','--quant_backend','auto','--oarfish_seq_tech','auto'])
run(['python','/review/workflow/support/validate_transcriptome_quant_outputs.py','--metadata',str(meta),'--quant-root',str(p/'quant')])
run(['amalgkit','merge','--out_dir',str(p),'--metadata',str(meta)])
print('AMALGKIT flow succeeded',file=log)
