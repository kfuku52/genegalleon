import json, os, resource, shutil, subprocess, time
from pathlib import Path
replicate=os.environ['REPLICATE']
directory=Path('/audit/measure-'+replicate)
directory.mkdir(exist_ok=True)
fixture=Path('/review/workflow/tests/data/csubst_scan_inference')
shutil.copyfile(fixture/'input.fa',directory/'full.fa')
shutil.copyfile(fixture/'tree.nwk',directory/'tree.nwk')
command=['csubst','inspect','--full_cds_alignment_file','full.fa','--rooted_tree_file','tree.nwk',
         '--nonsyn_recode','3di20','--iqtree_model','GY+FQ','--threads','2','--sa_cache','no',
         '--sa_state_cache','no','--sa_no_download','yes','--outdir','inspect']
env=dict(os.environ,CSUBST_CACHE_DIR='/audit/downloads-parallel',HF_HUB_OFFLINE='1',TRANSFORMERS_OFFLINE='1')
start=time.perf_counter()
with open(directory/'run.log','w') as log:
    proc=subprocess.run(command,cwd=directory,env=env,stdout=log,stderr=subprocess.STDOUT)
usage=resource.getrusage(resource.RUSAGE_CHILDREN)
result=dict(replicate=replicate,command=command,returncode=proc.returncode,wall_seconds=time.perf_counter()-start,
            peak_child_rss_kib=usage.ru_maxrss,user_seconds=usage.ru_utime,system_seconds=usage.ru_stime)
(directory/'measurement.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result),flush=True)
raise SystemExit(proc.returncode)
