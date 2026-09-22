import json, pathlib, sqlite3, subprocess, sys, tempfile
sys.path.insert(0, '/repo/workflow/support')
import generate_orthogroup_database as m
import pandas as pd
result={}
frame=pd.DataFrame({'OCNany2spe':[0.1,0.9]})
for value in ['OCNany2spe,0.8','OCNany2spe,typo','OCNany2sp,0.8','OCNany2spe,nan']:
    result[value]={'parsed':m.parse_cutoff_stat(value),'retained':m.apply_cutoff(frame,value)['OCNany2spe'].tolist()}
with tempfile.TemporaryDirectory() as td:
    root=pathlib.Path(td)
    tree=root/'stat_tree';tree.mkdir()
    branch=root/'stat_branch';branch.mkdir()
    (tree/'OG0001_stat.tree.tsv').write_text('num_branch\tnum_spe\tnum_dup\tnum_sp\n3\t1\t0\t2\n')
    (branch/'OG0001_stat.branch.tsv').write_text('branch_id\tnode_name\tnum_sp\tso_event\n0\tn0\t2\tS\n')
    db=root/'out.db'
    command=[sys.executable,'/repo/workflow/support/generate_orthogroup_database.py','--dbpath',str(db),'--dir_stat_tree',str(tree),'--dir_stat_branch',str(branch)]
    result['default_runs']=[]
    for i in range(2):
        proc=subprocess.run(command,cwd=root,text=True,capture_output=True)
        with sqlite3.connect(db) as c:
            result['default_runs'].append({'exit':proc.returncode,'tree_rows':c.execute('select count(*) from tree').fetchone()[0],'branch_rows':c.execute('select count(*) from branch').fetchone()[0]})
print(json.dumps(result,indent=2))
