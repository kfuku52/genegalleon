"""Actual spawned genomic replay must match the serial decision path exactly."""
import copy

import pytest
from Bio.Seq import Seq

from workflow.support import rescue_cached_validation as parallel
from workflow.support import rescue_gene_models as rescue
from workflow.support.gene_model_catalog import indexed_genome


def fixture(tmp_path, count):
    sequences = {'complete':'ATGAAATAA', 'stop':'ATGTGATAA', 'partial':'ATGAAA',
                 'minus':str(Seq('ATGAAATAA').reverse_complement()), 'splice':'ATGATAAAAGAAATAA'}
    path = tmp_path/'genome.fa'
    path.write_text(''.join('>'+name+'\n'+dna+'\n' for name,dna in sequences.items()))
    rows = []
    for i in range(count):
        seqid = list(sequences)[i % len(sequences)]
        row = dict(query='q'+str(i), seqid=seqid, strand='-' if seqid == 'minus' else '+',
            cds=[[0,len(sequences[seqid]),0]] if seqid != 'splice' else [[0,3,0],[10,16,0]],
            frameshift=i % 13 == 0, coverage=.6 if i % 7 == 0 else 1., identity=.8,
            search='genome_fallback', evidence={'target':'Target','donor':'Donor','query':'p','genome_only':True},
            unknown={'record':i,'nested':[i,{'preserve':True}]})
        rows.append(row)
    return path,rows


def evaluate(path, rows, cpus, code=1):
    arguments = ['plan', *[value for key in ('cds-dir','gff-dir','genome-dir','busco-dir','tree','output')
                          for value in ('--'+key, str(path.parent/key))]]
    defaults = rescue.parser().parse_args(arguments)
    params = {key:getattr(defaults,key) for key in rescue.PARAMETERS}
    proteins = {'Donor':{'p':'MK'}}
    audit = {}
    with indexed_genome(path) as genome:
        memo = rescue.RawValidationMemo(rescue.validate_model,genome,code,params,
                                        candidate_homology=rescue.homology_problems)
        def serial(model):
            return parallel.checked_prediction(model,memo,genome,code,params,'MK',engine=rescue)
        result = list(parallel.validate_cached(iter(rows),genome,code,params,proteins,cpus,serial,audit))
        assert genome.fetch('complete',0,9)=='ATGAAATAA'
    return result,audit


@pytest.mark.parametrize('code',[1,2])
def test_spawn_preserves_full_records_order_genetic_code_and_rejections(tmp_path,code):
    path,rows=fixture(tmp_path,577)
    original=copy.deepcopy(rows)
    serial,_=evaluate(path,copy.deepcopy(rows),1,code)
    spawned,audit=evaluate(path,rows,2,code)
    assert spawned==serial
    assert rows==original
    assert [r['query'] for r in spawned]==['q'+str(i) for i in range(577)]
    assert all(r['raw_prediction']['cds']==original[i]['cds'] for i,r in enumerate(spawned))
    assert audit['records']==577 and audit['workers']==2 and audit['max_pending_batches']<=4
    assert audit['worker_memo_counts']['calls']==577
    assert 'unanchored_genome_search' in spawned[1]['problems']
    spawned[0]['unknown']['nested'].append('changed')
    assert rows==original


def test_small_replay_keeps_in_process_callback_without_opening_genome_or_pool():
    rows=[{'q':i} for i in range(20)]
    seen=[]
    def serial(model):
        seen.append(model['q'])
        return {'same':model['q']}
    audit={}
    result=list(parallel.validate_cached(iter(rows),None,1,{},None,8,serial,audit))
    assert result==[{'same':i} for i in range(20)] and seen==list(range(20))
    assert audit['workers']==1 and audit['records']==20


def test_child_exception_is_not_suppressed_or_published_as_a_model(tmp_path):
    path,rows=fixture(tmp_path,513)
    rows[-1]['cds']=[[0,1000000,0]]
    for cpus in (1,2):
        with pytest.raises(ValueError,match='Predicted CDS outside genome'):
            evaluate(path,copy.deepcopy(rows),cpus)


@pytest.mark.parametrize('cpus',[0,-1,True,1.5])
def test_invalid_cpu_bound_is_rejected(cpus):
    with pytest.raises(ValueError,match='CPUs must be positive'):
        list(parallel.validate_cached(iter([]),None,1,{},None,cpus,lambda x:x,{}))
