"""Prove unchanged coding QC before extending completed native evidence.

Formatting, search-cache dispatch and UTR interpretation are not reused as
coding decisions. Their archived bytes remain bound, and changed native parsers
must reconstruct identical coding paths for the actual frozen prediction GFF.
"""
import ast
import copy
import hashlib
import json
import subprocess
import sys
from pathlib import Path


def _r():
    try:
        import rescue_gene_models
    except ImportError:
        from . import rescue_gene_models
    return rescue_gene_models


def _replay():
    try:
        import rescue_admission_replay
    except ImportError:
        from . import rescue_admission_replay
    return rescue_admission_replay


def _definitions(path, names):
    return _replay()._definitions(path, names)


class _CallerNames(ast.NodeTransformer):
    """Normalise only the extraction's explicit context arguments/module alias."""
    protein = ast.dump(ast.parse('proteins[evidence["donor"]][evidence["query"]]',mode='eval').body,
                       include_attributes=False)
    code = ast.dump(ast.parse('source["genetic_code"]',mode='eval').body,include_attributes=False)

    def visit_Name(self,node):
        return ast.copy_location(ast.Name(id={'raw_prediction':'raw','raw_validation':'memo'}.get(node.id,node.id),
                                         ctx=node.ctx),node)

    def visit_Attribute(self,node):
        if (isinstance(node.value,ast.Name) and node.value.id=='r'
                and node.attr in {'RAW_FIELDS','complete_terminals','revalidate_boundaries','validate_model','check_interval'}):
            return ast.copy_location(ast.Name(id=node.attr,ctx=node.ctx),node)
        return self.generic_visit(node)

    def visit_Subscript(self,node):
        signature=ast.dump(node,include_attributes=False)
        if signature in {self.protein,self.code}:
            return ast.copy_location(ast.Name(id='protein' if signature==self.protein else 'code',ctx=ast.Load()),node)
        return self.generic_visit(node)


def caller_contract(support):
    """Require the same complete checking recipe across the worker extraction."""
    support=Path(support)
    tree=ast.parse((support/'rescue_gene_models.py').read_text())
    nodes=[n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='checked_prediction']
    if len(nodes)!=1:
        raise ValueError('Coding QC requires one checked-prediction caller')
    function=nodes[0]
    if any(isinstance(n,ast.Attribute) and n.attr=='checked_prediction' for n in ast.walk(function)):
        expected=ast.parse('''def checked_prediction(model):
    evidence = model["evidence"]
    return rescue_cached_validation.checked_prediction(model, raw_validation, genome,
        source['genetic_code'], params, proteins[evidence['donor']][evidence['query']],
        engine=sys.modules[__name__])
''').body[0]
        if ast.dump(function,include_attributes=False)!=ast.dump(expected,include_attributes=False):
            raise ValueError('Coding QC worker context changed')
        helper=support/'rescue_cached_validation.py'
        functions=[n for n in ast.parse(helper.read_text()).body
                   if isinstance(n,ast.FunctionDef) and n.name=='checked_prediction']
        if len(functions)!=1:
            raise ValueError('Coding QC worker helper is unavailable')
        function=functions[0]
    ignored={ast.dump(ast.parse(text).body[0],include_attributes=False) for text in
             ('evidence = model["evidence"]','r = engine or _r()')}
    body=[_CallerNames().visit(copy.deepcopy(n)) for n in function.body
          if ast.dump(n,include_attributes=False) not in ignored]
    return ast.dump(ast.Module(body=body,type_ignores=[]),include_attributes=False)


def prove_contract(previous_plan,producer_source):
    r,replay=_r(),_replay()
    old=Path(producer_source).resolve(strict=True)/'workflow/support'
    current=Path(r.__file__).resolve().parent
    request=previous_plan['request']
    if not request.get('oriongeno') or request.get('annevo') or request.get('gemoma_jar'):
        raise ValueError('Validated native extension requires an Orion-only non-GeMoMa producer')
    before,after=request['tools'],r.identities()
    # These modules serve input provenance, cache replay or noncoding UTRs, not
    # the reused genomic/terminal/boundary/donor decisions. Actual source inputs
    # stay identical; a changed native coding parser is checked per target.
    separate={'rescue_prediction_cache.py','transcript_utrs.py','genome_sequence_identity.py'}
    strict={name:field for name,field in replay.MODULES.items() if name not in separate}
    functions={'rescue_gene_models.py':('validate_model','check_interval','iter_miniprot'),
               'annevo_gene_models.py':replay.FUNCTIONS['annevo_gene_models.py'],
               'rescue_coding_paths.py':replay.FUNCTIONS['rescue_coding_paths.py'],
               'input_generation_array_state.py':('digest','digest_paths','atomic_json')}
    bindings={**replay.MODULES,'rescue_gene_models.py':'implementation',
              'annevo_gene_models.py':'annevo_candidate_implementation',
              'rescue_coding_paths.py':'resolve_coding_paths_implementation',
              'input_generation_array_state.py':'state_implementation',
              'oriongeno_runtime.py':'oriongeno_candidate_implementation'}
    if 'cached_validation_implementation' in before:
        bindings['rescue_cached_validation.py']='cached_validation_implementation'
        if before['cached_validation_implementation']!=after['cached_validation_implementation']:
            raise ValueError('The producer cached-validation dispatcher changed')
    if request.get('oriongeno_genomewide'):
        bindings['rescue_native_genomewide.py']='native_genomewide_implementation'
    proof=r.FreshDigestBatch()
    paths=[*[old/name for name in bindings],*[current/name for name in strict.keys()|functions.keys()],
           current/'rescue_cached_validation.py',current/'transcript_utrs.py',current/'rescue_native_genomewide.py']
    observed=proof.read(paths)
    old_files={}
    for name,field in bindings.items():
        value=observed[str(old/name)]
        if before.get(field)!=value:
            raise ValueError('Native extension producer source differs from its plan: '+name)
        old_files[str(old/name)]=value
    for name in strict:
        if observed[str(old/name)]!=observed[str(current/name)]:
            raise ValueError('Native extension coding evidence kernel changed: '+name)
    for name,names in functions.items():
        if _definitions(old/name,names)!=_definitions(current/name,names):
            raise ValueError('Native extension coding decision function changed: '+name)
    if caller_contract(old)!=caller_contract(current):
        raise ValueError('Native extension checked-prediction recipe changed')
    if request.get('oriongeno_genomewide'):
        names=('nominate','transfer_structure','finish_audit')
        if _definitions(old/'rescue_native_genomewide.py',names)!=_definitions(current/'rescue_native_genomewide.py',names):
            raise ValueError('The completed genomewide search/structure decisions changed')
    permitted=set(replay.ALLOWED_TOOL_CHANGES)|{'verify_prediction_cache_implementation','transcript_utr_implementation',
        'annevo_genome_identity_implementation','native_genomewide_implementation',
        'cached_validation_implementation','genomewide_reuse_implementation',
        'oriongeno_candidate_implementation'}
    # The GPU launcher is not invoked here. Its old bytes remain bound, while
    # frozen prediction GFF/receipts and the coding reader/alignment are strict.
    if any(before.get(k)!=after.get(k) for k in before.keys()|after.keys() if k not in permitted):
        raise ValueError('Native extension coding runtime/tool identity changed')
    proof.check()
    return dict(schema=1,scope='native_genomewide_extension',producer_source=str(old.parent.parent),
        source_files=old_files,current_implementation=after['implementation'],
        checked_prediction_recipe_sha256=hashlib.sha256(caller_contract(current).encode()).hexdigest(),
        native_parser_equivalence_required=before['transcript_utr_implementation']!=after['transcript_utr_implementation'],
        strict_qc_unchanged=True,native_alignment_and_structure_unchanged=True,
        previous_searches_reexecuted=False,native_prediction_reexecuted=False)


def scientific_signature(plan):
    """Separate immutable scientific scope from declared producer/import paths."""
    result=copy.deepcopy(plan)
    result.pop('output',None)
    request=result['request']
    history=_replay()._legacy_annotation_history(plan)
    for name in ('tools','admission_cache','native_genomewide_parent','native_genomewide_worker_imports',
                 'integration_reclassification'):
        request.pop(name,None)
    request['frozen_annotation_history']=history
    return result


def freeze_completed(plan,output,producer_source):
    r=_r()
    previous=Path(output).resolve(strict=True)
    prior=r.load(previous,inspect_only=True)
    if not prior['request'].get('oriongeno_genomewide'):
        raise ValueError('Completed-worker reuse requires the same genomewide search scope')
    if scientific_signature(prior)!=scientific_signature(plan):
        raise ValueError('Completed native workers differ from the new scientific scope')
    contract=prove_contract(prior,producer_source)
    old=Path(contract['producer_source'])/'workflow/support'
    current=Path(r.__file__).parent
    for name in ('rescue_coding_paths.py','rescue_additional_candidates.py','rescue_model_quality.py'):
        if r.digest(old/name)!=r.digest(current/name):
            raise ValueError('Completed-worker admission rules changed: '+name)
    frozen=r.frozen_prediction_cache_key(previous)
    if not frozen['species']:
        raise ValueError('No completed native workers are available to reuse')
    return {'root':str(previous),'plan_sha256':frozen['plan_sha256'],'contract':contract,
            'workers':{name:value['receipt_sha256'] for name,value in frozen['species'].items()}}


def verify_completed(plan):
    frozen=plan['request'].get('native_genomewide_worker_imports')
    if frozen is None:
        return
    if frozen!=freeze_completed(plan,frozen['root'],frozen['contract']['producer_source']):
        raise ValueError('Completed native worker generations changed')


def import_worker(root,plan,name,cpus):
    """Publish independent copies with explicit original-producer provenance."""
    r=_r()
    frozen=plan['request']['native_genomewide_worker_imports']
    if name not in frozen['workers']:
        raise ValueError('Native worker is not an explicitly frozen completed import')
    previous=Path(frozen['root'])
    prior=r.load(previous,inspect_only=True)
    source=previous/'rescued'/name
    if r.digest(source/'receipt.json')!=frozen['workers'][name]:
        raise ValueError('Completed native worker receipt changed')
    old_key=r.rescue_key(previous,prior,name)
    new_key=r.rescue_key(root,plan,name)
    def guard():
        verify_completed(plan)
        r.require_same_key(new_key,r.rescue_key(root,plan,name))
        r.require_same_key(frozen['plan_sha256'],r.plan_digest(previous,prior))
        r.require_same_key(frozen['workers'][name],r.digest(source/'receipt.json'))
    def build(tmp):
        r.copy_verified_comparison(source,tmp,old_key,cpus)
        r.atomic_json(tmp/'completed_native_worker_import.json',{
            'schema':1,'source_worker':str(source),'source_plan_sha256':frozen['plan_sha256'],
            'source_receipt_sha256':frozen['workers'][name],'source_contract':frozen['contract'],
            'scientific_scope_identical':True,'new_scientific_computation':False})
    return r.stage(root,Path('rescued')/name,new_key,build,guard,hash_workers=cpus)


def parser_equivalence(plan,name,genome):
    contract=plan['request']['admission_cache']['contract']
    if contract['scope']!='native_genomewide_extension':
        raise ValueError('Native extension parser proof has the wrong scope')
    frozen=plan['request']['oriongeno']
    if name not in frozen['predictions']:
        return {'native_prediction_available':False}
    if not contract['native_parser_equivalence_required']:
        return {'native_prediction_available':True,'coding_parser_implementation_unchanged':True}
    import annevo_gene_models as native
    native.verify_inputs(frozen,[name])
    path=frozen['predictions'][name]['gff']
    rows=native.read_paths(path,genome,engine='oriongeno')
    def fingerprint(paths):
        paths=sorted(paths,key=lambda x:x['native_transcript_id'])
        return hashlib.sha256(json.dumps(paths,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    code='''import hashlib,json,pathlib,sys
sys.path.insert(0,sys.argv[1])
import annevo_gene_models as native
import transcript_utrs,pysam
assert pathlib.Path(native.__file__).resolve()==pathlib.Path(sys.argv[1])/'annevo_gene_models.py'
assert pathlib.Path(transcript_utrs.__file__).resolve()==pathlib.Path(sys.argv[1])/'transcript_utrs.py'
with pysam.FastaFile(sys.argv[3]) as genome:rows=native.read_paths(sys.argv[2],genome,engine='oriongeno')
rows=sorted(rows,key=lambda x:x['native_transcript_id'])
print(json.dumps({'count':len(rows),'sha256':hashlib.sha256(json.dumps(rows,sort_keys=True,separators=(',',':')).encode()).hexdigest()}))
'''
    import os
    old=Path(contract['producer_source'])/'workflow/support'
    result=subprocess.run([sys.executable,'-B','-c',code,str(old),path,os.fsdecode(genome.filename)],
                          capture_output=True,text=True,check=True)
    previous=json.loads(result.stdout)
    if previous!={'count':len(rows),'sha256':fingerprint(rows)}:
        raise ValueError('Frozen native coding paths changed under the current GFF parser')
    native.verify_inputs(frozen,[name])
    return {'native_prediction_available':True,'coding_paths_identical':True,
            'native_paths':len(rows),'coding_paths_sha256':previous['sha256']}


def fresh_candidates(root,plan,name,proteins,genome,ownership,params,directory,cpus,previous):
    r=_r()
    import annevo_gene_models
    import rescue_cached_validation
    import rescue_native_genomewide as native
    parser_proof=parser_equivalence(plan,name,genome)
    models,audit=native.nominate(plan['request']['oriongeno'],name,proteins,genome,ownership,params,directory,cpus=cpus)
    memo=r.RawValidationMemo(r.validate_model,genome,plan['request']['sources'][name]['genetic_code'],params,
                             candidate_homology=r.homology_problems)
    def check(model):
        evidence=model['evidence']
        return rescue_cached_validation.checked_prediction(model,memo,genome,
            plan['request']['sources'][name]['genetic_code'],params,proteins[evidence['donor']][evidence['query']],engine=r)
    checked=[check(model) for model in models]
    transferred=native.transfer_structure(checked,proteins,genome,
        plan['request']['sources'][name]['genetic_code'],params,directory,cpus=cpus)
    transfers=[check(model) for model in transferred]
    annevo_gene_models.confirm_transferred_structure(checked,[*previous,*transfers])
    audit['prior_native_coding_parser']=parser_proof
    audit['previous_evidence_reuse']=plan['request']['admission_cache']['contract']
    return [*transfers,*checked],audit
