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


class IncompatibleEvidence(ValueError):
    """An intact producer cannot provide the required equivalence proof."""


def archived_sources(output):
    """Discover archived code next to or inside its producer output."""
    output = Path(output).resolve(strict=True)
    for source in (output.parent / 'source', output / 'source'):
        if (source / 'workflow/support').is_dir():
            yield source


# The compact hit index preserves the old sorted set of scored protein pairs,
# length-rejected counts, and all downstream models. Qualify only the exact
# verified source transition; other changes still need the function proof.
_COMPACT_PROTEIN_PAIRS_REUSE = (
    '4f689c15a9ac56bd665dd6462dc2299a4ad14b28d6013d235e61136b06179931',
    '08ccc0a672ddc6a152bf95b1525fc7699c55f3323dee467bd0d31701d1fc49fa',
)


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
        raise IncompatibleEvidence('Coding QC requires one checked-prediction caller')
    function=nodes[0]
    if any(isinstance(n,ast.Attribute) and n.attr=='checked_prediction' for n in ast.walk(function)):
        expected=ast.parse('''def checked_prediction(model):
    evidence = model["evidence"]
    return rescue_cached_validation.checked_prediction(model, raw_validation, genome,
        source['genetic_code'], params, proteins[evidence['donor']][evidence['query']],
        engine=sys.modules[__name__])
''').body[0]
        if ast.dump(function,include_attributes=False)!=ast.dump(expected,include_attributes=False):
            raise IncompatibleEvidence('Coding QC worker context changed')
        helper=support/'rescue_cached_validation.py'
        functions=[n for n in ast.parse(helper.read_text()).body
                   if isinstance(n,ast.FunctionDef) and n.name=='checked_prediction']
        if len(functions)!=1:
            raise IncompatibleEvidence('Coding QC worker helper is unavailable')
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
        raise IncompatibleEvidence('Validated native extension requires an Orion-only non-GeMoMa producer')
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
    if 'source_verification_implementation' in before:
        bindings['rescue_source_verification.py']='source_verification_implementation'
    if 'cached_validation_implementation' in before:
        bindings['rescue_cached_validation.py']='cached_validation_implementation'
        if before['cached_validation_implementation']!=after['cached_validation_implementation']:
            raise IncompatibleEvidence('The producer cached-validation dispatcher changed')
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
            raise IncompatibleEvidence('Native extension producer source differs from its plan: '+name)
        old_files[str(old/name)]=value
    for name in strict:
        if observed[str(old/name)]!=observed[str(current/name)]:
            raise IncompatibleEvidence('Native extension coding evidence kernel changed: '+name)
    for name,names in functions.items():
        try:
            previous_functions = _definitions(old / name, names)
        except ValueError as exc:
            raise IncompatibleEvidence('Native extension coding decision function unavailable: ' + name) from exc
        if previous_functions != _definitions(current / name, names):
            raise IncompatibleEvidence('Native extension coding decision function changed: ' + name)
    if caller_contract(old)!=caller_contract(current):
        raise IncompatibleEvidence('Native extension checked-prediction recipe changed')
    if request.get('oriongeno_genomewide'):
        module = 'rescue_native_genomewide.py'
        old_hash, current_hash = observed[str(old / module)], observed[str(current / module)]
        if old_hash != current_hash and (old_hash, current_hash) != _COMPACT_PROTEIN_PAIRS_REUSE:
            names = ('nominate', '_candidate_protein_pairs', 'transfer_structure', 'finish_audit')
            try:
                unchanged = _definitions(old / module, names) == _definitions(current / module, names)
            except ValueError as exc:
                raise IncompatibleEvidence('Completed genomewide search/structure functions unavailable') from exc
            if not unchanged:
                raise IncompatibleEvidence('The completed genomewide search/structure decisions changed')
    permitted=set(replay.ALLOWED_TOOL_CHANGES)|{'verify_prediction_cache_implementation','transcript_utr_implementation',
        'annevo_genome_identity_implementation','native_genomewide_implementation',
        'cached_validation_implementation','genomewide_reuse_implementation',
        'oriongeno_candidate_implementation','nwkit'}
    if 'source_verification_implementation' not in before:
        permitted.add('source_verification_implementation')
    # The GPU launcher is not invoked here. Its old bytes remain bound, while
    # frozen prediction GFF/receipts and the coding reader/alignment are strict.
    # NWKit belongs to initial tree/donor selection. This route byte-imports the
    # frozen selection and requires the complete plan scope unchanged; no tree
    # command executes during validated extension or completed-worker import.
    if any(before.get(k)!=after.get(k) for k in before.keys()|after.keys() if k not in permitted):
        raise IncompatibleEvidence('Native extension coding runtime/tool identity changed')
    proof.check()
    return dict(schema=1,scope='native_genomewide_extension',producer_source=str(old.parent.parent),
        source_files=old_files,current_implementation=after['implementation'],
        checked_prediction_recipe_sha256=hashlib.sha256(caller_contract(current).encode()).hexdigest(),
        native_parser_equivalence_required=before['transcript_utr_implementation']!=after['transcript_utr_implementation'],
        frozen_tree_selection_reused=True,unused_nwkit_versions={'producer':before.get('nwkit'),'consumer':after.get('nwkit')},
        strict_qc_unchanged=True,native_alignment_and_structure_unchanged=True,
        previous_searches_reexecuted=False,native_prediction_reexecuted=False)


def scientific_signature(plan):
    """Separate immutable scientific scope from declared producer/import paths."""
    result=copy.deepcopy(plan)
    result.pop('output',None)
    result.pop('source_versions',None)  # The frozen request still compares every source SHA-256.
    request=result['request']
    history=_replay()._legacy_annotation_history(plan)
    for name in ('tools','admission_cache','native_genomewide_parent','native_genomewide_worker_imports',
                 'integration_reclassification'):
        request.pop(name,None)
    request['frozen_annotation_history']=history
    return result


def freeze_completed(plan, output, producer_source, *, names=None):
    r=_r()
    previous=Path(output).resolve(strict=True)
    prior=r.load(previous,inspect_only=True)
    if not prior['request'].get('oriongeno_genomewide'):
        raise ValueError('Completed-worker reuse requires the same genomewide search scope')
    if scientific_signature(prior)!=scientific_signature(plan):
        raise ValueError('Completed native workers differ from the new scientific scope')
    if producer_source is None:
        sources = list(archived_sources(previous))
        if not sources:
            raise ValueError('Completed native workers require their archived producer source')
        for source in sources:
            try:
                contract = prove_contract(prior, source)
                break
            except (IncompatibleEvidence, FileNotFoundError):
                continue
        else:
            raise ValueError('No archived source matches the completed native workers')
    else:
        contract = prove_contract(prior, producer_source)
    old=Path(contract['producer_source'])/'workflow/support'
    current=Path(r.__file__).parent
    for name in ('rescue_coding_paths.py', 'rescue_additional_candidates.py', 'rescue_model_quality.py'):
        if r.digest(old / name) != r.digest(current / name):
            raise ValueError('Completed-worker admission rules changed: ' + name)
    decisions = ('consolidate', 'revalidate_retained_coding_paths', 'split_model_candidates',
                 'automatic_model_repair_targets', 'validate_model_repair_predictions', 'model_repair_candidates')
    if _definitions(old / 'rescue_gene_models.py', decisions) != _definitions(current / 'rescue_gene_models.py', decisions):
        raise ValueError('Completed-worker admission or repair decisions changed')
    frozen = r.frozen_prediction_cache_key(previous, names)
    r.require_same_key(frozen['plan_sha256'], r.plan_digest(previous, prior))
    if not frozen['species']:
        raise ValueError('No completed native workers are available to reuse')
    return {'root':str(previous),'plan_sha256':frozen['plan_sha256'],'contract':contract,
            'workers':{name:value['receipt_sha256'] for name,value in frozen['species'].items()}}


def verify_completed(plan, *, name=None):
    frozen=plan['request'].get('native_genomewide_worker_imports')
    if frozen is None:
        return
    expected = frozen if name is None else {**frozen, 'workers': {name: frozen['workers'][name]}}
    current = freeze_completed(plan, frozen['root'], frozen['contract']['producer_source'],
                               names=expected['workers'])
    if expected != current:
        raise ValueError('Completed native worker generations changed')


def import_worker(root, plan, name, cpus, *, dependencies):
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
    old_key = r.rescue_key(previous, prior, name, dependency_snapshot=dependencies)
    new_key = r.rescue_key(root, plan, name, dependency_snapshot=dependencies)

    def guard():
        dependencies.check_all()
        verify_completed(plan, name=name)
        r.require_same_key(new_key, r.rescue_key(root, plan, name, dependency_snapshot=dependencies))
        r.require_same_key(frozen['plan_sha256'],r.plan_digest(previous,prior))
        r.require_same_key(frozen['workers'][name],r.digest(source/'receipt.json'))
    def build(tmp):
        contract = frozen['contract']
        audit = {'schema': 1, 'source_worker': str(source), 'source_plan_sha256': frozen['plan_sha256'],
                 'source_receipt_sha256': frozen['workers'][name], 'source_contract': contract,
                 'scientific_scope_identical': True, 'new_scientific_computation': False}
        if contract['native_parser_equivalence_required']:
            from gene_model_catalog import indexed_genome

            r.verify_sources(plan, [name], ['genome', 'fasta', 'gff'])
            with indexed_genome(plan['request']['sources'][name]['genome']) as genome:
                audit['native_coding_parser'] = parser_equivalence(plan, name, genome, contract=contract)
        r.copy_verified_comparison(source, tmp, old_key, cpus, receipt_sha256=frozen['workers'][name])
        r.atomic_json(tmp / 'completed_native_worker_import.json', audit)
    return r.stage(root,Path('rescued')/name,new_key,build,guard,hash_workers=cpus)


def parser_equivalence(plan, name, genome, *, contract=None):
    if contract is None:
        contract = plan['request']['admission_cache']['contract']
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
    native.verify_inputs(frozen,[name])
    if previous!={'count':len(rows),'sha256':fingerprint(rows)}:
        raise IncompatibleEvidence('Frozen native coding paths changed under the current GFF parser')
    return {'native_prediction_available':True,'coding_paths_identical':True,
            'native_paths':len(rows),'coding_paths_sha256':previous['sha256']}


def fresh_candidates(plan, name, proteins, genome, ownership, params, directory, cpus, previous, *, parser_proof=None):
    r=_r()
    import annevo_gene_models
    import rescue_cached_validation
    import rescue_native_genomewide as native
    if parser_proof is None:
        parser_proof = parser_equivalence(plan, name, genome)
    memo=r.RawValidationMemo(r.validate_model,genome,plan['request']['sources'][name]['genetic_code'],params,
                             candidate_homology=r.homology_problems)
    def check(model):
        evidence=model['evidence']
        return rescue_cached_validation.checked_prediction(model,memo,genome,
            plan['request']['sources'][name]['genetic_code'],params,proteins[evidence['donor']][evidence['query']],engine=r)
    checked, audit = native.nominate(plan['request']['oriongeno'], name, proteins, genome, ownership,
                                    params, directory, cpus=cpus, validate=check)
    transferred=native.transfer_structure(checked,proteins,genome,
        plan['request']['sources'][name]['genetic_code'],params,directory,cpus=cpus)
    transfers=[check(model) for model in transferred]
    annevo_gene_models.confirm_transferred_structure(checked,[*previous,*transfers])
    audit['prior_native_coding_parser']=parser_proof
    audit['previous_evidence_reuse']=plan['request']['admission_cache']['contract']
    return [*transfers,*checked],audit
