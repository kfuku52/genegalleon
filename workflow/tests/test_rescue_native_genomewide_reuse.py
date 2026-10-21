"""Source proof and actual coding-path equivalence guard validated reuse."""
import copy
import importlib
import shutil
import sys
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / 'support'
sys.path.insert(0, str(SUPPORT))
rescue = importlib.import_module('rescue_gene_models')
reuse = importlib.import_module('rescue_native_genomewide_reuse')
replay = importlib.import_module('rescue_admission_replay')
native_tests = importlib.import_module('test_rescue_native_genomewide')
catalog = importlib.import_module('gene_model_catalog')


@pytest.fixture
def producer(tmp_path):
    root = tmp_path / 'producer'
    shutil.copytree(SUPPORT, root / 'workflow/support')
    # Installed sources are immutable. Only this synthetic producer copy is
    # deliberately editable for corruption/semantic-change rejection tests.
    for path in (root / 'workflow/support').rglob('*'):
        if path.is_file():
            path.chmod(0o600)
    return root, {'request': {'oriongeno': {'predictions': {}}, 'tools': rescue.identities()}}


def test_source_bound_unchanged_science_is_proven(producer):
    source, plan = producer
    proof = reuse.prove_contract(plan, source)
    assert proof['strict_qc_unchanged'] and proof['native_alignment_and_structure_unchanged']
    assert proof['scope'] == 'native_genomewide_extension'
    assert not proof['native_parser_equivalence_required']
    assert str(source/'workflow/support/rescue_raw_validation.py') in proof['source_files']


@pytest.mark.parametrize('fault', ['raw-kernel', 'terminal', 'validator', 'caller-context',
                                  'native-confirmation', 'profile', 'runtime', 'unbound-source'])
def test_changed_or_unbound_scientific_evidence_is_rejected(producer, fault):
    source, plan = producer
    if fault == 'runtime':
        plan['request']['tools']['biopython'] = 'changed'
    else:
        filename, field = {
            'raw-kernel': ('rescue_raw_validation.py', 'raw_validation_implementation'),
            'terminal': ('rescue_terminal_completion.py', 'complete_terminals_implementation'),
            'validator': ('rescue_gene_models.py', 'implementation'),
            'caller-context': ('rescue_gene_models.py', 'implementation'),
            'native-confirmation': ('annevo_gene_models.py', 'annevo_candidate_implementation'),
            'profile': ('gene_model_species_profiles.py', 'species_profiles_implementation'),
            'unbound-source': ('rescue_model_quality.py', 'quality_implementation')}[fault]
        path = source/'workflow/support'/filename
        text = path.read_text()
        if fault in {'validator', 'native-confirmation'}:
            name = 'validate_model' if fault == 'validator' else 'confirm_transferred_structure'
            text = text.replace('def '+name+'(', 'def changed_'+name+'(')
        elif fault == 'caller-context':
            text = text.replace('checked_prediction(model, raw_validation, genome,',
                                'checked_prediction(model, raw_validation, None,')
        else:
            text += '\n# changed archived evidence\n'
        path.write_text(text)
        if fault != 'unbound-source':
            plan['request']['tools'][field] = rescue.digest(path)
    with pytest.raises((ValueError, KeyError)):
        reuse.prove_contract(plan, source)


@pytest.mark.parametrize('coding_change', [False, True])
def test_changed_utr_parser_must_preserve_actual_frozen_native_coding_paths(tmp_path, producer, coding_change):
    source, _ = producer
    genome, frozen, _, _ = native_tests.fixture(tmp_path)
    old_parser = source/'workflow/support/transcript_utrs.py'
    old_parser.write_text(old_parser.read_text() + ('\n_previous_load = load_models\n'
        'def load_models(path):\n'
        '    models = _previous_load(path)\n'
        '    for model in models.values():\n'
        '        a,b,p = model["cds"][0]\n'
        '        model["cds"][0] = (a+3,b,p)\n'
        '    return models\n' if coding_change else '\n# harmless UTR implementation change\n'))
    plan = {'request': {'oriongeno': frozen, 'admission_cache': {'contract': {
        'scope': 'native_genomewide_extension', 'producer_source': str(source),
        'native_parser_equivalence_required': True}}}}
    with catalog.indexed_genome(genome) as sequence:
        if coding_change:
            with pytest.raises(ValueError, match='coding paths changed'):
                reuse.parser_equivalence(plan, 'Target_species', sequence)
        else:
            proof = reuse.parser_equivalence(plan, 'Target_species', sequence)
            assert proof['coding_paths_identical'] and proof['native_paths'] == 1


def test_strict_failed_bridge_cannot_change_placement_or_final_admission():
    # Both placement and consolidation exclude strict failures before building
    # coding overlap components. Reusing their frozen rejection is equivalent.
    valid = dict(seqid='chr', strand='+', cds=[[0,9,0]], sequence='ATGAAATAA',
        coverage=1., identity=1., genetic_code=1, problems=['unanchored_genome_search'],
        evidence={'target':'Target', 'donor':'Donor', 'query':'p', 'genome_only':True})
    other = copy.deepcopy(valid)
    other['cds'] = [[18,27,0]]
    other['evidence']['query'] = 'q'
    rejected = copy.deepcopy(valid)
    rejected.update(cds=[[0,27,0]], problems=['internal_stop'])
    def decision(rows):
        rescue.reassess_unanchored_models(rows, 1)
        rescue.consolidate(rows, [], 'Target')
        return [(m['cds'],m['status'],m['problems']) for m in rows if m is not rows[-1] or 'internal_stop' not in m['problems']]
    assert decision(copy.deepcopy([valid,other])) == decision(copy.deepcopy([valid,other,rejected]))
    assert not replay.eligible({'problems':['internal_stop']})
