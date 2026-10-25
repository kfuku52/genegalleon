import gzip
import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
from Bio.SeqRecord import SeqRecord

ROOT = Path(__file__).resolve().parents[2]
HELPER = ROOT / 'workflow/support/species_tree_constraint.py'
CORE = ROOT / 'workflow/core/gg_genome_evolution_core.sh'


def fixture(tmp_path):
    alignment = tmp_path / 'input.fa'
    alignment.write_text(''.join(f'>{s}\nACDEFGHIKLMNPQRSTVWY\n' for s in 'ABCDEF'))
    constraint = tmp_path / 'constraint.nwk'
    constraint.write_text('((A,B),C,D,E,F);\n')
    return alignment, constraint


def run(alignment, constraint, *args):
    return subprocess.run([sys.executable, str(HELPER), '--alignment', str(alignment),
                           '--constraint', str(constraint), *map(str, args)], capture_output=True, text=True)


def load_module():
    spec = importlib.util.spec_from_file_location('constraint_contract', HELPER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def arguments(alignment, *, constraint='', taxa_list='', output=None):
    return SimpleNamespace(alignment=str(alignment), constraint=str(constraint), taxa_list=str(taxa_list),
                           scope='complete', tree='', output=str(output) if output is not None else None)


@pytest.mark.parametrize('compressed', [False, True], ids=['plain', 'gzip'])
def test_public_alignment_and_selection_apis_preserve_records_and_sets(tmp_path, compressed):
    module = load_module()
    path = tmp_path / ('alignment.fa.gz' if compressed else 'alignment.fa')
    payload = b''.join(b'>' + name + b' original description\nACD--X\n' for name in (b'D', b'C', b'B', b'A', b'E'))
    path.write_bytes(gzip.compress(payload, mtime=0) if compressed else payload)

    records = module.alignment(path)

    assert isinstance(records, list)
    assert all(isinstance(record, SeqRecord) for record in records)
    assert [record.id for record in records] == list('DCBAE')
    assert [record.description for record in records] == [name + ' original description' for name in 'DCBAE']
    assert [str(record.seq) for record in records] == ['ACD--X'] * 5
    assert module.selected_labels(records, '') == set('ABCDE')
    selection = tmp_path / 'taxa.txt'
    selection.write_text('A\nB\nC\nD\n')
    selected = module.selected_labels(records, selection)
    assert isinstance(selected, set) and selected == set('ABCD')


def test_alignment_preserves_existing_unique_empty_header_id_behavior(tmp_path, capsys):
    module = load_module()
    path = tmp_path / 'alignment.fa'
    path.write_text('>\nACD\n>A\nACD\n')

    # The existing error text says nonempty IDs, but only an empty alignment
    # and duplicate IDs are checked. RAM optimization must not tighten that API.
    records = module.alignment(path)

    assert [record.id for record in records] == ['', 'A']
    assert module.selected_labels(records, '') == {'', 'A'}
    module.run(arguments(path))
    assert json.loads(capsys.readouterr().out)['alignment_taxa'] == 2


@pytest.mark.parametrize('compressed', [False, True], ids=['plain', 'gzip'])
@pytest.mark.parametrize('with_output', [False, True], ids=['report-only', 'output'])
@pytest.mark.parametrize(
    ('payload', 'message'),
    [
        (b'', 'unique, nonempty sequence IDs'),
        (b'>A\nACD\n>B\nACD\n>A\nACD\n', 'unique, nonempty sequence IDs'),
        (b'>A\nACD\n>B\nACDE\n', 'same positive length'),
        (b'>A\n\n>B\n\n', 'same positive length'),
        (b'>A\nACD\n>B\nACD\n>A\nACDE\n', 'unique, nonempty sequence IDs'),
        (b'>\nACD\n>A\nACD\n>\nACD\n', 'unique, nonempty sequence IDs'),
    ],
    ids=['empty-alignment', 'late-duplicate', 'late-length', 'empty-sequences',
         'duplicate-before-length', 'duplicate-empty-id'],
)
def test_alignment_errors_precede_later_inputs_and_preserve_output(tmp_path, compressed, with_output, payload, message):
    module = load_module()
    path = tmp_path / ('alignment.fa.gz' if compressed else 'alignment.fa')
    path.write_bytes(gzip.compress(payload, mtime=0) if compressed else payload)
    output = tmp_path / 'previous.fa'
    output.write_text('previous output\n')
    args = arguments(path, constraint=tmp_path / 'missing.nwk', taxa_list=tmp_path / 'missing.txt',
                     output=output if with_output else None)

    with pytest.raises(ValueError, match=message):
        module.run(args)

    assert output.read_text() == 'previous output\n'


@pytest.mark.parametrize('with_output', [False, True], ids=['report-only', 'output'])
def test_truncated_gzip_tail_takes_priority_over_duplicate_and_length_errors(tmp_path, with_output):
    module = load_module()
    path = tmp_path / 'alignment.fa.gz'
    # More than one decoded read buffer puts the duplicate before the broken
    # trailer, so an implementation which stops scanning on a duplicate fails.
    payload = b'>A\nACD\n>B\nACD\n>A\nACD\n>late\n' + b'ACD' * 4000 + b'\n'
    path.write_bytes(gzip.compress(payload, mtime=0)[:-8])
    output = tmp_path / 'previous.fa'
    output.write_text('previous output\n')
    args = arguments(path, constraint=tmp_path / 'missing.nwk', taxa_list=tmp_path / 'missing.txt',
                     output=output if with_output else None)

    with pytest.raises(EOFError, match='end-of-stream'):
        module.run(args)

    assert output.read_text() == 'previous output\n'


@pytest.mark.parametrize('compressed', [False, True], ids=['plain', 'gzip'])
@pytest.mark.parametrize('subset', [False, True], ids=['all', 'subset'])
def test_report_only_and_output_modes_have_identical_reports(tmp_path, compressed, subset, capsys):
    module = load_module()
    alignment, constraint = fixture(tmp_path)
    payload = ''.join(f'>{name} retained description\nACD--X\n' for name in 'FEDCBA').encode()
    if compressed:
        alignment = tmp_path / 'input.fa.gz'
        alignment.write_bytes(gzip.compress(payload, mtime=0))
    else:
        alignment.write_bytes(payload)
    selection = ''
    if subset:
        selection = tmp_path / 'taxa.txt'
        selection.write_text('A\nB\nC\nD\n')
        constraint.write_text('((A,B),C,D);\n')
    output = tmp_path / 'output.fa'

    module.run(arguments(alignment, constraint=constraint, taxa_list=selection))
    report_only = json.loads(capsys.readouterr().out)
    module.run(arguments(alignment, constraint=constraint, taxa_list=selection, output=output))
    output_report = json.loads(capsys.readouterr().out)

    assert report_only == output_report
    assert report_only['alignment_taxa'] == 6
    assert report_only['alignment_columns'] == 6
    names = 'DCBA' if subset else 'FEDCBA'
    assert output.read_text() == ''.join(f'>{name} retained description\nACD--X\n' for name in names)


@pytest.mark.parametrize('constraint,message', [
    ('((A,A),C,D,E,F);', 'Duplicate'),
    ('((A,B),C,D,E,X);', 'absent'),
    ('((A,B),C,D);', 'omits'),
    ('(A,B,C,D,E,F);', 'nontrivial'),
    ('((A,B),C,D', 'constraint error'),
])
def test_invalid_constraint_retains_previous_output(tmp_path, constraint, message):
    aln, tree = fixture(tmp_path)
    tree.write_text(constraint)
    output = tmp_path / 'out.fa'
    output.write_text('previous')
    result = run(aln, tree, '--output', output)
    assert result.returncode == 2
    assert message in result.stderr
    assert output.read_text() == 'previous'


def test_selection_preserves_columns_and_complete_constraint(tmp_path):
    aln, tree = fixture(tmp_path)
    tree.write_text('((A,B),C,D);')
    taxa = tmp_path / 'taxa.txt'
    taxa.write_text('A\nB\nC\nD\n')
    output = tmp_path / 'out.fa'
    result = run(aln, tree, '--taxa-list', taxa, '--output', output)
    assert result.returncode == 0, result.stderr
    report = json.loads(result.stdout)
    assert report['alignment_taxa'] == 6 and report['selected_taxa'] == 4
    assert report['alignment_columns'] == 20
    assert len(report['assumed_unrooted_splits']) == 1
    assert output.read_text().count('>') == 4
    assert output.read_text().count('ACDEFGHIKLMNPQRSTVWY') == 4


def test_partial_constraint_checks_restricted_unrooted_splits(tmp_path):
    aln, tree = fixture(tmp_path)
    tree.write_text('((A,B),C,D);')
    result_tree = tmp_path / 'result.nwk'
    result_tree.write_text('(A,(B,E),(C,D,F));')
    assert run(aln, tree, '--scope', 'partial', '--tree', result_tree).returncode == 0
    result_tree.write_text('((A,C),(B,D),E,F);')
    result = run(aln, tree, '--scope', 'partial', '--tree', result_tree)
    assert result.returncode == 2 and 'violates' in result.stderr


@pytest.mark.parametrize('labels', ['A\nA\nC\nD\n', 'A\nB\nC\nX\n'])
def test_taxon_selection_rejects_duplicate_or_unknown_ids(tmp_path, labels):
    aln, tree = fixture(tmp_path)
    taxa = tmp_path / 'taxa.txt'
    taxa.write_text(labels)
    assert run(aln, tree, '--taxa-list', taxa).returncode == 2


def test_native_concat_stage_passes_constraint_and_invalidates_checkpoint(tmp_path):
    # Execute the real ordered core stage, with real IQ-TREE; only presentation
    # and artifact RPC plumbing are replaced by local fixture implementations.
    aln, constraint = fixture(tmp_path)
    # Distinct rows avoid a degenerate all-identical IQ-TREE input.
    aln.write_text('>A\nACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY\n'
                   '>B\nACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWF\n'
                   '>C\nACDEYGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWF\n'
                   '>D\nFCDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWF\n'
                   '>E\nACDEFGHIKLMNPQKSTVWYACDEFGHIKLMNPQRSTVWF\n'
                   '>F\nACDEFGHIKLMNPQKSTVWYACDEFGHIKLMNPQRSTVWY\n')
    core = CORE.read_text()
    start = core.index('task="IQ-TREE of the concatenated alignment with a Protein')
    end = core.index('task="Rooting of IQ-TREE\'s protein tree"', start)
    body = core[start:end]
    # Tiny fixture uses no bootstrap; production remains the original 1000/bnni.
    bin_dir = tmp_path / 'bin'
    bin_dir.mkdir()
    iqtree = subprocess.check_output(['which', 'iqtree'], text=True).strip()
    wrapper = bin_dir / 'iqtree'
    wrapper.write_text('#!/bin/bash\nargs=(); while (($#)); do\n'
                       'case "$1" in --ufboot) shift 2;; --bnni) shift;; *) args+=("$1"); shift;; esac; done\n'
                       f'{iqtree} "${{args[@]}}"\n'
                       'cp tmp.concat.pep.input.fasta.treefile tmp.concat.pep.input.fasta.contree\n')
    wrapper.chmod(0o755)
    prelude = '''set -euo pipefail
species_tree_iqtree_taxa=''
species_tree_topology_constraint_scope=complete
run_concat_iqtree_protein=1
GG_TASK_CPUS=1
GG_MEM_TOOL_GB=1
protein_model=LG
concat_iqtree_pep_provenance_args=()
gg_artifact_contract_init() { :; }
gg_artifact_prepare_stage() { printf -v "$1" 1; printf '%s\\n' "${@:3}" > "$dir_tmp/contract.txt"; }
gg_artifact_record() { :; }
gg_step_start() { :; }
gg_step_skip() { :; }
ensure_dir() { mkdir -p "$1"; }
gg_count_fasta_records() { grep -c '^>' "$1"; }
cp_out() { cp "$1" "$2"; }
'''
    env = dict(os.environ, PATH=f'{bin_dir}:{os.environ["PATH"]}',
               gg_support_dir=str(HELPER.parent), file_concat_pep=str(aln),
               species_tree_topology_constraint=str(constraint), dir_tmp=str(tmp_path),
               dir_concat_iqtree_pep=str(tmp_path / 'iqtree'), genome_evolution_provenance_dir=str(tmp_path),
               file_concat_iqtree_pep=str(tmp_path / 'published.nwk'))
    result = subprocess.run(['bash', '-c', prelude + body], env=env, cwd=tmp_path,
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    assert run(aln, constraint, '--tree', tmp_path / 'published.nwk').returncode == 0
    log = (tmp_path / 'iqtree/tmp.concat.pep.input.fasta.log').read_text()
    assert '-g ' + str(constraint) in log and '-redo' in log
    assert 'topology_constraint=' + str(constraint) in (tmp_path / 'contract.txt').read_text()
    assert json.loads((tmp_path / 'iqtree/topology_constraint.result.json').read_text())['result_constraint_verified']


def test_constraint_change_invalidates_tree_contract_only(tmp_path):
    aln, constraint = fixture(tmp_path)
    out = tmp_path / 'tree.nwk'
    out.write_text('((A,B),C,D,E,F);')
    script = ROOT / 'workflow/support/artifact_provenance.py'
    common = ['--workspace-root', str(tmp_path), '--logical-root', str(tmp_path / 'output/.gg_global_artifacts'),
              '--family-id', 'all_buscos', '--input', 'concat_peptide=' + str(aln)]
    tree_args = common + ['--manifest', str(tmp_path / 'tree.json'), '--step', 'species_tree_concat_iqtree_pep',
                          '--input', 'topology_constraint=' + str(constraint), '--output', 'tree=' + str(out),
                          '--parameter', 'constraint_scope=complete']
    msa_args = common + ['--manifest', str(tmp_path / 'msa.json'), '--step', 'species_tree_concatenate_alignment',
                         '--output', 'concat_peptide=' + str(aln)]
    def call(command, args):
        return subprocess.run([sys.executable, str(script), command, *args], capture_output=True, text=True)
    for args in [tree_args, msa_args]:
        assert call('record', args).returncode == 0
        assert call('needs-run', args).returncode == 1
    constraint.write_text('((A,C),B,D,E,F);')
    assert call('needs-run', tree_args).returncode == 3
    assert call('needs-run', msa_args).returncode == 1


@pytest.mark.parametrize('kind', ['pep', 'dna'])
def test_constraint_mode_never_auto_enables_stale_astral(tmp_path, kind):
    core = CORE.read_text()
    start = core.index('# Explicit IQ-TREE-only selection must also prevent stale-artifact auto-enabling.',
                       core.index(f'astral_{kind}_needs_update=0'))
    end = core.index('\nif [[ ${astral_', start)
    body = core[start:end]
    prelude = f'''set -euo pipefail
species_tree_topology_constraint=/fixture/constraint.nwk
species_tree_iqtree_taxa=''
astral_{kind}_provenance_args=()
gg_artifact_prepare_stage() {{ exit 23; }}
'''
    result = subprocess.run(['bash', '-c', prelude + body], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
