"""Run the actual inference stages and verify ML topology publication and reuse."""
import os
import random
import subprocess
import sys
from pathlib import Path

import ete4
import pytest

ROOT = Path(__file__).resolve().parents[2]
SUPPORT = ROOT / 'workflow/support'


def stage_body(stage):
    if stage == 'gene':
        core = (ROOT / 'workflow/core/gg_gene_evolution_core.sh').read_text()
        start = core.index('task="IQ-TREE"')
        end = core.index('gene_nwkit_identity=', start)
    else:
        core = (ROOT / 'workflow/core/gg_genome_evolution_core.sh').read_text()
        model, label = ('Protein', 'protein') if stage == 'protein' else ('DNA', 'DNA')
        start = core.index(f'task="IQ-TREE of the concatenated alignment with a {model}')
        end = core.index(f'task="Rooting of IQ-TREE\'s {label} tree"', start)
    return core[start:end]


PRELUDE = '''set -euo pipefail
run_iqtree=1
run_concat_iqtree_protein=1
run_concat_iqtree_dna=1
input_sequence_mode=protein
genetic_code=1
generax_model=LG
protein_model=LG
nucleotide_model=GTR+G
GG_TASK_CPUS=1
GG_MEM_TOOL_GB=1
og_id=OG0001
species_tree_iqtree_taxa=''
species_tree_topology_constraint=''
species_tree_topology_constraint_scope=complete
disable_if_no_input_file() { :; }
gg_step_start() { :; }
gg_step_skip() { :; }
gg_artifact_contract_init() {
  local -n contract=$1
  contract=(--step "$2" --family-id "$3" --manifest "$4"
    --logical-root "$dir_output_active" --workspace-root "$gg_workspace_dir")
}
gg_artifact_prepare_stage() {
  printf -v "$1" 1
  printf '%s\\n' "${@:3}" > "$dir_tmp/contract.txt"
  printf '%s\\n' "$artifact_legacy_policy" > "$dir_tmp/legacy-policy.txt"
}
gg_artifact_record() { :; }
gg_count_fasta_records() { grep -c '^>' "$1"; }
ensure_dir() { mkdir -p "$1"; }
assert_gene_evolution_aa_model_for_protein_mode() { :; }
build_iqtree_mem_args() { IQTREE_MEM_ARGS=(-mem 1G); }
cp_out() { printf '%s\\n' "$1" > "$dir_tmp/published-source.txt"; cp "$1" "$2"; }
'''


def write_alignment(path, stage, ntaxa):
    # Deterministic synthetic substitutions on three pairs of tips. No external
    # sequences or scientific fixture regeneration are involved.
    rng = random.Random(54)
    alphabet = 'ACGT' if stage == 'dna' else 'ACDEFGHIKLMNPQRSTVWY'
    ancestor = [rng.choice(alphabet) for _ in range(240)]

    def mutate(sequence, probability):
        return [rng.choice(alphabet.replace(base, '')) if rng.random() < probability else base
                for base in sequence]

    sequences = []
    for _ in range(3):
        pair = mutate(ancestor, 0.12)
        sequences.extend(mutate(pair, 0.04) for _ in range(2))
    path.write_text(''.join(f'>{name}\n{"".join(sequence)}\n'
                            for name, sequence in zip('ABCDEF'[:ntaxa], sequences[:ntaxa], strict=True)))


def assert_contract_requires_regeneration(tmp_path, parameter):
    args = (tmp_path / 'contract.txt').read_text().splitlines()
    assert parameter in args
    old_args = list(args)
    index = old_args.index(parameter)
    del old_args[index - 1:index + 1]

    def provenance(action, parameters):
        return subprocess.run([sys.executable, str(SUPPORT / 'artifact_provenance.py'), action, *parameters],
                              cwd=tmp_path, capture_output=True, text=True)

    recorded = provenance('record', old_args)
    assert recorded.returncode == 0, recorded.stderr
    assert provenance('needs-run', args).returncode == 3
    assert provenance('needs-run', args + ['--stale-policy', 'rebuild']).returncode == 0
    current = provenance('record', args)
    assert current.returncode == 0, current.stderr
    assert provenance('needs-run', args).returncode == 1


@pytest.mark.parametrize('stage,generax,ntaxa,fast_threshold,bootstrap', [
    ('gene', 0, 6, 2000, True),
    ('gene', 1, 6, 2000, False),
    ('gene', 0, 6, 5, False),
    ('gene', 0, 3, 2000, False),
    ('protein', 0, 6, 2000, True),
    ('dna', 0, 6, 2000, True),
])
def test_native_stages_publish_ml_tree_and_invalidate_consensus_contract(
        tmp_path, stage, generax, ntaxa, fast_threshold, bootstrap):
    alignment = tmp_path / 'input.fa'
    write_alignment(alignment, stage, ntaxa)
    output_dir = tmp_path / 'output'
    output_dir.mkdir()
    published = output_dir / 'iqtree_tree' / 'OG0001_iqtree.nwk'
    published.parent.mkdir()
    iqtree_dir = output_dir / 'iqtree'
    env = dict(os.environ, gg_support_dir=str(SUPPORT), gg_workspace_dir=str(tmp_path),
        dir_output_active=str(output_dir), dir_tmp=str(tmp_path),
        genome_evolution_provenance_dir=str(output_dir),
        file_og_trimmed_aln_analysis=str(alignment), file_og_iqtree_tree=str(published),
        file_concat_pep=str(alignment), file_concat_cds=str(alignment),
        dir_concat_iqtree_pep=str(iqtree_dir), dir_concat_iqtree_dna=str(iqtree_dir),
        file_concat_iqtree_pep=str(published), file_concat_iqtree_dna=str(published),
        run_generax=str(generax), iqtree_fast_mode_gt=str(fast_threshold), artifact_legacy_policy='adopt')
    completed = subprocess.run(['bash', '-c', PRELUDE + stage_body(stage)], env=env,
                               cwd=tmp_path, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    if stage != 'gene':
        # A requested rebuild must replace a completed IQ-TREE checkpoint too.
        rebuilt = subprocess.run(['bash', '-c', PRELUDE + stage_body(stage)], env=env,
                                 cwd=tmp_path, capture_output=True, text=True)
        assert rebuilt.returncode == 0, rebuilt.stdout + rebuilt.stderr
    source = Path((tmp_path / 'published-source.txt').read_text().strip())
    if not source.is_absolute():
        source = (tmp_path if stage == 'gene' else iqtree_dir) / source
    assert source.suffix == '.treefile'
    assert published.read_bytes() == source.read_bytes()
    assert (tmp_path / 'legacy-policy.txt').read_text().strip() == 'rebuild'
    tree = ete4.Tree(published.read_text(), parser=0)
    assert set(tree.leaf_names()) == set('ABCDEF'[:ntaxa])
    internal = [node for node in tree.traverse() if node.up is not None and not node.is_leaf]
    if bootstrap:
        assert internal and all(node.support is not None and 0 <= node.support <= 100 for node in internal)
        assert source.with_suffix('.contree').is_file()
    log = source.with_suffix('.log').read_text()
    assert ('--ufboot 1000' in log) == bootstrap
    assert ('--fast' in log) == (stage == 'gene' and ntaxa > fast_threshold)

    # Exercise the actual provenance comparison with the stage's declared
    # arguments: an old selection contract must stop or request a rebuild.
    previous = published.read_bytes()
    assert_contract_requires_regeneration(tmp_path, 'tree_output=treefile')
    assert published.read_bytes() == previous


def test_summary_contract_invalidates_old_sister_semantics(tmp_path):
    output_dir = tmp_path / 'output'
    branch = output_dir / 'stat_branch' / 'OG0001_stat.branch.tsv'
    summary = output_dir / 'stat_tree' / 'OG0001_stat.tree.tsv'
    for path in (branch, summary):
        path.parent.mkdir(parents=True)
        path.write_text('preserve prior statistics\n')
    core = (ROOT / 'workflow/core/gg_gene_evolution_core.sh').read_text()
    start = core.index('summary_provenance_args=(')
    end = core.index('disable_if_no_input_file "run_summary"', start)
    env = dict(os.environ, dir_output_active=str(output_dir), gg_workspace_dir=str(tmp_path),
        dir_tmp=str(tmp_path), file_og_stat_branch=str(branch), file_og_stat_tree=str(summary),
        file_og_summary_provenance=str(output_dir / 'artifact_provenance' / 'OG0001.summary_statistics.json'),
        treevis_clade_ortholog_prefix='singlecopy', mode_gene_evolution='orthogroup', treevis_query_marker='0',
        query_blast_coverage='0.5', artifact_legacy_policy='adopt')
    completed = subprocess.run(['bash', '-c', PRELUDE + 'summary_input_files=()\n' + core[start:end]],
                               env=env, cwd=tmp_path, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert (tmp_path / 'legacy-policy.txt').read_text().strip() == 'rebuild'
    assert_contract_requires_regeneration(tmp_path, 'sister_semantics=unique_or_sentinel')
    assert branch.read_text() == summary.read_text() == 'preserve prior statistics\n'
