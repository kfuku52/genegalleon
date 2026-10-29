"""External completed annotations must not redirect the scientific workspace."""
import shlex
import subprocess
from pathlib import Path

import pytest

CORE = Path(__file__).resolve().parents[1] / 'core/gg_genome_evolution_core.sh'


@pytest.mark.parametrize('override', ['', '/completed/Arabidopsis_thaliana_annotation.tsv'])
def test_native_go_annotation_input_selection_preserves_workspace(override):
    text = CORE.read_text()
    selected = text[text.index('annotation_species_resolved=""'):text.index('# Species tree\n', text.index('annotation_species_resolved=""'))]
    prefix = '''set -euo pipefail
annotation_species=Arabidopsis_thaliana
gg_workspace_output_dir=/original/workspace/output
dir_sp_cds=/original/workspace/input/species_cds
input_sequence_mode=cds
gg_species_names_from_annotation_dir() { echo Arabidopsis_thaliana; }
gg_resolve_annotation_species() { echo "$1"; }
'''
    result = subprocess.run(['bash', '-c', prefix + '\ngo_annotation_file=' + shlex.quote(override) + '\n' + selected +
                             '\nprintf "%s\\n%s\\n%s\\n" "$file_go_annotation" "$annotation_species_resolved" "$gg_workspace_output_dir"'],
                            text=True, capture_output=True, check=True)
    assert result.stdout.splitlines() == [override or '/original/workspace/output/species_cds_annotation/Arabidopsis_thaliana_annotation.tsv',
                                          'Arabidopsis_thaliana', '/original/workspace/output']


def test_relative_go_annotation_override_stops_before_analysis():
    text = CORE.read_text()
    selected = text[text.index('if [[ -n "${go_annotation_file:-}" ]]'):text.index('# Species tree\n', text.index('if [[ -n "${go_annotation_file:-}" ]]'))]
    result = subprocess.run(['bash', '-c', 'go_annotation_file=relative.tsv\n' + selected], text=True, capture_output=True)
    assert result.returncode == 1 and 'absolute annotation-table path' in result.stderr


@pytest.mark.parametrize('initial_run', ['0', '1'])
def test_fixed_go_universe_is_forwarded_when_stale_rebuild_enables_stage(tmp_path, initial_run):
    text = CORE.read_text()
    stage = text[text.index('task="GO enrichment analysis"'):text.index('if [[ ${delete_tmp_dir}', text.index('task="GO enrichment analysis"'))]
    frozen = tmp_path / 'frozen.tsv'
    frozen.write_text('fixed catalog')
    assignments = {
        'go_enrichment_method': 'event_all_go', 'run_go_enrichment': initial_run,
        'go_fixed_universe_file': str(frozen), 'genome_evolution_provenance_dir': str(tmp_path),
        'dir_cafe_output': str(tmp_path), 'file_gene_id': str(tmp_path / 'ids.tsv'),
        'file_go_annotation': str(tmp_path / 'ann.tsv'), 'file_go_enrichment_significant': str(tmp_path / 'out/hits.tsv'),
        'target_branch_go': '<37>', 'change_direction_go': 'decrease', 'go_category': 'BP,MF,CC',
        'gg_support_dir': str(tmp_path), 'go_family_alpha': '.05', 'annotation_species_resolved': 'Arabidopsis_thaliana',
        'dir_go_enrichment': str(tmp_path / 'out')}
    prelude = "set -euo pipefail\n" + "\n".join(k+'='+shlex.quote(v) for k,v in assignments.items()) + "\n"
    prelude += r'''
disable_if_no_input_file() { :; }
gg_artifact_contract_init() { go_enrichment_provenance_args=(); }
gg_artifact_add_input_if_present() { :; }
gg_artifact_prepare_stage() { printf -v "$1" '%s' 1; printf -v "$2" '%s' 1; }
gg_artifact_record() { :; }
gg_step_start() { :; }
gg_step_skip() { :; }
ensure_dir() { mkdir -p "$1"; }
Rscript() { [[ $# -eq 13 && "${!#}" == "$go_fixed_universe_file" ]] || return 77; echo FIXED_GO_ARGUMENT_PRESERVED; }
'''
    result = subprocess.run(['bash', '-c', prelude+stage], text=True, capture_output=True)
    assert result.returncode == 0, result.stdout+result.stderr
    assert 'FIXED_GO_ARGUMENT_PRESERVED' in result.stdout
