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
