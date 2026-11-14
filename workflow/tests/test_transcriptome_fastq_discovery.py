"""Exercise native input discovery without starting scientific processing."""

import json
import os
import subprocess
from pathlib import Path

import pytest

CORE = Path(__file__).resolve().parents[1] / "core/gg_transcriptome_generation_core.sh"


def discover(root, mode='files'):
    return subprocess.run(
        ["bash", "-c", 'GG_CORE_SOURCE_ONLY=1 source "$1"; discover_native_fastq_files "$2" "$3"',
         "test", str(CORE), str(root), mode],
        capture_output=True, text=True, check=False,
    )


def test_valid_file_links_preserve_input_names_and_sort_order(tmp_path):
    source = tmp_path / "outside"
    source.mkdir()
    reads = tmp_path / "species"
    reads.mkdir()
    (source / "raw-one").write_bytes(b"@read\nA\n+\nI\n")
    (source / "raw-two").write_bytes(b"@read\nT\n+\nI\n")
    (reads / "RUN_2.fastq.gz").symlink_to(source / "raw-two")
    (reads / "RUN_1.fastq.gz").symlink_to("../outside/raw-one")
    (reads / "RUN_3.fq").write_bytes(b"@read\nG\n+\nI\n")
    result = discover(reads)
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == [str(reads / name) for name in (
        "RUN_1.fastq.gz", "RUN_2.fastq.gz", "RUN_3.fq")]
    assert (reads / "RUN_1.fastq.gz").is_symlink()


def test_filters_hidden_nested_and_unsupported_entries(tmp_path):
    (tmp_path / "nested").mkdir()
    (tmp_path / "nested/ignored.fastq").write_bytes(b"nested")
    (tmp_path / "directory.fastq").mkdir()
    for name in (".hidden.fastq", "notes.txt", "upper.FASTQ"):
        (tmp_path / name).write_bytes(b"ignored")
    for name in ("a.fq", "b.fastq", "c.fq.gz", "d.fastq.gz"):
        (tmp_path / name).write_bytes(b"read")
    result = discover(tmp_path)
    assert result.returncode == 0, result.stderr
    assert [Path(path).name for path in result.stdout.splitlines()] == [
        "a.fq", "b.fastq", "c.fq.gz", "d.fastq.gz"]


@pytest.mark.parametrize("target_kind", ["missing", "directory", "loop"])
def test_invalid_fastq_link_fails_without_partial_inventory(tmp_path, target_kind):
    (tmp_path / "a.fastq").write_bytes(b"valid")
    target = tmp_path / "target"
    if target_kind == "directory":
        target.mkdir()
    elif target_kind == "loop":
        target.symlink_to("z.fastq.gz")
    (tmp_path / "z.fastq.gz").symlink_to(target)
    result = discover(tmp_path)
    assert result.returncode != 0
    assert result.stdout == ""
    assert "Native FASTQ input discovery failed:" in result.stderr


def test_species_directory_link_and_line_break_are_rejected(tmp_path):
    real = tmp_path / "real"
    real.mkdir()
    link = tmp_path / "link"
    link.symlink_to(real, target_is_directory=True)
    assert discover(link).returncode != 0
    (real / "bad\nname.fastq").write_bytes(b"read")
    assert "line break" in discover(real).stderr


def test_empty_inventory_is_empty_and_non_regular_file_is_ignored(tmp_path):
    os.mkfifo(tmp_path / "pipe.fastq")
    result = discover(tmp_path)
    assert result.returncode == 0, result.stderr
    assert result.stdout == ""


@pytest.mark.parametrize('change', ['none', 'retarget', 'add', 'remove'])
def test_metadata_stage_records_linked_fastqs_and_rejects_changed_bindings(tmp_path, change):
    """Run the owning metadata stage and real provenance writer, with offline integration."""
    source = tmp_path / 'source'
    source.mkdir()
    for name in ('first', 'second', 'other'):
        (source / name).write_text('@read\nA\n+\nI\n')
    reads = tmp_path / 'input/species_rnaseq/Plant_one'
    reads.mkdir(parents=True)
    for index, name in enumerate(('first', 'second'), 1):
        (reads / f'RUN_{index}.fastq').symlink_to(source / name)
    work = tmp_path / 'work'
    work.mkdir()
    output = tmp_path / 'output'
    output.mkdir()
    metadata_stage = CORE.read_text().split('task="amalgkit metadata/integrate"\n', 1)[1]
    metadata_stage = 'task="amalgkit metadata/integrate"\n' + metadata_stage.split(
        'ensure_parent_dir "${file_amalgkit_read_technology}"', 1)[0]
    shell = '''set -e
source "$support/gg_util.sh"
GG_CORE_SOURCE_ONLY=1 source "$core"
gg_support_dir=$support
gg_workspace_dir=$workspace
gg_workspace_output_dir=$workspace/output
dir_species_fastq=$workspace/input/species_rnaseq/Plant_one
transcriptome_provenance_dir=$workspace/output/provenance
file_amalgkit_metadata=$workspace/output/metadata.tsv
selected_transcriptome_mode=fastq
sp_ub=Plant_one
GG_TASK_CPUS=1
run_amalgkit_metadata_or_integrate=1
amalgkit_sra_strategy_query=''
dir_amalgkit_download_dir=$workspace/downloads
amalgkit() {
  printf 'run\\tscientific_name\\tlib_layout\\nRUN\\tPlant one\\tpaired\\n' > metadata_private_fastq.tsv
  case "$change" in
    retarget) ln -sfn "$workspace/source/other" "$dir_species_fastq/RUN_1.fastq" ;;
    add) ln -s "$workspace/source/other" "$dir_species_fastq/RUN_3.fastq" ;;
    remove) rm "$dir_species_fastq/RUN_2.fastq" ;;
  esac
}
'''
    env = dict(os.environ, support=str(CORE.parent.parent / 'support'), core=str(CORE),
               workspace=str(tmp_path), change=change)
    result = subprocess.run(['bash', '-c', shell + metadata_stage], cwd=work, env=env,
                            capture_output=True, text=True, timeout=30)
    manifest = output / 'provenance/Plant_one.metadata.json'
    if change != 'none':
        assert result.returncode != 0, result.stdout + result.stderr
        assert 'bindings changed' in result.stderr
        assert not manifest.exists()
        return
    assert result.returncode == 0, result.stdout + result.stderr
    assert manifest.is_file()
    payload = json.loads(manifest.read_text())
    assert {entry['label'] for entry in payload['inputs']} == {'fastq_1', 'fastq_2'}
    assert payload['parameters']['fastq_name_1'] == 'RUN_1.fastq'
    assert (reads / 'RUN_1.fastq').is_symlink()
    args = discover(reads, 'provenance').stdout.splitlines()
    script = CORE.parent.parent / 'support/artifact_provenance.py'
    command = ['python', str(script), 'needs-run', '--manifest', str(manifest),
               '--step', 'transcriptome_metadata', '--family-id', 'Plant_one',
               '--workspace-root', str(tmp_path), '--logical-root', str(output / '.gg_global_artifacts'),
               '--output', f'metadata={output / "metadata.tsv"}', '--parameter', 'input_mode=fastq',
               '--parameter', 'sra_strategy_query=', *args]
    assert subprocess.run(command, capture_output=True, text=True).returncode == 1
    (source / 'first').write_text('@read\nT\n+\nI\n')
    assert subprocess.run(command, capture_output=True, text=True).returncode == 3
