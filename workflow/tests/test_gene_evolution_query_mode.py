import subprocess
from pathlib import Path
from typing import Optional

CORE_PATH = Path(__file__).resolve().parents[1] / "core" / "gg_gene_evolution_core.sh"


def _query_fasta_stage_source() -> str:
    text = CORE_PATH.read_text(encoding="utf-8")
    start = text.index('task="Query fasta generation"')
    end = text.index('task="In-frame query BLAST', start)
    return text[start:end]


def _run_query_fasta_stage(
    mode: str,
    query_path: Optional[Path],
    capture_path: Path,
):
    script = "\n".join(
        [
            "set -eu",
            r'''
mode_gene_evolution=$1
query_path=$2
capture_path=$3
dir_output_active=/output
og_id=OG0001
gg_workspace_dir=/workspace
file_og_query_aa_fasta=/output/OG0001.query.aa.fasta
file_species_cds_store_manifest=/workspace/species_cds.json
genetic_code=1
run_extract_query_fasta=0
query_blast_method=diamond
GG_TASK_CPUS=1
if [[ -n "${query_path}" ]]; then
  file_query_gene=${query_path}
fi
gg_step_skip() {
  printf 'skipped=%s\n' "$1"
}
ensure_species_fasta_sequence_store() {
  return 0
}
gg_artifact_prepare_stage() {
  local needs_update_variable=$1
  shift 2
  printf -v "${needs_update_variable}" '%s' 0
  printf '%s\n' "$@" > "${capture_path}"
}
gg_artifact_record() {
  return 99
}
''',
            _query_fasta_stage_source(),
            "printf '%s\\n' later-orthogroup-stage-reached",
        ]
    )
    return subprocess.run(
        [
            "bash",
            "-s",
            "--",
            mode,
            "" if query_path is None else str(query_path),
            str(capture_path),
        ],
        input=script,
        text=True,
        capture_output=True,
        check=False,
    )


def test_orthogroup_mode_reaches_later_stage_without_query_variable(tmp_path):
    capture = tmp_path / "captured.txt"
    completed = _run_query_fasta_stage("orthogroup", None, capture)

    assert completed.returncode == 0, completed.stderr + completed.stdout
    assert "skipped=Query fasta generation" in completed.stdout
    assert "later-orthogroup-stage-reached" in completed.stdout
    assert not capture.exists()


def test_query2family_mode_binds_real_query_file_to_provenance(tmp_path):
    query = tmp_path / "query.faa"
    query.write_text(">query\nMPEPTIDE\n", encoding="utf-8")

    capture = tmp_path / "captured.txt"
    completed = _run_query_fasta_stage("query2family", query, capture)

    assert completed.returncode == 0, completed.stderr + completed.stdout
    captured = capture.read_text(encoding="utf-8")
    assert f"query={query}" in captured
    assert "later-orthogroup-stage-reached" in completed.stdout


def test_query2family_mode_keeps_query_file_required(tmp_path):
    completed = _run_query_fasta_stage(
        "query2family",
        None,
        tmp_path / "captured.txt",
    )

    assert completed.returncode != 0
    assert "file_query_gene" in completed.stderr
