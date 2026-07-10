import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW_DIR = REPO_ROOT / "workflow"
CORE_DIR = WORKFLOW_DIR / "core"


def _read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_species_label_parsing_does_not_use_legacy_two_token_sed():
    util = _read_text(WORKFLOW_DIR / "support" / "gg_util.sh")
    core = _read_text(CORE_DIR / "gg_genome_evolution_core.sh")
    legacy_sed = 'sed -e "s/_/|/" -e "s/_.*//" -e "s/|/_/"'

    assert legacy_sed not in util
    assert legacy_sed not in core
    assert "gg_fasta_relabel_headers_to_species |" in core
    assert "export -f gg_species_name_from_path_or_dot _gg_strip_species_terminal_suffixes" in util


def test_shell_species_label_parser_preserves_hybrid_binomials():
    script = (
        "source workflow/support/gg_util.sh; "
        "gg_species_name_from_path_or_dot "
        "Cenchrus_americanus_x_Cenchrus_purpureus.cds.fa.gz"
    )
    completed = subprocess.run(["bash", "-lc", script], cwd=REPO_ROOT, capture_output=True, text=True, check=False)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "Cenchrus_americanus_x_Cenchrus_purpureus"


def test_tmp_cleanup_is_limited_to_known_outputs():
    orthogroup_selection = _read_text(WORKFLOW_DIR / "support" / "orthogroup_selection.py")
    orthogroup_formatter = _read_text(WORKFLOW_DIR / "support" / "orthogroup_table_formatter.py")
    genome_evolution = _read_text(CORE_DIR / "gg_genome_evolution_core.sh")

    banned_tokens = [
        "os.listdir(os.getcwd()) if f.startswith('tmp.')",
        'rm -f -- tmp.*',
        "find . -maxdepth 1 -type f -name 'tmp.*' -delete",
    ]
    combined = "\n".join([orthogroup_selection, orthogroup_formatter, genome_evolution])
    for token in banned_tokens:
        assert token not in combined

    assert "ORTHOGROUP_SELECTION_TMP_FILES" in orthogroup_selection
