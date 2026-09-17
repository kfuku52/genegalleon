import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_shell_species_label_parser_preserves_hybrid_binomials():
    script = (
        "source workflow/support/gg_util.sh; "
        "gg_species_name_from_path_or_dot "
        "Cenchrus_americanus_x_Cenchrus_purpureus.cds.fa.gz"
    )
    completed = subprocess.run(["bash", "-lc", script], cwd=REPO_ROOT, capture_output=True, text=True, check=False)

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "Cenchrus_americanus_x_Cenchrus_purpureus"
