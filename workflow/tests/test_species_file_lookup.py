import shlex
import subprocess
from pathlib import Path


def test_species_lookup_follows_only_the_starting_directory_link(tmp_path):
    helper = Path(__file__).resolve().parents[1] / "support" / "gg_util" / "03_species_helpers.sh"
    real = tmp_path / "real"
    real.mkdir()
    target = real / "Species_example_cds.fa.gz"
    target.touch()
    (real / "Species_other_cds.fa.gz").touch()
    (real / ".Species_example_cds.fa.gz").touch()
    (real / "Species_example_link.fa.gz").symlink_to(target)
    nested = real / "nested"
    nested.mkdir()
    (nested / target.name).touch()
    link = tmp_path / "input"
    link.symlink_to(real, target_is_directory=True)
    for root in (real, link):
        result = subprocess.run(
            ["bash", "-c", f'source {shlex.quote(str(helper))}; '
             'gg_find_species_files_by_label "$1" Species_example', "test", str(root)],
            capture_output=True, text=True, check=True,
        )
        assert result.stdout.splitlines() == [str(root / target.name)]
