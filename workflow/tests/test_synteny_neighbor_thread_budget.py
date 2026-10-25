"""Every neighborhood subprocess must honor the requested CPU budget."""
import subprocess
from pathlib import Path

import pytest


@pytest.mark.parametrize("mode,representative", [("protein", ""), ("cds", "/representatives"), ("cds", "")])
def test_native_protein_families_share_annotation_cache_lock_namespace(mode, representative):
    core = (Path(__file__).resolve().parents[1] / "core/gg_gene_evolution_core.sh").read_text()
    start = core.index('  synteny_source_dir="${dir_sp_cds}"')
    setup = core[start:core.index('  run_synteny_generation=1', start)]
    script = ('set -eu\ndir_sp_cds=/cds\ndir_sp_protein_input=/proteins\n'
              'dir_synteny_gene_cache=/shared/cache\n'
              'species_protein_input_has_files() { return 0; }\n'
              'input_sequence_mode=$1\nrepresentative_inputs=$2\n'
              'for family in family_a family_b; do\n'
              ' file_og_parameters_dir=/parameters/$family\n' + setup +
              ' printf "%s\\n" "$synteny_lock_dir"\ndone\n')
    result = subprocess.run(["bash", "-c", script, "cache-lock", mode, representative],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert result.stdout.splitlines() == (["/shared/cache/.locks"] * 2 if mode == "protein" or representative else
                                         ["/parameters/family_a/synteny_locks", "/parameters/family_b/synteny_locks"])


@pytest.mark.parametrize("threads", [-2, 0, 1, 4])
@pytest.mark.parametrize("mode", ["cds", "protein"])
def test_all_neighbor_commands_keep_thread_budget_and_search_settings(tmp_path, monkeypatch, threads, mode):
    from workflow.support import synteny_neighbors as helper

    source = tmp_path / "input.fa"
    source.write_text(">Species_a_g1\nATGGCTGCT\n>Species_b_g1\nATGGCTGCT\n")
    commands = []

    def run(command):
        commands.append(command)
        if command[0] == "seqkit":
            Path(command[command.index("--out-file") + 1]).write_text(
                ">Species_a_g1\nMAA\n>Species_b_g1\nMAA\n"
            )
        if command[:2] == ["diamond", "blastp"]:
            Path(command[command.index("--out") + 1]).write_text(
                "Species_a_g1\tSpecies_b_g1\t1e-20\n"
            )

    monkeypatch.setattr(helper, "run_cmd", run)
    groups, sizes = helper.cluster_neighbors_by_similarity(
        str(source), mode, 0.01, 1, threads, str(tmp_path)
    )
    assert groups == {"Species_a_g1": "SG000001", "Species_b_g1": "SG000001"}
    assert sizes == {"SG000001": 2}
    assert [c[:2] for c in commands if c[0] == "diamond"] == [
        ["diamond", "makedb"], ["diamond", "blastp"]
    ]
    for command in commands:
        assert command[command.index("--threads") + 1] == str(max(1, threads))
    search = commands[-1]
    assert search[search.index("--evalue") + 1] == "0.01"
    assert search[search.index("--max-target-seqs") + 1] == "2000"
    assert not any(option in search for option in ("--block-size", "--index-chunks", "--fast", "--sensitive"))
