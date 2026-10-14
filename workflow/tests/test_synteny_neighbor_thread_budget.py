"""Every neighborhood subprocess must honor the requested CPU budget."""
from pathlib import Path

import pytest


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
