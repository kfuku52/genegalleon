"""Real neighborhood groups remain stable within the requested CPU budget."""
import random
import shutil

import pytest


@pytest.mark.parametrize("mode", ["cds", "protein"])
def test_native_diamond_keeps_identical_neighbor_groups_at_one_and_four_threads(tmp_path, mode):
    from workflow.support import synteny_neighbors as helper

    for tool in ("diamond", "seqkit"):
        assert shutil.which(tool), f"{tool} is required; run in the GeneGalleon container"
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    codons = dict(zip(amino_acids, [
        "GCT", "TGT", "GAT", "GAA", "TTT", "GGT", "CAT", "ATT", "AAA", "CTG",
        "ATG", "AAT", "CCT", "CAA", "CGT", "TCT", "ACT", "GTT", "TGG", "TAT",
    ], strict=True))
    rng = random.Random(314)
    shared = "M" + "".join(rng.choices(amino_acids, k=160))
    distinct = "M" + "".join(rng.choices(amino_acids, k=160))
    peptides = {"Species_a_g1": shared, "Species_b_g1": shared, "Species_c_g1": distinct}
    source = tmp_path / "input.fa"
    source.write_text("".join(
        f">{gene}\n{peptide if mode == 'protein' else ''.join(codons[aa] for aa in peptide)}\n"
        for gene, peptide in peptides.items()
    ))
    expected = ({"Species_a_g1": "SG000001", "Species_b_g1": "SG000001", "Species_c_g1": "SG000002"},
                {"SG000001": 2, "SG000002": 1})
    for threads in (1, 4):
        output = tmp_path / f"threads{threads}"
        output.mkdir()
        assert helper.cluster_neighbors_by_similarity(str(source), mode, 0.01, 1, threads, str(output)) == expected
