"""Real neighborhood groups remain stable within the requested CPU budget."""
import json
import random
import shutil
from pathlib import Path

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


@pytest.mark.parametrize("invalid_cds", [False, True])
def test_neighborhood_coordinates_report_utr_conflicts_and_reject_overlapping_cds(tmp_path, invalid_cds):
    import pandas

    from workflow.support import synteny_neighbors as helper

    gff = tmp_path / "Species_a.gff"
    rows = ["chr1\ts\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=g1.t1;gene_id=g1",
            "chr1\ts\tfive_prime_UTR\t2\t4\t.\t+\t.\tParent=g1.t1;gene_id=g1"]
    if invalid_cds:
        rows.append("chr1\ts\tCDS\t5\t12\t.\t+\t0\tID=cds2;Parent=g1.t1;gene_id=g1")
    gff.write_text("\n".join(rows) + "\n")
    fasta = tmp_path / "Species_a.fa"
    fasta.write_text(">Species_a_g1\nATGGCTGCT\n")
    args = dict(species_name="Species_a", species_cds_path=str(fasta), dir_sp_gff=str(tmp_path),
                cache_dir=str(tmp_path / "cache"), lock_dir=str(tmp_path / "locks"),
                gff2genestat_script=str(Path(helper.__file__).with_name("gff2genestat.py")), threads=1)
    if invalid_cds:
        with pytest.raises(RuntimeError, match="Overlapping GFF feature blocks"):
            helper.ensure_species_gene_cache(**args)
        assert not (tmp_path / "cache/Species_a.gff_info.tsv.provenance.json").exists()
        return
    path = helper.ensure_species_gene_cache(**args)
    table = pandas.read_csv(path, sep="\t", keep_default_na=False)
    assert len(table) == 1
    row = table.iloc[0]
    assert (row.start, row.end, row.strand, row.feature_size, row.feature_blocks) == (1, 9, "+", 9, "1-9")
    assert row.utr_blocks == "" and row.utr_status.startswith("conflicting:")
    manifest = json.loads(Path(path + ".provenance.json").read_text())
    assert manifest["parameters"]["structure_policy"] == "report"
    # An earlier strict-UTR cache must not be reused under the new policy.
    old = json.loads(json.dumps(manifest))
    old["parameters"].pop("structure_policy")
    Path(path + ".provenance.json").write_text(json.dumps(old))
    assert not helper.species_gene_cache_is_current(path, path + ".provenance.json",
                                                    helper.species_gene_cache_contract("Species_a", str(fasta), str(gff)))
