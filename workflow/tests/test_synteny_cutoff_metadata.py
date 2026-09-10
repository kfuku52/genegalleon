"""Keep the plotted cutoff tied to the search that produced the neighbors."""
import importlib.util
import sys
from pathlib import Path

import pandas as pd


def test_synteny_records_search_cutoff(tmp_path, monkeypatch):
    source = Path(__file__).parents[1] / "support" / "synteny_neighbors.py"
    spec = importlib.util.spec_from_file_location("synteny_neighbors", source)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    species = "Species_a"
    genes = [f"{species}_g{i}" for i in range(3)]
    fasta = tmp_path / f"{species}.fasta"
    fasta.write_text("".join(f">{g}\nATGATGATG\n" for g in genes))
    focal = tmp_path / "focal.fa"
    focal.write_text(f">{genes[1]}\nATGATGATG\n")
    cache = tmp_path / "genes.tsv"
    pd.DataFrame({"gene_id": genes, "chromosome": "chr1",
                  "start": [1, 11, 21], "end": [9, 19, 29]}).to_csv(cache, sep="\t", index=False)
    monkeypatch.setattr(mod, "ensure_species_gene_cache", lambda **kw: str(cache))
    observed = []

    def cluster(**kwargs):
        observed.append(kwargs["evalue_cutoff"])
        return {genes[0]: "SG1", genes[2]: "SG1"}, {"SG1": 2}

    monkeypatch.setattr(mod, "cluster_neighbors_by_similarity", cluster)
    outfile = tmp_path / "out.tsv"
    monkeypatch.setattr(sys, "argv", [str(source),
        "--focal_cds_fasta", str(focal), "--dir_sp_cds", str(tmp_path),
        "--dir_sp_gff", str(tmp_path), "--cache_dir", str(tmp_path),
        "--gff2genestat_script", "unused", "--evalue", "1e-10", "--outfile", str(outfile)])
    mod.main()
    result = pd.read_csv(outfile, sep="\t")
    assert observed == [1e-10]
    assert len(result) == 2
    assert set(result.evalue_cutoff) == set(observed)
    mod.write_empty_output(outfile)
    assert list(pd.read_csv(outfile, sep="\t").columns) == list(result.columns)
