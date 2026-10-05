"""Keep the plotted cutoff tied to the search that produced the neighbors."""
import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest


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


def test_neighbor_translation_uses_each_species_genetic_code(tmp_path, monkeypatch):
    from workflow.support import synteny_neighbors as mod

    genes = ['Species_a_g', 'Species_b_g']
    fasta = tmp_path / 'neighbors.fa'
    fasta.write_text('>Species_a_g\nATGTGAAAATAA\n>Species_b_g\nATGTGGAAATAA\n')
    table = tmp_path / 'codes.tsv'
    table.write_text('species\tgenetic_code\nSpecies_a\t4\nSpecies_b\t1\n')
    commands = []
    original = mod.run_cmd

    def command(args):
        commands.append(args)
        if args[:2] == ['diamond', 'makedb']:
            return None
        if args[:2] == ['diamond', 'blastp']:
            Path(args[args.index('--out') + 1]).write_text('Species_a_g\tSpecies_b_g\t1e-30\n')
            return None
        return original(args)

    monkeypatch.setattr(mod, 'run_cmd', command)
    groups, sizes = mod.cluster_neighbors_by_similarity(str(fasta), 'cds', 1e-5, 1, 1, str(tmp_path),
                                                        species_genetic_codes=mod.load_species_genetic_codes(table))
    assert groups[genes[0]] == groups[genes[1]] and sizes[groups[genes[0]]] == 2
    translated = mod.parse_fasta_subset(str(tmp_path / 'neighbors.pep.fasta'), set(genes))
    assert translated == {genes[0]: 'MWK*', genes[1]: 'MWK*'}
    assert {args[args.index('--transl-table') + 1] for args in commands if args[0] == 'seqkit'} == {'1', '4'}


@pytest.mark.parametrize('rows', ['Species_a\t999\n', 'Species_a\t4\nSpecies_a\t1\n'])
def test_neighbor_genetic_code_table_rejects_unknown_codes_and_duplicate_species(tmp_path, rows):
    from workflow.support.synteny_neighbors import load_species_genetic_codes

    table = tmp_path / 'codes.tsv'
    table.write_text('species\tgenetic_code\n' + rows)
    with pytest.raises(ValueError):
        load_species_genetic_codes(table)
