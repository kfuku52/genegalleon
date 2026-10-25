import gzip
import sys

import pandas
import pytest

from workflow.support import remove_contaminated_sequences as contamination
from workflow.support.remove_contaminated_sequences import (
    build_taxid_rank_lookup,
    resolve_lineage_rank_index,
    resolve_rank_taxid_from_lineage,
)


class FakeNcbi:
    def __init__(self, rank_map):
        self.rank_map = rank_map

    def get_rank(self, lineage_taxids):
        return {taxid: self.rank_map.get(taxid, "no rank") for taxid in lineage_taxids}


def test_resolve_rank_taxid_from_lineage_missing_returns_minus_one():
    ncbi = FakeNcbi({})
    assert resolve_rank_taxid_from_lineage(ncbi, "class", "nan") == -1
    assert resolve_rank_taxid_from_lineage(ncbi, "class", "") == -1


def test_resolve_rank_taxid_from_lineage_returns_rank_taxid_when_present():
    ncbi = FakeNcbi({2: "superkingdom", 2157: "class"})
    assert resolve_rank_taxid_from_lineage(ncbi, "class", "2;2157;9999") == 2157


def test_resolve_rank_taxid_from_lineage_returns_zero_when_rank_absent():
    ncbi = FakeNcbi({2: "superkingdom", 1224: "phylum"})
    assert resolve_rank_taxid_from_lineage(ncbi, "class", "2;1224;1236") == 0


def test_resolve_rank_taxid_from_lineage_treats_domain_and_superkingdom_as_aliases():
    ncbi = FakeNcbi({2: "domain", 33090: "kingdom"})
    assert resolve_rank_taxid_from_lineage(ncbi, "superkingdom", "1;2;33090") == 2


def test_build_taxid_rank_lookup_and_bulk_rank_resolution():
    ncbi = FakeNcbi({2: "superkingdom", 1224: "phylum", 28211: "class"})
    taxid_to_rank = build_taxid_rank_lookup(ncbi, ["2;1224;28211", "2;1224"])
    assert resolve_rank_taxid_from_lineage(ncbi, "class", "2;1224;28211", taxid_to_rank=taxid_to_rank) == 28211


def test_resolve_lineage_rank_index_treats_domain_and_superkingdom_as_aliases():
    df_lineage = pandas.DataFrame(
        {
            "rank": ["no rank", "domain", "kingdom", "species"],
            "name": ["root", "Eukaryota", "Viridiplantae", "Arabidopsis thaliana"],
            "taxid": [1, 2759, 33090, 3702],
        }
    )
    assert resolve_lineage_rank_index(df_lineage, "superkingdom") == 1


def test_resolve_lineage_rank_index_raises_clear_error_when_rank_absent():
    df_lineage = pandas.DataFrame(
        {
            "rank": ["no rank", "domain", "kingdom", "species"],
            "name": ["root", "Eukaryota", "Viridiplantae", "Arabidopsis thaliana"],
            "taxid": [1, 2759, 33090, 3702],
        }
    )
    try:
        resolve_lineage_rank_index(df_lineage, "class")
    except ValueError as exc:
        assert "Requested rank 'class' was not found" in str(exc)
    else:
        raise AssertionError("Expected ValueError for missing rank")


class MainNcbi(FakeNcbi):
    def __init__(self):
        super().__init__({1: "no rank", 2: "domain", 10: "class", 11: "species", 20: "class", 21: "species"})

    def get_name_translator(self, names):
        return {name: [11] for name in names if name == "Genus species"}

    def get_lineage(self, taxid):
        return [1, 2, 10, 11]

    def get_taxid_translator(self, taxids):
        names = {1: "root", 2: "Eukaryota", 10: "FocalClass", 11: "Genus species"}
        return {taxid: names[taxid] for taxid in taxids}


def prepare_main(tmp_path, monkeypatch, fasta, rows, *, extension="fa", rename=False, verbose=False):
    """Run production main with fixed taxonomy and ordinary on-disk inputs."""
    source = tmp_path / ("source." + extension)
    if extension in {"fa.gz", "fa.gzipped"}:
        source.write_bytes(gzip.compress(fasta.encode()))
    else:
        source.write_text(fasta)
    taxonomy = tmp_path / "taxonomy.tsv"
    fx = tmp_path / "fx.tsv"
    taxonomy_rows = []
    lengths = {}
    for identifier, length, status in rows:
        taxid, name, lineage = {"focal": (11, "Genus species", "1;2;10;11"),
                               "foreign": (21, "Other species", "1;2;20;21"),
                               "unknown": (0, "unclassified", "")}[status]
        taxonomy_rows.append(f"{identifier}\t{taxid}\t-\t{name}\t-\t-\t-\t-\t{lineage}\n")
        lengths[identifier] = length
    taxonomy.write_text("".join(taxonomy_rows))
    fx.write_text("#id\tlength\tGC\n" + "".join(f"{identifier}\t{length}\t50\n"
                                               for identifier, length in lengths.items()))
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("GG_TAXONOMY_DBFILE", raising=False)
    monkeypatch.setattr(contamination, "ncbi_taxa", lambda dbfile=None: MainNcbi())
    monkeypatch.setattr(sys, "argv", ["remove_contaminated_sequences.py", "--fasta_file", str(source),
        "--mmseqs2taxonomy_tsv", str(taxonomy), "--fx2tab_tsv", str(fx),
        "--species_name", "Genus_species", "--rank", "class", "--rename_prefix", "clean",
        "--rename_seq", "yes" if rename else "no", "--verbose", "yes" if verbose else "no"])
    return source


@pytest.mark.parametrize("extension", ["fa", "fa.gz", "fa.gzipped"])
@pytest.mark.parametrize("rename", [False, True])
@pytest.mark.parametrize("verbose", [False, True])
def test_main_fasta_order_renaming_and_unclassified(tmp_path, monkeypatch, capsys, extension, rename, verbose):
    prepare_main(tmp_path, monkeypatch,
        ">short original description\nACGT\n>foreign\nAAAAAAA\n>unknown\nNNNNNN\n>long description\nACGTACGTA\n",
        [("short", 4, "focal"), ("foreign", 7, "foreign"), ("unknown", 6, "unknown"), ("long", 9, "focal")],
        extension=extension, rename=rename, verbose=verbose)
    contamination.main()
    identifiers = ["clean1", "clean2", "clean3"] if rename else ["long", "unknown", "short"]
    assert (tmp_path / "clean_sequences.fa").read_text() == (
        f">{identifiers[0]}\nACGTACGTA\n>{identifiers[1]}\nNNNNNN\n>{identifiers[2]}\nACGT\n")
    lineage = pandas.read_csv(tmp_path / "lineage_compatibility.tsv", sep="\t")
    assert lineage["query"].tolist() == ["long", "foreign", "unknown", "short"]
    assert lineage["is_compatible_lineage"].tolist() == [True, False, True, True]
    if rename:
        assert lineage["new_seq_name"].fillna("").tolist() == ["clean1", "", "clean2", "clean3"]
    stdout = capsys.readouterr().out
    assert ("Renaming: long -> clean1 (9 bp)" in stdout) is (rename and verbose)


def test_main_duplicate_fasta_identifier_keeps_last_record(tmp_path, monkeypatch):
    prepare_main(tmp_path, monkeypatch, ">same first\nAAAAA\n>same last\nTTTTT\n",
                 [("same", 5, "focal")])
    contamination.main()
    assert (tmp_path / "clean_sequences.fa").read_text() == ">same\nTTTTT\n"


@pytest.mark.parametrize("rename", [False, True])
def test_main_duplicate_taxonomy_rows_preserve_record_mutation(tmp_path, monkeypatch, rename):
    prepare_main(tmp_path, monkeypatch, ">same\nACGT\n",
                 [("same", 4, "focal"), ("same", 4, "focal")], rename=rename)
    if rename:
        with pytest.raises(KeyError, match="clean1"):
            contamination.main()
        assert (tmp_path / "clean_sequences.fa").read_text() == ">clean1\nACGT\n"
    else:
        contamination.main()
        assert (tmp_path / "clean_sequences.fa").read_text() == ">same\nACGT\n>same\nACGT\n"
    assert len(pandas.read_csv(tmp_path / "lineage_compatibility.tsv", sep="\t")) == 2


def test_main_missing_fasta_identifier_preserves_keyerror(tmp_path, monkeypatch):
    prepare_main(tmp_path, monkeypatch, ">different\nACGT\n", [("wanted", 4, "focal")])
    with pytest.raises(KeyError, match="wanted"):
        contamination.main()
    assert (tmp_path / "clean_sequences.fa").read_bytes() == b""
    assert (tmp_path / "lineage_compatibility.tsv").is_file()


def test_main_corrupt_gzip_does_not_replace_existing_fasta(tmp_path, monkeypatch):
    source = prepare_main(tmp_path, monkeypatch, ">same\nACGT\n", [("same", 4, "focal")], extension="fa.gz")
    source.write_bytes(source.read_bytes()[:-6])
    (tmp_path / "clean_sequences.fa").write_text("existing FASTA\n")
    with pytest.raises(EOFError):
        contamination.main()
    assert (tmp_path / "clean_sequences.fa").read_text() == "existing FASTA\n"


def test_main_all_foreign_writes_empty_fasta(tmp_path, monkeypatch):
    prepare_main(tmp_path, monkeypatch, ">foreign\nACGT\n", [("foreign", 4, "foreign")])
    contamination.main()
    assert (tmp_path / "clean_sequences.fa").read_bytes() == b""
    lineage = pandas.read_csv(tmp_path / "lineage_compatibility.tsv", sep="\t")
    assert lineage["is_compatible_lineage"].tolist() == [False]
