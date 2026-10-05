import importlib.util
import sys
from pathlib import Path

import pytest

from workflow.support.representative_selection import RepresentativeMap, matching_identifier


def write_map(path, rows):
    path.write_text("species\tgene_id\tcandidate_id\tsource_transcript_id\tstatus\tscore\tmargin\treason\n"
                    + "".join("\t".join(row) + "\n" for row in rows))
    return path


def choice(gene="Plant_species_gene1", transcript="rna-tx1", candidate="c1", species="Plant_species"):
    return species, gene, candidate, transcript, "selected", "1.0", "0.2", "conserved"


def test_map_resolves_only_unique_structural_and_declared_species_aliases(tmp_path):
    path = write_map(tmp_path / "map.tsv", [choice()])
    selection = RepresentativeMap(path)
    assert selection.choice("gene1", "Plant_species")["source_transcript_id"] == "rna-tx1"
    assert selection.choice("Plant_species_gene1")["candidate_id"] == "c1"
    assert matching_identifier("tx1", ["rna-tx1", "tx2"], "Plant_species") == "rna-tx1"
    assert matching_identifier("tx1", ["tx1", "rna-tx1"], "Plant_species") == "tx1"
    with pytest.raises(ValueError, match="missing"):
        selection.choice("gene2", "Plant_species")
    with pytest.raises(ValueError, match="Ambiguous"):
        matching_identifier("Plant_species_tx1", ["tx1", "rna-tx1"], "Plant_species")


@pytest.mark.parametrize("rows,message", [
    ([choice(), choice()], "Duplicate"),
    ([choice(), choice("gene2", candidate="c2")], "ownership"),
    ([choice(), choice("gene2", transcript="tx2")], "ownership"),
    ([choice(transcript="")], "Invalid"),
    ([choice()[:5] + ("nan", "0.2", "conserved")], "Nonfinite"),
])
def test_map_rejects_duplicates_missing_ids_and_nonfinite_scores(tmp_path, rows, message):
    with pytest.raises(ValueError, match=message):
        RepresentativeMap(write_map(tmp_path / "map.tsv", rows))


def test_map_rejects_ambiguous_gene_aliases_without_erasing_exact_identity(tmp_path):
    path = write_map(tmp_path / "map.tsv", [choice("gene1"), choice("gene-gene1", "tx2", "c2")])
    selection = RepresentativeMap(path)
    assert selection.choice("gene1", "Plant_species")["candidate_id"] == "c1"
    with pytest.raises(ValueError, match="ambiguous"):
        selection.choice("Plant_species_gene1")


def test_unprefixed_gene_with_underscores_is_not_mistaken_for_species(tmp_path):
    selection = RepresentativeMap(write_map(tmp_path / "map.tsv", [choice("LOCUS_gene1")]))
    assert selection.choice("LOCUS_gene1")["candidate_id"] == "c1"
    selection = RepresentativeMap(write_map(tmp_path / "map.tsv", [choice("LOCUS_gene1"),
                                       choice("LOCUS_gene1", "tx2", "c2", "Other_species")]))
    with pytest.raises(ValueError, match="ambiguous"):
        selection.choice("LOCUS_gene1")


def test_validated_map_can_cross_script_and_package_import_boundaries(tmp_path, monkeypatch):
    source = Path(__file__).parents[1] / 'support/representative_selection.py'
    spec = importlib.util.spec_from_file_location('representative_selection_alternate', source)
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, spec.name, module)
    spec.loader.exec_module(module)
    path = write_map(tmp_path / 'map.tsv', [choice()])
    selection = RepresentativeMap(path)
    path.unlink()
    assert module.load_representative_map(selection) is selection
    assert selection.choice('gene1', 'Plant_species')['source_transcript_id'] == 'rna-tx1'
