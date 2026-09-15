import csv
import sys
from pathlib import Path
from types import SimpleNamespace

SUPPORT_DIR = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT_DIR))

import add_hgt_taxonomy_columns as migration  # noqa: E402


def _write_table(path: Path, fieldnames, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _read_table(path: Path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class _FakeNcbi:
    def __init__(self):
        self.lineage_ids = [1, 2, 3, 4, 5, 6, 7, 8]
        self.ranks = {
            1: "no rank",
            2: "superkingdom",
            3: "kingdom",
            4: "phylum",
            5: "subphylum",
            6: "class",
            7: "species",
            8: "clade",
        }
        self.names = {
            1: "root",
            2: "Eukaryota",
            3: "Viridiplantae",
            4: "Streptophyta",
            5: "Embryophyta",
            6: "Magnoliopsida",
            7: "Example species",
            8: "Example clade",
        }

    def get_name_translator(self, names):
        return {name: [7] for name in names if name == "Example species"}

    def get_lineage(self, taxid):
        return self.lineage_ids if int(taxid) == 7 else []

    def get_rank(self, taxids):
        return {int(taxid): self.ranks.get(int(taxid), "") for taxid in taxids}

    def get_taxid_translator(self, taxids):
        return {int(taxid): self.names[int(taxid)] for taxid in taxids if int(taxid) in self.names}


def test_resolve_taxonomy_annotation_keeps_standard_and_unusual_ranks():
    resolver = migration.TaxonomyResolver("")
    resolver.enabled = True
    resolver.ncbi = _FakeNcbi()

    annotation = migration.resolve_taxonomy_annotation("Example species", "", resolver)

    assert annotation["domain"] == "Eukaryota"
    assert annotation["superkingdom"] == "Eukaryota"
    assert annotation["kingdom"] == "Viridiplantae"
    assert annotation["phylum"] == "Streptophyta"
    assert annotation["subphylum"] == "Embryophyta"
    assert annotation["class"] == "Magnoliopsida"
    assert annotation["species"] == "Example species"
    assert "clade:Example clade" in annotation["taxonomy"]


def test_migration_appends_all_taxonomy_columns_and_preserves_existing_values(tmp_path, monkeypatch):
    branch_path = tmp_path / "hgt_branch_candidates.tsv"
    gene_path = tmp_path / "hgt_gene_candidates.tsv"
    orthogroup_path = tmp_path / "hgt_orthogroup_summary.tsv"

    _write_table(
        gene_path,
        [
            "orthogroup",
            "gene_id",
            "gene_taxon",
            "besthit_accession",
            "besthit_organism",
            "besthit_taxid",
            "old_score",
        ],
        [
            {
                "orthogroup": "OG0001",
                "gene_id": "geneA",
                "gene_taxon": "Arabidopsis thaliana",
                "besthit_accession": "P1",
                "besthit_organism": "Escherichia coli",
                "besthit_taxid": "562.0",
                "old_score": "0.75",
            },
            {
                "orthogroup": "OG0001",
                "gene_id": "geneB",
                "gene_taxon": "Arabidopsis thaliana",
                "besthit_accession": "P2",
                "besthit_organism": "Bacillus subtilis",
                "besthit_taxid": "1423.0",
                "old_score": "0.25",
            },
        ],
    )
    _write_table(
        branch_path,
        ["orthogroup", "branch_id", "candidate_genes", "old_branch_value"],
        [{"orthogroup": "OG0001", "branch_id": "3", "candidate_genes": "geneB; geneA", "old_branch_value": "keep"}],
    )
    _write_table(
        orthogroup_path,
        ["orthogroup", "hgt_branch_count", "old_og_value"],
        [{"orthogroup": "OG0001", "hgt_branch_count": "1", "old_og_value": "keep"}],
    )

    monkeypatch.setattr(migration, "TaxonomyResolver", lambda _dbfile: SimpleNamespace(enabled=True))

    def make_annotation(**values):
        annotation = {rank: "" for rank in migration.TAXONOMIC_RANKS}
        annotation["taxonomy"] = values.pop("taxonomy", "")
        values["class"] = values.pop("class_rank", "")
        annotation.update(values)
        return annotation

    def fake_resolve(name, taxid, _resolver):
        if "Arabidopsis" in str(name):
            return make_annotation(
                domain="Eukaryota",
                superkingdom="Eukaryota",
                kingdom="Viridiplantae",
                phylum="Streptophyta",
                class_rank="Magnoliopsida",
                order="Brassicales",
                family="Brassicaceae",
                genus="Arabidopsis",
                species="Arabidopsis thaliana",
                taxonomy="superkingdom:Eukaryota; kingdom:Viridiplantae; phylum:Streptophyta; class:Magnoliopsida; order:Brassicales; family:Brassicaceae; genus:Arabidopsis; species:Arabidopsis thaliana",
            )
        if "Escherichia" in str(name) or str(taxid).startswith("562"):
            return make_annotation(
                domain="Bacteria",
                superkingdom="Bacteria",
                phylum="Pseudomonadota",
                class_rank="Gammaproteobacteria",
                order="Enterobacterales",
                family="Enterobacteriaceae",
                genus="Escherichia",
                species="Escherichia coli",
                taxonomy="superkingdom:Bacteria; phylum:Pseudomonadota; class:Gammaproteobacteria; order:Enterobacterales; family:Enterobacteriaceae; genus:Escherichia; species:Escherichia coli",
            )
        if "Bacillus" in str(name) or str(taxid).startswith("1423"):
            return make_annotation(
                domain="Bacteria",
                superkingdom="Bacteria",
                phylum="Bacillota",
                class_rank="Bacilli",
                order="Bacillales",
                family="Bacillaceae",
                genus="Bacillus",
                species="Bacillus subtilis",
                taxonomy="superkingdom:Bacteria; phylum:Bacillota; class:Bacilli; order:Bacillales; family:Bacillaceae; genus:Bacillus; species:Bacillus subtilis",
            )
        return make_annotation()

    monkeypatch.setattr(migration, "resolve_taxonomy_annotation", fake_resolve)

    counts = migration.enrich_hgt_tables(branch_path, gene_path, orthogroup_path, "unused")

    assert counts == {"branch": 1, "gene": 2, "orthogroup": 1}
    genes = _read_table(gene_path)
    assert genes[0]["old_score"] == "0.75"
    assert genes[0]["recipient_phylum"] == "Streptophyta"
    assert genes[0]["donor_phylum"] == "Pseudomonadota"
    assert genes[1]["donor_phylum"] == "Bacillota"
    assert genes[0]["recipient_class"] == "Magnoliopsida"
    assert genes[0]["donor_family"] == "Enterobacteriaceae"
    assert "species:Arabidopsis thaliana" in genes[0]["recipient_taxonomy"]
    assert {"recipient_domain", "recipient_species", "donor_class", "donor_taxonomy"}.issubset(genes[0])

    branches = _read_table(branch_path)
    assert branches[0]["old_branch_value"] == "keep"
    assert branches[0]["recipient_phyla"] == "Streptophyta"
    assert branches[0]["donor_phyla"] == "Bacillota; Pseudomonadota"
    assert branches[0]["recipient_classes"] == "Magnoliopsida"
    assert branches[0]["donor_families"] == "Bacillaceae; Enterobacteriaceae"
    assert "class:Magnoliopsida" in branches[0]["recipient_taxonomies"]
    assert branches[0]["representative_gene_id"] == "geneA"
    assert branches[0]["representative_besthit_accession"] == "P1"
    assert branches[0]["representative_besthit_organism"] == "Escherichia coli"
    assert branches[0]["representative_recipient_phylum"] == "Streptophyta"
    assert branches[0]["representative_donor_phylum"] == "Pseudomonadota"
    assert branches[0]["representative_recipient_class"] == "Magnoliopsida"
    assert branches[0]["representative_donor_family"] == "Enterobacteriaceae"

    orthogroups = _read_table(orthogroup_path)
    assert orthogroups[0]["old_og_value"] == "keep"
    assert orthogroups[0]["recipient_phyla"] == "Streptophyta"
    assert orthogroups[0]["donor_phyla"] == "Pseudomonadota; Bacillota"
    assert orthogroups[0]["recipient_orders"] == "Brassicales"
    assert orthogroups[0]["donor_families"] == "Enterobacteriaceae; Bacillaceae"
    assert "species:Escherichia coli" in orthogroups[0]["donor_taxonomies"]
