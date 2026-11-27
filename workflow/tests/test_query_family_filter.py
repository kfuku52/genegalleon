import csv
from types import SimpleNamespace

import pytest

from workflow.support import query_family_filter as filtering


@pytest.fixture
def references(tmp_path):
    files = SimpleNamespace(**{key: tmp_path / key for key in ("entries", "tree", "mapping", "metadata")})
    files.entries.write_text(
        "ENTRY_AC\tENTRY_TYPE\nIPR1\tFamily\nIPR2\tFamily\nIPR3\tFamily\n"
        "IPR4\tFamily\nIPR5\tDomain\nIPR6\tFamily\n")
    files.tree.write_text(
        "IPR2::target::\n--IPR4::child::\n"
        "IPR1::parent::\n--IPR2::target::\n----IPR4::child::\n--IPR3::other::\n")
    files.mapping.write_text(
        "query\tinterpro_family_id\tbroader_interpro_family_ids\nFAMILY\tIPR2\tIPR6\n")
    files.metadata.write_text(
        "accession\tinterpro_ids\nTARGET\tIPR4;IPR1;IPR5\nOTHER\tIPR3;IPR1;IPR5\n"
        "BROAD\tIPR1;IPR5;IPR6\n")
    return files


def test_family_filter_prioritizes_supported_hits_before_hit_cap(tmp_path, references):
    genes = ("outside", "target", "conflict", "unknown", "fusion", "beyond_cap")
    args = SimpleNamespace(**vars(references), min_coverage=0, cap=4,
                           **{key: tmp_path / key for key in ("discovery", "genes", "selected", "audit")})
    args.genes.write_text("\n".join(genes) + "\n")
    args.discovery.write_text(
        "sacc\tqjointcov\tsstart\tsend\n"
        + "".join(f"{gene}\t0.05\t10;300\t109;399\n" for gene in genes)
    )
    hits = tmp_path / "hits.tsv"
    hits.write_text(
        "outside\tsp|OTHER|name\t100\t1e-30\t10\t109\n"
        "target\tTARGET\t100\t1e-30\t10\t109\n"
        "target\tOTHER\t89\t1e-25\t10\t109\n"
        "conflict\tTARGET\t100\t1e-30\t10\t109\n"
        "conflict\tOTHER\t90\t1e-24\t30\t129\n"
        "unknown\tBROAD\t100\t1e-30\t10\t109\n"
        "unknown\tOTHER\t101\t1e-30\t170\t200\n"
        "fusion\tTARGET\t100\t1e-27\t10\t109\n"
        "fusion\tOTHER\t150\t1e-37\t300\t399\n"
        "beyond_cap\tTARGET\t100\t1e-28\t10\t109\n"
    )
    context = filtering.read_family_context(args.entries, args.tree, args.mapping, "FAMILY")
    assert context == ("IPR2", {"IPR2", "IPR4"}, {"IPR3"})
    filtering.screen(args, filtering.ranked_genes(args.genes), hits, context)
    assert args.selected.read_text().splitlines() == ["target", "conflict", "fusion", "beyond_cap"]
    with args.audit.open() as handle:
        rows = {row["gene_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    assert [rows[gene]["status"] for gene in genes] == [
        "other_family_supported", "target_supported", "competing_families",
        "unresolved", "target_supported", "target_supported",
    ]
    assert rows["outside"]["other_accession"] == "OTHER"
    assert rows["target"]["target_family_ids"] == "IPR4"
    assert rows["fusion"]["other_accession"] == ""
    assert rows["unknown"]["retained"] == "0"
    assert rows["beyond_cap"]["retained"] == "1"


def test_mapped_accession_must_be_family_type(references):
    references.mapping.write_text("query\tinterpro_family_id\nFAMILY\tIPR5\n")
    with pytest.raises(ValueError, match="InterPro Family ID"):
        filtering.read_family_context(references.entries, references.tree, references.mapping, "FAMILY")


def test_unmapped_queries_keep_ranked_candidates_without_search(tmp_path, references, monkeypatch):
    references.mapping.write_text("query\tinterpro_family_id\nFAMILY\tunmapped\n")
    genes = tmp_path / "genes"
    genes.write_text("first\nsecond\n")
    audit, selected = tmp_path / "audit", tmp_path / "selected"
    monkeypatch.setattr(filtering.sys, "argv", [
        "query_family_filter.py", "--query", "FAMILY", "--min-coverage", "0", "--cap", "1",
        "--threads", "1", "--entries", str(references.entries), "--tree", str(references.tree),
        "--mapping", str(references.mapping), "--metadata", "/absent", "--genes", str(genes),
        "--discovery", "/absent", "--sequence-store", "/absent", "--swissprot-db", "/absent",
        "--audit", str(audit), "--selected", str(selected),
    ])
    filtering.main()
    assert selected.read_text() == "first\n"
    with audit.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [(row["target_family"], row["status"], row["retained"]) for row in rows] == [
        ("", "unmapped_query", "1"), ("", "unmapped_query", "0"),
    ]
