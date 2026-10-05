"""Exercise canonical selected IDs through the real nucleotide reader."""

import hashlib
import json
import subprocess
from pathlib import Path
from types import SimpleNamespace
from urllib.parse import quote

import pytest
from Bio.Seq import Seq
from kffractbias.jcvi import prepare_genome

from workflow.support import fractionation_selected_inputs as adapter
from workflow.support.fractionation_selected_inputs import prepare_selected_inputs
from workflow.support.representative_selection import REQUIRED_COLUMNS
from workflow.tests.test_gene_model_refinement_downstream import selected_view, tsv


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bundle(tmp_path, kind="normal", source_transcript="short"):
    rows, choices, original = [], [], {}
    for species, code in [("Target_species", 4), ("Query_species", 1)]:
        cds = "ATGTGAAAATAA" if code == 4 else "ATGTGGAAATAA"
        start, end, phase = 11, 22, 0
        if kind == "partial":
            cds, end, phase = "GA" + cds, 24, 2
        tid = source_transcript
        escaped = quote(tid, safe="._-:")
        gff = (f"chr1\ts\tgene\t1\t40\t.\t+\t.\tID=g\n"
               f"chr1\ts\tmRNA\t{start}\t{end}\t.\t+\t.\tID={escaped};Parent=g\n"
               f"chr1\ts\tCDS\t{start}\t{end}\t.\t+\t{phase}\tParent={escaped}\n"
               "chr1\ts\tmRNA\t1\t40\t.\t+\t.\tID=long;Parent=g\n"
               "chr1\ts\tCDS\t1\t39\t.\t+\t0\tParent=long\n")
        if kind == "gtf":
            gff = (f'chr1\ts\tCDS\t{start}\t{end}\t.\t+\t{phase}\tgene_id "g"; transcript_id "{tid}";\n'
                   'chr1\ts\tCDS\t1\t39\t.\t+\t0\tgene_id "g"; transcript_id "long";\n')
        if kind == "direct":
            gff = ("chr1\ts\tgene\t1\t40\t.\t+\t.\tID=g\n"
                   f"chr1\ts\tCDS\t{start}\t{end}\t.\t+\t{phase}\tID={escaped};Parent=g\n"
                   "chr1\ts\tCDS\t1\t39\t.\t+\t0\tID=long;Parent=g\n")
        if kind == "parentless":
            tid = "g"
            gff = f"chr1\ts\tCDS\t{start}\t{end}\t.\t+\t{phase}\tID=g\n"
        genome = "C" * 10 + cds + "C" * 40
        if kind == "minus_spliced":
            gff = ("chr1\ts\tgene\t1\t40\t.\t-\t.\tID=g\n"
                   f"chr1\ts\tmRNA\t11\t35\t.\t-\t.\tID={escaped};Parent=g\n"
                   f"chr1\ts\tCDS\t30\t35\t.\t-\t0\tParent={escaped}\n"
                   f"chr1\ts\tCDS\t11\t16\t.\t-\t0\tParent={escaped}\n")
            genome = "C" * 10 + str(Seq(cds[6:]).reverse_complement()) + "C" * 13 + str(Seq(cds[:6]).reverse_complement()) + "C" * 20
        paths = {role: tmp_path / (species + "." + role + ".fa") for role in ("cds", "protein", "genome")}
        paths["gff"] = tmp_path / (species + ".gff3")
        paths["cds"].write_text(">" + species + "_g\n" + cds + "\n")
        paths["protein"].write_text(">" + species + "_g\n" + str(Seq(cds[phase:]).translate(table=code)).removesuffix("*") + "\n")
        paths["genome"].write_text(">chr1\n" + genome + "\n")
        paths["gff"].write_text(gff)
        original[species] = digest(paths["gff"])
        rows.append({"species": species, "genetic_code": code, **paths})
        choices.append({"species": species, "gene_id": species + "_g", "candidate_id": species + "_candidate",
                        "source_transcript_id": tid, "status": "unchanged", "score": "", "margin": "", "reason": "source"})
    mapping = tmp_path / "choices.tsv"
    tsv(mapping, REQUIRED_COLUMNS, choices)
    for row in rows:
        row["representative_map"] = mapping
        for role in ("cds", "protein", "genome", "gff", "representative_map"):
            row[role + "_sha256"] = digest(row[role])
    manifest = tmp_path / "selected.tsv"
    tsv(manifest, list(rows[0]), rows)
    return manifest, original, rows


@pytest.mark.parametrize("kind,transcript", [
    ("normal", "short"), ("partial", "short"), ("direct", "short"), ("parentless", "g"), ("minus_spliced", "short"),
    ("normal", "t,1"), ("normal", "t%2C1"), ("normal", "t;1"),
    ("gtf", "t%2C1"), ("gtf", "t,1"),
])
def test_exact_selected_paths_map_through_real_kffractbias_and_official_cli(tmp_path, kind, transcript):
    manifest, originals, sources = bundle(tmp_path, kind, transcript)
    result = prepare_selected_inputs(manifest, ["Target_species", "Query_species"], tmp_path / "cache")
    for selected, source in zip(result, sources, strict=True):
        prepared = prepare_genome(selected["species"], selected["cds"], selected["gff"], tmp_path,
                                  feature="gene", attribute="ID")
        assert len(prepared.mapping.genes) == 1
        gene = prepared.mapping.genes[0]
        assert gene.gene_id == selected["species"] + "_g"
        assert gene.start == 10 and gene.end == (35 if kind == "minus_spliced" else 24 if kind == "partial" else 22)
        assert gene.strand == ("-" if kind == "minus_spliced" else "+")
        assert prepared.cds_path.read_text() == Path(source["cds"]).read_text()
        receipt = json.loads(Path(selected["receipt"]).read_text())
        assert receipt["mapping"][0]["source_transcript_id"] == transcript
        assert receipt["identity"]["inputs"]["gff"]["sha256"] == originals[selected["species"]]
        assert "gff2genestat.py" in receipt["identity"]["implementation"]
        assert "representative_selection.py" in receipt["identity"]["implementation"]
        assert "format_species_common.py" in receipt["identity"]["implementation"]
        assert "gff_attribute_syntax.py" in receipt["identity"]["implementation"]
        assert set(receipt["identity"]["runtime"]) == {"python", "numpy", "pandas"}
        assert digest(Path(source["gff"])) == originals[selected["species"]]
    target, query = result
    validation = subprocess.run(["kffractbias", "validate", "--target-cds", target["cds"], "--target-gff", target["gff"],
                                 "--query-cds", query["cds"], "--query-gff", query["gff"],
                                 "--target-feature", "gene", "--target-attribute", "ID",
                                 "--query-feature", "gene", "--query-attribute", "ID", "--pairwise"],
                                text=True, capture_output=True)
    assert validation.returncode == 0, validation.stderr
    assert json.loads(validation.stdout)["target"]["matched_gene_count"] == 1


def test_adapter_is_immutable_and_manifest_bound(tmp_path):
    manifest, _, _ = bundle(tmp_path)
    first = prepare_selected_inputs(manifest, ["Target_species"], tmp_path / "cache")[0]
    stamp = Path(first["receipt"]).stat().st_mtime_ns
    assert prepare_selected_inputs(manifest, ["Target_species"], tmp_path / "cache")[0] == first
    assert Path(first["receipt"]).stat().st_mtime_ns == stamp
    assert digest(manifest) in first["gff"]
    Path(first["gff"]).write_text(Path(first["gff"]).read_text() + "# tamper\n")
    with pytest.raises(ValueError, match="Changed fractionation"):
        prepare_selected_inputs(manifest, ["Target_species"], tmp_path / "cache")


def test_adapter_handles_a_real_published_refinement_view(tmp_path):
    effective, _, _ = selected_view.__wrapped__(tmp_path, SimpleNamespace(param={}))
    before = digest(effective / "species_gff/Target_species.gff3")
    rows = prepare_selected_inputs(effective / "inputs.tsv", ["Target_species", "Donor_one"], tmp_path / "reader_cache")
    for row in rows:
        prepared = prepare_genome(row["species"], row["cds"], row["gff"], tmp_path, feature="gene", attribute="ID")
        assert prepared.mapping.matched_gene_count == 1
    target = prepare_genome("check", rows[0]["cds"], rows[0]["gff"], tmp_path, feature="gene", attribute="ID")
    assert target.mapping.genes[0].start == 105
    assert json.loads(Path(rows[0]["receipt"]).read_text())["mapping"][0]["source_transcript_id"] == "short"
    assert digest(effective / "species_gff/Target_species.gff3") == before


@pytest.mark.parametrize("failure", ["missing_transcript", "length_mismatch"])
def test_adapter_rejects_inconsistent_external_bundle_without_longest_fallback(tmp_path, failure):
    manifest, _, rows = bundle(tmp_path)
    if failure == "missing_transcript":
        mapping = rows[0]["representative_map"]
        mapping.write_text(mapping.read_text().replace("\tshort\t", "\tabsent\t"))
        for row in rows:
            row["representative_map_sha256"] = digest(mapping)
    else:
        source = rows[0]["cds"]
        source.write_text(source.read_text().rstrip() + "AAA\n")
        rows[0]["cds_sha256"] = digest(source)
    tsv(manifest, list(rows[0]), rows)
    with pytest.raises(ValueError, match="absent|differs"):
        prepare_selected_inputs(manifest, ["Target_species"], tmp_path / "cache")
    assert not list((tmp_path / "cache").rglob("receipt.json"))


def test_adapter_detects_source_change_before_cached_view_delivery(tmp_path, monkeypatch):
    manifest, _, rows = bundle(tmp_path)
    first = prepare_selected_inputs(manifest, ["Target_species"], tmp_path / "cache")[0]
    real_digest = adapter.file_digest

    def change_after_cache_check(path):
        value = real_digest(path)
        if Path(path) == Path(first["gff"]):
            rows[0]["cds"].write_text(rows[0]["cds"].read_text() + "A\n")
        return value

    monkeypatch.setattr(adapter, "file_digest", change_after_cache_check)
    with pytest.raises(ValueError, match="changed before delivery"):
        prepare_selected_inputs(manifest, ["Target_species"], tmp_path / "cache")


def test_adapter_runtime_identity_changes_namespace(tmp_path, monkeypatch):
    manifest, _, _ = bundle(tmp_path)
    first = prepare_selected_inputs(manifest, ["Target_species"], tmp_path / "cache")[0]
    identity = adapter.runtime_identity()
    monkeypatch.setattr(adapter, "runtime_identity", lambda: {**identity, "pandas": "different-runtime"})
    second = prepare_selected_inputs(manifest, ["Target_species"], tmp_path / "cache")[0]
    assert first["gff"] != second["gff"]
    assert Path(first["gff"]).read_bytes() == Path(second["gff"]).read_bytes()
