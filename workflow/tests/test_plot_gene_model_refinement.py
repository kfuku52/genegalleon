"""Review figures use verified publications and distinguish paths from choices."""

import hashlib
import json
from importlib import import_module
from xml.etree import ElementTree

import pytest

from workflow.tests.test_gene_model_refinement import refinement, tiny_inputs

review = import_module("plot_gene_model_refinement")


def completed(tmp_path):
    inputs, edges, _ = tiny_inputs(tmp_path)
    root = tmp_path / "refinement"
    value = refinement.plan(root, inputs=inputs, edges=edges, mode="off")
    refinement.finalize(root, value)
    return root


def test_review_counts_complete_publication_and_renders_empty_change_summary(tmp_path):
    root = completed(tmp_path)
    data = review.collect(root, max_loci=2)
    assert len(data["species"]) == 3
    assert sum(s["source_candidates"] for s in data["species"].values()) == 6
    assert sum(s["changed_representatives"] for s in data["species"].values()) == 0
    assert data["details"] == []
    out = tmp_path / "report"
    out.mkdir()
    review.plot_summary(data, out)
    review.write_review(data, out)
    assert (out / "summary.png").read_bytes().startswith(b"\x89PNG")
    ElementTree.parse(out / "summary.svg")
    assert "data:image/png;base64," in (out / "review.html").read_text()
    # The source receipt is not replaced or extended by rendering.
    receipt = json.loads((root / "effective" / "receipt.json").read_text())
    assert all(not name.endswith(("png", "svg", "html")) for name in receipt["files"])


def test_review_rejects_prediction_content_tampering(tmp_path):
    root = completed(tmp_path)
    predictions = root / "predictions" / "Species_target" / "predictions.json"
    predictions.write_text('[{"status":"accepted"}]\n')
    with pytest.raises(ValueError, match="Review artifact changed"):
        review.collect(root)


def test_review_counts_interleaved_duplicates_and_uses_last_gallery_selection(tmp_path):
    root = completed(tmp_path)
    changes_path = root / "effective" / "changes.json"
    records = {row["species"]: row for row in json.loads(changes_path.read_text())}
    first = {**records["Species_target"], "status": "conserved", "reason": "earlier_selection"}
    last = {**records["Species_target"], "status": "sequence_repaired", "reason": "last_selection"}
    changes = [first, {**records["Species_donor2"], "status": "conserved"},
               records["Species_donor1"], last]
    changes_path.write_text(json.dumps(changes) + "\n")
    # Re-sign the changed fixture so collection exercises the full publication
    # proof instead of bypassing verification to inspect aggregation behavior.
    receipt_path = root / "effective" / "receipt.json"
    receipt = json.loads(receipt_path.read_text())
    receipt["files"]["changes.json"] = hashlib.sha256(changes_path.read_bytes()).hexdigest()
    receipt_path.write_text(json.dumps(receipt) + "\n")

    data = review.collect(root, max_loci=1, preferred_species="Species_target")

    target = data["species"]["Species_target"]
    assert target["changed_representatives"] == 2
    assert target["selection_status"] == {"conserved": 1, "sequence_repaired": 1}
    assert data["species"]["Species_donor2"]["changed_representatives"] == 1
    assert data["species"]["Species_donor1"]["changed_representatives"] == 0
    assert data["changed_loci_available"] == 2
    assert len(data["details"]) == 1
    assert data["details"][0]["species"] == "Species_target"
    assert data["details"][0]["selection"] == last


def test_locus_plot_escapes_labels_and_keeps_genomic_spacing_and_orientation():
    locus = {
        "species": "S", "gene_id": "g", "seqid": "chr&1", "strand": "-",
        "baseline_id": "S_t<1",
        "selection": {"candidate_id": "S_t<1"},
        "candidates": [{"candidate_id": "S_t<1", "blocks": [[10, 20, 2], [90, 100, 0]],
                        "origin": "original", "cds_length": 20}],
        "predictions": [],
    }
    svg = review.locus_svg(locus)
    node = ElementTree.fromstring(svg)
    assert node.attrib["aria-label"] == "Coding exon structure"
    assert "t&lt;1" in svg and "chr&amp;1:11–100 (-)" in svg
    assert 'width="72.22"' in svg
