"""Default stage figure and Swiss-Prot connections retain their data meaning."""
import argparse
from pathlib import Path

import pytest

from workflow.tests.test_gene_model_refinement_busco import busco, summary


def options(**overrides):
    value = dict(no_rescue_swissprot=False, swissprot_db_prefix=None,
                 swissprot_metadata=None, swissprot_cache=None, cpus=2, swissprot_memory_gb=8)
    return argparse.Namespace(**(value | overrides))


def test_default_audit_discovery_verifies_the_exact_rescue_and_rejects_partial_audits(tmp_path, monkeypatch):
    anchor = tmp_path / "rescue"
    changes = {"rescue_reference_selection": {"rescue_output": str(anchor)}}
    called = []
    monkeypatch.setattr(busco, "collect_rescue_swissprot_evidence", lambda c, d: called.append(d) or c)
    with pytest.raises(ValueError, match="audit is required"):
        busco.connect_rescue_swissprot(changes, report=tmp_path / "report", options=options())
    assert busco.connect_rescue_swissprot(changes, options=options(no_rescue_swissprot=True)) is changes
    assert not called
    audit = Path(str(anchor) + ".swissprot")
    audit.mkdir()
    with pytest.raises(ValueError, match="Incomplete Swiss-Prot"):
        busco.connect_rescue_swissprot(changes, options=options())
    (audit / "receipt.json").write_text("{}")
    assert busco.connect_rescue_swissprot(changes, options=options()) is changes
    assert called == [audit]
    monkeypatch.setattr(busco, "collect_rescue_swissprot_evidence", lambda *a: (_ for _ in ()).throw(ValueError("different rescue inputs")))
    with pytest.raises(ValueError, match="different rescue inputs"):
        busco.connect_rescue_swissprot(changes, options=options())


def test_missing_audit_search_receives_shared_db_cache_and_resource_budget(tmp_path, monkeypatch):
    import rescue_swissprot_evidence as swiss
    changes = {"rescue_reference_selection": {"rescue_output": str(tmp_path / "rescue")}}
    requests, collections = [], []
    monkeypatch.setattr(swiss, "audit", lambda request: requests.append(request))
    monkeypatch.setattr(busco, "collect_rescue_swissprot_evidence", lambda c, d: collections.append(d) or c)
    report = tmp_path / "report"
    args = options(swissprot_db_prefix=tmp_path / "sprot", swissprot_metadata=tmp_path / "metadata",
                   swissprot_cache=tmp_path / "cache", cpus=3, swissprot_memory_gb=5)
    busco.connect_rescue_swissprot(changes, report=report, options=args)
    assert len(requests) == 1
    request = requests[0]
    assert (request.cpus, request.memory_gb, request.cache, request.db_prefix) == (3, 5, args.swissprot_cache, args.swissprot_db_prefix)
    assert request.rescue_output == tmp_path / "rescue"
    assert request.evalue == swiss.DEFAULTS["evalue"]
    assert request.minimum_alignment == swiss.DEFAULTS["minimum_alignment"]
    assert request.query_coverage == swiss.DEFAULTS["query_coverage"]
    assert request.search_evalue is None
    assert collections == [report / "rescue_swissprot"]


def test_search_paths_require_an_annotation_binding_before_redraw(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr("sys.argv", ["busco", "--plot-only", "--report", str(tmp_path),
                                    "--swissprot-db-prefix", "/db/sprot", "--swissprot-metadata", "/db/metadata",
                                    "--swissprot-cache", "/cache"])
    with pytest.raises(SystemExit) as error:
        busco.main()
    assert error.value.code == 2
    assert "requires --output" in capsys.readouterr().err


def test_default_figure_reorders_equal_panels_and_preserves_signed_counts_and_unavailable_rows(tmp_path, monkeypatch):
    from matplotlib.figure import Figure
    from matplotlib.text import Text
    results = []
    for i, complete in enumerate((7, 9, 8)):
        path = tmp_path / f"result{i}"
        path.write_text(summary(complete))
        results.append(busco.read_result(path))
    rows = [busco.staged_result({"species": name, "refinement_status": status}, results[1], results[2], results[0])
            for name, status in (("Ancistrocladus_abbreviatus", "analysed"), ("Drosophyllum_lusitanicum", "not_analysed"))]
    changes = {"species": {
        rows[0]["species"]: {"refinement_status": "analysed", "prior_rescued_loci": 3,
                             "accepted_repair_paths": 1, "accepted_isoform_paths": 2,
                             "rescue_repeat_groups": dict(zip(busco.REPEAT_GROUPS, (3, 0, 0, 0), strict=True))},
        rows[1]["species"]: {"refinement_status": "not_analysed", "prior_rescued_loci": None,
                             "accepted_repair_paths": None, "accepted_isoform_paths": None}}}
    original = Figure.savefig
    saved = []
    def save(fig, path, **kwargs):
        fig.canvas.draw()
        assert fig.get_size_inches()[0] == pytest.approx(7.2)
        assert [a.get_title().replace("\n", " ") for a in fig.axes] == [
            "Missing-gene rescue", "Gene-model refinement", "BUSCO: before / after rescue / after refinement", "Change in complete groups"]
        assert len({round(a.get_position().width, 8) for a in fig.axes}) == 1
        assert all(t.get_fontsize() == 8 for t in fig.findobj(match=Text))
        labels = fig.axes[0].get_yticklabels()
        assert labels[0].get_text() == "Ancistrocladus\nabbreviatus"
        assert all(t.get_fontstyle() == "italic" for t in labels)
        busco_axis, delta_axis = fig.axes[2:]
        assert len(busco_axis.patches) == 17  # one shaded row + 12 analysed + 4 pre-rescue
        assert not any(t.get_text() == "Not analysed" for t in busco_axis.texts)
        bars = [p for p in delta_axis.patches if p.get_height() == pytest.approx(.65)]
        assert [(p.get_x(), p.get_width()) for p in bars] == [(0, 20), (0, -10)]
        assert delta_axis.get_xlim()[0] < -10 and delta_axis.get_xlim()[1] > 20
        assert all(len(a.child_axes) == 1 for a in (busco_axis, delta_axis))
        for axis in (busco_axis, delta_axis):
            forward, inverse = axis.child_axes[0]._functions
            assert forward(20) == pytest.approx(2)
            assert inverse(2) == pytest.approx(20)
        saved.append(Path(path).suffix)
        return original(fig, path, **kwargs)
    monkeypatch.setattr(Figure, "savefig", save)
    busco.plot_comparison(rows, tmp_path, changes)
    assert set(saved) == {".png", ".svg", ".pdf"}
    svg = (tmp_path / "busco_comparison.svg").read_text()
    assert "Swiss-Prot audit not supplied" in svg
    assert "Repeat annotation" not in svg
    assert (tmp_path / "busco_comparison.pdf").read_bytes().startswith(b"%PDF")
