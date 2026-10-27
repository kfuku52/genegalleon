import json
import subprocess
import sys
from pathlib import Path

import pytest
from PIL import Image

from workflow.support import synteny_multitrack as multi
from workflow.tests.test_pairwise_synteny_runtime import run_core, write_genome

REPO = Path(__file__).resolve().parents[2]


@pytest.fixture
def saved_panel(tmp_path):
    tracks = []
    for species, prefix in (("Anc_species", "a"), ("Dio_species", "d"), ("Tri_species", "t")):
        bed = tmp_path / (species + ".bed")
        bed.write_text("".join(f"chr1\t{i * 10}\t{i * 10 + 5}\t{prefix}{i}\n" for i in range(8)))
        tracks.append({"species": species, "bed": str(bed)})
    tracks[0]["fixed_orientation"] = {"chr1": "+"}
    tracks[2]["fixed_orientation"] = {"chr1": "+"}
    ab, bc = tmp_path / "ad.simple", tmp_path / "dt.simple"
    ab.write_text("a0 a7 d0 d7 8 -\n")
    bc.write_text("d0 d1 t0 t1 2 +\n")
    return {"tracks": tracks, "edges": [{"tracks": [0, 1], "blocks": str(ab)}, {"tracks": [1, 2], "blocks": str(bc)}]}


@pytest.mark.parametrize("scale", ["shared", "independent"])
def test_native_jcvi_uses_one_reversed_middle_track_and_preserves_all_inputs(saved_panel, tmp_path, monkeypatch, scale):
    from jcvi.graphics.karyotype import Track
    from pypdf import PdfReader

    sources = [Path(t["bed"]) for t in saved_panel["tracks"]] + [Path(e["blocks"]) for e in saved_panel["edges"]]
    original = {str(p): p.read_bytes() for p in sources}
    observed = []
    real = Track.get_coords
    def coordinates(track, gene):
        result = real(track, gene)
        observed.append((str(track), gene, result, set(track.rev)))
        return result
    monkeypatch.setattr(Track, "get_coords", coordinates)
    output = tmp_path / "plot"
    report = multi.render_panel(saved_panel, output, scale_mode=scale)
    assert {str(p): p.read_bytes() for p in sources} == original
    assert report["species_order"] == ["Anc_species", "Dio_species", "Tri_species"]
    assert report["display_orientation"]["Dio_species"] == {"chr1": "-"}
    assert report["displayed_block_counts"] == [1, 1]
    assert report["gaps"][1] * 518.4 == pytest.approx(8)
    assert report["scale_bar_count"] == (1 if scale == "shared" else 3)
    assert (output / "seqids").read_text() == "chr1\nchr1-\nchr1\n"
    middle = [(gene, point) for species, gene, point, reverse in observed if species == "Dio species" and reverse == {"chr1"}]
    assert {g for g, _ in middle} == {"d0", "d1", "d7"}
    for gene, point in middle:
        assert point[0] == pytest.approx(report["xstart"] + report["track_ratios"][1] * (7 - int(gene[1:])))
    pdf = PdfReader(output / "karyotype.pdf")
    assert len(pdf.pages) == 1
    assert float(pdf.pages[0].mediabox.width) == pytest.approx(518.4)
    for species in ("Anc species", "Dio species", "Tri species"):
        assert pdf.pages[0].extract_text().count(species) == 1
    assert "(rev)" in pdf.pages[0].extract_text()
    with Image.open(output / "karyotype.png") as image:
        image.verify()
    manifest = json.loads((output / "manifest.json").read_text())
    assert all((output / t["bed"]).is_file() for t in manifest["tracks"])
    assert all((output / e["blocks"]).is_file() for e in manifest["edges"])


def test_standalone_saved_bundle_redraws_from_another_working_directory(saved_panel, tmp_path):
    one = tmp_path / "one"
    before = multi.render_panel(saved_panel, one, formats=("svg",))
    destination = tmp_path / "two"
    result = subprocess.run([sys.executable, str(REPO / "workflow/support/synteny_multitrack.py"),
                             "--manifest", str(one / "manifest.json"), "--output", str(destination), "--formats", "svg"],
                            cwd=tmp_path.parent, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    after = json.loads((destination / "karyotype_layout.json").read_text())
    for key in ("display_order", "display_orientation", "objective_after", "species_order", "displayed_block_counts"):
        assert before[key] == after[key]
    assert (destination / "karyotype.svg").is_file()


def test_source_mutation_during_save_does_not_publish_a_bundle(saved_panel, tmp_path, monkeypatch):
    import matplotlib.pyplot as plt
    from matplotlib.figure import Figure

    original = Figure.savefig
    def changed(fig, path, **kwargs):
        original(fig, path, **kwargs)
        source = Path(saved_panel["tracks"][1]["bed"])
        replacement = source.with_suffix(".replacement")
        replacement.write_bytes(source.read_bytes())
        replacement.replace(source)
    monkeypatch.setattr(Figure, "savefig", changed)
    figures = plt.get_fignums()
    with pytest.raises(OSError, match="changed"):
        multi.render_panel(saved_panel, tmp_path / "new", formats=("pdf",))
    assert not (tmp_path / "new").exists()
    assert plt.get_fignums() == figures


def test_existing_bundle_is_not_overwritten(saved_panel, tmp_path):
    output = tmp_path / "old"
    output.mkdir()
    (output / "user.txt").write_text("keep")
    with pytest.raises(ValueError, match="must be new"):
        multi.render_panel(saved_panel, output)
    assert (output / "user.txt").read_text() == "keep"


def test_custom_labels_preserve_gene_and_sequence_identifiers(saved_panel, tmp_path):
    saved_panel["tracks"][1]["labels"] = {"chr1": "1"}
    output = tmp_path / "short-labels"
    report = multi.render_panel(saved_panel, output, formats=("svg",))
    assert report["display_order"]["Dio_species"] == ["chr1"]
    assert report["display_labels"]["Dio_species"] == {"chr1": "1"}
    assert "1 (rev)" in (output / "karyotype.svg").read_text()
    assert "chr1-" in (output / "seqids").read_text()
    assert (output / "inputs/track1.bed").read_bytes() == Path(saved_panel["tracks"][1]["bed"]).read_bytes()
    saved_panel["tracks"][1]["labels"] = {"unknown": "1"}
    with pytest.raises(ValueError, match="labels"):
        multi.render_panel(saved_panel, tmp_path / "invalid-labels", formats=("svg",))


def test_pair_renderer_accepts_generated_reverse_token(saved_panel, tmp_path):
    from workflow.support.pairwise_synteny_karyotype import render_karyotype

    directory = tmp_path / "pair"
    directory.mkdir()
    for side, track in zip(("target", "query"), saved_panel["tracks"], strict=False):
        (directory / (side + ".bed")).write_bytes(Path(track["bed"]).read_bytes())
    (directory / "seqids").write_text("chr1\nchr1-\n")
    (directory / "layout").write_text("0.7,0.12,0.92,0,,Anc,top,target.bed,top\n0.3,0.12,0.92,0,,Dio,bottom,query.bed,bottom\ne,0,1,blocks.simple\n")
    (directory / "blocks.simple").write_text("#88afc4*a0 a7 d0 d7 8 -\n")
    colors = {"chromosomes": {s: {"chr1": "#88afc4"} for s in ("target", "query")}}
    style = render_karyotype(directory, {"target_species": "Anc_species", "query_species": "Dio_species"}, "svg", colors)
    assert style["display_orientation"] == {"target": {"chr1": "+"}, "query": {"chr1": "-"}}
    assert "chr1 (rev)" in (directory / "karyotype.svg").read_text()


def test_three_species_workflow_panels_are_plot_only_and_reject_inconsistent_shared_beds(tmp_path):
    workspace = tmp_path / "workspace"
    for species, reverse in (("Anc_species", False), ("Dio_species", True), ("Tri_species", False)):
        write_genome(workspace, species, "protein", (("Chr1", reverse),))
    (workspace / "input/synteny_pairs.tsv").write_text("analysis_id\ttarget_species\tquery_species\nad\tAnc_species\tDio_species\ndt\tDio_species\tTri_species\n")
    first = run_core(workspace, synteny_plot_formats="svg", synteny_karyotype_sort="none")
    assert first.returncode == 0, first.stdout + first.stderr
    root = workspace / "output/genome_evolution/synteny"
    analysis_before = {str(p.relative_to(root)): (p.read_bytes(), p.stat().st_mtime_ns) for p in (root / "analysis").rglob("*") if p.is_file()}
    panels = workspace / "input/panels.json"
    config = {"panels": [{"panel_id": "anc_dio_tri", "tracks": [{"species": s} for s in ("Anc_species", "Dio_species", "Tri_species")],
                           "comparisons": ["ad", "dt"]}]}
    panels.write_text(json.dumps(config))
    result = run_core(workspace, synteny_plot_only="1", synteny_plot_formats="svg", synteny_karyotype_sort="order_orientation",
                      synteny_karyotype_panels_file=str(panels), artifact_stale_policy="rebuild")
    assert result.returncode == 0, result.stdout + result.stderr
    assert analysis_before == {str(p.relative_to(root)): (p.read_bytes(), p.stat().st_mtime_ns) for p in (root / "analysis").rglob("*") if p.is_file()}
    output = root / "plots/panels/anc_dio_tri"
    report = json.loads((output / "karyotype_layout.json").read_text())
    assert report["species_order"] == ["Anc_species", "Dio_species", "Tri_species"]
    assert report["orientation_changed"] is True
    assert (output / "karyotype.svg").is_file()
    pair_report = json.loads((root / "plots/ad/karyotype_order.json").read_text())
    assert isinstance(pair_report["objective_after"], float)
    assert pair_report["objective_after"] == pair_report["objective_components_after"]["total"]
    unchanged = (output / "karyotype.svg").stat().st_mtime_ns
    again = run_core(workspace, synteny_plot_only="1", synteny_plot_formats="svg", synteny_karyotype_sort="order_orientation",
                     synteny_karyotype_panels_file=str(panels))
    assert again.returncode == 0, again.stdout + again.stderr
    assert (output / "karyotype.svg").stat().st_mtime_ns == unchanged
    config["panels"][0]["tracks"][1]["fixed_orientation"] = {"Chr1": "+"}
    panels.write_text(json.dumps(config))
    changed = run_core(workspace, synteny_plot_only="1", synteny_plot_formats="svg", synteny_karyotype_sort="order_orientation",
                       synteny_karyotype_panels_file=str(panels), artifact_stale_policy="rebuild")
    assert changed.returncode == 0, changed.stdout + changed.stderr
    updated = json.loads((output / "karyotype_layout.json").read_text())
    assert updated["display_orientation"]["Dio_species"] == {"Chr1": "+"}
    assert (output / "karyotype.svg").stat().st_mtime_ns != unchanged
    assert analysis_before == {str(p.relative_to(root)): (p.read_bytes(), p.stat().st_mtime_ns) for p in (root / "analysis").rglob("*") if p.is_file()}
    # The generic assembler checks the repeated middle species before combining
    # native inputs; do not infer a coordinate conversion or silently mix them.
    (root / "analysis/dt/target.bed").write_text((root / "analysis/dt/target.bed").read_text().replace("\t0\t", "\t1\t", 1))
    plan = {"pairs": [{"analysis_id": "ad", "target": {"species": "Anc_species"}, "query": {"species": "Dio_species"}},
                      {"analysis_id": "dt", "target": {"species": "Dio_species"}, "query": {"species": "Tri_species"}}],
            "karyotype_panels": config["panels"], "formats": ["svg"]}
    with pytest.raises(ValueError, match="different gene coordinates"):
        multi.render_workflow_panels(plan, root / "analysis", tmp_path / "inconsistent")


def test_plot_contract_config_does_not_invalidate_analysis_or_ds(tmp_path):
    from workflow.support.pairwise_synteny import contract_args

    plan = {"workspace": str(tmp_path), "pairs": [], "parameters": {}, "tools": {"jcvi": "test"}, "formats": ["svg"],
            "ds_tools": {"source_hashes": {}}, "karyotype_panels_file": str(tmp_path / "panels.json"),
            "karyotype_panels": [{"panel_id": "example"}]}
    for phase in ("analysis", "ds", "plots"):
        before = contract_args(plan, phase)
        plan["karyotype_panels"][0]["seed"] = 42
        after = contract_args(plan, phase)
        assert (before == after) == (phase != "plots")
        plan["karyotype_panels"][0].pop("seed")
