"""Method-based figures retain counting units, provenance and physical geometry."""

import copy
import hashlib
import importlib
import json
import sys
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"
if str(SUPPORT) not in sys.path:
    sys.path.insert(0, str(SUPPORT))

figure = importlib.import_module("gene_model_improvement_summary")
methods = importlib.import_module("refinement_path_methods")


def fixture(groups=100):
    rows = []
    for index in range(2):
        rows.append(
            dict(
                species=f"Genus{index}_species",
                analysed=True,
                busco={
                    phase: dict(single=groups - 3, duplicated=1, fragmented=1, missing=1) for phase in figure.PHASES
                },
                rescue=dict(
                    method_counts=dict(synteny_based=3, miniprot=3, oriongeno=2),
                    unique_loci=5,
                    loci_supported_by_multiple_methods=2,
                    exclusive_intersections={
                        "synteny_based": 1,
                        "miniprot": 1,
                        "oriongeno": 1,
                        "miniprot+synteny_based": 1,
                        "miniprot+oriongeno+synteny_based": 1,
                    },
                ),
                coding_paths=dict(repairs=1, isoforms=2, counts=[1, 1, 1], colours=list(figure.ORIGIN_COLOURS)),
                utr=dict(
                    final_genes=10,
                    sides=dict(five_only=2, three_only=1, both=4, none=3),
                    origins=dict(original_gff_only=3, oriongeno_only=2, both=2, none=3),
                ),
            )
        )
    rows.append(
        dict(
            species="Excluded_species",
            analysed=False,
            busco=copy.deepcopy(rows[0]["busco"]),
            rescue=None,
            coding_paths=None,
            utr=None,
        )
    )
    return dict(
        schema_version=1,
        mode="proteins",
        groups=groups,
        species=rows,
        coding_path_legend=dict(labels=list(methods.LABELS), colours=list(methods.COLOURS)),
    )


def accepted_models():
    models = []
    identity = dict(GFF_sha256="a" * 64, receipt_sha256="b" * 64)
    for index, origins in enumerate((("miniprot",), ("oriongeno",), ("miniprot", "oriongeno"))):
        alignments = [
            dict(
                donor_species="donor",
                supports_path=True,
                **(
                    dict(
                        oriongeno={
                            **identity,
                            "native_transcript_id": "native_t",
                            "native_gene_id": "native_g",
                            "genomic_structure_transfer": {"confirmed": True},
                        }
                    )
                    if origin == "oriongeno"
                    else {}
                ),
            )
            for origin in origins
        ]
        models.append(
            dict(
                status="accepted",
                gene_id="one_locus",
                change_type="model_revision" if index == 0 else "isoform_addition",
                donors=["donor"],
                candidate=dict(candidate_id=f"path{index}", support=dict(donors=["donor"])),
                alignments=alignments,
            )
        )
    return models, identity


def publication(tmp_path):
    root = tmp_path / "refinement"
    root.mkdir()
    models, identity = accepted_models()
    plan = dict(
        request=dict(
            sources={f"Genus{i}_species": {} for i in range(2)},
            oriongeno=dict(predictions={f"Genus{i}_species": dict(identity=identity) for i in range(2)}),
        )
    )
    raw = json.dumps(plan).encode()
    (root / "plan.json").write_bytes(raw)
    plan_sha = hashlib.sha256(raw).hexdigest()
    for name in plan["request"]["sources"]:
        directory = root / "predictions" / name
        directory.mkdir(parents=True)
        raw = json.dumps(models).encode()
        (directory / "predictions.json").write_bytes(raw)
        (directory / "receipt.json").write_text(
            json.dumps(dict(key=dict(plan=plan_sha), files={"predictions.json": hashlib.sha256(raw).hexdigest()}))
        )
    return root


def bound_fixture(root):
    data = fixture()
    data["refinement_plan_sha256"] = hashlib.sha256((root / "plan.json").read_bytes()).hexdigest()
    return data


def test_methods_count_paths_not_loci_and_qc_is_not_a_second_nomination():
    models, native = accepted_models()
    result = methods.classify_path_methods(models, native)
    assert result["counts"] == dict(miniprot_only=1, oriongeno_only=1, both=1)
    assert len(result["paths"]) == 3
    assert {p["gene_id"] for p in result["paths"].values()} == {"one_locus"}
    assert result["paths"]["path1"]["category"] == "oriongeno_only"
    models[1]["alignments"].append(dict(donor_species="other", supports_path=False))
    assert methods.classify_path_methods(models, native) == result


@pytest.mark.parametrize("fault", ["duplicate", "imported", "annevo", "native_identity", "donors"])
def test_unresolved_nomination_evidence_is_not_invented(fault):
    models, native = accepted_models()
    if fault == "duplicate":
        models.append(models[0])
    if fault == "imported":
        models[0]["alignments"][0]["rescue_revision_model_id"] = "old_path"
    if fault == "annevo":
        models[0]["alignments"][0]["annevo"] = {"native": True}
    if fault == "native_identity":
        native["GFF_sha256"] = "c" * 64
    if fault == "donors":
        models[0]["donors"] = ["other"]
    with pytest.raises(ValueError):
        methods.classify_path_methods(models, native)


def test_receipt_bound_reclassification_does_not_modify_source_models(tmp_path):
    root = publication(tmp_path)
    data = bound_fixture(root)
    original = copy.deepcopy(data)
    before = {p: p.read_bytes() for p in root.rglob("*") if p.is_file()}
    derived, proof = methods.collect_path_methods(data, root)
    assert data == original
    assert derived["species"][0]["coding_paths"]["counts"] == [1, 1, 1]
    assert proof["predictors_executed"] is False and len(proof["species"]) == 2
    assert all(p.read_bytes() == raw for p, raw in before.items())


@pytest.mark.parametrize("fault", ["content", "plan", "species", "totals", "binding"])
def test_reclassification_rejects_wrong_publication_and_counts(tmp_path, fault):
    root = publication(tmp_path)
    data = bound_fixture(root)
    directory = root / "predictions" / "Genus0_species"
    if fault == "content":
        (directory / "predictions.json").write_text("[]")
    if fault == "plan":
        (root / "plan.json").write_text((root / "plan.json").read_text() + "\n")
    if fault == "species":
        data["species"][0]["species"] = "Other_species"
    if fault == "totals":
        data["species"][0]["coding_paths"]["repairs"] += 1
    if fault == "binding":
        data["refinement_plan_sha256"] = "0" * 64
    with pytest.raises(ValueError):
        methods.collect_path_methods(data, root)


@pytest.mark.parametrize("groups", [10, 2026, 5000])
def test_figure_preserves_units_dual_axes_palettes_and_bar_geometry(tmp_path, monkeypatch, groups):
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_hex
    from matplotlib.figure import Figure
    from matplotlib.text import Text

    data = fixture(groups)
    data["species"][0]["busco"]["before"]["single"] += 1
    data["species"][0]["busco"]["before"]["missing"] -= 1
    data["species"][0]["busco"]["after"] = copy.deepcopy(data["species"][0]["busco"]["before"])
    original = Figure.savefig
    inspected = []

    def save(fig, path, **kwargs):
        fig.canvas.draw()
        assert fig.get_size_inches()[0] == pytest.approx(7.2)
        assert len(fig.axes) == 5 and len({round(a.get_position().width, 9) for a in fig.axes}) == 1
        assert all(t.get_fontsize() == 8 for t in fig.findobj(match=Text))
        assert fig.axes[1].get_title(loc="left") == "Refined coding\npaths"
        assert fig.axes[2].get_title(loc="left") == "BUSCO\ngroups"
        assert fig.axes[3].get_xlim()[0] == 0
        heights = []
        for ax in fig.axes:
            for bar in ax.patches:
                if bar.get_height() < 0.5:
                    a, z = (
                        ax.transData.transform((0, bar.get_y())),
                        ax.transData.transform((0, bar.get_y() + bar.get_height())),
                    )
                    heights.append(abs(z[1] - a[1]) * 72 / fig.dpi)
        assert heights and min(heights) == pytest.approx(4.4) and max(heights) == pytest.approx(4.4)
        bars = fig.axes[1].patches[:5]
        assert [to_hex(p.get_facecolor()) for p in bars] == [*figure.PATH_CHANGE_COLOURS, *figure.ORIGIN_COLOURS]
        assert [p.get_width() for p in bars] == [1, 2, 1, 1, 1]
        for axis in fig.axes[2:4]:
            assert len(axis.child_axes) == 1
            for percent in (0, 0.5, 1):
                assert axis.transData.transform((percent, 0))[0] == pytest.approx(
                    axis.child_axes[0].transData.transform((percent * groups / 100, 0))[0]
                )
        assert sum("Not analysed" == t.get_text() for ax in fig.axes for t in ax.texts) == 4
        assert fig.axes[4].texts[0].get_text() == "n=10"
        renderer = fig.canvas.get_renderer()
        footer = fig.texts[0].get_window_extent(renderer)
        assert footer.y0 == pytest.approx(0.08 * fig.dpi)
        assert all(not legend.get_window_extent(renderer).overlaps(footer) for legend in fig.legends)
        inspected.append(Path(path).suffix)
        return original(fig, path, **kwargs)

    monkeypatch.setattr(Figure, "savefig", save)
    figure.plot_method_summary(data, tmp_path)
    assert inspected == [".pdf", ".png", ".svg"] and not plt.get_fignums()
    assert (tmp_path / "busco_comparison.pdf").read_bytes().startswith(b"%PDF")


def test_negative_deltas_and_baseline_only_species_are_supported(tmp_path, monkeypatch):
    from matplotlib.figure import Figure

    data = fixture(10)
    data["species"][0]["busco"]["after"]["single"] -= 1
    data["species"][0]["busco"]["after"]["missing"] += 1

    def negative(fig, *args, **kwargs):
        assert fig.axes[3].get_xlim()[0] < 0

    monkeypatch.setattr(Figure, "savefig", negative)
    figure.plot_method_summary(data, tmp_path / "negative")
    data["species"] = data["species"][-1:]
    monkeypatch.setattr(Figure, "savefig", lambda *args, **kwargs: None)
    figure.plot_method_summary(data, tmp_path / "baseline_only")


def test_all_zero_deltas_keep_compact_ticks_and_separate_footer(tmp_path, monkeypatch):
    from matplotlib.figure import Figure

    def inspect(fig, *args, **kwargs):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        delta, utr = fig.axes[3:]
        visible = [
            t for t in delta.get_xticklabels() if delta.get_xlim()[0] <= t.get_position()[0] <= delta.get_xlim()[1]
        ]
        assert [t.get_text() for t in visible] == ["0", "0.5", "1"]
        assert visible[-1].get_window_extent(renderer).x1 < utr.get_xticklabels()[0].get_window_extent(renderer).x0
        assert all(
            not legend.get_window_extent(renderer).overlaps(fig.texts[0].get_window_extent(renderer))
            for legend in fig.legends
        )

    monkeypatch.setattr(Figure, "savefig", inspect)
    figure.plot_method_summary(fixture(2026), tmp_path)


@pytest.mark.parametrize("fault", ["busco", "origins", "paths", "loci", "intersection", "unanalysed", "donor_groups"])
def test_inconsistent_denominators_and_groupings_do_not_reach_render(fault):
    data = fixture()
    row = data["species"][0]
    if fault == "busco":
        row["busco"]["after"]["single"] += 1
    if fault == "origins":
        row["utr"]["origins"]["both"] += 1
    if fault == "paths":
        row["coding_paths"]["counts"][0] += 1
    if fault == "loci":
        row["rescue"]["unique_loci"] += 1
    if fault == "intersection":
        row["rescue"]["exclusive_intersections"]["miniprot+miniprot"] = 0
    if fault == "unanalysed":
        row["analysed"] = False
    if fault == "donor_groups":
        data["coding_path_legend"]["labels"][0] = "Nearest relatives only"
    with pytest.raises(ValueError):
        figure.validate(data)


def test_report_only_cli_uses_verified_predictions_and_writes_provenance(tmp_path, monkeypatch):
    plot = importlib.import_module("gene_model_improvement_plot")
    root = publication(tmp_path)
    data = tmp_path / "counts.json"
    data.write_text(json.dumps(bound_fixture(root)))
    output = tmp_path / "report"
    monkeypatch.setattr(
        "sys.argv", ["plot", "--figure-data", str(data), "--output", str(output), "--refinement-output", str(root)]
    )
    plot.main()
    proof = json.loads((output / "improvement_method_provenance.json").read_text())
    assert not proof["predictors_executed"] and not proof["busco_executed"]
    assert proof["figure_data_sha256"] == hashlib.sha256(data.read_bytes()).hexdigest()
    assert len(proof["path_methods"]["species"]) == 2


def test_cli_keeps_reports_outside_immutable_refinement(tmp_path, monkeypatch):
    plot = importlib.import_module("gene_model_improvement_plot")
    root = publication(tmp_path)
    data = tmp_path / "counts.json"
    data.write_text(json.dumps(bound_fixture(root)))
    monkeypatch.setattr(
        "sys.argv",
        ["plot", "--figure-data", str(data), "--output", str(root / "figures"), "--refinement-output", str(root)],
    )
    with pytest.raises(SystemExit) as error:
        plot.main()
    assert error.value.code == 2 and not (root / "figures").exists()


def test_cli_does_not_overwrite_its_source_or_accept_a_racing_input(tmp_path, monkeypatch):
    plot = importlib.import_module("gene_model_improvement_plot")
    source = tmp_path / "busco_comparison.svg"
    source.write_text(json.dumps(fixture()))
    monkeypatch.setattr("sys.argv", ["plot", "--figure-data", str(source), "--output", str(tmp_path)])
    before = source.read_bytes()
    with pytest.raises(SystemExit):
        plot.main()
    assert source.read_bytes() == before
    monkeypatch.setattr("sys.argv", ["plot", "--figure-data", str(source), "--output", str(tmp_path / "new")])
    monkeypatch.setattr(plot, "plot_method_summary", lambda *args: source.write_text("{}"))
    with pytest.raises(OSError, match="changed"):
        plot.main()
    assert not (tmp_path / "new" / "improvement_method_provenance.json").exists()
