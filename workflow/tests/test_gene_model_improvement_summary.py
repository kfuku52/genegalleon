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
def test_figure_preserves_absolute_counts_and_matching_label_colors(tmp_path, monkeypatch, groups):
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_rgba
    from matplotlib.figure import Figure

    data = fixture(groups)
    row = data["species"][0]
    row["coding_paths"]["architecture"] = {"merges": 4, "splits": 3}
    species = row["species"]
    omark = {species: {"input": ((10, 12), (9, 10))}}
    original = Figure.savefig
    inspected = []

    def save(fig, path, **kwargs):
        assert len(fig.axes) == 6
        _, refinement, utr, busco, hogs, proteins = fig.axes
        assert busco.get_title(loc="left") == "BUSCO\ngroups"
        bars = refinement.patches[:4]
        assert [bar.get_width() for bar in bars] == [1, 2, 4, 3]
        assert [bar.get_x() for bar in bars] == [0, 1, 3, 7]
        assert len({bar.get_y() for bar in bars}) == 1
        assert [t.get_text() for t in refinement.texts[:4]] == ["1", "2", "4", "3"]
        for bar, text in zip(bars, refinement.texts, strict=False):
            assert to_rgba(text.get_color()) == bar.get_facecolor()
            assert text.get_position()[1] > bar.get_y() + bar.get_height()
        assert [bar.get_width() for bar in busco.patches[:4]] == [groups - 3, 1, 1, 1]
        assert [bar.get_width() for bar in hogs.patches] == [10, 2]
        assert [bar.get_width() for bar in proteins.patches] == [9, 1]
        assert proteins.texts[0].get_text() == "90.0% (n=12)"
        assert "70.0%" in [t.get_text() for t in utr.texts]
        assert sum(t.get_text() == "Not analyzed" for ax in fig.axes for t in ax.texts) == 3
        assert not fig.texts
        inspected.append(Path(path).suffix)
        return original(fig, path, **kwargs)

    monkeypatch.setattr(Figure, "savefig", save)
    figure.plot_method_summary(data, tmp_path, omark=omark, gene_totals={(species, "input"): 12})
    assert inspected == [".pdf", ".png", ".svg"] and not plt.get_fignums()
    assert (tmp_path / "busco_comparison.pdf").read_bytes().startswith(b"%PDF")


def test_baseline_only_species_are_supported(tmp_path, monkeypatch):
    from matplotlib.figure import Figure

    data = fixture(10)
    data["species"] = data["species"][-1:]

    def inspect(fig, *args, **kwargs):
        assert len(fig.axes) == 4
        assert sum(t.get_text() == "Not analyzed" for ax in fig.axes for t in ax.texts) == 3
        assert len(fig.axes[-1].patches) == 4

    monkeypatch.setattr(Figure, "savefig", inspect)
    figure.plot_method_summary(data, tmp_path)


@pytest.mark.parametrize('column_scale', [1.0, .8])
def test_measured_legend_space_preserves_bar_geometry_in_a_compact_layout(tmp_path, monkeypatch, column_scale):
    from matplotlib.figure import Figure

    original = figure.layout
    def compact(rows, **kwargs):
        geometry = original(rows, **kwargs)
        return {**geometry, 'bottom': geometry['bottom'] - 1.0, 'height': geometry['height'] - 1.0,
                'positions': [(x, width * column_scale) for x, width in geometry['positions']]}
    monkeypatch.setattr(figure, 'layout', compact)
    inspected = []
    def inspect(fig, path, **kwargs):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        footer = next(t for t in fig.texts if t.get_text().startswith('Miniprot ')).get_window_extent(renderer)
        assert fig.get_size_inches()[0] == pytest.approx(5.99)
        assert fig.get_size_inches()[1] > compact(3)['height']
        assert footer.y0 == pytest.approx(.15 * fig.dpi)
        assert len(fig.legends) == 7 and len(fig.texts) == 2
        titles = [''.join(legend.get_title().get_text().split()) for legend in fig.legends]
        assert 'Swiss-ProtannotationoftheMiniprot/OrionGenounion' in titles
        assert all(not legend.get_window_extent(renderer).overlaps(footer) for legend in fig.legends)
        for ax in fig.axes:
            for bar in ax.patches:
                if bar.get_height() < .5:
                    y0 = ax.transData.transform((0, bar.get_y()))[1]
                    y1 = ax.transData.transform((0, bar.get_y() + bar.get_height()))[1]
                    assert abs(y1 - y0) * 72 / fig.dpi == pytest.approx(4.4)
        inspected.append(Path(path).suffix)
    monkeypatch.setattr(Figure, 'savefig', inspect)
    figure.plot_method_summary(fixture(2026), tmp_path, show_method_notes=True)
    assert inspected == ['.pdf', '.png', '.svg']


def test_swissprot_union_baseline_slot_and_column_legends(tmp_path, monkeypatch):
    from matplotlib.colors import to_hex
    from matplotlib.figure import Figure

    data = fixture(2026)
    data["busco_dataset"] = "embryophyta_odb12"
    data["species"][0]["rescue"]["swissprot_support"] = dict(zip(figure.SWISS_KEYS, (1, 2, 1, 1, 0), strict=True))

    def inspect(fig, *args, **kwargs):
        fig.canvas.draw()
        assert "embryophyta_odb12" in fig.axes[3].get_title(loc="left")
        bars = fig.axes[0].patches[:9]
        assert [p.get_width() for p in bars[-5:]] == [1, 2, 1, 1, 0]
        assert [to_hex(p.get_facecolor()) for p in bars] == [
            *figure.MINIPROT_STACK_COLOURS, figure.METHOD_COLOURS[2], *figure.SWISS_COLOURS,
        ]
        # Historical inputs lacking support counts must not become no-hit calls.
        assert [p.get_width() for p in fig.axes[0].patches[13:18]] == [0, 0, 0, 0, 5]
        baseline = fig.axes[3].patches[-4:]
        assert all(p.get_y() + p.get_height() / 2 == pytest.approx(2 + figure.bar_offsets(3)[0]) for p in baseline)
        titles = [legend.get_title().get_text().replace("\n", " ") for legend in fig.legends]
        assert titles.count("BUSCO categories") == 1
        assert {"UTR ends", "UTR sources", "Model changes", "Miniprot", "OrionGeno"} <= set(titles)
        # Width-based wrapping may split a long word. Keep every character of
        # the complete title while allowing only layout whitespace to vary.
        assert 'Swiss-ProtannotationoftheMiniprot/OrionGenounion' in [''.join(title.split()) for title in titles]
        assert not fig.texts
        for column, counts in enumerate((3, 1, 2, 1)):
            start = sum((3, 1, 2, 1)[:column])
            x, width = figure.layout(3)["positions"][column]
            legends = fig.legends[start:start + counts]
            for legend in legends:
                bbox = legend.get_window_extent(fig.canvas.get_renderer())
                assert bbox.x0 >= x * fig.dpi - .1 and bbox.x1 <= (x + width) * fig.dpi + .1
            for upper, lower in zip(legends, legends[1:], strict=False):
                assert upper.get_window_extent().y0 - lower.get_window_extent().y1 >= .09 * fig.dpi - .1

    monkeypatch.setattr(Figure, "savefig", inspect)
    figure.plot_method_summary(data, tmp_path)


def test_reported_thresholds_are_explicit_and_never_invented():
    data = fixture()
    assert "not recorded" in " ".join(figure.method_notes(data))
    data["rescue_parameters"] = dict(minimum_coverage=.7, minimum_identity=.5, minimum_broad_coverage=.9,
                                     minimum_lower_coverage_identity=.7, min_anchors=4, cscore=.7, max_interval=200000)
    data["rescue_swissprot_evidence"] = dict(parameters=dict(evalue=1e-5, minimum_alignment=50,
                                                           query_coverage=.5, target_coverage=.6, score_fraction=.9))
    text = " ".join(figure.method_notes(data))
    for expected in ("coverage >= 70%", "identity >= 70% at 70-<90%", ">= 50% at >= 90%",
                     "blocks >= 4", "C-score >= 0.7", "200 kb", "E <= 1e-05", "50%/60%", "90% of best"):
        assert expected in text
    for key, value in (("minimum_coverage", float("nan")), ("min_anchors", True), ("cscore", 2)):
        changed = copy.deepcopy(data)
        changed["rescue_parameters"][key] = value
        with pytest.raises(ValueError):
            figure.validate(changed)


@pytest.mark.parametrize("fault", ["total", "negative", "missing", "bool", "dataset", "swiss_threshold"])
def test_swissprot_reporting_rejects_invalid_counts_or_metadata(fault):
    data = fixture()
    support = dict(zip(figure.SWISS_KEYS, (1, 2, 1, 1, 0), strict=True))
    data["species"][0]["rescue"]["swissprot_support"] = support
    if fault == "total":
        support["both"] += 1
    if fault == "negative":
        support["not_assessed"] = -1
    if fault == "missing":
        del support["both"]
    if fault == "bool":
        support["te_only"] = True
    if fault == "dataset":
        data["busco_dataset"] = "embryophyta\nodb12"
    if fault == "swiss_threshold":
        data["rescue_swissprot_evidence"] = dict(parameters=dict(evalue=0))
    with pytest.raises(ValueError):
        figure.validate(data)


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


def test_staged_output_uses_same_renderer_with_bound_omark(tmp_path, monkeypatch):
    plot = importlib.import_module("gene_model_improvement_plot")
    prepared = fixture(2026)
    changes = {"plan_sha256": "plan", "effective_receipt_sha256": "receipt", "species": {}}
    rows = []
    for entry in prepared["species"]:
        species = entry["species"]
        rows.append({"species": species, "refinement_status": "analysed" if entry["analysed"] else "not_analysed",
                     **{f"{phase}_result": {**result, "total": 2026, "mode": "proteins", "lineage": "test"}
                        for phase, result in entry["busco"].items()}})
        if entry["analysed"]:
            changes["species"][species] = {
                "prior_rescued_loci": entry["rescue"]["unique_loci"],
                "rescue_method_intersections": entry["rescue"]["exclusive_intersections"],
                "figure_annotation": {"utr": entry["utr"], "repairs": (1, 2),
                                      "architecture": {"merges": 3, "splits": 4}},
            }
    output = tmp_path / "busco"
    name = rows[0]["species"]
    omark = tmp_path / "omark_comparison.tsv"
    omark.write_text(
        "species\tstage\trepresentative_proteins\tconsistent_full\tconserved_hogs\tmissing_hogs\t"
        "annotated_gene_loci\trefinement_plan_sha256\teffective_receipt_sha256\n"
        f"{name}\trefined\t10\t9\t12\t2\t13\tplan\treceipt\n"
    )
    calls = []
    monkeypatch.setattr(plot, "plot_method_summary", lambda data, dest, **opts: calls.append((data, dest, opts)))
    plot.plot_staged_comparison(rows, output, changes)
    data, destination, opts = calls.pop()
    assert destination == output
    assert data["species"][0]["coding_paths"]["architecture"] == {"merges": 3, "splits": 4}
    assert opts["omark"] == {name: {"refined": ((10, 12), (9, 10))}}
    assert opts["gene_totals"] == {(name, "refined"): 13}
    omark.write_text(omark.read_text().replace("\tplan\treceipt", "\tplan\tstale"))
    plot.plot_staged_comparison(rows, output, changes)
    assert "omark" not in calls.pop()[2]


def test_report_cli_accepts_relative_quality_paths(tmp_path, monkeypatch):
    plot = importlib.import_module("gene_model_improvement_plot")
    species = fixture()["species"][0]["species"]
    (tmp_path / "figure.json").write_text(json.dumps(fixture()))
    (tmp_path / "quality.tsv").write_text(
        "species\tstage\trepresentative_proteins\tconsistent_full\tconserved_hogs\tmissing_hogs\n"
        f"{species}\trefined\t10\t9\t12\t2\n"
    )
    (tmp_path / "totals.tsv").write_text(f"species\tstage\ttotal_genes\n{species}\trefined\t13\n")
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr("sys.argv", ["plot", "--figure-data", "figure.json", "--output", "report",
                                     "--omark-comparison", "quality.tsv", "--gene-totals", "totals.tsv"])
    calls = []
    monkeypatch.setattr(plot, "plot_method_summary", lambda data, dest, **opts: calls.append(opts))
    plot.main()
    assert calls[0]["gene_totals"] == {(species, "refined"): 13}


def test_report_uses_zero_for_unpublished_structure_events(tmp_path):
    import gzip

    summaries = importlib.import_module("refinement_figure_inputs")
    name = "Example_species"
    audit = f"utr_audit/{name}.json.gz"
    inputs = {
        "utr_summary.tsv": ("species\texported_transcripts\twith_5prime_utr\twith_3prime_utr\twith_both_utrs\n"
                            f"{name}\t1\t1\t1\t1\n").encode(),
        audit: gzip.compress(json.dumps([{"utr5": {"length": 15, "source": "selected"},
                                          "utr3": {"length": 22, "source": "oriongeno"}}]).encode()),
        "split_model_merges.json": json.dumps([{"species": name}]).encode(),
    }
    (tmp_path / "utr_audit").mkdir()
    for relative, content in inputs.items():
        (tmp_path / relative).write_bytes(content)
    publication = {"files": {path: hashlib.sha256(content).hexdigest() for path, content in inputs.items()}}
    annotation = summaries.refinement_annotations(tmp_path, publication, {name})[name]
    assert annotation["architecture"] == {"merges": 1, "splits": 0}
    assert annotation["repairs"] == (0, 0)
    assert annotation["utr"]["origins"]["both"] == 1


def test_gene_total_counts_coding_loci_once(tmp_path):
    artifact = importlib.import_module("refinement_figure_inputs")
    gff = tmp_path / "source.gff3"
    gff.write_text(
        "chr\tref\tgene\t1\t50\t.\t+\t.\tID=g1\n"
        "chr\tref\tmRNA\t1\t50\t.\t+\t.\tID=t1;Parent=g1\n"
        "chr\tref\tmRNA\t1\t50\t.\t+\t.\tID=t2;Parent=g1\n"
        "chr\tref\tCDS\t4\t30\t.\t+\t0\tParent=t1,t2\n"
        "chr\tref\tgene\t60\t90\t.\t+\t.\tID=g2;gene_biotype=lncRNA\n"
        "chr\tref\tgene\t91\t99\t.\t+\t.\tID=g3;gene_biotype=protein_coding\n"
    )
    assert artifact._coding_loci(gff) == 2
