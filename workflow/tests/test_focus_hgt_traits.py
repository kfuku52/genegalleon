import json
import sys
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT))

from focus_hgt_gene_trees import annotate, background_supported  # noqa: E402
from focus_hgt_traits import generate, read_tsv, write_tsv  # noqa: E402
from species_trait_schema import schema_path, schema_payload  # noqa: E402


def supported_link(event_id, side, gene):
    return dict(event_id=event_id, orthogroup="OG1", side=side, gene_id=gene,
                eligible_for_context="True", host_scaffold_status="measured", host_scaffold_id="scaffold1",
                host_scaffold_background_class_total_count="20", host_scaffold_background_class_compatible_count="9",
                host_scaffold_background_class_incompatible_count="1", host_scaffold_background_class_unresolved_count="10",
                host_scaffold_background_class_classified_fraction="0.5",
                host_scaffold_background_class_compatible_fraction="0.9")


def focused_node_source():
    event = dict(event_id="OG1:3:1", orthogroup="OG1", gene_tree_branch_id="3", gene_tree_node="n3",
                 event_index="1", generax_transfer="Y@D@A", support_used="90")
    stat = [dict(branch_id="0", node_name="A_gene", child1="-999", child2="-999", generax_transfer="N",
                 support_generax_ufboot=""),
            dict(branch_id="3", node_name="n3", child1="0", child2="1", generax_transfer="Y@D@A;Y@D@B",
                 support_generax_ufboot="90"),
            dict(branch_id="1", node_name="D_gene", child1="-999", child2="-999", generax_transfer="N",
                 support_generax_ufboot="")]
    links = [supported_link(event["event_id"], side, gene) for side, gene in
             (("donor", "D_gene"), ("recipient", "A_gene"))]
    return stat, [event], links


def test_focused_gene_nodes_use_exact_event_branch_and_inclusive_bilateral_profile():
    stat, events, links = focused_node_source()
    output, audit = annotate(stat, events, links)
    assert audit[0]["status"] == "selected" and audit[0]["support_generax_ufboot"] == 90
    assert output[0]["hgtfocus_recipient_flag"] == 1
    assert output[1]["hgtfocus_event_ids"] == "OG1:3:1"
    assert output[1]["hgtfocus_node_label"] == "HGT1 UF=90"
    assert stat[1].get("hgtfocus_event_ids") is None  # Original table is preserved.
    # Two events sharing a gene-tree node are evaluated separately.
    other = dict(events[0], event_id="OG1:3:2", event_index="2", generax_transfer="Y@D@B")
    extra = supported_link(other["event_id"], "recipient", "B_gene")
    output, audit = annotate(stat, events + [other], links + [extra])
    assert [row["status"] for row in audit] == ["selected", "withheld"]
    assert output[1]["hgtfocus_event_count"] == 1  # Missing donor evidence cannot borrow from another event.


@pytest.mark.parametrize("alteration,reason", [
    ("node", "gene_tree_branch_node_unmapped"),
    ("terminal", "terminal_gene_tree_branch"),
    ("token", "transfer_token_or_event_index_unmapped"),
    ("index", "transfer_token_or_event_index_unmapped"),
    ("missing_support", "generax_ufboot_unavailable"),
    ("low_support", "generax_ufboot_below_threshold"),
    ("recipient_background", "bilateral_class_background_not_supported_or_unavailable"),
    ("context_gene", "supported_context_gene_unmapped_to_family_tip"),
])
def test_focused_gene_nodes_withhold_missing_or_unmapped_evidence(alteration, reason):
    stat, events, links = focused_node_source()
    if alteration == "node":
        events[0]["gene_tree_node"] = "wrong"
    elif alteration == "terminal":
        stat[1]["child1"] = "-999"
    elif alteration == "token":
        events[0]["generax_transfer"] = "Y@X@A"
    elif alteration == "index":
        events[0]["event_index"] = "2"
    elif alteration.endswith("support"):
        stat[1]["support_generax_ufboot"] = "" if alteration == "missing_support" else "89.9"
    elif alteration == "context_gene":
        links[0]['gene_id'] = 'other_family_gene'
    else:
        links[1]["host_scaffold_background_class_total_count"] = ""
    output, audit = annotate(stat, events, links)
    assert audit[0]["reason"] == reason
    assert all(row["hgtfocus_event_count"] == 0 for row in output)


def test_focused_gene_nodes_reject_conflicting_support_and_background_fractions():
    stat, events, links = focused_node_source()
    events[0]["support_used"] = "100"
    with pytest.raises(ValueError, match="support disagrees"):
        annotate(stat, events, links)
    links[0]["host_scaffold_background_class_compatible_fraction"] = "1"
    with pytest.raises(ValueError, match="background fraction"):
        background_supported(links[0])


def test_focused_gene_tree_folder_audits_missing_families_and_preserves_inputs(tmp_path):
    from focus_hgt_gene_trees import export_gene_trees

    stat, events, links = focused_node_source()
    root = tmp_path / "families"
    table = root / "stat_branch/OG1_stat.branch.tsv"
    write_tsv(table, list(stat[0]), stat)
    original = table.read_bytes()
    def renderer(annotated, pdf):
        assert read_tsv(annotated)[1][1]["hgtfocus_event_count"] == "1"
        pdf.write_bytes(b"%PDF-native-renderer-fixture")
    output = tmp_path / "tree_plot"
    report = export_gene_trees(output, events, links, root, renderer=renderer)
    assert report["rendered_family_count"] == report["selected_event_count"] == 1
    assert table.read_bytes() == original
    assert (output / "OG1_focused_hgt_tree_plot.pdf").is_file()
    missing = dict(events[0], orthogroup="OG2")
    report = export_gene_trees(tmp_path / "missing", [missing], [], root, renderer=renderer)
    assert report["selected_event_count"] == 0
    assert read_tsv(tmp_path / "missing/event_node_audit.tsv")[1][0]["reason"] == "stat_branch_unavailable"


@pytest.fixture
def source(tmp_path):
    tree = tmp_path / "tree.nwk"
    tree.write_text("((A:1,B:1)0042:1,(C:1,D:1)mixed:1)root;")
    trait = tmp_path / "traits.tsv"
    trait.write_text("species\tbinary\tcategory\tcontinuous\nA\t1\t1\t1\nB\t1\t2\t0\nC\t0\t1\t0.5\nD\tNA\tNA\t2\n")
    schema_path(trait).write_bytes(schema_payload(trait.read_bytes(),
                                                {"binary": "binary", "category": "categorical", "continuous": "numeric"}))
    events = []
    for number, recipient in enumerate(("A", "0042", "mixed", "B", "C"), 1):
        events.append(dict(event_id=f"OG1:3:{number}", orthogroup="OG1", generax_transfer=f"Y@D@{recipient}",
                           generax_donor_node="D", generax_recipient_node=recipient, mapping_status="matched",
                           support_generax_ufboot="90", existing_annotation="Product name", sequence_quality_status="NA"))
    event = tmp_path / "events.tsv"
    write_tsv(event, list(events[0]), events)
    links = []
    for row in events:
        for species, side in (("D", "donor"), ("A", "recipient"), ("B", "recipient")):
            links.append(dict(event_id=row["event_id"], orthogroup="OG1", gene_id=species + "_gene", gene_species=species,
                              side=side, eligible_for_context="True", product_name="Protein " + species, synteny_support_score="NA"))
    link = tmp_path / "links.tsv"
    write_tsv(link, list(links[0]), links)
    return event, link, tree, trait, tmp_path / "output"


def test_category_and_binary_focus_preserve_cohort_fields_and_event_identity(source):
    manifest = generate(*source, plots=False)
    root = source[-1]
    _, binary = read_tsv(root / "traits/binary/all_category1/events.tsv")
    assert {row["event_id"] for row in binary} == {"OG1:3:1", "OG1:3:2", "OG1:3:4"}
    _, category = read_tsv(root / "traits/category/all_category1/events.tsv")
    assert {row["event_id"] for row in category} == {"OG1:3:1", "OG1:3:5"}
    assert all(row["support_generax_ufboot"] == "90" and row["existing_annotation"] == "Product name" for row in binary)
    assert next(row for row in binary if row["generax_recipient_node"] == "0042")["recipient_clade_tip_labels"] == "A; B"
    assert all(row["focus_ancestral_state_inferred"] == "0" for row in binary)
    _, direct_a = read_tsv(root / "traits/binary/tips/A/direct_events.tsv")
    _, ancestors_a = read_tsv(root / "traits/binary/tips/A/ancestral_recipient_events.tsv")
    _, ancestors_b = read_tsv(root / "traits/binary/tips/B/ancestral_recipient_events.tsv")
    assert [row["event_id"] for row in direct_a] == ["OG1:3:1"]
    assert [row["event_id"] for row in ancestors_a] == [row["event_id"] for row in ancestors_b] == ["OG1:3:2"]
    assert not (root / "traits/binary/tips/D").exists()
    assert not (root / "traits/continuous").exists()
    assert manifest["source_event_count"] == 5
    assert next(row for row in manifest["trait_selection"] if row["trait"] == "continuous")["status"] == "excluded"
    _, recipient = read_tsv(root / "traits/binary/tips/A/recipient_genes.tsv")
    assert len(recipient) == 1 and recipient[0]["gene_species"] == "A"
    assert recipient[0]["product_name"] == "Protein A" and recipient[0]["synteny_support_score"] == "NA"


def test_category1_gene_tree_export_receives_only_aggregate_events(source, monkeypatch):
    import focus_hgt_gene_trees
    import focus_hgt_traits

    calls = []
    real_export = focus_hgt_traits.export_bundle
    def export_without_species_pdf(*args, **kwargs):
        # Keep this wiring test independent of the species-tree renderer.
        positional = list(args)
        positional[8] = False
        return real_export(*positional, **kwargs)
    def capture(directory, events, links, root):
        calls.append((directory, {e['event_id'] for e in events}, root))
        return dict(rendered_family_count=1, selected_event_count=len(events))
    monkeypatch.setattr(focus_hgt_traits, 'export_bundle', export_without_species_pdf)
    monkeypatch.setattr(focus_hgt_gene_trees, 'export_gene_trees', capture)
    report = generate(*source, plots=True, gene_family_root='existing-families')
    assert len(calls) == 2  # One aggregate per binary/categorical trait, no per-tip rendering.
    assert all(path.name == 'tree_plot' and root == 'existing-families' for path, _, root in calls)
    assert calls[0][1] == {'OG1:3:1', 'OG1:3:2', 'OG1:3:4'}
    assert calls[1][1] == {'OG1:3:1', 'OG1:3:5'}
    assert report['gene_tree_plots']['binary']['selected_event_count'] == 3


def test_empty_category1_targets_are_reported(source):
    source[0].write_text(source[0].read_text().replace("Y@D@A", "Y@D@D").replace("\tD\tA\t", "\tD\tD\t"))
    generate(*source, plots=False)
    summary = json.loads((source[-1] / "traits/category/tips/A/summary.json").read_text())
    assert summary["event_count"] == 0
    assert read_tsv(source[-1] / "traits/category/tips/A/events.tsv")[1] == []


@pytest.mark.parametrize("alteration,match", [
    ("duplicate_event", "Duplicate or empty transfer"),
    ("duplicate_link", "Duplicate or invalid"),
    ("family_mismatch", "identity mismatch"),
    ("tree_mismatch", "clade tips disagree"),
    ("tip_count_mismatch", "clade tip count disagrees"),
    ("token_mismatch", "transfer token disagrees"),
    ("schema_mismatch", "schema does not match"),
    ("bad_binary", "Invalid binary"),
    ("alias_duplicate", "unique nonempty species"),
    ("positive_outside_tree", "absent from analysis"),
])
def test_invalid_inputs_do_not_publish_partial_bundle(source, alteration, match):
    event, link, _, trait, out = source
    if alteration == "duplicate_event":
        with event.open("a") as handle:
            handle.write(event.read_text().splitlines()[1] + "\n")
    elif alteration == "duplicate_link":
        with link.open("a") as handle:
            handle.write(link.read_text().splitlines()[1] + "\n")
    elif alteration == "family_mismatch":
        link.write_text(link.read_text().replace("\tOG1\t", "\tOTHER\t", 1))
    elif alteration in {"tree_mismatch", "tip_count_mismatch"}:
        fields, rows = read_tsv(event)
        for row in rows:
            row["recipient_clade_tip_labels"] = "Z" if alteration == "tree_mismatch" else {"A": "A", "B": "B", "C": "C", "0042": "A; B", "mixed": "C; D"}[row["generax_recipient_node"]]
            row["recipient_clade_tip_count"] = "99"
        write_tsv(event, fields + ["recipient_clade_tip_labels", "recipient_clade_tip_count"], rows)
    elif alteration == "token_mismatch":
        event.write_text(event.read_text().replace("Y@D@A", "Y@C@A"))
    else:
        replacement = {"schema_mismatch": ("A\t1", "A\t0"), "bad_binary": ("A\t1", "A\t2"),
                       "alias_duplicate": ("B\t1", "A\t1"), "positive_outside_tree": ("A\t1", "Unknown\t1")}[alteration]
        trait.write_text(trait.read_text().replace(*replacement))
        if alteration != "schema_mismatch":
            schema = json.loads(schema_path(trait).read_text())["traits"]
            schema_path(trait).write_bytes(schema_payload(trait.read_bytes(), schema))
    with pytest.raises(ValueError, match=match):
        generate(*source, plots=False)
    assert not out.exists()
    assert not list(out.parent.glob(".hgt-trait-focus-*"))


def test_schema_free_binary_only_and_missing_retained(source):
    schema_path(source[3]).unlink()
    generate(*source, plots=False)
    root = source[-1]
    assert (root / "traits/binary").exists()
    assert not (root / "traits/category").exists()
    assert not (root / "traits/continuous").exists()
    _, values = read_tsv(root / "traits/binary/species_trait_category1.tsv")
    assert next(row for row in values if row["species"] == "D")["binary"] == ""


def test_native_plot_exports_all_selected_arrows(source, monkeypatch):
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.patches import FancyArrowPatch

    alphas = []
    original = PdfPages.savefig

    def capture(self, fig, **kwargs):
        alphas.extend(p.get_alpha() for ax in fig.axes for p in ax.patches if isinstance(p, FancyArrowPatch))
        return original(self, fig, **kwargs)

    monkeypatch.setattr(PdfPages, "savefig", capture)
    manifest = generate(*source, plots=True, arrow_alpha=0.4)
    assert manifest["transfer_arrow_alpha"] == 0.4
    assert manifest["plot_scope"] == "trait_aggregate_only"
    assert alphas and set(alphas) == {0.4}
    root = source[-1] / "traits/binary/all_category1"
    _, edges = read_tsv(root / "transfer_edges.tsv")
    assert sum(int(row["hgt_event_count"]) for row in edges) == 3
    assert (root / "transfer_tree.pdf").read_bytes().startswith(b"%PDF")
    assert {str(path.relative_to(source[-1])) for path in source[-1].rglob("*.pdf")} == {
        "traits/binary/all_category1/transfer_tree.pdf", "traits/category/all_category1/transfer_tree.pdf"}
    for target in manifest["result_index"]:
        if target["target_type"] != "aggregate":
            directory = source[-1] / target["relative_path"]
            _, target_edges = read_tsv(directory / "transfer_edges.tsv")
            assert sum(int(row["hgt_event_count"]) for row in target_edges) == target["event_count"]
            assert (directory / "events.tsv").is_file()
    from plot_hgt_summary import read_transfer_traits
    assert list(read_transfer_traits(source[3]).columns) == ["binary", "continuous"]
    # Republishing a managed legacy bundle removes old per-recipient PDFs.
    legacy = source[-1] / "traits/binary/tips/A/transfer_tree.pdf"
    legacy.write_bytes(b"previous recipient figure")
    generate(*source, plots=True, arrow_alpha=0.4)
    assert not legacy.exists()
    assert (legacy.parent / "events.tsv").is_file()
    assert (legacy.parent / "transfer_edges.tsv").is_file()


def test_republication_replaces_managed_bundle_and_failure_preserves_it(source, monkeypatch):
    import focus_hgt_traits
    from focus_hgt_traits import build_focus
    generate(*source, plots=False)
    first = (source[-1] / "manifest.json").read_bytes()
    generate(*source, plots=False)
    assert (source[-1] / "manifest.json").read_bytes() == first
    def fail(*args, **kwargs):
        build_focus(*args, **kwargs)
        raise ValueError("Intentional failure")
    monkeypatch.setattr(focus_hgt_traits, "build_focus", fail)
    with pytest.raises(ValueError, match="Intentional failure"):
        generate(*source, plots=False)
    assert (source[-1] / "manifest.json").read_bytes() == first


def test_unmanaged_output_and_input_containment_are_rejected(source):
    source[-1].mkdir()
    with pytest.raises(ValueError, match="unmanaged"):
        generate(*source, plots=False)
    with pytest.raises(ValueError, match="contain an input"):
        generate(*source[:-1], source[0].parent, plots=False)


def test_observation_columns_are_not_focus_traits(source):
    import hashlib

    from species_trait_contract import CONTRACT_VERSION
    trait = source[3]
    metadata = {"schema_version": CONTRACT_VERSION, "table_sha256": hashlib.sha256(trait.read_bytes()).hexdigest(),
                "traits": {"binary": {"role": "quality", "source": "user"}}}
    Path(str(trait) + ".metadata.json").write_text(json.dumps(metadata))
    generate(*source, plots=False)
    assert not (source[-1] / "traits/binary").exists()


def test_categorical_only_table_can_render_overview_and_focused_native_tree(source):
    from plot_hgt_summary import read_transfer_traits
    trait = source[3]
    trait.write_text("species\tcategory\nA\t1\nB\t2\nC\t1\nD\tNA\n")
    schema_path(trait).write_bytes(schema_payload(trait.read_bytes(), {"category": "categorical"}))
    assert read_transfer_traits(trait).shape == (4, 0)
    generate(*source, plots=True)
    assert (source[-1] / "traits/category/all_category1/transfer_tree.pdf").is_file()


def test_changed_input_preserves_previous_bundle(source, monkeypatch):
    import focus_hgt_traits
    generate(*source, plots=False)
    previous = (source[-1] / "manifest.json").read_bytes()
    build = focus_hgt_traits.build_focus
    def change_source(*args, **kwargs):
        result = build(*args, **kwargs)
        with source[0].open("a") as handle:
            handle.write("\n")
        return result
    monkeypatch.setattr(focus_hgt_traits, "build_focus", change_source)
    with pytest.raises(ValueError, match="inputs changed"):
        generate(*source, plots=False)
    assert (source[-1] / "manifest.json").read_bytes() == previous


def test_failed_publication_restores_previous_bundle(source, monkeypatch):
    import focus_hgt_traits
    generate(*source, plots=False)
    previous = (source[-1] / "manifest.json").read_bytes()
    replace = focus_hgt_traits.os.replace
    def fail_stage(path, destination):
        if Path(path).name.startswith(".hgt-trait-focus-") and not Path(path).name.startswith(".hgt-trait-focus-backup-"):
            raise OSError("publication failed")
        return replace(path, destination)
    monkeypatch.setattr(focus_hgt_traits.os, "replace", fail_stage)
    with pytest.raises(OSError, match="publication failed"):
        generate(*source, plots=False)
    assert (source[-1] / "manifest.json").read_bytes() == previous
