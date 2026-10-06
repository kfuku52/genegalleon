import json
import sys
from pathlib import Path

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT))

from focus_hgt_traits import generate, read_tsv, write_tsv  # noqa: E402
from species_trait_schema import schema_path, schema_payload  # noqa: E402


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
    assert alphas and set(alphas) == {0.4}
    root = source[-1] / "traits/binary/all_category1"
    _, edges = read_tsv(root / "transfer_edges.tsv")
    assert sum(int(row["hgt_event_count"]) for row in edges) == 3
    assert (root / "transfer_tree.pdf").read_bytes().startswith(b"%PDF")
    from plot_hgt_summary import read_transfer_traits
    assert list(read_transfer_traits(source[3]).columns) == ["binary", "continuous"]


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
