"""Bounded input and evidence regressions, independent of the full WGD scan."""

import copy
import csv
import json
from io import StringIO
from pathlib import Path
from types import SimpleNamespace

import pytest

from workflow.support import wgd_ssd
from workflow.support.wgd_evidence import branch_for_ks, combine_node, read_table, summarize_events


def write_rows(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def candidate(taxa="A,B"):
    return {"species_event_id": "event_AB", "descendant_taxa": taxa,
            "count_support": "count_supported_conditional", "nuisance_bound_reached": "False",
            "background_nuisance_bound_reached": "False", "branch_burst_nuisance_bound_reached": "False"}


def anchor_rows(species=("A", "B"), blocks=3):
    return [{"species_event_id": "event_AB", "species": name, "block_id": f"block{i}",
             "gene_a": f"{name}_a{i}", "gene_b": f"{name}_b{i}", "ks": "0.5", "ks_status": "ok",
             "placement_status": "interval_supported"}
            for name in species for i in range(blocks)]


def pair_positions(tandem=False):
    return {"A_a": {"species": "A", "seqid": "chr1", "rank": 1, "locus_id": "locus_a", "start": 0, "end": 80},
            "A_b": {"species": "A", "seqid": "chr1" if tandem else "chr2", "rank": 2,
                    "locus_id": "locus_b", "start": 100, "end": 180}}


def anchor_positions(rows):
    return {gene: {"species": row["species"], "seqid": arm, "rank": i + 1, "locus_id": gene,
                   "start": i * 100, "end": i * 100 + 80}
            for i, row in enumerate(rows) for arm in ("gene_a", "gene_b") for gene in (row[arm],)}


@pytest.mark.parametrize("text", [
    "\n", "species\tspecies\nA\tA\n", "species\tcount\nA\n", "species\nA\t1\n",
])
def test_tsv_rejects_ambiguous_headers_and_row_widths(tmp_path, text):
    path = tmp_path / "bad.tsv"
    path.write_text(text, encoding="utf-8")
    with pytest.raises(ValueError, match="Invalid TSV"):
        read_table(path)


def test_event_support_is_an_evidence_rule_not_a_probability():
    anchors = anchor_rows()
    row = summarize_events([candidate()], anchors, {name: {"num_annotated_loci": 20} for name in ("A", "B")},
                           positions=anchor_positions(anchors))[0]
    assert row["event_support"] == "WGD-supported"
    assert row["support_meaning"] == "experimental_evidence_rule_not_posterior"
    assert row["synteny_coverage_definition"] == "unique_branch_matched_anchor_loci_over_annotated_loci"
    assert not {"posterior", "probability", "ploidy"}.intersection(row)


def test_event_summary_does_not_join_unrelated_species_to_candidate_branch():
    row = summarize_events([candidate()], anchor_rows(("C", "D")),
                           {name: {"num_annotated_loci": 20} for name in ("C", "D")}, positions=anchor_positions(anchor_rows(("C", "D"))))[0]
    assert row["event_support"] == "unresolved"
    assert row["num_branch_matched_synteny_species"] == 0
    assert row["branch_matched_synteny_species"] == ""


def test_an_outside_species_cannot_supply_the_second_descendant_replication():
    row = summarize_events([candidate()], anchor_rows(("A", "C")),
                           {name: {"num_annotated_loci": 20} for name in ("A", "C")}, positions=anchor_positions(anchor_rows(("A", "C"))))[0]
    assert row["event_support"] == "unresolved"
    assert row["branch_matched_synteny_species"] == "A"


@pytest.mark.parametrize("summary", [None, {}, {"num_genes": 0}])
def test_missing_or_zero_prepared_gene_denominator_cannot_support_wgd(summary):
    summaries = {} if summary is None else {"A": summary}
    row = summarize_events([candidate("A")], anchor_rows(("A",)), summaries, positions=anchor_positions(anchor_rows(("A",))))[0]
    assert row["event_support"] == "unresolved"
    assert row["num_branch_matched_synteny_species"] == 0


def test_repeated_anchor_rows_do_not_inflate_blocks_or_gene_coverage():
    anchors = anchor_rows(("A",), blocks=1) * 30
    row = summarize_events([candidate("A")], anchors, {"A": {"num_annotated_loci": 20}}, positions=anchor_positions(anchors))[0]
    assert row["event_support"] == "unresolved"
    assert row["num_branch_matched_synteny_species"] == 0


def test_block_span_coverage_cannot_replace_unique_anchor_gene_coverage():
    row = summarize_events([candidate("A")], anchor_rows(("A",)),
                           {"A": {"num_annotated_loci": 100, "gene_span_coverage": 1.0}}, positions=anchor_positions(anchor_rows(("A",))))[0]
    assert row["event_support"] == "unresolved"


@pytest.mark.parametrize("field,value", [
    ("count_support", "not_calibrated"), ("count_support", "branch_burst_preferred"),
    ("nuisance_bound_reached", "True"), ("nuisance_bound_reached", "1"),
    ("background_nuisance_bound_reached", "True"), ("background_nuisance_bound_reached", "1"),
    ("background_nuisance_bound_reached", "NA"),
    ("branch_burst_nuisance_bound_reached", "True"), ("branch_burst_nuisance_bound_reached", "1"),
    ("branch_burst_nuisance_bound_reached", "NA"),
])
def test_uncalibrated_burst_or_boundary_counts_cannot_support_wgd(field, value):
    row = summarize_events([{**candidate("A"), field: value}], anchor_rows(("A",)),
                           {"A": {"num_annotated_loci": 20}}, positions=anchor_positions(anchor_rows(("A",))))[0]
    assert row["event_support"] == "unresolved"


@pytest.mark.parametrize("field", ["background_nuisance_bound_reached", "branch_burst_nuisance_bound_reached"])
def test_missing_competing_fit_diagnostics_cannot_support_wgd(field):
    value = candidate("A")
    value.pop(field)
    anchors = anchor_rows(("A",))
    row = summarize_events([value], anchors, {"A": {"num_annotated_loci": 20}},
                           positions=anchor_positions(anchors))[0]
    assert row["event_support"] == "unresolved"


def test_node_requires_the_same_species_branch_even_with_positive_wgd_events():
    row = {**anchor_rows(("A",), 1)[0], "gene_a": "A_a", "gene_b": "A_b"}
    anchors = {("A_a", "A_b"): [row]}
    events = {"event_AB": {"event_support": "WGD-supported"}}
    assert combine_node("event_A", [("A_a", "A_b")], pair_positions(), anchors, events)[0] == "unresolved"
    assert combine_node("event_AB", [("A_a", "A_b")], pair_positions(), anchors, events)[0] == "WGD-supported"


@pytest.mark.parametrize("status", ["boundary_overlap", "missing_ks"])
def test_node_cannot_trust_branch_id_when_anchor_age_placement_is_unresolved(status):
    row = {**anchor_rows(("A",), 1)[0], "gene_a": "A_a", "gene_b": "A_b", "placement_status": status}
    result = combine_node("event_AB", [("A_a", "A_b")], pair_positions(),
                          {("A_a", "A_b"): [row]}, {"event_AB": {"event_support": "WGD-supported"}})
    assert result[0] == "unresolved"


def test_node_anchor_species_must_match_the_positioned_pair():
    row = {**anchor_rows(("B",), 1)[0], "gene_a": "A_a", "gene_b": "A_b"}
    result = combine_node("event_AB", [("A_a", "A_b")], pair_positions(),
                          {("A_a", "A_b"): [row]}, {"event_AB": {"event_support": "WGD-supported"}})
    assert result[0] == "unresolved"


def test_terminal_tandem_and_matched_wgd_anchor_conflict_is_not_resolved_by_precedence():
    row = {**anchor_rows(("A",), 1)[0], "gene_a": "A_a", "gene_b": "A_b"}
    result = combine_node("event_AB", [("A_a", "A_b")], pair_positions(tandem=True),
                          {("A_a", "A_b"): [row]}, {"event_AB": {"event_support": "WGD-supported"}})
    assert result[:2] == ("unresolved", "tandem_and_collinearity_conflict")


def test_unknown_positions_and_empty_pairs_are_not_ssd_evidence():
    assert combine_node("event_AB", [], {}, {}, {})[0] == "unresolved"
    assert combine_node("event_AB", [("missing_a", "missing_b")], {}, {}, {})[0] == "unresolved"


@pytest.fixture
def boundary_tree(tmp_path):
    from nwkit.clade_index import CladeIndex

    path = tmp_path / "species.nwk"
    path.write_text("(((A:1,B:1):1,C:2):1,D:3);\n", encoding="utf-8")
    tree = wgd_ssd.species_tree(path)
    index = CladeIndex(tree)
    leaf = next(node for node in tree.leaves() if node.name == "A")
    ids = {"A": index.clade_id_for_node(leaf), "AB": index.clade_id_for_node(leaf.up),
           "ABC": index.clade_id_for_node(leaf.up.up)}
    bounds = {("A", ids[name]): {"status": "ok", "interval_status": "ok", "ci_lower": low,
                                 "ci_upper": high, "monotone_from_younger_node": monotone}
              for name, low, high, monotone in (("AB", 0.9, 1.1, "not_comparable"), ("ABC", 1.9, 2.1, "yes"))}
    return tree, ids, bounds


def test_ks_boundaries_map_to_child_branches_not_the_boundary_node(boundary_tree):
    tree, ids, bounds = boundary_tree
    assert branch_for_ks("A", (0.4, 0.6), tree, bounds) == (ids["A"], "interval_supported")
    assert branch_for_ks("A", (1.4, 1.6), tree, bounds) == (ids["AB"], "interval_supported")


@pytest.mark.parametrize("interval", [(0.8, 1.0), (1.0, 1.2), (0.0, 0.2)])
def test_ks_interval_touching_a_boundary_is_unresolved(boundary_tree, interval):
    tree, _, bounds = boundary_tree
    assert branch_for_ks("A", interval, tree, bounds)[0] is None


def test_another_focals_boundary_cannot_be_reused(boundary_tree):
    tree, _, bounds = boundary_tree
    other_focal_bounds = {("B", key[1]): row for key, row in bounds.items()}
    assert branch_for_ks("A", 0.5, tree, other_focal_bounds)[0] is None


def test_reversed_boundary_confidence_limits_cannot_place_an_anchor(boundary_tree):
    tree, ids, bounds = boundary_tree
    bounds[("A", ids["AB"])] = {**bounds[("A", ids["AB"])], "ci_upper": 0.1}
    assert branch_for_ks("A", 0.5, tree, bounds)[0] is None


@pytest.fixture
def count_plan(tmp_path):
    tree = tmp_path / "species.nwk"
    tree.write_text("((A:1,B:1):1,C:2);\n", encoding="utf-8")
    counts = tmp_path / "counts.tsv"
    output = tmp_path / "normalized"
    output.mkdir()
    return {"counts": str(counts), "species": ["A", "B", "C"], "species_tree": str(tree),
            "parameters": {"max_count_families": 100, "seed": 1}}, output


@pytest.mark.parametrize("missing", ["NA", "nan", "", "."])
def test_count_normalization_preserves_unknown_counts_and_observed_zero(count_plan, missing):
    plan, output = count_plan
    write_rows(Path(plan["counts"]), ["Orthogroup", "A", "B", "C"], [
        {"Orthogroup": "OG1", "A": missing, "B": "1", "C": "1"},
        {"Orthogroup": "OG2", "A": "1", "B": "0", "C": "1"},
        {"Orthogroup": "absent", "A": "0", "B": "0", "C": "1"},
    ])
    wgd_ssd.normalized_counts(plan, output)
    rows = {row["family_id"]: row for row in read_table(output / "counts.tsv")}
    assert rows["OG1"]["A"] == "NA"
    assert rows["OG2"]["B"] == "0"
    assert "absent" not in rows


@pytest.mark.parametrize("invalid", ["-1", "1.5", "inf", "not_a_count"])
def test_invalid_counts_are_not_silently_converted_to_missing(count_plan, invalid):
    plan, output = count_plan
    write_rows(Path(plan["counts"]), ["Orthogroup", "A", "B", "C"], [
        {"Orthogroup": "OG1", "A": invalid, "B": "1", "C": "1"},
        {"Orthogroup": "OG2", "A": "1", "B": "1", "C": "1"},
    ])
    with pytest.raises(ValueError, match="count"):
        wgd_ssd.normalized_counts(plan, output)


def test_filtered_sequence_set_cannot_create_false_terminal_tandem(tmp_path, monkeypatch):
    fasta, gff = tmp_path / "A.fa", tmp_path / "A.gff3"
    fasta.write_text(">g1\nMALW\n>g3\nMALW\n", encoding="utf-8")
    gff.write_text("##gff-version 3\n"
                   "chr1\ttest\tgene\t1\t90\t.\t+\t.\tID=g1\n"
                   "chr1\ttest\tgene\t101\t190\t.\t+\t.\tID=g2\n"
                   "chr1\ttest\tgene\t201\t290\t.\t+\t.\tID=g3\n", encoding="utf-8")
    source = {"species": "A", "mode": "protein", "fasta": str(fasta), "gff": str(gff), "cds": "",
              "genetic_code": 1, "feature": "gene", "attribute": "ID"}

    def empty_self_scan(_argv, directory, label):
        if label == "align":
            (directory / "self.self.last.inverse.filtered").write_text("", encoding="utf-8")
        else:
            (directory / "self.self.lifted.anchors").write_text("", encoding="utf-8")
            (directory / "self.self.raw.summary.json").write_text(json.dumps({"num_genes": 2}), encoding="utf-8")

    monkeypatch.setattr(wgd_ssd, "run_tool", empty_self_scan)
    output = tmp_path / "evidence"
    output.mkdir()
    plan = {"genomes": [source], "parameters": {"minimum_mapping_fraction": 1, "cscore": 0.7,
                                               "self_hit_percent": 98, "diagonal_bound": 300}}
    positions, _, _, _ = wgd_ssd.prepare_sources(plan, output, cpus=1)
    positions = {row["gene_id"]: row for row in positions}
    assert {"A_g1", "A_g3"}.issubset(positions)
    result = combine_node("event_A", [("A_g1", "A_g3")], positions, {}, {})
    assert result[0] == "unresolved", "The annotation contains intervening locus g2 even though it has no supplied sequence"


def plan_args(workspace):
    values = {"workspace": workspace, "outfile": workspace / "plan.json", "species_tree": "species.nwk", "counts": "counts.tsv",
              "members": "members.tsv", "genomes": "", "sequence_mode": "protein", "genetic_code": 1,
              "count_bootstrap": 0, "ks_bootstrap": 20, "seed": 1, "max_pairs": 100,
              "max_ks_families": 100, "max_count_families": 100, "diagonal_bound": 300, "min_blocks": 3,
              "max_states": 128, "max_iterations": 20, "multiplicity": 2, "cscore": 0.7,
              "self_hit_percent": 98, "minimum_mapping_fraction": 1, "alpha": 0.05, "min_coverage": 0.2}
    return SimpleNamespace(**values)


def test_previously_absent_genetic_code_input_is_part_of_the_plan_contract(tmp_path, monkeypatch):
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    (workspace / "species.nwk").write_text("(A:1,B:1);\n", encoding="utf-8")
    (workspace / "counts.tsv").write_text("Orthogroup\tA\tB\nOG1\t1\t1\nOG2\t1\t1\n", encoding="utf-8")
    (workspace / "members.tsv").write_text("Orthogroup\tA\tB\nOG1\tA_g1\tB_g1\n", encoding="utf-8")
    for name in ("A", "B"):
        for directory, suffix, text in (("species_protein", ".fa", ">g1\nMALW\n"),
                                        ("species_gff", ".gff3", "##gff-version 3\n")):
            root = workspace / "input" / directory
            root.mkdir(parents=True, exist_ok=True)
            (root / f"{name}{suffix}").write_text(text, encoding="utf-8")
    write_rows(workspace / "input/wgd_genomes.tsv", ["species", "mode", "fasta", "gff"], [
        {"species": name, "mode": "protein", "fasta": f"input/species_protein/{name}.fa",
         "gff": f"input/species_gff/{name}.gff3"} for name in ("A", "B")])
    monkeypatch.setattr(wgd_ssd, "owned_identity", lambda: {"source_hashes": {}, "ds": {"source_hashes": {}}})
    plan = wgd_ssd.make_plan(plan_args(workspace))
    wgd_ssd.verify(plan)
    code = workspace / "input/species_genetic_code/species_genetic_code.tsv"
    code.parent.mkdir(parents=True)
    code.write_text("species\tgenetic_code\nA\t4\nB\t1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="appeared|changed"):
        wgd_ssd.verify(plan)


def test_verify_rejects_changed_hashed_input(tmp_path):
    path = tmp_path / "input.tsv"
    path.write_text("original\n", encoding="utf-8")
    plan = {"input_hashes": {str(path): wgd_ssd.digest(path)}, "absent_inputs": []}
    path.write_text("changed\n", encoding="utf-8")
    with pytest.raises(ValueError, match="input changed"):
        wgd_ssd.verify(plan)


def test_classification_rejects_evidence_whose_genomic_input_has_changed(tmp_path, monkeypatch):
    from nwkit.clade_index import CladeIndex

    species_tree, gene_tree = tmp_path / "species.nwk", tmp_path / "gene.nwk"
    species_tree.write_text("(A:1,B:1);\n", encoding="utf-8")
    gene_tree.write_text("(A_a:1,A_b:1);\n", encoding="utf-8")
    genomic_input = tmp_path / "A.gff3"
    genomic_input.write_text("original annotation\n", encoding="utf-8")
    plan = {"species_tree": str(species_tree), "input_hashes": {
        str(species_tree): wgd_ssd.digest(species_tree), str(genomic_input): wgd_ssd.digest(genomic_input)},
        "absent_inputs": []}
    evidence = tmp_path / "evidence"
    evidence.mkdir()
    (evidence / "summary.json").write_text(json.dumps({"plan": plan}), encoding="utf-8")
    positions = [{"gene_id": name, **data} for name, data in pair_positions(tandem=True).items()]
    write_rows(evidence / "gene_positions.tsv", ["gene_id", "species", "seqid", "rank", "locus_id", "start", "end"], positions)
    write_rows(evidence / "wgd_events.tsv", ["species_event_id", "event_support"], [])
    write_rows(evidence / "anchor_evidence.tsv", ["gene_a", "gene_b", "species_event_id"], [])
    genomic_input.write_text("new annotation with an intervening gene\n", encoding="utf-8")

    def reconciliation_stub(_argv, directory, label):
        assert label == "reconcile"
        tree = wgd_ssd.species_tree(gene_tree)
        fields = ["event_type", "gene_clade_id", "species_event_id", "mapping_status", "event_status"]
        write_rows(directory / "reconciliation.tsv", fields, [{
            "event_type": "duplication", "gene_clade_id": CladeIndex(tree).clade_id_for_node(tree),
            "species_event_id": "event_A", "mapping_status": "mapped", "event_status": "resolved"}])

    monkeypatch.setattr(wgd_ssd, "run_tool", reconciliation_stub)
    args = SimpleNamespace(evidence=evidence, output=tmp_path / "classification", species_tree=species_tree,
                           gene_tree=gene_tree, family_id="OG1", species_parser="taxonomic", species_regex="",
                           species_map="", proximal_distance=10)
    with pytest.raises(ValueError, match="input changed"):
        wgd_ssd.classify(args)


@pytest.mark.parametrize("defect", ["count", "duplicate_gene", "duplicate_family", "missing_family", "extra_species"])
def test_membership_contract_rejects_inconsistent_or_ambiguous_data(tmp_path, defect):
    counts, members = tmp_path / "counts.tsv", tmp_path / "members.tsv"
    write_rows(counts, ["Orthogroup", "A", "B"], [
        {"Orthogroup": family, "A": 1, "B": 1} for family in ("OG1", "OG2")])
    rows = [{"Orthogroup": "OG1", "A": "A_a", "B": "B_a"},
            {"Orthogroup": "OG2", "A": "A_b", "B": "B_b"}]
    fields = ["Orthogroup", "A", "B"]
    if defect == "count":
        rows[0]["A"] = "A_a, A_c"
    elif defect == "duplicate_gene":
        rows[1]["A"] = "A_a"
    elif defect == "duplicate_family":
        rows[1]["Orthogroup"] = "OG1"
    elif defect == "missing_family":
        rows.pop()
    else:
        fields.append("C")
    write_rows(members, fields, rows)
    with pytest.raises(ValueError):
        wgd_ssd.validated_members({"members": str(members), "counts": str(counts), "species": ["A", "B"]})


def test_output_protection_preserves_curated_inputs_and_input_ancestors(tmp_path):
    workspace = tmp_path / "workspace"
    source = workspace / "input/A.gff"
    for output in (source, workspace, workspace / "input/new.json"):
        with pytest.raises(ValueError, match="output"):
            wgd_ssd.protect_output(output, [source], workspace)
    wgd_ssd.protect_output(workspace / "output/new", [source], workspace)


def test_classification_rejects_edited_evidence_table_before_running_tools(tmp_path, monkeypatch):
    species = tmp_path / "species.nwk"
    species.write_text("(A:1,B:1);\n")
    evidence = tmp_path / "evidence"
    evidence.mkdir()
    table = evidence / "gene_positions.tsv"
    table.write_text("original\n")
    for name in ("wgd_events.tsv", "anchor_evidence.tsv"):
        (evidence / name).write_text("empty\n")
    summary = {"plan": {"species_tree": str(species), "input_hashes": {str(species): wgd_ssd.digest(species)}},
               "output_hashes": {path.name: wgd_ssd.digest(path) for path in evidence.iterdir()}}
    (evidence / "summary.json").write_text(json.dumps(summary))
    table.write_text("changed\n")
    monkeypatch.setattr(wgd_ssd, "run_tool", lambda *_: pytest.fail("stale evidence must fail before execution"))
    with pytest.raises(ValueError, match="evidence output changed"):
        wgd_ssd.classify(SimpleNamespace(evidence=evidence, species_tree=species, output=tmp_path / "results", native_tree_likelihood=0))


def test_long_labels_and_deep_topology_plots(tmp_path):
    from ete4 import Tree

    rows = [{"branch_id": i, "descendant_taxa": "Alpha_species_with_long_name,Beta_species_with_long_name,Gamma_species",
             "likelihood_ratio_statistic": i * 1.5, "event_support": "WGD-supported" if i == 2 else "unresolved"}
            for i in range(1, 7)]
    wgd_ssd.plot_candidates(rows, tmp_path)
    names = [f"Alpha_species_with_long_name_gene_identifier_{i:03d}_additional_description" for i in range(32)]
    topology = names[0]
    for name in names[1:]:
        topology = f"({topology},{name})"
    tree = Tree(topology + ";", parser=1)
    for index, node in enumerate(tree.traverse()):
        if node.children:
            node.add_prop("duplication_origin", ("WGD-supported", "SSD-supported", "unresolved")[index % 3])
    wgd_ssd.plot_origins(tree, tmp_path)
    for name in ("wgd_candidates.pdf", "duplication_origins.pdf"):
        assert (tmp_path / name).stat().st_size > 1000


def test_missing_gene_support_is_na_but_explicit_zero_and_one_are_preserved(tmp_path):
    species, gene = tmp_path / "species.nwk", tmp_path / "gene.nwk"
    species.write_text("(A:1,B:1);\n")
    gene.write_text("((A_a:1,A_b:1)1:1,(A_c:1,A_d:1)0:1);\n")
    evidence = tmp_path / "evidence"
    evidence.mkdir()
    positions = [{"gene_id": name, "species": "A", "locus_id": name, "seqid": "chr1", "rank": rank,
                  "start": rank * 100, "end": rank * 100 + 80}
                 for name, rank in zip(("A_a", "A_b", "A_c", "A_d"), (1, 2, 4, 5), strict=True)]
    write_rows(evidence / "gene_positions.tsv", ["gene_id", "species", "locus_id", "seqid", "rank", "start", "end"], positions)
    write_rows(evidence / "anchor_evidence.tsv", ["gene_a", "gene_b", "species_event_id"], [])
    write_rows(evidence / "wgd_events.tsv", ["species_event_id", "event_support"], [])
    summary = {"plan": {"species_tree": str(species), "input_hashes": {str(species): wgd_ssd.digest(species)}},
               "output_hashes": {path.name: wgd_ssd.digest(path) for path in evidence.iterdir()}}
    (evidence / "summary.json").write_text(json.dumps(summary))
    output = tmp_path / "classification"
    wgd_ssd.classify(SimpleNamespace(evidence=evidence, species_tree=species, gene_tree=gene, output=output,
                                   family_id="OGsupport", species_parser="legacy", species_regex="^([AB])_",
                                   species_map="", proximal_distance=10, native_tree_likelihood=0))
    support = [row["input_tree_support"] for row in read_table(output / "duplication_origins.tsv")]
    assert support.count("NA") == 1
    assert {float(value) for value in support if value != "NA"} == {0, 1}


def test_duplicate_count_families_are_not_overwritten_by_membership_validation(tmp_path):
    counts, members = tmp_path / "counts.tsv", tmp_path / "members.tsv"
    write_rows(counts, ["family_id", "A", "B"], [{"family_id": "OG1", "A": 1, "B": 1}] * 2)
    write_rows(members, ["family_id", "A", "B"], [{"family_id": "OG1", "A": "A_a", "B": "B_a"}])
    with pytest.raises(ValueError, match="Duplicate"):
        wgd_ssd.validated_members({"members": str(members), "counts": str(counts), "species": ["A", "B"]})


@pytest.mark.parametrize("gene", ["B_wrong", "unprefixed", "A_", "A_gene|other"])
def test_membership_ids_must_belong_to_the_declared_species(tmp_path, gene):
    counts, members = tmp_path / "counts.tsv", tmp_path / "members.tsv"
    write_rows(counts, ["family_id", "A", "B"], [{"family_id": "OG1", "A": 1, "B": 1}])
    write_rows(members, ["family_id", "A", "B"], [{"family_id": "OG1", "A": gene, "B": "B_a"}])
    with pytest.raises(ValueError, match="identifier|prefix"):
        wgd_ssd.validated_members({"members": str(members), "counts": str(counts), "species": ["A", "B"]})


def test_output_cannot_be_nested_inside_the_evidence_input_directory(tmp_path):
    evidence = tmp_path / "evidence"
    evidence.mkdir()
    with pytest.raises(ValueError, match="output"):
        wgd_ssd.protect_output(evidence / "classification", [evidence])


def test_ks_numeric_strings_are_normalized_before_interval_comparison(boundary_tree):
    tree, ids, bounds = boundary_tree
    assert branch_for_ks("A", ("1.4", "1.6"), tree, bounds) == (ids["AB"], "interval_supported")


def classification_fixture(tmp_path, gene_text, positions, anchors=(), events=()):
    species, gene = tmp_path / "species.nwk", tmp_path / "gene.nwk"
    species.write_text("(((A:1,B:1)AB:1,C:2)ABC:1,D:3)ROOT;\n")
    gene.write_text(gene_text + "\n")
    evidence = tmp_path / "evidence"
    evidence.mkdir()
    write_rows(evidence / "gene_positions.tsv", ["gene_id", "species", "locus_id", "seqid", "rank", "start", "end"], positions)
    write_rows(evidence / "anchor_evidence.tsv", ["species", "block_id", "gene_a", "gene_b", "ks", "ks_status", "species_event_id", "placement_status"], anchors)
    write_rows(evidence / "wgd_events.tsv", ["species_event_id", "event_support"], events)
    summary = {"plan": {"species_tree": str(species), "input_hashes": {str(species): wgd_ssd.digest(species)}},
               "output_hashes": {path.name: wgd_ssd.digest(path) for path in evidence.iterdir()}}
    (evidence / "summary.json").write_text(json.dumps(summary))
    return SimpleNamespace(evidence=evidence, species_tree=species, gene_tree=gene, output=tmp_path / "classification",
                           family_id="OGedge", species_parser="legacy", species_regex="^([ABCD])_",
                           species_map="", proximal_distance=10, native_tree_likelihood=0)


def positioned(name, chromosome="chr1", rank=1, species=None, locus=None):
    return {"gene_id": name, "species": species or name[0], "locus_id": locus or name, "seqid": chromosome,
            "rank": rank, "start": rank * 100, "end": rank * 100 + 80}


UNROOTED_SPECIES_INPUTS = (
    "[&U]((A:1,B:1):1,C:2);",
    "((A:1,B:1):1,C:2)[&&NHX:nwkit_rooted=no];",
    "(A:1,B:1,C:1);",
    "(A:1,B:1,C:1)[&&NHX:nwkit_rooted=unknown];",
)


@pytest.mark.parametrize("text", UNROOTED_SPECIES_INPUTS)
def test_species_reader_rejects_unrooted_or_unknown_orientation(tmp_path, text):
    path = tmp_path / "species.nwk"
    path.write_text(text)
    with pytest.raises(ValueError, match="rooted species tree"):
        wgd_ssd.species_tree(path)


@pytest.mark.parametrize("text", (
    "((A:1,B:1):1,C:2);",
    "[&R](A:1,B:1,C:1);",
    "(A:1,B:1,C:1)[&&NHX:nwkit_rooted=yes];",
))
def test_declared_root_polytomies_and_legacy_binary_roots_keep_orientation(tmp_path, text):
    from nwkit.rooting_state import is_rooted

    path = tmp_path / "species.nwk"
    path.write_text(text)
    tree = wgd_ssd.species_tree(path)
    assert is_rooted(tree)
    assert len(tree.children) == (2 if text.startswith("((") else 3)


@pytest.mark.parametrize("text", UNROOTED_SPECIES_INPUTS)
def test_count_normalization_rejects_unrooted_tree_without_writing(count_plan, text):
    plan, output = count_plan
    Path(plan["species_tree"]).write_text(text)
    write_rows(Path(plan["counts"]), ["family_id", "A", "B", "C"], [
        {"family_id": family, "A": 1, "B": 1, "C": 1} for family in ("OG1", "OG2")])
    with pytest.raises(ValueError, match="rooted species tree"):
        wgd_ssd.normalized_counts(plan, output)
    assert not list(output.iterdir())


@pytest.mark.parametrize("text", (
    "[&R](A:1,B:1,C:1);",
    "(A:1,B:1,C:1)[&&NHX:nwkit_rooted=yes];",
))
def test_count_normalization_uses_all_declared_root_polytomy_clades(count_plan, text):
    plan, output = count_plan
    Path(plan["species_tree"]).write_text(text)
    write_rows(Path(plan["counts"]), ["family_id", "A", "B", "C"], [
        {"family_id": family, "A": 1, "B": 1, "C": 1} for family in ("OG1", "OG2")]
        + [{"family_id": "absent_C", "A": 1, "B": 1, "C": 0}])
    wgd_ssd.normalized_counts(plan, output)
    assert [row["family_id"] for row in read_table(output / "counts.tsv")] == ["OG1", "OG2"]
    assert read_table(output / "count_family_selection.tsv")[-1]["selection"] == "excluded_root_clade_absence"


@pytest.mark.parametrize("text", UNROOTED_SPECIES_INPUTS)
def test_plan_rejects_unrooted_species_tree_before_tool_identity(tmp_path, monkeypatch, text):
    path = tmp_path / "species.nwk"
    path.write_text(text)
    counts, members = tmp_path / "counts.tsv", tmp_path / "members.tsv"
    write_rows(counts, ["family_id", "A", "B", "C"], [])
    write_rows(members, ["family_id", "A", "B", "C"], [])
    args = SimpleNamespace(workspace=tmp_path, species_tree=str(path), counts=str(counts), members=str(members),
                           genomes="", sequence_mode="cds", genetic_code=1, count_bootstrap=199, ks_bootstrap=199,
                           seed=1, max_pairs=5000, max_ks_families=500, max_count_families=10000,
                           diagonal_bound=300, cscore=0.7, self_hit_percent=98, minimum_mapping_fraction=1,
                           alpha=0.05, min_coverage=0.2, min_blocks=3, max_states=256, max_iterations=200,
                           multiplicity=2, outfile=tmp_path / "plan.json")
    monkeypatch.setattr(wgd_ssd, "owned_identity", lambda: pytest.fail("invalid rooting must fail before tools"))
    with pytest.raises(ValueError, match="rooted species tree"):
        wgd_ssd.make_plan(args)


@pytest.mark.parametrize("text", (
    "[&U](A_a:1,A_b:1);",
    "(A_a:1,A_b:1)[&&NHX:nwkit_rooted=no];",
    "(A_a:1,A_b:1,A_c:1);",
    "(A_a:1,A_b:1)[&&NHX:nwkit_rooted=unknown];",
))
def test_classification_rejects_unrooted_gene_before_outputs(tmp_path, monkeypatch, text):
    args = classification_fixture(tmp_path, text, [])
    monkeypatch.setattr(wgd_ssd, "run_tool", lambda *_: pytest.fail("invalid rooting must fail before tools"))
    with pytest.raises(ValueError, match="rooted gene tree"):
        wgd_ssd.classify(args)
    assert not args.output.exists()


@pytest.mark.parametrize("text", UNROOTED_SPECIES_INPUTS)
def test_classification_rejects_unrooted_species_before_outputs(tmp_path, monkeypatch, text):
    args = classification_fixture(tmp_path, "(A_a:1,A_b:1);", [])
    args.species_tree.write_text(text)
    summary_path = args.evidence / "summary.json"
    summary = json.loads(summary_path.read_text())
    summary["plan"]["input_hashes"][str(args.species_tree)] = wgd_ssd.digest(args.species_tree)
    summary_path.write_text(json.dumps(summary))
    monkeypatch.setattr(wgd_ssd, "run_tool", lambda *_: pytest.fail("invalid rooting must fail before tools"))
    with pytest.raises(ValueError, match="rooted species tree"):
        wgd_ssd.classify(args)
    assert not args.output.exists()


@pytest.mark.parametrize("species_text", (
    "[&R]((A:1,B:1):1,C:2);",
    "((A:1,B:1):1,C:2)[&&NHX:nwkit_rooted=yes];",
))
@pytest.mark.parametrize("gene_text", (
    "[&R](A_a:1,A_b:1);",
    "(A_a:1,A_b:1)[&&NHX:nwkit_rooted=yes];",
))
def test_classification_accepts_explicit_roots_without_reorientation(tmp_path, species_text, gene_text):
    args = classification_fixture(tmp_path, gene_text, [positioned("A_a"), positioned("A_b", rank=2)])
    args.species_tree.write_text(species_text)
    summary_path = args.evidence / "summary.json"
    summary = json.loads(summary_path.read_text())
    summary["plan"]["input_hashes"][str(args.species_tree)] = wgd_ssd.digest(args.species_tree)
    summary_path.write_text(json.dumps(summary))
    wgd_ssd.classify(args)
    rows = read_table(args.output / "duplication_origins.tsv")
    assert len(rows) == 1 and rows[0]["classification"] == "SSD-supported"
    annotated = (args.output / "classified_gene_tree.nhx").read_text()
    from nwkit.rooting_state import require_rooted
    from nwkit.util import read_tree
    # Rooting tokens are off by default; the rooted topology and annotations
    # must survive the actual writer/reader round trip for both input forms.
    assert not annotated.startswith("[&R]")
    tree = read_tree(str(args.output / "classified_gene_tree.nhx"), "auto", True, rooted="auto")
    require_rooted(tree, "classification round trip")
    assert [(leaf.name, leaf.dist) for leaf in tree.leaves()] == [("A_a", 1.0), ("A_b", 1.0)]
    assert tree.props["duplication_origin"] == "SSD-supported"


@pytest.mark.parametrize("text", (
    "[&R](A:1,B:1,C:1);",
    "(A:1,B:1,C:1)[&&NHX:nwkit_rooted=yes];",
))
def test_declared_species_polytomy_is_not_arbitrarily_resolved_for_reconciliation(tmp_path, text):
    args = classification_fixture(tmp_path, "[&R](A_a:1,A_b:1);", [])
    args.species_tree.write_text(text)
    summary_path = args.evidence / "summary.json"
    summary = json.loads(summary_path.read_text())
    summary["plan"]["input_hashes"][str(args.species_tree)] = wgd_ssd.digest(args.species_tree)
    summary_path.write_text(json.dumps(summary))
    with pytest.raises(RuntimeError, match="strictly bifurcating"):
        wgd_ssd.classify(args)
    assert not (args.output / "summary.json").exists()


def test_species_parser_and_coordinate_species_disagreement_remains_unresolved(tmp_path):
    args = classification_fixture(tmp_path, "(A_a:1,A_b:1);",
                                  [positioned("A_a", species="B"), positioned("A_b", rank=2, species="B")])
    wgd_ssd.classify(args)
    row = read_table(args.output / "duplication_origins.tsv")[0]
    assert row["classification"] == "unresolved"
    assert row["reason"] == "genomic_species_mapping_conflict"


def test_generated_nhx_origins_are_recomputed_not_inherited(tmp_path):
    args = classification_fixture(tmp_path,
        "(A_a:1[&&NHX:duplication_origin=WGD-supported:mul_mapping_status=ambiguous],B_a:1)[&&NHX:D=N:duplication_origin=WGD-supported:conditional_wgd_probability=0.99:mul_mapping_status=consistent:mul_dl_duplication=all:mul_gene_node=42:mul_best_hypotheses=99:mul_optimal_mappings=999:custom=keep];",
        [positioned("A_a"), positioned("B_a")])
    wgd_ssd.classify(args)
    annotated = (args.output / "classified_gene_tree.nhx").read_text()
    assert "duplication_origin=" not in annotated
    assert "conditional_wgd_probability=" not in annotated
    assert "mul_" not in annotated
    assert "custom=keep" in annotated and "D=N" in annotated


def test_classification_rejects_inputs_changed_during_reconciliation(tmp_path, monkeypatch):
    args = classification_fixture(tmp_path, "(A_a:1,A_b:1);", [positioned("A_a"), positioned("A_b", rank=2)])
    real_run = wgd_ssd.run_tool

    def mutate(argv, directory, label):
        real_run(argv, directory, label)
        args.gene_tree.write_text("(A_c:1,A_d:1);\n")

    monkeypatch.setattr(wgd_ssd, "run_tool", mutate)
    with pytest.raises(ValueError, match="changed"):
        wgd_ssd.classify(args)
    assert not (args.output / "summary.json").exists()


def test_nhx_ancestral_branch_mapping_cannot_get_terminal_tandem_label(tmp_path):
    args = classification_fixture(tmp_path, "(A_a:1,A_b:1)[&&NHX:D=Y:S=AB];",
                                  [positioned("A_a"), positioned("A_b", rank=2)])
    wgd_ssd.classify(args)
    assert read_table(args.output / "duplication_origins.tsv")[0]["classification"] == "unresolved"


def test_internal_wgd_pair_must_cross_the_reconciled_duplication_children(tmp_path):
    from nwkit.clade_index import CladeIndex

    tree_file = tmp_path / "reference.nwk"
    tree_file.write_text("(((A:1,B:1)AB:1,C:2)ABC:1,D:3)ROOT;\n")
    tree = wgd_ssd.species_tree(tree_file)
    event = CladeIndex(tree).clade_id_for_node(tree.children[0].children[0])
    rows = [positioned("A_a"), positioned("A_b", chromosome="chr2"), positioned("B_a"), positioned("B_b", chromosome="chr2")]
    anchor = {"species": "A", "block_id": "block1", "gene_a": "A_a", "gene_b": "A_b", "ks": "1.5",
              "ks_status": "ok", "species_event_id": event, "placement_status": "interval_supported"}
    args = classification_fixture(tmp_path, "((A_a:1,B_a:1):1,(A_b:1,B_b:1):1);", rows,
                                  [anchor], [{"species_event_id": event, "event_support": "WGD-supported"}])
    wgd_ssd.classify(args)
    origins = read_table(args.output / "duplication_origins.tsv")
    assert len(origins) == 1 and origins[0]["species_event_id"] == event
    assert origins[0]["classification"] == "WGD-supported"


@pytest.mark.parametrize("defect", ["bad_ks", "same_locus", "redundant_blocks", "unknown_gene", "missing_denominator"])
def test_event_support_rejects_invalid_or_redundant_genomic_evidence(defect):
    anchors = anchor_rows(("A",))
    positions = {gene: positioned(gene, chromosome="chr1" if arm == "gene_a" else "chr2", rank=i + 1)
                 for i, row in enumerate(anchors) for arm in ("gene_a", "gene_b") for gene in (row[arm],)}
    summaries = {"A": {"num_annotated_loci": 20, "num_genes": 6}}
    if defect == "bad_ks":
        for row in anchors:
            row["ks_status"] = "saturated"
    elif defect == "same_locus":
        for row in positions.values():
            row["locus_id"] = "one_locus"
    elif defect == "redundant_blocks":
        anchors = [{**row, "block_id": f"redundant{i}"} for i in range(3) for row in anchors]
    elif defect == "unknown_gene":
        positions = {}
    else:
        summaries["A"].pop("num_annotated_loci")
    row = summarize_events([candidate("A")], anchors, summaries, positions=positions)[0]
    assert row["event_support"] == "unresolved"


def ks_fixture(tmp_path, monkeypatch, counts=(1, 1, 1), ds_defect=None):
    members, count_file = tmp_path / "members.tsv", tmp_path / "counts.tsv"
    genes = {name: [f"{name}_g{i}" for i in range(2 if count == 2 else 1)]
             for name, count in zip(("A", "B", "C"), counts, strict=True)}
    write_rows(members, ["family_id", "A", "B", "C"], [{"family_id": "OG1", **{name: ", ".join(ids) for name, ids in genes.items()}}])
    write_rows(count_file, ["family_id", "A", "B", "C"], [{"family_id": "OG1", **dict(zip(("A", "B", "C"), counts, strict=True))}])
    plan = {"species": ["A", "B", "C"], "members": str(members), "counts": str(count_file),
            "genomes": [{"species": name, "genetic_code": 1} for name in genes],
            "parameters": {"seed": 7, "max_pairs": 2, "max_ks_families": 100}}
    sequences = {name: {gene: ("ATGGCT", "MA") for gene in ids} for name, ids in genes.items()}

    def aligned(pair, _sequences, _code):
        return {"pair_id": "|".join(pair), "target_gene": pair[0], "query_gene": pair[1], "sequence_1": "ATGGCT", "sequence_2": "ATGGCT"}

    def ds(argv, _directory, _label):
        pairs = read_table(Path(argv[argv.index("--pairs_file") + 1]))
        rows = [{"pair_id": pair["pair_id"], "status": "ok", "dS": "0.5"} for pair in pairs]
        if ds_defect == "duplicate":
            rows.append(rows[0])
        elif ds_defect == "missing":
            rows.pop()
        elif ds_defect == "extra":
            rows.append({"pair_id": "unknown|pair", "status": "ok", "dS": "0.5"})
        elif ds_defect == "invalid_ok":
            rows[0]["dS"] = "NA"
        elif ds_defect == "saturated":
            rows[0].update(status="saturated", dS="NA")
        write_rows(Path(argv[argv.index("--out_file") + 1]), ["pair_id", "status", "dS"], rows)

    monkeypatch.setattr(wgd_ssd, "align_pair", aligned)
    monkeypatch.setattr(wgd_ssd, "run_tool", ds)
    return plan, sequences


def test_known_multicopy_family_is_not_a_single_copy_ks_calibration(tmp_path, monkeypatch):
    plan, sequences = ks_fixture(tmp_path, monkeypatch, counts=(2, 1, 1))
    _, contrasts = wgd_ssd.estimate_ks(plan, tmp_path, sequences, [], 1)
    assert contrasts == []
    assert read_table(tmp_path / "ks_family_selection.tsv")[0]["selection"] == "known_multicopy_family"


def test_unknown_count_is_not_assumed_single_copy_for_ks(tmp_path, monkeypatch):
    plan, sequences = ks_fixture(tmp_path, monkeypatch, counts=("NA", 1, 1))
    _, contrasts = wgd_ssd.estimate_ks(plan, tmp_path, sequences, [], 1)
    assert {(row["species_a"], row["species_b"]) for row in contrasts} == {("B", "C")}


@pytest.mark.parametrize("defect", ["duplicate", "missing", "extra", "invalid_ok"])
def test_ds_report_must_cover_requested_pairs_exactly_without_invalid_ok_values(tmp_path, monkeypatch, defect):
    plan, sequences = ks_fixture(tmp_path, monkeypatch, ds_defect=defect)
    with pytest.raises(ValueError, match="dS|Ks"):
        wgd_ssd.estimate_ks(plan, tmp_path, sequences, [], 1)


def test_saturated_ks_pair_is_audited_not_used_as_a_divergence_contrast(tmp_path, monkeypatch):
    plan, sequences = ks_fixture(tmp_path, monkeypatch, ds_defect="saturated")
    reports, contrasts = wgd_ssd.estimate_ks(plan, tmp_path, sequences, [], 1)
    assert reports["A_g0|B_g0"]["status"] == "saturated"
    assert len(contrasts) == 2
    assert read_table(tmp_path / "pair_selection.tsv")[0]["status"] == "saturated"


def test_self_pair_sampling_audit_includes_unsampled_and_unavailable_cds(tmp_path, monkeypatch):
    plan, sequences = ks_fixture(tmp_path, monkeypatch)
    for i in range(5):
        sequences["A"][f"A_self{i}"] = ("ATGGCT", "MA")
    anchors = [{"species": "A", "gene_a": "A_g0", "gene_b": f"A_self{i}"} for i in range(6)]
    wgd_ssd.estimate_ks(plan, tmp_path, sequences, anchors, 1)
    audit = [row for row in read_table(tmp_path / "pair_selection.tsv") if row["pair_role"] == "self_anchor"]
    assert len(audit) == 6
    assert sum(row["status"] == "ok" for row in audit) == 2
    assert sum(row["status"] == "not_sampled" for row in audit) == 3
    assert sum(row["status"] == "cds_unavailable" for row in audit) == 1


def test_selected_isoform_bounds_do_not_hide_full_locus_overlap(tmp_path, monkeypatch):
    fasta, gff = tmp_path / "A.fa", tmp_path / "A.gff3"
    fasta.write_text(">t1\nMAAAAA\n>t1b\nMAA\n>t2\nMAAAA\n")
    gff.write_text("##gff-version 3\n"
                   "chr1\tt\tgene\t1\t300\t.\t+\t.\tID=g1\n"
                   "chr1\tt\tmRNA\t1\t60\t.\t+\t.\tID=t1;Parent=g1\n"
                   "chr1\tt\tmRNA\t151\t300\t.\t+\t.\tID=t1b;Parent=g1\n"
                   "chr1\tt\tgene\t101\t140\t.\t+\t.\tID=g2\n"
                   "chr1\tt\tmRNA\t101\t140\t.\t+\t.\tID=t2;Parent=g2\n")
    source = {"species": "A", "mode": "protein", "fasta": str(fasta), "gff": str(gff), "cds": "",
              "genetic_code": 1, "feature": "mRNA", "attribute": "ID"}

    def empty_scan(_argv, directory, label):
        if label == "align":
            (directory / "self.self.last.inverse.filtered").write_text("")
        else:
            (directory / "self.self.lifted.anchors").write_text("")
            (directory / "self.self.raw.summary.json").write_text(json.dumps({"num_genes": 2}))

    monkeypatch.setattr(wgd_ssd, "run_tool", empty_scan)
    output = tmp_path / "results"
    output.mkdir()
    plan = {"genomes": [source], "parameters": {"minimum_mapping_fraction": 1, "cscore": 0.7,
                                               "self_hit_percent": 98, "diagonal_bound": 300}}
    rows, _, summaries, _ = wgd_ssd.prepare_sources(plan, output, 1)
    positions = {row["gene_id"]: row for row in rows}
    assert set(positions) == {"A_t1", "A_t2"}
    assert (positions["A_t1"]["start"], positions["A_t1"]["end"]) == (0, 300)
    assert summaries["A"]["num_annotated_loci"] == 2
    assert combine_node("event_A", [("A_t1", "A_t2")], positions, {}, {})[:2] == (
        "unresolved", "overlapping_or_ambiguous_loci")


@pytest.mark.parametrize("annotation", ["D=N", "D=Y:H=Y", "D=Y:S=does_not_exist"])
def test_unsupported_or_nonduplication_nhx_events_never_receive_confident_origin_labels(tmp_path, annotation):
    args = classification_fixture(tmp_path, f"(A_a:1,A_b:1)[&&NHX:{annotation}:duplication_origin=WGD-supported];",
                                  [positioned("A_a"), positioned("A_b", rank=2)])
    wgd_ssd.classify(args)
    rows = read_table(args.output / "duplication_origins.tsv")
    if annotation == "D=Y:S=does_not_exist":
        assert len(rows) == 1 and rows[0]["classification"] == "unresolved"
        assert rows[0]["reason"] == "unresolved_reconciliation"
    else:
        assert rows == []
    annotated = (args.output / "classified_gene_tree.nhx").read_text()
    assert "duplication_origin=WGD-supported" not in annotated and "duplication_origin=SSD-supported" not in annotated


@pytest.mark.parametrize("total", ["-1", "2.5", "1", "3", "invalid"])
def test_inconsistent_or_invalid_count_total_is_not_discarded(tmp_path, total):
    path = tmp_path / "counts.tsv"
    write_rows(path, ["family_id", "A", "B", "Total"], [{"family_id": "OG1", "A": 1, "B": 1, "Total": total}])
    with pytest.raises(ValueError, match="Total"):
        wgd_ssd.validated_counts({"species": ["A", "B"], "counts": str(path)})


@pytest.mark.parametrize("method", ["family-bootstrap-percentile", "pair-median-bonferroni"])
def test_ks_placement_uses_primary_bounds_independently_of_ci_method_or_bootstrap_diagnostics(boundary_tree, method):
    tree, ids, bounds = boundary_tree
    bounds = {key: {**row, "ci_method": method, "num_bootstrap": 0,
                    "bootstrap_ci_lower": "0.4", "bootstrap_ci_upper": "0.6"} for key, row in bounds.items()}
    assert branch_for_ks("A", 0.5, tree, bounds) == (ids["A"], "interval_supported")
    bounds[("A", ids["AB"])].update(ci_lower="0.4", ci_upper="0.6")
    assert branch_for_ks("A", 0.5, tree, bounds) == (None, "boundary_overlap")


@pytest.mark.parametrize("taxa", ["", "A,A"])
def test_empty_or_duplicated_candidate_taxa_cannot_support_an_event(taxa):
    anchors = anchor_rows(("A",))
    assert summarize_events([candidate(taxa)], anchors, {"A": {"num_annotated_loci": 20}},
                            positions=anchor_positions(anchors))[0]["event_support"] == "unresolved"


def test_unknown_nuisance_boundary_status_is_not_treated_as_false():
    anchors = anchor_rows(("A",))
    value = candidate("A")
    value.pop("nuisance_bound_reached")
    assert summarize_events([value], anchors, {"A": {"num_annotated_loci": 20}},
                            positions=anchor_positions(anchors))[0]["event_support"] == "unresolved"


@pytest.mark.parametrize("defect", ["duplicate_gene", "duplicate_event"])
def test_ambiguous_evidence_indexes_do_not_silently_overwrite_rows(tmp_path, defect):
    positions = [positioned("A_a"), positioned("A_b", rank=2)]
    events = [{"species_event_id": "event", "event_support": "unresolved"}]
    if defect == "duplicate_gene":
        positions.append(positions[0])
    else:
        events.append(events[0])
    args = classification_fixture(tmp_path, "(A_a:1,A_b:1);", positions, events=events)
    with pytest.raises(ValueError, match="Duplicate"):
        wgd_ssd.classify(args)


def test_isoform_ambiguity_within_one_child_blocks_ancestral_wgd_label(tmp_path):
    from nwkit.clade_index import CladeIndex

    reference = tmp_path / "reference.nwk"
    reference.write_text("(((A:1,B:1)AB:1,C:2)ABC:1,D:3)ROOT;\n")
    tree = wgd_ssd.species_tree(reference)
    event = CladeIndex(tree).clade_id_for_node(next(leaf for leaf in tree.leaves() if leaf.name == "A"))
    positions = [positioned("A_a", locus="g1"), positioned("A_iso", rank=2, locus="g1"),
                 positioned("A_b", chromosome="chr2")]
    anchor = {"species": "A", "block_id": "block1", "gene_a": "A_a", "gene_b": "A_b", "ks": "0.5",
              "ks_status": "ok", "species_event_id": event, "placement_status": "interval_supported"}
    args = classification_fixture(tmp_path, "((A_a:1,A_iso:1):1,A_b:2);", positions, [anchor],
                                  [{"species_event_id": event, "event_support": "WGD-supported"}])
    wgd_ssd.classify(args)
    rows = read_table(args.output / "duplication_origins.tsv")
    assert len(rows) == 2
    assert all(row["classification"] == "unresolved" and row["reason"] == "same_locus_annotation_ambiguity" for row in rows)


@pytest.mark.parametrize("different_families", [False, True])
def test_known_annotation_isoforms_cannot_inflate_family_gene_counts(tmp_path, different_families):
    members, counts = tmp_path / "members.tsv", tmp_path / "counts.tsv"
    rows = ([{"family_id": "OG1", "A": "A_t1, A_t2", "B": "B_t1"}] if not different_families else
            [{"family_id": "OG1", "A": "A_t1", "B": "B_t1"}, {"family_id": "OG2", "A": "A_t2", "B": ""}])
    write_rows(members, ["family_id", "A", "B"], rows)
    write_rows(counts, ["family_id", "A", "B"], [{"family_id": row["family_id"], "A": len(row["A"].split(",")),
               "B": 1 if row["B"] else 0} for row in rows])
    gff = tmp_path / "A.gff3"
    gff.write_text("chr1\tt\tmRNA\t1\t80\t.\t+\t.\tID=t1;Parent=locus1\n"
                   "chr1\tt\tmRNA\t1\t100\t.\t+\t.\tID=A_t2;Parent=locus1\n")
    plan = {"species": ["A", "B"], "counts": str(counts), "members": str(members),
            "genomes": [{"species": "A", "gff": str(gff), "feature": "mRNA", "attribute": "ID"}]}
    wgd_ssd.validated_members(plan)
    with pytest.raises(ValueError, match="isoform|locus"):
        wgd_ssd.validate_annotation_counts(plan, rows)


def mul_join_fixture(gene_text="(x1_X,x2_X);"):
    from nwkit.mul_reconcile import run_search, topology_text
    from nwkit.mul_reconcile_model import hypotheses
    from nwkit.mul_reconcile_nodes import write_node_diagnostics
    from nwkit.reconcile import build_reconciliation_table
    from nwkit.species_parser import get_species_parser
    from nwkit.util import read_tree

    gene = read_tree(gene_text, "auto", True, quiet=True)
    species = read_tree("((A,X),B);", "auto", True, quiet=True)
    parser = get_species_parser(species_parser="legacy", species_regex=r".*_([^_]+)$")
    candidates = hypotheses(species, "X", "A B X")
    scores = {c: s for c, s, _ in run_search(candidates, [gene], parser)}
    best = [c for c in candidates if scores[c.id] == min(scores.values())]
    handle = StringIO()
    write_node_diagnostics(handle, best, [gene], parser, species, max_state_pairs=10000000, max_maps=100000)
    rows = list(csv.DictReader(StringIO(handle.getvalue()), delimiter="\t"))
    model = {"method": "exact-MUL-LCA-DL-parsimony-v1", "num_gene_trees": 1,
             "node_diagnostics": {"schema": "nwkit-mul-node-assignments-v1"},
             "best_hypotheses": [c.id for c in best],
             "scores": [{"mul.tree": c.id, "score": scores[c.id], "h1.node": c.h1, "h2.node": c.h2,
                         "hypothesis.kind": c.kind, "labeled.tree": topology_text(c.tree)} for c in candidates]}
    reconciliation = build_reconciliation_table(gene, species, {n.name: parser.parse(n.name).species_label for n in gene.leaves()}).to_dict("records")
    origins = [{"gene_clade_id": r["gene_clade_id"], "classification": "SSD-supported", "reason": "fixture_evidence"}
               for r in reconciliation if r["event_type"] == "duplication"]
    return gene, species, rows, model, reconciliation, origins


def test_mul_node_join_keeps_ties_and_does_not_override_origin_evidence():
    from workflow.support.wgd_mul_diagnostics import join_nodes

    fixture = mul_join_fixture()
    origins = copy.deepcopy(fixture[-1])
    result = join_nodes(*fixture, "OGtest")
    assert fixture[-1] == origins
    assert len(result) == 1 and result[0]["classification"] == "SSD-supported"
    assert result[0]["mul_best_hypotheses"] == 3
    assert result[0]["mul_optimal_mappings"] == 6
    assert result[0]["mul_mapping_status"] == "ambiguous"
    assert result[0]["mul_dl_duplication"] == "none"
    assert "probability" not in result[0]


def test_mul_node_join_is_clade_based_not_postorder_based():
    from workflow.support.wgd_mul_diagnostics import join_nodes

    fixture = mul_join_fixture("((a_A,x1_X),(b_B,x2_X));")
    for node in fixture[0].traverse():
        node.children.reverse()
        node.dist = 77
    result = join_nodes(*fixture, "OGtest")
    assert len(result) == 3
    assert {r["mul_mapping_status"] for r in result} == {"consistent"}


@pytest.mark.parametrize("kind", ["gene_topology", "species_topology", "missing_node", "duplicate_node",
                                  "missing_candidate", "count", "duplication", "mapped_tips", "wrong_model",
                                  "candidate_metadata", "wrong_clade", "reconciliation", "origins",
                                  "leaf_species", "missing_column", "nonbest_candidate", "repeat_mapping"])
def test_mul_node_join_rejects_mismatched_or_incomplete_records(kind):
    from workflow.support.wgd_mul_diagnostics import join_nodes

    fixture = list(mul_join_fixture())
    rows, model = fixture[2], fixture[3]
    if kind == "gene_topology":
        rows[0]["gene_topology_id"] = "other tree"
    elif kind == "species_topology":
        rows[0]["species_topology_id"] = "other tree"
    elif kind == "missing_node":
        rows.pop()
    elif kind == "duplicate_node":
        rows.append(dict(rows[0]))
    elif kind == "missing_candidate":
        fixture[2] = [r for r in rows if r["mul.tree"] != "3"]
    elif kind == "count":
        rows[0]["optimal.mappings"] = "99"
    elif kind == "duplication":
        rows[0]["duplication"] = "2"
    elif kind == "mapped_tips":
        rows[0]["mul_descendant_tips"] = '["unknown"]'
    elif kind == "wrong_model":
        model["method"] = "locus-mc"
    elif kind == "candidate_metadata":
        rows[0]["h2.node"] = "wrong"
    elif kind == "wrong_clade":
        rows[0]["gene_clade_id"] = "unknown"
    elif kind == "reconciliation":
        fixture[4].pop()
    elif kind == "origins":
        fixture[5].clear()
    elif kind == "leaf_species":
        next(r for r in fixture[4] if r["event_type"] == "leaf")["species_name"] = "A"
    elif kind == "missing_column":
        rows[0].pop("gene_node")
    elif kind == "nonbest_candidate":
        model["best_hypotheses"].append(0)
    elif kind == "repeat_mapping":
        copies = [dict(r) for r in rows if r["mul.tree"] == "1" and r["mapping.id"] == "1"]
        for r in rows:
            if r["mul.tree"] == "1":
                r["optimal.mappings"] = "3"
        for r in copies:
            r["mapping.id"], r["optimal.mappings"] = "3", "3"
        rows.extend(copies)
    with pytest.raises(ValueError):
        join_nodes(*fixture, "OGtest")
    assert all("mul_mapping_status" not in n.props for n in fixture[0].traverse())


@pytest.mark.parametrize("supported", [False, True])
def test_real_classification_mul_diagnostics_preserve_wgd_unresolved_and_nhx(tmp_path, supported):
    from nwkit.clade_index import CladeIndex

    reference = tmp_path / "reference.nwk"
    reference.write_text("(((A:1,B:1)AB:1,C:2)ABC:1,D:3)ROOT;")
    species = wgd_ssd.species_tree(reference)
    branch = CladeIndex(species).clade_id_for_node(next(n for n in species.traverse() if n.name == "AB"))
    positions = [positioned("A_a"), positioned("A_c", chromosome="chr2"),
                 positioned("B_b"), positioned("B_d", chromosome="chr2")]
    anchor = {"species": "A", "block_id": "block1", "gene_a": "A_a", "gene_b": "A_c", "ks": "0.5",
              "ks_status": "ok", "species_event_id": branch, "placement_status": "interval_supported"}
    args = classification_fixture(tmp_path, "((A_a:1,B_b:1):1,(A_c:1,B_d:1):1)[&&NHX:D=Y:annotation=retained];",
                                  positions, [anchor] if supported else [],
                                  [{"species_event_id": branch, "event_support": "WGD-supported"}] if supported else [])
    wgd_ssd.classify(args)
    previous = (args.output / "duplication_origins.tsv").read_bytes()
    args.output = tmp_path / "with_mul"
    args.mul_diagnostics, args.mul_h1, args.mul_h2 = 1, "A,B", "A,B"
    args.mul_max_candidates, args.mul_max_state_pairs, args.mul_max_maps = 100, 100000, 10000
    wgd_ssd.classify(args)
    assert (args.output / "duplication_origins.tsv").read_bytes() == previous
    assert (args.output / "duplication_origins_mul.pdf").stat().st_size > 1000
    diagnostic = read_table(args.output / "node_diagnostics.tsv")
    assert len(diagnostic) == 3
    classes = {r["classification"] for r in diagnostic}
    assert ("WGD-supported" if supported else "unresolved") in classes
    nhx = (args.output / "classified_gene_tree.nhx").read_text()
    assert "annotation=retained" in nhx and "D=Y" in nhx and "mul_mapping_status=" in nhx
    assert json.loads((args.output / "summary.json").read_text())["mul_diagnostics"]["origin_rule_changed"] is False


@pytest.mark.parametrize("tips, label_width", [(80, 0), (8, 200)])
def test_mul_node_plot_deep_tree_labels_do_not_overlap_or_clip(tmp_path, monkeypatch, tips, label_width):
    import matplotlib.figure
    from nwkit.util import read_tree

    from workflow.support.wgd_mul_diagnostics import plot_diagnostics

    label = "W" * label_width
    text = f"({label}tip0,{label}tip1)"
    for number in range(2, tips):
        text = f"({text},{label}tip{number})"
    tree = read_tree(text + ";", "auto", True, quiet=True)
    for number, node in enumerate(tree.traverse("postorder")):
        if not node.is_leaf:
            node.add_prop("mul_mapping_status", "ambiguous")
            node.add_prop("mul_gene_node", number)
    save = matplotlib.figure.Figure.savefig

    def checked_save(fig, *args, **kwargs):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        boxes = [t.get_window_extent(renderer) for t in fig.axes[0].texts]
        for index, box in enumerate(boxes):
            assert fig.bbox.contains(box.x0, box.y0) and fig.bbox.contains(box.x1, box.y1)
            assert not any(box.overlaps(other) for other in boxes[index + 1:])
        return save(fig, *args, **kwargs)

    monkeypatch.setattr(matplotlib.figure.Figure, "savefig", checked_save)
    plot_diagnostics(tree, tmp_path / "deep.pdf")


@pytest.mark.parametrize("mul_enabled", [0, 1])
def test_mul_node_diagnostics_preserve_single_tip_rejection(tmp_path, mul_enabled):
    args = classification_fixture(tmp_path, "A_a:1[&&NHX:nwkit_rooted=yes];", [positioned("A_a")])
    args.mul_diagnostics, args.mul_h1, args.mul_h2 = mul_enabled, "A", "A"
    args.mul_max_candidates, args.mul_max_state_pairs, args.mul_max_maps = 100, 100000, 10000
    with pytest.raises(RuntimeError, match="at least two tips"):
        wgd_ssd.classify(args)
    assert not (args.output / "summary.json").exists()
    assert not (args.output / "node_diagnostics.tsv").exists()
