import argparse
import csv
import itertools
import json
import random
from types import SimpleNamespace

import numpy as np
import pytest

from workflow.support import pairwise_synteny as synteny
from workflow.support import pairwise_synteny_dotplot as dotplot
from workflow.support import pairwise_synteny_karyotype as karyotype
from workflow.support import pairwise_synteny_layout as layout


def write_inputs(tmp_path, fasta=">Species_name_t1\nMPEPTIDE\n>Species_name_t2\nMPEPTIDEAAA\n"):
    source = tmp_path / "Species_name.protein.fa"
    source.write_text(fasta)
    gff = tmp_path / "Species_name.gff3"
    gff.write_text("##gff-version 3\nchr1\ttest\tgene\t1\t100\t.\t+\t.\tID=g1\n"
                   "chr1\ttest\tmRNA\t1\t40\t.\t+\t.\tID=t1;Parent=g1\n"
                   "chr1\ttest\tmRNA\t1\t100\t.\t+\t.\tID=t2;Parent=g1\n")
    return {"species": "Species_name", "fasta": str(source), "gff": str(gff), "mode": "protein",
            "feature": "", "attribute": "", "genetic_code": None}


def test_prepare_selects_one_isoform_and_preserves_id_mapping(tmp_path):
    source = write_inputs(tmp_path)
    genes, metadata = synteny.prepare_genome(source, tmp_path, "target", 1)
    assert [(g.gene_id, g.start, g.end) for g in genes] == [("Species_name_t2", 0, 100)]
    assert metadata["collapsed_isoform_count"] == 1
    assert (tmp_path / "target.pep").read_text() == ">Species_name_t2\nMPEPTIDEAAA\n"
    with (tmp_path / "target.id_map.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["status"] for row in rows] == ["isoform_excluded", "selected"]
    assert all(row["locus_id"] == "g1" for row in rows)


@pytest.mark.parametrize("fasta,message", [
    (">t1\nMPEPTIDE\n>t1\nMPEPTIDE\n", "Duplicate"),
    (">Species_name_t1\nMPEPTIDE\n>t1\nMPEPTIDE\n", "Ambiguous"),
    (">unmapped\nMPEPTIDE\n", "GFF"),
    (">t1\nMP*EPTIDE\n", "internal stop"),
    (">t1\nMPſEPTIDE\n", "Invalid protein"),
])
def test_invalid_fasta_annotation_mapping_fails(tmp_path, fasta, message):
    source = write_inputs(tmp_path, fasta)
    with pytest.raises(ValueError, match=message):
        synteny.prepare_genome(source, tmp_path, "target", 1)


def test_cds_translation_uses_selected_genetic_code(tmp_path):
    source = write_inputs(tmp_path, ">t1\nATGTAGTGA\n")
    source.update(mode="cds", genetic_code=6)
    genes, _ = synteny.prepare_genome(source, tmp_path, "query", 1)
    assert genes[0].gene_id == "Species_name_t1"
    assert (tmp_path / "query.pep").read_text() == ">Species_name_t1\nMQ\n"


def test_unknown_display_chromosome_is_rejected(tmp_path):
    source = write_inputs(tmp_path)
    genes, _ = synteny.prepare_genome(source, tmp_path, "target", 1)
    with pytest.raises(ValueError, match="annotated chromosomes"):
        synteny.selected_seqids("Chr1", genes)


def test_source_discovery_rejects_multiple_annotation_releases(tmp_path):
    directory = tmp_path / "input/species_gff"
    directory.mkdir(parents=True)
    for release in (1, 2):
        (directory / f"Species_name.release{release}.gff3").write_text("##gff-version 3\n")
    with pytest.raises(ValueError, match="exactly one"):
        synteny.source_file(tmp_path, "", "species_gff", "Species_name", (".gff3",))


def recorded_connection_analysis(tmp_path):
    (tmp_path / "logs").mkdir()
    synteny.write_json(tmp_path / "summary.json", {"parameters": {
        "cscore": 0.85, "min_anchors": 6, "distance": 11, "quota": None}})
    (tmp_path / "logs/01.mcscan.log").write_text(
        "diamond blastp --evalue 2e-9 --outfmt 6\n"
        "running the local dups filter (tandem_Nmax=7)\n"
        "0 new pairs found (dist=5).\n")
    return tmp_path


def test_connection_legend_uses_recorded_nondefault_analysis_values(tmp_path):
    criteria = karyotype.connection_criteria(recorded_connection_analysis(tmp_path))
    assert criteria["protein_evalue"] == 2e-9
    assert criteria["seed_cscore"] == 0.85
    assert criteria["seed_min_unique_genes_per_genome"] == 6
    assert criteria["chaining_max_gene_rank_gap_per_genome"] == 11
    assert criteria["tandem_gene_rank_distance"] == 7
    assert criteria["liftover_gene_rank_distance"] == 5
    assert criteria["liftover_metric"] == "Manhattan" and criteria["liftover_bound"] == "strict"
    assert criteria["quota"] is None and criteria["ds_filter"] is False
    assert criteria["source_hashes"] == {str(path.resolve()): synteny.digest(path)
                                         for path in (tmp_path / "summary.json", tmp_path / "logs/01.mcscan.log")}
    text = karyotype.connection_legend_text(criteria)
    assert "E <= 2e-9" in text and "C-score >= 0.85" in text
    assert ">= 6 unique genes / genome" in text and "gap <= 11 gene ranks / genome" in text
    assert "tandem distance = 7 gene ranks" in text and "|dx| + |dy| < 5 gene ranks" in text
    assert "no quota; no dS filter" in text
    criteria["quota"] = "2:2"
    assert "quota = 2:2" in karyotype.connection_legend_text(criteria)


@pytest.mark.parametrize("defect", ["evalue", "tandem", "liftover", "ambiguous", "invalid"])
def test_connection_legend_refuses_unrecorded_or_invalid_thresholds(tmp_path, defect):
    recorded_connection_analysis(tmp_path)
    log_file = tmp_path / "logs/01.mcscan.log"
    if defect == "ambiguous":
        log_file.write_text(log_file.read_text() + "diamond blastp --evalue 1e-5\n")
    elif defect == "invalid":
        summary = json.loads((tmp_path / "summary.json").read_text())
        summary["parameters"]["cscore"] = 2
        synteny.write_json(tmp_path / "summary.json", summary)
    else:
        patterns = {"evalue": "--evalue", "tandem": "tandem_Nmax", "liftover": "new pairs found"}
        log_file.write_text("\n".join(line for line in log_file.read_text().splitlines() if patterns[defect] not in line))
    with pytest.raises(ValueError, match="Cannot annotate connections"):
        karyotype.connection_criteria(tmp_path)


def test_plan_rejects_duplicate_pair_ids_before_running_tools(tmp_path):
    pairs = tmp_path / "pairs.tsv"
    pairs.write_text("analysis_id\ttarget_species\tquery_species\na\tTarget_species\tQuery_species\na\tTarget_species\tQuery_species\n")
    args = argparse.Namespace(workspace=tmp_path, pairs=pairs, sequence_mode="auto", genetic_code=1,
                              cscore=0.7, min_anchors=4, distance=20, minimum_mapping_fraction=1,
                              formats="pdf", karyotype_sort="none")
    with pytest.raises(ValueError, match="unique"):
        synteny.build_plan(args)


def test_plot_settings_do_not_change_analysis_contract(tmp_path):
    plan = {"workspace": str(tmp_path), "parameters": {}, "tools": {"jcvi": "example"}, "formats": ["pdf"],
            "ds_tools": {"source_hashes": {}},
            "pairs": [{"analysis_id": "pair", "target_species": "Target_species", "query_species": "Query_species",
                       "target": {"fasta": "/input/t.fa", "gff": "/input/t.gff"},
                       "query": {"fasta": "/input/q.fa", "gff": "/input/q.gff"},
                       "target_seqids": "", "query_seqids": "",
                       "ds": {"target": {"fasta": "/input/t.cds.fa", "genetic_code": 1},
                              "query": {"fasta": "/input/q.cds.fa", "genetic_code": 1}}}]}
    before_analysis = synteny.contract_args(plan, "analysis")
    before_ds = synteny.contract_args(plan, "ds")
    before_plots = synteny.contract_args(plan, "plots")
    plan["formats"] = ["png"]
    plan["pairs"][0]["query_seqids"] = "chr2,chr1"
    assert synteny.contract_args(plan, "analysis") == before_analysis
    assert synteny.contract_args(plan, "plots") != before_plots
    before_plots = synteny.contract_args(plan, "plots")
    plan["karyotype_color"] = "homoeolog"
    assert synteny.contract_args(plan, "analysis") == before_analysis
    assert synteny.contract_args(plan, "plots") != before_plots
    before_plots = synteny.contract_args(plan, "plots")
    plan["karyotype_scale"] = "independent"
    assert synteny.contract_args(plan, "analysis") == before_analysis
    assert synteny.contract_args(plan, "ds") == before_ds
    assert synteny.contract_args(plan, "plots") != before_plots
    before_plots = synteny.contract_args(plan, "plots")
    plan["karyotype_track_order"] = "query-target"
    assert synteny.contract_args(plan, "analysis") == before_analysis
    assert synteny.contract_args(plan, "ds") == before_ds
    assert synteny.contract_args(plan, "plots") != before_plots
    before_plots = synteny.contract_args(plan, "plots")
    plan.update(dotplot_sort="none", dotplot_min_length=2000000)
    plan["pairs"][0]["dotplot_lengths"] = {"target": {"kind": "genome", "path": "/input/genome.fa"}}
    assert synteny.contract_args(plan, "analysis") == before_analysis
    assert synteny.contract_args(plan, "plots") != before_plots
    for mode in ("target", "query", "target_length", "query_length", "both_length"):
        before_analysis = synteny.contract_args(plan, "analysis")
        before_plots = synteny.contract_args(plan, "plots")
        plan["karyotype_sort"] = mode
        assert synteny.contract_args(plan, "analysis") == before_analysis
        assert synteny.contract_args(plan, "plots") != before_plots


@pytest.mark.parametrize("scale_mode", ["shared", "independent"])
@pytest.mark.parametrize("unequal", [False, True])
def test_both_tracks_solver_matches_independent_joint_permutations(tmp_path, scale_mode, unequal):
    selected = [["m1", "m2", "m3"], ["q1", "q2", "q3"]]
    totals = ((6, 8, 10), (12, 18, 4)) if unequal else ((10, 10, 10), (10, 10, 10))
    genomes = [display_genes(prefix, tuple(zip(seqids, counts, strict=True)))
               for prefix, seqids, counts in zip(("m", "q"), selected, totals, strict=True)]
    partners = [(selected[0][i], selected[1][j], totals[0][i] - 1, totals[1][j] - 1)
                for i, j in ((0, 2), (1, 0), (2, 1))]
    blocks = [(f"m{a}_0", f"m{a}_{count_a}", f"q{b}_0", f"q{b}_{count_b}")
              for a, b, count_a, count_b in partners]
    simple = tmp_path / "blocks.simple"
    simple.write_text("\n".join(" ".join(block) + " 10 +" for block in blocks) + "\n")
    ordered, metadata = synteny.order_karyotype(selected, genomes, tmp_path / "unused", "both_length", simple, scale_mode)
    geometries = layout.pair_track_geometry(selected, genomes, scale_mode)
    def independent(orders):
        positions = [layout.offsets(order, geometry[0], geometry[3]) for order, geometry in zip(orders, geometries, strict=True)]
        total = 0
        for a, b, count_a, count_b in partners:
            width = (geometries[0][2] * count_a + geometries[1][2] * count_b) / 2
            center_a = positions[0][a] + geometries[0][2] * count_a / 2
            center_b = positions[1][b] + geometries[1][2] * count_b / 2
            total += width * ((20 * (center_a - center_b)) ** 2 + 3.2 ** 2) ** 0.5
        return total
    optimum = min(independent((a, b)) for a, b in itertools.product(itertools.permutations(selected[0]), itertools.permutations(selected[1])))
    assert independent(ordered) == pytest.approx(optimum, abs=1e-12)
    assert metadata["objective_after"] == pytest.approx(independent(ordered))
    assert metadata["globally_optimal"] is True
    assert metadata["fixed_side"] is None
    if not unequal:
        assert ordered[1] != selected[1]
    assert metadata["scale_mode"] == scale_mode
    assert metadata["track_ratios"] == pytest.approx([g[2] for g in geometries])


@pytest.mark.parametrize("scale_mode", ["shared", "independent"])
def test_pair_geometry_uses_one_width_per_gene_only_in_shared_mode(scale_mode):
    selected = [["t1", "t2"], [f"q{i}" for i in range(18)]]
    genomes = [display_genes("t", (("t1", 8), ("t2", 8))),
               display_genes("q", tuple((sid, 4) for sid in selected[1]))]
    geometry = layout.pair_track_geometry(selected, genomes, scale_mode)
    spans = [sum(widths.values()) + (len(widths) - 1) * gap for widths, _, _, gap in geometry]
    assert max(spans) == pytest.approx(layout.XEND - layout.XSTART)
    if scale_mode == "shared":
        assert geometry[0][2] == geometry[1][2]
        assert geometry[0][0]["t1"] == pytest.approx(2 * geometry[1][0]["q0"])
        assert spans[0] < spans[1]
    else:
        assert spans == pytest.approx([0.8, 0.8])
        assert geometry[0][2] != geometry[1][2]
    with pytest.raises(ValueError, match="karyotype-scale"):
        layout.pair_track_geometry(selected, genomes, "bad")


def test_large_joint_search_is_bounded_deterministic_and_not_claimed_global(tmp_path, monkeypatch):
    monkeypatch.setattr(layout, "JOINT_EXACT_WORK_LIMIT", 0)
    monkeypatch.setattr(layout, "JOINT_RESTARTS", 3)
    monkeypatch.setattr(layout, "JOINT_MAX_ROUNDS", 4)
    selected = [["m1", "m2"], ["q1", "q2"]]
    genomes = [display_genes("m", (("m1", 10), ("m2", 10))), display_genes("q", (("q1", 10), ("q2", 10)))]
    simple = tmp_path / "blocks.simple"
    simple.write_text("mm1_0 mm1_9 qq2_0 qq2_9 10 +\nmm2_0 mm2_9 qq1_0 qq1_9 10 +\n")
    first = layout.order_by_ribbon_length(selected, genomes, simple, None)
    assert first == layout.order_by_ribbon_length(selected, genomes, simple, None)
    metadata = first[1]
    assert metadata["globally_optimal"] is False
    assert metadata["objective_after"] <= metadata["objective_before"]
    assert metadata["track_solves"] <= 1 + 3 * 4 * 2


def display_genes(prefix, chromosomes):
    return [SimpleNamespace(seqid=seqid, gene_id=f"{prefix}{seqid}_{rank}", start=rank * 100, end=rank * 100 + 50)
            for seqid, count in chromosomes for rank in range(count)]


@pytest.mark.parametrize("mode", ["target", "query"])
def test_karyotype_sort_fixes_opposite_track_and_uses_dominant_partner(tmp_path, mode):
    moving = display_genes("m", (("m1", 4), ("m2", 1), ("m3", 1), ("m4", 1)))
    fixed = display_genes("f", (("q1", 4), ("q2", 2), ("hidden", 1)))
    # m1's dominant partner is q2, even though it also has secondary support on q1.
    pairs = [("mm1_0", "fq2_0"), ("mm1_1", "fq2_1"), ("mm1_2", "fq1_0"),
             ("mm1_3", "fhidden_0"), ("mm2_0", "fq1_3"), ("mm3_0", "fq1_0")]
    if mode == "target":
        genomes = (moving, fixed)
        selected = (["m1", "m2", "m3", "m4"], ["q1", "q2"])
    else:
        genomes = (fixed, moving)
        selected = (["q1", "q2"], ["m1", "m2", "m3", "m4"])
        pairs = [(b, a) for a, b in pairs]
    anchors = tmp_path / "lifted.anchors"
    anchors.write_text("###\n" + "\n".join(f"{a}\t{b}\t100" for a, b in pairs) + "\n")
    ordered, metadata = synteny.order_karyotype(selected, genomes, anchors, mode)
    moving_index = 0 if mode == "target" else 1
    assert ordered[moving_index] == ["m3", "m2", "m1", "m4"]
    assert ordered[1 - moving_index] == ["q1", "q2"]
    assert selected[moving_index] == ["m1", "m2", "m3", "m4"]  # caller input is not mutated
    assert metadata["orientation_changed"] is False
    m1 = next(row for row in metadata["chromosomes"] if row["seqid"] == "m1")
    assert m1["dominant_partner"] == "q2"
    assert m1["dominant_anchor_count"] == 2
    assert m1["anchor_count"] == 3  # hidden chromosomes do not influence display sorting
    assert metadata["display_order"][mode] == ordered[moving_index]


def test_sort_deduplicates_anchors_and_preserves_exact_ties_and_unsupported_order(tmp_path):
    moving = display_genes("m", (("m1", 2), ("m2", 1), ("u1", 1), ("u2", 1)))
    fixed = display_genes("f", (("q1", 3), ("q2", 1)))
    anchors = tmp_path / "lifted.anchors"
    anchors.write_text("###\nmm1_0 fq1_0 100\nmm1_1 fq1_2 100\nmm2_0 fq1_1 100\n"
                       "###\nmm1_0 fq1_0 100\nmm1_0 fq1_0 100\n")
    selected = (["u2", "m2", "m1", "u1"], ["q1", "q2"])
    ordered, metadata = synteny.order_karyotype(selected, (moving, fixed), anchors, "target")
    assert ordered == [["m2", "m1", "u2", "u1"], ["q1", "q2"]]
    m1 = next(row for row in metadata["chromosomes"] if row["seqid"] == "m1")
    assert m1["anchor_count"] == 2
    assert m1["partner_gene_rank_mean"] == 1


def test_dominant_partner_tie_uses_explicit_fixed_order_not_natural_order(tmp_path):
    moving = display_genes("m", (("m1", 2), ("m2", 1)))
    fixed = display_genes("f", (("q1", 1), ("q2", 1)))
    anchors = tmp_path / "lifted.anchors"
    anchors.write_text("mm1_0 fq1_0 100\nmm1_1 fq2_0 100\nmm2_0 fq1_0 100\n")
    selected = (["m2", "m1"], ["q2", "q1"])
    ordered, metadata = synteny.order_karyotype(selected, (moving, fixed), anchors, "target")
    assert ordered == [["m1", "m2"], ["q2", "q1"]]
    assert metadata["chromosomes"][1]["dominant_partner"] == "q2"


def test_disabled_sort_preserves_supplied_order_without_reading_anchors(tmp_path):
    selected = (["chr10", "chr2"], ["scaffold2", "scaffold1"])
    ordered, metadata = synteny.order_karyotype(selected, ((), ()), tmp_path / "unused", "none")
    assert ordered == list(selected)
    assert metadata["method"] == "input_order"


@pytest.mark.parametrize("line", ["missing q1 100", "m1 missing 100", "m1"])
def test_sort_rejects_invalid_anchors(tmp_path, line):
    anchors = tmp_path / "anchors"
    anchors.write_text(line + "\n")
    moving = [SimpleNamespace(seqid="m", gene_id="m1", start=0, end=10)]
    fixed = [SimpleNamespace(seqid="q", gene_id="q1", start=0, end=10)]
    with pytest.raises(ValueError, match="Invalid JCVI anchor"):
        synteny.order_karyotype((["m"], ["q"]), (moving, fixed), anchors, "target")


def test_sort_rejects_unknown_mode(tmp_path):
    with pytest.raises(ValueError, match="karyotype-sort"):
        synteny.order_karyotype(([], []), ((), ()), tmp_path / "unused", "both")


@pytest.mark.parametrize("seed", range(8))
def test_weighted_length_solver_matches_independent_all_permutations(seed):
    rng = random.Random(seed)
    order = [str(i) for i in range(5)]
    widths = {sid: rng.uniform(0.025, 0.12) for sid in order}
    targets = {sid: [(rng.uniform(0.12, 0.92), rng.uniform(0.001, 1)) for _ in range(3)] for sid in order}

    def cost(sid, position):
        return sum(weight * np.hypot((position + widths[sid] / 2 - fixed) * 20, 3.2)
                   for fixed, weight in targets[sid])

    def independent(permutation):
        value, position = 0, 0.12
        for sid in permutation:
            for fixed, weight in targets[sid]:
                value += weight * ((20 * (position + widths[sid] / 2 - fixed)) ** 2 + 3.2 ** 2) ** 0.5
            position += widths[sid] + 0.01
        return value

    result, metadata = layout.minimize_order(order, widths, 0.01, cost)
    assert independent(result) == pytest.approx(min(map(independent, itertools.permutations(order))), abs=1e-12)
    assert metadata["objective_after"] == pytest.approx(independent(result))
    assert metadata["objective_after"] <= metadata["objective_before"] + 1e-12
    assert metadata["globally_optimal"] is True


def test_length_ties_preserve_input_and_large_track_reports_non_exact_solver():
    order = ["b", "c", "a", "d"]
    result, metadata = layout.minimize_order(order, dict.fromkeys(order, 0.1), 0.01, lambda _sid, start: np.zeros_like(start))
    assert result == order
    assert metadata["objective_after"] == 0
    _, metadata = layout.minimize_order(order, dict.fromkeys(order, 0.1), 0.01,
                                       lambda _sid, start: np.zeros_like(start), exact_limit=3)
    assert metadata["globally_optimal"] is False
    assert metadata["solver"] == "bounded_adjacent_swap"


@pytest.mark.parametrize("mode", ["target_length", "query_length"])
def test_ribbon_width_weights_and_complete_block_positions_choose_optimal_track(tmp_path, mode):
    moving = display_genes("m", (("m1", 10), ("m2", 10)))
    fixed = display_genes("f", (("q1", 10), ("q2", 10)))
    # m1 connects broadly to q2 but also narrowly to q1. Counting the two
    # connections equally would lose the width information.
    blocks = [("mm1_0", "mm1_9", "fq2_0", "fq2_9"),
              ("mm1_0", "mm1_1", "fq1_0", "fq1_1"),
              ("mm2_0", "mm2_9", "fq1_0", "fq1_9")]
    if mode == "target_length":
        genomes, selected = (moving, fixed), (["m1", "m2"], ["q1", "q2"])
    else:
        genomes, selected = (fixed, moving), (["q1", "q2"], ["m1", "m2"])
        blocks = [(c, d, a, b) for a, b, c, d in blocks]
    simple = tmp_path / "blocks.simple"
    simple.write_text("\n".join(" ".join((*row, "10", "+")) for row in blocks) + "\n")
    ordered, metadata = synteny.order_karyotype(selected, genomes, tmp_path / "unused", mode, simple)
    moving_index = 0 if mode == "target_length" else 1
    assert ordered[moving_index] == ["m2", "m1"]
    assert ordered[1 - moving_index] == ["q1", "q2"]
    assert metadata["displayed_block_count"] == 3
    assert metadata["objective_after"] < metadata["objective_before"]
    assert metadata["weight"] == "mean_endpoint_width"
    assert metadata["globally_optimal"] is True
    assert metadata["orientation_changed"] is False


def test_layout_gap_and_ranks_follow_jcvi_gene_rank_geometry():
    selected = [f"chr{i}" for i in range(18)]
    genes = display_genes("g", [(sid, 10 + i) for i, sid in enumerate(selected)])
    widths, ranks, ratio, gap = layout.track_geometry(selected, genes)
    assert gap == pytest.approx(min(0.01, 0.01 * 16 / 18 + 0.001))
    assert sum(widths.values()) + 17 * gap == pytest.approx(0.8)
    assert ranks["gchr17_26"] == 26
    assert widths["chr17"] == pytest.approx(27 * ratio)


def test_layout_shared_start_ranks_use_literal_identifier_not_end_or_natural_id():
    genes = [SimpleNamespace(seqid="chr1", gene_id="g2", start=0, end=10),
             SimpleNamespace(seqid="chr1", gene_id="g10", start=0, end=20)]
    assert layout.track_geometry(["chr1"], genes)[1] == {"g10": 0, "g2": 1}


def test_empty_anchors_are_not_a_completed_analysis(tmp_path):
    (tmp_path / "target.query.lifted.anchors").write_text("###\n")
    with pytest.raises(ValueError, match="No syntenic blocks"):
        synteny.summarize_anchors(tmp_path, ((), ()))


def test_changed_source_is_rejected_before_publication(tmp_path):
    source = tmp_path / "source.fa"
    source.write_text(">gene\nMPEPTIDE\n")
    plan = {"input_hashes": {str(source): synteny.digest(source)}}
    synteny.verify_inputs(plan)
    source.write_text(">gene\nMPEPTIDEX\n")
    with pytest.raises(ValueError, match="Input changed"):
        synteny.verify_inputs(plan)


def test_dotplot_physical_length_boundary_and_explicit_karyotype_order(tmp_path):
    genomes = [display_genes("t", (("short", 1), ("exact", 2), ("long", 2))),
               display_genes("q", (("query", 2),))]
    target, query = tmp_path / "target.sizes", tmp_path / "query.gff"
    target.write_text("seqid\tlength\nshort\t999999\nexact\t1000000\nlong\t2000000\n")
    query.write_text("##sequence-region query 1 1000000\n")
    pair = {"target_species": "Target_species", "query_species": "Query_species",
            "target": {"gff": str(query)}, "query": {"gff": str(query)},
            "dotplot_lengths": {"target": {"kind": "sizes", "path": str(target)},
                                "query": {"kind": "gff", "path": str(query)}}}
    anchors = "###\ntshort_0 qquery_0 10\n###\ntexact_0 qquery_0 10\ntlong_1 qquery_1 10\n###\ntshort_0 qquery_1 10\n"
    (tmp_path / "target.query.lifted.anchors").write_text(anchors)
    metadata = dotplot.prepare_dotplot(tmp_path, pair, genomes, [["long", "exact", "short"], ["query"]], 1000000, "karyotype")
    assert metadata["display_order"] == {"target": ["long", "exact"], "query": ["query"]}
    assert metadata["anchor_count"] == 2
    assert metadata["excluded_anchor_count"] == 2
    assert metadata["displayed_block_count"] == 1
    assert (tmp_path / "dotplot.anchors").read_text() == "###\ntexact_0 qquery_0 10\ntlong_1 qquery_1 10\n"
    assert (tmp_path / "target.query.lifted.anchors").read_text() == anchors
    assert (tmp_path / "dotplot.target.bed").read_text().splitlines()[0].startswith("long\t")
    assert "short" not in (tmp_path / "dotplot.anchors").read_text()
    with (tmp_path / "dotplot_display.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["length_bp"] == "999999"
    assert rows[0]["selected_for_dotplot"] == "0"
    assert rows[1]["display_rank"] == "2"


@pytest.mark.parametrize("kind,text,message", [
    ("gff", "chr1\t.\tgene\t1\t50\t.\t+\t.\tID=g\n", "Missing physical"),
    ("gff", "##sequence-region chr1 10 2000000\n", "Missing physical"),
    ("gff", "##sequence-region chr1 1 2000000\n##sequence-region chr1 1 2000001\n", "Duplicate"),
    ("sizes", "chr1\t0\n", "invalid"),
    ("sizes", "chr1\t10\n", "exceeds"),
    ("genome", ">chr1\nACGT\n>chr1\nACGT\n", "Duplicate"),
])
def test_dotplot_length_metadata_must_be_trustworthy(tmp_path, kind, text, message):
    path = tmp_path / "lengths"
    path.write_text(text)
    with pytest.raises(ValueError, match=message):
        dotplot.chromosome_lengths({"kind": kind, "path": str(path)}, display_genes("t", (("chr1", 1),)), 1000000)


def test_dotplot_genome_lengths_use_full_sequence_not_gene_extent(tmp_path):
    path = tmp_path / "genome.fa"
    path.write_text(">chr1 comment\n" + "N" * 999999 + "\nA\n")
    lengths, _ = dotplot.chromosome_lengths({"kind": "genome", "path": str(path)}, display_genes("t", (("chr1", 1),)), 1000000)
    assert lengths == {"chr1": 1000000}


def test_dotplot_2x2_sort_groups_both_tracks_preserves_weak_connections_and_deduplicates(tmp_path):
    selected = [["t1", "t2", "t3", "t4"], ["q1", "q2", "q3", "q4"]]
    genomes = [display_genes("T", tuple((x, 8) for x in selected[0])),
               display_genes("Q", tuple((x, 8) for x in selected[1]))]
    pairs = [(a, b) for a in ("t1", "t3") for b in ("q2", "q4")]
    pairs += [(a, b) for a in ("t2", "t4") for b in ("q1", "q3")]
    lines = [f"T{a}_{i} Q{b}_{i} 10" for a, b in pairs for i in range(8)]
    lines += ["Tt1_0 Qq1_0 10"] * 50  # Duplicate anchors cannot increase grouping support.
    path = tmp_path / "anchors"
    original = "###\n" + "\n".join(lines) + "\n"
    path.write_text(original)
    ordered, metadata = dotplot.homoeolog_order(selected, genomes, path)
    assert ordered == [["t1", "t3", "t2", "t4"], ["q2", "q4", "q1", "q3"]]
    assert all(group["anchor_counts_2x2"] == [8, 8, 8, 8] for group in metadata["groups"])
    assert (ordered, metadata) == dotplot.homoeolog_order(selected, genomes, path)
    assert metadata["globally_optimal"] is False
    assert path.read_text() == original


def test_dotplot_2x2_sort_does_not_invent_groups_from_two_diagonal_matches(tmp_path):
    selected = [["t1", "t2"], ["q1", "q2"]]
    genomes = [display_genes("T", (("t1", 1), ("t2", 1))), display_genes("Q", (("q1", 1), ("q2", 1)))]
    path = tmp_path / "anchors"
    path.write_text("###\nTt1_0 Qq1_0 10\nTt2_0 Qq2_0 10\n")
    ordered, metadata = dotplot.homoeolog_order(selected, genomes, path)
    assert ordered == selected
    assert metadata["groups"] == []


def test_dotplot_sort_preserves_original_within_chromosome_bed_ranks(tmp_path):
    genes = [SimpleNamespace(seqid="chr", gene_id="g2", start=0, end=10),
             SimpleNamespace(seqid="chr", gene_id="g10", start=0, end=20)]
    gff = tmp_path / "gff"
    gff.write_text("##sequence-region chr 1 1000000\n")
    pair = {"target_species": "Target_species", "query_species": "Query_species",
            "target": {"gff": str(gff)}, "query": {"gff": str(gff)}}
    (tmp_path / "target.query.lifted.anchors").write_text("###\ng2 g10 1\n")
    dotplot.prepare_dotplot(tmp_path, pair, [genes, genes], [["chr"], ["chr"]], 1000000, "karyotype")
    assert [line.split()[3] for line in (tmp_path / "dotplot.target.bed").read_text().splitlines()] == ["g2", "g10"]


@pytest.mark.parametrize("text", ["species\twrong_column\nSpecies_name\t1\n", "species\tgenetic_code\nSpecies_name\n"])
def test_plan_rejects_malformed_genetic_code_metadata_cleanly(tmp_path, text):
    directory = tmp_path / "input/species_genetic_code"
    directory.mkdir(parents=True)
    (directory / "species_genetic_code.tsv").write_text(text)
    pairs = tmp_path / "pairs.tsv"
    pairs.write_text("analysis_id\ttarget_species\tquery_species\npair\tTarget_species\tQuery_species\n")
    args = argparse.Namespace(workspace=tmp_path, pairs=pairs, sequence_mode="auto", genetic_code=1,
                              cscore=0.7, min_anchors=4, distance=20, minimum_mapping_fraction=1,
                              formats="png", karyotype_sort="none")
    with pytest.raises(ValueError, match="[Gg]enetic-code"):
        synteny.build_plan(args)


def test_genetic_code_metadata_is_rechecked_before_publication(tmp_path, monkeypatch):
    for directory in ("species_cds", "species_gff", "species_genetic_code"):
        (tmp_path / "input" / directory).mkdir(parents=True)
    for species in ("Target_species", "Query_species"):
        (tmp_path / f"input/species_cds/{species}.fa").write_text(">gene\nATGGCT\n")
        (tmp_path / f"input/species_gff/{species}.gff").write_text("##gff-version 3\n")
    code_file = tmp_path / "input/species_genetic_code/species_genetic_code.tsv"
    code_file.write_text("species\tgenetic_code\nTarget_species\t1\nQuery_species\t1\n")
    pairs = tmp_path / "pairs.tsv"
    pairs.write_text("analysis_id\ttarget_species\tquery_species\npair\tTarget_species\tQuery_species\n")
    args = argparse.Namespace(workspace=tmp_path, pairs=pairs, sequence_mode="cds", genetic_code=1,
                              cscore=0.7, min_anchors=4, distance=20, minimum_mapping_fraction=1,
                              formats="png", karyotype_sort="none")
    monkeypatch.setattr(synteny, "tool_identity", lambda: {"jcvi": "test"})
    plan = synteny.build_plan(args)
    synteny.verify_inputs(plan)
    code_file.write_text("species\tgenetic_code\nTarget_species\t2\nQuery_species\t1\n")
    with pytest.raises(ValueError, match="Input changed"):
        synteny.verify_inputs(plan)


def test_pair_table_is_rechecked_before_publication_without_invalidating_display_only_changes(tmp_path, monkeypatch):
    for directory in ("species_cds", "species_gff"):
        (tmp_path / "input" / directory).mkdir(parents=True)
    for species in ("Target_species", "Query_species"):
        (tmp_path / f"input/species_cds/{species}.fa").write_text(">gene\nATGGCT\n")
        (tmp_path / f"input/species_gff/{species}.gff").write_text("##gff-version 3\n")
    pairs = tmp_path / "pairs.tsv"
    pairs.write_text("analysis_id\ttarget_species\tquery_species\npair\tTarget_species\tQuery_species\n")
    args = argparse.Namespace(workspace=tmp_path, pairs=pairs, sequence_mode="cds", genetic_code=1,
                              cscore=0.7, min_anchors=4, distance=20, minimum_mapping_fraction=1,
                              formats="png", karyotype_sort="none")
    monkeypatch.setattr(synteny, "tool_identity", lambda: {"jcvi": "test"})
    plan = synteny.build_plan(args)
    before = synteny.contract_args(plan, "analysis")
    synteny.verify_inputs(plan)
    pairs.write_text("analysis_id\ttarget_species\tquery_species\ttarget_seqids\n"
                     "pair\tTarget_species\tQuery_species\tchr1\n")
    with pytest.raises(ValueError, match="Input changed"):
        synteny.verify_inputs(plan)
    # The next complete run may still reuse its scientific analysis: a pair-table
    # byte fingerprint belongs to the run guard, not every phase's cache key.
    updated = synteny.build_plan(args)
    assert updated["karyotype_scale"] == "shared"
    assert synteny.contract_args(updated, "analysis") == before


def test_shared_homoeolog_colors_are_opt_in_and_cover_unlisted_chromosomes(tmp_path):
    selected = [["t1", "t2"], ["q1", "q2"]]
    genomes = [display_genes("T", (("t1", 1), ("t2", 1), ("hidden", 1))),
               display_genes("Q", (("q1", 1), ("q2", 1)))]
    anchors = tmp_path / "anchors"
    anchors.write_text("###\n" + "".join(f"T{a}_0 Q{b}_0 10\n" for a in selected[0] for b in selected[1]))
    default = karyotype.chromosome_colors(selected, genomes, anchors)
    assert default["mode"] == "chromosome"
    assert default["groups"] == []
    assert default["chromosomes"]["target"]["t1"] != default["chromosomes"]["target"]["t2"]
    assert "hidden" in default["chromosomes"]["target"]
    paired = karyotype.chromosome_colors(selected, genomes, anchors, "homoeolog")
    assert len(paired["groups"]) == 1
    colors = {paired["chromosomes"][side][seqid] for side, track in zip(("target", "query"), selected, strict=True) for seqid in track}
    assert len(colors) == 1
    assert paired["chromosomes"]["target"]["hidden"] not in colors
    assert default == karyotype.chromosome_colors([track[::-1] for track in selected], genomes, anchors)
    with pytest.raises(ValueError, match="karyotype-color"):
        karyotype.chromosome_colors(selected, genomes, anchors, "unknown")


@pytest.mark.parametrize("totals,expected", [([3, 4], 1), ([60, 500], 10), ([600, 400], 50), ([18000, 23000], 2000)])
def test_gene_scale_is_integral_and_uses_the_smaller_track(totals, expected):
    assert karyotype.gene_scale(totals) == expected


def test_dS_blue_palette_darkens_monotonically_and_stays_visible_on_white():
    from matplotlib.colors import to_rgb

    from workflow.support.pairwise_synteny_style import MISSING_COLOR, ds_colormap

    cmap = ds_colormap()
    rgb = np.asarray([cmap(x)[:3] for x in np.linspace(0, 1, 256)])
    linear = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    luminance = linear @ np.asarray((0.2126, 0.7152, 0.0722))
    assert np.all(np.diff(luminance) < 0)
    assert np.all(1.05 / (luminance + 0.05) >= 3)
    assert np.all(rgb[:, 2] > rgb[:, 1]) and np.all(rgb[:, 1] > rgb[:, 0])
    assert tuple(cmap(np.nan)[:3]) == to_rgb(MISSING_COLOR)
