"""Same-coordinate source-DNA repairs retain ownership and independent support."""
import copy
import json
import sqlite3
import sys
import tempfile
import unittest
from importlib import import_module
from pathlib import Path
from unittest.mock import patch
from xml.etree import ElementTree

import pytest

SUPPORT = Path(__file__).resolve().parents[1] / "support"
sys.path.insert(0, str(SUPPORT))
refine = import_module("gene_model_refinement")
review = import_module("plot_gene_model_refinement")
SPECIES = "Species_target"
GENE = SPECIES + "_g1"
DONORS = ("Species_relative", "Species_balanced")
PARAMS = {**refine.DEFAULTS, "isoform_adoption": "conservation_supported"}
GENOMIC_SEQUENCE = "ATG" + "AAA" * 194 + "TGA"
SUPPLIED_SEQUENCE = GENOMIC_SEQUENCE[:-3] + "NNN" + "AAA" * 18 + "NNN"


def source_fixture():
    candidate = {
        "candidate_id": SPECIES + "_t1", "source_transcript_id": "t1",
        "source_gene_id": "g1", "gene_token": "g1", "seqid": "chr1",
        "strand": "+", "blocks": [[1000, 1588, 0]], "source_blocks": [[1000, 1588, 0]],
        "cds": GENOMIC_SEQUENCE, "origin": "original", "junctions": [],
        "source_fasta_ids": ["g1"],
        "source_cds": [{"fasta_id": "g1", "cds": SUPPLIED_SEQUENCE, "sequence_agreement": False}],
    }
    candidate["quality"] = refine.validate_candidate(candidate, 1)
    candidate["quality"].update(sequence_mismatch=True, source_sequence_agreement=False,
                                valid_orf=False, usable=False)
    gene = {
        "species": SPECIES, "gene_id": GENE, "source_gene_id": "g1",
        "gene_token": "g1", "seqid": "chr1", "strand": "+",
        "source_baseline_candidate_id": candidate["candidate_id"],
    }
    return {"gene": gene, "candidate": candidate,
            "supplied_sequence": SUPPLIED_SEQUENCE, "genomic_sequence": GENOMIC_SEQUENCE}


FIXTURE = source_fixture()


def data():
    gene = {**copy.deepcopy(FIXTURE["gene"]), "candidates": [copy.deepcopy(FIXTURE["candidate"])],
            "source_baseline_candidate_id": FIXTURE["candidate"]["candidate_id"]}
    cat = {"species": SPECIES, "genetic_code": 1, "loci": [gene]}
    edges = [{"species_a": SPECIES, "gene_a": GENE, "species_b": donor,
              "gene_b": donor + "_synthetic", "ambiguous": False} for donor in DONORS]
    models = [{"gene_id": GENE, "seqid": gene["seqid"], "strand": gene["strand"],
               "cds": copy.deepcopy(gene["candidates"][0]["blocks"]),
               "sequence": FIXTURE["genomic_sequence"], "problems": [],
               "coverage": 1.0, "identity": 1.0, "donor_species": donor,
               "donor_candidate": donor + "_synthetic",
               "synthetic_evidence": True,
               "support_group": "nearest" if donor == DONORS[0] else "phylogenetic_balanced"}
              for donor in DONORS]
    return cat, edges, models

def classify(cat=None, edges=None, models=None):
    original = data()
    cat, edges, models = (cat or original[0], edges or original[1], models or original[2])
    with patch.object(refine.rescue, "run", side_effect=AssertionError("External dispatch forbidden")):
        return refine.classify_predictions(models, cat, edges, PARAMS, [], "0" * 64)

def decision():
    return {"selections": [{"species": SPECIES, "gene_id": GENE,
                            "candidate_id": FIXTURE["candidate"]["candidate_id"],
                            "source_transcript_id": FIXTURE["candidate"]["source_transcript_id"],
                            "status": "ambiguous_correspondence", "reason": "unresolved_paralog_copy"}]}

class SameCoordinateRepair(unittest.TestCase):
    def test_two_distinct_external_species_repair_without_mutating_original(self):
        cat, edges, models = data()
        frozen = copy.deepcopy(cat)
        row, = classify(cat, edges, models)
        self.assertEqual(row["status"], "accepted")
        self.assertEqual(row["change_type"], "model_revision")
        self.assertEqual(row["same_coordinate_source_repair"]["independent_external_species"], sorted(DONORS))
        self.assertEqual(row["candidate"]["source_transcript_id"], FIXTURE["candidate"]["source_transcript_id"])
        self.assertEqual(row["candidate"]["gene_id"], GENE)
        self.assertFalse(row["candidate"]["quality"]["rna_supported"])
        self.assertEqual(cat, frozen)
        self.assertEqual(len(cat["loci"][0]["candidates"][0]["source_cds"][0]["cds"]), 645)

    def test_selection_metrics_and_flags_follow_final_representative_with_bounded_graph_audit(self):
        row, = classify()
        old_id = row["same_coordinate_source_repair"]["original_candidate_id"]
        repair_id = row["candidate"]["candidate_id"]
        before = decision()
        before["metrics"] = {"changed_representatives": 0, "pair_evaluations": 7}
        before["scores"] = [
            {"species": SPECIES, "gene_id": GENE, "candidate_id": old_id, "selected": True, "score": None},
            {"species": SPECIES, "gene_id": GENE, "candidate_id": repair_id, "selected": False, "score": 3},
            {"species": "Other_species", "gene_id": "other", "candidate_id": "other.1", "selected": True, "score": 1}]
        original_decision = copy.deepcopy(before["selections"][0])
        original_metrics = copy.deepcopy(before["metrics"])
        scores_identity = before["scores"]
        unrelated_score = before["scores"][2]
        after = refine.apply_same_coordinate_repairs(before, [row])
        self.assertIs(after["scores"], scores_identity)
        self.assertIs(after["scores"][2], unrelated_score)
        self.assertEqual([r["selected"] for r in after["scores"]], [False, True, True])
        self.assertEqual(after["metrics"]["changed_representatives"], 1)
        self.assertEqual(after["metrics"]["same_coordinate_sequence_repairs"], 1)
        self.assertEqual(after["source_correspondence_metrics"], original_metrics)
        self.assertEqual(after["selections"][0]["source_correspondence_decision"], original_decision)
        audit, = after["same_coordinate_repair_audit"]
        self.assertEqual(audit["source_selected_flags"],
                         [{"candidate_id": old_id, "selected": True},
                          {"candidate_id": repair_id, "selected": False}])
        self.assertEqual(audit["change_count_delta"], 1)

    def test_zero_repair_selection_is_byte_identical_and_same_scores_object(self):
        row, = classify()
        row["status"] = "proposal"
        selection = decision()
        selection["metrics"] = {"changed_representatives": 0, "pair_evaluations": 7}
        selection["scores"] = [{"species": SPECIES, "gene_id": GENE,
                                "candidate_id": FIXTURE["candidate"]["candidate_id"], "selected": True}]
        before = json.dumps(selection, indent=2).encode()
        score_identity = selection["scores"]
        after = refine.apply_same_coordinate_repairs(selection, [row])
        self.assertIs(after, selection)
        self.assertIs(after["scores"], score_identity)
        self.assertEqual(json.dumps(after, indent=2).encode(), before)

    def test_already_graph_selected_repair_does_not_double_count_changed_representative(self):
        row, = classify()
        before = decision()
        before["selections"][0].update(candidate_id=row["candidate"]["candidate_id"], status="conserved")
        before["metrics"] = {"changed_representatives": 1}
        after = refine.apply_same_coordinate_repairs(before, [row])
        self.assertEqual(after["metrics"]["changed_representatives"], 1)
        self.assertEqual(after["same_coordinate_repair_audit"][0]["change_count_delta"], 0)

    def test_invalid_isoform_repair_is_revision_and_does_not_force_normal_representative(self):
        cat, edges, models = data()
        normal = copy.deepcopy(cat["loci"][0]["candidates"][0])
        normal["candidate_id"] += "_normal_isoform"
        normal["source_transcript_id"] += "_normal_isoform"
        normal["blocks"] = [[a+3, b+3, phase] for a, b, phase in normal["blocks"]]
        normal["quality"].update(sequence_mismatch=False, valid_orf=True, usable=True,
                                 source_sequence_agreement=True, representative_eligible=True)
        normal["source_cds"][0].update(cds=normal["cds"], sequence_agreement=True)
        cat["loci"][0]["candidates"].append(normal)
        row, = classify(cat, edges, models)
        self.assertEqual(row["status"], "accepted")
        self.assertEqual(row["change_type"], "model_revision")
        selected = decision()
        selected["selections"][0].update(candidate_id=normal["candidate_id"],
                                         source_transcript_id=normal["source_transcript_id"],
                                         status="conserved", reason="normal_supported_isoform")
        self.assertEqual(refine.apply_same_coordinate_repairs(selected, [row]), selected)

    def test_select_streamed_predictions_are_full_byte_receipt_fenced(self):
        for mutation in ("none", "before_guard", "during_build"):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory(prefix="same_coordinate_selection_toy_") as td:
                root = Path(td)
                (root / "plan.json").write_text("{}")
                paths = {
                    "catalog": root / "catalog" / SPECIES,
                    "index": root / "catalog_index_final",
                    "correspondence": root / "correspondence",
                    "predictions": root / "predictions" / SPECIES,
                }
                row, = classify()
                for kind, directory in paths.items():
                    directory.mkdir(parents=True)
                    member = ("predictions.json" if kind == "predictions" else
                              "edges.json" if kind == "correspondence" else "small.json")
                    payload = json.dumps([row] if kind == "predictions" else [])
                    (directory / member).write_text(payload)
                    key = {"fixture": kind}
                    (directory / "receipt.json").write_text(json.dumps({
                        "key": key, "files": {member: refine.digest(directory / member)}}))
                predictions = paths["predictions"] / "predictions.json"
                value = {"species": [SPECIES], "request": {"parameters": PARAMS}}
                if mutation == "before_guard":
                    predictions.write_text("[]")
                original_stream = refine.stream_json_array

                def stream(path, mutation=mutation, original_stream=original_stream):
                    yield from original_stream(path)
                    if mutation == "during_build":
                        path.write_text("[]")

                with patch.object(refine, "catalog_index", return_value=paths["index"] / "loci.sqlite3"), \
                        patch.object(refine, "correspondence", return_value=paths["correspondence"]), \
                        patch.object(refine, "select_from_store", return_value=decision()), \
                        patch.object(refine, "load", return_value=value), \
                        patch.object(refine, "stream_json_array", side_effect=stream), \
                        patch.object(refine.rescue, "run", side_effect=AssertionError("No scientific dispatch")):
                    if mutation == "none":
                        with refine.invocation_context():
                            result = refine.select(root, value, predictions=True)
                        selected = json.loads((result / "selection.json").read_text())["selections"][0]
                        self.assertEqual(selected["status"], "sequence_repaired")
                        self.assertTrue(refine.rescue.verified(result, json.loads((result / "receipt.json").read_text())["key"]))
                    else:
                        with self.assertRaises(ValueError), refine.invocation_context():
                            refine.select(root, value, predictions=True)
                        self.assertFalse((root / "selection_final" / "receipt.json").exists())

    def test_single_species_two_support_groups_never_count_as_two(self):
        cat, edges, models = data()
        models[1]["donor_species"] = models[0]["donor_species"]
        models[1]["donor_candidate"] = models[0]["donor_candidate"]
        row, = classify(cat, edges, models)
        self.assertEqual(row["status"], "proposal")
        self.assertIn("existing_coding_path", row["problems"])
        self.assertIn("unresolved_source_sequence_mismatch", row["problems"])
        self.assertNotIn("same_coordinate_source_repair", row)

    def test_self_does_not_count_as_independent_external_species(self):
        cat, edges, models = data()
        models[1]["donor_species"] = SPECIES
        edges.append({"species_a": SPECIES, "gene_a": GENE, "species_b": SPECIES, "gene_b": GENE, "ambiguous": False})
        row, = classify(cat, edges, models)
        self.assertEqual(row["status"], "proposal")

    def test_untrusted_or_low_quality_donor_cannot_establish_repair(self):
        for variant in ("untrusted", "coverage", "identity", "sequence"):
            cat, edges, models = data()
            if variant == "untrusted":
                edges[1]["ambiguous"] = True
            elif variant == "coverage":
                models[1]["coverage"] = 0.2
            elif variant == "identity":
                models[1]["identity"] = 0.2
            else:
                models[1]["sequence"] = models[1]["sequence"][:30] + "AAG" + models[1]["sequence"][33:]
            with self.subTest(variant=variant):
                row, = classify(cat, edges, models)
                self.assertEqual(row["status"], "proposal")
                self.assertNotIn("same_coordinate_source_repair", row)

    def test_grouped_support_never_substitutes_different_DNA_at_same_coordinates(self):
        cat, edges, models = data()
        different = copy.deepcopy(models[0])
        different["sequence"] = different["sequence"][:30] + "AAG" + different["sequence"][33:]
        row, = classify(cat, edges, [different, *models])
        self.assertEqual(row["status"], "proposal")
        self.assertNotIn("same_coordinate_source_repair", row)
        self.assertIn("existing_coding_path", row["problems"])

    def test_normal_raw_model_keeps_duplicate_path_gate_and_selection(self):
        cat, edges, models = data()
        c = cat["loci"][0]["candidates"][0]
        c["quality"].update(sequence_mismatch=False, valid_orf=True, usable=True, source_sequence_agreement=True)
        c["source_cds"][0].update(cds=c["cds"], sequence_agreement=True)
        row, = classify(cat, edges, models)
        self.assertEqual(row["status"], "proposal")
        self.assertIn("existing_coding_path", row["problems"])
        before = decision()
        self.assertEqual(refine.apply_same_coordinate_repairs(before, [row]), before)

    def test_other_coordinates_remain_held_even_in_repair_only_search(self):
        cat, edges, models = data()
        for m in models:
            m["cds"] = [[a+3, b+3, p] for a, b, p in m["cds"]]
            m["same_coordinate_repair_only"] = True
        row, = classify(cat, edges, models)
        self.assertEqual(row["status"], "proposal")
        self.assertIn("same_coordinate_repair_required", row["problems"])

    def test_original_partial_exception_unknown_phase_or_ambiguous_ownership_never_repairs(self):
        for flag in ("annotated_partial", "annotated_exception", "phase_unknown", "translation_exception",
                     "structure_problem", "ambiguous_coordinates", "frameshift"):
            cat, edges, models = data()
            if flag == "ambiguous_coordinates":
                cat["loci"][0][flag] = True
            else:
                cat["loci"][0]["candidates"][0]["quality"][flag] = "annotated_translation_exception" if "exception" in flag else True
            with self.subTest(flag=flag):
                row, = classify(cat, edges, models)
                self.assertEqual(row["status"], "proposal")
                self.assertNotIn("same_coordinate_source_repair", row)

    def test_noncanonical_splice_or_other_locus_collision_never_repairs(self):
        cat, edges, models = data()
        for m in models:
            m["problems"] = ["noncanonical_splice"]
        row, = classify(cat, edges, models)
        self.assertEqual(row["status"], "proposal")
        self.assertIn("noncanonical_splice", row["problems"])
        cat, edges, models = data()
        other = copy.deepcopy(cat["loci"][0])
        other["gene_id"] += "_other"
        cat["loci"].append(other)
        row, = classify(cat, edges, models)
        self.assertEqual(row["status"], "proposal")
        self.assertIn("overlap_other_locus", row["problems"])

    def test_partial_orf_cannot_be_promoted(self):
        cat, edges, models = data()
        for m in models:
            m["sequence"] = "TTT" + m["sequence"][3:]
        cat["loci"][0]["candidates"][0]["cds"] = models[0]["sequence"]
        row, = classify(cat, edges, models)
        self.assertEqual(row["status"], "proposal")
        self.assertIn("invalid_predicted_coding_path", row["problems"])

    def test_unknown_and_absent_genetic_code_fail_closed(self):
        cat, edges, models = data()
        cat["genetic_code"] = 999
        with self.assertRaisesRegex(ValueError, "Unsupported genetic code"):
            classify(cat, edges, models)
        cat, edges, models = data()
        del cat["genetic_code"]
        row, = classify(cat, edges, models)
        self.assertEqual(row["status"], "proposal")

    def test_copy_ambiguity_is_audited_separately_from_exact_locus_DNA_repair(self):
        cat, edges, models = data()
        edges.append({"species_a": SPECIES, "gene_a": GENE, "species_b": "Other_species",
                      "gene_b": "other_synthetic", "ambiguous": True})
        row, = classify(cat, edges, models)
        self.assertEqual(row["status"], "accepted")
        self.assertTrue(row["same_coordinate_source_repair"]["correspondence_ambiguity"])
        after = refine.apply_same_coordinate_repairs(decision(), [row])["selections"][0]
        self.assertEqual(after["status"], "sequence_repaired")
        self.assertEqual(after["source_correspondence_decision"], decision()["selections"][0])
        self.assertEqual(after["orthology"], "unassigned")
        self.assertEqual(after["expected_copy"], "unassigned")

    def test_mutated_candidate_context_not_reused_as_repair(self):
        row, = classify()
        for variant in ("coordinates", "DNA", "code", "species_support"):
            bad = copy.deepcopy(row)
            if variant == "coordinates":
                bad["candidate"]["blocks"][0][0] += 3
            elif variant == "DNA":
                bad["candidate"]["cds"] = "TTT" + bad["candidate"]["cds"][3:]
            elif variant == "code":
                bad["same_coordinate_source_repair"]["genetic_code"] = 999
            else:
                bad["same_coordinate_source_repair"]["independent_external_species"] = [DONORS[0], DONORS[0]]
            with self.subTest(variant=variant):
                before = decision()
                self.assertEqual(refine.apply_same_coordinate_repairs(before, [bad]), before)

    def test_unqualified_same_coordinate_mismatch_retains_original_two_gates(self):
        cat, edges, models = data()
        row, = classify(cat, edges, models[:1])
        self.assertEqual(row["status"], "proposal")
        self.assertIn("existing_coding_path", row["problems"])
        self.assertIn("unresolved_source_sequence_mismatch", row["problems"])
        self.assertEqual(row["donors"], [DONORS[0]])
        self.assertNotIn("same_coordinate_source_repair", row)

    def test_masked_terminal_stop_alone_is_existing_convention_but_extension_is_not(self):
        import gene_model_catalog
        c = copy.deepcopy(FIXTURE["candidate"])
        self.assertEqual(gene_model_catalog._source_convention(c["cds"][:-3] + "NNN", c), "masked_terminal_stop")
        self.assertEqual(gene_model_catalog._source_convention(FIXTURE["supplied_sequence"], c), "")

    def test_finalize_preserves_gene_transcript_original_cds_and_GFF(self):
        row, = classify()
        chosen = refine.apply_same_coordinate_repairs(decision(), [row])["selections"]
        with tempfile.TemporaryDirectory(prefix="same_coordinate_finalize_toy_") as td:
            root = Path(td)
            final = root / "selection_final"
            final.mkdir()
            (final / "selection.json").write_text(json.dumps({"selections": chosen}))
            (final / "representative_map.tsv").write_text("species\tgene_id\tcandidate_id\n")
            (final / "receipt.json").write_text('{"key":{},"files":{}}')
            db = root / "catalog_index_final/loci.sqlite3"
            db.parent.mkdir()
            (db.parent / "receipt.json").write_text('{"key":{},"files":{}}')
            catalog = root / "catalog" / SPECIES
            catalog.mkdir(parents=True)
            (catalog / "catalog_metadata.json").write_text(json.dumps({"species": SPECIES, "genetic_code": 1}))
            (catalog / "receipt.json").write_text('{"key":{},"files":{}}')
            gene = data()[0]["loci"][0]
            gene["candidates"].append(row["candidate"])
            raw = root / "raw.fa"
            raw.write_text(">" + GENE + "\n" + FIXTURE["supplied_sequence"] + "\n")
            gff = root / "original.gff3"
            annotation = "##gff-version 3\n" + refine.candidate_gff(gene, FIXTURE["candidate"], SPECIES,
                                             gene_id=gene["source_gene_id"], gene_token=GENE)
            gff.write_text(annotation)
            genome = root / "genome.fa"
            genome.write_text(">chr1\n" + "N" * 1000 + FIXTURE["genomic_sequence"] + "N" * 10 + "\n")
            value = {"species": [SPECIES], "request": {"sources": {SPECIES: {
                "fasta": str(raw), "gff": str(gff), "genome": str(genome), "genetic_code": 1}}}}
            tmp = root / "effective"
            tmp.mkdir()
            def build_only(unused_root, name, key, builder, **kwargs):
                self.assertEqual(name, "effective")
                builder(tmp)
                return tmp
            with patch.object(refine, "select", return_value=final), \
                    patch.object(refine, "catalog_index", return_value=db), \
                    patch.object(refine, "iter_loci", return_value=iter([gene])), \
                    patch.object(refine, "stage", side_effect=build_only), \
                    patch.object(refine.rescue, "run", side_effect=AssertionError("No scientific dispatch")):
                refine.finalize(root, value)
            self.assertEqual((tmp / "source_cds" / (SPECIES + ".fa")).read_bytes(), raw.read_bytes())
            self.assertEqual((tmp / "source_annotation" / (SPECIES + ".gff3")).read_bytes(), gff.read_bytes())
            self.assertEqual((tmp / "full_annotation" / (SPECIES + ".gff3")).read_text(), annotation)
            self.assertEqual((tmp / "species_gff" / (SPECIES + ".gff3")).read_text().count("\tmRNA\t"), 1)
            self.assertIn("ID=" + FIXTURE["candidate"]["source_transcript_id"], (tmp / "species_gff" / (SPECIES + ".gff3")).read_text())
            self.assertEqual((tmp / "species_cds" / (SPECIES + ".fa")).read_text(),
                             ">" + GENE + "\n" + FIXTURE["genomic_sequence"] + "\n")
            self.assertEqual((tmp / "all_candidates" / (SPECIES + ".fa")).read_text().count(">"), 2)
            self.assertEqual((tmp / "effective_exclusions.tsv").read_text().count("\n"), 1)
            summary = json.loads((tmp / "summary.json").read_text())
            self.assertEqual(summary["same_coordinate_sequence_repairs"], 1)
            self.assertEqual(summary["changed_representatives"], 1)



TARGET = "Species_target"


def plot_fixtures():
    return import_module("workflow.tests.test_plot_gene_model_refinement")


def replace_member(directory, name, value):
    path = directory / name
    path.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n")
    receipt_path = directory / "receipt.json"
    receipt = json.loads(receipt_path.read_text())
    receipt["files"][name] = review.digest(path)
    receipt_path.write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n")


def refresh_dependency(directory, kind, name, member_directory):
    path = directory / "receipt.json"
    receipt = json.loads(path.read_text())
    dependencies = receipt["key"].get("dependencies", {})
    if kind not in dependencies:
        return
    if name is None:
        dependencies[kind] = review.digest(member_directory / "receipt.json")
    else:
        dependencies[kind][name] = review.digest(member_directory / "receipt.json")
    path.write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n")


def future_plot_schema(root, *, selected_repair=True, include_accepted=True):
    """Adapt only an owned toy publication to exercise new report schema.

    This is a plot-schema fixture, not a biological admission.
    The source-DNA repair tests above verify the independent support gates.
    """
    changes = json.loads((root / "effective/changes.json").read_text())
    chosen = next(r for r in changes if r["species"] == TARGET)
    old_id = chosen["candidate_id"]
    db = root / "catalog_index_final/loci.sqlite3"
    with sqlite3.connect(db) as connection:
        locus = json.loads(connection.execute(
            "SELECT json FROM loci WHERE species=? AND gene_id=?",
            (TARGET, chosen["gene_id"])).fetchone()[0])
        original = next(c for c in locus["candidates"] if c["candidate_id"] == old_id)
        template = original if selected_repair else next(
            c for c in locus["candidates"] if c["candidate_id"] != old_id)
        candidate = copy.deepcopy(template)
        candidate["candidate_id"] += "_same_coordinate_fixture"
        candidate["origin"] = "predicted"
        candidate["support"] = {
            "donors": ["Species_donor1", "Species_donor2"], "rna_paths": [],
            "class": "homology_supported_same_coordinate_repair"}
        candidate["same_coordinate_source_repair"] = {"schema": 1, "fixture_only": True}
        locus["candidates"].append(candidate)
        connection.execute("UPDATE loci SET json=? WHERE species=? AND gene_id=?",
                           (json.dumps(locus, sort_keys=True), TARGET, chosen["gene_id"]))
        connection.execute(
            "INSERT INTO candidate_owners(species,candidate_id,gene_id) VALUES(?,?,?)",
            (TARGET, candidate["candidate_id"], chosen["gene_id"]))
    receipt_path = db.parent / "receipt.json"
    receipt = json.loads(receipt_path.read_text())
    receipt["files"][db.name] = review.digest(db)
    receipt_path.write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n")
    model = {
        "gene_id": chosen["gene_id"], "status": "accepted", "change_type": "model_revision",
        "evidence_class": "homology_supported_same_coordinate_repair",
        "donors": candidate["support"]["donors"], "rna_paths": [], "problems": [],
        "alignments": [{"donor_species": donor, "donor_candidate": donor + "_t",
                        "coverage": 1.0, "identity": 1.0, "supports_path": True, "problems": []}
                       for donor in candidate["support"]["donors"]],
        "candidate": candidate,
    }
    prediction_dir = root / "predictions" / TARGET
    replace_member(prediction_dir, "predictions.json", [model] if include_accepted else [])
    refresh_dependency(db.parent, "predictions", TARGET, prediction_dir)
    if selected_repair:
        previous = copy.deepcopy(chosen)
        chosen.update(candidate_id=candidate["candidate_id"], status="sequence_repaired",
                      selected_origin="predicted", reason="trusted_same_coordinate_genome_cds",
                      source_correspondence_decision=previous,
                      orthology="unassigned", expected_copy="unassigned")
    final = root / "selection_final"
    selection = json.loads((final / "selection.json").read_text())
    for row in selection["selections"]:
        if row["species"] == TARGET and selected_repair:
            row.update({k: v for k, v in chosen.items() if k != "selected_origin"})
    replace_member(final, "selection.json", selection)
    refresh_dependency(final, "predictions", TARGET, prediction_dir)
    refresh_dependency(final, "index", db.parent.name, db.parent)
    replace_member(root / "effective", "changes.json", changes)
    refresh_dependency(root / "effective", "selection", None, final)
    refresh_dependency(root / "effective", "index", db.parent.name, db.parent)
    return old_id, candidate["candidate_id"]


@pytest.mark.parametrize("case_name", [
    "test_review_counts_complete_publication_and_renders_empty_change_summary",
    "test_review_rejects_prediction_content_tampering",
    "test_locus_plot_escapes_labels_and_keeps_genomic_spacing_and_orientation",
])
def test_existing_plot_fixtures_use_repair_consumer(case_name, tmp_path, monkeypatch):
    fixtures = plot_fixtures()
    monkeypatch.setattr(fixtures, "review", review)
    case = getattr(fixtures, case_name)
    if "locus_plot" in case_name:
        case()
    else:
        case(tmp_path)


def test_repair_count_changed_representative_and_gallery_agree(tmp_path):
    root = plot_fixtures().completed(tmp_path)
    old_id, repair_id = future_plot_schema(root)
    data = review.collect(root, max_loci=2, preferred_species=TARGET)
    stats = data["species"][TARGET]
    assert stats["accepted_repair_paths"] == 1
    assert stats["accepted_isoform_paths"] == 0
    assert stats["changed_representatives"] == 1
    assert stats["selection_status"] == {"sequence_repaired": 1}
    detail, = data["details"]
    assert detail["selection"]["candidate_id"] == repair_id
    assert detail["selection"]["source_correspondence_decision"]["candidate_id"] == old_id
    assert detail["selection"]["orthology"] == "unassigned"
    assert data["changed_loci_available"] == 1
    out = tmp_path / "repaired_report"
    out.mkdir()
    review.plot_summary(data, out)
    review.write_review(data, out)
    ElementTree.parse(out / "summary.svg")
    html = (out / "review.html").read_text()
    assert "sequence_repaired" in html and repair_id in html
    assert "Repair coding paths" in html and "Changed representatives" in html


def test_sequence_repaired_status_enters_gallery_without_prediction_fallback(tmp_path):
    root = plot_fixtures().completed(tmp_path)
    _, repair_id = future_plot_schema(root, include_accepted=False)
    data = review.collect(root, max_loci=2)
    assert data["species"][TARGET]["accepted_repair_paths"] == 0
    assert data["species"][TARGET]["changed_representatives"] == 1
    assert data["changed_loci_available"] == 1
    assert data["details"][0]["selection"]["candidate_id"] == repair_id


def assert_same_render(baseline_data, candidate_data, output):
    assert json.dumps(baseline_data).encode() == json.dumps(candidate_data).encode()
    old = output / "baseline"
    new = output / "candidate"
    old.mkdir(parents=True)
    new.mkdir()
    review.plot_summary(baseline_data, old)
    review.write_review(baseline_data, old)
    review.plot_summary(candidate_data, new)
    review.write_review(candidate_data, new)
    assert (old / "summary.png").read_bytes() == (new / "summary.png").read_bytes()
    assert (old / "review.html").read_bytes() == (new / "review.html").read_bytes()


def test_zero_repair_preserves_data_and_PNG_HTML_bytes(tmp_path):
    root = plot_fixtures().completed(tmp_path)
    assert_same_render(review.collect(root), review.collect(root), tmp_path / "zero_report")


def test_normal_representative_stays_unchanged_with_accepted_other_path_repair(tmp_path):
    root = plot_fixtures().completed(tmp_path)
    old_id, repair_id = future_plot_schema(root, selected_repair=False)
    candidate = review.collect(root)
    stats = candidate["species"][TARGET]
    assert stats["accepted_repair_paths"] == 1
    assert stats["changed_representatives"] == 0
    detail, = candidate["details"]
    assert detail["selection"]["candidate_id"] == old_id
    assert detail["selection"]["candidate_id"] != repair_id
    assert_same_render(review.collect(root), candidate, tmp_path / "normal_report")


def test_same_coordinate_repair_import_selection_export_end_to_end(tmp_path, monkeypatch):
    """Real ownership import and export; only external prediction is synthetic."""
    input_root = tmp_path / "input"
    input_root.mkdir()
    rows = []
    annotation = ("##gff-version 3\n"
                  "chr1\tfixture\tgene\t1\t588\t.\t+\t.\tID=g1\n"
                  "chr1\tfixture\tmRNA\t1\t588\t.\t+\t.\tID=t1;Parent=g1\n"
                  "chr1\tfixture\tCDS\t1\t588\t.\t+\t0\tID=c1;Parent=t1\n")
    for name in (SPECIES, *DONORS):
        fasta, gff, genome = [input_root / (name + suffix)
                              for suffix in (".cds.fa", ".gff3", ".genome.fa")]
        fasta.write_text(">g1\n" + (SUPPLIED_SEQUENCE if name == SPECIES else GENOMIC_SEQUENCE) + "\n")
        genome.write_text(">chr1\n" + GENOMIC_SEQUENCE + "\n")
        gff.write_text(annotation)
        rows.append({"species": name, "cds": str(fasta), "gff": str(gff),
                     "genome": str(genome), "genetic_code": 1})
    inputs, links = input_root / "inputs.tsv", input_root / "edges.tsv"
    refine.rescue.write_tsv(inputs, list(rows[0]), [list(r.values()) for r in rows])
    edges = [{"species_a": SPECIES, "gene_a": GENE, "species_b": donor,
              "gene_b": donor + "_g1"} for donor in DONORS]
    refine.rescue.write_tsv(links, list(edges[0]), [list(r.values()) for r in edges])
    # Freeze a synthetic predictor identity: fast CI has no external miniprot.
    # Prediction/search dispatch below remains forbidden; genomic QC is real.
    original_which = refine.shutil.which
    monkeypatch.setattr(
        refine.shutil, "which",
        lambda name, *args, **kwargs: (
            None if name == "miniprot" else original_which(name, *args, **kwargs)))
    with pytest.raises(ValueError, match="miniprot is required for prediction"):
        refine.plan(tmp_path / "missing_predictor", inputs=inputs, edges=links,
                    mode="conservative", isoform_adoption="conservation_supported")
    predictor = tmp_path / "fixture_miniprot"
    predictor.write_text("#!/bin/sh\nexit 97\n")
    predictor.chmod(0o755)
    monkeypatch.setattr(
        refine.shutil, "which",
        lambda name, *args, **kwargs: (
            str(predictor) if name == "miniprot" else original_which(name, *args, **kwargs)))
    root = tmp_path / "refinement"
    value = refine.plan(root, inputs=inputs, edges=links, mode="conservative",
                        isoform_adoption="conservation_supported")
    assert value["request"]["miniprot"] == {
        "path": str(predictor), "sha256": refine.digest(predictor)}
    initial = refine.select(root, value, predictions=False)
    correspondence = root / "correspondence"
    graph = json.loads((correspondence / "edges.json").read_text())
    catalog_dir = root / "catalog" / SPECIES
    catalog = json.loads((catalog_dir / "catalog_metadata.json").read_text())
    catalog["loci"] = list(refine.iter_loci(root / "catalog_index/loci.sqlite3", SPECIES))
    original, = catalog["loci"][0]["candidates"]
    assert original["quality"]["sequence_mismatch"]
    assert original["cds"] == GENOMIC_SEQUENCE
    assert original["source_cds"][0]["cds"] == SUPPLIED_SEQUENCE
    raw_models = []
    with refine.indexed_genome(rows[0]["genome"]) as genome:
        for donor in DONORS:
            model = {"gene_id": GENE, "seqid": "chr1", "strand": "+",
                     "cds": [[0, 588, 0]], "frameshift": False,
                     "coverage": 1.0, "identity": 1.0,
                     "donor_species": donor, "donor_candidate": donor + "_g1",
                     "support_group": "nearest" if donor == DONORS[0] else "phylogenetic_balanced"}
            checked = refine.rescue.validate_model(model, genome, 1, value["request"]["parameters"])
            assert checked["problems"] == []
            raw_models.append(checked)
    accepted, = refine.classify_predictions(
        raw_models, catalog, graph, value["request"]["parameters"], [],
        refine.digest(Path(rows[0]["genome"])))
    assert accepted["status"] == "accepted"
    assert accepted["change_type"] == "model_revision"
    frozen_predictions = {}
    for name in value["species"]:
        dependencies = {"catalog": {n: refine.digest(root / "catalog" / n / "receipt.json")
                                    for n in value["species"]},
                        "initial": refine.digest(initial / "receipt.json"),
                        "correspondence": refine.digest(correspondence / "receipt.json")}
        models = [accepted] if name == SPECIES else []

        def builder(directory, models=models):
            refine.atomic_json(directory / "predictions.json", models)

        frozen_predictions[name] = refine.stage(
            root, "predictions/" + name,
            {"dependencies": dependencies, "synthetic_external_prediction_fixture": True}, builder)

    def synthetic_prediction(unused_root, unused_value, name, *args, **kwargs):
        assert unused_root == root
        path = frozen_predictions[name]
        receipt = json.loads((path / "receipt.json").read_text())
        assert refine.rescue.verified(path, receipt["key"])
        return path

    monkeypatch.setattr(refine, "predict_species", synthetic_prediction)
    monkeypatch.setattr(refine.rescue, "run",
                        lambda *args, **kwargs: pytest.fail("External prediction/search forbidden"))
    with refine.invocation_context():
        final_db = refine.catalog_index(root, value, predictions=True)
        with sqlite3.connect(final_db) as connection:
            imported = json.loads(connection.execute(
                "SELECT json FROM loci WHERE species=? AND gene_id=?", (SPECIES, GENE)).fetchone()[0])
            owner, = connection.execute(
                "SELECT gene_id FROM candidate_owners WHERE species=? AND candidate_id=?",
                (SPECIES, accepted["candidate"]["candidate_id"])).fetchone()
        assert owner == GENE
        assert len(imported["candidates"]) == 2
        assert imported["candidates"][0]["candidate_id"] == original["candidate_id"]
        assert imported["candidates"][0]["source_cds"][0]["cds"] == SUPPLIED_SEQUENCE
        selected = refine.select(root, value, predictions=True)
        selection = json.loads((selected / "selection.json").read_text())
        representative = next(r for r in selection["selections"] if r["species"] == SPECIES)
        assert representative["status"] == "sequence_repaired"
        assert representative["candidate_id"] == accepted["candidate"]["candidate_id"]
        assert representative["source_transcript_id"] == "t1"
        assert representative["gene_id"] == GENE
        effective = refine.finalize(root, value)
        refine.verify_inputs(effective / "inputs.tsv")
    assert (effective / "source_cds" / (SPECIES + ".fa")).read_bytes() == Path(rows[0]["cds"]).read_bytes()
    assert (effective / "source_annotation" / (SPECIES + ".gff3")).read_text() == annotation
    assert (effective / "full_annotation" / (SPECIES + ".gff3")).read_text() == annotation
    assert (effective / "species_cds" / (SPECIES + ".fa")).read_text() == ">" + GENE + "\n" + GENOMIC_SEQUENCE + "\n"
    assert (effective / "analysis_cds" / (SPECIES + ".fa")).read_text() == ">" + GENE + "\n" + GENOMIC_SEQUENCE + "\n"
    assert (effective / "species_protein" / (SPECIES + ".fa")).read_text() == ">" + GENE + "\n" + "M" + "K" * 194 + "\n"
    selected_gff = (effective / "species_gff" / (SPECIES + ".gff3")).read_text()
    assert selected_gff.count("\tmRNA\t") == 1
    features = [line.split("\t") for line in selected_gff.splitlines()
                if line and not line.startswith("#")]
    gene_ids = [refine.parse_gff_attributes(fields[8])["ID"]
                for fields in features if fields[2] == "gene"]
    transcript_ids = [refine.parse_gff_attributes(fields[8])["ID"]
                      for fields in features if fields[2] == "mRNA"]
    assert gene_ids == [("g1",)]
    assert transcript_ids == [("t1",)]
    assert "Parent=t1" in selected_gff
    assert (effective / "all_candidates" / (SPECIES + ".fa")).read_text().count(">") == 2
    changes = json.loads((effective / "changes.json").read_text())
    change = next(r for r in changes if r["species"] == SPECIES)
    assert change["status"] == "sequence_repaired"
    assert change["selected_origin"] == "predicted"
    assert change["gene_id"] == GENE
    summary = json.loads((effective / "summary.json").read_text())
    assert summary["same_coordinate_sequence_repairs"] == 1
    assert summary["changed_representatives"] == 1
