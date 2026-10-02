#!/usr/bin/env python3
"""Verify the native reconciliation exports required by GeneGalleon."""

import csv
import json
import subprocess
import tempfile
from pathlib import Path


def check_mul_reconciliation(work):
    gene, species = work / "mul_gene.nwk", work / "mul_species.nwk"
    scores, detail, model = work / "mul_scores.tsv", work / "mul_detail.tsv", work / "mul_model.json"
    gene.write_text("((a_A,x1_X),(b_B,x2_X));\n")
    species.write_text("((A,X),B);\n")
    subprocess.run([
        "nwkit", "mul-reconcile", "--infile", str(gene), "--species-tree", str(species),
        "--species-regex", ".*_([^_]+)$", "--h1", "X", "--h2", "B",
        "--outfile", str(scores), "--report", str(detail), "--model-out", str(model),
    ], check=True, stdout=subprocess.DEVNULL)
    with scores.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if int(rows[0]["score"]) != 0 or rows[0]["hypothesis.kind"] != "allopolyploid":
        raise ValueError("NWKIT MUL reconciliation failed the known allopolyploid control.")
    with detail.open() as handle:
        maps = list(csv.DictReader(handle, delimiter="\t"))
    if not maps or any(int(row["total.score"]) != 0 for row in maps):
        raise ValueError("NWKIT MUL reconciliation did not export its optimal mappings.")
    for row in maps:
        nodes = json.loads(row["node.maps"])
        if len(nodes) != 7 or sum(
            node["duplication"] + node["child_edge_losses"] + node["root_losses"]
            for node in nodes
        ) != int(row["total.score"]):
            raise ValueError("NWKIT MUL node mappings do not reconstruct D+L scores.")
    if json.loads(model.read_text())["method"] != "exact-MUL-LCA-DL-parsimony-v1":
        raise ValueError("NWKIT MUL reconciliation lacks the required model contract.")


def main():
    from nwkit.util import read_tree_strings

    with tempfile.TemporaryDirectory(prefix="gg-nwkit-reconciliation-") as temporary:
        work = Path(temporary)
        gene, species = work / "gene.nwk", work / "species.nwk"
        selected, candidates, table = work / "selected.nwk", work / "roots.nwk", work / "events.tsv"
        gene.write_text("((A_a_1:1,A_a_2:2):3,(A_a_3:4,A_a_4:5):6);")
        species.write_text("(A_a:1,B_b:1);")
        subprocess.run([
            "nwkit", "root", "--method", "reconciliation", "--infile", str(gene),
            "--species-tree", str(species), "--outfile", str(selected),
            "--candidates-out", str(candidates), "--duplication-cost", "1.5", "--loss-cost", "1",
        ], check=True, stdout=subprocess.DEVNULL)
        if len(read_tree_strings(str(candidates))) != 5:
            raise ValueError("NWKIT did not export all five tied root edges.")
        gene.write_text("((A_a_1:1,B_b_1:1):1,A_a_2:2);")
        species.write_text("(A_a:1,(B_b:1,C_c:1):1);")
        subprocess.run([
            "nwkit", "reconcile", "--infile", str(gene), "--species-tree", str(species),
            "--event-source", "lca", "--unmatched", "error", "--outfile", str(table),
        ], check=True, stdout=subprocess.DEVNULL)
        with table.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        if sum(row["event_type"] == "duplication" for row in rows) != 1:
            raise ValueError("NWKIT returned an incorrect duplication count.")
        if sum(int(row["implied_losses"]) for row in rows) != 2:
            raise ValueError("NWKIT returned an incorrect implied-loss count.")
        check_mul_reconciliation(work)
    print("NWKIT optimal roots, LCA losses and exact MUL reconciliation verified.")


if __name__ == "__main__":
    main()
