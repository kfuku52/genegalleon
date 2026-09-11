#!/usr/bin/env python3
"""Verify the native reconciliation exports required by GeneGalleon."""

import csv
import subprocess
import tempfile
from pathlib import Path


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
    print("NWKIT optimal roots and LCA loss exports verified.")


if __name__ == "__main__":
    main()
