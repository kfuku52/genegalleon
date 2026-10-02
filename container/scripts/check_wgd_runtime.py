#!/usr/bin/env python3
"""Verify the owning dependencies expose the required native WGD contracts."""

import csv
import importlib
import inspect
import json
import subprocess
import tempfile
from pathlib import Path


def main():
    from kffractbias.io import Gene, annotation_locus_order
    from kffractbias.selfevidence import write_self_evidence
    from kffractbias.selfscan import align_self, scan_self

    for module in ("ksrate", "ksrate_cli", "ksrate_model", "wgd_count", "wgd_count_cli",
                   "wgd_count_fit", "wgd_count_model", "wgd_tree", "wgd_tree_cli", "wgd_tree_model"):
        importlib.import_module("nwkit." + module)
    if "sequence_type" not in inspect.signature(align_self).parameters:
        raise ValueError("kfFractBias lacks protein self-alignment support.")
    if not {"screening", "allow_empty"}.issubset(inspect.signature(scan_self).parameters):
        raise ValueError("kfFractBias lacks explicit unquota/empty self evidence.")
    with tempfile.TemporaryDirectory(prefix="gg-wgd-contract-") as temporary:
        root = Path(temporary)
        gff = root / "loci.gff3"
        gff.write_text("chr1\tt\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1\n"
                       "chr1\tt\tgene\t11\t19\t.\t+\t.\tID=intervening\n"
                       "chr1\tt\tmRNA\t21\t29\t.\t+\t.\tID=t2;Parent=g2\n")
        if [gene.gene_id for gene in annotation_locus_order(gff, feature="mRNA", attribute="ID")] != ["g1", "intervening", "g2"]:
            raise ValueError("kfFractBias lost unsequenced annotation loci.")
        genes = tuple(Gene(chrom, 0, 10, name) for chrom, name in (("chr1", "a"), ("chr2", "b")))
        for name in ("self.self.anchors", "self.self.lifted.anchors"):
            (root / name).write_text("###\na\tb\t10\n")
        paths = write_self_evidence([[(0, 1, 10)]], genes, root)
        with paths["depth"].open() as handle:
            depths = list(csv.DictReader(handle, delimiter="\t"))
        if [int(row["block_arm_depth"]) for row in depths] != [1, 1]:
            raise ValueError("kfFractBias raw block-arm coverage is incompatible.")
        if json.loads(paths["summary"].read_text())["gene_span_coverage"] != 1:
            raise ValueError("kfFractBias raw coverage summary is incompatible.")
    for command in ("ksrate", "wgd-count", "wgd-tree"):
        subprocess.run(["nwkit", command, "--help"], check=True, stdout=subprocess.DEVNULL)
    print("Native WGD imports, CLI commands, full locus ordering and raw self evidence verified.")


if __name__ == "__main__":
    main()
