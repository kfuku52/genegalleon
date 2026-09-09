#!/usr/bin/env python3
"""Require resident IQ-TREE likelihoods to agree with the installed CLI."""

import json
import tempfile
from pathlib import Path

import numpy as np
from ete4 import Tree
from nwkit.iqtree_library import select_worker
from nwkit.radte_inputs import build_chronology
from nwkit.radte_iqtree import IQTreeLikelihood
from nwkit.reconcile import build_reconciliation_table


def main():
    identity = select_worker(interface="library")
    gene = Tree("((A_1:0.1,B_1:0.1)S1:0.1,(A_2:0.2,B_2:0.2)S2:0.2)D;", parser=1)
    species = Tree("(A:10,B:10)AB;", parser=1)
    table = build_reconciliation_table(
        gene, species, {node.name: node.name.split("_")[0] for node in gene.leaves()}
    )
    chronology = build_chronology(gene, species, table, None, 30)
    rng = np.random.default_rng(83)
    ancestor = rng.choice(list("ACGT"), 200)
    with tempfile.TemporaryDirectory(prefix="gg-iqtree-smoke-") as directory:
        alignment = Path(directory) / "alignment.fa"
        records = []
        for name in ("A_1", "A_2", "B_1", "B_2"):
            sequence = ancestor.copy()
            sequence[rng.choice(200, 30, replace=False)] = rng.choice(list("ACGT"), 30)
            records.append(f">{name}\n{''.join(sequence)}\n")
        alignment.write_text("".join(records))
        library = IQTreeLikelihood(chronology, alignment, "JC", interface="auto")
        cli = None
        try:
            if library.interface != "library-worker-v1":
                raise RuntimeError("Default IQ-TREE interface did not select the installed worker")
            cli = IQTreeLikelihood(chronology, alignment, "JC", interface="cli")
            pid = library.worker.process.pid
            library.alignment.unlink()
            base = np.array([edge.dist for edge in chronology.edges])
            for scale in (1.0, 1.13, 0.87):
                actual = library.evaluate(base * scale)
                expected = cli.evaluate(base * scale)
                np.testing.assert_allclose(actual[1], expected[1], rtol=0, atol=2e-5)
                np.testing.assert_allclose(
                    actual[2][actual[4]], expected[2][expected[4]], rtol=1e-5, atol=3e-4
                )
                np.testing.assert_allclose(
                    actual[3][actual[4]], np.diag(expected[3])[expected[4]],
                    rtol=1e-5, atol=3e-3,
                )
                if library.worker.process.pid != pid or library.worker.process.poll() is not None:
                    raise RuntimeError("IQ-TREE worker did not remain resident")
        finally:
            library.close()
            if cli is not None:
                cli.close()
    print(json.dumps({"worker": identity, "evaluations": 3, "cli_agreement": True}, sort_keys=True))


if __name__ == "__main__":
    main()
