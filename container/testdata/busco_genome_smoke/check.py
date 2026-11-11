#!/usr/bin/env python3
"""Run BUSCO genome mode offline against a synthetic single-marker lineage."""

import gzip
import json
import random
import subprocess
import tempfile
from pathlib import Path


def main():
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    codons = "GCT TGT GAT GAA TTT GGT CAT ATT AAA CTG ATG AAT CCT CAA CGT TCT ACT GTT TGG TAT".split()
    genetic_code = dict(zip(amino_acids, codons, strict=True))
    generator = random.Random(743)
    protein = "M" + "".join(generator.choices(amino_acids, k=219))
    scaffold = "ACGT" * 50 + "".join(genetic_code[aa] for aa in protein) + "TAA" + "ACGT" * 50

    with tempfile.TemporaryDirectory(prefix="gg-busco-smoke-") as directory:
        root = Path(directory)
        lineage = root / "smoke_odb12"
        (lineage / "hmms").mkdir(parents=True)
        (lineage / "dataset.cfg").write_text(
            "name=smoke_odb12\ndomain=eukaryota\ncreation_date=2026-01-01\n"
            "number_of_BUSCOs=1\nnumber_of_species=1\n", encoding="ascii",
        )
        (lineage / "scores_cutoff").write_text("1at1 30\n", encoding="ascii")
        (lineage / "lengths_cutoff").write_text(f"1at1 {len(protein)} 10\n", encoding="ascii")
        marker = root / "marker.faa"
        marker.write_text(f">1at1\n{protein}\n", encoding="ascii")
        with gzip.open(lineage / "refseq_db.faa.gz", "wt", encoding="ascii") as handle:
            handle.write(f">1at1_smoke\n{protein}\n")
        genome = root / "genome.fa"
        genome.write_text(f">scaffold\n{scaffold}\n", encoding="ascii")

        subprocess.run(
            ["hmmbuild", "-n", "1at1", str(lineage / "hmms/1at1.hmm"), str(marker)],
            cwd=root, check=True, timeout=30,
        )
        subprocess.run(
            ["busco", "-i", str(genome), "-l", str(lineage), "-o", "result",
             "-m", "genome", "-c", "1", "--offline", "--opt-out-run-stats"],
            cwd=root, check=True, timeout=120,
        )
        summary_file, = (root / "result").glob("short_summary.specific.*.json")
        summary = json.loads(summary_file.read_text(encoding="utf-8"))
        metrics = summary["metrics"]
        if int(metrics["Number of scaffolds"]) != 1 or int(metrics["Scaffold N50"]) != len(scaffold):
            raise ValueError(f"BUSCO did not record the expected assembly statistics: {metrics}")
        if summary["results"]["Complete BUSCOs"] != 1:
            raise ValueError(f"BUSCO did not recover the synthetic marker: {summary['results']}")


if __name__ == "__main__":
    main()
