import csv
import gzip
import random
import subprocess
import sys
from pathlib import Path

import pytest
from Bio.Data.CodonTable import unambiguous_dna_by_id
from Bio.SeqIO.FastaIO import SimpleFastaParser

SUPPORT = Path(__file__).resolve().parents[1] / "support"


@pytest.mark.parametrize(("assembler", "input_ids", "group"), [
    ("trinity", ("TRINITY_DN1_c0_g1_i1", "DN1_c0_g1_i2"), "DN1-c0-g1"),
    ("rnaspades", ("NODE_1_length_360_cov_10_g1_i1", "NODE_2_length_330_cov_10_g1_i2"), "g1"),
])
def test_catalog_rescue_adds_missing_family_and_preserves_existing_cds(tmp_path, assembler, input_ids, group):
    species = "Test_species"
    group = f"{species}_{group}"
    normalized = (f"{group}-i1", f"{group}-i2")
    random_numbers = random.Random(123)
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    peptides = ["M" + "".join(random_numbers.choices(alphabet, k=n)) for n in (119, 109)]
    codon_for = {aa: codon for codon, aa in unambiguous_dna_by_id[1].forward_table.items()}
    cds = ["".join(codon_for[aa] for aa in peptide) + "TAA" for peptide in peptides]
    assembly = tmp_path / "assembly.fa.gz"
    with gzip.open(assembly, "wt") as stream:
        for name, seq in zip(input_ids, cds, strict=True):
            stream.write(f">{name} annotation\n{seq}\n")
    baseline = tmp_path / "baseline.fa"
    baseline.write_text("".join(f">{name}\n{seq}\n" for name, seq in zip(normalized, cds, strict=True)))
    alternates = tmp_path / "alternates.fa"
    alternates.touch()
    targets = tmp_path / "targets.faa"
    targets.write_text("".join(f">family{i}\n{seq}\n" for i, seq in enumerate(peptides)))
    db = tmp_path / "targets.dmnd"
    subprocess.run(
        ["diamond", "makedb", "--in", str(targets), "--db", str(db), "--threads", "1"],
        check=True, capture_output=True,
    )
    output, selection = tmp_path / "selected.fa", tmp_path / "selection.tsv"
    rescue, groups = tmp_path / "rescue.tsv", tmp_path / "groups.tsv"
    subprocess.run(
        [sys.executable, str(SUPPORT / "transcriptome_homology.py"),
         "--baseline", str(baseline), "--alternates", str(alternates),
         "--database", str(db), "--expression", r"\-i[0-9].*", "--genetic-code", "1", "--cpus", "1",
         "--output", str(output), "--report", str(selection),
         "--split-transcripts", str(assembly), "--split-species", species,
         "--split-assembler", assembler, "--split-targets", str(targets),
         "--orf-transcripts", str(baseline),
         "--split-extra-baseline", str(baseline),
         "--rescue-report", str(rescue), "--rescue-groups", str(groups)],
        check=True, capture_output=True,
    )
    with rescue.open(newline="") as stream:
        added, = (row for row in csv.DictReader(stream, delimiter="\t") if row["decision"] == "add")
    assert added["parent_group_id"] == group
    assert added["transcript_id"] == normalized[1]
    rescued_id = added["candidate_cds_id"]
    assert rescued_id.startswith(f"{group}.split")
    published = subprocess.run(
        [sys.executable, str(SUPPORT / "rnaspades_group_ids.py"),
         "rename-selected", str(groups), str(output)],
        check=True, capture_output=True, text=True,
    )
    assert list(SimpleFastaParser(iter(published.stdout.splitlines()))) == [(group, cds[0]), (rescued_id, cds[1])]
