"""Exercise domain ranking and the native MAFFT stage through both sequence modes."""

import csv
import gzip
import json
import os
import subprocess
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Data.CodonTable import unambiguous_dna_by_id
from Bio.Seq import Seq
from cdskit.alignment_io import write_fasta

from workflow.support.filter_empty_protein_alignment import filter_alignment

SUPPORT = Path(__file__).resolve().parents[1] / "support"
QUERY, DUPLICATE, PARTIAL, REPRESENTATIVE, CONFLICT = (
    "Homo_sapiens_query", "Homo_sapiens_copy", "Danio_rerio_partial",
    "Xenopus_laevis_full", "Arabidopsis_thaliana_alternate",
)


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_protein_trimming_filter_preserves_coordinates_and_partial_ambiguities(tmp_path):
    source = tmp_path / "clipkit.fa"
    output = tmp_path / "trimmed.fa"
    excluded = tmp_path / "excluded.tsv"
    source.write_text(">missing description\n-XU?O*\n>partial description\n-bZ.J-\n>one-residue\n-----A\n")

    assert filter_alignment(source, output, excluded) == (2, 1)
    assert output.read_text() == ">partial description\n-bZ.J-\n>one-residue\n-----A\n"
    assert read_tsv(excluded) == [{"gene_id": "missing", "reason": "no_usable_residues_after_trimming"}]


@pytest.mark.parametrize("mode", ("protein", "cds"))
def test_native_backbone_stage_recovers_all_residues_and_reuses_verified_alignment(tmp_path, mode):
    output = tmp_path / "output"
    common = "MKTFFVALKQLISEAISQLDRDLAETFGWVKL"
    sequences = {
        PARTIAL: common[:15] + "GGWV" + common[15:],
        QUERY: common,
        REPRESENTATIVE: common.replace("TFF", "TFY"),
        DUPLICATE: common + "A",
    }
    protein = tmp_path / "proteins.fa.gz"
    query = tmp_path / "query.fa"
    with gzip.open(protein, "wt") as handle:
        for gene, sequence in sequences.items():
            handle.write(f">{gene}\n{sequence}\n")
    write_fasta(query, sequences, [QUERY])
    inverse = {aa: codon for codon, aa in sorted(unambiguous_dna_by_id[1].forward_table.items())}
    cds = tmp_path / "cds.fa"
    write_fasta(cds, {g: "".join(inverse[a] for a in s) for g, s in sequences.items()}, sequences)
    rps = tmp_path / "rps.tsv"
    with rps.open("w") as h:
        w = csv.writer(h, delimiter="\t", lineterminator="\n")
        w.writerow(("qacc", "sacc", "qstart", "qend", "sstart", "send", "bitscore"))
        for gene in sequences:
            end = 4 if gene == PARTIAL else 10
            w.writerow((gene, "A", 1, end, 1, end, 50))
    screen = tmp_path / "screen.tsv"
    screen.write_text("gene_id\tstatus\tretained\n" + "".join(
        f"{gene}\ttarget_supported\t1\n" for gene in (REPRESENTATIVE, DUPLICATE)
    ))
    alignment = output / "mafft/TEST_cds.aln.fa.gz"
    variables = dict(os.environ, gg_support_dir=str(SUPPORT), gg_workspace_dir=str(tmp_path),
                     dir_output_active=str(output), og_id="TEST", mode_gene_evolution="query2family",
                     input_sequence_mode=mode, file_og_primary_fasta=str(protein if mode == "protein" else cds),
                     file_og_pep_fasta=str(protein), file_og_cds_fasta=str(cds), file_og_mafft=str(alignment),
                     file_og_query_aa_fasta=str(query), file_og_rpsblast=str(rps),
                     genetic_code="1", run_mafft="1", mafft_strategy="query_backbone", mafft_backbone_size="2",
                     query_family_filter="swissprot_interpro", query_family_filter_report=str(screen),
                     species_label_parser="extended-taxonomic", species_label_regex="", species_label_map_tsv="",
                     GG_TASK_CPUS="1", artifact_stale_policy="stop")
    source = (SUPPORT.parent / "core/gg_gene_evolution_core.sh").read_text()
    start = source.index('task="In-frame mafft alignment"')
    end = source.index('task="cdskit stats for original alignment"', start)
    script = r"""
set -euo pipefail
source "${gg_support_dir}/gg_util.sh"
gg_artifact_server_request() { local program=$1; shift; python "${program}" "$@"; }
gg_step_start() { printf 'started=%s\n' "$1"; }
gg_step_skip() { printf 'skipped=%s\n' "$1"; }
""" + source[source.index("prepare_gene_family_rpsblast_options()"):source.index('task="Protein RPS-BLAST"')] + source[start:end]
    first = subprocess.run(["bash", "-c", script], env=variables, cwd=tmp_path,
                           text=True, capture_output=True, timeout=90)
    assert first.returncode == 0, first.stdout + first.stderr

    with gzip.open(alignment, "rt") as h:
        rows = {record.id: str(record.seq).upper() for record in SeqIO.parse(h, "fasta")}
    assert list(rows) == list(sequences)
    for gene, aligned in rows.items():
        protein_row = str(Seq(aligned).translate(gap="-")) if mode == "cds" else aligned
        assert protein_row.replace("-", "") == sequences[gene]
    audit = read_tsv(output / "mafft_backbone/TEST_mafft_backbone.tsv")
    assert {row["gene_id"] for row in audit if row["backbone_reason"] != "added"} == {QUERY, REPRESENTATIVE}
    columns = read_tsv(output / "mafft_backbone/TEST_mafft_backbone_columns.tsv")
    assert audit[1]["family_screen"] == ""
    assert not (tmp_path / "downloads/pfam/cdskit").exists()
    assert len(columns) == len(rows[QUERY]) // (3 if mode == "cds" else 1)
    assert any(not row["backbone_column"] for row in columns)
    manifest = output / "artifact_provenance/TEST.mafft.json"
    receipt = json.loads(manifest.read_text())
    assert receipt["parameters"]["mafft_strategy"] == "query_backbone"
    assert {row["label"] for row in receipt["outputs"]} == {"mafft", "backbone_selection", "backbone_columns"}
    before = alignment.read_bytes()
    again = subprocess.run(["bash", "-c", script], env=variables, cwd=tmp_path,
                           text=True, capture_output=True, timeout=90)
    assert again.returncode == 0, again.stdout + again.stderr
    assert "skipped=In-frame mafft alignment" in again.stdout
    assert alignment.read_bytes() == before

    incompatible = subprocess.run(["bash", "-c", script], env=dict(variables, mafft_strategy="standard"),
                                  cwd=tmp_path, text=True, capture_output=True, timeout=90)
    assert incompatible.returncode != 0
    assert alignment.read_bytes() == before
