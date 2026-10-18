import copy
import gzip
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pytest

from workflow.support.curate_gff_utr_coordinates import curate
from workflow.support.format_species_annotation.common import reverse_complement


def fixture(tmp_path, strand="+"):
    dna = "CCCATGAAGTCCAGACCCTAAGGG"
    # First CDS block has five bases, so the second requires phase 1.
    features = [("gene", 1, 24, ".", "ID=gene1"),
                ("mRNA", 1, 24, ".", "ID=tx1;Parent=gene1"),
                ("exon", 1, 8, ".", "ID=e1;Parent=tx1"),
                ("CDS", 4, 8, "0", "ID=c1;Parent=tx1"),
                ("exon", 15, 24, ".", "ID=e2;Parent=tx1"),
                ("CDS", 15, 21, "1", "ID=c2;Parent=tx1"),
                ("five_prime_UTR", 1, 3, ".", "ID=u5;Parent=tx1"),
                ("three_prime_UTR", 2, 24, ".", "ID=u3;Parent=tx1")]
    if strand == "-":
        dna = reverse_complement(dna)
        features = [(kind, 25 - end, 25 - start, phase, attrs)
                    for kind, start, end, phase, attrs in features]
    gff = tmp_path / "original.gff"
    gff.write_text("".join(f"chr1\tsource\t{kind}\t{start}\t{end}\t.\t{strand}\t{phase}\t{attrs}\n"
                           for kind, start, end, phase, attrs in features))
    genome, cds, transcripts = [tmp_path / name for name in ("genome.fa", "cds.fa", "rna.fa")]
    genome.write_text(">chr1\n" + dna + "\n")
    cds.write_text(">tx1 published CDS\nATGAAACCCTAA\n")
    transcripts.write_text(">tx1 independent RNA\nCCCATGAAACCCTAAGGG\n")
    proposal = tmp_path / "proposal.json"
    change = dict(line=8, transcript_id="tx1", feature_id="u3",
                  old_start=2, old_end=24, new_start=22, new_end=24)
    if strand == "-":
        change.update(old_start=1, old_end=23, new_start=1, new_end=3)
    proposal.write_text(json.dumps(dict(schema_version=1, gff_sha256=hashlib.sha256(gff.read_bytes()).hexdigest(),
                                        changes=[change])))
    outdir = tmp_path / "curated"
    outdir.mkdir()
    return dict(gff=gff, genome=genome, cds=cds, transcripts=transcripts,
                proposal=proposal, outfile=outdir / "Species_a.gff.gz", audit=outdir / "proof.json")


def change_source(args, key, transform):
    path = args[key]
    path.write_text(transform(path.read_text()))
    if key == "gff":
        proposal = json.loads(args["proposal"].read_text())
        proposal["gff_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
        args["proposal"].write_text(json.dumps(proposal))


@pytest.mark.parametrize("strand", ["+", "-"])
def test_sequence_supported_utr_correction_retains_cds_and_passes_native_strict_reader(tmp_path, strand):
    args = fixture(tmp_path, strand)
    before = {key: path.read_bytes() for key, path in args.items() if key not in {"outfile", "audit"}}
    result = curate(**args)
    assert result["cds_file_unchanged"] and result["cds_exon_rna_gene_features_unchanged"]
    assert result["models"]["tx1"]["cds_rna_offset"] == 3
    assert result["models"]["tx1"]["cds_length"] == 12
    assert len(result["implementation_sha256"]) == 64
    assert result == json.loads(args["audit"].read_text())
    assert all(path.read_bytes() == before[key] for key, path in args.items() if key in before)
    with gzip.open(args["outfile"], "rt") as source:
        lines = source.readlines()
    original = args["gff"].read_text().splitlines(keepends=True)
    assert lines[:7] == original[:7]
    assert lines[7].split("\t")[:3] == original[7].split("\t")[:3]
    assert lines[7].split("\t")[5:] == original[7].split("\t")[5:]
    assert hashlib.sha256(args["outfile"].read_bytes()).hexdigest() == result["output_sha256"]
    formatted = tmp_path / "formatted.fa"
    formatted.write_text(">Species_a_gene1\nATGAAACCCTAA\n")
    completed = subprocess.run([sys.executable, str(Path(__file__).resolve().parents[1] / "support/gff2genestat.py"),
                                "--dir_gff", str(args["outfile"].parent), "--seqfile", str(formatted),
                                "--outfile", str(tmp_path / "native.tsv"), "--validate-cds-length", "--require-matches"],
                               capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr + completed.stdout
    assert "available" in (tmp_path / "native.tsv").read_text()


@pytest.mark.parametrize(("key", "before", "after", "reason"), [
    ("transcripts", "CCCATGAA", "CCCCTGAA", "Published RNA differs"),
    ("cds", "ATGAAACCCTAA", "ATGAAACCCTAG", "Supplied CDS differs"),
    ("gff", "CDS\t15\t21\t.\t+\t1", "CDS\t15\t21\t.\t+\t2", "CDS phases"),
    ("gff", "mRNA\t1\t24", "mRNA\t1\t23", "RNA bounds"),
    ("gff", "CDS\t15\t21", "CDS\t16\t21", "CDS is not contiguous"),
    ("cds", ">tx1 published", ">different_id published", "Missing exact transcript"),
])
def test_inconsistent_source_evidence_publishes_nothing(tmp_path, key, before, after, reason):
    args = fixture(tmp_path)
    change_source(args, key, lambda text: text.replace(before, after))
    with pytest.raises(ValueError, match=reason):
        curate(**args)
    assert not args["outfile"].exists() and not args["audit"].exists()


@pytest.mark.parametrize("defect", ["wrong_interval", "wrong_original", "wrong_hash", "duplicate_change",
                                     "shared_parent", "duplicate_rna", "duplicate_sequence"])
def test_ambiguous_or_changed_proposals_remain_rejected(tmp_path, defect):
    args = fixture(tmp_path)
    if defect in {"shared_parent", "duplicate_rna"}:
        change_source(args, "gff", lambda text: text.replace("ID=u3;Parent=tx1", "ID=u3;Parent=tx1,tx2")
                      if defect == "shared_parent" else text + text.splitlines(keepends=True)[1])
    elif defect == "duplicate_sequence":
        change_source(args, "cds", lambda text: text + text)
    else:
        proposal = json.loads(args["proposal"].read_text())
        if defect == "wrong_interval":
            proposal["changes"][0]["new_start"] = 23
        elif defect == "wrong_original":
            proposal["changes"][0]["old_start"] = 3
        elif defect == "wrong_hash":
            proposal["gff_sha256"] = "0" * 64
        else:
            proposal["changes"].append(copy.deepcopy(proposal["changes"][0]))
        args["proposal"].write_text(json.dumps(proposal))
    with pytest.raises(ValueError):
        curate(**args)
    assert not args["outfile"].exists() and not args["audit"].exists()


def test_output_replay_and_protected_input_overlap_are_rejected(tmp_path):
    args = fixture(tmp_path)
    invalid = {**args, "outfile": args["cds"]}
    original = args["cds"].read_bytes()
    with pytest.raises(ValueError, match="new output paths"):
        curate(**invalid)
    assert args["cds"].read_bytes() == original
    curate(**args)
    with pytest.raises(ValueError, match="new output paths"):
        curate(**args)


def test_internal_uncoded_exon_segment_does_not_become_an_invented_utr(tmp_path):
    args = fixture(tmp_path)
    args["genome"].write_text(">chr1\nCCCATGCCCAAGTCCAGACCCTAAGGG\n")
    args["transcripts"].write_text(">tx1\nCCCATGCCCAAACCCTAAGGG\n")
    change_source(args, "gff", lambda text: (
        "chr1\tsource\tgene\t1\t27\t.\t+\t.\tID=gene1\n"
        "chr1\tsource\tmRNA\t1\t27\t.\t+\t.\tID=tx1;Parent=gene1\n"
        "chr1\tsource\texon\t1\t11\t.\t+\t.\tID=e1;Parent=tx1\n"
        "chr1\tsource\tCDS\t4\t6\t.\t+\t0\tID=c1a;Parent=tx1\n"
        "chr1\tsource\tCDS\t10\t11\t.\t+\t0\tID=c1b;Parent=tx1\n"
        "chr1\tsource\texon\t18\t27\t.\t+\t.\tID=e2;Parent=tx1\n"
        "chr1\tsource\tCDS\t18\t24\t.\t+\t1\tID=c2;Parent=tx1\n"
        "chr1\tsource\tfive_prime_UTR\t1\t3\t.\t+\t.\tID=u5;Parent=tx1\n"
        "chr1\tsource\tthree_prime_UTR\t2\t27\t.\t+\t.\tID=u3;Parent=tx1\n"))
    request = json.loads(args["proposal"].read_text())
    request["changes"][0].update(line=9, old_end=27, new_start=25, new_end=27)
    args["proposal"].write_text(json.dumps(request))
    with pytest.raises(ValueError, match="CDS is not contiguous"):
        curate(**args)
    assert not args["outfile"].exists() and not args["audit"].exists()


def test_changed_source_after_sequence_proof_does_not_publish(tmp_path, monkeypatch):
    import workflow.support.curate_gff_utr_coordinates as module

    args = fixture(tmp_path)
    original = module.model_proof

    def change_after_proof(*positional):
        proof = original(*positional)
        args["cds"].write_text(args["cds"].read_text() + "\n")
        return proof

    monkeypatch.setattr(module, "model_proof", change_after_proof)
    with pytest.raises(ValueError, match="input changed"):
        curate(**args)
    assert not args["outfile"].exists() and not args["audit"].exists()


def test_cli_requires_exact_proposal_and_retains_failure(tmp_path):
    args = fixture(tmp_path)
    completed = subprocess.run([sys.executable, str(Path(__file__).resolve().parents[1] / "support/curate_gff_utr_coordinates.py"),
                                *[word for key, value in args.items() for word in ("--" + key, str(value))]],
                               capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    assert json.loads(completed.stdout)["status"] == "sequence_verified_utr_curation"
