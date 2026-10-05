import csv
import gzip
import io
import json
import os
import shutil
import sqlite3
import ssl
import subprocess
import sys
import tarfile
from http.client import RemoteDisconnected
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from urllib.error import URLError

import pytest

SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "format_species_inputs.py"
VALIDATE_MAPPING_SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "validate_cds_gff_mapping.py"
VALIDATE_LONGEST_SCRIPT_PATH = Path(__file__).resolve().parents[1] / "support" / "validate_longest_cds_selection.py"
SMALL_DATASET_ROOT = Path(__file__).resolve().parent / "data" / "small_gfe_dataset"


@pytest.mark.parametrize("strand,expected", [("+", "ATGAAATAA"), ("-", "TTATTTCAT")])
def test_normalisation_indexes_reference_lengths_once(tmp_path, monkeypatch, strand, expected):
    load_module()
    import cds_model_normalisation
    import pysam

    genome, gff = tmp_path / "genome.fa", tmp_path / "source.gff"
    genome.write_text(">chr1\nATGAAATAA\n>chr2\nACGT\n")
    gff.write_text("##gff-version 3\n")
    real_fasta = pysam.FastaFile
    calls = {"references": 0, "lengths": 0}

    class CountingFasta:
        def __init__(self, path):
            self.delegate = real_fasta(path)

        @property
        def references(self):
            calls["references"] += 1
            return self.delegate.references

        @property
        def lengths(self):
            calls["lengths"] += 1
            return self.delegate.lengths

        def fetch(self, *args):
            return self.delegate.fetch(*args)

        def close(self):
            self.delegate.close()

    monkeypatch.setattr(pysam, "FastaFile", CountingFasta)
    normaliser = cds_model_normalisation.CdsModelNormaliser(
        {"genome": str(genome), "gff": str(gff)}, tmp_path, "test")
    row = {"seqid": "chr1", "strand": strand, "start": 0, "end": 9}
    try:
        for _ in range(3):
            assert normaliser.fetch([row]) == expected
        assert calls == {"references": 1, "lengths": 1}
        for changed in ({"seqid": "missing"}, {"start": -1}, {"end": 10}, {"end": 0}):
            with pytest.raises(ValueError, match="outside anchor genome"):
                normaliser.fetch([{**row, **changed}])
    finally:
        normaliser.close()
    assert not list(tmp_path.glob(".anchor-genome-*"))


@pytest.mark.parametrize("provider", ["direct", "local", "coge"])
def test_structured_coge_headers_keep_gene_identity_under_any_transport(tmp_path, provider):
    module = load_module()
    raw = tmp_path / "Species_one.cds.fa"
    raw.write_text(
        ">Species one||chr1||1||6||model-1-mRNA-1||1||CDS||101||1\nATGAAA\n"
        ">Species one||chr1||1||9||model-1-mRNA-2||1||CDS||102||2\nATGAAATTT\n"
        ">Species one||chr1||20||25||model-2-mRNA-1||-1||CDS||103||3\nATGCCC\n")
    task = dict(provider=provider, species_key="Species_one", species_prefix="Species_one",
                cds_path=raw, gff_path=None, genome_path=None)
    out = tmp_path / "out"
    out.mkdir()
    result = module.format_cds(task, out, overwrite=False, dry_run=False)
    assert result["before_count"] == 3
    assert result["after_count"] == 2
    with gzip.open(result["output_path"], "rt") as handle:
        assert handle.read() == ">Species_one_model-1\nATGAAATTT\n>Species_one_model-2\nATGCCC\n"


@pytest.mark.parametrize("provider", ["direct", "local", "coge"])
@pytest.mark.parametrize("model,start", [("", "1"), ("model1", "bad"), ("model1", "0")])
def test_direct_coge_malformed_identity_is_never_collapsed_to_species(tmp_path, model, start, provider):
    module = load_module()
    raw = tmp_path / "Species_one.cds.fa"
    raw.write_text(f">Species one||chr1||{start}||9||{model}||1||CDS||101||1\nATGAAATTT\n")
    task = dict(provider=provider, species_key="Species_one", species_prefix="Species_one",
                cds_path=raw, gff_path=None, genome_path=None)
    with pytest.raises(ValueError, match="Malformed CoGe CDS header"):
        module.format_cds(task, tmp_path / "out", overwrite=False, dry_run=False)
    assert not (tmp_path / "out").exists()


def test_header_only_cds_identifier_sanitization_collision_fails(tmp_path):
    module = load_module()
    raw = tmp_path / "Species_one.cds.fa"
    raw.write_text(
        ">Species one||chr1||1||6||gene:a-mRNA-1||1||CDS||101||1\nATGAAA\n"
        ">Species one||chr2||1||6||gene_a-mRNA-1||1||CDS||102||2\nATGCCC\n")
    task = dict(provider="direct", species_key="Species_one", species_prefix="Species_one",
                cds_path=raw, gff_path=None, genome_path=None)
    with pytest.raises(ValueError, match="Distinct CDS source IDs collide"):
        module.format_cds(task, tmp_path / "out", overwrite=False, dry_run=False)


@pytest.mark.parametrize("provider", ["direct", "local", "coge"])
def test_structured_coge_gene_sanitization_collision_fails_across_different_isoform_ids(tmp_path, provider):
    module = load_module()
    raw = tmp_path / "Species_one.cds.fa"
    raw.write_text(">Species one||chr1||1||6||gene:a-mRNA-1||1||CDS||101||1\nATGAAA\n"
                   ">Species one||chr2||1||6||gene_a-mRNA-2||1||CDS||102||2\nATGCCC\n")
    task = dict(provider=provider, species_key="Species_one", species_prefix="Species_one",
                cds_path=raw, gff_path=None, genome_path=None)
    with pytest.raises(ValueError, match="Distinct CoGe genes collide"):
        module.format_cds(task, tmp_path / "out", overwrite=False, dry_run=False)


@pytest.mark.parametrize("mode", ["strict", "rescue_overlap"])
def test_supplied_cds_uses_explicit_rna_parent_without_inventing_cds_coordinates(tmp_path, mode):
    module = load_module()
    gff = tmp_path / "source.gff"
    gff.write_text(
        "chr1\ts\tgene\t1\t20\t.\t+\t.\tID=g1\n"
        "chr1\ts\tmRNA\t1\t20\t.\t+\t.\tID=t1;Parent=g1;Accession=RNA1\n"
        "chr1\ts\tCDS\t1\t6\t.\t+\t0\tParent=t1\n"
        "chr1\ts\tmRNA\t1\t20\t.\t+\t.\tID=t2;Parent=g1;Accession=RNA2\n")
    raw = tmp_path / "Species_one.cds.fa"
    raw.write_text(">RNA1\nATGAAA\n>RNA2\nATGAAATTT\n")
    task = dict(provider="direct", species_key="Species_one", species_prefix="Species_one",
                cds_path=raw, gff_path=gff, genome_path=None, gene_grouping_mode=mode)
    index = module.build_gff_cds_grouping_index(task)
    assert index["transcript_gene_tokens"]["t2"] == "g1"
    assert len(index["location_to_gene_tokens"]) == 1
    (tmp_path / "out").mkdir()
    result = module.format_cds(task, tmp_path / "out", overwrite=False, dry_run=False)
    assert result["after_count"] == 1
    with gzip.open(result["output_path"], "rt") as handle:
        assert handle.read() == ">Species_one_g1\nATGAAATTT\n"


def normalisation_bundle(tmp_path, models, *, provided=True, code=1):
    from Bio.Seq import Seq

    cds, gff, genome = tmp_path / "source.cds.fa", tmp_path / "source.gff3", tmp_path / "source.genome.fa"
    features, sequences, supplied = ["##gff-version 3"], [], []
    for number, model in enumerate(models, 1):
        key, strand = f"t{number}", model.get("strand", "+")
        blocks = model.get("blocks", [(model.get("cds", "ATGAAATAA"), 0)])
        utr5, utr3 = model.get("utr5", ""), model.get("utr3", "")
        parts = [value for value, _ in blocks]
        parts[0] = utr5 + parts[0]
        parts[-1] += utr3
        dna = "N" * 7
        intervals = []
        for index, (part, (coding, phase)) in enumerate(zip(parts, blocks, strict=True)):
            start = len(dna)
            intervals.append((start, start + len(part), start + (len(utr5) if index == 0 else 0),
                              start + (len(utr5) if index == 0 else 0) + len(coding), phase))
            dna += part + "N" * 7
        if strand == "-":
            dna = str(Seq(dna).reverse_complement())
        sequences.append(f">chr{number}\n{dna}\n")
        def row(kind, start, end, phase, attrs, dna=dna, strand=strand, number=number):
            if strand == "-":
                start, end = len(dna) - end, len(dna) - start
            return f"chr{number}\tsynthetic\t{kind}\t{start + 1}\t{end}\t.\t{strand}\t{phase}\t{attrs}"
        features += [row("gene", 0, len(dna), ".", f"ID=g{number}"),
                     row("mRNA", 0, len(dna), ".", f"ID={key};Parent=g{number}" + model.get("attributes", ""))]
        for start, end, left, right, phase in intervals:
            features += [row("CDS", left, right, phase, f"Parent={key}"), row("exon", start, end, ".", f"Parent={key}")]
        raw = model.get("supplied", "".join(value for value, _ in blocks))
        supplied.append(f">{key}\n{raw}\n")
    cds.write_text("".join(supplied))
    gff.write_text("\n".join(features) + "\n")
    genome.write_text("".join(sequences))
    return dict(provider="direct", species_key="Test_species", species_prefix="Test_species",
                cds_path=cds if provided else None, gff_path=gff, genome_path=genome,
                gene_grouping_mode="strict", genetic_code=code)


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("provided", [False, True])
def test_formatting_normalises_partial_cds_and_paired_gff_without_rescue_flag(tmp_path, strand, provided, monkeypatch):
    from Bio.Seq import Seq

    monkeypatch.setenv("GG_INPUT_RUN_GENE_MODEL_RESCUE", "0")
    module = load_module()
    task = normalisation_bundle(tmp_path, [{}, {"blocks": [("TAAT", 1), ("GAAACCCTAAC", 2)], "strand": strand}], provided=provided)
    before = {name: Path(task[name]).read_bytes() for name in ("gff_path", "genome_path")}
    if provided:
        before["cds_path"] = task["cds_path"].read_bytes()
    output = tmp_path / "output"
    output.mkdir()
    result = module.format_cds(task, output, False, False)
    assert result["cds_normalised_records"] == 1
    records = dict(module.iter_fasta_records(result["output_path"]))
    assert records["Test_species_g2"] == "ATGAAACCCTAA"
    paired = module.format_gff(task, output, False, False, formatted_cds_path=result["output_path"])
    with gzip.open(paired["output_path"], "rt") as handle:
        rows = [line.rstrip().split("\t") for line in handle if "\tCDS\t" in line and "Parent=t2" in line]
    rows.sort(key=lambda row: int(row[3]), reverse=strand == "-")
    assert [int(row[7]) for row in rows] == [0, 1]
    genome = {header: sequence for header, sequence in module.iter_fasta_records(task["genome_path"])}
    pieces = [genome[row[0]][int(row[3]) - 1:int(row[4])] for row in rows]
    genomic = "".join(str(Seq(part).reverse_complement()) if strand == "-" else part for part in pieces)
    assert genomic == records["Test_species_g2"]
    assert all("gg_cds_normalisation=annotated_partial_frame" in row[8] for row in rows)
    assert {name: Path(task[name]).read_bytes() for name in before} == before
    audit = json.loads(Path(str(result["output_path"]) + ".cds-normalisation.json").read_text())
    assert audit["records"][0]["selected_evidence"][0]["trailing_bases_omitted"] == 1
    assert not list(tmp_path.glob("*.fai"))


@pytest.mark.parametrize("explicit", [False, True])
def test_cds_normalisation_uses_task_scratch_for_genome_reconstruction(tmp_path, monkeypatch, explicit):
    import tempfile

    # Load process-wide annotation caches before changing this task's scratch.
    # The test must also work in isolation, without an earlier module import.
    module = load_module()
    import cds_model_normalisation

    runtime = tmp_path / "runtime"
    runtime.mkdir()
    override = tmp_path / "override"
    override.mkdir()
    monkeypatch.setenv("TMPDIR", str(runtime))
    monkeypatch.setattr(tempfile, "tempdir", None)
    selected = override if explicit else runtime
    original = tempfile.TemporaryDirectory
    created = []

    def temporary_directory(*args, **kwargs):
        if kwargs.get("prefix") == ".anchor-genome-":
            assert Path(kwargs["dir"]) == selected
            created.append(Path(kwargs["dir"]))
        return original(*args, **kwargs)

    monkeypatch.setattr(cds_model_normalisation.tempfile, "TemporaryDirectory", temporary_directory)
    task = normalisation_bundle(tmp_path, [{"blocks": [("TAAT", 1), ("GAAACCCTAAC", 2)]}])
    if explicit:
        task["_normalisation_scratch"] = override
    output = tmp_path / "output"
    output.mkdir()
    result = module.format_cds(task, output, False, False)
    assert dict(module.iter_fasta_records(result["output_path"])) == {"Test_species_g1": "ATGAAACCCTAA"}
    assert created == [selected]
    assert not list(selected.iterdir())


@pytest.mark.parametrize("strand", ["+", "-"])
def test_formatting_removes_proven_utr_and_retains_unresolved_originals(tmp_path, strand):
    module = load_module()
    task = normalisation_bundle(tmp_path, [
        {"cds": "ATGAAATAA", "utr5": "TAA", "supplied": "TAAATGAAATAA", "strand": strand},
        {"cds": "ATGTGACCCTAA", "attributes": ";transl_except=(pos:4..6,aa:OTHER)"},
        {"blocks": [("TAATGAAATAA", 1)]}])
    result = module.format_cds(task, tmp_path, False, False)
    records = dict(module.iter_fasta_records(result["output_path"]))
    assert records == {"Test_species_g1": "ATGAAATAA", "Test_species_g2": "ATGTGACCCTAA", "Test_species_g3": "TAATGAAATAAN"}
    assert (result["cds_normalised_records"], result["cds_unresolved_records"]) == (1, 2)
    audit = json.loads(Path(str(result["output_path"]) + ".cds-normalisation.json").read_text())
    assert audit["counts"] == {"normalised": 1, "retained_unresolved": 2}


def test_formatting_chooses_longest_after_utr_correction(tmp_path):
    module = load_module()
    task = normalisation_bundle(tmp_path, [{"cds": "ATGAAATAA", "utr5": "TAATAATAA", "supplied": "TAATAATAAATGAAATAA"}])
    with task["cds_path"].open("a") as handle:
        handle.write(">t2\nATGCCCCCCTAA\n")
    genome = task["genome_path"].read_text().splitlines()[1]
    start = len(genome) + 1
    task["genome_path"].write_text(">chr1\n" + genome + "ATGCCCCCCTAA\n")
    lines = task["gff_path"].read_text().splitlines()
    gene = lines[1].split("\t")
    gene[4] = str(len(genome) + 12)
    lines[1] = "\t".join(gene)
    lines += [f"chr1\tsynthetic\tmRNA\t{start}\t{start + 11}\t.\t+\t.\tID=t2;Parent=g1",
              f"chr1\tsynthetic\tCDS\t{start}\t{start + 11}\t.\t+\t0\tParent=t2"]
    task["gff_path"].write_text("\n".join(lines) + "\n")
    result = module.format_cds(task, tmp_path, False, False)
    assert dict(module.iter_fasta_records(result["output_path"])) == {"Test_species_g1": "ATGCCCCCCTAA"}
    with Path(str(result["output_path"]) + ".gff-grouping.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["raw_sequence_length"] == "18" and rows[0]["effective_cds_length"] == "9"
    assert rows[1]["selected_longest"] == "1"
    from validate_longest_cds_selection import collect_expected_longest_records
    assert collect_expected_longest_records(task)[0]["Test_species_g1"]["sequence"] == "ATGCCCCCCTAA"


def test_formatting_normalisation_audit_controls_reuse_and_genetic_code(tmp_path):
    module = load_module()
    task = normalisation_bundle(tmp_path, [{"cds": "ATGTGACCCTAA"}], code=4)
    first = module.format_cds(task, tmp_path, False, False)
    assert first["cds_unresolved_records"] == 0
    assert module.format_cds(task, tmp_path, False, False)["status"] == "skip"
    task["genetic_code"] = 1
    assert module.format_cds(task, tmp_path, False, False)["cds_unresolved_records"] == 1
    audit_path = Path(str(first["output_path"]) + ".cds-normalisation.json")
    audit = json.loads(audit_path.read_text())
    assert audit["contract"]["genetic_code"] == 1
    task["genome_path"].write_text(task["genome_path"].read_text().replace("ATGTGA", "ATGAGA"))
    assert module.format_cds(task, tmp_path, False, False)["status"] == "write"


def test_formatting_dry_run_does_not_publish_corrections_or_indexes(tmp_path):
    module = load_module()
    task = normalisation_bundle(tmp_path, [{"blocks": [("TAAT", 1), ("GAAATAA", 2)]}])
    output = tmp_path / "output"
    output.mkdir()
    assert module.format_cds(task, output, False, True)["status"] == "dry-run"
    assert list(output.iterdir()) == [] and not list(tmp_path.glob("*.fai"))


@pytest.mark.parametrize("provided", [True, False])
def test_formatting_refuses_unverified_paired_corrections_and_rebuilds(tmp_path, provided):
    module = load_module()
    task = normalisation_bundle(tmp_path, [{"blocks": [("TAAT", 1), ("GAAACCCTAAC", 2)]}], provided=provided)
    first = module.format_cds(task, tmp_path, False, False)
    paired = module.format_gff(task, tmp_path, False, False, formatted_cds_path=first["output_path"])
    output_before = paired["output_path"].read_bytes()
    audit = Path(str(first["output_path"]) + ".cds-normalisation.json")
    audit_before = audit.stat()
    audit.write_text("{broken")
    with pytest.raises(ValueError, match="Missing or stale CDS normalisation audit"):
        module.format_gff(task, tmp_path, False, False, formatted_cds_path=first["output_path"])
    assert paired["output_path"].read_bytes() == output_before
    assert module.format_cds(task, tmp_path, False, False)["status"] == "write"
    # Filesystems with coarse timestamps can recreate the identical audited pair
    # within one tick. Explicitly invalidate the recorded sidecar fingerprint.
    os.utime(audit, ns=(audit_before.st_atime_ns, audit_before.st_mtime_ns + 2_000_000_000))
    rebuilt = module.format_gff(task, tmp_path, False, False, formatted_cds_path=first["output_path"])
    assert rebuilt["status"] == "write"
    with gzip.open(rebuilt["output_path"], "rt") as handle:
        rebuilt_lines = handle.read()
    assert rebuilt_lines == gzip.decompress(output_before).decode()
    audit.unlink()
    with pytest.raises(ValueError, match="Missing or stale CDS normalisation audit"):
        module.format_gff(task, tmp_path, True, False, formatted_cds_path=first["output_path"])


def test_formatting_normalises_against_archived_genome_with_original_seqid(tmp_path):
    module = load_module()
    task = normalisation_bundle(tmp_path, [{"blocks": [("TAAT", 1), ("GAAACCCTAAC", 2)]}])
    dna = task["genome_path"].read_text().splitlines()[1]
    contents = f">acc1 OriSeqID=chr1 Len={len(dna)}\n{dna}\n".encode()
    archive = tmp_path / "genome.fa.tar.gz"
    with tarfile.open(archive, "w:gz") as handle:
        member = tarfile.TarInfo("genome.fa")
        member.size = len(contents)
        handle.addfile(member, io.BytesIO(contents))
    task["genome_path"] = archive
    result = module.format_cds(task, tmp_path, False, False)
    paired = module.format_gff(task, tmp_path, False, False, formatted_cds_path=result["output_path"])
    assert result["cds_normalised_records"] == 1
    with gzip.open(paired["output_path"], "rt") as handle:
        assert all(line.split("\t")[0] == "acc1" for line in handle if not line.startswith("#"))
    assert not list(tmp_path.glob("*.fai"))


def test_formatting_bad_genome_coordinates_fail_before_cds_publication(tmp_path):
    module = load_module()
    task = normalisation_bundle(tmp_path, [{"blocks": [("TAAT", 1), ("GAAACCCTAAC", 2)]}])
    task["genome_path"].write_text(">chr1\nA\n")
    sources = {name: Path(task[name]).read_bytes() for name in ("cds_path", "gff_path", "genome_path")}
    output = tmp_path / "output"
    output.mkdir()
    with pytest.raises(ValueError, match="outside anchor genome|exceeds genome"):
        module.format_cds(task, output, False, False)
    assert list(output.iterdir()) == []
    assert {name: Path(task[name]).read_bytes() for name in sources} == sources


def test_formatting_drops_unused_missing_sequence_region_declarations(tmp_path):
    module = load_module()
    task = normalisation_bundle(tmp_path, [{}])
    task["gff_path"].write_text("##sequence-region removed_plastid 1 150718\n" + task["gff_path"].read_text())
    result = module.format_cds(task, tmp_path, False, False)
    paired = module.format_gff(task, tmp_path, False, False, formatted_cds_path=result["output_path"])
    with gzip.open(paired["output_path"], "rt") as handle:
        assert "removed_plastid" not in handle.read()
    task["gff_path"].write_text(task["gff_path"].read_text().replace("chr1", "removed_plastid"))
    with pytest.raises(ValueError, match="absent from genome"):
        module.format_gff(task, tmp_path, True, False, formatted_cds_path=result["output_path"])


@pytest.mark.parametrize("header,source_id,canonical", [
    ("acc1 OriSeqID=Chr1 Len=9", "Chr1", "acc1"),
    ("lcl|chr1", "chr1", "lcl|chr1"),
    ("evm.model.chr1", "evm.model.chr1", "chr1"),
])
@pytest.mark.parametrize("with_cds", [False, True])
def test_formatted_gff_references_match_exported_genome(tmp_path, header, source_id, canonical, with_cds):
    module = load_module()
    genome = tmp_path / "genome.fa"
    genome.write_text(">" + header + "\nATGAAATTT\n")
    gff = tmp_path / "source.gff"
    gff.write_text("##sequence-region " + source_id + " 1 9\n" + source_id + "\tsrc\tCDS\t1\t9\t.\t+\t0\tID=gene1;Parent=gene1\n")
    cds = tmp_path / "selected.fa"
    cds.write_text(">Test_species_gene1\nATGAAATTT\n")
    task = dict(provider="direct", species_prefix="Test_species", gff_path=gff, genome_path=genome,
                gene_grouping_mode="strict", gff_repair_mode="safe")
    result = module.format_gff(task, tmp_path, False, False, formatted_cds_path=cds if with_cds else None)
    with gzip.open(result["output_path"], "rt") as handle:
        lines = handle.read().splitlines()
    assert lines[0] == "##sequence-region " + canonical + " 1 9"
    assert lines[1].split("\t")[0] == canonical
    from format_species_annotation.reference import validate_gff_genome_references
    output_genome = module.format_genome(task, tmp_path, False, False)["output_path"]
    assert validate_gff_genome_references(result["output_path"], output_genome) == 1


@pytest.mark.parametrize("headers,message", [
    (">chr2\nATG\n", "absent from genome"),
    (">acc1 OriSeqID=chr1 Len=3\nATG\n>acc2 OriSeqID=chr1 Len=3\nATG\n", "Ambiguous"),
    (">chr1\nATG\n>chr1\nATG\n", "duplicate"),
    (">acc1 OriSeqID=chr1 Len=4\nATG\n", "length disagrees"),
])
def test_reference_normalization_rejects_unresolved_or_ambiguous_genomes(tmp_path, headers, message):
    load_module()
    from format_species_annotation.reference import gff_reference_mapping
    genome = tmp_path / "genome.fa"
    genome.write_text(headers)
    gff = tmp_path / "source.gff"
    gff.write_text("chr1\tsrc\tCDS\t1\t3\t.\t+\t0\tParent=gene1\n")
    with pytest.raises(ValueError, match=message):
        gff_reference_mapping(gff, genome)


@pytest.mark.parametrize("provided_cds", [False, True])
def test_overlap_rescue_preserves_distinct_declared_gene_parents_without_gene_rows(tmp_path, provided_cds):
    module = load_module()
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\n" + "ATG" * 10 + "\n")
    gff = tmp_path / "source.gff"
    models = [("499.g6.t1_499.g7.t1", "499.g6_499.g7", 18), ("499.g7.t1.1.hash", "499.g7", 12)]
    gff.write_text("".join(
        f"chr1\tsrc\tmRNA\t1\t{end}\t.\t+\t.\tID={tid};Parent={parent}\n"
        f"chr1\tsrc\tCDS\t1\t{end}\t.\t+\t0\tParent={tid}\n" for tid, parent, end in models))
    task = dict(provider="direct", species_prefix="Test_species", species_key="Test_species",
                gff_path=gff, genome_path=genome, gene_grouping_mode="rescue_overlap", gff_repair_mode="safe")
    if provided_cds:
        cds = tmp_path / "source.cds"
        cds.write_text("".join(f">{tid}\n" + "ATG" * (end // 3) + "\n" for tid, parent, end in models))
        task["cds_path"] = cds
    result = module.format_cds(task, tmp_path, False, False, strict=True)
    assert dict(module.iter_fasta_records(result["output_path"])) == {
        "Test_species_499.g6_499.g7": "ATG" * 6,
        "Test_species_499.g7": "ATG" * 4,
    }


def test_genome_reference_index_reads_multiple_archived_fasta_members(tmp_path):
    load_module()
    from format_species_annotation.reference import genome_reference_index
    archive_path = tmp_path / "genome.fa.tar.bz2"
    with tarfile.open(archive_path, "w:bz2") as archive:
        for name, text in (("nested/one.fa", ">acc1 OriSeqID=Chr1 Len=3\nATG\n"),
                           ("two.fna", ">acc2\nAAACCC\n"), ("README.txt", "metadata\n")):
            data = text.encode()
            member = tarfile.TarInfo(name)
            member.size = len(data)
            archive.addfile(member, io.BytesIO(data))
    index = genome_reference_index(archive_path)
    assert index["Chr1"] == index["acc1"] == 3
    assert index["acc2"] == 6
    assert index.canonical_ids["Chr1"] == "acc1"


def load_module():
    spec = spec_from_file_location("format_species_inputs_module", SCRIPT_PATH)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_script(*args, env=None):
    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH), *args],
        capture_output=True,
        text=True,
        check=False,
        env=env,
    )


def run_validate_mapping_script(*args):
    return subprocess.run(
        [sys.executable, str(VALIDATE_MAPPING_SCRIPT_PATH), *args],
        capture_output=True,
        text=True,
        check=False,
    )


def run_validate_longest_script(*args):
    return subprocess.run(
        [sys.executable, str(VALIDATE_LONGEST_SCRIPT_PATH), *args],
        capture_output=True,
        text=True,
        check=False,
    )


def test_cli_default_outputs_honor_isolated_output_root(monkeypatch, tmp_path):
    module = load_module()
    output_root = tmp_path / "isolated-output"
    monkeypatch.setenv("GG_INPUT_GENERATION_OUTPUT_ROOT", str(output_root))

    args = module.build_arg_parser().parse_args(["--provider", "local"])

    assert Path(args.species_cds_dir) == output_root / "species_cds"
    assert Path(args.species_gff_dir) == output_root / "species_gff"
    assert Path(args.species_genome_dir) == output_root / "species_genome"
    assert Path(args.species_summary_output) == output_root / "gg_input_generation_species.tsv"


def test_provider_resolvers_share_one_registry_interface():
    module = load_module()

    assert "coge" in module.PROVIDER_ID_RESOLVERS
    assert "figshare" in module.PROVIDER_ID_RESOLVERS
    assert all(hasattr(resolver, "resolve") for resolver in module.PROVIDER_ID_RESOLVERS.values())


def test_transient_network_error_treats_ssl_eof_as_retryable():
    module = load_module()

    error = URLError(ssl.SSLEOFError(8, "UNEXPECTED_EOF_WHILE_READING"))

    assert module.is_transient_network_error(error)


def test_cds_extension_is_treated_as_fasta():
    module = load_module()

    assert module.is_fasta_filename("GZX_Primary.cds")
    assert module.is_fasta_filename("GZX_Primary.cds.gz")
    assert not module.is_fasta_filename("GZX_Primary.gff")
    assert (
        module.normalize_cds_output_basename("GZX_Primary.cds", "Fakus_species")
        == "Fakus_species_GZX_Primary.cds.fa.gz"
    )
    assert (
        module.normalize_cds_output_basename("GZX_Primary.cds.gz", "Fakus_species")
        == "Fakus_species_GZX_Primary.cds.fa.gz"
    )


def test_coge_cds_header_uses_delimited_feature_id():
    module = load_module()
    headers = (
        "Gilia yorkii||Gy1||47683||49946||GY000001-RA||-1||CDS||3560718418||1",
        "Gilia yorkii||Gy1||92931||101228||GY000002-RA||-1||CDS||3560718421||2",
        (
            "Sarracenia purpurea ||CM117427.1_RagTag||149705||152539||"
            "evm.model.CM117427.1_RagTag.1||-1||CDS||4559620697||1"
        ),
    )

    assert [module.extract_provider_id("coge", header) for header in headers] == [
        "GY000001-RA",
        "GY000002-RA",
        "evm.model.CM117427.1_RagTag.1",
    ]
    gilia_task = {"provider": "coge", "species_prefix": "Gilia_yorkii"}
    assert module.build_formatted_cds_id(gilia_task, headers[0]) == "Gilia_yorkii_GY000001-RA"
    assert module.build_formatted_cds_id(gilia_task, headers[1]) == "Gilia_yorkii_GY000002-RA"


def test_coge_coding_export_requires_cds_features(tmp_path):
    module = load_module()
    gff_path = tmp_path / "gene_only.gff"
    gff_path.write_text("chr1\tCoGe\tgene\t1\t9\t.\t+\t.\tID=g1\n", encoding="utf-8")

    with pytest.raises(ValueError, match="coding-only GFF export with no CDS records"):
        module.validate_coge_export_gff_file(gff_path, gid="69349")


def test_format_coge_delimited_cds_headers_map_to_gff_genes(tmp_path):
    input_dir = tmp_path / "CoGe" / "species_wise_original"
    species_dir = input_dir / "Gilia_yorkii"
    species_dir.mkdir(parents=True)
    (species_dir / "Gilia_yorkii.coge.gid62042.cds.fa").write_text(
        (
            ">Gilia yorkii||Gy1||1||9||GY000001-RA||1||CDS||101||1\nATGAAATTT\n"
            ">Gilia yorkii||Gy1||20||28||GY000002-RA||1||CDS||102||2\nATGCCCTTT\n"
        ),
        encoding="utf-8",
    )
    sarracenia_dir = input_dir / "Sarracenia_purpurea"
    sarracenia_dir.mkdir()
    (sarracenia_dir / "Sarracenia_purpurea.coge.gid69349.cds.fa").write_text(
        (
            ">Sarracenia purpurea ||CM117427.1_RagTag||1||9||"
            "evm.model.CM117427.1_RagTag.1||1||CDS||201||1\nATGAAATTT\n"
        ),
        encoding="utf-8",
    )
    (sarracenia_dir / "Sarracenia_purpurea.gid69349.gff").write_text(
        "\n".join(
            [
                "##gff-version\t3",
                (
                    "CM117427.1_RagTag\tCoGe\tgene\t1\t9\t.\t+\t.\t"
                    "ID=EVM%20prediction%20CM117427.1_RagTag.1;"
                    "Name=EVM%20prediction%20CM117427.1_RagTag.1;"
                    "gene=EVM%20prediction%20CM117427.1_RagTag.1"
                ),
                (
                    "CM117427.1_RagTag\tCoGe\tmRNA\t1\t9\t.\t+\t.\t"
                    "ID=evm.model.CM117427.1_RagTag.1;"
                    "Parent=EVM%20prediction%20CM117427.1_RagTag.1"
                ),
                (
                    "CM117427.1_RagTag\tCoGe\tCDS\t1\t9\t.\t+\t0\t"
                    "ID=evm.model.CM117427.1_RagTag.1.cds;"
                    "Parent=evm.model.CM117427.1_RagTag.1;"
                    "Alias=evm.model.CM117427.1_RagTag.1"
                ),
                "",
            ]
        ),
        encoding="utf-8",
    )
    (species_dir / "Gilia_yorkii.gid62042.gff").write_text(
        "\n".join(
            [
                "##gff-version\t3",
                "Gy1\tCoGe\tgene\t1\t9\t.\t+\t.\tID=GY000001;Name=GY000001;gene=GY000001",
                "Gy1\tCoGe\tmRNA\t1\t9\t.\t+\t.\tID=GY000001.mRNA1;Parent=GY000001;Alias=GY000001-RA",
                "Gy1\tCoGe\tCDS\t1\t9\t.\t+\t0\tID=GY000001-RA;Parent=GY000001.mRNA1;CDS=GY000001-RA",
                "Gy1\tCoGe\tgene\t20\t28\t.\t+\t.\tID=GY000002;Name=GY000002;gene=GY000002",
                "Gy1\tCoGe\tmRNA\t20\t28\t.\t+\t.\tID=GY000002.mRNA1;Parent=GY000002;Alias=GY000002-RA",
                "Gy1\tCoGe\tCDS\t20\t28\t.\t+\t0\tID=GY000002-RA;Parent=GY000002.mRNA1;CDS=GY000002-RA",
                "",
            ]
        ),
        encoding="utf-8",
    )

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    completed = run_script(
        "--provider",
        "coge",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(tmp_path / "species_genome"),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Gilia_yorkii_coge.gid62042.cds.fa.gz"
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        headers = [line.strip() for line in handle if line.startswith(">")]
    assert headers == [">Gilia_yorkii_GY000001", ">Gilia_yorkii_GY000002"]
    sarracenia_cds = out_cds / "Sarracenia_purpurea_coge.gid69349.cds.fa.gz"
    with gzip.open(sarracenia_cds, "rt", encoding="utf-8") as handle:
        sarracenia_headers = [line.strip() for line in handle if line.startswith(">")]
    assert sarracenia_headers == [">Sarracenia_purpurea_EVM_prediction_CM117427.1_RagTag.1"]

    mapping = run_validate_mapping_script(
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
    )
    assert mapping.returncode == 0, mapping.stderr + "\n" + mapping.stdout
    assert "[Gilia_yorkii] CDS-to-GFF mapping OK: 2/2 IDs" in mapping.stdout
    assert "[Sarracenia_purpurea] CDS-to-GFF mapping OK: 1/1 IDs" in mapping.stdout


def test_local_gff_genome_lcl_prefix_is_resolved_when_deriving_cds(tmp_path):
    input_dir = tmp_path / "Local" / "species_wise_original"
    species_dir = input_dir / "Gilia_yorkii"
    species_dir.mkdir(parents=True)
    (species_dir / "Gilia_yorkii.genome.fa").write_text(
        ">lcl|Gy1\nATGAAATTT\n>lcl|Gy2\nATGCCCTTT\n",
        encoding="utf-8",
    )
    (species_dir / "Gilia_yorkii.gff").write_text(
        "Gy1\tCoGe\tgene\t1\t9\t.\t+\t.\tID=GY000001\n"
        "Gy1\tCoGe\tmRNA\t1\t9\t.\t+\t.\tID=GY000001.mRNA1;Parent=GY000001\n"
        "Gy1\tCoGe\tCDS\t1\t9\t.\t+\t0\tID=GY000001-RA;Parent=GY000001.mRNA1\n"
        "Gy2\tCoGe\tgene\t1\t9\t.\t+\t.\tID=GY000002\n"
        "Gy2\tCoGe\tmRNA\t1\t9\t.\t+\t.\tID=GY000002.mRNA1;Parent=GY000002\n"
        "Gy2\tCoGe\tCDS\t1\t9\t.\t+\t0\tID=GY000002-RA;Parent=GY000002.mRNA1\n",
        encoding="utf-8",
    )

    out_cds = tmp_path / "species_cds"
    completed = run_script(
        "--provider",
        "local",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(tmp_path / "species_gff"),
        "--species-genome-dir",
        str(tmp_path / "species_genome"),
    )

    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    formatted_cds = next(out_cds.glob("*.fa.gz"))
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert ">Gilia_yorkii_GY000001" in text
    assert ">Gilia_yorkii_GY000002" in text
    assert "ATGAAATTT" in text
    assert "ATGCCCTTT" in text


def test_local_coge_gff_percent_encoded_gene_ids_remain_distinct(tmp_path):
    input_dir = tmp_path / "Local" / "species_wise_original"
    species_dir = input_dir / "Sarracenia_purpurea"
    species_dir.mkdir(parents=True)
    genome_path = species_dir / "Sarracenia_purpurea.genome.fa"
    genome_path.write_text(
        ">lcl|CM117427.1_RagTag\nATGAAATTTATGCCCTTT\n",
        encoding="utf-8",
    )
    gff_path = species_dir / "Sarracenia_purpurea.gff"
    gff_path.write_text(
        "\n".join(
            [
                (
                    "CM117427.1_RagTag\tCoGe\tgene\t1\t9\t.\t+\t.\t"
                    "ID=EVM%20prediction%20CM117427.1_RagTag.1;"
                    "Name=EVM%20prediction%20CM117427.1_RagTag.1;"
                    "gene=EVM%20prediction%20CM117427.1_RagTag.1"
                ),
                (
                    "CM117427.1_RagTag\tCoGe\tmRNA\t1\t9\t.\t+\t.\t"
                    "ID=evm.model.CM117427.1_RagTag.1;"
                    "Parent=EVM%20prediction%20CM117427.1_RagTag.1"
                ),
                (
                    "CM117427.1_RagTag\tCoGe\tCDS\t1\t9\t.\t+\t0\t"
                    "ID=evm.model.CM117427.1_RagTag.1.cds;"
                    "Parent=evm.model.CM117427.1_RagTag.1"
                ),
                (
                    "CM117427.1_RagTag\tCoGe\tgene\t10\t18\t.\t+\t.\t"
                    "ID=EVM%20prediction%20CM117427.1_RagTag.2;"
                    "Name=EVM%20prediction%20CM117427.1_RagTag.2;"
                    "gene=EVM%20prediction%20CM117427.1_RagTag.2"
                ),
                (
                    "CM117427.1_RagTag\tCoGe\tmRNA\t10\t18\t.\t+\t.\t"
                    "ID=evm.model.CM117427.1_RagTag.2;"
                    "Parent=EVM%20prediction%20CM117427.1_RagTag.2"
                ),
                (
                    "CM117427.1_RagTag\tCoGe\tCDS\t10\t18\t.\t+\t0\t"
                    "ID=evm.model.CM117427.1_RagTag.2.cds;"
                    "Parent=evm.model.CM117427.1_RagTag.2"
                ),
                "",
            ]
        ),
        encoding="utf-8",
    )

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    summary_path = tmp_path / "gg_input_generation_species.tsv"
    completed = run_script(
        "--provider",
        "local",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(tmp_path / "species_genome"),
        "--species-summary-output",
        str(summary_path),
    )

    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    formatted_cds = next(out_cds.glob("*.fa.gz"))
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        headers = [line.strip() for line in handle if line.startswith(">")]
    assert headers == [
        ">Sarracenia_purpurea_EVM_prediction_CM117427.1_RagTag.1",
        ">Sarracenia_purpurea_EVM_prediction_CM117427.1_RagTag.2",
    ]

    mapping = run_validate_mapping_script(
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
    )
    assert mapping.returncode == 0, mapping.stderr + "\n" + mapping.stdout
    assert "[Sarracenia_purpurea] CDS-to-GFF mapping OK: 2/2 IDs" in mapping.stdout

    longest = run_validate_longest_script(
        "--species-cds-dir",
        str(out_cds),
        "--species-summary",
        str(summary_path),
    )
    assert longest.returncode == 0, longest.stderr + "\n" + longest.stdout
    assert "[Sarracenia_purpurea] Longest CDS validation OK:" in longest.stdout


def test_gff_genome_seqid_map_rejects_many_to_one_lcl_alias():
    module = load_module()

    with pytest.raises(ValueError, match="both resolve to FASTA sequence ID 'lcl\\|Gy1'"):
        module.build_gff_genome_seqid_map(
            {"lcl|Gy1": "ATGAAATTT"},
            {"Gy1", "lcl|Gy1"},
        )


def test_gff_genome_seqid_map_reports_unresolved_seqids():
    module = load_module()

    seqid_map, missing_seqids = module.build_gff_genome_seqid_map(
        {"lcl|Gy1": "ATGAAATTT", "Gy2": "ATGCCCTTT"},
        {"Gy1", "Gy2", "Gy3"},
    )

    assert seqid_map == {"Gy1": "lcl|Gy1", "Gy2": "Gy2"}
    assert missing_seqids == ("Gy3",)


def test_gff_genome_seqid_map_uses_verified_original_sequence_ids(tmp_path):
    module = load_module()
    genome = tmp_path / "genome.fa"
    genome.write_text(
        ">GWHEUWF00000001.1 Chromosome 1a OriSeqID=Chr1A Len=9\nATGAAATTT\n",
        encoding="utf-8",
    )
    sequences = module.load_genome_sequences(genome)

    seqid_map, missing = module.build_gff_genome_seqid_map(sequences, {"Chr1A"})

    assert missing == ()
    assert seqid_map == {"Chr1A": "GWHEUWF00000001.1"}
    with pytest.raises(ValueError, match="both resolve"):
        module.build_gff_genome_seqid_map(sequences, {"Chr1A", "GWHEUWF00000001.1"})


def test_gff_genome_seqid_map_rejects_incorrect_or_colliding_original_ids(tmp_path):
    module = load_module()
    genome = tmp_path / "genome.fa"
    genome.write_text(">acc1 OriSeqID=Chr1 Len=8\nATGAAATTT\n", encoding="utf-8")
    with pytest.raises(ValueError, match="length disagrees"):
        module.load_genome_sequences(genome)

    genome.write_text(
        ">acc1 OriSeqID=Chr1 Len=9\nATGAAATTT\n"
        ">acc2 OriSeqID=Chr1 Len=9\nATGCCCTTT\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="alias 'Chr1' is shared"):
        module.load_genome_sequences(genome)


def test_direct_ncbi_like_cds_header_uses_locus_tag_for_gene_grouping():
    module = load_module()
    task = {
        "provider": "direct",
        "species_prefix": "Vanilla_planifolia",
    }
    header = (
        "lcl|CM028150.1_cds_KAG0495310.1_1 "
        "[locus_tag=HPP92_000001] [protein_id=KAG0495310.1] [db_xref=NCBI_GP:KAG0495310.1]"
    )

    gene_id = module.build_gene_aggregate_id(
        task,
        header,
        "Vanilla_planifolia_lcl_CM028150.1_cds_KAG0495310.1_1",
    )

    assert gene_id == "Vanilla_planifolia_HPP92_000001"


def test_direct_species_discovery_treats_lone_chr_fasta_as_genome_when_gff_exists(tmp_path):
    module = load_module()
    input_dir = tmp_path / "Direct" / "species_wise_original"
    species_dir = input_dir / "Santalum_album"
    species_dir.mkdir(parents=True)
    (species_dir / "tanxiang.FINAL.chr.fa").write_text(">lcl|Chr01\nATGAAATTT\n", encoding="utf-8")
    (species_dir / "tanxiang.FINAL.chr_modified.gff").write_text(
        "Chr01\tGnomon\tmRNA\t1\t9\t.\t+\t.\tID=SA1G00001\n"
        "Chr01\tGnomon\tCDS\t1\t9\t.\t+\t0\tID=SA1G00001.cds1;Parent=SA1G00001\n",
        encoding="utf-8",
    )

    tasks, warnings, errors = module.discover_tasks("direct", input_dir)

    assert errors == []
    assert warnings == []
    assert len(tasks) == 1
    assert tasks[0]["cds_path"] is None
    assert tasks[0]["genome_path"].name == "tanxiang.FINAL.chr.fa"
    assert tasks[0]["gff_path"].name == "tanxiang.FINAL.chr_modified.gff"


def test_direct_species_discovery_keeps_lone_nonmatching_fasta_as_cds_when_gff_exists(tmp_path):
    module = load_module()
    input_dir = tmp_path / "Direct" / "species_wise_original"
    species_dir = input_dir / "Fakus_species"
    species_dir.mkdir(parents=True)
    (species_dir / "sequences.fa").write_text(">gene1\nATGAAA\n", encoding="utf-8")
    (species_dir / "annotation.gff3").write_text(
        "chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tID=cds1;Parent=gene1\n",
        encoding="utf-8",
    )

    tasks, warnings, errors = module.discover_tasks("direct", input_dir)

    assert errors == []
    assert warnings == []
    assert len(tasks) == 1
    assert tasks[0]["cds_path"].name == "sequences.fa"
    assert tasks[0]["genome_path"] is None
    assert tasks[0]["gff_path"].name == "annotation.gff3"


class FakeTextPipe:
    def __init__(self):
        self.parts = []

    def write(self, text):
        self.parts.append(text)
        return len(text)

    def close(self):
        return None

    def getvalue(self):
        return "".join(self.parts)


class FakeBinaryResponse:
    def __init__(self, payload):
        self._payload = payload

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        return False

    def read(self):
        return self._payload


def write_test_taxonomy_fixture(tmp_path):
    db_path = tmp_path / "taxa.sqlite"
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute(
        "CREATE TABLE species (taxid INT PRIMARY KEY, parent INT, spname VARCHAR(50) COLLATE NOCASE, common VARCHAR(50) COLLATE NOCASE, rank VARCHAR(50), track TEXT)"
    )
    cur.execute("CREATE TABLE synonym (taxid INT,spname VARCHAR(50) COLLATE NOCASE, PRIMARY KEY (spname, taxid))")
    species_rows = [
        (1, 1, "root", "", "no rank", "1"),
        (131567, 1, "cellular organisms", "", "cellular root", "131567,1"),
        (2759, 131567, "Eukaryota", "", "domain", "2759,131567,1"),
        (33090, 2759, "Viridiplantae", "", "kingdom", "33090,2759,131567,1"),
        (242159, 33090, "Ostreococcus lucimarinus", "", "species", "242159,33090,2759,131567,1"),
        (5786, 2759, "Dictyostelium discoideum", "", "species", "5786,2759,131567,1"),
        (5911, 2759, "Tetrahymena thermophila", "", "species", "5911,2759,131567,1"),
    ]
    cur.executemany(
        "INSERT INTO species (taxid, parent, spname, common, rank, track) VALUES (?, ?, ?, ?, ?, ?)", species_rows
    )
    cur.execute("INSERT INTO synonym (taxid, spname) VALUES (?, ?)", (242159, "Ostreococcus_lucimarinus"))
    conn.commit()
    conn.close()

    def nodes_line(taxid, parent, rank, gc_id, mito_gc_id):
        return "{}\t|\t{}\t|\t{}\t|\t\t|\t0\t|\t0\t|\t{}\t|\t0\t|\t{}\t|\t0\t|\t0\t|\t0\t|\t\t|\n".format(
            taxid, parent, rank, gc_id, mito_gc_id
        )

    gencode_text = (
        "1\t|\tSGC0\t|\tStandard\t|\t\t|\t\t|\n"
        "4\t|\tSGC4\t|\tMold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma\t|\t\t|\t\t|\n"
        "6\t|\tSGC6\t|\tCiliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear\t|\t\t|\t\t|\n"
        "11\t|\tSGC11\t|\tBacterial, Archaeal and Plant Plastid\t|\t\t|\t\t|\n"
    )
    nodes_text = "".join(
        [
            nodes_line(1, 1, "no rank", 1, 0),
            nodes_line(131567, 1, "cellular root", 1, 0),
            nodes_line(2759, 131567, "domain", 1, 0),
            nodes_line(33090, 2759, "kingdom", 1, 0),
            nodes_line(242159, 33090, "species", 1, 1),
            nodes_line(5786, 2759, "species", 1, 4),
            nodes_line(5911, 2759, "species", 6, 4),
        ]
    )

    taxdump_path = tmp_path / "taxdump.tar.gz"
    with tarfile.open(taxdump_path, "w:gz") as archive:
        for name, text in {
            "gencode.dmp": gencode_text,
            "nodes.dmp": nodes_text,
            "readme.txt": "test fixture\n",
        }.items():
            payload = text.encode("utf-8")
            info = tarfile.TarInfo(name=name)
            info.size = len(payload)
            archive.addfile(info, io.BytesIO(payload))

    return db_path, taxdump_path


def test_format_species_inputs_with_small_fixture_all_providers(tmp_path):
    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    stats_json = tmp_path / "stats.json"

    completed = run_script(
        "--provider",
        "all",
        "--input-dir",
        str(SMALL_DATASET_ROOT),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--stats-output",
        str(stats_json),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    ensembl_cds = out_cds / "Ostreococcus_lucimarinus_ASM9206v1.cds.all.fa.gz"
    phycocosm_cds = out_cds / "Microglena_spYARC_MicrYARC1_GeneCatalog_CDS_20220803.fa.gz"
    phytozome_cds = out_cds / "Hydrocotyle_leucocephala_HleucocephalaHAP1_768_v2.1.cds_primaryTranscriptOnly.fa.gz"
    assert ensembl_cds.exists()
    assert phycocosm_cds.exists()
    assert phytozome_cds.exists()

    with gzip.open(ensembl_cds, "rt", encoding="utf-8") as handle:
        ensembl_text = handle.read()
    assert ensembl_text.count(">Ostreococcus_lucimarinus_OSTLU_25062") == 1
    assert ">Ostreococcus_lucimarinus_OSTLU_99999" in ensembl_text
    assert "ATGAAN" in ensembl_text

    with gzip.open(phycocosm_cds, "rt", encoding="utf-8") as handle:
        phycocosm_text = handle.read()
    assert ">Microglena_spYARC_mRNA.MigICE15955" in phycocosm_text
    assert ">Microglena_spYARC_mRNA.MigICE_00468" in phycocosm_text

    with gzip.open(phytozome_cds, "rt", encoding="utf-8") as handle:
        phytozome_text = handle.read()
    assert ">Hydrocotyle_leucocephala_HyleuH1.06G006800" in phytozome_text
    assert "ATGATGAN" in phytozome_text

    ensembl_gff = out_gff / "Ostreococcus_lucimarinus_ASM9206v1.56.gff.gz"
    phytozome_gff = out_gff / "Hydrocotyle_leucocephala_HleucocephalaHAP1_768_v2.1.gene.gff.gz"
    phytozome_exons = out_gff / "Hydrocotyle_leucocephala_HleucocephalaHAP1_768_v2.1.gene_exons.gff.gz"
    assert ensembl_gff.exists()
    assert phytozome_gff.exists()
    assert not phytozome_exons.exists()

    with gzip.open(ensembl_gff, "rt", encoding="utf-8") as handle:
        gff_text = handle.read()
    assert "evm.model." not in gff_text
    assert "Oropetium_20150105_" not in gff_text

    with gzip.open(phytozome_gff, "rt", encoding="utf-8") as handle:
        phytozome_gff_text = handle.read()
    assert "evm_27.model." not in phytozome_gff_text
    stats = json.loads(stats_json.read_text(encoding="utf-8"))
    assert stats["species_processed"] == 3
    assert stats["num_species_cds_files"] == 3
    assert stats["num_species_gff_files"] == 3
    assert stats["cds_sequences_before"] >= stats["cds_sequences_after"]
    assert stats["cds_first_sequence_name"] != ""


def test_species_taxonomy_metadata_resolver_supports_nonstandard_nuclear_codes(tmp_path):
    mod = load_module()
    db_path, taxdump_path = write_test_taxonomy_fixture(tmp_path)
    resolver = mod.SpeciesTaxonomyMetadataResolver(str(db_path), str(taxdump_path))

    tetrahymena = resolver.resolve("Tetrahymena thermophila")
    assert tetrahymena["taxid"] == "5911"
    assert tetrahymena["nuclear_genetic_code_id"] == "6"
    assert tetrahymena["nuclear_genetic_code_name"] == "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear"
    assert tetrahymena["mitochondrial_genetic_code_id"] == "4"
    assert tetrahymena["plastid_genetic_code_id"] == ""


def test_parse_species_key_candidate_preserves_taxonomic_qualifiers():
    mod = load_module()

    assert mod.parse_species_key_candidate("Dictyostelium cf. discoideum") == "Dictyostelium_cf_discoideum"
    assert mod.parse_species_key_candidate("Bacillus subtilis subsp. subtilis") == "Bacillus_subtilis_subsp_subtilis"
    assert mod.parse_species_key_candidate("Amoeba sp. JDS-Ruffled") == "Amoeba_sp_JDSRuffled"
    assert mod.parse_species_key_candidate("Amoeba sp.") == "Amoeba_sp_unknown"
    assert (
        mod.parse_species_key_candidate("Solanum lycopersicum cultivar Heinz 1706")
        == "Solanum_lycopersicum_cultivar_Heinz1706"
    )
    assert mod.parse_species_key_candidate("Escherichia coli serovar O157") == "Escherichia_coli_serovar_O157"
    assert mod.parse_species_key_candidate("Citrus x limon") == "Citrus_x_limon"
    assert (
        mod.parse_species_key_candidate("Cenchrus americanus x Cenchrus purpureus")
        == "Cenchrus_americanus_x_Cenchrus_purpureus"
    )
    assert (
        mod.parse_species_key_candidate("Cenchrus americanus \u00d7 Cenchrus purpureus")
        == "Cenchrus_americanus_x_Cenchrus_purpureus"
    )


def test_source_id_candidates_add_ensemblplants_species_and_accession_variants():
    mod = load_module()

    candidates = mod.source_id_candidates(
        "ensemblplants",
        "GCA_910589775.1",
        "Avena_eriantha",
    )
    assert "Avena_eriantha" in candidates
    assert "Avena_eriantha_gca910589775v1cm" in candidates

    candidates = mod.source_id_candidates(
        "ensemblplants",
        "GCA_000695525.1",
        "Brassica_oleracea_var._oleracea",
    )
    assert "Brassica_oleracea" in candidates

    candidates = mod.source_id_candidates(
        "ensemblplants",
        "GCA_001952365.2",
        "Oryza_sativa_aus_subgroup",
    )
    assert "Oryza_sativa_aus" in candidates


def test_source_id_candidates_add_ensemblmetazoa_accession_suffix_variants():
    mod = load_module()

    candidates = mod.source_id_candidates(
        "ensemblmetazoa",
        "GCA_003676215.3",
        "Rhopalosiphum_maidis",
    )

    assert "Rhopalosiphum_maidis" in candidates
    assert "Rhopalosiphum_maidis_gca003676215v3" in candidates


def test_expand_ensemblgenomes_id_candidates_prefers_matching_directory_suffix(monkeypatch):
    mod = load_module()
    provider_module = sys.modules[mod.expand_ensemblgenomes_id_candidates.__module__]

    monkeypatch.setattr(
        provider_module,
        "fetch_ensemblgenomes_dir_ids",
        lambda provider, timeout, headers: ["rhopalosiphum_maidis_gca003676215v3"],
    )

    candidates = mod.expand_ensemblgenomes_id_candidates(
        "ensemblmetazoa",
        ["Rhopalosiphum_maidis"],
        timeout=1,
        headers={},
    )

    assert candidates[0] == "rhopalosiphum_maidis_gca003676215v3"
    assert "Rhopalosiphum_maidis" in candidates


def test_infer_provider_from_ensemblgenomes_prefix_and_url():
    mod = load_module()

    assert mod.infer_provider_from_id("ensemblmetazoa:anopheles_gambiae") == "ensemblmetazoa"
    assert mod.infer_provider_from_id("ensemblprotists:phytophthora_parasitica") == "ensemblprotists"
    assert (
        mod.infer_provider_from_id(
            "https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/metazoa/fasta/anopheles_gambiae/cds/"
        )
        == "ensemblmetazoa"
    )
    assert (
        mod.infer_provider_from_id(
            "https://ftp.ensemblgenomes.ebi.ac.uk/pub/current/protists/gff3/phytophthora_parasitica/"
        )
        == "ensemblprotists"
    )


def test_discover_ncbi_like_tasks_normalizes_bare_sp_species_key(tmp_path):
    mod = load_module()
    input_dir = tmp_path / "NCBI_Genome" / "species_wise_original"
    species_dir = input_dir / "Amoeba_sp"
    species_dir.mkdir(parents=True, exist_ok=True)
    (species_dir / "GCA_000000001.1_demo_cds_from_genomic.fna.gz").write_text("", encoding="utf-8")
    (species_dir / "GCA_000000001.1_demo_genomic.gff.gz").write_text("", encoding="utf-8")
    (species_dir / "GCA_000000001.1_demo_genomic.fna.gz").write_text("", encoding="utf-8")

    tasks, warnings, errors = mod.discover_ncbi_like_tasks(input_dir, "ncbi")

    assert warnings == []
    assert errors == []
    assert len(tasks) == 1
    assert tasks[0]["species_key"] == "Amoeba_sp_unknown"
    assert tasks[0]["species_prefix"] == "Amoeba_sp_unknown"


def test_species_prefix_keeps_hybrid_marker_with_following_epithet():
    mod = load_module()
    assert mod.species_prefix_from_value("Citrus_x_limon") == "Citrus_x_limon"
    assert mod.species_prefix_from_value("Petunia_x_hybrida") == "Petunia_x_hybrida"
    assert (
        mod.species_prefix_from_value("Cenchrus_americanus_x_Cenchrus_purpureus.cds.fa.gz")
        == "Cenchrus_americanus_x_Cenchrus_purpureus"
    )


def test_species_prefix_preserves_dotted_taxonomic_qualifiers():
    mod = load_module()

    assert (
        mod.species_prefix_from_value("Asimitellaria_furusei_var._furusei_demo.fa.gz")
        == "Asimitellaria_furusei_var._furusei"
    )
    assert (
        mod.species_prefix_from_value("Asimitellaria_furusei_var._subramosa.cds.fa.gz")
        == "Asimitellaria_furusei_var._subramosa"
    )
    assert mod.species_prefix_from_value("Arisaema_sp._aooni_demo.fa.gz") == "Arisaema_sp._aooni"
    assert (
        mod.species_prefix_from_value("Bacillus_subtilis_subsp._subtilis_demo.gff.gz")
        == "Bacillus_subtilis_subsp._subtilis"
    )
    assert mod.species_prefix_from_value("homo_sapiens.GRCh38.cds.all.fa.gz") == "homo_sapiens"


def test_discover_ensembl_like_tasks_preserves_dotted_taxonomic_qualifiers(tmp_path):
    mod = load_module()
    input_dir = tmp_path / "Ensembl" / "original_files"
    input_dir.mkdir(parents=True)
    for species_key in (
        "Asimitellaria_furusei_var._furusei",
        "Asimitellaria_furusei_var._subramosa",
    ):
        (input_dir / f"{species_key}.ASM.cds.all.fa.gz").write_text("", encoding="utf-8")
        (input_dir / f"{species_key}.ASM.1.gff3.gz").write_text("", encoding="utf-8")
        (input_dir / f"{species_key}.ASM.dna.primary_assembly.fa.gz").write_text("", encoding="utf-8")

    tasks, warnings, errors = mod.discover_ensembl_like_tasks(input_dir, "ensembl")

    assert warnings == []
    assert errors == []
    assert [task["species_key"] for task in tasks] == [
        "Asimitellaria_furusei_var._furusei",
        "Asimitellaria_furusei_var._subramosa",
    ]
    assert [task["species_prefix"] for task in tasks] == [
        "Asimitellaria_furusei_var._furusei",
        "Asimitellaria_furusei_var._subramosa",
    ]


def test_discover_ncbi_like_tasks_rejects_incomplete_taxonomic_qualifier_species_key(tmp_path):
    mod = load_module()
    input_dir = tmp_path / "NCBI_Genome" / "species_wise_original"
    species_dir = input_dir / "Dictyostelium_cf"
    species_dir.mkdir(parents=True, exist_ok=True)
    (species_dir / "GCA_054859205.1_ASM5485920v1_cds_from_genomic.fna.gz").write_text("", encoding="utf-8")
    (species_dir / "GCA_054859205.1_ASM5485920v1_genomic.gff.gz").write_text("", encoding="utf-8")
    (species_dir / "GCA_054859205.1_ASM5485920v1_genomic.fna.gz").write_text("", encoding="utf-8")

    tasks, warnings, errors = mod.discover_ncbi_like_tasks(input_dir, "ncbi")

    assert tasks == []
    assert warnings == []
    assert len(errors) == 1
    assert "requires a following epithet" in errors[0]
    assert "Dictyostelium_cf" in errors[0]


def test_species_taxonomy_metadata_resolver_falls_back_from_qualified_name_to_base_species(tmp_path):
    mod = load_module()
    db_path, taxdump_path = write_test_taxonomy_fixture(tmp_path)
    resolver = mod.SpeciesTaxonomyMetadataResolver(str(db_path), str(taxdump_path))

    dicty = resolver.resolve("Dictyostelium_cf_discoideum")
    assert dicty["taxid"] == "5786"
    assert dicty["nuclear_genetic_code_id"] == "1"
    assert dicty["mitochondrial_genetic_code_id"] == "4"


def test_manifest_declared_providers_preserves_manifest_order_and_skips_unlisted_rows():
    mod = load_module()

    providers = mod.manifest_declared_providers(
        [
            {"provider": "ncbi", "id": "GCA_000000001.1"},
            {"provider": "direct", "id": "sample_direct"},
            {"provider": "NCBI", "id": "GCA_000000002.1"},
            {"provider": "oryza_minuta", "id": "gramene_tetraploids"},
            {"provider": "", "id": "missing_provider"},
            {"provider": "unsupported", "id": "unsupported_provider"},
            {"provider": "local", "id": "sample_local"},
        ],
        provider_filter="all",
    )

    assert providers == ["ncbi", "direct", "oryza_minuta", "local"]


def test_species_summary_includes_taxid_and_genetic_codes_when_taxonomy_cache_is_available(tmp_path):
    db_path, taxdump_path = write_test_taxonomy_fixture(tmp_path)
    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    species_summary = tmp_path / "gg_input_generation_species.tsv"
    env = dict(os.environ)
    env["GG_TAXONOMY_DBFILE"] = str(db_path)
    env["GG_TAXONOMY_TAXDUMPFILE"] = str(taxdump_path)

    completed = run_script(
        "--provider",
        "ensemblplants",
        "--input-dir",
        str(SMALL_DATASET_ROOT / "20230216_EnsemblPlants" / "original_files"),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
        "--species-summary-output",
        str(species_summary),
        env=env,
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    with open(species_summary, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    row = rows[0]
    assert row["species_prefix"] == "Ostreococcus_lucimarinus"
    assert row["taxid"] == "242159"
    assert row["nuclear_genetic_code_id"] == "1"
    assert row["nuclear_genetic_code_name"] == "Standard"
    assert row["mitochondrial_genetic_code_id"] == "1"
    assert row["mitochondrial_genetic_code_name"] == "Standard"
    assert row["plastid_genetic_code_id"] == "11"
    assert row["plastid_genetic_code_name"] == "Bacterial, Archaeal and Plant Plastid"


def test_format_species_inputs_strict_mode_accepts_cds_only_inputs(tmp_path):
    dataset_copy = tmp_path / "small_gfe_dataset_copy"
    shutil.copytree(SMALL_DATASET_ROOT, dataset_copy)
    missing_gff = (
        dataset_copy / "20230216_EnsemblPlants" / "original_files" / "Ostreococcus_lucimarinus.ASM9206v1.56.gff3"
    )
    missing_gff.unlink()

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    species_summary = tmp_path / "gg_input_generation_species.tsv"
    completed = run_script(
        "--provider",
        "ensemblplants",
        "--input-dir",
        str(dataset_copy / "20230216_EnsemblPlants" / "original_files"),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-summary-output",
        str(species_summary),
        "--strict",
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    assert any(path.name.endswith(".fa.gz") for path in out_cds.iterdir())
    assert not any(out_gff.iterdir())
    with open(species_summary, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["gff_status"] == "missing"
    assert rows[0]["genome_status"] == "missing"


def test_format_species_inputs_derives_cds_from_gff_and_genome_when_cds_is_missing(tmp_path):
    input_dir = tmp_path / "Direct" / "species_wise_original"
    species_dir = input_dir / "Arabidopsis_thaliana"
    species_dir.mkdir(parents=True, exist_ok=True)
    gff_path = species_dir / "Arabidopsis_thaliana.annotation.gff3"
    genome_path = species_dir / "Arabidopsis_thaliana.genome.fa"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene1",
                "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=gene1.t1;Parent=gene1",
                "chr1\tsrc\tCDS\t1\t3\t.\t+\t0\tID=cds1;Parent=gene1.t1",
                "chr1\tsrc\tCDS\t7\t9\t.\t+\t0\tID=cds2;Parent=gene1.t1",
                "chr2\tsrc\tgene\t1\t9\t.\t-\t.\tID=gene2",
                "chr2\tsrc\tmRNA\t1\t9\t.\t-\t.\tID=gene2.t1;Parent=gene2",
                "chr2\tsrc\tCDS\t1\t3\t.\t-\t0\tID=cds3;Parent=gene2.t1",
                "chr2\tsrc\tCDS\t7\t9\t.\t-\t0\tID=cds4;Parent=gene2.t1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    genome_path.write_text(
        ">chr1\nATGAAATTT\n>chr2\nTTTAAACAT\n",
        encoding="utf-8",
    )

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    species_summary = tmp_path / "gg_input_generation_species.tsv"
    completed = run_script(
        "--provider",
        "direct",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
        "--species-summary-output",
        str(species_summary),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Arabidopsis_thaliana_annotation.derived.cds.fa.gz"
    assert formatted_cds.exists()
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert text.count(">Arabidopsis_thaliana_gene1") == 1
    assert text.count(">Arabidopsis_thaliana_gene2") == 1
    assert "ATGTTT" in text
    assert "ATGAAA" in text

    with open(species_summary, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    row = rows[0]
    assert str(gff_path) in row["cds_input_path"]
    assert str(genome_path) in row["cds_input_path"]
    assert "derived CDS" in row["cds_input_path"]
    assert row["gff_input_path"] == str(gff_path)
    assert row["genome_input_path"] == str(genome_path)


def test_format_species_inputs_skips_mixed_strand_transcript_when_deriving_cds(tmp_path):
    input_dir = tmp_path / "Direct" / "species_wise_original"
    species_dir = input_dir / "Ophrys_sphegodes"
    species_dir.mkdir(parents=True, exist_ok=True)
    gff_path = species_dir / "Ophrys_sphegodes.annotation.gff3"
    genome_path = species_dir / "Ophrys_sphegodes.genome.fa"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tID=gene1",
                "chr1\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=gene1.t1;Parent=gene1",
                "chr1\tsrc\tCDS\t1\t3\t.\t+\t0\tID=cds1;Parent=gene1.t1",
                "chr1\tsrc\tCDS\t10\t12\t.\t-\t0\tID=cds2;Parent=gene1.t1",
                "chr1\tsrc\tgene\t20\t25\t.\t+\t.\tID=gene2",
                "chr1\tsrc\tmRNA\t20\t25\t.\t+\t.\tID=gene2.t1;Parent=gene2",
                "chr1\tsrc\tCDS\t20\t25\t.\t+\t0\tID=cds3;Parent=gene2.t1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    genome_path.write_text(">chr1\nATGAAATTTCCCGGGAAATGAAATAA\n", encoding="utf-8")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "direct",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    assert "skipping transcript 'gene1.t1'" in completed.stderr

    formatted_cds = out_cds / "Ophrys_sphegodes_annotation.derived.cds.fa.gz"
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert ">Ophrys_sphegodes_gene1" not in text
    assert ">Ophrys_sphegodes_gene2" in text


def test_format_species_inputs_skips_mixed_strand_transcript_before_utr_trimming(tmp_path):
    input_dir = tmp_path / "Direct" / "species_wise_original"
    species_dir = input_dir / "Ophrys_sphegodes"
    species_dir.mkdir(parents=True, exist_ok=True)
    gff_path = species_dir / "Ophrys_sphegodes.annotation.gff3"
    genome_path = species_dir / "Ophrys_sphegodes.genome.fa"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tID=gene1",
                "chr1\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=gene1.t1;Parent=gene1",
                "chr1\tsrc\tfive_prime_UTR\t1\t3\t.\t+\t.\tParent=gene1.t1",
                "chr1\tsrc\tCDS\t1\t3\t.\t+\t0\tID=cds1;Parent=gene1.t1",
                "chr1\tsrc\tCDS\t10\t12\t.\t-\t0\tID=cds2;Parent=gene1.t1",
                "chr1\tsrc\tgene\t20\t25\t.\t+\t.\tID=gene2",
                "chr1\tsrc\tmRNA\t20\t25\t.\t+\t.\tID=gene2.t1;Parent=gene2",
                "chr1\tsrc\tCDS\t20\t25\t.\t+\t0\tID=cds3;Parent=gene2.t1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    genome_path.write_text(">chr1\nATGAAATTTCCCGGGAAATGAAATAA\n", encoding="utf-8")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "direct",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    assert "skipping transcript 'gene1.t1'" in completed.stderr

    formatted_cds = out_cds / "Ophrys_sphegodes_annotation.derived.cds.fa.gz"
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert ">Ophrys_sphegodes_gene1" not in text
    assert ">Ophrys_sphegodes_gene2" in text


def test_format_species_inputs_treats_softmasked_direct_fasta_as_genome(tmp_path):
    input_dir = tmp_path / "Direct" / "species_wise_original"
    species_dir = input_dir / "Arabidopsis_thaliana"
    species_dir.mkdir(parents=True, exist_ok=True)
    gff_path = species_dir / "Arabidopsis_thaliana.genes.gff3"
    genome_path = species_dir / "Arabidopsis_thaliana.softmasked.fa"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene1",
                "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=gene1.t1;Parent=gene1",
                "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=gene1.t1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    genome_path.write_text(">chr1\nATGAAATTT\n", encoding="utf-8")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "direct",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Arabidopsis_thaliana_genes.derived.cds.fa.gz"
    formatted_genome = out_genome / "Arabidopsis_thaliana_softmasked.fa.gz"
    assert formatted_cds.exists()
    assert formatted_genome.exists()
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert ">Arabidopsis_thaliana_gene1" in text
    assert "ATGAAATTT" in text


def test_coge_export_reconstructs_complete_cds_from_fragment_ids(tmp_path):
    mod = load_module()
    genome = tmp_path / "genome.fa"
    genome.write_text(">lcl|chr1\nATGCCCTAA\n>lcl|chr2\nTTACCCCAT\n")
    gff = tmp_path / "export.gff"
    gff.write_text("\n".join([
        "chr1\tCoGe\tCDS\t1\t3\t.\t+\t.\tID=gene1;Name=gene1;CDS=gene1;coge_fid=101",
        "chr1\tCoGe\tCDS\t7\t9\t.\t+\t.\tID=gene1.CDS2;Name=gene1;CDS=gene1;coge_fid=101",
        "chr2\tCoGe\tCDS\t1\t3\t.\t-\t.\tID=gene2;Name=gene2;CDS=gene2;coge_fid=102",
        "chr2\tCoGe\tCDS\t7\t9\t.\t-\t.\tID=gene2.CDS2;Name=gene2;CDS=gene2;coge_fid=102",
        "",
    ]))
    records = list(mod.derive_cds_records_from_gff_and_genome({
        "provider": "coge", "species_key": "Pinguicula_agnata",
        "gff_path": gff, "genome_path": genome,
    }))
    assert records == [("gene1 [gene=gene1]", "ATGTAA"), ("gene2 [gene=gene2]", "ATGTAA")]


@pytest.mark.parametrize("second", [
    "chr2\tCoGe\tCDS\t7\t9\t.\t+\t.\tID=gene1.CDS2;Name=gene1;CDS=gene1;coge_fid=101",
    "chr1\tCoGe\tCDS\t7\t9\t.\t-\t.\tID=gene1.CDS2;Name=gene1;CDS=gene1;coge_fid=101",
    "chr1\tCoGe\tCDS\t7\t9\t.\t+\t.\tID=gene1.CDS2;Name=gene1;CDS=gene1;coge_fid=102",
    "chr1\tCoGe\tCDS\t7\t9\t.\t+\t.\tID=gene1.CDS2;Name=gene1;CDS=gene2;coge_fid=101",
])
def test_coge_export_rejects_conflicting_source_models(tmp_path, second):
    mod = load_module()
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\nATGCCCTAA\n>chr2\nATGCCCTAA\n")
    gff = tmp_path / "export.gff"
    gff.write_text(
        "chr1\tCoGe\tCDS\t1\t3\t.\t+\t.\tID=gene1;Name=gene1;CDS=gene1;coge_fid=101\n"
        + second + "\n"
    )
    with pytest.raises(ValueError, match="CoGe CDS feature identity"):
        list(mod.derive_cds_records_from_gff_and_genome({
            "provider": "coge", "species_key": "Pinguicula_agnata",
            "gff_path": gff, "genome_path": genome,
        }))


@pytest.mark.parametrize("provider,suffix", [("direct", ".t"), ("coge", ".mRNA")])
@pytest.mark.parametrize("mode", ["strict", "rescue_overlap"])
@pytest.mark.parametrize("provided_cds", [False, True])
def test_orphan_transcript_suffixes_select_one_longest_at_shared_coding_locus(
    tmp_path, provider, suffix, mode, provided_cds,
):
    mod = load_module()
    genome, gff = tmp_path / "genome.fa", tmp_path / "annotation.gff"
    genome.write_text(">chr1\nATGAAACCC\n")
    t1, t2 = "locusX" + suffix + "1", "locusX" + suffix + "2"
    if provider == "direct":
        rows = [
            f"chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID={t1}",
            f"chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent={t1}",
            f"chr1\tsrc\tmRNA\t1\t6\t.\t+\t.\tID={t2}",
            f"chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tID=cds2;Parent={t2}",
        ]
    else:
        rows = [
            f"chr1\tCoGe\tCDS\t1\t9\t.\t+\t.\tID={t1};Name={t1};coge_fid=1",
            f"chr1\tCoGe\tCDS\t1\t6\t.\t+\t.\tID={t2};Name={t2};coge_fid=2",
        ]
    gff.write_text("\n".join(rows) + "\n")
    task = dict(provider=provider, species_key="Test_species", species_prefix="Test_species",
                gff_path=gff, genome_path=genome, gene_grouping_mode=mode, gff_repair_mode="safe")
    if provided_cds:
        cds = tmp_path / "raw.fa"
        cds.write_text(f">{t1}\nATGAAACCC\n>{t2}\nATGAAA\n")
        task["cds_path"] = cds
    output = tmp_path / "cds"
    output.mkdir()
    result = mod.format_cds(task, output, False, False)
    with gzip.open(result["output_path"], "rt") as handle:
        assert handle.read() == ">Test_species_locusX\nATGAAACCC\n"
    gff_output = tmp_path / "gff"
    gff_output.mkdir()
    formatted_gff = mod.format_gff(task, gff_output, False, False, formatted_cds_path=result["output_path"])
    from validate_cds_gff_mapping import validate_single_species
    mapping = validate_single_species(dict(index=1, species_prefix="Test_species",
        cds_file=result["output_path"], gff_file=formatted_gff["output_path"], strict=True), 10)
    assert mapping["ok"], json.dumps(mapping)


@pytest.mark.parametrize("mode", ["strict", "rescue_overlap"])
@pytest.mark.parametrize("provided_cds", [False, True])
def test_orphan_suffixes_preserve_disjoint_or_opposite_strand_coding_loci(tmp_path, mode, provided_cds):
    mod = load_module()
    genome, gff = tmp_path / "genome.fa", tmp_path / "annotation.gff"
    genome.write_text(">chr1\nATGAAACCCGGGTTTAAACCC\n>chr2\nATGAAACCC\n")
    # Shared transcript span alone is insufficient; these CDSs are disjoint.
    models = [("locusX.t1", "chr1", "+", 1, 6), ("locusX.t2", "chr1", "+", 10, 15),
              ("locusX.t3", "chr1", "-", 1, 6), ("locusX.t4", "chr2", "+", 1, 6)]
    rows = []
    for tid, axis, strand, start, end in models:
        rows.extend([
            f"{axis}\tsrc\tmRNA\t1\t21\t.\t{strand}\t.\tID={tid}",
            f"{axis}\tsrc\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID=cds-{tid};Parent={tid}",
        ])
    gff.write_text("\n".join(rows) + "\n")
    task = dict(provider="direct", species_key="Test_species", species_prefix="Test_species",
                gff_path=gff, genome_path=genome, gene_grouping_mode=mode)
    if provided_cds:
        cds = tmp_path / "raw.fa"
        cds.write_text("".join(f">{tid}\nATGAAA\n" for tid, *_ in models))
        task["cds_path"] = cds
    output = tmp_path / "cds"
    output.mkdir()
    result = mod.format_cds(task, output, False, False)
    records = dict(mod.iter_fasta_records(result["output_path"]))
    assert set(records) == {"Test_species_" + tid for tid, *_ in models}


@pytest.mark.parametrize("mode", ["strict", "rescue_overlap"])
@pytest.mark.parametrize("provided_cds", [False, True])
@pytest.mark.parametrize("strand", ["+", "-"])
def test_missing_parent_siblings_keep_shared_owner_with_export_suffix(
    tmp_path, mode, provided_cds, strand,
):
    from Bio.Seq import Seq
    mod = load_module()
    genome, gff = tmp_path / "genome.fa", tmp_path / "annotation.gff"
    sequence = "ATG" * 9
    genome.write_text(">chr1\n" + sequence + "\n")
    tids = ["locusX.t1", "locusX.t1.1.5db15392"]
    parts = [[(1, 9), (19, 27)], [(1, 15)]]
    rows = []
    for tid, blocks in zip(tids, parts, strict=True):
        rows.append(f"chr1\tsrc\tmRNA\t1\t27\t.\t{strand}\t.\tID={tid};Parent=locusX;Name=locusX.t1")
        rows.extend(f"chr1\tsrc\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID=cds.{tid};Parent={tid}"
                    for start, end in blocks)
    gff.write_text("\n".join(rows) + "\n")
    task = dict(provider="direct", species_key="Test_species", species_prefix="Test_species",
                gff_path=gff, genome_path=genome, gene_grouping_mode=mode, gff_repair_mode="safe")
    coding = ["".join(sequence[start - 1:end] for start, end in blocks) for blocks in parts]
    if strand == "-":
        coding = [str(Seq(value).reverse_complement()) for value in coding]
    if provided_cds:
        cds = tmp_path / "raw.fa"
        cds.write_text("".join(f">{tid}\n{value}\n" for tid, value in zip(tids, coding, strict=True)))
        task["cds_path"] = cds
    output = tmp_path / "output"
    output.mkdir()
    result = mod.format_cds(task, output, False, False, strict=True)
    assert result["before_count"] == 2 and result["after_count"] == 1
    assert dict(mod.iter_fasta_records(result["output_path"])) == {"Test_species_locusX": coding[0]}
    formatted = mod.format_gff(task, output, False, False, formatted_cds_path=result["output_path"])
    from validate_cds_gff_mapping import validate_single_species
    mapping = validate_single_species(dict(index=1, species_prefix="Test_species",
        cds_file=result["output_path"], gff_file=formatted["output_path"], strict=True), 10)
    assert mapping["ok"], json.dumps(mapping)
    from gff2genestat import process_single_gff
    traits = process_single_gff(Path(formatted["output_path"]).name, str(output),
        ["Test_species_locusX"], "CDS", "longest",
        ["sequence", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"],
        ["gene_id", "feature_size"])
    assert traits.gene_id.tolist() == ["Test_species_locusX"]
    assert traits.feature_size.tolist() == [18]


@pytest.mark.parametrize("mode", ["strict", "rescue_overlap"])
def test_missing_parent_export_siblings_preserve_declared_owner(tmp_path, mode):
    mod = load_module()
    gff = tmp_path / "annotation.gff"
    gff.write_text(
        "chr1\tsrc\tmRNA\t1\t6\t.\t+\t.\tID=locusX.t1;Parent=locusX\n"
        "chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tID=c1;Parent=locusX.t1\n"
        "chr1\tsrc\tmRNA\t20\t25\t.\t+\t.\tID=locusX.t1.1.5db15392;Parent=locusX\n"
        "chr1\tsrc\tCDS\t20\t25\t.\t+\t0\tID=c2;Parent=locusX.t1.1.5db15392\n"
    )
    index = mod.build_gff_cds_grouping_index(dict(provider="direct", species_key="Test_species",
        gff_path=gff, gene_grouping_mode=mode))
    assert set(index["transcript_gene_tokens"].values()) == {"locusX"}


@pytest.mark.parametrize("mode", ["strict", "rescue_overlap"])
def test_overlapping_missing_parent_genes_remain_distinct(tmp_path, mode):
    mod = load_module()
    gff = tmp_path / "annotation.gff"
    gff.write_text(
        "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=gene%3Ag1\n"
        "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=t1\n"
        "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=t2;Parent=g2\n"
        "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=c2;Parent=t2\n")
    index = mod.build_gff_cds_grouping_index(dict(provider="direct", species_key="Test_species",
        gff_path=gff, gene_grouping_mode=mode))
    assert index["transcript_gene_tokens"] == {"t1": "g1", "t2": "g2"}


@pytest.mark.parametrize("provided_cds", [False, True])
@pytest.mark.parametrize("prefix,encoded", [("", False), ("evm.model.", False), ("evm_27.model.", False), ("evm.model.", True)])
def test_coge_overlap_rescue_uses_complete_export_models(tmp_path, provided_cds, prefix, encoded):
    mod = load_module()
    genome, gff = tmp_path / "genome.fa", tmp_path / "annotation.gff"
    genome.write_text(">chr1\n" + "ATG" * 9 + "\n")
    gff.write_text(
        "chr1\tCoGe\tCDS\t1\t6\t.\t+\t0\tID=modelA;Name=modelA;coge_fid=101\n"
        "chr1\tCoGe\tCDS\t19\t27\t.\t+\t0\tID=modelA.CDS2;Name=modelA;coge_fid=101\n"
        "chr1\tCoGe\tCDS\t19\t27\t.\t+\t0\tID=modelB;Name=modelB;coge_fid=102\n"
    )
    gff.write_text(gff.read_text().replace("modelA", prefix + "modelA").replace("modelB", prefix + "modelB"))
    if encoded:
        gff.write_text(gff.read_text().replace(prefix, prefix.replace('.', '%2E')))
    task = dict(provider="coge", species_key="Test_species", species_prefix="Test_species",
                gff_path=gff, genome_path=genome, gene_grouping_mode="rescue_overlap", gff_repair_mode="safe")
    if provided_cds:
        cds = tmp_path / "raw.fa"
        cds.write_text(">" + prefix + "modelA\n" + "ATG" * 5 + "\n>" + prefix + "modelB\n" + "ATG" * 3 + "\n")
        task["cds_path"] = cds
    index = mod.build_gff_cds_grouping_index(task)
    assert index["transcripts_total"] == 2
    assert index["transcript_gene_tokens"] == {prefix + "modelA": prefix + "modelA", prefix + "modelB": prefix + "modelA"}
    output = tmp_path / "output"
    output.mkdir()
    result = mod.format_cds(task, output, False, False, strict=True)
    assert dict(mod.iter_fasta_records(result["output_path"])) == {"Test_species_modelA": "ATG" * 5}
    formatted = mod.format_gff(task, output, False, False, formatted_cds_path=result["output_path"])
    from validate_cds_gff_mapping import validate_single_species
    mapping = validate_single_species(dict(index=1, species_prefix="Test_species",
        cds_file=result["output_path"], gff_file=formatted["output_path"], strict=True), 10)
    assert mapping["ok"], json.dumps(mapping)
    from gff2genestat import process_single_gff
    traits = process_single_gff(Path(formatted["output_path"]).name, str(output),
        ["Test_species_modelA"], "CDS", "longest",
        ["sequence", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"],
        ["gene_id", "feature_size"])
    assert traits.feature_size.tolist() == [15]


@pytest.mark.parametrize("mode", ["strict", "rescue_overlap"])
@pytest.mark.parametrize("provided_cds", [False, True])
@pytest.mark.parametrize("strand", ["+", "-"])
def test_disconnected_missing_parent_numeric_models_project_complete_gene_to_gff(
    tmp_path, mode, provided_cds, strand,
):
    from Bio.Seq import Seq
    from gff2genestat import process_single_gff
    mod = load_module()
    genome, gff = tmp_path / "genome.fa", tmp_path / "annotation.gff"
    sequence = "ATG" * 15
    genome.write_text(">chr1\n" + sequence + "\n")
    # The source explicitly assigns all transcripts to the same missing gene.
    # Its ownership takes precedence over suffix or coding-locus inference.
    root = "Calam.05G208800"
    models = [(root + ".1", 1, 6), (root + ".6", 1, 12),
              (root + ".7", 31, 39)]
    rows, coding = [], {}
    for tid, start, end in models:
        rows.extend([
            f"chr1\tsrc\tmRNA\t{start}\t{end}\t.\t{strand}\t.\tID={tid};Parent={root}",
            f"chr1\tsrc\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID=cds.{tid};Parent={tid}",
        ])
        coding[tid] = sequence[start - 1:end]
        if strand == "-":
            coding[tid] = str(Seq(coding[tid]).reverse_complement())
    gff.write_text("\n".join(rows) + "\n")
    task = dict(provider="direct", species_key="Test_species", species_prefix="Test_species",
                gff_path=gff, genome_path=genome, gene_grouping_mode=mode, gff_repair_mode="safe")
    if provided_cds:
        cds = tmp_path / "raw.fa"
        cds.write_text("".join(f">{tid}\n{value}\n" for tid, value in coding.items()))
        task["cds_path"] = cds
    output = tmp_path / "output"
    output.mkdir()
    result = mod.format_cds(task, output, False, False, strict=True)
    ids = ["Test_species_" + root]
    assert dict(mod.iter_fasta_records(result["output_path"])) == {
        ids[0]: coding[root + ".6"],
    }
    formatted = mod.format_gff(task, output, False, False, formatted_cds_path=result["output_path"])
    with gzip.open(formatted["output_path"], "rt") as handle:
        emitted = [line.rstrip().split("\t") for line in handle]
    assert [parts[:8] for parts in emitted] == [line.split("\t")[:8] for line in rows]
    assert all(parts[8].startswith(source.split("\t")[8] + ";gene_id=")
               for parts, source in zip(emitted, rows, strict=True))
    traits = process_single_gff(Path(formatted["output_path"]).name, str(output), ids,
        "CDS", "longest",
        ["sequence", "source", "feature", "start", "end", "score", "strand", "phase", "attributes"],
        ["gene_id", "feature_size"])
    assert dict(zip(traits.gene_id, traits.feature_size, strict=True)) == {ids[0]: 12}


def test_reverse_complement_preserves_complete_iupac_dna_alphabet():
    from Bio.Seq import Seq
    from format_species_annotation.common import reverse_complement
    sequence = "ACGTRYSWKMBDHVNacgtryswkmbdhvn"
    assert reverse_complement(sequence) == str(Seq(sequence).reverse_complement())
    assert reverse_complement(reverse_complement(sequence)) == sequence


def test_minus_strand_multiexon_iupac_cds_matches_independent_genome_reconstruction(tmp_path):
    from Bio.Seq import Seq
    mod = load_module()
    sequence = "ACGTRYSWKMBD" + "NNNNNN" + "HVNRYSWKMBDH"
    genome, gff = tmp_path / "genome.fa", tmp_path / "source.gff"
    genome.write_text(">chr1\n" + sequence + "\n")
    gff.write_text(
        "chr1\ts\tgene\t1\t30\t.\t-\t.\tID=g1\n"
        "chr1\ts\tmRNA\t1\t30\t.\t-\t.\tID=t1;Parent=g1\n"
        "chr1\ts\tCDS\t1\t12\t.\t-\t0\tID=c1;Parent=t1\n"
        "chr1\ts\tCDS\t19\t30\t.\t-\t0\tID=c2;Parent=t1\n")
    task = dict(provider="direct", species_key="Test_species", species_prefix="Test_species",
                genome_path=genome, gff_path=gff, gene_grouping_mode="rescue_overlap")
    expected = str(Seq(sequence[:12] + sequence[18:30]).reverse_complement())
    out = tmp_path / "output"
    out.mkdir()
    result = mod.format_cds(task, out, False, False, strict=True)
    assert dict(mod.iter_fasta_records(result["output_path"])) == {"Test_species_g1": expected}


@pytest.mark.parametrize("axis", [("chr2", "+"), ("chr1", "-")])
def test_formatter_missing_parent_rejects_conflicting_axes(tmp_path, axis):
    mod = load_module()
    gff = tmp_path / "source.gff"
    gff.write_text("chr1\ts\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1\n"
                   f"{axis[0]}\ts\tmRNA\t1\t9\t.\t{axis[1]}\t.\tID=t2;Parent=g1\n")
    task = dict(provider="direct", species_key="Test_species", species_prefix="Test_species",
                gff_path=gff, gene_grouping_mode="rescue_overlap")
    with pytest.raises(ValueError, match="Conflicting axes"):
        mod.build_gff_cds_grouping_index(task)


def test_suffix_grouping_keeps_explicit_gene_ids_and_connected_components(tmp_path):
    mod = load_module()
    genome, gff = tmp_path / "genome.fa", tmp_path / "annotation.gff"
    genome.write_text(">chr1\n" + "ATG" * 12 + "\n")
    gff.write_text(
        "chr1\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=X.t1;gene_id=first\n"
        "chr1\tsrc\tCDS\t1\t12\t.\t+\t0\tID=c1;Parent=X.t1\n"
        "chr1\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=X.t2;gene_id=second\n"
        "chr1\tsrc\tCDS\t1\t12\t.\t+\t0\tID=c2;Parent=X.t2\n"
        "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=Y.t1\n"
        "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=c3;Parent=Y.t1\n"
        "chr1\tsrc\tmRNA\t7\t15\t.\t+\t.\tID=Y.t2\n"
        "chr1\tsrc\tCDS\t7\t15\t.\t+\t0\tID=c4;Parent=Y.t2\n"
        "chr1\tsrc\tmRNA\t13\t21\t.\t+\t.\tID=Y.t3\n"
        "chr1\tsrc\tCDS\t13\t21\t.\t+\t0\tID=c5;Parent=Y.t3\n"
    )
    task = dict(provider="direct", species_key="Test_species", species_prefix="Test_species",
                gff_path=gff, genome_path=genome, gene_grouping_mode="rescue_overlap")
    records = dict(mod.derive_cds_records_from_gff_and_genome(task))
    assert set(records) == {"X.t1 [gene=first]", "X.t2 [gene=second]",
                            "Y.t1 [gene=Y]", "Y.t2 [gene=Y]", "Y.t3 [gene=Y]"}


def test_orphan_suffix_does_not_collide_with_explicit_gene_or_late_parent(tmp_path):
    mod = load_module()
    genome, gff = tmp_path / "genome.fa", tmp_path / "annotation.gff"
    genome.write_text(">chr1\n" + "ATG" * 10 + "\n")
    gff.write_text(
        "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=late.t1;Parent=actual_gene\n"
        "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=late.t1\n"
        "chr1\tsrc\tmRNA\t1\t6\t.\t+\t.\tID=X.t1\n"
        "chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tID=c2;Parent=X.t1\n"
        "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=X.t2\n"
        "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=c3;Parent=X.t2\n"
        "chr1\tsrc\tCDS\t20\t28\t.\t+\t0\tID=explicit;gene_id=X\n"
        "chr1\tsrc\tCDS\t10\t18\t.\t+\t0\tID=Z.t1\n"
        "chr1\tsrc\tCDS\t22\t30\t.\t+\t0\tID=explicit_z;gene_id=Z\n"
        "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=actual_gene\n"
    )
    task = dict(provider="direct", species_key="Test_species", gff_path=gff, genome_path=genome)
    records = dict(mod.derive_cds_records_from_gff_and_genome(task))
    assert set(records) == {"late.t1 [gene=actual_gene]", "explicit [gene=X]",
                            "X.t1 [gene=X.t1]", "X.t2 [gene=X.t1]",
                            "Z.t1 [gene=Z.t1]", "explicit_z [gene=Z]"}


def test_non_coge_parentless_cds_preserves_distinct_feature_ids(tmp_path):
    mod = load_module()
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\nATGCCCTAA\n")
    gff = tmp_path / "annotation.gff"
    gff.write_text(
        "chr1\tCoGe\tCDS\t1\t3\t.\t+\t.\tID=gene1;Name=gene1;coge_fid=101\n"
        "chr1\tCoGe\tCDS\t7\t9\t.\t+\t.\tID=gene1.CDS2;Name=gene1;coge_fid=101\n"
    )
    records = list(mod.derive_cds_records_from_gff_and_genome({
        "provider": "direct", "species_key": "Pinguicula_agnata",
        "gff_path": gff, "genome_path": genome,
    }))
    assert len(records) == 2
    assert [sequence for _, sequence in records] == ["ATG", "TAA"]


def test_coge_canonical_source_names_match_formatted_cds_and_gff(tmp_path):
    mod = load_module()
    from validate_cds_gff_mapping import validate_single_species

    genome, gff = tmp_path / "genome.fa", tmp_path / "export.gff"
    genome.write_text(">lcl|chr1\nATGCCCTAA\n")
    gff.write_text(
        "chr1\tCoGe\tCDS\t1\t3\t.\t+\t.\tID=model__1;Name=model__1;CDS=model__1;coge_fid=101\n"
        "chr1\tCoGe\tCDS\t7\t9\t.\t+\t.\tID=model__1.CDS2;Name=model__1;CDS=model__1;coge_fid=101\n"
    )
    task = dict(provider="coge", species_key="Pinguicula_agnata", species_prefix="Pinguicula_agnata",
                gff_path=gff, genome_path=genome, gene_grouping_mode="rescue_overlap", gff_repair_mode="safe")
    cds_dir, gff_dir = tmp_path / "cds", tmp_path / "gff"
    cds_dir.mkdir()
    gff_dir.mkdir()
    cds = mod.format_cds(task, cds_dir, False, False)
    result = mod.format_gff(task, gff_dir, False, False, formatted_cds_path=cds["output_path"])
    with gzip.open(cds["output_path"], "rt") as handle:
        assert handle.read() == ">Pinguicula_agnata_model_1\nATGTAA\n"
    with gzip.open(result["output_path"], "rt") as handle:
        text = handle.read()
    assert "ID=model__1.CDS2;Name=model_1;CDS=model_1;coge_fid=101" in text
    assert text.count("\tCDS\t") == 2
    mapping = validate_single_species(dict(index=1, species_prefix="Pinguicula_agnata",
                                          cds_file=cds["output_path"], gff_file=result["output_path"], strict=True), 10)
    assert mapping["ok"], mapping
    audit = json.loads(Path(str(result["output_path"]) + ".repair.json").read_text())
    assert audit["status"] == "repaired" and audit["changed_values"] == 4


def test_coge_source_names_reject_canonical_collisions(tmp_path):
    mod = load_module()
    genome, gff = tmp_path / "genome.fa", tmp_path / "export.gff"
    genome.write_text(">chr1\nATGCCCTAA\n")
    gff.write_text(
        "chr1\tCoGe\tCDS\t1\t3\t.\t+\t.\tID=model__1;Name=model__1;coge_fid=101\n"
        "chr1\tCoGe\tCDS\t7\t9\t.\t+\t.\tID=model_1;Name=model_1;coge_fid=102\n"
    )
    task = dict(provider="coge", species_key="Pinguicula_agnata", gff_path=gff, genome_path=genome)
    with pytest.raises(ValueError, match="Conflicting CoGe CDS feature identity"):
        list(mod.derive_cds_records_from_gff_and_genome(task))


def test_format_species_inputs_derives_cds_without_trimming_nonzero_phase(tmp_path):
    input_dir = tmp_path / "Direct" / "species_wise_original"
    species_dir = input_dir / "Arabidopsis_thaliana"
    species_dir.mkdir(parents=True, exist_ok=True)
    gff_path = species_dir / "Arabidopsis_thaliana.annotation.gff3"
    genome_path = species_dir / "Arabidopsis_thaliana.genome.fa"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t15\t.\t+\t.\tID=gene1",
                "chr1\tsrc\tmRNA\t1\t15\t.\t+\t.\tID=gene1.t1;Parent=gene1",
                "chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tID=cds1;Parent=gene1.t1",
                "chr1\tsrc\tCDS\t10\t15\t.\t+\t2\tID=cds2;Parent=gene1.t1",
                "chr2\tsrc\tgene\t1\t15\t.\t-\t.\tID=gene2",
                "chr2\tsrc\tmRNA\t1\t15\t.\t-\t.\tID=gene2.t1;Parent=gene2",
                "chr2\tsrc\tCDS\t1\t6\t.\t-\t1\tID=cds3;Parent=gene2.t1",
                "chr2\tsrc\tCDS\t10\t15\t.\t-\t0\tID=cds4;Parent=gene2.t1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    genome_path.write_text(
        ">chr1\nATGAAACCCGGGTTT\n>chr2\nAAACCCGGGTTTCAT\n",
        encoding="utf-8",
    )

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "direct",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Arabidopsis_thaliana_annotation.derived.cds.fa.gz"
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert ">Arabidopsis_thaliana_gene1" in text
    assert ">Arabidopsis_thaliana_gene2" in text
    assert "ATGAAAGGGTTT" in text
    assert "ATGAAATTT" not in text


def test_build_formatted_cds_id_ensembl_uses_transcript_token_for_tie_breaks():
    module = load_module()
    task = {
        "provider": "ensemblplants",
        "species_prefix": "Arabidopsis_thaliana",
    }
    raw_header = "AT1G01520.5 cds chromosome:TAIR10:1:190478:192436:1 gene:AT1G01520 gene_symbol:ASG4"
    derived_header = "transcript:AT1G01520.1 gene=AT1G01520"

    assert module.build_formatted_cds_id(task, raw_header) == "Arabidopsis_thaliana_AT1G01520.5"
    assert module.build_formatted_cds_id(task, derived_header) == "Arabidopsis_thaliana_AT1G01520.1"
    assert (
        module.build_gene_aggregate_id(task, raw_header, "Arabidopsis_thaliana_AT1G01520.5")
        == "Arabidopsis_thaliana_AT1G01520"
    )
    assert (
        module.build_gene_aggregate_id(task, derived_header, "Arabidopsis_thaliana_AT1G01520.1")
        == "Arabidopsis_thaliana_AT1G01520"
    )


def test_format_species_inputs_uses_gff_hierarchy_for_provided_cds_longest_selection(tmp_path):
    input_dir = tmp_path / "Direct" / "species_wise_original"
    species_dir = input_dir / "Arabidopsis_thaliana"
    species_dir.mkdir(parents=True, exist_ok=True)
    cds_path = species_dir / "models.cds.fa"
    gff_path = species_dir / "models.gff3"
    species_summary = tmp_path / "gg_input_generation_species.tsv"
    cds_path.write_text(
        ">transcript_alpha\nATGAAATTT\n>transcript_beta\nATGAAACCCGGGTTT\n",
        encoding="utf-8",
    )
    gff_path.write_text(
        "\n".join(
            [
                "##gff-version 3",
                "chr1\tsrc\tgene\t1\t30\t.\t+\t.\tID=gene_from_gff",
                "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=transcript_alpha;Parent=gene_from_gff;longest=1",
                "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=cds_alpha;Parent=transcript_alpha",
                "chr1\tsrc\tmRNA\t16\t30\t.\t+\t.\tID=transcript_beta;Parent=gene_from_gff;longest=0",
                "chr1\tsrc\tCDS\t16\t30\t.\t+\t0\tID=cds_beta;Parent=transcript_beta",
                "",
            ]
        ),
        encoding="utf-8",
    )

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "direct",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
        "--species-summary-output",
        str(species_summary),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = next(out_cds.glob("*.fa.gz"))
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert text == ">Arabidopsis_thaliana_gene_from_gff\nATGAAACCCGGGTTT\n"

    with open(species_summary, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["cds_grouping_source"] == "gff"
    assert rows[0]["cds_gff_records_mapped"] == "2"
    assert rows[0]["cds_gff_records_unmapped"] == "0"
    audit_path = Path(rows[0]["cds_gff_grouping_audit_path"])
    assert audit_path.exists()
    with open(audit_path, "rt", encoding="utf-8", newline="") as handle:
        audit_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["mapping_status"] for row in audit_rows] == ["mapped", "mapped"]
    assert [row["selected_longest"] for row in audit_rows] == ["0", "1"]

    original_stat = formatted_cds.stat()
    corrupted_bytes = bytearray(formatted_cds.read_bytes())
    corrupted_bytes[len(corrupted_bytes) // 2] ^= 1
    formatted_cds.write_bytes(corrupted_bytes)
    os.utime(formatted_cds, ns=(original_stat.st_atime_ns, original_stat.st_mtime_ns))
    repaired = run_script(
        "--provider",
        "direct",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
        "--species-summary-output",
        str(species_summary),
    )
    assert repaired.returncode == 0, repaired.stderr + "\n" + repaired.stdout
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        assert handle.read() == ">Arabidopsis_thaliana_gene_from_gff\nATGAAACCCGGGTTT\n"

    cds_input_stat = cds_path.stat()
    cds_path.write_text(
        ">transcript_alpha\nATGAAATTT\n>transcript_beta\nATGCCCAAAGGGTTT\n",
        encoding="utf-8",
    )
    os.utime(cds_path, ns=(cds_input_stat.st_atime_ns, cds_input_stat.st_mtime_ns))
    changed_cds = run_script(
        "--provider",
        "direct",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
        "--species-summary-output",
        str(species_summary),
    )
    assert changed_cds.returncode == 0, changed_cds.stderr + "\n" + changed_cds.stdout
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        assert handle.read() == ">Arabidopsis_thaliana_gene_from_gff\nATGCCCAAAGGGTTT\n"

    gff_input_stat = gff_path.stat()
    changed_gff_text = gff_path.read_text(encoding="utf-8").replace("gene_from_gff", "gene_from_xff")
    assert len(changed_gff_text.encode("utf-8")) == gff_input_stat.st_size
    gff_path.write_text(changed_gff_text, encoding="utf-8")
    os.utime(gff_path, ns=(gff_input_stat.st_atime_ns, gff_input_stat.st_mtime_ns))
    next(out_gff.glob("*.gff.gz")).unlink()
    changed_gff = run_script(
        "--provider",
        "direct",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
        "--species-summary-output",
        str(species_summary),
    )
    assert changed_gff.returncode == 0, changed_gff.stderr + "\n" + changed_gff.stdout
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        assert handle.read() == ">Arabidopsis_thaliana_gene_from_xff\nATGCCCAAAGGGTTT\n"
    with open(str(formatted_cds) + ".gff-grouping.json", "rt", encoding="utf-8") as handle:
        audit = json.load(handle)
    assert audit["version"] == 15
    assert len(audit["cds_input"]["sha256"]) == 64
    assert len(audit["gff_input"]["sha256"]) == 64

    longest = run_validate_longest_script(
        "--species-cds-dir",
        str(out_cds),
        "--species-summary",
        str(species_summary),
    )
    assert longest.returncode == 0, longest.stderr + "\n" + longest.stdout


@pytest.mark.parametrize(
    ("species_prefix", "gene_id", "isoforms", "selected_protein"),
    (
        (
            "Erythranthe_tilingii",
            "ACP275_01G000100",
            (
                ("KAL7121685.1", "rna-Mitil.01G000100.2", 9),
                ("KAL7121688.1", "rna-Mitil.01G000100.3", 9),
                ("KAL7121687.1", "rna-Mitil.01G000100.1", 10),
                ("KAL7121686.1", "rna-Mitil.01G000100.4", 6),
            ),
            "KAL7121687.1",
        ),
        (
            "Phlomoides_rotata",
            "ACS0TY_000073",
            (
                ("KAL8550850.1", "rna-Phrot01aG0007500.1", 7),
                ("KAL8550851.1", "rna-Phrot01aG0007500.2", 5),
                ("KAL8550852.1", "rna-Phrot01aG0007500.3", 5),
                ("KAL8550853.1", "rna-Phrot01aG0007500.4", 5),
            ),
            "KAL8550850.1",
        ),
        (
            "Saponaria_officinalis",
            "RND81_02G052900",
            tuple(
                ("KAK97483{}.1".format(suffix), "rna-Sapof.02G052900.{}".format(index), 8)
                for index, suffix in enumerate(range(69, 75), 1)
            ),
            "KAK9748369.1",
        ),
        (
            "Rhododendron_molle",
            "RHMOL_Rhmol01G0002200",
            (
                ("KAI8570035.1", "rna-Rhmol01G0002200.1", 4),
                ("KAI8570036.1", "rna-Rhmol01G0002200.2", 4),
                ("KAI8570037.1", "rna-Rhmol01G0002200.3", 3),
                ("KAI8570038.1", "rna-Rhmol01G0002200.4", 5),
            ),
            "KAI8570038.1",
        ),
    ),
)
def test_issue_26_ncbi_isoforms_group_by_gff_gene(
    tmp_path,
    species_prefix,
    gene_id,
    isoforms,
    selected_protein,
):
    module = load_module()
    cds_path = tmp_path / "models_cds_from_genomic.fna.gz"
    gff_path = tmp_path / "models_genomic.gff.gz"
    with gzip.open(cds_path, "wt", encoding="utf-8") as handle:
        for index, (protein_id, _transcript_id, codon_length) in enumerate(isoforms, 1):
            handle.write(
                ">lcl|CM000001.1_cds_{}_{} [protein_id={}]\n{}\n".format(
                    protein_id,
                    index,
                    protein_id,
                    "ATG" * codon_length,
                )
            )
    with gzip.open(gff_path, "wt", encoding="utf-8") as handle:
        handle.write("chr1\tsrc\tgene\t1\t1000\t.\t+\t.\tID=gene-{};locus_tag={}\n".format(gene_id, gene_id))
        for index, (protein_id, transcript_id, codon_length) in enumerate(isoforms, 1):
            handle.write(
                "chr1\tsrc\tmRNA\t{}\t{}\t.\t+\t.\tID={};Parent=gene-{}\n".format(
                    index,
                    index + codon_length * 3 - 1,
                    transcript_id,
                    gene_id,
                )
            )
            handle.write(
                "chr1\tsrc\tCDS\t{}\t{}\t.\t+\t0\tID=cds-{};Parent={};protein_id={};locus_tag={}\n".format(
                    index,
                    index + codon_length * 3 - 1,
                    protein_id,
                    transcript_id,
                    protein_id,
                    gene_id,
                )
            )

    task = {
        "provider": "ncbi",
        "species_key": species_prefix,
        "species_prefix": species_prefix,
        "cds_path": cds_path,
        "gff_path": gff_path,
        "gene_grouping_mode": "strict",
        "format_strict": True,
    }
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    result = module.format_cds(task, output_dir, overwrite=False, dry_run=False)

    assert result["before_count"] == len(isoforms)
    assert result["after_count"] == 1
    assert result["gff_records_mapped"] == len(isoforms)
    assert result["gff_records_unmapped"] == 0
    with gzip.open(result["output_path"], "rt", encoding="utf-8") as handle:
        assert handle.readline().strip() == ">{}_{}".format(species_prefix, gene_id)
    with open(result["gff_grouping_audit_path"], "rt", encoding="utf-8", newline="") as handle:
        audit_rows = list(csv.DictReader(handle, delimiter="\t"))
    selected_rows = [row for row in audit_rows if row["selected_longest"] == "1"]
    assert len(selected_rows) == 1
    assert selected_protein in selected_rows[0]["raw_cds_id"]


def test_issue_26_figshare_prefers_gene_gff_and_groups_isoforms(tmp_path):
    module = load_module()
    input_dir = tmp_path / "Figshare" / "species_wise_original"
    species_dir = input_dir / "Euryodendron_excelsum"
    species_dir.mkdir(parents=True)
    (species_dir / "Euryodendron_excelsum.cds.fa").write_text(
        "".join(
            ">FUN_001415-T{}\n{}\n".format(index, "ATG" * codon_length)
            for index, codon_length in enumerate((9, 8, 7, 6, 5, 4, 3), 1)
        ),
        encoding="utf-8",
    )
    canonical_gff = species_dir / "Euryodendron_excelsum.gff3"
    canonical_gff.write_text(
        "chr1\tfunannotate\tgene\t1\t100\t.\t+\t.\tID=FUN_001415;\n"
        + "".join(
            "chr1\tfunannotate\tmRNA\t1\t100\t.\t+\t.\tID=FUN_001415-T{0};Parent=FUN_001415;\n"
            "chr1\tfunannotate\tCDS\t{0}\t{1}\t.\t+\t0\tID=FUN_001415-T{0}.cds;Parent=FUN_001415-T{0};\n".format(
                index,
                index + codon_length * 3 - 1,
            )
            for index, codon_length in enumerate((9, 8, 7, 6, 5, 4, 3), 1)
        ),
        encoding="utf-8",
    )
    (species_dir / "Euryodendron_excelsum.EDTA.gff3").write_text(
        "chr1\tEDTA\trepeat_region\t1\t100\t.\t+\t.\tID=repeat1\n",
        encoding="utf-8",
    )
    (species_dir / "Euryodendron_excelsum.only_long-transcripts.gff3").write_text(
        "chr1\tfunannotate\tgene\t1\t100\t.\t+\t.\tID=FUN_001415;\n"
        "chr1\tfunannotate\tmRNA\t1\t100\t.\t+\t.\tID=FUN_001415-T1;Parent=FUN_001415;\n"
        "chr1\tfunannotate\tCDS\t1\t27\t.\t+\t0\tID=FUN_001415-T1.cds;Parent=FUN_001415-T1;\n",
        encoding="utf-8",
    )

    tasks, warnings, errors = module.discover_tasks("figshare", input_dir)
    assert errors == []
    assert tasks[0]["gff_path"] == canonical_gff
    assert tasks[0]["gff_auto_selected_from_multiple"] is True
    assert len(tasks[0]["gff_selection_candidates"]) == 3
    assert any("Using 'Euryodendron_excelsum.gff3'" in warning for warning in warnings)

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    summary = tmp_path / "summary.tsv"
    completed = run_script(
        "--provider",
        "figshare",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
        "--species-summary-output",
        str(summary),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    formatted_cds = next(out_cds.glob("*.fa.gz"))
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        assert handle.read() == ">Euryodendron_excelsum_FUN_001415\n{}\n".format("ATG" * 9)
    with open(summary, "rt", encoding="utf-8", newline="") as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["cds_gff_records_mapped"] == "7"
    assert row["cds_gff_records_unmapped"] == "0"


def test_figshare_rejects_all_unmapped_auto_selected_gff(tmp_path):
    module = load_module()
    cds_path = tmp_path / "Euryodendron_excelsum.cds.fa"
    gff_path = tmp_path / "Euryodendron_excelsum.EDTA.gff3"
    cds_path.write_text(">FUN_001415-T1\nATG\n>FUN_001415-T2\nATGAAA\n", encoding="utf-8")
    gff_path.write_text(
        "chr1\tEDTA\trepeat_region\t1\t100\t.\t+\t.\tID=repeat1\n",
        encoding="utf-8",
    )
    task = {
        "provider": "figshare",
        "species_key": "Euryodendron_excelsum",
        "species_prefix": "Euryodendron_excelsum",
        "cds_path": cds_path,
        "gff_path": gff_path,
        "gene_grouping_mode": "strict",
        "gff_auto_selected_from_multiple": True,
        "gff_selection_candidates": (
            "Euryodendron_excelsum.EDTA.gff3",
            "Euryodendron_excelsum.gff3",
        ),
    }
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    with pytest.raises(ValueError, match=r"unexpected_unmapped=2 ambiguous=0"):
        module.format_cds(task, output_dir, overwrite=False, dry_run=False)
    audit_json, audit_tsv = module.cds_gff_grouping_audit_paths(
        output_dir / module.normalize_cds_output_basename(cds_path.name, task["species_prefix"])
    )
    assert not list(output_dir.glob("*.fa.gz"))
    assert json.loads(audit_json.read_text(encoding="utf-8"))["status"] == "failed"
    assert len(audit_tsv.read_text(encoding="utf-8").splitlines()) == 3


def test_empty_raw_cds_does_not_create_formatted_output(tmp_path):
    module = load_module()
    cds_path = tmp_path / "empty.cds.fa"
    cds_path.write_text("", encoding="utf-8")
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    task = {
        "provider": "direct",
        "species_key": "Example_species",
        "species_prefix": "Example_species",
        "cds_path": cds_path,
        "gff_path": None,
    }

    result = module.format_cds(task, output_dir, overwrite=False, dry_run=False)

    assert result["status"] == "empty"
    assert result["output_path"] is None
    assert not list(output_dir.iterdir())


def test_provided_cds_gff_grouping_rejects_unmapped_records_by_default(tmp_path):
    input_dir = tmp_path / "Direct" / "species_wise_original"
    species_dir = input_dir / "Arabidopsis_thaliana"
    species_dir.mkdir(parents=True, exist_ok=True)
    (species_dir / "models.cds.fa").write_text(
        ">known_transcript\nATGAAATTT\n>missing_transcript\nATGCCCTTT\n",
        encoding="utf-8",
    )
    (species_dir / "models.gff3").write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=known_gene",
                "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=known_transcript;Parent=known_gene",
                "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tParent=known_transcript",
                "",
            ]
        ),
        encoding="utf-8",
    )

    normal_root = tmp_path / "normal"
    normal_summary = normal_root / "summary.tsv"
    normal = run_script(
        "--provider",
        "direct",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(normal_root / "cds"),
        "--species-gff-dir",
        str(normal_root / "gff"),
        "--species-genome-dir",
        str(normal_root / "genome"),
        "--species-summary-output",
        str(normal_summary),
    )
    assert normal.returncode != 0
    assert "unexpected_unmapped=1 ambiguous=0" in normal.stderr

    strict_root = tmp_path / "strict"
    strict = run_script(
        "--provider",
        "direct",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(strict_root / "cds"),
        "--species-gff-dir",
        str(strict_root / "gff"),
        "--species-genome-dir",
        str(strict_root / "genome"),
        "--strict",
    )
    assert strict.returncode != 0
    assert "unexpected_unmapped=1 ambiguous=0" in strict.stderr


def test_provided_cds_gff_grouping_tolerates_only_low_rate_residual_mismatch(tmp_path):
    module = load_module()
    cds_path = tmp_path / "models.cds.fa"
    gff_path = tmp_path / "models.gff3"
    total_records = 1000
    mapped_records = total_records - 1
    cds_path.write_text(
        "".join(">tx{}\nATG\n".format(index) for index in range(total_records)),
        encoding="utf-8",
    )
    gff_path.write_text(
        "".join(
            "chr1\tsrc\tCDS\t{}\t{}\t.\t+\t0\tID=tx{};gene_id=gene{}\n".format(
                index * 3 + 1,
                index * 3 + 3,
                index,
                index,
            )
            for index in range(mapped_records)
        ),
        encoding="utf-8",
    )
    task = {
        "provider": "direct",
        "species_key": "Adiantum_capillus-veneris",
        "species_prefix": "Adiantum_capillus-veneris",
        "cds_path": cds_path,
        "gff_path": gff_path,
        "gene_grouping_mode": "strict",
        "format_strict": False,
    }
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    result = module.format_cds(task, output_dir, overwrite=False, dry_run=False)

    assert result["status"] == "write"
    assert result["gff_records_unmapped"] == 1
    audit_path = Path(str(result["output_path"]) + ".gff-grouping.json")
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    assert audit["stats"]["mapping_fallback_tolerated"] == 1
    assert audit["stats"]["unexpected_mapping_records"] == 1

    strict_task = dict(task, format_strict=True)
    strict_output_dir = tmp_path / "strict-out"
    strict_output_dir.mkdir()
    with pytest.raises(ValueError, match=r"unexpected_unmapped=1 ambiguous=0"):
        module.format_cds(strict_task, strict_output_dir, overwrite=False, dry_run=False)


def test_gfacs_bare_column9_ids_preserve_gene_boundaries_and_are_normalized(tmp_path):
    input_dir = tmp_path / "Direct" / "species_wise_original"
    species_dir = input_dir / "Sequoiadendron_giganteum"
    species_dir.mkdir(parents=True, exist_ok=True)
    (species_dir / "Segi.2_0.cds.fa").write_text(
        ">SEGI_00001\nATGAAATTT\n>SEGI_00002\nATGCCCTTT\n",
        encoding="utf-8",
    )
    (species_dir / "Segi.2_0.gtf").write_text(
        (
            "chr3\tGFACS\tgene\t1\t100\t.\t+\t.\tSEGI_00001\n"
            "chr3\tGFACS\tCDS\t1\t100\t.\t+\t0\tSEGI_00001\n"
            "chr3\tGFACS\tgene\t5\t100\t.\t+\t.\tSEGI_00002\n"
            "chr3\tGFACS\tCDS\t5\t100\t.\t+\t0\tSEGI_00002\n"
        ),
        encoding="utf-8",
    )
    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    completed = run_script(
        "--provider",
        "direct",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(tmp_path / "species_genome"),
        "--gene-grouping-mode",
        "rescue_overlap",
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = next(out_cds.glob("*.fa.gz"))
    audit = json.loads(Path(str(formatted_cds) + ".gff-grouping.json").read_text(encoding="utf-8"))
    assert audit["stats"]["mapped"] == 2
    assert audit["stats"]["ambiguous"] == 0
    assert audit["stats"]["coordinate_rescued_transcripts"] == 0
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        cds_text = handle.read()
    assert cds_text.count(">Sequoiadendron_giganteum_SEGI_") == 2

    formatted_gff = next(out_gff.glob("*.gff.gz"))
    with gzip.open(formatted_gff, "rt", encoding="utf-8") as handle:
        gff_text = handle.read()
    assert "\tgene\t1\t100\t.\t+\t.\tID=SEGI_00001;gene_id=SEGI_00001" in gff_text
    assert "\tCDS\t1\t100\t.\t+\t0\tParent=SEGI_00001;gene_id=SEGI_00001" in gff_text
    repair_audit = json.loads(Path(str(formatted_gff) + ".repair.json").read_text(encoding="utf-8"))
    assert repair_audit["normalized_bare_attribute_lines"] == 4
    validation = run_validate_mapping_script(
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
    )
    assert validation.returncode == 0, validation.stderr + "\n" + validation.stdout
    assert "CDS-to-GFF mapping OK: 2/2 IDs" in validation.stdout


def test_invalid_utf8_in_gff_attributes_is_replaced_and_audited(tmp_path):
    input_dir = tmp_path / "Direct" / "species_wise_original"
    species_dir = input_dir / "Pinus_tabuliformis"
    species_dir.mkdir(parents=True, exist_ok=True)
    (species_dir / "models.cds.fa").write_text(">tx1\nATGAAATTT\n", encoding="utf-8")
    (species_dir / "models.gff3").write_bytes(
        b"chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene1;Name=PtGT\xa6\xc3-N\n"
        b"chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=tx1;Parent=gene1\n"
        b"chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tParent=tx1\n"
    )
    out_gff = tmp_path / "species_gff"
    completed = run_script(
        "--provider",
        "direct",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(tmp_path / "species_cds"),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(tmp_path / "species_genome"),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    assert "replaced 2 invalid UTF-8 byte(s)" in completed.stderr

    formatted_gff = next(out_gff.glob("*.gff.gz"))
    with gzip.open(formatted_gff, "rt", encoding="utf-8") as handle:
        assert "Name=PtGT��-N" in handle.read()
    audit = json.loads(Path(str(formatted_gff) + ".repair.json").read_text(encoding="utf-8"))
    assert audit["invalid_utf8_bytes"] == 2
    assert audit["invalid_utf8_line_count"] == 1
    assert audit["invalid_utf8_lines"] == [1]


def test_gff_grouping_rescue_preserves_declared_parent_genes(tmp_path):
    module = load_module()
    gff_path = tmp_path / "models.gff3"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t18\t.\t+\t.\tID=badGeneA",
                "chr1\tsrc\tmRNA\t1\t18\t.\t+\t.\tID=locusX.t1;Parent=badGeneA",
                "chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tID=cds1;Parent=locusX.t1",
                "chr1\tsrc\tCDS\t13\t18\t.\t+\t0\tID=cds2;Parent=locusX.t1",
                "chr1\tsrc\tgene\t1\t18\t.\t+\t.\tID=badGeneB",
                "chr1\tsrc\tmRNA\t1\t18\t.\t+\t.\tID=locusX.t2;Parent=badGeneB",
                "chr1\tsrc\tCDS\t1\t3\t.\t+\t0\tID=cds3;Parent=locusX.t2",
                "chr1\tsrc\tCDS\t13\t18\t.\t+\t0\tID=cds4;Parent=locusX.t2",
                "",
            ]
        ),
        encoding="utf-8",
    )
    base_task = {
        "provider": "direct",
        "species_prefix": "Arabidopsis_thaliana",
        "gff_path": gff_path,
    }

    strict_index = module.build_gff_cds_grouping_index({**base_task, "gene_grouping_mode": "strict"})
    rescue_index = module.build_gff_cds_grouping_index({**base_task, "gene_grouping_mode": "rescue_overlap"})
    strict_a = module.resolve_cds_header_gff_gene(base_task, "locusX.t1", strict_index)
    strict_b = module.resolve_cds_header_gff_gene(base_task, "locusX.t2", strict_index)
    rescue_a = module.resolve_cds_header_gff_gene(base_task, "locusX.t1", rescue_index)
    rescue_b = module.resolve_cds_header_gff_gene(base_task, "locusX.t2", rescue_index)

    assert {strict_a["gene_token"], strict_b["gene_token"]} == {"badGeneA", "badGeneB"}
    assert {rescue_a["gene_token"], rescue_b["gene_token"]} == {"badGeneA", "badGeneB"}
    assert rescue_index["coordinate_rescued_transcripts"] == 0
    assert rescue_index["coordinate_rescued_groups"] == 0


def test_gff_grouping_preserves_unmatched_terminal_quotes_in_feature_ids(tmp_path):
    module = load_module()
    attrs = module.parse_gff_attributes(
        "ID=gene-TIL';Name=gene-TIL%27;Alias=\"balanced value\";gene_id='quoted-gene'"
    )
    assert attrs["ID"] == ("gene-TIL'",)
    assert attrs["Name"] == ("gene-TIL'",)
    assert attrs["Alias"] == ("balanced value",)
    assert attrs["gene_id"] == ("quoted-gene",)

    gff_path = tmp_path / "literal-apostrophe.gff3"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene-TIL;gene=gene-TIL",
                "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=tx-plain;Parent=gene-TIL",
                "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tParent=tx-plain",
                "chr1\tsrc\tgene\t20\t28\t.\t+\t.\tID=gene-TIL';gene=gene-TIL'",
                "chr1\tsrc\tmRNA\t20\t28\t.\t+\t.\tID=tx-apostrophe;Parent=gene-TIL'",
                "chr1\tsrc\tCDS\t20\t28\t.\t+\t0\tParent=tx-apostrophe",
                "",
            ]
        ),
        encoding="utf-8",
    )
    task = {
        "provider": "direct",
        "species_prefix": "Test_species",
        "gff_path": gff_path,
        "gene_grouping_mode": "strict",
    }

    index = module.build_gff_cds_grouping_index(task)

    plain = module.resolve_cds_header_gff_gene(task, "tx-plain", index)
    apostrophe = module.resolve_cds_header_gff_gene(
        task,
        "tx-apostrophe",
        index,
    )
    assert plain["gene_token"] == "TIL"
    assert apostrophe["gene_token"] == "TIL'"


def test_gff_grouping_rescue_preserves_distinct_authoritative_loci(tmp_path):
    module = load_module()
    gff_path = tmp_path / "distinct-loci.gff3"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t18\t.\t+\t.\tID=gene-L1;locus_tag=L1",
                "chr1\tsrc\tmRNA\t1\t18\t.\t+\t.\tID=L1.t1;Parent=gene-L1",
                "chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tParent=L1.t1",
                "chr1\tsrc\tCDS\t13\t18\t.\t+\t0\tParent=L1.t1",
                "chr1\tsrc\tgene\t1\t18\t.\t+\t.\tID=gene-L2;locus_tag=L2",
                "chr1\tsrc\tmRNA\t1\t18\t.\t+\t.\tID=L2.t1;Parent=gene-L2",
                "chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tParent=L2.t1",
                "chr1\tsrc\tCDS\t13\t18\t.\t+\t0\tParent=L2.t1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    task = {
        "provider": "direct",
        "species_prefix": "Test_species",
        "gff_path": gff_path,
        "gene_grouping_mode": "rescue_overlap",
    }

    index = module.build_gff_cds_grouping_index(task)

    left = module.resolve_cds_header_gff_gene(task, "L1.t1", index)
    right = module.resolve_cds_header_gff_gene(task, "L2.t1", index)
    assert {left["gene_token"], right["gene_token"]} == {"L1", "L2"}
    assert index["coordinate_rescued_transcripts"] == 0
    assert index["coordinate_rescued_groups"] == 0


def test_gff_grouping_rescue_does_not_bridge_authoritative_loci_through_unlabeled_model(tmp_path):
    module = load_module()
    gff_path = tmp_path / "bridged-loci.gff3"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t18\t.\t+\t.\tID=gene-L1;locus_tag=L1",
                "chr1\tsrc\tmRNA\t1\t18\t.\t+\t.\tID=L1.t1;Parent=gene-L1",
                "chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tParent=L1.t1",
                "chr1\tsrc\tCDS\t13\t18\t.\t+\t0\tParent=L1.t1",
                "chr1\tsrc\tgene\t1\t18\t.\t+\t.\tID=unlabeledGene",
                "chr1\tsrc\tmRNA\t1\t18\t.\t+\t.\tID=bridge.t1;Parent=unlabeledGene",
                "chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tParent=bridge.t1",
                "chr1\tsrc\tCDS\t13\t18\t.\t+\t0\tParent=bridge.t1",
                "chr1\tsrc\tgene\t1\t18\t.\t+\t.\tID=gene-L2;locus_tag=L2",
                "chr1\tsrc\tmRNA\t1\t18\t.\t+\t.\tID=L2.t1;Parent=gene-L2",
                "chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tParent=L2.t1",
                "chr1\tsrc\tCDS\t13\t18\t.\t+\t0\tParent=L2.t1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    task = {
        "provider": "direct",
        "species_prefix": "Test_species",
        "gff_path": gff_path,
        "gene_grouping_mode": "rescue_overlap",
    }

    index = module.build_gff_cds_grouping_index(task)

    gene_tokens = {
        module.resolve_cds_header_gff_gene(task, transcript_id, index)["gene_token"]
        for transcript_id in ("L1.t1", "bridge.t1", "L2.t1")
    }
    assert gene_tokens == {"L1", "unlabeledGene", "L2"}
    assert index["coordinate_rescued_transcripts"] == 0
    assert index["coordinate_rescued_groups"] == 0


def test_gff_grouping_does_not_merge_overlapping_opposite_strands(tmp_path):
    module = load_module()
    gff_path = tmp_path / "models.gff3"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t18\t.\t+\t.\tID=geneA",
                "chr1\tsrc\tmRNA\t1\t18\t.\t+\t.\tID=locusY.t1;Parent=geneA",
                "chr1\tsrc\tCDS\t1\t18\t.\t+\t0\tParent=locusY.t1",
                "chr1\tsrc\tgene\t1\t18\t.\t-\t.\tID=geneB",
                "chr1\tsrc\tmRNA\t1\t18\t.\t-\t.\tID=locusY.t2;Parent=geneB",
                "chr1\tsrc\tCDS\t1\t18\t.\t-\t0\tParent=locusY.t2",
                "",
            ]
        ),
        encoding="utf-8",
    )
    task = {
        "provider": "direct",
        "species_prefix": "Arabidopsis_thaliana",
        "gff_path": gff_path,
        "gene_grouping_mode": "rescue_overlap",
    }
    index = module.build_gff_cds_grouping_index(task)

    left = module.resolve_cds_header_gff_gene(task, "locusY.t1", index)
    right = module.resolve_cds_header_gff_gene(task, "locusY.t2", index)
    assert {left["gene_token"], right["gene_token"]} == {"geneA", "geneB"}
    assert index["coordinate_rescued_transcripts"] == 0


def test_gff_grouping_prefers_locus_identity_over_shared_gene_symbol(tmp_path):
    module = load_module()
    gff_path = tmp_path / "shared-symbol.gff3"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene-L1;gene=SHARED;locus_tag=L1",
                "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=tx1;Parent=gene-L1;locus_tag=L1",
                "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=cds1;Parent=tx1;gene=SHARED;locus_tag=L1;protein_id=P1",
                "chr1\tsrc\tgene\t20\t28\t.\t+\t.\tID=gene-L2;gene=SHARED;locus_tag=L2",
                "chr1\tsrc\tmRNA\t20\t28\t.\t+\t.\tID=tx2;Parent=gene-L2;locus_tag=L2",
                "chr1\tsrc\tCDS\t20\t28\t.\t+\t0\tID=cds2;Parent=tx2;gene=SHARED;locus_tag=L2;protein_id=P2",
                "",
            ]
        ),
        encoding="utf-8",
    )
    task = {
        "provider": "direct",
        "species_prefix": "Test_species",
        "gff_path": gff_path,
        "gene_grouping_mode": "strict",
    }
    index = module.build_gff_cds_grouping_index(task)

    assert module.resolve_cds_header_gff_gene(task, "tx1", index)["gene_token"] == "L1"
    assert module.resolve_cds_header_gff_gene(task, "P1", index)["gene_token"] == "L1"
    assert module.resolve_cds_header_gff_gene(task, "tx2", index)["gene_token"] == "L2"
    assert module.resolve_cds_header_gff_gene(task, "P2", index)["gene_token"] == "L2"


def test_gff_grouping_normalizes_separators_without_hiding_collisions(tmp_path):
    module = load_module()
    gff_path = tmp_path / "separator-aliases.gff3"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=riceGene",
                "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=LOC_Os10g36420.1;Parent=riceGene",
                "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tParent=LOC_Os10g36420.1",
                "chr1\tsrc\tgene\t20\t28\t.\t+\t.\tID=cephalotusGene",
                "chr1\tsrc\tmRNA\t20\t28\t.\t+\t.\tID=Cfol_v3_15267;Parent=cephalotusGene",
                "chr1\tsrc\tCDS\t20\t28\t.\t+\t0\tParent=Cfol_v3_15267",
                "chr1\tsrc\tgene\t40\t48\t.\t+\t.\tID=collisionA",
                "chr1\tsrc\tmRNA\t40\t48\t.\t+\t.\tID=X-A;Parent=collisionA",
                "chr1\tsrc\tCDS\t40\t48\t.\t+\t0\tParent=X-A",
                "chr1\tsrc\tgene\t60\t68\t.\t+\t.\tID=collisionB",
                "chr1\tsrc\tmRNA\t60\t68\t.\t+\t.\tID=X_A;Parent=collisionB",
                "chr1\tsrc\tCDS\t60\t68\t.\t+\t0\tParent=X_A",
                "",
            ]
        ),
        encoding="utf-8",
    )
    task = {
        "provider": "direct",
        "species_prefix": "Test_species",
        "gff_path": gff_path,
        "gene_grouping_mode": "strict",
    }
    index = module.build_gff_cds_grouping_index(task)

    assert module.resolve_cds_header_gff_gene(task, "LOC-Os10g36420.1", index)["gene_token"] == "riceGene"
    assert module.resolve_cds_header_gff_gene(task, "Cfol-v3-15267", index)["gene_token"] == "cephalotusGene"
    assert module.resolve_cds_header_gff_gene(task, "X-A", index)["gene_token"] == "collisionA"
    collision = module.resolve_cds_header_gff_gene(task, "X__A", index)
    assert collision["status"] == "ambiguous"
    assert set(collision["candidate_gene_tokens"]) == {"collisionA", "collisionB"}


def test_gff_grouping_uses_gene_alias_to_disambiguate_shared_protein_id(tmp_path):
    module = load_module()
    gff_path = tmp_path / "shared-protein.gff3"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene-L1;locus_tag=L1",
                "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=tx1;Parent=gene-L1;locus_tag=L1",
                "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tParent=tx1;locus_tag=L1;protein_id=P_SHARED",
                "chr1\tsrc\tgene\t20\t28\t.\t+\t.\tID=gene-L2;locus_tag=L2",
                "chr1\tsrc\tmRNA\t20\t28\t.\t+\t.\tID=tx2;Parent=gene-L2;locus_tag=L2",
                "chr1\tsrc\tCDS\t20\t28\t.\t+\t0\tParent=tx2;locus_tag=L2;protein_id=P_SHARED",
                "",
            ]
        ),
        encoding="utf-8",
    )
    task = {
        "provider": "ncbi",
        "species_prefix": "Test_species",
        "gff_path": gff_path,
        "gene_grouping_mode": "strict",
    }
    index = module.build_gff_cds_grouping_index(task)
    hit = module.resolve_cds_header_gff_gene(
        task,
        "record [protein_id=P_SHARED] [locus_tag=L2]",
        index,
    )

    assert hit["status"] == "mapped"
    assert hit["gene_token"] == "L2"


def test_gff_grouping_prefers_gtf_gene_id_over_transcript_dbxref(tmp_path):
    module = load_module()
    gtf_path = tmp_path / "models.gtf"
    gtf_path.write_text(
        "\n".join(
            [
                'chr1\tsrc\ttranscript\t1\t9\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; db_xref "Ensembl:ENST000001";',
                'chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tgene_id "G1"; transcript_id "T1"; protein_id "P1";',
                (
                    'chr1\tsrc\ttranscript\t20\t31\t.\t+\t.\tgene_id "G1"; '
                    'transcript_id "T2"; db_xref "Ensembl:ENST000002";'
                ),
                'chr1\tsrc\tCDS\t20\t31\t.\t+\t0\tgene_id "G1"; transcript_id "T2"; protein_id "P2";',
                "",
            ]
        ),
        encoding="utf-8",
    )
    task = {
        "provider": "ensemblplants",
        "species_prefix": "Test_species",
        "gff_path": gtf_path,
        "gene_grouping_mode": "strict",
    }

    index = module.build_gff_cds_grouping_index(task)

    assert module.resolve_cds_header_gff_gene(task, "P1", index)["gene_token"] == "G1"
    assert module.resolve_cds_header_gff_gene(task, "P2", index)["gene_token"] == "G1"
    assert module.resolve_cds_header_gff_gene(task, "ENST000001", index)["gene_token"] == "G1"


def test_gff_grouping_maps_polypeptide_and_protein_derives_from_aliases(tmp_path):
    module = load_module()
    gff_path = tmp_path / "protein-aliases.gff3"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t30\t.\t+\t.\tID=G1",
                "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=T1;Parent=G1",
                "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tParent=T1",
                "chr1\tsrc\tpolypeptide\t1\t9\t.\t+\t.\tID=P1;Derives_from=T1",
                "chr1\tsrc\tmRNA\t16\t30\t.\t+\t.\tID=T2;Parent=G1",
                "chr1\tsrc\tCDS\t16\t30\t.\t+\t0\tParent=T2",
                "chr1\tsrc\tprotein\t16\t30\t.\t+\t.\tID=P2;Derives_from=T2",
                "",
            ]
        ),
        encoding="utf-8",
    )
    task = {
        "provider": "direct",
        "species_prefix": "Test_species",
        "gff_path": gff_path,
        "gene_grouping_mode": "strict",
    }

    index = module.build_gff_cds_grouping_index(task)

    assert module.resolve_cds_header_gff_gene(task, "P1", index)["gene_token"] == "G1"
    assert module.resolve_cds_header_gff_gene(task, "P2", index)["gene_token"] == "G1"


def test_gff_grouping_rejects_conflicting_unique_header_aliases(tmp_path):
    module = load_module()
    cds_path = tmp_path / "conflicting-aliases.cds.fa"
    gff_path = tmp_path / "conflicting-aliases.gff3"
    cds_path.write_text(">PA [locus_tag=L2]\nATGAAATTT\n", encoding="utf-8")
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene-L1;locus_tag=L1",
                "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=T1;Parent=gene-L1",
                "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tParent=T1;protein_id=PA",
                "chr1\tsrc\tgene\t20\t28\t.\t+\t.\tID=gene-L2;locus_tag=L2",
                "chr1\tsrc\tmRNA\t20\t28\t.\t+\t.\tID=T2;Parent=gene-L2",
                "chr1\tsrc\tCDS\t20\t28\t.\t+\t0\tParent=T2;protein_id=PB",
                "",
            ]
        ),
        encoding="utf-8",
    )
    task = {
        "provider": "ncbi",
        "species_key": "Test_species",
        "species_prefix": "Test_species",
        "cds_path": cds_path,
        "gff_path": gff_path,
        "gene_grouping_mode": "strict",
    }

    index = module.build_gff_cds_grouping_index(task)
    hit = module.resolve_cds_header_gff_gene(task, "PA [locus_tag=L2]", index)

    assert hit["status"] == "ambiguous"
    assert set(hit["candidate_gene_tokens"]) == {"L1", "L2"}
    with pytest.raises(ValueError, match="unexpected_unmapped=0 ambiguous=1"):
        module.format_cds(task, tmp_path / "out", overwrite=False, dry_run=False, strict=True)


def test_gff_grouping_rejects_sanitized_gene_id_collisions(tmp_path):
    module = load_module()
    cds_path = tmp_path / "models.cds.fa"
    gff_path = tmp_path / "models.gff3"
    cds_path.write_text(">T1\nATGAAATTT\n>T2\nATGCCCAAATTT\n", encoding="utf-8")
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=A:B",
                "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=T1;Parent=A:B",
                "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tParent=T1",
                "chr1\tsrc\tgene\t20\t31\t.\t+\t.\tID=A%2FB",
                "chr1\tsrc\tmRNA\t20\t31\t.\t+\t.\tID=T2;Parent=A%2FB",
                "chr1\tsrc\tCDS\t20\t31\t.\t+\t0\tParent=T2",
                "",
            ]
        ),
        encoding="utf-8",
    )
    task = {
        "provider": "direct",
        "species_key": "Test_species",
        "species_prefix": "Test_species",
        "cds_path": cds_path,
        "gff_path": gff_path,
        "gene_grouping_mode": "strict",
    }

    with pytest.raises(ValueError, match="collide after identifier sanitization"):
        module.format_cds(task, tmp_path / "out", overwrite=False, dry_run=False)


@pytest.mark.parametrize("gene_grouping_mode", ("strict", "rescue_overlap"))
def test_gff_grouping_rejects_gene_feature_prefix_collisions(tmp_path, gene_grouping_mode):
    module = load_module()
    gff_path = tmp_path / "prefix-collision.gff3"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene-A",
                "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=T1;Parent=gene-A",
                "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tParent=T1",
                "chr1\tsrc\tgene\t20\t28\t.\t+\t.\tID=A",
                "chr1\tsrc\tmRNA\t20\t28\t.\t+\t.\tID=T2;Parent=A",
                "chr1\tsrc\tCDS\t20\t28\t.\t+\t0\tParent=T2",
                "",
            ]
        ),
        encoding="utf-8",
    )
    task = {
        "provider": "direct",
        "species_prefix": "Test_species",
        "gff_path": gff_path,
        "gene_grouping_mode": gene_grouping_mode,
    }

    with pytest.raises(ValueError, match="collapse after prefix normalization"):
        module.build_gff_cds_grouping_index(task)


def test_gff_grouping_allows_prefix_variants_confirmed_by_same_locus_tag(tmp_path):
    module = load_module()
    gff_path = tmp_path / "confirmed-prefix-variants.gff3"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene-A;locus_tag=A",
                "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=T1;Parent=gene-A",
                "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tParent=T1",
                "chr1\tsrc\tgene\t20\t28\t.\t+\t.\tID=A;locus_tag=A",
                "chr1\tsrc\tmRNA\t20\t28\t.\t+\t.\tID=T2;Parent=A",
                "chr1\tsrc\tCDS\t20\t28\t.\t+\t0\tParent=T2",
                "",
            ]
        ),
        encoding="utf-8",
    )
    task = {
        "provider": "direct",
        "species_prefix": "Test_species",
        "gff_path": gff_path,
        "gene_grouping_mode": "strict",
    }

    index = module.build_gff_cds_grouping_index(task)

    assert module.resolve_cds_header_gff_gene(task, "T1", index)["gene_token"] == "A"
    assert module.resolve_cds_header_gff_gene(task, "T2", index)["gene_token"] == "A"


@pytest.mark.parametrize("gene_grouping_mode", ("strict", "rescue_overlap"))
@pytest.mark.parametrize("parent_order", ("gene-L1,gene-L2", "gene-L2,gene-L1"))
def test_gff_grouping_marks_distinct_multi_parent_genes_ambiguous(
    tmp_path,
    gene_grouping_mode,
    parent_order,
):
    module = load_module()
    gff_path = tmp_path / "multi-parent.gff3"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t30\t.\t+\t.\tID=gene-L1;locus_tag=L1",
                "chr1\tsrc\tgene\t1\t30\t.\t+\t.\tID=gene-L2;locus_tag=L2",
                f"chr1\tsrc\tmRNA\t1\t30\t.\t+\t.\tID=T1;Parent={parent_order}",
                "chr1\tsrc\tCDS\t1\t30\t.\t+\t0\tParent=T1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    task = {
        "provider": "direct",
        "species_prefix": "Test_species",
        "gff_path": gff_path,
        "gene_grouping_mode": gene_grouping_mode,
    }

    index = module.build_gff_cds_grouping_index(task)
    resolved = module.resolve_cds_header_gff_gene(task, "T1", index)

    assert resolved["status"] == "ambiguous"
    assert resolved["gene_token"] == ""
    assert resolved["candidate_gene_tokens"] == ("L1", "L2")
    assert index["ambiguous_transcript_gene_tokens"] == {"T1": ("L1", "L2")}


def test_provided_cds_strict_format_rejects_distinct_multi_parent_genes(tmp_path):
    module = load_module()
    cds_path = tmp_path / "multi-parent.cds.fa"
    gff_path = tmp_path / "multi-parent.gff3"
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    cds_path.write_text(">T1\nATGAAAAAA\n", encoding="utf-8")
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene-L1;locus_tag=L1",
                "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene-L2;locus_tag=L2",
                "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=T1;Parent=gene-L1,gene-L2",
                "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tParent=T1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    task = {
        "provider": "direct",
        "species_key": "Test_species",
        "species_prefix": "Test_species",
        "cds_path": cds_path,
        "gff_path": gff_path,
        "gene_grouping_mode": "strict",
    }

    with pytest.raises(ValueError, match="ambiguous=1"):
        module.format_cds(task, output_dir, overwrite=False, dry_run=False, strict=True)


@pytest.mark.parametrize("reverse_definitions", (False, True))
def test_gff_grouping_rejects_conflicting_duplicate_feature_ids(tmp_path, reverse_definitions):
    module = load_module()
    gff_path = tmp_path / "duplicate-id.gff3"
    transcript_definitions = [
        "chr1\tsrc\tmRNA\t1\t30\t.\t+\t.\tID=T1;Parent=gene-L1",
        "chr1\tsrc\tmRNA\t1\t30\t.\t+\t.\tID=T1;Parent=gene-L2",
    ]
    if reverse_definitions:
        transcript_definitions.reverse()
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t30\t.\t+\t.\tID=gene-L1;locus_tag=L1",
                "chr1\tsrc\tgene\t1\t30\t.\t+\t.\tID=gene-L2;locus_tag=L2",
                *transcript_definitions,
                "chr1\tsrc\tCDS\t1\t30\t.\t+\t0\tParent=T1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    task = {
        "provider": "direct",
        "species_prefix": "Test_species",
        "gff_path": gff_path,
        "gene_grouping_mode": "strict",
    }

    with pytest.raises(ValueError, match="conflicting definitions for feature ID T1.*parents"):
        module.build_gff_cds_grouping_index(task)


def test_gff_grouping_merges_compatible_duplicate_feature_ids(tmp_path):
    module = load_module()
    gff_path = tmp_path / "compatible-duplicate-id.gff3"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t30\t.\t+\t.\tID=gene-A;locus_tag=A",
                "chr1\tsrc\tgene\t1\t30\t.\t+\t.\tID=alternate-A;locus_tag=A",
                "chr1\tsrc\tmRNA\t1\t30\t.\t+\t.\tID=T1;Parent=gene-A,alternate-A;Alias=T1a",
                "chr1\tsrc\tmRNA\t1\t30\t.\t+\t.\tID=T1;Parent=alternate-A,gene-A;Alias=T1b",
                "chr1\tsrc\tCDS\t1\t30\t.\t+\t0\tParent=T1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    task = {
        "provider": "direct",
        "species_prefix": "Test_species",
        "gff_path": gff_path,
        "gene_grouping_mode": "strict",
    }

    index = module.build_gff_cds_grouping_index(task)

    assert module.resolve_cds_header_gff_gene(task, "T1a", index)["gene_token"] == "A"
    assert module.resolve_cds_header_gff_gene(task, "T1b", index)["gene_token"] == "A"


@pytest.mark.parametrize('strand', ['+', '-'])
def test_ordered_multipart_gene_keeps_one_longest_cds(tmp_path, strand):
    module = load_module()
    gff = tmp_path / 'ordered.gff'
    gff.write_text(
        f'chr1\tDDBJ\tgene\t1\t9\t.\t{strand}\t.\tID=gene-A;locus_tag=A;is_ordered=true;partial=true\n'
        f'chr1\tDDBJ\tgene\t31\t39\t.\t{strand}\t.\tID=gene-A;locus_tag=A;is_ordered=true;partial=true\n'
        f'chr1\tDDBJ\tCDS\t1\t9\t.\t{strand}\t0\tID=cds-P1;Parent=gene-A;protein_id=P1\n'
        f'chr1\tDDBJ\tCDS\t31\t39\t.\t{strand}\t0\tID=cds-P1;Parent=gene-A;protein_id=P1\n'
        f'chr1\tDDBJ\tCDS\t31\t36\t.\t{strand}\t0\tID=cds-P2;Parent=gene-A;protein_id=P2\n')
    cds = tmp_path / 'cds.fa'
    cds.write_text('>P1\nATGATGATGATGATGATG\n>P2\nATGATG\n')
    task = dict(provider='ncbi', species_key='Test_species', species_prefix='Test_species',
                cds_path=cds, gff_path=gff)
    result = module.format_cds(task, tmp_path, False, False)
    assert gzip.open(result['output_path'], 'rt').read() == '>Test_species_A\nATGATGATGATGATGATG\n'
    audit = json.loads(Path(str(result['output_path']) + '.gff-grouping.json').read_text())
    assert audit['stats']['mapped'] == 2
    assert audit['after_count'] == 1


@pytest.mark.parametrize('second', [
    ('chr1', '+', '', 'A'), ('chr2', '+', ';is_ordered=true', 'A'),
    ('chr1', '-', ';is_ordered=true', 'A'), ('chr1', '+', ';is_ordered=true', 'B'),
])
def test_ordered_multipart_gene_rejects_conflicting_identity(tmp_path, second):
    module = load_module()
    seqid, strand, marker, token = second
    gff = tmp_path / 'conflicting.gff'
    gff.write_text('chr1\tDDBJ\tgene\t1\t9\t.\t+\t.\tID=gene-A;locus_tag=A;is_ordered=true\n' +
        f'{seqid}\tDDBJ\tgene\t31\t39\t.\t{strand}\t.\tID=gene-A;locus_tag={token}{marker}\n')
    with pytest.raises(ValueError, match='conflicting definitions for feature ID gene-A'):
        module.build_gff_cds_grouping_index(dict(provider='ncbi', species_prefix='Test_species', gff_path=gff))


@pytest.mark.parametrize("transcript_first", (False, True))
def test_gff_grouping_accepts_maker_gene_transcript_shared_id(tmp_path, transcript_first):
    module = load_module()
    gff_path = tmp_path / "maker-shared-id.gff3"
    shared_id_rows = [
        "chr1\tmaker\tgene\t1\t30\t.\t+\t.\tID=Dm-0001;Name=Dm-0001;Alias=gene-alias",
        "chr1\tmaker\tmRNA\t1\t30\t.\t+\t.\tID=Dm-0001;Parent=Dm-0001;Alias=mrna-alias",
    ]
    if transcript_first:
        shared_id_rows.reverse()
    gff_path.write_text(
        "\n".join(
            [
                *shared_id_rows,
                "chr1\tmaker\tCDS\t1\t30\t.\t+\t0\tID=Dm-0001:cds;Parent=Dm-0001",
                "",
            ]
        ),
        encoding="utf-8",
    )
    task = {
        "provider": "direct",
        "species_prefix": "Test_species",
        "gff_path": gff_path,
        "gene_grouping_mode": "strict",
    }

    index = module.build_gff_cds_grouping_index(task)

    assert module.resolve_cds_header_gff_gene(task, "Dm-0001", index)["gene_token"] == "Dm-0001"
    assert module.resolve_cds_header_gff_gene(task, "gene-alias", index)["gene_token"] == "Dm-0001"
    assert module.resolve_cds_header_gff_gene(task, "mrna-alias", index)["gene_token"] == "Dm-0001"


def test_gff_grouping_handles_deep_parent_graph_without_recursion(tmp_path):
    module = load_module()
    gff_path = tmp_path / "deep-parent-graph.gff3"
    depth = 1500
    rows = ["chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=gene-G;locus_tag=G"]
    for index in range(depth):
        parent = f"N{index + 1}" if index + 1 < depth else "gene-G"
        rows.append(f"chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=N{index};Parent={parent}")
    rows.extend(("chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tParent=N0", ""))
    gff_path.write_text("\n".join(rows), encoding="utf-8")
    task = {
        "provider": "direct",
        "species_prefix": "Test_species",
        "gff_path": gff_path,
        "gene_grouping_mode": "strict",
    }

    index = module.build_gff_cds_grouping_index(task)
    resolved = module.resolve_cds_header_gff_gene(task, "N0", index)

    assert resolved["status"] == "mapped"
    assert resolved["gene_token"] == "G"


@pytest.mark.parametrize("reverse_cds_order", (False, True))
def test_gff_grouping_cycle_fallback_is_input_order_independent(tmp_path, reverse_cds_order):
    module = load_module()
    gff_path = tmp_path / "cyclic-parent-graph.gff3"
    cds_rows = [
        "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tParent=T1",
        "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tParent=T2",
    ]
    if reverse_cds_order:
        cds_rows.reverse()
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=T1;Parent=T2",
                "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=T2;Parent=T1",
                *cds_rows,
                "",
            ]
        ),
        encoding="utf-8",
    )
    task = {
        "provider": "direct",
        "species_prefix": "Test_species",
        "gff_path": gff_path,
        "gene_grouping_mode": "strict",
    }

    index = module.build_gff_cds_grouping_index(task)

    assert module.resolve_cds_header_gff_gene(task, "T1", index)["gene_token"] == "T1"
    assert module.resolve_cds_header_gff_gene(task, "T2", index)["gene_token"] == "T1"


def test_gff_grouping_cycle_with_conflicting_gene_tokens_is_ambiguous(tmp_path):
    module = load_module()
    gff_path = tmp_path / "cyclic-conflicting-identity.gff3"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=T1;Parent=T2;gene=A",
                "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=T2;Parent=T1;gene=B",
                "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tParent=T1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    task = {
        "provider": "direct",
        "species_prefix": "Test_species",
        "gff_path": gff_path,
        "gene_grouping_mode": "strict",
    }

    index = module.build_gff_cds_grouping_index(task)
    resolved = module.resolve_cds_header_gff_gene(task, "T1", index)

    assert resolved["status"] == "ambiguous"
    assert resolved["candidate_gene_tokens"] == ("A", "B")


def test_gff_rescue_compares_matching_authoritative_loci():
    module = load_module()
    features = {
        "T1": [{"seqid": "chr1", "start": 1, "end": 30, "strand": "+", "gene_token": "badA"}],
        "T2": [{"seqid": "chr1", "start": 1, "end": 30, "strand": "+", "gene_token": "badB"}],
    }
    authoritative = {"T1": ("L1",), "T2": ("L1",)}
    task = {
        "provider": "direct",
        "species_prefix": "Test_species",
        "gene_grouping_mode": "rescue_overlap",
    }

    resolved = module.build_rescued_gene_tokens_for_transcripts(task, features, authoritative)

    assert resolved["T1"] == resolved["T2"] == "badA"


def test_gff_rescue_skips_pairwise_checks_between_distinct_authoritative_loci(monkeypatch):
    module = load_module()
    call_count = 0
    original = module.build_rescued_gene_tokens_for_transcripts.__globals__[
        "should_rescue_overlapping_transcripts"
    ]

    def count_calls(left, right):
        nonlocal call_count
        call_count += 1
        return original(left, right)

    monkeypatch.setitem(
        module.build_rescued_gene_tokens_for_transcripts.__globals__,
        "should_rescue_overlapping_transcripts",
        count_calls,
    )
    transcript_count = 500
    features = {
        f"T{index}": [
            {
                "seqid": "chr1",
                "start": 1,
                "end": 1000,
                "strand": "+",
                "gene_token": f"G{index}",
            }
        ]
        for index in range(transcript_count)
    }
    authoritative = {f"T{index}": (f"G{index}",) for index in range(transcript_count)}
    task = {
        "provider": "direct",
        "species_prefix": "Test_species",
        "gene_grouping_mode": "rescue_overlap",
    }

    resolved = module.build_rescued_gene_tokens_for_transcripts(task, features, authoritative)

    assert resolved == {f"T{index}": f"G{index}" for index in range(transcript_count)}
    assert call_count == 0


def test_provided_cds_longest_selection_compares_lengths_before_padding(tmp_path):
    module = load_module()
    cds_path = tmp_path / "models.cds.fa"
    gff_path = tmp_path / "models.gff3"
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    cds_path.write_text(
        ">A_short\nATGAAAAA\n>B_long\nATGAAAAAA\n",
        encoding="utf-8",
    )
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t20\t.\t+\t.\tID=G1",
                "chr1\tsrc\tmRNA\t1\t8\t.\t+\t.\tID=A_short;Parent=G1",
                "chr1\tsrc\tCDS\t1\t8\t.\t+\t0\tParent=A_short",
                "chr1\tsrc\tmRNA\t12\t20\t.\t+\t.\tID=B_long;Parent=G1",
                "chr1\tsrc\tCDS\t12\t20\t.\t+\t0\tParent=B_long",
                "",
            ]
        ),
        encoding="utf-8",
    )
    task = {
        "provider": "direct",
        "species_key": "Test_species",
        "species_prefix": "Test_species",
        "cds_path": cds_path,
        "gff_path": gff_path,
        "gene_grouping_mode": "strict",
    }

    result = module.format_cds(task, output_dir, overwrite=False, dry_run=False)

    with gzip.open(result["output_path"], "rt", encoding="utf-8") as handle:
        assert handle.read() == ">Test_species_G1\nATGAAAAAA\n"
    audit_json_path = Path(str(result["output_path"]) + ".gff-grouping.json")
    audit_tsv_path = Path(str(result["output_path"]) + ".gff-grouping.tsv")
    with open(audit_json_path, "rt", encoding="utf-8") as handle:
        audit = json.load(handle)
    with open(audit_tsv_path, "rt", encoding="utf-8", newline="") as handle:
        audit_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert audit["version"] == 15
    assert [row["raw_sequence_length"] for row in audit_rows] == ["8", "9"]
    assert [row["sequence_length"] for row in audit_rows] == ["9", "9"]
    assert [row["selected_longest"] for row in audit_rows] == ["0", "1"]


def test_provided_cds_gff_grouping_regenerates_older_audit_version(tmp_path):
    module = load_module()
    cds_path = tmp_path / "models.cds.fa"
    gff_path = tmp_path / "models.gff3"
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    cds_path.write_text(">T1\nATGAAAAAA\n", encoding="utf-8")
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=G1",
                "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=T1;Parent=G1",
                "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tParent=T1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    task = {
        "provider": "direct",
        "species_key": "Test_species",
        "species_prefix": "Test_species",
        "cds_path": cds_path,
        "gff_path": gff_path,
        "gene_grouping_mode": "strict",
    }

    first = module.format_cds(task, output_dir, overwrite=False, dry_run=False)
    audit_path = Path(str(first["output_path"]) + ".gff-grouping.json")
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    audit["version"] = 4
    audit_path.write_text(json.dumps(audit), encoding="utf-8")

    regenerated = module.format_cds(task, output_dir, overwrite=False, dry_run=False)
    skipped = module.format_cds(task, output_dir, overwrite=False, dry_run=False)

    assert regenerated["status"] == "write"
    assert json.loads(audit_path.read_text(encoding="utf-8"))["version"] == 15
    assert skipped["status"] == "skip"


def test_format_species_inputs_preserves_cds_despite_overlapping_utrs(tmp_path):
    input_dir = tmp_path / "Direct" / "species_wise_original"
    species_dir = input_dir / "Arabidopsis_thaliana"
    species_dir.mkdir(parents=True, exist_ok=True)
    gff_path = species_dir / "Arabidopsis_thaliana.annotation.gff3"
    genome_path = species_dir / "Arabidopsis_thaliana.genome.fa"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t15\t.\t+\t.\tID=gene1",
                "chr1\tsrc\tmRNA\t1\t15\t.\t+\t.\tID=gene1.t1;Parent=gene1",
                "chr1\tsrc\tfive_prime_UTR\t1\t6\t.\t+\t.\tParent=gene1.t1",
                "chr1\tsrc\tCDS\t1\t15\t.\t+\t0\tID=cds1;Parent=gene1.t1",
                "chr2\tsrc\tgene\t1\t15\t.\t-\t.\tID=gene2",
                "chr2\tsrc\tmRNA\t1\t15\t.\t-\t.\tID=gene2.t1;Parent=gene2",
                "chr2\tsrc\tfive_prime_UTR\t10\t15\t.\t-\t.\tParent=gene2.t1",
                "chr2\tsrc\tCDS\t1\t15\t.\t-\t0\tID=cds2;Parent=gene2.t1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    genome_path.write_text(
        ">chr1\nATGAAACCCGGGTTT\n>chr2\nAAACCCGGGTTTCAT\n",
        encoding="utf-8",
    )

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "direct",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Arabidopsis_thaliana_annotation.derived.cds.fa.gz"
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert ">Arabidopsis_thaliana_gene1" in text
    assert ">Arabidopsis_thaliana_gene2" in text
    assert "CCCGGGTTT" in text
    assert text.count("ATGAAACCCGGGTTT") == 2


def test_format_species_inputs_rescue_preserves_different_parent_genes(tmp_path):
    input_dir = tmp_path / "Direct" / "species_wise_original"
    species_dir = input_dir / "Arabidopsis_thaliana"
    species_dir.mkdir(parents=True, exist_ok=True)
    gff_path = species_dir / "Arabidopsis_thaliana.annotation.gff3"
    genome_path = species_dir / "Arabidopsis_thaliana.genome.fa"
    species_summary = tmp_path / "gg_input_generation_species.tsv"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t18\t.\t+\t.\tID=badGeneA",
                "chr1\tsrc\tmRNA\t1\t18\t.\t+\t.\tID=locusX.t1;Parent=badGeneA",
                "chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tID=cds1;Parent=locusX.t1",
                "chr1\tsrc\tCDS\t13\t18\t.\t+\t0\tID=cds2;Parent=locusX.t1",
                "chr1\tsrc\tgene\t1\t18\t.\t+\t.\tID=badGeneB",
                "chr1\tsrc\tmRNA\t1\t18\t.\t+\t.\tID=locusX.t2;Parent=badGeneB",
                "chr1\tsrc\tCDS\t1\t3\t.\t+\t0\tID=cds3;Parent=locusX.t2",
                "chr1\tsrc\tCDS\t13\t18\t.\t+\t0\tID=cds4;Parent=locusX.t2",
                "",
            ]
        ),
        encoding="utf-8",
    )
    genome_path.write_text(">chr1\nATGAAACCCGGGTTTCCC\n", encoding="utf-8")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "direct",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
        "--species-summary-output",
        str(species_summary),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Arabidopsis_thaliana_annotation.derived.cds.fa.gz"
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        headers = [line.strip() for line in handle if line.startswith(">")]
    assert headers == [">Arabidopsis_thaliana_badGeneA", ">Arabidopsis_thaliana_badGeneB"]

    with open(species_summary, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["gene_grouping_mode"] == "rescue_overlap"
    assert rows[0]["aggregated_cds_removed"] == "0"


def test_format_species_inputs_uses_locus_tag_for_genbank_style_ncbi_cds(tmp_path):
    input_dir = tmp_path / "NCBI_Genome" / "species_wise_original"
    species_dir = input_dir / "Dictyostelium_cf_discoideum"
    species_dir.mkdir(parents=True, exist_ok=True)
    cds_path = species_dir / "GCA_054859205.1_ASM5485920v1_cds_from_genomic.fna.gz"
    gff_path = species_dir / "GCA_054859205.1_ASM5485920v1_genomic.gff.gz"
    genome_path = species_dir / "GCA_054859205.1_ASM5485920v1_genomic.fna.gz"
    with gzip.open(cds_path, "wt", encoding="utf-8") as handle:
        handle.write(
            (
                ">lcl|JBTAPH010000036.1_cds_KAM9986187.1_1 [locus_tag=ACTFIY_010592] [protein=hypothetical protein] [protein_id=KAM9986187.1] [location=complement(join(9..158,229..966))] [gbkey=CDS]\n"
                "ATGTCTACCACTGTTAACAATAATGATGCCTCTAGTAGTAGTAGCTCTGCCTCTAATAACGATGAATCCTTTGATTTAAGAATGAAATCAATGGAGGATCAAATCAATAACCTTTCATTAGCCTTTACCAGATTCATGAAAGAACCTATGTTCTCTTCTAATACCAAATCACGTAGCCAACCTTCTCATGATAACTCTGACACTGAGAATGAACAAAGTGATGACGAATCAAGTAACAAT\n"
                ">lcl|JBTAPH010000036.1_cds_KAM9986188.1_2 [locus_tag=ACTFIY_010593] [protein=hypothetical protein] [protein_id=KAM9986188.1] [location=complement(join(2114..2192,2246..2381,2478..2691))] [gbkey=CDS]\n"
                "ATGGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATC\n"
            )
        )
    # Negative-strand block lengths 214, 136, 79 imply phases 0, 2, 1.
    with gzip.open(gff_path, "wt", encoding="utf-8") as handle:
        handle.write(
            "\n".join(
                [
                    "##gff-version 3",
                    "JBTAPH010000036.1\tGenbank\tgene\t9\t966\t.\t-\t.\tID=gene-ACTFIY_010592;Name=ACTFIY_010592;gbkey=Gene;gene_biotype=protein_coding;locus_tag=ACTFIY_010592",
                    "JBTAPH010000036.1\tGenbank\tmRNA\t9\t966\t.\t-\t.\tID=rna-mrna.DD_M4_00007442-RA:cds;Parent=gene-ACTFIY_010592;gbkey=mRNA;locus_tag=ACTFIY_010592;orig_protein_id=gnl|WGS:JBTAPH|DD_M4_00007442-RA:cds;orig_transcript_id=gnl|WGS:JBTAPH|mrna.DD_M4_00007442-RA:cds;product=hypothetical protein",
                    "JBTAPH010000036.1\tGenbank\tCDS\t229\t966\t.\t-\t0\tID=cds-KAM9986187.1;Parent=rna-mrna.DD_M4_00007442-RA:cds;Dbxref=NCBI_GP:KAM9986187.1;Name=KAM9986187.1;gbkey=CDS;locus_tag=ACTFIY_010592;orig_transcript_id=gnl|WGS:JBTAPH|mrna.DD_M4_00007442-RA:cds;product=hypothetical protein;protein_id=KAM9986187.1",
                    "JBTAPH010000036.1\tGenbank\tCDS\t9\t158\t.\t-\t0\tID=cds-KAM9986187.1;Parent=rna-mrna.DD_M4_00007442-RA:cds;Dbxref=NCBI_GP:KAM9986187.1;Name=KAM9986187.1;gbkey=CDS;locus_tag=ACTFIY_010592;orig_transcript_id=gnl|WGS:JBTAPH|mrna.DD_M4_00007442-RA:cds;product=hypothetical protein;protein_id=KAM9986187.1",
                    "JBTAPH010000036.1\tGenbank\tgene\t2114\t2691\t.\t-\t.\tID=gene-ACTFIY_010593;Name=ACTFIY_010593;gbkey=Gene;gene_biotype=protein_coding;locus_tag=ACTFIY_010593",
                    "JBTAPH010000036.1\tGenbank\tmRNA\t2114\t2691\t.\t-\t.\tID=rna-mrna.DD_M4_00007443-RA:cds;Parent=gene-ACTFIY_010593;gbkey=mRNA;locus_tag=ACTFIY_010593;orig_protein_id=gnl|WGS:JBTAPH|DD_M4_00007443-RA:cds;orig_transcript_id=gnl|WGS:JBTAPH|mrna.DD_M4_00007443-RA:cds;product=hypothetical protein",
                    "JBTAPH010000036.1\tGenbank\tCDS\t2478\t2691\t.\t-\t0\tID=cds-KAM9986188.1;Parent=rna-mrna.DD_M4_00007443-RA:cds;Dbxref=NCBI_GP:KAM9986188.1;Name=KAM9986188.1;gbkey=CDS;locus_tag=ACTFIY_010593;orig_transcript_id=gnl|WGS:JBTAPH|mrna.DD_M4_00007443-RA:cds;product=hypothetical protein;protein_id=KAM9986188.1",
                    "JBTAPH010000036.1\tGenbank\tCDS\t2246\t2381\t.\t-\t2\tID=cds-KAM9986188.1;Parent=rna-mrna.DD_M4_00007443-RA:cds;Dbxref=NCBI_GP:KAM9986188.1;Name=KAM9986188.1;gbkey=CDS;locus_tag=ACTFIY_010593;orig_transcript_id=gnl|WGS:JBTAPH|mrna.DD_M4_00007443-RA:cds;product=hypothetical protein;protein_id=KAM9986188.1",
                    "JBTAPH010000036.1\tGenbank\tCDS\t2114\t2192\t.\t-\t1\tID=cds-KAM9986188.1;Parent=rna-mrna.DD_M4_00007443-RA:cds;Dbxref=NCBI_GP:KAM9986188.1;Name=KAM9986188.1;gbkey=CDS;locus_tag=ACTFIY_010593;orig_transcript_id=gnl|WGS:JBTAPH|mrna.DD_M4_00007443-RA:cds;product=hypothetical protein;protein_id=KAM9986188.1",
                    "",
                ]
            )
        )
    with gzip.open(genome_path, "wt", encoding="utf-8") as handle:
        handle.write(">JBTAPH010000036.1\n" + ("A" * 3000) + "\n")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "ncbi",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Dictyostelium_cf_discoideum_GCA_054859205.1_ASM5485920v1_cds_from_genomic.fa.gz"
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        headers = [line.strip() for line in handle if line.startswith(">")]
    assert headers == [
        ">Dictyostelium_cf_discoideum_ACTFIY_010592",
        ">Dictyostelium_cf_discoideum_ACTFIY_010593",
    ]

    mapping = run_validate_mapping_script(
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
    )
    assert mapping.returncode == 0, mapping.stderr + "\n" + mapping.stdout
    assert "[Dictyostelium_cf_discoideum] CDS-to-GFF mapping OK: 2/2 IDs" in mapping.stdout


def test_format_species_inputs_excludes_only_unlinkable_anonymous_ncbi_cds(tmp_path):
    input_dir = tmp_path / "NCBI_Genome" / "species_wise_original"
    species_dir = input_dir / "Example_species"
    species_dir.mkdir(parents=True)
    cds_path = species_dir / "GCA_000000001.1_demo_cds_from_genomic.fna.gz"
    gff_path = species_dir / "GCA_000000001.1_demo_genomic.gff.gz"
    genome_path = species_dir / "GCA_000000001.1_demo_genomic.fna.gz"
    total_mappable = 999
    with gzip.open(cds_path, "wt", encoding="utf-8") as handle:
        for index in range(total_mappable):
            handle.write(
                ">lcl|NC_000001.1_cds_XP{}.1_{} [locus_tag=LOC{}] "
                "[protein_id=XP{}.1] [gbkey=CDS]\nATGAAA\n".format(
                    index, index + 1, index, index
                )
            )
        handle.write(
            ">lcl|NC_000001.1_cds_1000 [location=9001..9006] [gbkey=CDS]\nATGAAA\n"
        )
    with gzip.open(gff_path, "wt", encoding="utf-8") as handle:
        for index in range(total_mappable):
            handle.write(
                "NC_000001.1\tsrc\tCDS\t{}\t{}\t.\t+\t0\t"
                "ID=cds-XP{}.1;locus_tag=LOC{};protein_id=XP{}.1\n".format(
                    index * 6 + 1,
                    index * 6 + 6,
                    index,
                    index,
                    index,
                )
            )
    with gzip.open(genome_path, "wt", encoding="utf-8") as handle:
        handle.write(">NC_000001.1\n" + ("A" * 10000) + "\n")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    species_summary = tmp_path / "summary.tsv"
    completed = run_script(
        "--provider",
        "ncbi",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
        "--species-summary-output",
        str(species_summary),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = next(out_cds.glob("*.fa.gz"))
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        headers = [line.strip() for line in handle if line.startswith(">")]
    assert len(headers) == total_mappable
    assert all("NC_000001.1_cds_1000" not in header for header in headers)

    audit_path = Path(str(formatted_cds) + ".gff-grouping.tsv")
    with open(audit_path, "rt", encoding="utf-8", newline="") as handle:
        audit_rows = list(csv.DictReader(handle, delimiter="\t"))
    excluded = [row for row in audit_rows if row["mapping_status"] == "excluded_anonymous_unmapped"]
    assert len(excluded) == 1
    assert excluded[0]["exclusion_reason"] == "anonymous_ncbi_cds_without_gff_link"

    mapping = run_validate_mapping_script(
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--strict",
    )
    assert mapping.returncode == 0, mapping.stderr + "\n" + mapping.stdout
    assert "CDS-to-GFF mapping OK: 999/999 IDs" in mapping.stdout

    with open(species_summary, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["cds_gff_records_excluded_anonymous"] == "1"


@pytest.mark.parametrize("provider", ("ncbi", "direct"))
def test_format_species_inputs_maps_anonymous_ncbi_cds_by_exact_location(tmp_path, provider):
    input_dir = tmp_path / "NCBI_Genome" / "species_wise_original"
    species_dir = input_dir / "Clitoria_ternatea"
    species_dir.mkdir(parents=True)
    cds_path = species_dir / "GCA_037962975.1_Cter_1_cds_from_genomic.fna.gz"
    gff_path = species_dir / "GCA_037962975.1_Cter_1_genomic.gff.gz"
    genome_path = species_dir / "GCA_037962975.1_Cter_1_genomic.fna.gz"
    with gzip.open(cds_path, "wt", encoding="utf-8") as handle:
        handle.write(
            ">lcl|CM075669.1_cds_1 "
            "[db_xref=InterPro:IPR000504,InterPro:IPR035979] "
            "[protein=hypothetical protein] [pseudo=true] "
            "[location=complement(join(101..109,201..209))] [gbkey=CDS]\n"
            "ATGAAATTTATGAAATTT\n"
            ">lcl|CM075669.1_cds_2 [location=join(301..309,401..409)] [gbkey=CDS]\n"
            "ATGCCCTTTATGCCCTTT\n"
        )
    with gzip.open(gff_path, "wt", encoding="utf-8") as handle:
        handle.write(
            "\n".join(
                [
                    "##gff-version 3",
                    "CM075669.1\tGenbank\tgene\t101\t209\t.\t-\t.\tID=gene-RJT34_00001;Name=RJT34_00001;locus_tag=RJT34_00001;pseudo=true",
                    "CM075669.1\tGenbank\tmRNA\t101\t209\t.\t-\t.\tID=rna-RJT34_mrna00001;Parent=gene-RJT34_00001;locus_tag=RJT34_00001;pseudo=true",
                    "CM075669.1\tGenbank\tCDS\t201\t209\t.\t-\t0\tID=cds-RJT34_00001;Parent=rna-RJT34_mrna00001;Dbxref=InterPro:IPR000504,InterPro:IPR035979;locus_tag=RJT34_00001;pseudo=true",
                    "CM075669.1\tGenbank\tCDS\t101\t109\t.\t-\t0\tID=cds-RJT34_00001;Parent=rna-RJT34_mrna00001;Dbxref=InterPro:IPR000504,InterPro:IPR035979;locus_tag=RJT34_00001;pseudo=true",
                    "CM075669.1\tGenbank\tgene\t301\t409\t.\t+\t.\tID=gene-RJT34_00002;Name=RJT34_00002;locus_tag=RJT34_00002",
                    "CM075669.1\tGenbank\tmRNA\t301\t409\t.\t+\t.\tID=rna-RJT34_mrna00002;Parent=gene-RJT34_00002;locus_tag=RJT34_00002",
                    "CM075669.1\tGenbank\tCDS\t301\t309\t.\t+\t0\tID=cds-RJT34_00002;Parent=rna-RJT34_mrna00002;locus_tag=RJT34_00002;pseudo=true",
                    "CM075669.1\tGenbank\tCDS\t401\t409\t.\t+\t0\tID=cds-RJT34_00002;Parent=rna-RJT34_mrna00002;locus_tag=RJT34_00002;pseudo=true",
                    "",
                ]
            )
        )
    with gzip.open(genome_path, "wt", encoding="utf-8") as handle:
        handle.write(">CM075669.1\n" + ("A" * 500) + "\n")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        provider,
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
        "--strict",
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = next(out_cds.glob("*.fa.gz"))
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        headers = [line.strip() for line in handle if line.startswith(">")]
    assert headers == [
        ">Clitoria_ternatea_RJT34_00001",
        ">Clitoria_ternatea_RJT34_00002",
    ]

    with open(str(formatted_cds) + ".gff-grouping.tsv", "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["mapping_status"] for row in rows] == ["mapped", "mapped"]
    assert [row["matched_aliases"] for row in rows] == ["location", "location"]
    assert all("InterPro" not in row["matched_aliases"] for row in rows)

    mapping = run_validate_mapping_script(
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--strict",
    )
    assert mapping.returncode == 0, mapping.stderr + "\n" + mapping.stdout
    assert "CDS-to-GFF mapping OK: 2/2 IDs" in mapping.stdout


def test_direct_ncbi_pseudogene_location_match_requires_unique_exact_coordinates(tmp_path):
    module = load_module()
    gff_path = tmp_path / "pseudogenes.gff"
    gff_path.write_text(
        "\n".join(
            [
                "CM000001.1\tGenbank\tpseudogene\t1\t9\t.\t+\t.\tID=gene-L1;locus_tag=L1;pseudo=true",
                "CM000001.1\tGenbank\tmRNA\t1\t9\t.\t+\t.\tID=rna-L1;Parent=gene-L1;locus_tag=L1;pseudo=true",
                "CM000001.1\tGenbank\tCDS\t1\t9\t.\t+\t0\tID=cds-L1;Parent=rna-L1;pseudo=true",
                "CM000001.1\tGenbank\tpseudogene\t20\t28\t.\t-\t.\tID=gene-L2;locus_tag=L2;pseudo=true",
                "CM000001.1\tGenbank\tmRNA\t20\t28\t.\t-\t.\tID=rna-L2;Parent=gene-L2;locus_tag=L2;pseudo=true",
                "CM000001.1\tGenbank\tCDS\t20\t28\t.\t-\t0\tID=cds-L2;Parent=rna-L2;pseudo=true",
                "",
            ]
        ),
        encoding="utf-8",
    )
    cds_path = tmp_path / "pseudogenes_cds.fna"
    cds_path.write_text(
        ">lcl|CM000001.1_cds_1 [location=1..9] [gbkey=CDS]\nATGAAATTT\n",
        encoding="utf-8",
    )
    task = {
        "provider": "direct",
        "species_key": "Example_species",
        "species_prefix": "Example_species",
        "cds_path": cds_path,
        "gff_path": gff_path,
    }
    index = module.build_gff_cds_grouping_index(task)
    mapped = module.resolve_cds_header_gff_gene(
        task, "lcl|CM000001.1_cds_1 [location=1..9] [gbkey=CDS]", index
    )
    assert (mapped["status"], mapped["gene_token"], mapped["matched_aliases"]) == (
        "mapped", "L1", ("location",)
    )
    assert module.resolve_cds_header_gff_gene(
        {**task, "provider": "ensembl"},
        "lcl|CM000001.1_cds_1 [location=1..9]",
        index,
    )["status"] == "unmapped"
    output_dir = tmp_path / "formatted"
    output_dir.mkdir()
    formatted = module.format_cds(task, output_dir, overwrite=False, dry_run=False, strict=True)
    with gzip.open(formatted["output_path"], "rt", encoding="utf-8") as handle:
        assert handle.read() == ">Example_species_L1\nATGAAATTT\n"
    for header in (
        "lcl|CM000001.1_cds_2 [location=complement(1..9)] [gbkey=CDS]",
        "lcl|CM000001.1_cds_3 [location=1..10] [gbkey=CDS]",
        "custom_cds_4 [location=1..9]",
    ):
        assert module.resolve_cds_header_gff_gene(task, header, index)["status"] == "unmapped"
        assert not module.is_unlinkable_anonymous_ncbi_cds(task, header, "unmapped")

    with gff_path.open("a", encoding="utf-8") as handle:
        handle.write(
            "CM000001.1\tGenbank\tpseudogene\t1\t9\t.\t+\t.\tID=gene-L3;locus_tag=L3;pseudo=true\n"
            "CM000001.1\tGenbank\tmRNA\t1\t9\t.\t+\t.\tID=rna-L3;Parent=gene-L3;locus_tag=L3;pseudo=true\n"
            "CM000001.1\tGenbank\tCDS\t1\t9\t.\t+\t0\tID=cds-L3;Parent=rna-L3;pseudo=true\n"
        )
    ambiguous_index = module.build_gff_cds_grouping_index(task)
    ambiguous = module.resolve_cds_header_gff_gene(
        task, "lcl|CM000001.1_cds_1 [location=1..9] [gbkey=CDS]", ambiguous_index
    )
    assert ambiguous["status"] == "ambiguous"
    assert ambiguous["candidate_gene_tokens"] == ("L1", "L3")


def test_ncbi_conflicting_header_locus_requires_matching_protein_and_exact_location(tmp_path):
    module = load_module()
    gff_path = tmp_path / "models.gff3"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tGenbank\tgene\t1\t9\t.\t+\t.\tID=gene-L1;locus_tag=L1",
                "chr1\tGenbank\tmRNA\t1\t9\t.\t+\t.\tID=rna-L1;Parent=gene-L1;locus_tag=L1",
                "chr1\tGenbank\tCDS\t1\t9\t.\t+\t0\tID=cds-P1;Parent=rna-L1;protein_id=P1;locus_tag=L1",
                "chr1\tGenbank\tgene\t20\t28\t.\t+\t.\tID=gene-L2;locus_tag=L2",
                "chr1\tGenbank\tmRNA\t20\t28\t.\t+\t.\tID=rna-L2;Parent=gene-L2;locus_tag=L2",
                "chr1\tGenbank\tCDS\t20\t28\t.\t+\t0\tID=cds-P2;Parent=rna-L2;protein_id=P2;locus_tag=L2",
                "",
            ]
        ),
        encoding="utf-8",
    )
    cds_path = tmp_path / "models.fna"
    header = "lcl|chr1_cds_P1_1 [locus_tag=L2] [protein_id=P1] [location=1..9]"
    cds_path.write_text(">" + header + "\nATGAAATTT\n", encoding="utf-8")
    task = {
        "provider": "direct",
        "species_key": "Example_species",
        "species_prefix": "Example_species",
        "cds_path": cds_path,
        "gff_path": gff_path,
        "gene_grouping_mode": "strict",
    }
    index = module.build_gff_cds_grouping_index(task)
    mapped = module.resolve_cds_header_gff_gene(task, header, index)
    assert mapped["status"] == "mapped"
    assert mapped["gene_token"] == "L1"
    assert mapped["ignored_conflicting_locus_tag"] == "L2"
    anonymous = "lcl|chr1_cds_L2_1 [locus_tag=L2] [location=1..9] [gbkey=CDS]"
    assert module.resolve_cds_header_gff_gene(task, anonymous, index)["gene_token"] == "L1"
    output_dir = tmp_path / "formatted"
    output_dir.mkdir()
    result = module.format_cds(task, output_dir, overwrite=False, dry_run=False, strict=True)
    with gzip.open(result["output_path"], "rt", encoding="utf-8") as handle:
        assert handle.read() == ">Example_species_L1\nATGAAATTT\n"
    with open(result["gff_grouping_audit_path"], newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["ignored_conflicting_locus_tag"] == "L2"

    for rejected_header in (
        header.replace("[protein_id=P1]", "[protein_id=P2]"),
        header.replace("[location=1..9]", "[location=2..9]"),
        header + " [transcript_id=rna-L2]",
        anonymous.replace(" [gbkey=CDS]", ""),
    ):
        assert module.resolve_cds_header_gff_gene(task, rejected_header, index)["status"] != "mapped"

    gff_path.write_text(
        gff_path.read_text(encoding="utf-8").replace("protein_id=P2", "protein_id=P1"),
        encoding="utf-8",
    )
    duplicate_protein_index = module.build_gff_cds_grouping_index(task)
    assert module.resolve_cds_header_gff_gene(task, header, duplicate_protein_index)["status"] == "ambiguous"


def test_format_species_inputs_does_not_exclude_named_ncbi_cds_from_unrelated_gff(tmp_path):
    module = load_module()
    task = {
        "provider": "ncbi",
        "species_prefix": "Example_species",
    }
    header = (
        "lcl|NC_000001.1_cds_XP_123.1_1 [protein_id=XP_123.1] "
        "[location=1..6] [gbkey=CDS]"
    )

    assert not module.is_unlinkable_anonymous_ncbi_cds(task, header, "unmapped")


def test_format_species_inputs_maps_protein_id_through_gff_gene_hierarchy(tmp_path):
    input_dir = tmp_path / "NCBI_Genome" / "species_wise_original"
    species_dir = input_dir / "Dictyostelium_firmibasis"
    species_dir.mkdir(parents=True, exist_ok=True)
    cds_path = species_dir / "GCA_036169595.1_ASM3616959v1_cds_from_genomic.fna.gz"
    gff_path = species_dir / "GCA_036169595.1_ASM3616959v1_genomic.gff.gz"
    genome_path = species_dir / "GCA_036169595.1_ASM3616959v1_genomic.fna.gz"
    with gzip.open(cds_path, "wt", encoding="utf-8") as handle:
        handle.write(
            (
                ">lcl|CM069765.1_cds_KAK5581746.1_1 [protein=hypothetical protein] [protein_id=KAK5581746.1] [location=join(5022..5093,5192..5424)] [gbkey=CDS]\n"
                "ATGCAAACAAATACATTTAGCAATGTACCTGGCTCACTTAATATTGAAGACCTATTAAATAAAATAGAAACTGTAGTATT\n"
                ">lcl|CM069765.1_cds_KAK5581747.1_2 [protein=hypothetical protein] [protein_id=KAK5581747.1] [location=8575..9000] [gbkey=CDS]\n"
                "ATGGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAAATCGCCGAA\n"
            )
        )
    with gzip.open(gff_path, "wt", encoding="utf-8") as handle:
        handle.write(
            "\n".join(
                [
                    "##gff-version 3",
                    "CM069765.1\tGenbank\tgene\t5010\t5424\t.\t+\t.\tID=gene-RB653_003324;Name=RB653_003324;gbkey=Gene;gene_biotype=protein_coding;locus_tag=RB653_003324",
                    "CM069765.1\tGenbank\tmRNA\t5010\t5424\t.\t+\t.\tID=rna-DFI_00002-RA;Parent=gene-RB653_003324;gbkey=mRNA;locus_tag=RB653_003324;orig_protein_id=gnl|WGS:JAVFKY|DFI_00002-RA:cds;orig_transcript_id=gnl|WGS:JAVFKY|DFI_00002-RA;product=hypothetical protein",
                    "CM069765.1\tGenbank\tCDS\t5022\t5093\t.\t+\t0\tID=cds-KAK5581746.1;Parent=rna-DFI_00002-RA;Dbxref=NCBI_GP:KAK5581746.1;Name=KAK5581746.1;gbkey=CDS;locus_tag=RB653_003324;orig_transcript_id=gnl|WGS:JAVFKY|DFI_00002-RA;product=hypothetical protein;protein_id=KAK5581746.1",
                    "CM069765.1\tGenbank\tCDS\t5192\t5424\t.\t+\t0\tID=cds-KAK5581746.1;Parent=rna-DFI_00002-RA;Dbxref=NCBI_GP:KAK5581746.1;Name=KAK5581746.1;gbkey=CDS;locus_tag=RB653_003324;orig_transcript_id=gnl|WGS:JAVFKY|DFI_00002-RA;product=hypothetical protein;protein_id=KAK5581746.1",
                    "CM069765.1\tGenbank\tgene\t8575\t9000\t.\t+\t.\tID=gene-RB653_003325;Name=RB653_003325;gbkey=Gene;gene_biotype=protein_coding;locus_tag=RB653_003325",
                    "CM069765.1\tGenbank\tmRNA\t8575\t9000\t.\t+\t.\tID=rna-DFI_00003-RA;Parent=gene-RB653_003325;gbkey=mRNA;locus_tag=RB653_003325;orig_protein_id=gnl|WGS:JAVFKY|DFI_00003-RA:cds;orig_transcript_id=gnl|WGS:JAVFKY|DFI_00003-RA;product=hypothetical protein",
                    "CM069765.1\tGenbank\tCDS\t8575\t9000\t.\t+\t0\tID=cds-KAK5581747.1;Parent=rna-DFI_00003-RA;Dbxref=NCBI_GP:KAK5581747.1;Name=KAK5581747.1;gbkey=CDS;locus_tag=RB653_003325;orig_transcript_id=gnl|WGS:JAVFKY|DFI_00003-RA;product=hypothetical protein;protein_id=KAK5581747.1",
                    "",
                ]
            )
        )
    with gzip.open(genome_path, "wt", encoding="utf-8") as handle:
        handle.write(">CM069765.1\n" + ("A" * 12000) + "\n")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "ncbi",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Dictyostelium_firmibasis_GCA_036169595.1_ASM3616959v1_cds_from_genomic.fa.gz"
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        headers = [line.strip() for line in handle if line.startswith(">")]
    assert headers == [
        ">Dictyostelium_firmibasis_RB653_003324",
        ">Dictyostelium_firmibasis_RB653_003325",
    ]

    mapping = run_validate_mapping_script(
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
    )
    assert mapping.returncode == 0, mapping.stderr + "\n" + mapping.stdout
    assert "[Dictyostelium_firmibasis] CDS-to-GFF mapping OK: 2/2 IDs" in mapping.stdout


def test_format_species_inputs_strict_gene_grouping_keeps_misassigned_gene_ids_separate(tmp_path):
    input_dir = tmp_path / "Direct" / "species_wise_original"
    species_dir = input_dir / "Arabidopsis_thaliana"
    species_dir.mkdir(parents=True, exist_ok=True)
    gff_path = species_dir / "Arabidopsis_thaliana.annotation.gff3"
    genome_path = species_dir / "Arabidopsis_thaliana.genome.fa"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\tsrc\tgene\t1\t18\t.\t+\t.\tID=badGeneA",
                "chr1\tsrc\tmRNA\t1\t18\t.\t+\t.\tID=locusX.t1;Parent=badGeneA",
                "chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tID=cds1;Parent=locusX.t1",
                "chr1\tsrc\tCDS\t13\t18\t.\t+\t0\tID=cds2;Parent=locusX.t1",
                "chr1\tsrc\tgene\t1\t18\t.\t+\t.\tID=badGeneB",
                "chr1\tsrc\tmRNA\t1\t18\t.\t+\t.\tID=locusX.t2;Parent=badGeneB",
                "chr1\tsrc\tCDS\t1\t3\t.\t+\t0\tID=cds3;Parent=locusX.t2",
                "chr1\tsrc\tCDS\t13\t18\t.\t+\t0\tID=cds4;Parent=locusX.t2",
                "",
            ]
        ),
        encoding="utf-8",
    )
    genome_path.write_text(">chr1\nATGAAACCCGGGTTTCCC\n", encoding="utf-8")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "direct",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
        "--gene-grouping-mode",
        "strict",
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Arabidopsis_thaliana_annotation.derived.cds.fa.gz"
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        headers = [line.strip() for line in handle if line.startswith(">")]
    assert headers == [
        ">Arabidopsis_thaliana_badGeneA",
        ">Arabidopsis_thaliana_badGeneB",
    ]


def test_format_species_inputs_fernbase_prefers_primary_annotation_files(tmp_path):
    input_dir = tmp_path / "FernBase" / "species_wise_original"
    species_dir = input_dir / "Azolla_filiculoides"
    species_dir.mkdir(parents=True, exist_ok=True)
    (species_dir / "Azolla_filiculoides.CDS.lowconfidence_v1.1.fasta").write_text(
        ">Azfi_g1.t1\nATG\n", encoding="utf-8"
    )
    (species_dir / "Azolla_filiculoides.CDS.highconfidence_v1.1.fasta").write_text(
        ">Azfi_g1.t1\nATGAA\n", encoding="utf-8"
    )
    (species_dir / "Azolla_filiculoides.transcript.highconfidence_v1.1.fasta").write_text(
        ">Azfi_g1.t1\nATGAAA\n", encoding="utf-8"
    )
    (species_dir / "Azolla_filiculoides.gene_models.lowconfidence_v1.1.gff").write_text(
        "chr1\tsrc\tgene\t1\t3\t.\t+\t.\tID=Azfi_g1\n",
        encoding="utf-8",
    )
    (species_dir / "Azolla_filiculoides.gene_models.highconfidence_v1.1.gff").write_text(
        "chr1\tsrc\tgene\t1\t5\t.\t+\t.\tID=Azfi_g1\n",
        encoding="utf-8",
    )
    (species_dir / "Azolla_filiculoides.genome_v1.2.fasta").write_text(">chr1\nATGCATGC\n", encoding="utf-8")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "fernbase",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Azolla_filiculoides_CDS.highconfidence_v1.1.fa.gz"
    formatted_gff = out_gff / "Azolla_filiculoides_gene_models.highconfidence_v1.1.gff.gz"
    formatted_genome = out_genome / "Azolla_filiculoides_genome_v1.2.fa.gz"
    assert formatted_cds.exists()
    assert formatted_gff.exists()
    assert formatted_genome.exists()

    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        cds_text = handle.read()
    assert ">Azolla_filiculoides_Azfi_g1" in cds_text
    assert "ATGAAN" in cds_text
    assert "lowconfidence" not in formatted_gff.name


def test_format_species_inputs_fernbase_accepts_markerless_genome_fasta(tmp_path):
    input_dir = tmp_path / "FernBase" / "species_wise_original"
    species_dir = input_dir / "Ceratopteris_richardii"
    species_dir.mkdir(parents=True, exist_ok=True)
    (species_dir / "Crichardii_676_v2.0_cds.fa").write_text(">Crich_g1.t1\nATGAA\n", encoding="utf-8")
    (species_dir / "Crichardii_676_v2.1.gene.gff3").write_text(
        "chr1\tsrc\tgene\t1\t5\t.\t+\t.\tID=Crich_g1\n",
        encoding="utf-8",
    )
    (species_dir / "Crichardii_676_v2.0.fa").write_text(">chr1\nATGCATGC\n", encoding="utf-8")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "fernbase",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_genome = out_genome / "Ceratopteris_richardii_Crichardii_676_v2.0.fa.gz"
    assert formatted_genome.exists()


def test_format_species_inputs_fernbase_uses_plain_gene_tag_for_aggregation(tmp_path):
    input_dir = tmp_path / "FernBase" / "species_wise_original"
    species_dir = input_dir / "Pteridium_aquilinum"
    species_dir.mkdir(parents=True, exist_ok=True)
    (species_dir / "pt_aq.cds.fa").write_text(
        (
            ">pteridium_mrna-253 gene=pteridium_gene-252 seq_id=scaff1 type=cds\n"
            "ATGAA\n"
            ">pteridium_mrna-303 gene=pteridium_gene-295 seq_id=scaff1 type=cds\n"
            "ATGAAA\n"
            ">pteridium_mrna-308 gene=pteridium_gene-295 seq_id=scaff1 type=cds\n"
            "ATGAAATAG\n"
        ),
        encoding="utf-8",
    )
    (species_dir / "pt_aq.gff3").write_text(
        (
            "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=pteridium_gene-252\n"
            "chr1\tsrc\tgene\t20\t40\t.\t+\t.\tID=pteridium_gene-295\n"
        ),
        encoding="utf-8",
    )
    (species_dir / "pt_aq_final_genome.fasta").write_text(">chr1\nATGCATGCATGC\n", encoding="utf-8")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    species_summary = tmp_path / "gg_input_generation_species.tsv"
    completed = run_script(
        "--provider",
        "fernbase",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
        "--species-summary-output",
        str(species_summary),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Pteridium_aquilinum_pt_aq.cds.fa.gz"
    assert formatted_cds.exists()
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        headers = [line.strip() for line in handle if line.startswith(">")]
    assert headers == [
        ">Pteridium_aquilinum_pteridium_gene-252",
        ">Pteridium_aquilinum_pteridium_gene-295",
    ]

    with open(species_summary, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["species_prefix"] == "Pteridium_aquilinum"
    assert rows[0]["aggregated_cds_removed"] == "1"
    assert rows[0]["cds_sequences_before"] == "3"
    assert rows[0]["cds_sequences_after"] == "2"


def test_reuse_gff_grouped_cds_without_audit_reports_unverified_grouping(tmp_path):
    mod = load_module()
    cds_path = tmp_path / "input.cds.fa"
    gff_path = tmp_path / "input.gff3"
    cds_path.write_text(">transcript1\nATG\n>transcript2\nATGA\n", encoding="utf-8")
    gff_path.write_text("chr1\tsrc\tgene\t1\t4\t.\t+\t.\tID=gene1\n", encoding="utf-8")
    task = {
        "provider": "local",
        "species_key": "Species_a",
        "species_prefix": "Species_a",
        "cds_path": cds_path,
        "gff_path": gff_path,
    }
    output_dir = tmp_path / "output"
    output_dir.mkdir()
    output_path = output_dir / mod.normalize_cds_output_basename(cds_path.name, "Species_a")
    with gzip.open(output_path, "wt", encoding="utf-8") as handle:
        handle.write(">Species_a_gene1\nATG\n")

    result = mod.format_cds(task, output_dir, overwrite=False, dry_run=False, reuse_existing=True)

    assert result["status"] == "skip"
    assert result["output_path"] == output_path
    assert result["after_count"] == 1
    assert result["first_sequence_name"] == "Species_a_gene1"
    assert result["grouping_source"] == "reused_without_audit"
    assert result["gff_grouping_audit_path"] == ""
    assert not Path(str(output_path) + ".gff-grouping.json").exists()


def test_write_fasta_records_gzip_prefers_seqkit(monkeypatch, tmp_path):
    mod = load_module()
    output_path = tmp_path / "species.fa.gz"
    calls = {}

    class FakeSeqkitProcess:
        def __init__(self, command, **kwargs):
            calls["command"] = command
            calls["kwargs"] = kwargs
            self.command = command
            self._stdin_pipe = FakeTextPipe()
            self.stdin = self._stdin_pipe
            self.stderr = io.StringIO("")

        def wait(self):
            output_arg = self.command[self.command.index("-o") + 1]
            with gzip.open(output_arg, "wt", encoding="utf-8") as handle:
                handle.write(self._stdin_pipe.getvalue())
            return 0

        def kill(self):
            return None

    monkeypatch.setenv("GG_TASK_CPUS", "3")
    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/{}".format(name))
    monkeypatch.setattr(subprocess, "Popen", FakeSeqkitProcess)

    mod.write_fasta_records_gzip(output_path, [("seq1", "ATG"), ("seq2", "ATGA")])

    with gzip.open(output_path, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert ">seq1" in text
    assert ">seq2" in text
    assert calls["command"][:4] == ["/usr/bin/seqkit", "seq", "--threads", "3"]
    assert calls["command"][-1] == "-"


def test_fallback_fasta_writer_preserves_existing_output_on_record_error(monkeypatch, tmp_path):
    mod = load_module()
    output_path = tmp_path / "species.fa.gz"
    with gzip.open(output_path, "wt", encoding="utf-8") as handle:
        handle.write(">original\nACGT\n")
    monkeypatch.setattr(shutil, "which", lambda _name: None)

    def broken_records():
        yield "new", "ATG"
        raise RuntimeError("source read failed")

    with pytest.raises(RuntimeError, match="source read failed"):
        mod.write_fasta_records_gzip(output_path, broken_records())

    with gzip.open(output_path, "rt", encoding="utf-8") as handle:
        assert handle.read() == ">original\nACGT\n"
    assert not list(tmp_path.glob(".species.fa.tmp.*"))


def test_write_gff_gzip_prefers_pigz(monkeypatch, tmp_path):
    mod = load_module()
    input_path = tmp_path / "input.gff3"
    input_path.write_text("chr1\tsrc\tgene\t1\t3\t.\t+\t.\tID=evm.model.g1\n", encoding="utf-8")
    output_path = tmp_path / "output.gff.gz"
    calls = {}

    class FakePigzProcess:
        def __init__(self, command, **kwargs):
            calls["command"] = command
            calls["kwargs"] = kwargs
            self.stdout = kwargs["stdout"]
            self._stdin_pipe = FakeTextPipe()
            self.stdin = self._stdin_pipe
            self.stderr = io.StringIO("")

        def wait(self):
            with gzip.GzipFile(fileobj=self.stdout, mode="wb") as handle:
                handle.write(self._stdin_pipe.getvalue().encode("utf-8"))
            return 0

        def kill(self):
            return None

    monkeypatch.setenv("GG_TASK_CPUS", "4")
    monkeypatch.setattr(shutil, "which", lambda name: "/usr/bin/{}".format(name))
    monkeypatch.setattr(subprocess, "Popen", FakePigzProcess)

    line_count = mod.write_gff_gzip(input_path, output_path)

    with gzip.open(output_path, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert line_count == 1
    assert "evm.model." not in text
    assert calls["command"] == ["/usr/bin/pigz", "-p", "4", "-c"]


def test_fallback_gff_writers_preserve_existing_output_on_error(monkeypatch, tmp_path):
    mod = load_module()
    monkeypatch.setattr(shutil, "which", lambda _name: None)
    output_path = tmp_path / "species.gff.gz"
    with gzip.open(output_path, "wt", encoding="utf-8") as handle:
        handle.write("original\n")

    def broken_lines():
        yield "chr1\tsrc\tgene\t1\t3\t.\t+\t.\tID=new\n"
        raise RuntimeError("GFF generation failed")

    with pytest.raises(RuntimeError, match="GFF generation failed"):
        mod.write_gff_lines_gzip(output_path, broken_lines())
    with gzip.open(output_path, "rt", encoding="utf-8") as handle:
        assert handle.read() == "original\n"

    input_path = tmp_path / "truncated.gff.gz"
    input_path.write_bytes(gzip.compress(b"chr1\tsrc\tgene\t1\t3\t.\t+\t.\tID=new\n")[:-4])
    with pytest.raises((EOFError, OSError)):
        mod.write_gff_gzip(input_path, output_path)
    with gzip.open(output_path, "rt", encoding="utf-8") as handle:
        assert handle.read() == "original\n"
    assert not list(tmp_path.glob(".species.gff.tmp.*"))


def test_resolve_provider_download_limits_keeps_fernbase_and_insectbase_default_caps_at_two(monkeypatch):
    mod = load_module()
    monkeypatch.delenv("GG_INPUT_MAX_CONCURRENT_DOWNLOADS_FERNBASE", raising=False)
    monkeypatch.delenv("GG_INPUT_MAX_CONCURRENT_DOWNLOADS_INSECTBASE", raising=False)
    limits = mod.resolve_provider_download_limits(8)
    assert limits["fernbase"] == 2
    assert limits["insectbase"] == 2


def test_resolve_ncbi_download_urls_from_id_retries_transient_remote_disconnect(monkeypatch):
    mod = load_module()
    calls = {"esearch": 0, "esummary": 0}

    def fake_urlopen(request, timeout):
        url = request.full_url
        if "esearch.fcgi" in url:
            calls["esearch"] += 1
            if calls["esearch"] == 1:
                raise RemoteDisconnected("Remote end closed connection without response")
            payload = {
                "header": {"type": "esearch", "version": "0.3"},
                "esearchresult": {"idlist": ["12345"]},
            }
            return FakeBinaryResponse(json.dumps(payload).encode("utf-8"))
        if "esummary.fcgi" in url:
            calls["esummary"] += 1
            payload = {
                "header": {"type": "esummary", "version": "0.3"},
                "result": {
                    "uids": ["12345"],
                    "12345": {
                        "ftppath_genbank": (
                            "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/036/169/595/GCA_036169595.1_ASM3616959v1"
                        ),
                        "organism": "Dictyostelium firmibasis",
                        "speciesname": "Dictyostelium firmibasis",
                    },
                },
            }
            return FakeBinaryResponse(json.dumps(payload).encode("utf-8"))
        raise AssertionError(url)

    provider_module = sys.modules[mod.resolve_ncbi_download_urls_from_id.__module__]
    monkeypatch.setattr(provider_module, "urlopen", fake_urlopen)
    monkeypatch.setattr(provider_module, "throttle_ncbi_eutils_request", lambda: None)
    monkeypatch.setattr(provider_module.time, "sleep", lambda _seconds: None)

    resolved = mod.resolve_ncbi_download_urls_from_id("GCA_036169595.1", timeout=1.0)
    assert calls["esearch"] == 2
    assert calls["esummary"] == 1
    assert resolved["species_key"] == "Dictyostelium_firmibasis"
    assert resolved["cds_filename"] == "GCA_036169595.1_ASM3616959v1_cds_from_genomic.fna.gz"
    assert resolved["ncbi_source_db"] == "genbank"


def test_resolve_ncbi_download_urls_from_id_uses_datasets_when_ftp_path_is_missing(monkeypatch):
    mod = load_module()

    def fake_urlopen(request, timeout):
        url = request.full_url
        if "esearch.fcgi" in url:
            payload = {
                "header": {"type": "esearch", "version": "0.3"},
                "esearchresult": {"idlist": ["34905841"]},
            }
            return FakeBinaryResponse(json.dumps(payload).encode("utf-8"))
        if "esummary.fcgi" in url:
            payload = {
                "header": {"type": "esummary", "version": "0.3"},
                "result": {
                    "uids": ["34905841"],
                    "34905841": {
                        "assemblyaccession": "GCA_059696495.1",
                        "assemblyname": "Tcin_v1.0",
                        "organism": "Tanacetum cinerariifolium (pyrethrum)",
                        "speciesname": "Tanacetum cinerariifolium",
                        "ftppath_refseq": "",
                        "ftppath_genbank": "",
                    },
                },
            }
            return FakeBinaryResponse(json.dumps(payload).encode("utf-8"))
        raise AssertionError(url)

    provider_module = sys.modules[mod.resolve_ncbi_download_urls_from_id.__module__]
    monkeypatch.setattr(provider_module, "urlopen", fake_urlopen)
    monkeypatch.setattr(provider_module, "throttle_ncbi_eutils_request", lambda: None)

    resolved = mod.resolve_ncbi_download_urls_from_id("GCA_059696495.1", timeout=1.0)

    assert resolved["species_key"] == "Tanacetum_cinerariifolium"
    assert resolved["gbff_url"] == "ncbi-datasets://GCA_059696495.1/gbff"
    assert resolved["genome_url"] == "ncbi-datasets://GCA_059696495.1/genome"
    assert resolved["gbff_filename"] == "GCA_059696495.1_Tcin_v1.0_genomic.gbff.gz"
    assert resolved["genome_filename"] == "GCA_059696495.1_Tcin_v1.0_genomic.fna.gz"
    assert resolved["ncbi_source_db"] == "datasets_api"


def test_iter_fasta_records_reads_tar_bz2_archive(tmp_path):
    mod = load_module()
    archive_path = tmp_path / "example.genome.fa.tar.bz2"
    payload = ">chr1 description\nATGC\nATGC\n".encode("utf-8")
    with tarfile.open(archive_path, "w:bz2") as archive:
        info = tarfile.TarInfo(name="nested/example.genome.fa")
        info.size = len(payload)
        archive.addfile(info, io.BytesIO(payload))
        note = tarfile.TarInfo(name="README.txt")
        note_payload = b"fixture\n"
        note.size = len(note_payload)
        archive.addfile(note, io.BytesIO(note_payload))

    records = list(mod.iter_fasta_records(archive_path))
    assert records == [("chr1 description", "ATGCATGC")]


def test_format_species_inputs_fernbase_prefers_namespaced_transcript_gene_id_when_gene_tag_is_short(tmp_path):
    input_dir = tmp_path / "FernBase" / "species_wise_original"
    species_dir = input_dir / "Salvinia_molesta"
    species_dir.mkdir(parents=True, exist_ok=True)
    (species_dir / "sg2.cds").write_text(
        (">sg2.g13251.t1 gene=g13251\nATGAAATAG\n>sg2.g9.t1 gene=g9\nATGAAATAA\n"),
        encoding="utf-8",
    )
    (species_dir / "sg2.gff3").write_text(
        ("Chr_1\tgmst\tgene\t1\t9\t.\t+\t.\tID=sg2.g13251\nChr_1\tgmst\tgene\t20\t28\t.\t+\t.\tID=sg2.g9\n"),
        encoding="utf-8",
    )
    (species_dir / "sg2_genome.fasta").write_text(">Chr_1\nATGCATGCATGC\n", encoding="utf-8")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "fernbase",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Salvinia_molesta_sg2.cds.fa.gz"
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        headers = [line.strip() for line in handle if line.startswith(">")]
    assert headers == [
        ">Salvinia_molesta_sg2.g13251",
        ">Salvinia_molesta_sg2.g9",
    ]


def test_format_species_inputs_fernbase_strips_amt_suffix_when_gff_uses_base_gene_id(tmp_path):
    input_dir = tmp_path / "FernBase" / "species_wise_original"
    species_dir = input_dir / "Azolla_filiculoides"
    species_dir.mkdir(parents=True, exist_ok=True)
    (species_dir / "Azolla_filiculoides.CDS.highconfidence_v1.1.fasta").write_text(
        (">Azfi_s0034.g025227.AMT2\nATGAAATAG\n>Azfi_s0093.g043301.AMT2\nATGAAATAA\n"),
        encoding="utf-8",
    )
    (species_dir / "Azolla_filiculoides.gene_models.highconfidence_v1.1.gff").write_text(
        (
            "SCAF_1\tAUGUSTUS\tgene\t1\t9\t.\t+\t.\tID=Azfi_s0034.g025227\n"
            "SCAF_1\tAUGUSTUS\tgene\t20\t28\t.\t+\t.\tID=Azfi_s0093.g043301\n"
        ),
        encoding="utf-8",
    )
    (species_dir / "Azolla_filiculoides.genome_v1.2.fasta").write_text(">SCAF_1\nATGCATGCATGC\n", encoding="utf-8")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "fernbase",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Azolla_filiculoides_CDS.highconfidence_v1.1.fa.gz"
    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        headers = [line.strip() for line in handle if line.startswith(">")]
    assert headers == [
        ">Azolla_filiculoides_Azfi_s0034.g025227",
        ">Azolla_filiculoides_Azfi_s0093.g043301",
    ]


def test_generic_discovery_prefers_cds_and_gene_model_gff_candidates(tmp_path):
    module = load_module()
    input_dir = tmp_path / "DDBJ" / "species_wise_original"
    species_dir = input_dir / "Test_species"
    species_dir.mkdir(parents=True)
    (species_dir / "a.cdna.fa").write_text(">T1\nATG\n", encoding="utf-8")
    cds_path = species_dir / "z.cds.fa"
    cds_path.write_text(">T1\nATG\n", encoding="utf-8")
    (species_dir / "a.repeat.gff3").write_text(
        "chr1\trepeat\trepeat_region\t1\t3\t.\t+\t.\tID=R1\n", encoding="utf-8"
    )
    gff_path = species_dir / "z.genes.gff3"
    gff_path.write_text("chr1\tsrc\tgene\t1\t3\t.\t+\t.\tID=T1\n", encoding="utf-8")

    tasks, _warnings, errors = module.discover_tasks("ddbj", input_dir)

    assert errors == []
    assert tasks[0]["cds_path"] == cds_path
    assert tasks[0]["gff_path"] == gff_path


def test_plantgarden_discovery_accepts_provider_native_filenames(tmp_path):
    module = load_module()
    input_dir = tmp_path / "PlantGARDEN" / "species_wise_original"
    species_dir = input_dir / "Actinidia_polygama"
    species_dir.mkdir(parents=True)
    cds_path = species_dir / "APO1.1.cds.fasta.gz"
    cds_path.write_text(">G1\nATG\n", encoding="utf-8")
    gff_path = species_dir / "APO1.1.genes.gff.gz"
    gff_path.write_text("chr1\tsrc\tgene\t1\t3\t.\t+\t.\tID=G1\n", encoding="utf-8")
    genome_path = species_dir / "APO_r1.1.pmol.fasta.gz"
    genome_path.write_text(">chr1\nATG\n", encoding="utf-8")

    tasks, _warnings, errors = module.discover_tasks("plantgarden", input_dir)

    assert errors == []
    assert len(tasks) == 1
    assert tasks[0]["cds_path"] == cds_path
    assert tasks[0]["gff_path"] == gff_path
    assert tasks[0]["genome_path"] == genome_path


def test_phycocosm_discovery_never_uses_assembly_as_cds(tmp_path):
    module = load_module()
    input_dir = tmp_path / "PhycoCosm" / "species_wise_original"
    species_dir = input_dir / "Microglena_spYARC"
    species_dir.mkdir(parents=True)
    cds_path = species_dir / "MicrYARC1_GeneCatalog_CDS.fasta"
    cds_path.write_text(">T1\nATG\n", encoding="utf-8")
    genome_path = species_dir / "MicrYARC1_genome_assembly.fasta"
    genome_path.write_text(">chr1\nATG\n", encoding="utf-8")
    (species_dir / "MicrYARC1_genes.gff3").write_text(
        "chr1\tsrc\tgene\t1\t3\t.\t+\t.\tID=T1\n", encoding="utf-8"
    )

    tasks, _warnings, errors = module.discover_tasks("phycocosm", input_dir)

    assert errors == []
    assert tasks[0]["cds_path"] == cds_path
    assert tasks[0]["genome_path"] == genome_path


def test_ncbi_discovery_keeps_assembly_accession_bundle_coherent(tmp_path):
    module = load_module()
    input_dir = tmp_path / "NCBI" / "species_wise_original"
    species_dir = input_dir / "Test_species"
    species_dir.mkdir(parents=True)
    (species_dir / "GCA_000001.1_old_cds_from_genomic.fna").write_text(">T1\nATG\n", encoding="utf-8")
    gff_path = species_dir / "GCA_000002.1_new_genomic.gff"
    gff_path.write_text("chr1\tsrc\tCDS\t1\t3\t.\t+\t0\tID=T2\n", encoding="utf-8")
    genome_path = species_dir / "GCA_000002.1_new_genomic.fna"
    genome_path.write_text(">chr1\nATG\n", encoding="utf-8")

    tasks, warnings, errors = module.discover_tasks("ncbi", input_dir)

    assert errors == []
    assert len(tasks) == 1
    assert tasks[0]["source_bundle_id"] == "GCA_000002.1"
    assert tasks[0]["cds_path"] is None
    assert tasks[0]["gff_path"] == gff_path
    assert tasks[0]["genome_path"] == genome_path
    assert any("coherent bundle 'GCA_000002.1'" in warning for warning in warnings)


def test_single_wrong_gff_is_rejected_without_strict_mode(tmp_path):
    module = load_module()
    cds_path = tmp_path / "models.cds.fa"
    gff_path = tmp_path / "models.gff3"
    cds_path.write_text(">T1\nATG\n", encoding="utf-8")
    gff_path.write_text("chr1\trepeat\trepeat_region\t1\t3\t.\t+\t.\tID=R1\n", encoding="utf-8")
    task = {
        "provider": "direct",
        "species_key": "Test_species",
        "species_prefix": "Test_species",
        "cds_path": cds_path,
        "gff_path": gff_path,
        "gene_grouping_mode": "strict",
    }
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    with pytest.raises(ValueError, match="unexpected_unmapped=1 ambiguous=0"):
        module.format_cds(task, output_dir, overwrite=False, dry_run=False)


def test_species_summary_is_incremental_and_persistent_across_runs(tmp_path):
    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    species_summary = tmp_path / "gg_input_generation_species.tsv"

    first = run_script(
        "--provider",
        "ensemblplants",
        "--input-dir",
        str(SMALL_DATASET_ROOT / "20230216_EnsemblPlants" / "original_files"),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
        "--species-summary-output",
        str(species_summary),
    )
    assert first.returncode == 0, first.stderr + "\n" + first.stdout
    with open(species_summary, "rt", encoding="utf-8", newline="") as handle:
        rows_first = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows_first) == 1
    assert rows_first[0]["provider"] == "ensemblplants"
    assert rows_first[0]["species_prefix"] == "Ostreococcus_lucimarinus"
    assert int(rows_first[0]["cds_sequences_before"]) >= int(rows_first[0]["cds_sequences_after"])
    assert rows_first[0]["cds_first_sequence_name"] != ""

    second = run_script(
        "--provider",
        "phycocosm",
        "--input-dir",
        str(SMALL_DATASET_ROOT / "PhycoCosm" / "species_wise_original"),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
        "--species-summary-output",
        str(species_summary),
    )
    assert second.returncode == 0, second.stderr + "\n" + second.stdout
    with open(species_summary, "rt", encoding="utf-8", newline="") as handle:
        rows_second = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows_second) == 2
    by_key = {(row["provider"], row["species_prefix"]): row for row in rows_second}
    assert ("ensemblplants", "Ostreococcus_lucimarinus") in by_key
    assert ("phycocosm", "Microglena_spYARC") in by_key
    assert by_key[("phycocosm", "Microglena_spYARC")]["cds_first_sequence_name"] != ""


def test_format_species_inputs_derives_from_gbff_and_genome_when_gff_and_cds_are_missing(tmp_path):
    input_dir = tmp_path / "Direct" / "species_wise_original"
    species_dir = input_dir / "Arabidopsis_thaliana"
    species_dir.mkdir(parents=True, exist_ok=True)
    gbff_path = species_dir / "Arabidopsis_thaliana.genomic.gbff"
    genome_path = species_dir / "Arabidopsis_thaliana.genome.fa"
    gbff_path.write_text(
        "\n".join(
            [
                "LOCUS       chr1               9 bp    DNA     linear   PLN 01-JAN-2000",
                "DEFINITION  test.",
                "ACCESSION   chr1",
                "VERSION     chr1",
                "FEATURES             Location/Qualifiers",
                "     gene            1..9",
                '                     /locus_tag="gene1"',
                '                     /gene="gene1"',
                "     CDS             join(1..3,7..9)",
                '                     /locus_tag="gene1"',
                '                     /gene="gene1"',
                '                     /protein_id="gene1.t1"',
                "ORIGIN",
                "        1 atgaaattt",
                "//",
                "",
            ]
        ),
        encoding="utf-8",
    )
    genome_path.write_text(">chr1\nATGAAATTT\n", encoding="utf-8")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    species_summary = tmp_path / "gg_input_generation_species.tsv"
    completed = run_script(
        "--provider",
        "direct",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
        "--species-summary-output",
        str(species_summary),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout

    formatted_cds = out_cds / "Arabidopsis_thaliana_genomic.derived.cds.fa.gz"
    formatted_gff = out_gff / "Arabidopsis_thaliana_genomic.derived.gff.gz"
    formatted_genome = out_genome / "Arabidopsis_thaliana_genome.fa.gz"
    assert formatted_cds.exists()
    assert formatted_gff.exists()
    assert formatted_genome.exists()

    with gzip.open(formatted_cds, "rt", encoding="utf-8") as handle:
        cds_text = handle.read()
    assert cds_text.count(">Arabidopsis_thaliana_gene1") == 1
    assert "ATGTTT" in cds_text

    with gzip.open(formatted_gff, "rt", encoding="utf-8") as handle:
        gff_text = handle.read()
    assert "##gff-version 3" in gff_text
    assert "\tgene\t" in gff_text
    assert "\tmRNA\t" in gff_text
    assert "\tCDS\t" in gff_text

    with open(species_summary, "rt", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    row = rows[0]
    assert str(gbff_path) in row["cds_input_path"]
    assert "derived CDS" in row["cds_input_path"]
    assert str(gbff_path) in row["gff_input_path"]
    assert "derived GFF" in row["gff_input_path"]


def test_format_species_inputs_derives_cds_gff_and_genome_from_embl(tmp_path):
    input_dir = tmp_path / "Direct" / "species_wise_original"
    species_dir = input_dir / "Picea_abies"
    species_dir.mkdir(parents=True, exist_ok=True)
    embl_path = species_dir / "Picea_abies.embl"
    embl_path.write_text(
        """ID   chr1; SV 1; linear; genomic DNA; STD; PLN; 9 BP.
XX
AC   chr1;
XX
DE   test
XX
OS   Picea abies
XX
FH   Key             Location/Qualifiers
FH
FT   source          1..9
FT   gene            1..9
FT                   /locus_tag="gene1"
FT                   /gene="gene1"
FT   CDS             join(1..3,7..9)
FT                   /locus_tag="gene1"
FT                   /gene="gene1"
FT                   /protein_id="gene1.t1"
XX
SQ   Sequence 9 BP; 3 A; 1 C; 1 G; 4 T; 0 other;
     atgaaattt                                                               9
//
""",
        encoding="utf-8",
    )
    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    completed = run_script(
        "--provider",
        "direct",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
    )
    assert completed.returncode == 0, completed.stderr + "\n" + completed.stdout
    assert len(list(out_cds.glob("*.fa.gz"))) == 1
    assert len(list(out_gff.glob("*.gff.gz"))) == 1
    assert len(list(out_genome.glob("*.fa.gz"))) == 1
    with gzip.open(next(out_cds.glob("*.fa.gz")), "rt", encoding="utf-8") as handle:
        cds_text = handle.read()
    assert ">Picea_abies_gene1" in cds_text
    assert "ATGTTT" in cds_text
    with gzip.open(next(out_genome.glob("*.fa.gz")), "rt", encoding="utf-8") as handle:
        assert "ATGAAATTT" in handle.read()
    validation = run_validate_mapping_script(
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
    )
    assert validation.returncode == 0, validation.stderr + "\n" + validation.stdout


def test_annotation_only_gff_does_not_load_genome(tmp_path, monkeypatch):
    from format_species_annotation import genbank

    gff = tmp_path / "assembly.gff"
    gff.write_text("##gff-version 3\nchr1\t.\tregion\t1\t9\t.\t+\t.\tID=chr1\n")

    def unexpected_load(path):
        pytest.fail("A GFF without nuclear CDS must not load the genome")

    monkeypatch.setattr(genbank, "load_genome_sequences", unexpected_load)
    task = {"provider": "direct", "species_key": "Species_a", "species_prefix": "Species_a",
            "gff_path": gff, "genome_path": tmp_path / "huge.fa"}
    assert list(genbank.derive_cds_records_from_gff_and_genome(task)) == []


@pytest.mark.parametrize("phase", [".", "0"])
def test_coge_duplicate_complete_models_are_read_once(tmp_path, phase):
    module = load_module()
    genome = tmp_path / "sequence.fa"
    genome.write_text(">chr1\nATGCCCTAA\n")
    gff = tmp_path / "duplicate.gff"
    gff.write_text("\n".join(
        f"chr1\tCoGe\tCDS\t{start}\t{end}\t.\t+\t{phase}\t"
        f"ID=gene1{suffix};Name=gene1;CDS=gene1;coge_fid={fid}"
        for fid in (101, 102)
        for start, end, suffix in ((1, 3, ""), (7, 9, ".CDS2"))
    ) + "\n")
    task = {"provider": "coge", "species_key": "Alpha_alba", "species_prefix": "Alpha_alba", "gff_path": gff, "genome_path": genome}
    records = list(module.derive_cds_records_from_gff_and_genome(task))
    assert len(records) == 1
    assert records[0][1] == "ATGTAA"
    from format_species_annotation.gff_repair import iter_repaired_gff_lines
    counters = {"changed_lines": 0, "changed_values": 0, "changed_references": 0, "normalized_bare_attribute_lines": 0}
    lines = list(iter_repaired_gff_lines(gff, {}, counters, coge=True))
    assert len(lines) == 2
    assert counters["duplicate_coge_cds_blocks_removed"] == 2


def test_coge_duplicate_model_phase_conflict_still_fails(tmp_path):
    module = load_module()
    gff = tmp_path / "duplicate.gff"
    gff.write_text("\n".join(
        f"chr1\tCoGe\tCDS\t1\t3\t.\t+\t{phase}\tID=gene1;Name=gene1;CDS=gene1;coge_fid={fid}"
        for fid, phase in ((101, "0"), (102, "1"))
    ) + "\n")
    with pytest.raises(ValueError, match="Conflicting CoGe CDS feature identity"):
        list(module.derive_cds_records_from_gff_and_genome({"provider": "coge", "species_key": "Alpha_alba", "gff_path": gff, "genome_path": tmp_path / "missing.fa"}))
    from format_species_annotation.gff_repair import iter_repaired_gff_lines
    with pytest.raises(ValueError, match="Conflicting CoGe CDS feature identity"):
        list(iter_repaired_gff_lines(gff, {}, {}, coge=True))


def test_format_species_inputs_does_not_write_empty_gbff_derived_outputs(tmp_path):
    input_dir = tmp_path / "Direct" / "species_wise_original"
    species_dir = input_dir / "Fakus_emptyus"
    species_dir.mkdir(parents=True, exist_ok=True)
    gbff_path = species_dir / "Fakus_emptyus.genomic.gbff"
    genome_path = species_dir / "Fakus_emptyus.genome.fa"
    gbff_path.write_text(
        "\n".join(
            [
                "LOCUS       chr1               9 bp    DNA     linear   PLN 01-JAN-2000",
                "DEFINITION  test.",
                "ACCESSION   chr1",
                "VERSION     chr1",
                "FEATURES             Location/Qualifiers",
                "     gene            1..9",
                '                     /locus_tag="gene1"',
                '                     /gene="gene1"',
                "ORIGIN",
                "        1 atgaaattt",
                "//",
                "",
            ]
        ),
        encoding="utf-8",
    )
    genome_path.write_text(">chr1\nATGAAATTT\n", encoding="utf-8")

    out_cds = tmp_path / "species_cds"
    out_gff = tmp_path / "species_gff"
    out_genome = tmp_path / "species_genome"
    species_summary = tmp_path / "gg_input_generation_species.tsv"
    completed = run_script(
        "--provider",
        "direct",
        "--input-dir",
        str(input_dir),
        "--species-cds-dir",
        str(out_cds),
        "--species-gff-dir",
        str(out_gff),
        "--species-genome-dir",
        str(out_genome),
        "--species-summary-output",
        str(species_summary),
    )
    assert completed.returncode == 2
    assert "derived CDS contained no records" in completed.stderr
    assert "derived GFF contained no feature rows" in completed.stderr
    assert "CDS=Fakus_emptyus.genome.fa (derived CDS) (empty, NA" in completed.stdout
    assert "GFF=Fakus_emptyus.genomic.gbff (derived GFF) (empty, lines=1)" in completed.stdout
    assert not (out_cds / "Fakus_emptyus_genomic.derived.cds.fa.gz").exists()
    assert not (out_gff / "Fakus_emptyus_genomic.derived.gff.gz").exists()
    assert (out_genome / "Fakus_emptyus_genome.fa.gz").exists()
    assert species_summary.exists()
    with open(species_summary, "rt", encoding="utf-8", newline="") as handle:
        assert len(list(csv.DictReader(handle, delimiter="\t"))) == 0


def test_ensembl_like_gff_selection_prefers_full_annotation_over_partial_files():
    module = load_module()
    candidates = [
        "https://ftp.example/Saccharomyces_cerevisiae.R64-1-1.115.abinitio.gff3.gz",
        "https://ftp.example/Saccharomyces_cerevisiae.R64-1-1.115.chromosome.I.gff3.gz",
        "https://ftp.example/Saccharomyces_cerevisiae.R64-1-1.115.gff3.gz",
    ]
    assert module.select_best_url_for_label("ensembl", "GFF", candidates).endswith(".115.gff3.gz")

    candidates = [
        "https://ftp.example/Penaeus_japonicus_gca017312705v1.Mj_TUMSAT_v1.0.62.chr.gff3.gz",
        "https://ftp.example/Penaeus_japonicus_gca017312705v1.Mj_TUMSAT_v1.0.62.primary_assembly.NC_007010.1.gff3.gz",
        "https://ftp.example/Penaeus_japonicus_gca017312705v1.Mj_TUMSAT_v1.0.62.gff3.gz",
    ]
    assert module.select_best_url_for_label("ensemblmetazoa", "GFF", candidates).endswith(".62.gff3.gz")


def test_ensembl_like_explicit_partial_gff_url_is_replaced_with_full_annotation(monkeypatch):
    module = load_module()
    provider_module = sys.modules[module.resolve_preferred_ensembl_like_gff_url.__module__]

    partial_url = "https://ftp.example/current_gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.115.abinitio.gff3.gz"
    full_url = "https://ftp.example/current_gff3/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.115.gff3.gz"

    def fake_resolve_urls_from_index_url(provider, index_url, timeout, headers):
        assert provider == "ensembl"
        assert index_url == "https://ftp.example/current_gff3/saccharomyces_cerevisiae/"
        return {"gff_url": full_url}

    monkeypatch.setattr(provider_module, "resolve_urls_from_index_url", fake_resolve_urls_from_index_url)

    assert module.resolve_preferred_ensembl_like_gff_url("ensembl", partial_url, 1, {}) == full_url


def test_remove_stale_ensembl_like_partial_gff_outputs_keeps_full_annotation(tmp_path):
    module = load_module()
    gff_dir = tmp_path / "species_gff"
    gff_dir.mkdir()
    full = gff_dir / "Saccharomyces_cerevisiae_R64-1-1.115.gff.gz"
    abinitio = gff_dir / "Saccharomyces_cerevisiae_R64-1-1.115.abinitio.gff.gz"
    chromosome = gff_dir / "Saccharomyces_cerevisiae_R64-1-1.115.chromosome.I.gff.gz"
    unrelated = gff_dir / "Drosophila_melanogaster_BDGP6.54.115.abinitio.gff.gz"
    longer_species = gff_dir / "Saccharomyces_cerevisiae_subsp_x_R64-1-1.115.abinitio.gff.gz"
    for path in (full, abinitio, chromosome, unrelated, longer_species):
        path.write_bytes(b"dummy")
    abinitio_audit = Path(str(abinitio) + ".repair.json")
    chromosome_audit = Path(str(chromosome) + ".repair.json")
    abinitio_audit.write_text("{}", encoding="utf-8")
    chromosome_audit.write_text("{}", encoding="utf-8")

    removed = module.remove_stale_ensembl_like_partial_gff_outputs(
        "ensembl",
        "Saccharomyces_cerevisiae",
        gff_dir,
        full,
    )

    assert sorted(removed) == sorted([abinitio.name, chromosome.name])
    assert full.exists()
    assert unrelated.exists()
    assert longer_species.exists()
    assert not abinitio.exists()
    assert not chromosome.exists()
    assert not abinitio_audit.exists()
    assert not chromosome_audit.exists()


@pytest.mark.parametrize('strands,expected', [(['-', '-'], 'CATGGG'), (['+', '-'], 'ATGGGG')])
def test_derive_explicit_trans_splicing_without_dropping_or_reordering(tmp_path, monkeypatch, strands, expected):
    monkeypatch.syspath_prepend(str(SCRIPT_PATH.parent))
    from format_species_annotation.genbank import derive_cds_records_from_gff_and_genome
    gff = tmp_path/'annotation.gff3'
    genome = tmp_path/'genome.fa'
    genome.write_text('>chr\nATGAAACCC\n')
    gff.write_text('chr\tsrc\tgene\t1\t9\t.\t?\t.\tID=g\n' +
        'chr\tsrc\tmRNA\t1\t9\t.\t?\t.\tID=t;Parent=g\n' +
        f'chr\tsrc\tCDS\t7\t9\t.\t{strands[1]}\t0\tParent=t;exception=trans-splicing;part=2\n' +
        f'chr\tsrc\tCDS\t1\t3\t.\t{strands[0]}\t0\tParent=t;exception=trans-splicing;part=1\n')
    records = list(derive_cds_records_from_gff_and_genome(dict(
        provider='direct',species_key='Species_a',gff_path=gff,genome_path=genome,gene_grouping_mode='strict')))
    assert len(records)==1 and records[0][0].split()[0]=='t' and records[0][1]==expected


@pytest.mark.parametrize('provided', [False, True])
@pytest.mark.parametrize('mode', ['strict', 'rescue_overlap'])
@pytest.mark.parametrize('stable', [False, True])
def test_declared_gene_boundaries_survive_overlap_and_reused_labels(tmp_path, provided, mode, stable):
    mod = load_module()
    gff, genome, cds = tmp_path/'source.gff', tmp_path/'genome.fa', tmp_path/'cds.fa'
    genome.write_text('>chr1\nATGATGATG\n')
    label = ';locus_tag=shared' if stable else ''
    gff.write_text('chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=g1'+label+'\n'
                  'chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1\n'
                  'chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=t1\n'
                  'chr1\tsrc\tgene\t4\t9\t.\t+\t.\tID=g2'+label+'\n'
                  'chr1\tsrc\tmRNA\t4\t9\t.\t+\t.\tID=t2;Parent=g2\n'
                  'chr1\tsrc\tCDS\t4\t9\t.\t+\t0\tID=c2;Parent=t2\n')
    task=dict(provider='direct', species_key='Test_species', species_prefix='Test_species',
              gff_path=gff, genome_path=genome, gene_grouping_mode=mode, format_strict=True)
    if provided:
        cds.write_text('>t1\nATGATGATG\n>t2\nATGATG\n')
        task['cds_path']=cds
    result=mod.format_cds(task,tmp_path,False,False)
    assert gzip.open(result['output_path'],'rt').read()=='>Test_species_g1\nATGATGATG\n>Test_species_g2\nATGATG\n'


@pytest.mark.parametrize('provided', [False, True])
@pytest.mark.parametrize('root', ['Lavan.20G002400','Lavan.S003640'])
def test_numeric_author_model_suffix_requires_coding_locus(tmp_path, provided, root):
    mod=load_module()
    gff=tmp_path/'source.gff'
    genome=tmp_path/'genome.fa'
    genome.write_text('>chr1\nATGATGATG\n')
    gff.write_text(f'chr1\ts\tmRNA\t1\t9\t.\t+\t.\tID={root}.1\n'
                  f'chr1\ts\tCDS\t1\t9\t.\t+\t0\tParent={root}.1\n'
                  f'chr1\ts\tmRNA\t1\t6\t.\t+\t.\tID={root}.2\n'
                  f'chr1\ts\tCDS\t1\t6\t.\t+\t0\tParent={root}.2\n')
    task=dict(provider='direct',species_key='Test_species',species_prefix='Test_species',gff_path=gff,genome_path=genome)
    if provided:
        cds=tmp_path/'cds.fa'
        cds.write_text(f'>{root}.1\nATGATGATG\n>{root}.2\nATGATG\n')
        task['cds_path']=cds
    result=mod.format_cds(task,tmp_path,False,False)
    assert gzip.open(result['output_path'],'rt').read()==f'>Test_species_{root}\nATGATGATG\n'
    assert mod.collapse_transcript_suffix('direct','Pn1.1301')=='Pn1.1301'
    assert mod.collapse_transcript_suffix('direct',root+'.1')==root+'.1'


def test_derived_cds_uses_parent_gene_identity_before_display_symbol(tmp_path):
    mod=load_module()
    gff=tmp_path/'source.gff'
    genome=tmp_path/'genome.fa'
    genome.write_text('>chr1\nATGATGATGATG\n')
    gff.write_text('chr1\ts\tgene\t1\t12\t.\t+\t.\tID=g1;Name=display\n'
                  'chr1\ts\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1\n'
                  'chr1\ts\tCDS\t1\t9\t.\t+\t0\tParent=t1;gene=first\n'
                  'chr1\ts\tmRNA\t1\t12\t.\t+\t.\tID=t2;Parent=g1\n'
                  'chr1\ts\tCDS\t1\t12\t.\t+\t0\tParent=t2;gene=second\n')
    result=mod.format_cds(dict(provider='direct',species_key='Test_species',species_prefix='Test_species',
        gff_path=gff,genome_path=genome),tmp_path,False,False)
    assert gzip.open(result['output_path'],'rt').read()=='>Test_species_g1\nATGATGATGATG\n'


def test_transcript_cds_spans_choose_longest_cds_and_record_utr_removal(tmp_path):
    mod=load_module()
    gff=tmp_path/'source.gff'
    cds=tmp_path/'mrna.fa'
    gff.write_text('chr1\ts\tgene\t1\t39\t.\t+\t.\tID=g1\n'
                  'chr1\ts\tmRNA\t1\t21\t.\t+\t.\tID=t1;Parent=g1\n'
                  'chr1\ts\tCDS\t4\t12\t.\t+\t0\tParent=t1\n'
                  'chr1\ts\tmRNA\t22\t39\t.\t+\t.\tID=t2;Parent=g1\n'
                  'chr1\ts\tCDS\t22\t36\t.\t+\t0\tParent=t2\n')
    cds.write_text('>t1 CDS=4-12\nAAAATGATGATGTTTTTTTTT\n>t2 CDS=1-15\nATGATGATGATGATGTTT\n')
    result=mod.format_cds(dict(provider='direct',species_key='Test_species',species_prefix='Test_species',
        gff_path=gff,cds_path=cds),tmp_path,False,False)
    assert gzip.open(result['output_path'],'rt').read()=='>Test_species_g1\nATGATGATGATGATG\n'
    audit=json.loads(Path(str(result['output_path'])+'.gff-grouping.json').read_text())
    assert len(audit['rna_conversion']['trimmed'])==2


@pytest.mark.parametrize('strand,expected', [('+','CCCGGTTTA'),('-','ACCCGGTTT')])
def test_gwh_rna_cds_extraction_respects_splicing_and_orientation(tmp_path,strand,expected):
    from format_species_annotation.rna import extract_input_cds
    gff=tmp_path/'source.gff'
    gff.write_text(f'chr1\ts\tmRNA\t1\t18\t.\t{strand}\t.\tID=r1;Accession=GWHT1\n'
                  f'chr1\ts\texon\t1\t6\t.\t{strand}\t.\tParent=r1\n'
                  f'chr1\ts\texon\t10\t18\t.\t{strand}\t.\tParent=r1\n'
                  f'chr1\ts\tCDS\t4\t6\t.\t{strand}\t0\tParent=r1\n'
                  f'chr1\ts\tCDS\t11\t16\t.\t{strand}\t0\tParent=r1\n')
    assert extract_input_cds(dict(gff_path=gff),'GWHT1 Type=mRNA','AAACCCGGGTTTAAA')==expected


def test_noncoding_rna_is_recorded_and_not_retained_as_an_unmapped_gene(tmp_path):
    mod=load_module()
    gff=tmp_path/'source.gff'
    cds=tmp_path/'mrna.fa'
    gff.write_text('chr1\ts\tgene\t1\t15\t.\t+\t.\tID=g1\n'
                  'chr1\ts\tmRNA\t1\t9\t.\t+\t.\tID=r1;Accession=GWHT1;Parent=g1\n'
                  'chr1\ts\texon\t1\t9\t.\t+\t.\tParent=r1\n'
                  'chr1\ts\tCDS\t1\t9\t.\t+\t0\tParent=r1\n'
                  'chr1\ts\tmRNA\t1\t15\t.\t+\t.\tID=r2;Accession=GWHT2;Parent=g1\n'
                  'chr1\ts\texon\t1\t15\t.\t+\t.\tParent=r2\n')
    cds.write_text('>GWHT1 Type=mRNA\nATGATGATG\n>GWHT2 Type=mRNA\nATGATGATGATGATG\n')
    result=mod.format_cds(dict(provider='direct',species_key='Test_species',species_prefix='Test_species',
        gff_path=gff,cds_path=cds,format_strict=True),tmp_path,False,False)
    assert gzip.open(result['output_path'],'rt').read()=='>Test_species_g1\nATGATGATG\n'
    audit=json.loads(Path(str(result['output_path'])+'.gff-grouping.json').read_text())
    assert audit['rna_conversion']['excluded_noncoding']==['GWHT2']


@pytest.mark.parametrize("provided", [False, True])
def test_btu_author_notes_restore_genes_and_transcript_parents(tmp_path, provided):
    mod = load_module()
    gff, genome, cds = (tmp_path / name for name in ("models.gff", "genome.fa", "cds.fa"))
    genome.write_text(">chr1\nATGATGATGATGATGATG\n")
    rows = []
    for index, (start, end, author) in enumerate(((1, 9, "Btu.g00001t01"),
                                                (1, 6, "Btu.g00001t02"),
                                                (10, 18, "Btu.g00002t01")), 1):
        rows.extend((
            f"chr1\tEMBL\tgene\t{start}\t{end}\t.\t+\t.\tID=gene-BTU1;locus_tag=BTU1",
            f"chr1\tEMBL\tmRNA\t{start}\t{end}\t.\t+\t.\tID=rna-BTU1;Parent=gene-BTU1;Note=ID:{author}%3B~source:EVM;locus_tag=BTU1",
            f"chr1\tEMBL\texon\t{start}\t{end}\t.\t+\t.\tParent=rna-BTU1;Note=ID:{author};locus_tag=BTU1",
            f"chr1\tEMBL\tCDS\t{start}\t{end}\t.\t+\t0\tID=cds-P{index};protein_id=P{index};Parent=gene-BTU1;Note=ID:{author}.CDS;locus_tag=BTU1",
        ))
    gff.write_text("\n".join(rows) + "\n")
    original = gff.read_bytes()
    task = dict(provider="direct", species_key="Test_species", species_prefix="Test_species",
                gff_path=gff, genome_path=genome, format_strict=True)
    if provided:
        cds.write_text(">P1\nATGATGATG\n>P2\nATGATG\n>P3\nATGATGATG\n")
        task["cds_path"] = cds
    result = mod.format_cds(task, tmp_path, False, False)
    assert gzip.open(result["output_path"], "rt").read() == (
        ">Test_species_Btu.g00001\nATGATGATG\n>Test_species_Btu.g00002\nATGATGATG\n")
    repaired = mod.format_gff(task, tmp_path, False, False, formatted_cds_path=result["output_path"])
    text = gzip.open(repaired["output_path"], "rt").read()
    assert "Parent=Btu.g00001t01" in text
    assert "orig_export_parent=gene-BTU1" in text
    assert "##genegalleon-original-id-normalization" in text
    assert gff.read_bytes() == original


@pytest.mark.parametrize("failure", ["missing_note", "outside_rna", "conflicting_rna"])
def test_btu_author_normalization_rejects_unproven_models(tmp_path, failure):
    from format_species_annotation.source_identity import source_annotation_path
    gff = tmp_path / "source.gff"
    rna = "chr1\tEMBL\tmRNA\t1\t9\t.\t+\t.\tID=r;Note=ID:Btu.g00001t01\n"
    note = "" if failure == "missing_note" else ";Note=ID:Btu.g00001t01.CDS"
    end = 12 if failure == "outside_rna" else 9
    gff.write_text(rna + f"chr1\tEMBL\tCDS\t1\t{end}\t.\t+\t0\tID=c;Parent=r{note}\n" +
                   (rna.replace("\t1\t9\t", "\t1\t12\t") if failure == "conflicting_rna" else ""))
    with pytest.raises(ValueError):
        source_annotation_path(gff)


def test_derived_cds_rejects_ambiguous_gene_parents(tmp_path):
    mod = load_module()
    gff, genome = tmp_path / "source.gff", tmp_path / "genome.fa"
    genome.write_text(">chr1\nATGATGATG\n")
    gff.write_text("chr1\ts\tgene\t1\t9\t.\t+\t.\tID=g1\n"
                   "chr1\ts\tgene\t1\t9\t.\t+\t.\tID=g2\n"
                   "chr1\ts\tmRNA\t1\t9\t.\t+\t.\tID=t;Parent=g1,g2\n"
                   "chr1\ts\tCDS\t1\t9\t.\t+\t0\tParent=t\n")
    with pytest.raises(ValueError, match="Ambiguous GFF gene parents"):
        list(mod.derive_cds_records_from_gff_and_genome(dict(
            provider="direct", species_key="Test_species", gff_path=gff, genome_path=genome)))


@pytest.mark.parametrize("header", ["t CDS=0-9", "t CDS=1-12", "t Type=mRNA"])
def test_rna_conversion_rejects_invalid_cds_or_missing_identity(header):
    from format_species_annotation.rna import extract_input_cds
    with pytest.raises(ValueError):
        extract_input_cds({}, header, "ATGATGATG")


@pytest.mark.parametrize("author", ["Lavan.01G000100", "Other.01G000100"])
def test_lavan_original_gene_number_restores_disjoint_isoforms_only_in_known_export(tmp_path, author):
    mod = load_module()
    gff, cds = tmp_path / "source.gff", tmp_path / "cds.fa"
    gff.write_text(f"chr1\ts\tmRNA\t1\t9\t.\t+\t.\tID={author}.1\n"
                   f"chr1\ts\tCDS\t1\t9\t.\t+\t0\tParent={author}.1\n"
                   f"chr1\ts\tmRNA\t20\t25\t.\t+\t.\tID={author}.2\n"
                   f"chr1\ts\tCDS\t20\t25\t.\t+\t0\tParent={author}.2\n")
    cds.write_text(f">{author}.1\nATGATGATG\n>{author}.2\nATGATG\n")
    result = mod.format_cds(dict(provider="direct", species_key="Test_species", species_prefix="Test_species",
                                 gff_path=gff, cds_path=cds), tmp_path, False, False)
    assert result["after_count"] == (1 if author.startswith("Lavan.") else 2)


def test_source_author_normalization_is_idempotent(tmp_path):
    from format_species_annotation.source_identity import source_annotation_path
    gff = tmp_path / "source.gff"
    gff.write_text("chr1\ts\tmRNA\t1\t9\t.\t+\t.\tID=Lavan.S000100.1\n"
                   "chr1\ts\tCDS\t1\t9\t.\t+\t0\tParent=Lavan.S000100.1\n")
    normalized = source_annotation_path(gff)
    assert source_annotation_path(normalized) == normalized


def test_embedded_gff_fasta_never_enters_feature_or_encoding_parsing(tmp_path):
    mod = load_module()
    gff, cds = tmp_path / "source.gff", tmp_path / "cds.fa"
    annotation = ("##gff-version 3\nchr1\ts\tgene\t1\t9\t.\t+\t.\tID=g1\n"
                  "chr1\ts\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1\n"
                  "chr1\ts\tCDS\t1\t9\t.\t+\t0\tParent=t1\n")
    gff.write_bytes(annotation.encode() + b"##FASTA\n>chr1\n" + b"A" * 1_000_000 + b"\xff\n")
    cds.write_text(">t1\nATGATGATG\n")
    task = dict(provider="direct", species_key="Test_species", species_prefix="Test_species",
                gff_path=gff, cds_path=cds)
    result = mod.format_cds(task, tmp_path, False, False)
    repaired = mod.format_gff(task, tmp_path, False, False, formatted_cds_path=result["output_path"])
    text = gzip.open(repaired["output_path"], "rt").read()
    assert "##FASTA" not in text and ">chr1" not in text
    assert repaired["invalid_utf8_bytes"] == 0
    assert text.count("\tCDS\t") == 1


def test_gff_gene_repair_detects_alignment_collisions_before_and_after_genes(tmp_path):
    from workflow.support.format_species_annotation.gff_repair import choose_gene_id_repairs

    gff = tmp_path / "source.gff"
    gff.write_text(
        "chr1\tsrc\tmatch_part\t1\t3\t.\t+\t.\tID=targetTaken\n"
        "chr1\tsrc\tgene\t1\t3\t.\t+\t.\tID=rawTaken;Alias=targetTaken\n"
        "chr1\tsrc\tgene\t4\t6\t.\t+\t.\tID=rawShared;Alias=targetShared\n"
        "chr1\tsrc\tgene\t7\t9\t.\t+\t.\tID=rawSafe;Alias=targetSafe\n"
        "chr1\tsrc\tprotein_match\t4\t6\t.\t+\t.\tID=rawShared\n"
        + "".join(f"chr1\tsrc\tmatch_part\t1\t3\t.\t+\t.\tID=unrelated{i}\n" for i in range(1000))
    )
    plan = choose_gene_id_repairs(gff, {"targetTaken", "targetShared", "targetSafe"})
    assert plan["id_mapping"] == {"rawSafe": "targetSafe"}
    assert {row["reason"] for row in plan["collisions"]} == {
        "target_id_already_exists", "source_id_is_shared_with_non_gene_feature"}


@pytest.mark.parametrize('provided', [False, True])
@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('prefix,encoded', [('', False), ('evm.model.', False), ('evm_27.model.', False), ('evm.model.', True)])
def test_overlap_rescue_keeps_existing_isoform_group_atomic(tmp_path, provided, strand, prefix, encoded):
    mod = load_module()
    gff, genome, cds = tmp_path/'source.gff', tmp_path/'genome.fa', tmp_path/'cds.fa'
    genome.write_text('>chr1\n'+'ATG'*30+'\n')
    parts = {'a.t1': [(1, 9), (31, 39)],
             'b.t1': [(1, 9), (21, 23), (31, 39), (61, 69)],
             'b.t2': [(1, 9), (31, 39), (61, 69)]}
    parts = {prefix + name: intervals for name, intervals in parts.items()}
    lines = []
    for transcript, intervals in parts.items():
        lines.append(f'chr1\ts\ttranscript\t1\t69\t.\t{strand}\t.\tID={transcript}')
        for start, end in intervals:
            lines.append(f'chr1\ts\tCDS\t{start}\t{end}\t.\t{strand}\t0\tParent={transcript}')
    gff.write_text('\n'.join(lines)+'\n')
    if encoded:
        gff.write_text(gff.read_text().replace(prefix, prefix.replace('.', '%2E')))
    task = dict(provider='direct', species_key='Test_species', species_prefix='Test_species',
                gff_path=gff, genome_path=genome, gene_grouping_mode='rescue_overlap', format_strict=True)
    if provided:
        unit = 'ATG' if strand == '+' else 'CAT'
        cds.write_text(''.join(f'>{name}\n'+unit*(sum(e-s+1 for s,e in intervals)//3)+'\n'
                               for name, intervals in parts.items()))
        task['cds_path'] = cds
    result = mod.format_cds(task, tmp_path, False, False)
    text = gzip.open(result['output_path'], 'rt').read()
    assert result['after_count'] == 1
    assert text.splitlines()[0] == '>Test_species_a'
    assert len(text.splitlines()[1]) == 30
    formatted = mod.format_gff(task, tmp_path, False, False, formatted_cds_path=result['output_path'])
    sys.path.insert(0, str(SCRIPT_PATH.parent))
    import gff2genestat as reader
    rows = reader.process_single_gff(formatted['output_path'].name, str(tmp_path), ['Test_species_a'],
                                    'CDS', 'longest',
                                    ['sequence','source','feature','start','end','score','strand','phase','attributes'],
                                    ['gene_id','feature_size','feature_blocks','chromosome','start','end','strand','feature_type'])
    assert rows.gene_id.tolist() == ['Test_species_a']
    assert rows.feature_size.tolist() == [30]
    original_rows = [line.split('\t') for line in gff.read_text().splitlines()]
    output_rows = [line.split('\t') for line in gzip.open(formatted['output_path'], 'rt').read().splitlines()]
    assert [r[:8] for r in output_rows] == [r[:8] for r in original_rows]
    for original, output in zip(original_rows, output_rows, strict=True):
        assert output[8].startswith(mod.apply_common_replacements(original[8]))


def test_gff_rescued_owner_normalisation_rejects_distinct_owner_collision(tmp_path, monkeypatch):
    load_module()
    from format_species_annotation import gff_repair
    gff, cds, output = tmp_path/'source.gff', tmp_path/'cds.fa', tmp_path/'output.gff.gz'
    gff.write_text('chr1\ts\tCDS\t1\t6\t.\t+\t0\tParent=evm.model.t\n'
                   'chr1\ts\tCDS\t20\t25\t.\t+\t0\tParent=t\n')
    cds.write_text('>Test_species_gA\nATGAAA\n>Test_species_gB\nATGCCC\n')
    # A legacy paired CDS can predate normalization; conflicting source owners
    # must be rejected before their transcript IDs become identical in the GFF.
    monkeypatch.setattr(gff_repair, 'build_gff_cds_grouping_index', lambda _task: {
        'rescued_transcript_gene_tokens': {'evm.model.t': 'gA', 't': 'gB'},
        'suffix_inferred_transcript_gene_tokens': {},
        'explicit_missing_parent_gene_tokens': (),
        'transcript_gene_tokens': {'evm.model.t': 'gA', 't': 'gB'},
    })
    task = dict(provider='direct', species_key='Test_species', species_prefix='Test_species', gff_path=gff)
    with pytest.raises(ValueError, match='Conflicting normalized rescued GFF gene owners'):
        gff_repair.write_repaired_gff(gff, cds, output, 'Test_species', 'safe', source_task=task)
    assert not output.exists()


@pytest.mark.parametrize('strict', [False, True])
def test_ncbi_unique_protein_owner_resolves_wrong_locus_at_shared_cds_coordinates(tmp_path, strict):
    mod = load_module()
    gff, cds = tmp_path/'source.gff', tmp_path/'cds.fa'
    lines = []
    for gene, protein in [('g1', 'P1.1'), ('g2', 'P2.1')]:
        lines.extend([f'chr1\ts\tgene\t1\t9\t.\t+\t.\tID=gene-{gene};locus_tag={gene}',
                      f'chr1\ts\tmRNA\t1\t9\t.\t+\t.\tID=rna-{gene};Parent=gene-{gene}',
                      f'chr1\ts\tCDS\t1\t9\t.\t+\t0\tID=cds-{protein};Parent=rna-{gene};protein_id={protein}'])
    gff.write_text('\n'.join(lines)+'\n')
    cds.write_text(''.join(f'>lcl|chr1_cds_{protein}_{i} [locus_tag=g1] [protein_id={protein}] [location=1..9] [gbkey=CDS]\nATGATGATG\n'
                           for i, protein in enumerate(['P1.1', 'P2.1'], 1)))
    task = dict(provider='ncbi', species_key='Test_species', species_prefix='Test_species',
                gff_path=gff, cds_path=cds, gene_grouping_mode='rescue_overlap', format_strict=strict)
    result = mod.format_cds(task, tmp_path, False, False)
    assert gzip.open(result['output_path'], 'rt').read() == '>Test_species_g1\nATGATGATG\n>Test_species_g2\nATGATGATG\n'
    audit = json.loads(Path(str(result['output_path'])+'.gff-grouping.json').read_text())
    assert audit['stats']['ambiguous'] == 0


@pytest.mark.parametrize("mode", ["off", "safe", "strict"])
@pytest.mark.parametrize("paired", [False, True])
def test_augustus_metadata_colon_separators_preserve_literal_values(tmp_path, mode, paired):
    module = load_module()
    from gff_attribute_syntax import validate_gff

    value = "Protein 1, other OS=Plant|Ontology_id GO:1"
    raw = tmp_path / "raw.gff3"
    rows = ["chr1\tAUGUSTUS\tgene\t1\t12\t.\t+\t.\tID=g1;Name:" + value + ";\n",
            "chr1\tAUGUSTUS\tmRNA\t1\t12\t.\t+\t.\tID=t1;Parent=g1;Blast2Go:" + value + ";\n",
            "chr1\tAUGUSTUS\tCDS\t1\t12\t.\t+\t0\tParent=t1; 5_prime_partial=true;Note=literal%253B;\n"]
    raw.write_text("".join(rows))
    cds = tmp_path / "Species_one.cds.fa"
    cds.write_text(">Species_one_g1\nATGAAACCCTAA\n")
    task = {"provider": "local", "species_key": "Species_one", "species_prefix": "Species_one",
            "gff_path": raw, "gff_repair_mode": mode}
    output = tmp_path / "out"
    output.mkdir()
    result = module.format_gff(task, output, overwrite=False, dry_run=False,
                              formatted_cds_path=cds if paired else None)
    validate_gff(result["output_path"])
    expected = "".join(rows).replace("Name:" + value, "Name=" + value.replace(",", "%2C").replace("=", "%3D"))
    expected = expected.replace("Blast2Go:" + value, "Blast2Go=" + value.replace(",", "%2C").replace("=", "%3D"))
    with gzip.open(result["output_path"], "rt") as handle:
        assert handle.read() == expected
    assert raw.read_text() == "".join(rows)


@pytest.mark.parametrize("mode", ["off", "safe", "strict"])
@pytest.mark.parametrize("paired", [False, True])
def test_funannotate_metadata_semicolons_are_repaired_before_consumption(tmp_path, mode, paired):
    module = load_module()
    from gff_attribute_syntax import file_sha256, validate_gff

    raw = tmp_path / "raw.gff3"
    rows = ["##gff-version 3\n", "# preserved\n",
            "chr1\tfunannotate\tgene\t1\t12\t.\t+\t.\tID=g1;Name=SULTR4;1_1;\n",
            "chr1\tfunannotate\tmRNA\t1\t12\t.\t+\t.\tID=t1;Parent=g1;product=Nucleosome assembly protein 1;3, variant 2;Dbxref=PFAM:PF00956;note=COG:B,COG:D;\n",
            "chr1\tfunannotate\tCDS\t1\t12\t.\t+\t0\tParent=t1,t2;Name=already%3Bescaped;Note=literal%253B;\n"]
    raw.write_text("".join(rows))
    original = raw.read_bytes()
    cds = tmp_path / "Species_one.cds.fa"
    cds.write_text(">Species_one_g1\nATGAAACCCTAA\n")
    output = tmp_path / "out"
    output.mkdir()
    task = {"provider": "local", "species_key": "Species_one", "gff_path": raw,
            "species_prefix": "Species_one", "gff_repair_mode": mode}
    result = module.format_gff(task, output, overwrite=False, dry_run=False,
                               formatted_cds_path=cds if paired else None)
    validate_gff(result["output_path"])
    with gzip.open(result["output_path"], "rt") as handle:
        actual = handle.readlines()
    expected = rows.copy()
    expected[2] = expected[2].replace("Name=SULTR4;1_1", "Name=SULTR4%3B1_1")
    expected[3] = expected[3].replace("protein 1;3,", "protein 1%3B3%2C")
    assert actual == expected
    audit = json.loads(Path(str(result["output_path"]) + ".repair.json").read_text())["attribute_syntax"]
    assert audit["changed_rows"] == 2
    assert [r["source_line"] for r in audit["changes"]] == [3, 4]
    assert audit["source_sha256"] == file_sha256(raw)
    assert audit["output_sha256"] == file_sha256(result["output_path"])
    assert raw.read_bytes() == original
    # Formatting the output again must not double escape it.
    task["gff_path"] = result["output_path"]
    another = tmp_path / "again"
    another.mkdir()
    rerun = module.format_gff(task, another, overwrite=False, dry_run=False,
                              formatted_cds_path=cds if paired else None)
    with gzip.open(rerun["output_path"], "rt") as handle:
        assert handle.readlines() == actual


@pytest.mark.parametrize("source,feature,attrs", [
    ("funannotate", "gene", "ID=g1;1;Name=ABC"),
    ("funannotate", "mRNA", "ID=t1;Parent=g1;2;product=ABC"),
    ("other", "gene", "ID=g1;Name=SULTR4;1"),
    ("funannotate", "gene", "ID=g1;Name=SULTR4;unknown"),
    ("funannotate", "gene", "ID=g1;Name=SULTR4;1;2"),
])
def test_unrecoverable_gff_attributes_fail_without_publishing(tmp_path, source, feature, attrs):
    module = load_module()
    raw = tmp_path / "invalid.gff3"
    raw.write_text(f"##gff-version 3\nchr1\t{source}\t{feature}\t1\t9\t.\t+\t.\t{attrs}\n")
    original = raw.read_bytes()
    out = tmp_path / "out"
    out.mkdir()
    task = {"gff_path": raw, "species_prefix": "Species_one"}
    with pytest.raises(ValueError, match=r"invalid.gff3:2: Unrecoverable"):
        module.format_gff(task, out, overwrite=False, dry_run=False)
    assert not list(out.iterdir()) and raw.read_bytes() == original


def test_attribute_syntax_preserves_gtf_quotes_and_ignores_embedded_fasta(tmp_path):
    load_module()
    from gff_attribute_syntax import normalise_attributes, validate_gff

    text = 'gene_id "g1"; transcript_id "t1"; product "protein 1;2";'
    assert normalise_attributes(text, "funannotate", "mRNA") == text
    gff = tmp_path / "quoted.gtf"
    gff.write_text("chr1\tsynthetic\tCDS\t1\t9\t.\t+\t0\t" + text + "\n##FASTA\n>chr1\nATGAAATAA\n")
    validate_gff(gff)


@pytest.mark.parametrize("reuse", [False, True])
def test_invalid_legacy_gff_cache_is_never_silently_reused(tmp_path, reuse):
    module = load_module()
    raw = tmp_path / "raw.gff3"
    raw.write_text("chr1\tfunannotate\tgene\t1\t9\t.\t+\t.\tID=g1;Name=SULTR4;1;\n")
    output = tmp_path / "out"
    output.mkdir()
    task = {"gff_path": raw, "species_prefix": "Species_one", "gff_repair_mode": "off"}
    cached = output / module.normalize_gff_output_basename(raw.name, "Species_one")
    with gzip.open(cached, "wt") as handle:
        handle.write(raw.read_text())
    before = cached.read_bytes()
    if reuse:
        with pytest.raises(ValueError, match="regenerate with gg_input_generation"):
            module.format_gff(task, output, overwrite=False, dry_run=False, reuse_existing=True)
        assert cached.read_bytes() == before
    else:
        result = module.format_gff(task, output, overwrite=False, dry_run=False)
        assert result["status"] == "write"
        with gzip.open(cached, "rt") as handle:
            assert "Name=SULTR4%3B1" in handle.read()


@pytest.mark.parametrize("transport", ["plain", "gz", "bz2", "tar", "tar-single"])
@pytest.mark.parametrize("compression", ["seqkit", "python"])
def test_streamed_genome_matches_record_writer_across_transport_and_edge_cases(tmp_path, monkeypatch, transport, compression):
    import bz2
    import re
    module = load_module()
    from format_species_annotation.common import extract_header_tag_value, first_token, iter_fasta_records
    from format_species_annotation.organelle import gff_organelle_seqids
    from format_species_writers import apply_common_replacements, write_fasta_records_gzip
    if compression == "python":
        import format_species_writers
        monkeypatch.setattr(format_species_writers.shutil, "which", lambda _: None)
    else:
        assert shutil.which("seqkit"), "Qualified runtime must contain real seqkit"
    contents = ("ignored before header\r\n>chr:1 description\r\na c\tgt\u2003nß\r\n\r\n"
                ">cp chloroplast\nacgt\n>other OriSeqID=cp;\nACGT\n>lcl|cp\nACGT\n"
                ">empty\n>\nacg\n>duplicate\nACT\n>duplicate\nGTT\n"
                ">long " + "header " * 150000 + "\n" + "a cgt\t" * 400000 + "\n>last\ntt")
    raw = contents.encode()
    path = tmp_path / ("Species_one.genome.fa" + {"plain": "", "gz": ".gz", "bz2": ".bz2", "tar": ".tar.gz", "tar-single": ".tar.gz"}[transport])
    if transport.startswith("tar"):
        with tarfile.open(path, "w:gz") as archive:
            member = tarfile.TarInfo("source.fa" if transport == "tar" else "single.data")
            member.size = len(raw)
            archive.addfile(member, io.BytesIO(raw))
            if transport == "tar":
                member = tarfile.TarInfo("ignored.txt")
                member.size = 1
                archive.addfile(member, io.BytesIO(b"x"))
                member = tarfile.TarInfo("second.fasta")
                second = b">second\nacgt\n"
                member.size = len(second)
                archive.addfile(member, io.BytesIO(second))
    else:
        path.write_bytes(gzip.compress(raw) if transport == "gz" else bz2.compress(raw) if transport == "bz2" else raw)
    gff = tmp_path / "source.gff"
    gff.write_text("cp\tsrc\tregion\t1\t4\t.\t+\t.\tID=cp;genome=chloroplast\n")
    organelles = gff_organelle_seqids(gff)
    expected_records = []
    for header, sequence in iter_fasta_records(path):
        record_id = first_token(apply_common_replacements(header)) or "unnamed"
        original = extract_header_tag_value(header, "OriSeqID").rstrip(";")
        if record_id in organelles or record_id.removeprefix("lcl|") in organelles or apply_common_replacements(original) in organelles:
            continue
        expected_records.append((record_id, re.sub(r"\s+", "", sequence).upper()))
    out = tmp_path / "out"
    out.mkdir()
    expected = out / "expected.fa.gz"
    write_fasta_records_gzip(expected, expected_records)
    result = module.format_genome({"species_prefix": "Species_one", "genome_path": path, "gff_path": gff}, out, True, False)
    with gzip.open(expected, "rb") as before, gzip.open(result["output_path"], "rb") as after:
        assert before.read() == after.read()
    assert result["written"] == len(expected_records)


@pytest.mark.parametrize("compression", ["seqkit", "python"])
def test_streamed_genome_read_failure_preserves_published_output(tmp_path, monkeypatch, compression):
    module = load_module()
    if compression == "python":
        import format_species_writers
        monkeypatch.setattr(format_species_writers.shutil, "which", lambda _: None)
    path = tmp_path / "Species_one.genome.fa"
    path.write_text(">chr1\nATGAAA\n")
    out = tmp_path / "out"
    out.mkdir()
    task = {"species_prefix": "Species_one", "genome_path": path}
    result = module.format_genome(task, out, True, False)
    published = result["output_path"].read_bytes()
    path.write_bytes(b">chr1\n" + b"ACGT\n" * 300000 + b"\xff")
    with pytest.raises(UnicodeDecodeError):
        module.format_genome(task, out, True, False)
    assert result["output_path"].read_bytes() == published
    assert list(out.iterdir()) == [result["output_path"]]
