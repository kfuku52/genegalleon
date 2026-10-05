"""Bounded regressions for native input retries; raw scientific inputs stay intact."""

import gzip
import json
import sys
from pathlib import Path

import pandas
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "support"))
import input_generation_stage_resume as resume
from cds_model_normalisation import CdsModelNormaliser
from format_species_annotation.gff_repair import gene_reference_repairs
from gff2genestat import read_gff_table
from gff_attribute_syntax import normalise_attributes, normalise_line


def test_description_linkage_separator_is_encoded_and_audited():
    raw = "ID=Sb01g07380;description=PREDICTED: endo-1,3;1,4-beta-D-glucanase-like [Sesamum indicum];transl_table=1"
    expected = raw.replace("1,3;1,4-", "1%2C3%3B1%2C4-")
    assert normalise_attributes(raw, "EVM", "gene") == expected
    assert normalise_attributes(expected, "EVM", "gene") == expected
    changes = []
    line = "chr1\tEVM\tgene\t1\t9\t.\t-\t.\t" + raw + "\n"
    assert normalise_line(line, "source.gff", 23978, changes).endswith(expected + "\n")
    assert changes[0]["reason"] == "escaped_description_linkage_semicolon"


@pytest.mark.parametrize(
    "attrs",
    [
        "ID=g;1,4-beta-D-glucanase",
        "ID=g;Parent=t;1,4-beta-D-glucanase",
        "description=endo-1,3;Parent:t",
        "description=endo-1,3;1,4-beta;extra",
        "description=endo;1,4-beta",
        "description=endo-1,3;ID:t",
    ],
)
def test_description_repair_never_consumes_structural_or_ambiguous_fragments(attrs):
    with pytest.raises(ValueError, match="Unrecoverable"):
        normalise_attributes(attrs, "EVM", "gene")


@pytest.mark.parametrize("mismatch", ["", "strand", "span", "mixed", "duplicate", "in_bounds"])
def test_out_of_bounds_parent_gene_requires_exact_unanimous_child_span(tmp_path, mismatch):
    (genome, gff) = (tmp_path / "genome.fa", tmp_path / "source.gff")
    genome.write_text(">wrong\n" + "A" * (90 if mismatch == "in_bounds" else 6) + "\n>right\n" + "A" * 90 + "\n")
    gene = "wrong\tEVM\tgene\t10\t30\t.\t-\t.\tID=g\n"
    rna = "right\tEVM\tmRNA\t10\t30\t.\t-\t.\tID=t;Parent=g\n"
    if mismatch == "strand":
        rna = rna.replace("\t-\t", "\t+\t")
    if mismatch == "span":
        rna = rna.replace("\t30\t", "\t29\t")
    text = gene + rna + "right\tEVM\tCDS\t10\t30\t.\t-\t0\tParent=t\n"
    if mismatch == "mixed":
        text += rna.replace("right", "wrong").replace("ID=t", "ID=u")
    if mismatch == "duplicate":
        text += gene
    gff.write_text(text)
    repairs = gene_reference_repairs(gff, genome)
    assert bool(repairs) == (mismatch == "")
    if repairs:
        assert repairs["g"]["to_seqid"] == "right"
    assert gff.read_text() == text


def annotation(alignment_count=0):
    return (
        "chr1\ts\tgene\t1\t9\t.\t+\t.\tID=g;pseudo=true\nchr1\ts\tmatch\t1\t9\t.\t+\t.\tID=bridge;Parent=g\nchr1\ts\tmRNA\t1\t9\t.\t+\t.\tID=t;Parent=bridge\nchr1\ts\texon\t1\t9\t.\t+\t.\tParent=t\nchr1\ts\tCDS\t1\t9\t.\t+\t0\tParent=t\n"
        + "".join(
            (f"chr1\ts\tmatch_part\t1\t9\t.\t+\t.\tID=alignment{i};Parent=other{i}\n" for i in range(alignment_count))
        )
    )


def test_coding_normalisation_retains_complete_ancestor_closure_and_same_decision(tmp_path):
    gff = tmp_path / "source.gff"
    gff.write_text(annotation(3000))
    args = dict(
        source=dict(gff=str(gff), genome="unused", species="S", genetic_code=1), directory=tmp_path, side="format"
    )
    (full, bounded) = (CdsModelNormaliser(**args), CdsModelNormaliser(**args, coding_only=True))
    assert len(full.features) == 3005 and len(bounded.features) == 5
    assert full.ancestry(["t"]) == bounded.ancestry(["t"]) == {"t", "bridge", "g"}
    assert full.models == bounded.models and full.exons == bounded.exons
    full.normalise("seq", "ATGAAATAA", {"t"})
    bounded.normalise("seq", "ATGAAATAA", {"t"})
    assert full.rows == bounded.rows and full.corrected == bounded.corrected


@pytest.mark.parametrize("compressed", [False, True])
def test_coding_gff_table_retains_referenced_alignment_ancestors(tmp_path, compressed):
    gff = tmp_path / ("source.gff" + (".gz" if compressed else ""))
    text = annotation(3000)
    if compressed:
        with gzip.open(gff, "wt") as handle:
            handle.write(text)
    else:
        gff.write_text(text)
    full = read_gff_table(str(gff))
    bounded = read_gff_table(str(gff), coding_only=True)
    assert len(bounded) == 5
    pandas.testing.assert_frame_equal(
        full.iloc[:5].reset_index(drop=True), bounded.sort_index().iloc[[0, 4, 1, 2, 3]].reset_index(drop=True)
    )


@pytest.fixture
def busco_pair(tmp_path):
    (donor, target) = (tmp_path / "old/output/input_generation", tmp_path / "new/output/input_generation")
    settings = []
    for root in (donor, target):
        (root / "tmp").mkdir(parents=True)
        (root / "tmp/busco_lineage.resolved.txt").write_text("embryophyta_odb12\n")
        item = dict(run_species_busco="1", busco_lineage="embryophyta_odb12")
        for kind in ("full", "short"):
            directory = root / ("busco_" + kind)
            directory.mkdir()
            item["species_busco_" + kind + "_dir"] = str(directory)
        settings.append(item)
    (source_cds, target_cds) = (donor / "cds.fa", target / "cds.fa")
    for f in (source_cds, target_cds):
        f.write_text(">S_g\nATGAAATAA\n")
    outputs = []
    for kind in ("full", "short"):
        f = Path(settings[0]["species_busco_" + kind + "_dir"]) / (
            "S.busco." + ("full.tsv" if kind == "full" else "short.txt")
        )
        f.write_text("BUSCO validated result " + kind + "\n")
        outputs.append(dict(label="busco_" + kind, path=str(f), sha256=resume.digest(f)))
    payload = dict(
        schema_version=1,
        step="input_generation_species_busco",
        family_id="S",
        parameters=dict(
            busco_lineage_request="embryophyta_odb12",
            busco_lineage_resolved="embryophyta_odb12",
            busco_mode="transcriptome",
            evalue="1e-03",
            limit="20",
        ),
        inputs=[dict(label="species_cds", path=str(source_cds), sha256=resume.digest(source_cds))],
        outputs=outputs,
    )
    (donor / "artifact_provenance").mkdir()
    (donor / "artifact_provenance/busco.S.json").write_text(json.dumps(payload))
    return (donor, target, source_cds, target_cds, *settings)


def test_identical_busco_import_writes_native_target_provenance(busco_pair):
    assert resume.import_busco(*busco_pair[:2], "S", *busco_pair[2:])
    payload = json.loads((busco_pair[1] / "artifact_provenance/busco.S.json").read_text())
    assert payload["inputs"][0]["path"] == "output/input_generation/cds.fa"
    assert all((x["path"].startswith("output/input_generation/busco_") for x in payload["outputs"]))


@pytest.mark.parametrize("change", ["lineage", "resolved", "cds", "output", "parameter"])
def test_busco_import_rejects_changed_scientific_contract(busco_pair, change):
    (donor, target, source_cds, target_cds, source_settings, target_settings) = busco_pair
    if change == "lineage":
        target_settings["busco_lineage"] = "viridiplantae_odb12"
    elif change == "resolved":
        (target / "tmp/busco_lineage.resolved.txt").write_text("viridiplantae_odb12\n")
    elif change == "cds":
        target_cds.write_text("changed")
    elif change == "output":
        next(Path(source_settings["species_busco_full_dir"]).iterdir()).write_text("changed")
    else:
        path = donor / "artifact_provenance/busco.S.json"
        payload = json.loads(path.read_text())
        payload["parameters"]["limit"] = "21"
        path.write_text(json.dumps(payload))
    if change in ("lineage", "resolved"):
        assert not resume.import_busco(donor, target, "S", source_cds, target_cds, source_settings, target_settings)
    else:
        with pytest.raises(ValueError):
            resume.import_busco(donor, target, "S", source_cds, target_cds, source_settings, target_settings)
    assert not (target / "artifact_provenance/busco.S.json").exists()


def test_parent_reference_repair_runs_through_standard_formatter_and_strict_validator(tmp_path):
    from format_species_annotation.reference import validate_gff_genome_references
    from format_species_inputs import format_cds, format_gff

    (genome, gff) = (tmp_path / "genome.fa", tmp_path / "source.gff")
    genome.write_text(">wrong\nATGAAA\n>right\nATGAAATAA\n")
    raw = "wrong\tEVM\tgene\t1\t9\t.\t+\t.\tID=g\nright\tEVM\tmRNA\t1\t9\t.\t+\t.\tID=t;Parent=g\nright\tEVM\tCDS\t1\t9\t.\t+\t0\tParent=t\n"
    gff.write_text(raw)
    task = dict(
        provider="local",
        species_prefix="Species_one",
        species_key="Species_one",
        gff_path=gff,
        genome_path=genome,
        gff_repair_mode="safe",
        gene_grouping_mode="rescue_overlap",
        format_strict=True,
    )
    cds = format_cds(task, tmp_path, False, False)
    formatted = format_gff(task, tmp_path, False, False, formatted_cds_path=cds["output_path"])
    assert validate_gff_genome_references(formatted["output_path"], genome) == 3
    with gzip.open(formatted["output_path"], "rt") as h:
        assert h.read() == raw.replace("wrong\t", "right\t")
    audit = json.loads(Path(str(formatted["output_path"]) + ".repair.json").read_text())
    assert audit["gene_reference_repairs"]["g"]["from_seqid"] == "wrong"
    assert gff.read_text() == raw


def test_busco_import_rejects_a_source_change_during_copy(busco_pair, monkeypatch):
    copy = resume.copy_atomic
    (donor, target, source_cds, target_cds, source_settings, target_settings) = busco_pair

    def racing_copy(source, destination, **kwargs):
        copy(source, destination, **kwargs)
        source_cds.write_text("raced source")

    monkeypatch.setattr(resume, "copy_atomic", racing_copy)
    with pytest.raises(OSError, match="File changed"):
        resume.import_busco(donor, target, "S", source_cds, target_cds, source_settings, target_settings)
    assert not (target / "artifact_provenance/busco.S.json").exists()
